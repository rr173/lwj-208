import uuid
import hashlib
import json
import threading
from typing import List, Dict, Optional, Tuple
from datetime import datetime
from sqlalchemy.orm import Session
from sqlalchemy import func

from app import models, schemas
from app.database import SessionLocal
from app.services.batch_service import task_manager
from app.phylogeny.distance import (
    build_sequence_matrix,
    compute_distance_matrix,
)
from app.phylogeny.tree_building import (
    build_tree,
    TreeNode,
    to_newick_format,
)
from app.phylogeny.bootstrap import compute_bootstrap_support
from app.phylogeny.molecular_clock import apply_molecular_clock
from app.phylogeny.tree_compare import compare_trees


def _get_cache_key(
    sample_ids: List[int],
    reference_name: str,
    tree_method: str,
    bootstrap_replicates: int,
    with_molecular_clock: bool
) -> str:
    sorted_ids = sorted(sample_ids)
    key_str = f"{sorted_ids}:{reference_name}:{tree_method}:{bootstrap_replicates}:{with_molecular_clock}"
    return hashlib.md5(key_str.encode()).hexdigest()


def _compute_data_hash(
    db: Session,
    sample_ids: List[int],
    reference_id: int
) -> str:
    sorted_ids = sorted(sample_ids)
    
    all_variants = db.query(
        models.SampleVariantSpectrum.sample_id,
        models.SampleVariantSpectrum.ref_pos,
        models.SampleVariantSpectrum.variant_type,
        models.SampleVariantSpectrum.ref_base,
        models.SampleVariantSpectrum.alt_base,
    ).filter(
        models.SampleVariantSpectrum.sample_id.in_(sorted_ids),
        models.SampleVariantSpectrum.reference_id == reference_id,
        models.SampleVariantSpectrum.variant_type == "SNP",
    ).order_by(
        models.SampleVariantSpectrum.sample_id,
        models.SampleVariantSpectrum.ref_pos,
    ).all()
    
    hash_input = json.dumps([
        (v.sample_id, v.ref_pos, v.variant_type, v.ref_base, v.alt_base)
        for v in all_variants
    ], sort_keys=True)
    
    return hashlib.md5(hash_input.encode()).hexdigest()


def _get_variants_for_samples(
    db: Session,
    sample_ids: List[int],
    reference_id: int
) -> Tuple[Dict[int, Dict[int, str]], List[int], Dict[int, str]]:
    sorted_ids = sorted(sample_ids)
    
    all_variants = db.query(models.SampleVariantSpectrum).filter(
        models.SampleVariantSpectrum.sample_id.in_(sorted_ids),
        models.SampleVariantSpectrum.reference_id == reference_id,
        models.SampleVariantSpectrum.variant_type == "SNP",
    ).all()
    
    variants_by_sample: Dict[int, Dict[int, str]] = {sid: {} for sid in sorted_ids}
    all_positions: set = set()
    ref_base_by_pos: Dict[int, str] = {}
    
    for v in all_variants:
        variants_by_sample[v.sample_id][v.ref_pos] = v.alt_base
        all_positions.add(v.ref_pos)
        ref_base_by_pos[v.ref_pos] = v.ref_base
    
    sorted_positions = sorted(all_positions)
    return variants_by_sample, sorted_positions, ref_base_by_pos


def _dict_to_tree_node(data: Dict) -> TreeNode:
    node = TreeNode(data.get("name"))
    node.branch_length = data.get("branch_length")
    node.bootstrap_support = data.get("bootstrap_support")
    node.divergence_time = data.get("divergence_time")
    node.is_leaf = data.get("is_leaf", False)
    for child_data in data.get("children", []):
        child = _dict_to_tree_node(child_data)
        child.parent = node
        node.children.append(child)
    return node


def _tree_node_to_schema(node: TreeNode) -> schemas.PhyloTreeNode:
    return schemas.PhyloTreeNode(
        name=node.name,
        children=[_tree_node_to_schema(c) for c in node.children],
        branch_length=node.branch_length,
        bootstrap_support=node.bootstrap_support,
        divergence_time=node.divergence_time,
        is_leaf=node.is_leaf,
    )


def _get_sample_names(db: Session, sample_ids: List[int]) -> Dict[int, str]:
    samples = db.query(models.Sample).filter(models.Sample.id.in_(sample_ids)).all()
    return {s.id: s.name for s in samples}


def _get_sample_dates(db: Session, sample_ids: List[int]) -> Dict[str, datetime]:
    samples = db.query(models.Sample).filter(models.Sample.id.in_(sample_ids)).all()
    dates = {}
    for s in samples:
        if s.collection_date:
            dates[s.name] = s.collection_date
    return dates


def compute_distance_matrix_endpoint(
    db: Session,
    request: schemas.DistanceMatrixRequest
) -> schemas.DistanceMatrixOut:
    ref = db.query(models.ReferenceSequence).filter(
        models.ReferenceSequence.name == request.reference_name
    ).first()
    if not ref:
        raise ValueError(f"Reference sequence '{request.reference_name}' not found")
    
    sorted_ids = sorted(request.sample_ids)
    for sid in sorted_ids:
        sample = db.query(models.Sample).filter(models.Sample.id == sid).first()
        if not sample:
            raise ValueError(f"Sample {sid} not found")
    
    variants_by_sample, all_positions, ref_base_by_pos = _get_variants_for_samples(
        db, sorted_ids, ref.id
    )
    
    seq_matrix = build_sequence_matrix(
        variants_by_sample, all_positions, sorted_ids, ref_base_by_pos
    )
    
    dist_matrix, warnings, _ = compute_distance_matrix(seq_matrix, sorted_ids)
    
    sample_names_map = _get_sample_names(db, sorted_ids)
    sample_names = [sample_names_map[sid] for sid in sorted_ids]
    
    return schemas.DistanceMatrixOut(
        sample_ids=sorted_ids,
        sample_names=sample_names,
        distance_matrix=dist_matrix,
        warnings=warnings,
    )


def _check_cached_result(
    db: Session,
    cache_key: str,
    data_hash: str
) -> Optional[models.PhyloTreeResult]:
    result = db.query(models.PhyloTreeResult).filter(
        models.PhyloTreeResult.cache_key == cache_key
    ).first()
    
    if not result:
        return None
    
    result.last_checked_at = func.now()
    
    if result.is_stale or result.data_hash != data_hash:
        result.is_stale = True
        db.commit()
        return None
    
    db.commit()
    return result


def _save_result(
    db: Session,
    cache_key: str,
    sample_ids: List[int],
    reference_name: str,
    tree_method: str,
    bootstrap_replicates: int,
    with_molecular_clock: bool,
    root: TreeNode,
    newick_tree: str,
    distance_matrix: List[List[float]],
    sample_names: List[str],
    distance_warnings: List[Dict],
    molecular_clock_info: Optional[Dict],
    data_hash: str
) -> models.PhyloTreeResult:
    tree_json = root.to_dict()
    
    result = models.PhyloTreeResult(
        cache_key=cache_key,
        sample_ids=sorted(sample_ids),
        reference_name=reference_name,
        tree_method=tree_method,
        bootstrap_replicates=bootstrap_replicates,
        with_molecular_clock=with_molecular_clock,
        newick_tree=newick_tree,
        tree_json=tree_json,
        distance_matrix=distance_matrix,
        sample_names=sample_names,
        distance_warnings=distance_warnings,
        molecular_clock_info=molecular_clock_info,
        data_hash=data_hash,
        is_stale=False,
    )
    db.add(result)
    db.commit()
    db.refresh(result)
    return result


def _run_phylogeny_task(task_id: str):
    db = SessionLocal()
    try:
        task = db.query(models.PhyloTreeTask).filter(
            models.PhyloTreeTask.task_id == task_id
        ).first()
        
        if not task:
            return
        
        task.status = "running"
        db.commit()
        
        try:
            ref = db.query(models.ReferenceSequence).filter(
                models.ReferenceSequence.name == task.reference_name
            ).first()
            if not ref:
                raise ValueError(f"Reference sequence '{task.reference_name}' not found")
            
            sorted_ids = sorted(task.sample_ids)
            for sid in sorted_ids:
                sample = db.query(models.Sample).filter(models.Sample.id == sid).first()
                if not sample:
                    raise ValueError(f"Sample {sid} not found")
            
            task.progress_message = "Building sequence matrix..."
            db.commit()
            
            variants_by_sample, all_positions, ref_base_by_pos = _get_variants_for_samples(
                db, sorted_ids, ref.id
            )
            
            seq_matrix = build_sequence_matrix(
                variants_by_sample, all_positions, sorted_ids, ref_base_by_pos
            )
            
            task.progress_message = "Computing distance matrix..."
            db.commit()
            
            dist_matrix, dist_warnings, _ = compute_distance_matrix(seq_matrix, sorted_ids)
            
            sample_names_map = _get_sample_names(db, sorted_ids)
            sample_names = [sample_names_map[sid] for sid in sorted_ids]
            
            task.progress_message = f"Building {task.tree_method.upper()} tree..."
            db.commit()
            
            root, newick = build_tree(dist_matrix, sample_names, task.tree_method)
            
            if task.bootstrap_replicates > 0:
                task.progress_message = "Running bootstrap..."
                task.total_steps = task.bootstrap_replicates
                task.completed_steps = 0
                db.commit()
                
                def bootstrap_progress(completed: int):
                    task_db = SessionLocal()
                    try:
                        t = task_db.query(models.PhyloTreeTask).filter(
                            models.PhyloTreeTask.task_id == task_id
                        ).first()
                        if t:
                            t.completed_steps = completed
                            t.progress_message = f"Bootstrap: {completed}/{task.bootstrap_replicates}"
                            task_db.commit()
                            task_manager.notify_progress_sync(
                                task_id, completed, task.bootstrap_replicates, "running"
                            )
                    finally:
                        task_db.close()
                
                root = compute_bootstrap_support(
                    root,
                    seq_matrix,
                    sorted_ids,
                    sample_names,
                    n_replicates=task.bootstrap_replicates,
                    method=task.tree_method,
                    progress_callback=bootstrap_progress
                )
                
                newick = to_newick_format(root)
            
            clock_info = None
            if task.with_molecular_clock:
                task.progress_message = "Fitting molecular clock..."
                db.commit()
                
                sample_dates = _get_sample_dates(db, sorted_ids)
                root, clock_info, clock_warning = apply_molecular_clock(
                    root, sample_dates, sample_names
                )
                newick = to_newick_format(root)
            
            task.progress_message = "Saving result..."
            db.commit()
            
            cache_key = _get_cache_key(
                task.sample_ids,
                task.reference_name,
                task.tree_method,
                task.bootstrap_replicates,
                task.with_molecular_clock,
            )
            
            data_hash = _compute_data_hash(db, sorted_ids, ref.id)
            
            result = _save_result(
                db,
                cache_key,
                sorted_ids,
                task.reference_name,
                task.tree_method,
                task.bootstrap_replicates,
                task.with_molecular_clock,
                root,
                newick,
                dist_matrix,
                sample_names,
                dist_warnings,
                clock_info,
                data_hash,
            )
            
            task.result_id = result.id
            task.status = "completed"
            task.completed_at = datetime.utcnow()
            task.completed_steps = task.total_steps
            task.progress_message = "Completed"
            db.commit()
            
            task_manager.notify_progress_sync(
                task_id, task.total_steps, task.total_steps, "completed"
            )
            
        except Exception as e:
            task.status = "failed"
            task.error_message = str(e)
            task.completed_at = datetime.utcnow()
            db.commit()
            task_manager.notify_progress_sync(task_id, 0, 1, "failed")
            
    finally:
        db.close()


def create_phylogeny_task(
    db: Session,
    request: schemas.PhyloTreeRequest
) -> schemas.PhyloTreeOut:
    ref = db.query(models.ReferenceSequence).filter(
        models.ReferenceSequence.name == request.reference_name
    ).first()
    if not ref:
        raise ValueError(f"Reference sequence '{request.reference_name}' not found")
    
    sorted_ids = sorted(request.sample_ids)
    for sid in sorted_ids:
        sample = db.query(models.Sample).filter(models.Sample.id == sid).first()
        if not sample:
            raise ValueError(f"Sample {sid} not found")
    
    cache_key = _get_cache_key(
        sorted_ids,
        request.reference_name,
        request.tree_method,
        request.bootstrap_replicates,
        request.with_molecular_clock,
    )
    
    data_hash = _compute_data_hash(db, sorted_ids, ref.id)
    
    cached_result = _check_cached_result(db, cache_key, data_hash)
    
    if cached_result:
        root_node = _dict_to_tree_node(cached_result.tree_json)
        tree_schema = _tree_node_to_schema(root_node)
        
        clock_info = None
        if cached_result.molecular_clock_info:
            clock_info = schemas.MolecularClockInfo(**cached_result.molecular_clock_info)
        
        return schemas.PhyloTreeOut(
            task_id="cached",
            status="completed",
            newick_tree=cached_result.newick_tree,
            tree_json=tree_schema,
            sample_ids=cached_result.sample_ids,
            sample_names=cached_result.sample_names,
            tree_method=cached_result.tree_method,
            bootstrap_replicates=cached_result.bootstrap_replicates,
            distance_matrix=cached_result.distance_matrix,
            distance_warnings=cached_result.distance_warnings or [],
            molecular_clock=clock_info,
            is_stale=cached_result.is_stale,
            created_at=cached_result.created_at,
        )
    
    task_uuid = str(uuid.uuid4())
    
    is_async = len(sorted_ids) > 20 or request.bootstrap_replicates > 0
    
    total_steps = request.bootstrap_replicates if request.bootstrap_replicates > 0 else 100
    
    task = models.PhyloTreeTask(
        task_id=task_uuid,
        status="pending",
        sample_ids=sorted_ids,
        reference_name=request.reference_name,
        tree_method=request.tree_method,
        bootstrap_replicates=request.bootstrap_replicates,
        with_molecular_clock=request.with_molecular_clock,
        total_steps=total_steps,
        completed_steps=0,
        progress_message="Task created",
        is_async=is_async,
    )
    db.add(task)
    db.commit()
    db.refresh(task)
    
    if is_async:
        thread = threading.Thread(target=_run_phylogeny_task, args=(task_uuid,), daemon=True)
        thread.start()
        
        return schemas.PhyloTreeOut(
            task_id=task_uuid,
            status="running",
            sample_ids=sorted_ids,
            sample_names=[],
            tree_method=request.tree_method,
            bootstrap_replicates=request.bootstrap_replicates,
            created_at=task.created_at,
            progress=0,
            total=total_steps,
            progress_message="Task queued",
        )
    else:
        _run_phylogeny_task(task_uuid)
        
        db.refresh(task)
        return get_task_result(db, task_uuid)


def get_task_result(db: Session, task_id: str) -> schemas.PhyloTreeOut:
    if task_id == "cached":
        raise ValueError("Invalid task ID")
    
    task = db.query(models.PhyloTreeTask).filter(
        models.PhyloTreeTask.task_id == task_id
    ).first()
    
    if not task:
        raise ValueError(f"Task {task_id} not found")
    
    if task.status != "completed" or not task.result_id:
        return schemas.PhyloTreeOut(
            task_id=task_id,
            status=task.status,
            sample_ids=task.sample_ids,
            sample_names=[],
            tree_method=task.tree_method,
            bootstrap_replicates=task.bootstrap_replicates,
            created_at=task.created_at,
            progress=task.completed_steps,
            total=task.total_steps,
            progress_message=task.progress_message,
            error_message=task.error_message,
        )
    
    result = db.query(models.PhyloTreeResult).filter(
        models.PhyloTreeResult.id == task.result_id
    ).first()
    
    if not result:
        return schemas.PhyloTreeOut(
            task_id=task_id,
            status="failed",
            sample_ids=task.sample_ids,
            sample_names=[],
            tree_method=task.tree_method,
            bootstrap_replicates=task.bootstrap_replicates,
            created_at=task.created_at,
            error_message="Result not found",
        )
    
    root_node = _dict_to_tree_node(result.tree_json)
    tree_schema = _tree_node_to_schema(root_node)
    
    clock_info = None
    if result.molecular_clock_info:
        clock_info = schemas.MolecularClockInfo(**result.molecular_clock_info)
    
    return schemas.PhyloTreeOut(
        task_id=task_id,
        status=task.status,
        newick_tree=result.newick_tree,
        tree_json=tree_schema,
        sample_ids=result.sample_ids,
        sample_names=result.sample_names,
        tree_method=result.tree_method,
        bootstrap_replicates=result.bootstrap_replicates,
        distance_matrix=result.distance_matrix,
        distance_warnings=result.distance_warnings or [],
        molecular_clock=clock_info,
        is_stale=result.is_stale,
        created_at=result.created_at,
    )


def get_task_status(db: Session, task_id: str) -> schemas.PhyloTaskStatusOut:
    task = db.query(models.PhyloTreeTask).filter(
        models.PhyloTreeTask.task_id == task_id
    ).first()
    
    if not task:
        raise ValueError(f"Task {task_id} not found")
    
    return schemas.PhyloTaskStatusOut(
        task_id=task.task_id,
        status=task.status,
        progress=task.completed_steps,
        total=task.total_steps,
        progress_message=task.progress_message,
        is_async=task.is_async,
        created_at=task.created_at,
        completed_at=task.completed_at,
        error_message=task.error_message,
    )


def compare_trees_endpoint(
    db: Session,
    request: schemas.TreeCompareRequest
) -> schemas.TreeCompareOut:
    task_a = db.query(models.PhyloTreeTask).filter(
        models.PhyloTreeTask.task_id == request.task_id_a
    ).first()
    task_b = db.query(models.PhyloTreeTask).filter(
        models.PhyloTreeTask.task_id == request.task_id_b
    ).first()
    
    if not task_a or not task_a.result_id:
        raise ValueError(f"Task {request.task_id_a} not found or not completed")
    if not task_b or not task_b.result_id:
        raise ValueError(f"Task {request.task_id_b} not found or not completed")
    
    result_a = db.query(models.PhyloTreeResult).filter(
        models.PhyloTreeResult.id == task_a.result_id
    ).first()
    result_b = db.query(models.PhyloTreeResult).filter(
        models.PhyloTreeResult.id == task_b.result_id
    ).first()
    
    if not result_a or not result_b:
        raise ValueError("Result not found")
    
    cache_key = hashlib.md5(
        f"{min(request.task_id_a, request.task_id_b)}:{max(request.task_id_a, request.task_id_b)}".encode()
    ).hexdigest()
    
    cached = db.query(models.TreeComparisonResult).filter(
        models.TreeComparisonResult.cache_key == cache_key
    ).first()
    
    if cached:
        return schemas.TreeCompareOut(
            rf_distance=cached.rf_distance,
            normalized_rf_distance=cached.normalized_rf_distance,
            matching_splits=cached.matching_splits,
            total_splits=cached.total_splits,
            weighted_branch_length_distance=cached.weight_dist,
            inconsistent_branches=[
                schemas.InconsistentBranch(**b)
                for b in (cached.inconsistent_branches or [])
            ],
            sample_names_a=result_a.sample_names,
            sample_names_b=result_b.sample_names,
        )
    
    tree_a = _dict_to_tree_node(result_a.tree_json)
    tree_b = _dict_to_tree_node(result_b.tree_json)
    
    comp_result = compare_trees(
        tree_a, tree_b, result_a.sample_names, result_b.sample_names
    )
    
    comp_record = models.TreeComparisonResult(
        cache_key=cache_key,
        task_id_a=request.task_id_a,
        task_id_b=request.task_id_b,
        rf_distance=comp_result["rf_distance"],
        normalized_rf_distance=comp_result["normalized_rf_distance"],
        matching_splits=comp_result["matching_splits"],
        total_splits=comp_result["total_splits"],
        weight_dist=comp_result["weight_dist"],
        inconsistent_branches=comp_result["inconsistent_branches"],
    )
    db.add(comp_record)
    db.commit()
    
    return schemas.TreeCompareOut(
        rf_distance=comp_result["rf_distance"],
        normalized_rf_distance=comp_result["normalized_rf_distance"],
        matching_splits=comp_result["matching_splits"],
        total_splits=comp_result["total_splits"],
        weighted_branch_length_distance=comp_result["weight_dist"],
        inconsistent_branches=[
            schemas.InconsistentBranch(**b)
            for b in comp_result["inconsistent_branches"]
        ],
        sample_names_a=result_a.sample_names,
        sample_names_b=result_b.sample_names,
    )


def invalidate_tree_results_for_samples(
    db: Session,
    sample_ids: List[int]
):
    affected = db.query(models.PhyloTreeResult).filter(
        models.PhyloTreeResult.is_stale == False
    ).all()
    
    sample_id_set = set(sample_ids)
    for result in affected:
        result_sample_ids = set(result.sample_ids or [])
        if result_sample_ids & sample_id_set:
            result.is_stale = True
    
    db.commit()
