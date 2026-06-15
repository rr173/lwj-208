from fastapi import APIRouter, Depends, HTTPException
from sqlalchemy.orm import Session
from app import schemas
from app.database import get_db
from app.services.genome_service import get_reference_by_name
from app.primer.primer_design import (
    design_primers,
    get_cached_primer_result,
    save_primer_result,
    compute_tm_consistency,
    get_primer_design_history,
)

router = APIRouter(prefix="/api/primer", tags=["primer"])


@router.post("/design", response_model=schemas.PrimerDesignResult)
def design_primer_pair(request: schemas.PrimerDesignRequest, db: Session = Depends(get_db)):
    if request.target_start >= request.target_end:
        raise HTTPException(
            status_code=400,
            detail="target_start must be less than target_end",
        )

    reference = get_reference_by_name(db, request.reference_name)
    if not reference:
        raise HTTPException(
            status_code=404,
            detail=f"Reference sequence '{request.reference_name}' not found",
        )

    if request.target_end >= reference.length:
        raise HTTPException(
            status_code=400,
            detail=f"target_end ({request.target_end}) exceeds reference length ({reference.length})",
        )

    cached = get_cached_primer_result(db, request.reference_name, request.target_start, request.target_end)
    if cached is not None:
        cached["cached"] = True
        return cached

    result = design_primers(
        reference_sequence=reference.sequence,
        ref_name=request.reference_name,
        target_start=request.target_start,
        target_end=request.target_end,
    )

    if "error" in result:
        raise HTTPException(status_code=422, detail=result["error"])

    result["cached"] = False

    save_primer_result(db, request.reference_name, request.target_start, request.target_end, result)

    return result


@router.post("/batch-design", response_model=schemas.BatchPrimerDesignResult)
def batch_design_primers(request: schemas.BatchPrimerDesignRequest, db: Session = Depends(get_db)):
    reference = get_reference_by_name(db, request.reference_name)
    if not reference:
        raise HTTPException(
            status_code=404,
            detail=f"Reference sequence '{request.reference_name}' not found",
        )

    best_pairs = []
    batch_results = []
    success_count = 0
    failure_count = 0

    for idx, target in enumerate(request.targets):
        target_start = target.target_start
        target_end = target.target_end

        item_result = {
            "target_index": idx,
            "target_start": target_start,
            "target_end": target_end,
            "best_primer_pair": None,
            "success": False,
            "error_message": None,
            "cached": False,
        }

        if target_start >= target_end:
            item_result["error_message"] = "target_start must be less than target_end"
            batch_results.append(item_result)
            best_pairs.append(None)
            failure_count += 1
            continue

        if target_end >= reference.length:
            item_result["error_message"] = (
                f"target_end ({target_end}) exceeds reference length ({reference.length})"
            )
            batch_results.append(item_result)
            best_pairs.append(None)
            failure_count += 1
            continue

        try:
            cached = get_cached_primer_result(db, request.reference_name, target_start, target_end)
            if cached is not None:
                item_result["cached"] = True
                pairs = cached.get("primer_pairs", [])
                if pairs:
                    item_result["best_primer_pair"] = pairs[0]
                    best_pairs.append(pairs[0])
                else:
                    best_pairs.append(None)
                item_result["success"] = True
                success_count += 1
                batch_results.append(item_result)
                continue

            result = design_primers(
                reference_sequence=reference.sequence,
                ref_name=request.reference_name,
                target_start=target_start,
                target_end=target_end,
            )

            if "error" in result:
                item_result["error_message"] = result["error"]
                best_pairs.append(None)
                failure_count += 1
                batch_results.append(item_result)
                continue

            save_primer_result(db, request.reference_name, target_start, target_end, result)

            pairs = result.get("primer_pairs", [])
            if pairs:
                item_result["best_primer_pair"] = pairs[0]
                best_pairs.append(pairs[0])
            else:
                best_pairs.append(None)

            item_result["success"] = True
            success_count += 1
            batch_results.append(item_result)

        except Exception as e:
            item_result["error_message"] = str(e)
            best_pairs.append(None)
            failure_count += 1
            batch_results.append(item_result)

    tm_consistency = compute_tm_consistency(best_pairs)

    return {
        "reference_name": request.reference_name,
        "total_targets": len(request.targets),
        "success_count": success_count,
        "failure_count": failure_count,
        "results": batch_results,
        "tm_consistency": tm_consistency,
    }


@router.get("/history/{reference_name}", response_model=schemas.PrimerDesignHistoryList)
def get_primer_history(reference_name: str, db: Session = Depends(get_db)):
    reference = get_reference_by_name(db, reference_name)
    if not reference:
        raise HTTPException(
            status_code=404,
            detail=f"Reference sequence '{reference_name}' not found",
        )

    history = get_primer_design_history(db, reference_name)

    return {
        "reference_name": reference_name,
        "total_records": len(history),
        "history": history,
    }
