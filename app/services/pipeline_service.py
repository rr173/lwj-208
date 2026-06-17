import re
import time
import signal
import copy
import hashlib
import json
from typing import List, Optional, Dict, Any, Tuple
from datetime import datetime
from sqlalchemy.orm import Session

from app import models, schemas
from app.services import mutational_signature_service, scoring_service, ld_service

PIPELINE_TIMEOUT_SECONDS = 300

STEP_SERVICE_MAP = {
    "spectrum": {
        "service": mutational_signature_service.build_mutational_spectrum,
        "required_params": ["sample_ids", "reference_name"],
        "optional_params": [],
    },
    "nmf_extract": {
        "service": mutational_signature_service.extract_signatures,
        "required_params": ["sample_ids", "reference_name", "k_value"],
        "optional_params": [],
    },
    "scoring": {
        "service": scoring_service.score_sample,
        "required_params": ["sample_id"],
        "optional_params": [],
    },
    "ld_analysis": {
        "service": ld_service.compute_ld_matrix_endpoint,
        "required_params": ["reference_name"],
        "optional_params": ["start", "end", "min_maf"],
    },
    "temporal_analysis": {
        "service": mutational_signature_service.analyze_signature_temporal_trends,
        "required_params": ["sample_ids", "reference_name"],
        "optional_params": ["k_value", "p_value_threshold"],
    },
}


class PipelineTimeoutError(Exception):
    pass


class PipelineResumeError(Exception):
    pass


def _timeout_handler(signum, frame):
    raise PipelineTimeoutError(f"Pipeline execution timed out after {PIPELINE_TIMEOUT_SECONDS} seconds")


def _validate_step_type(step_type: str) -> None:
    if step_type not in STEP_SERVICE_MAP:
        raise ValueError(
            f"Invalid step type: {step_type}. "
            f"Valid types are: {list(STEP_SERVICE_MAP.keys())}"
        )


def _get_nested_value(data: Dict[str, Any], path: str) -> Any:
    keys = path.split(".")
    current = data
    for key in keys:
        if isinstance(current, dict) and key in current:
            current = current[key]
        else:
            raise ValueError(f"Cannot resolve path '{path}' in data")
    return current


def _resolve_from_step(from_step_expr: str, step_outputs: Dict[str, Any]) -> Any:
    match = re.match(r"^(step_\d+)\.(.+)$", from_step_expr)
    if not match:
        raise ValueError(f"Invalid from_step format: {from_step_expr}")
    step_ref, field_path = match.groups()
    if step_ref not in step_outputs:
        raise ValueError(f"Step {step_ref} output not available")
    return _get_nested_value(step_outputs[step_ref], field_path)


def _resolve_step_params(
    step_config: Dict[str, Any],
    step_outputs: Dict[str, Any],
    initial_params: Dict[str, Any]
) -> Dict[str, Any]:
    resolved_params = {}
    for param in step_config.get("input_params", []):
        param_name = param["name"]
        if param.get("value") is not None:
            resolved_params[param_name] = param["value"]
        elif param.get("from_step"):
            from_step = param["from_step"]
            if from_step.startswith("initial."):
                field_path = from_step[len("initial."):]
                resolved_params[param_name] = _get_nested_value(initial_params, field_path)
            else:
                resolved_params[param_name] = _resolve_from_step(from_step, step_outputs)
    return resolved_params


def _evaluate_condition(condition_expr: str, step_outputs: Dict[str, Any]) -> bool:
    allowed_pattern = r'^(step_\d+)\.([a-zA-Z0-9_.]+)\s*([><=!]+)\s*(-?\d+\.?\d*|true|false|"[^"]*")$'
    match = re.match(allowed_pattern, condition_expr)
    if not match:
        raise ValueError(f"Invalid condition expression: {condition_expr}")

    step_ref, field_path, operator, value_str = match.groups()

    if step_ref not in step_outputs:
        raise ValueError(f"Step {step_ref} output not available for condition evaluation")

    try:
        actual_value = _get_nested_value(step_outputs[step_ref], field_path)
    except ValueError:
        raise ValueError(f"Cannot resolve field '{field_path}' in {step_ref} output")

    if value_str.lower() == "true":
        compare_value = True
    elif value_str.lower() == "false":
        compare_value = False
    elif value_str.startswith('"') and value_str.endswith('"'):
        compare_value = value_str[1:-1]
    elif "." in value_str:
        compare_value = float(value_str)
    else:
        compare_value = int(value_str)

    operators = {
        ">": lambda a, b: a > b,
        ">=": lambda a, b: a >= b,
        "<": lambda a, b: a < b,
        "<=": lambda a, b: a <= b,
        "==": lambda a, b: a == b,
        "!=": lambda a, b: a != b,
    }

    if operator not in operators:
        raise ValueError(f"Invalid operator: {operator}")

    return operators[operator](actual_value, compare_value)


def _summarize_output(step_type: str, output: Any) -> Dict[str, Any]:
    summary = {}
    if step_type == "spectrum":
        if hasattr(output, "sample_ids"):
            summary["sample_count"] = len(output.sample_ids)
            summary["sample_ids"] = output.sample_ids
            summary["reference_name"] = output.reference_name
            summary["mutation_type_count"] = len(output.mutation_types) if hasattr(output, "mutation_types") else 0
        elif isinstance(output, dict):
            summary["sample_count"] = len(output.get("sample_ids", []))
            summary["sample_ids"] = output.get("sample_ids", [])
            summary["reference_name"] = output.get("reference_name")
            summary["mutation_type_count"] = len(output.get("mutation_types", []))
    elif step_type == "nmf_extract":
        if hasattr(output, "sample_ids"):
            summary["sample_count"] = len(output.sample_ids)
            summary["k_value"] = output.k_value if hasattr(output, "k_value") else None
            summary["signature_count"] = len(output.signatures) if hasattr(output, "signatures") else 0
            summary["reconstruction_error"] = output.reconstruction_error if hasattr(output, "reconstruction_error") else None
        elif isinstance(output, dict):
            summary["sample_count"] = len(output.get("sample_ids", []))
            summary["k_value"] = output.get("k_value")
            summary["signature_count"] = len(output.get("signatures", []))
            summary["reconstruction_error"] = output.get("reconstruction_error")
    elif step_type == "scoring":
        if hasattr(output, "total_score"):
            summary["total_score"] = output.total_score
            summary["classification"] = output.classification if hasattr(output, "classification") else None
            summary["total_variants"] = output.total_variants if hasattr(output, "total_variants") else None
        elif isinstance(output, dict):
            summary["total_score"] = output.get("total_score")
            summary["classification"] = output.get("classification")
            summary["total_variants"] = output.get("total_variants")
    elif step_type == "ld_analysis":
        if hasattr(output, "total_samples"):
            summary["total_samples"] = output.total_samples
            summary["snp_count"] = output.snp_count if hasattr(output, "snp_count") else 0
            summary["reference_name"] = output.reference_name if hasattr(output, "reference_name") else None
        elif isinstance(output, dict):
            summary["total_samples"] = output.get("total_samples")
            summary["snp_count"] = output.get("snp_count", 0)
            summary["reference_name"] = output.get("reference_name")
    elif step_type == "temporal_analysis":
        if hasattr(output, "sample_ids"):
            summary["sample_count"] = len(output.sample_ids)
            summary["total_signatures"] = output.total_signatures if hasattr(output, "total_signatures") else 0
            summary["significant_rising"] = output.significant_rising if hasattr(output, "significant_rising") else []
            summary["significant_falling"] = output.significant_falling if hasattr(output, "significant_falling") else []
        elif isinstance(output, dict):
            summary["sample_count"] = len(output.get("sample_ids", []))
            summary["total_signatures"] = output.get("total_signatures", 0)
            summary["significant_rising"] = output.get("significant_rising", [])
            summary["significant_falling"] = output.get("significant_falling", [])

    return summary


def _make_json_safe(obj: Any) -> Any:
    if isinstance(obj, (str, int, float, bool)) or obj is None:
        return obj
    if isinstance(obj, (list, tuple)):
        return [_make_json_safe(item) for item in obj]
    if isinstance(obj, dict):
        return {str(k): _make_json_safe(v) for k, v in obj.items()}
    if hasattr(obj, "model_dump"):
        try:
            return _make_json_safe(obj.model_dump())
        except Exception:
            return str(obj)
    if hasattr(obj, "__dict__"):
        try:
            return _make_json_safe(vars(obj))
        except Exception:
            return str(obj)
    try:
        json.dumps(obj)
        return obj
    except (TypeError, ValueError):
        return str(obj)


def _output_to_dict(output: Any) -> Dict[str, Any]:
    if isinstance(output, dict):
        return _make_json_safe(output)
    result = {}
    for attr in dir(output):
        if not attr.startswith("_") and not callable(getattr(output, attr)):
            value = getattr(output, attr)
            result[attr] = _make_json_safe(value)
    return result


def _execute_step(
    db: Session,
    step_type: str,
    params: Dict[str, Any]
) -> Any:
    if step_type not in STEP_SERVICE_MAP:
        raise ValueError(f"Unknown step type: {step_type}")

    service_info = STEP_SERVICE_MAP[step_type]
    service_func = service_info["service"]
    required_params = service_info["required_params"]

    for param in required_params:
        if param not in params:
            raise ValueError(f"Missing required parameter '{param}' for step type '{step_type}'")

    return service_func(db, **params)


def _get_template_signature(steps: List[Dict[str, Any]]) -> str:
    signature_data = [
        {
            "step_index": i + 1,
            "step_type": s.get("step_type"),
            "step_name": s.get("name"),
        }
        for i, s in enumerate(steps)
    ]
    raw = json.dumps(signature_data, sort_keys=True, ensure_ascii=False)
    return hashlib.sha256(raw.encode("utf-8")).hexdigest()


def _verify_template_consistency(
    snapshot_steps: List[Dict[str, Any]],
    current_steps: List[Dict[str, Any]]
) -> None:
    if len(snapshot_steps) != len(current_steps):
        raise PipelineResumeError(
            f"Template has been modified: step count changed "
            f"({len(snapshot_steps)} -> {len(current_steps)}). Resume aborted."
        )
    for i, (snap_step, curr_step) in enumerate(zip(snapshot_steps, current_steps)):
        idx = i + 1
        if snap_step.get("step_type") != curr_step.get("step_type"):
            raise PipelineResumeError(
                f"Template has been modified: step {idx} type changed "
                f"({snap_step.get('step_type')} -> {curr_step.get('step_type')}). Resume aborted."
            )
        if snap_step.get("name") != curr_step.get("name"):
            raise PipelineResumeError(
                f"Template has been modified: step {idx} name changed "
                f"('{snap_step.get('name')}' -> '{curr_step.get('name')}'). Resume aborted."
            )


def _run_single_step_with_retry(
    db: Session,
    step_config: Dict[str, Any],
    step_execution: models.PipelineStepExecution,
    step_outputs: Dict[str, Any],
    initial_params: Dict[str, Any],
    execution: models.PipelineExecution,
    pipeline_start_time: float,
) -> Tuple[bool, Optional[Dict[str, Any]], Optional[float]]:
    step_index = step_execution.step_index
    step_key = f"step_{step_index}"
    step_name = step_config.get("name", f"Step {step_index}")
    step_type = step_config["step_type"]
    max_retries = step_config.get("retry_count", 0) or 0
    retry_interval = step_config.get("retry_interval_seconds", 0) or 0

    attempt = 0
    last_error = None
    step_duration = None

    while attempt <= max_retries:
        if attempt > 0 and retry_interval > 0:
            time.sleep(retry_interval)

        step_execution.status = "running"
        step_execution.started_at = datetime.utcnow()
        step_execution.retry_count = attempt
        db.commit()

        condition = step_config.get("condition")
        if condition and attempt == 0:
            try:
                should_execute = _evaluate_condition(condition["expression"], step_outputs)
                if not should_execute:
                    step_execution.status = "skipped"
                    step_execution.completed_at = datetime.utcnow()
                    step_execution.duration_seconds = 0
                    db.commit()
                    return True, None, 0.0
            except Exception as e:
                last_error = f"Condition evaluation failed: {str(e)}"
                step_execution.last_error = last_error
                step_execution.status = "failed"
                step_execution.error_message = last_error
                step_execution.completed_at = datetime.utcnow()
                step_execution.duration_seconds = (
                    step_execution.completed_at - step_execution.started_at
                ).total_seconds() if step_execution.started_at else 0
                execution.status = "failed"
                execution.error_message = f"Step {step_index} ({step_name}) condition failed: {str(e)}"
                execution.completed_at = datetime.utcnow()
                execution.total_duration_seconds = time.time() - pipeline_start_time
                db.commit()
                return False, None, step_execution.duration_seconds

        try:
            resolved_params = _resolve_step_params(step_config, step_outputs, initial_params)
            step_execution.input_params = resolved_params
            db.commit()

            step_start_time = time.time()
            output = _execute_step(db, step_type, resolved_params)
            step_duration = time.time() - step_start_time

            output_dict = _output_to_dict(output)
            output_summary = _summarize_output(step_type, output)
            step_result = {
                "output": output_dict,
                "summary": output_summary
            }

            step_execution.status = "completed"
            step_execution.output_summary = output_summary
            step_execution.output_data = output_dict
            step_execution.completed_at = datetime.utcnow()
            step_execution.duration_seconds = step_duration
            step_execution.last_error = None
            db.commit()

            return True, step_result, step_duration

        except PipelineTimeoutError as e:
            last_error = str(e)
            step_execution.last_error = last_error
            step_execution.status = "failed"
            step_execution.error_message = last_error
            step_execution.completed_at = datetime.utcnow()
            step_execution.duration_seconds = (
                step_execution.completed_at - step_execution.started_at
            ).total_seconds() if step_execution.started_at else 0
            execution.status = "failed"
            execution.error_message = last_error
            execution.completed_at = datetime.utcnow()
            execution.total_duration_seconds = time.time() - pipeline_start_time
            db.commit()
            return False, None, step_execution.duration_seconds

        except Exception as e:
            last_error = str(e)
            step_execution.last_error = last_error
            attempt += 1
            db.commit()

    step_execution.status = "failed"
    step_execution.error_message = last_error
    step_execution.completed_at = datetime.utcnow()
    step_execution.duration_seconds = (
        step_execution.completed_at - step_execution.started_at
    ).total_seconds() if step_execution.started_at else 0
    execution.status = "failed"
    retry_info = f" (after {max_retries} retries)" if max_retries > 0 else ""
    execution.error_message = f"Step {step_index} ({step_name}) failed{retry_info}: {last_error}"
    execution.completed_at = datetime.utcnow()
    execution.total_duration_seconds = time.time() - pipeline_start_time
    db.commit()
    return False, None, step_execution.duration_seconds


def create_pipeline_template(
    db: Session,
    request: schemas.PipelineTemplateCreate
) -> models.PipelineTemplate:
    for step in request.steps:
        _validate_step_type(step.step_type)

    existing = get_pipeline_template_by_name(db, request.name)
    if existing:
        raise ValueError(f"Pipeline template name '{request.name}' already exists")

    steps_data = [step.model_dump() for step in request.steps]

    template = models.PipelineTemplate(
        name=request.name,
        description=request.description or "",
        steps=steps_data,
    )
    db.add(template)
    db.commit()
    db.refresh(template)
    return template


def get_pipeline_template(
    db: Session,
    template_id: int
) -> Optional[models.PipelineTemplate]:
    return db.query(models.PipelineTemplate).filter(
        models.PipelineTemplate.id == template_id
    ).first()


def get_pipeline_template_by_name(
    db: Session,
    name: str
) -> Optional[models.PipelineTemplate]:
    return db.query(models.PipelineTemplate).filter(
        models.PipelineTemplate.name == name
    ).first()


def list_pipeline_templates(
    db: Session,
    skip: int = 0,
    limit: int = 100
) -> List[models.PipelineTemplate]:
    return db.query(models.PipelineTemplate).offset(skip).limit(limit).all()


def update_pipeline_template(
    db: Session,
    template: models.PipelineTemplate,
    request: schemas.PipelineTemplateUpdate
) -> models.PipelineTemplate:
    if request.name is not None:
        existing = get_pipeline_template_by_name(db, request.name)
        if existing and existing.id != template.id:
            raise ValueError(f"Pipeline template name '{request.name}' already exists")
        template.name = request.name

    if request.description is not None:
        template.description = request.description

    if request.steps is not None:
        for step in request.steps:
            _validate_step_type(step.step_type)
        template.steps = [step.model_dump() for step in request.steps]

    db.commit()
    db.refresh(template)
    return template


def delete_pipeline_template(
    db: Session,
    template: models.PipelineTemplate
) -> None:
    db.delete(template)
    db.commit()


def execute_pipeline(
    db: Session,
    template: models.PipelineTemplate,
    initial_params: Dict[str, Any]
) -> models.PipelineExecution:
    steps = template.steps
    step_outputs: Dict[str, Any] = {}
    step_executions = []

    template_snapshot = copy.deepcopy(steps)

    execution = models.PipelineExecution(
        template_id=template.id,
        initial_params=initial_params,
        status="running",
        template_snapshot=template_snapshot,
        resume_count=0,
    )
    db.add(execution)
    db.flush()

    for i, step_config in enumerate(steps):
        step_index = i + 1
        step_name = step_config.get("name", f"Step {step_index}")
        step_type = step_config["step_type"]

        step_execution = models.PipelineStepExecution(
            execution_id=execution.id,
            step_index=step_index,
            step_type=step_type,
            step_name=step_name,
            input_params={},
            status="pending",
            retry_count=0,
            last_error=None,
        )
        db.add(step_execution)
        db.flush()
        step_executions.append(step_execution)

    db.commit()

    signal.signal(signal.SIGALRM, _timeout_handler)
    signal.alarm(PIPELINE_TIMEOUT_SECONDS)

    pipeline_start_time = time.time()

    try:
        for i, step_config in enumerate(steps):
            step_index = i + 1
            step_key = f"step_{step_index}"
            step_execution = step_executions[i]

            if execution.status != "running":
                break

            success, step_result, _ = _run_single_step_with_retry(
                db=db,
                step_config=step_config,
                step_execution=step_execution,
                step_outputs=step_outputs,
                initial_params=initial_params,
                execution=execution,
                pipeline_start_time=pipeline_start_time,
            )

            if not success:
                break

            if step_result is not None:
                step_outputs[step_key] = step_result

        if execution.status == "running":
            execution.status = "completed"
            execution.completed_at = datetime.utcnow()
            execution.total_duration_seconds = time.time() - pipeline_start_time
            db.commit()

    except PipelineTimeoutError as e:
        execution.status = "failed"
        execution.error_message = f"Pipeline timed out after {PIPELINE_TIMEOUT_SECONDS}s: {str(e)}"
        execution.completed_at = datetime.utcnow()
        execution.total_duration_seconds = time.time() - pipeline_start_time
        db.commit()

    finally:
        signal.alarm(0)

    db.refresh(execution)
    return execution


def resume_pipeline_execution(
    db: Session,
    execution: models.PipelineExecution,
) -> Tuple[models.PipelineExecution, int, Optional[int]]:
    template = execution.template
    if template is None:
        raise PipelineResumeError("Associated pipeline template no longer exists")

    if execution.template_snapshot is None:
        raise PipelineResumeError(
            "Original template snapshot is missing. Cannot verify template consistency. "
            "Please re-run the pipeline from scratch."
        )

    _verify_template_consistency(execution.template_snapshot, template.steps)

    if execution.status not in ("failed", "pending"):
        raise PipelineResumeError(
            f"Cannot resume execution with status '{execution.status}'. "
            f"Only 'failed' or 'pending' executions can be resumed."
        )

    snapshot_steps = execution.template_snapshot
    step_executions = list(execution.step_executions)

    if len(step_executions) != len(snapshot_steps):
        raise PipelineResumeError(
            f"Step execution count mismatch with template snapshot "
            f"({len(step_executions)} vs {len(snapshot_steps)})."
        )

    step_outputs: Dict[str, Any] = {}
    start_step_index = None
    skipped_count = 0

    for i, (step_config, step_exec) in enumerate(zip(snapshot_steps, step_executions)):
        idx = i + 1
        step_key = f"step_{idx}"
        if step_exec.status == "completed":
            output_to_use = step_exec.output_data if step_exec.output_data is not None else step_exec.output_summary
            if output_to_use is not None:
                step_outputs[step_key] = {
                    "output": output_to_use,
                    "summary": step_exec.output_summary if step_exec.output_summary is not None else output_to_use,
                }
            skipped_count += 1
            continue
        if step_exec.status == "skipped":
            skipped_count += 1
            continue
        if step_exec.status in ("failed", "pending"):
            if start_step_index is None:
                start_step_index = idx
            break

    if start_step_index is None:
        if execution.status == "completed":
            return execution, skipped_count, None
        all_completed = all(se.status in ("completed", "skipped") for se in step_executions)
        if all_completed:
            execution.status = "completed"
            execution.error_message = None
            db.commit()
            db.refresh(execution)
            return execution, skipped_count, None
        start_step_index = 1
        for i, se in enumerate(step_executions):
            if se.status not in ("completed", "skipped"):
                start_step_index = i + 1
                break

    execution.status = "running"
    execution.error_message = None
    execution.resume_count = (execution.resume_count or 0) + 1
    db.commit()

    signal.signal(signal.SIGALRM, _timeout_handler)
    signal.alarm(PIPELINE_TIMEOUT_SECONDS)

    pipeline_start_time = time.time()

    try:
        for i in range(start_step_index - 1, len(snapshot_steps)):
            step_config = snapshot_steps[i]
            step_execution = step_executions[i]
            step_index = i + 1
            step_key = f"step_{step_index}"

            if execution.status != "running":
                break

            if step_execution.status == "completed":
                if step_key not in step_outputs:
                    output_to_use = step_execution.output_data if step_execution.output_data is not None else step_execution.output_summary
                    if output_to_use is not None:
                        step_outputs[step_key] = {
                            "output": output_to_use,
                            "summary": step_execution.output_summary if step_execution.output_summary is not None else output_to_use,
                        }
                continue
            if step_execution.status == "skipped":
                continue

            success, step_result, _ = _run_single_step_with_retry(
                db=db,
                step_config=step_config,
                step_execution=step_execution,
                step_outputs=step_outputs,
                initial_params=execution.initial_params,
                execution=execution,
                pipeline_start_time=pipeline_start_time,
            )

            if not success:
                break

            if step_result is not None:
                step_outputs[step_key] = step_result

        if execution.status == "running":
            execution.status = "completed"
            execution.completed_at = datetime.utcnow()
            if execution.total_duration_seconds is None:
                execution.total_duration_seconds = time.time() - pipeline_start_time
            else:
                execution.total_duration_seconds += time.time() - pipeline_start_time
            db.commit()
        else:
            if execution.total_duration_seconds is None:
                execution.total_duration_seconds = time.time() - pipeline_start_time
            else:
                execution.total_duration_seconds += time.time() - pipeline_start_time
            db.commit()

    except PipelineTimeoutError as e:
        execution.status = "failed"
        execution.error_message = f"Pipeline timed out after {PIPELINE_TIMEOUT_SECONDS}s during resume: {str(e)}"
        execution.completed_at = datetime.utcnow()
        if execution.total_duration_seconds is None:
            execution.total_duration_seconds = time.time() - pipeline_start_time
        else:
            execution.total_duration_seconds += time.time() - pipeline_start_time
        db.commit()

    finally:
        signal.alarm(0)

    db.refresh(execution)
    return execution, skipped_count, start_step_index


def get_pipeline_executions(
    db: Session,
    template_id: int,
    status: Optional[str] = None,
    skip: int = 0,
    limit: int = 100
) -> Tuple[List[models.PipelineExecution], int]:
    query = db.query(models.PipelineExecution).filter(
        models.PipelineExecution.template_id == template_id
    )

    if status:
        if status not in models.PIPELINE_STATUSES:
            raise ValueError(
                f"Invalid status: {status}. "
                f"Valid statuses are: {models.PIPELINE_STATUSES}"
            )
        query = query.filter(models.PipelineExecution.status == status)

    total = query.count()
    executions = query.order_by(
        models.PipelineExecution.triggered_at.desc()
    ).offset(skip).limit(limit).all()

    return executions, total


def get_pipeline_execution(
    db: Session,
    execution_id: int
) -> Optional[models.PipelineExecution]:
    return db.query(models.PipelineExecution).filter(
        models.PipelineExecution.id == execution_id
    ).first()
