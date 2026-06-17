import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from datetime import datetime
from app.database import SessionLocal, Base, engine
from app import models
from app.services import genome_service, sample_service, pipeline_service
from app import schemas
from app.sample_data import (
    build_sample_reference,
    SAMPLE_REFERENCE_NAME,
    SAMPLE_ANNOTATIONS,
    get_sample_queries,
)


def separator(title):
    print(f"\n{'='*60}")
    print(f"  {title}")
    print(f"{'='*60}\n")


def init_database():
    Base.metadata.create_all(bind=engine)
    db = SessionLocal()
    try:
        from sqlalchemy import text, inspect
        inspector = inspect(engine)

        exec_cols = [c["name"] for c in inspector.get_columns("pipeline_executions")]
        if "template_snapshot" not in exec_cols:
            db.execute(text("ALTER TABLE pipeline_executions ADD COLUMN template_snapshot JSON"))
        if "resume_count" not in exec_cols:
            db.execute(text("ALTER TABLE pipeline_executions ADD COLUMN resume_count INTEGER DEFAULT 0"))

        step_cols = [c["name"] for c in inspector.get_columns("pipeline_step_executions")]
        if "retry_count" not in step_cols:
            db.execute(text("ALTER TABLE pipeline_step_executions ADD COLUMN retry_count INTEGER DEFAULT 0"))
        if "last_error" not in step_cols:
            db.execute(text("ALTER TABLE pipeline_step_executions ADD COLUMN last_error TEXT"))
        if "output_data" not in step_cols:
            db.execute(text("ALTER TABLE pipeline_step_executions ADD COLUMN output_data JSON"))

        db.commit()

        existing = genome_service.get_reference_by_name(db, SAMPLE_REFERENCE_NAME)
        if not existing:
            ref_seq = build_sample_reference()
            ref = genome_service.create_reference_sequence(db, SAMPLE_REFERENCE_NAME, ref_seq)
            genome_service.add_gene_annotations(db, ref.id, SAMPLE_ANNOTATIONS)
            queries = get_sample_queries(ref_seq)
            for q in queries:
                genome_service.align_query_all_references(db, q["sequence"])
            db.commit()
            print("Sample reference and alignment data initialized.")
        else:
            print("Sample reference already exists.")

        samples = sample_service.list_samples(db)
        if len(samples) < 10:
            from app.services.batch_service import create_batch_task
            ref = genome_service.get_reference_by_name(db, SAMPLE_REFERENCE_NAME)
            queries = get_sample_queries(ref.sequence)

            for i in range(10):
                sample_data = schemas.SampleCreate(
                    name=f"PipelineTest_Sample_{i+1:03d}",
                    species="Homo sapiens",
                    collection_date=datetime(2025, 1, 1 + i),
                    notes=f"Test sample {i+1} for pipeline testing"
                )
                sample = sample_service.create_sample(db, sample_data)

                batch = create_batch_task(db, [q["sequence"] for q in queries[:5]])
                sample_service.link_sample_to_alignments(
                    db, sample.id, [item.alignment_id for item in batch.items if item.alignment_id]
                )
            db.commit()
            print("Test samples created and linked.")
        else:
            print("Test samples already exist.")
    finally:
        db.close()


def test_step_type_validation():
    separator("TEST 1: Step Type Validation")
    db = SessionLocal()
    try:
        invalid_request = schemas.PipelineTemplateCreate(
            name="test_invalid",
            description="Test with invalid step type",
            steps=[
                schemas.PipelineStepCreate(
                    name="invalid_step",
                    step_type="invalid_type",
                    input_params=[]
                )
            ]
        )
        print("✗ Should have failed for invalid step type")
    except ValueError as e:
        print(f"✓ Correctly rejected invalid step type: {e}")
    finally:
        db.close()


def test_template_crud():
    separator("TEST 2: Pipeline Template CRUD Operations")
    db = SessionLocal()
    try:
        create_req = schemas.PipelineTemplateCreate(
            name="test_pipeline_crud",
            description="Test pipeline for CRUD operations",
            steps=[
                schemas.PipelineStepCreate(
                    name="build_spectrum",
                    step_type="spectrum",
                    input_params=[
                        schemas.PipelineStepInputParam(
                            name="sample_ids",
                            from_step="initial.sample_ids"
                        ),
                        schemas.PipelineStepInputParam(
                            name="reference_name",
                            from_step="initial.reference_name"
                        )
                    ]
                ),
                schemas.PipelineStepCreate(
                    name="extract_signatures",
                    step_type="nmf_extract",
                    input_params=[
                        schemas.PipelineStepInputParam(
                            name="sample_ids",
                            from_step="step_1.output.sample_ids"
                        ),
                        schemas.PipelineStepInputParam(
                            name="reference_name",
                            from_step="step_1.output.reference_name"
                        ),
                        schemas.PipelineStepInputParam(
                            name="k_value",
                            value=3
                        )
                    ],
                    condition=schemas.PipelineStepCondition(
                        expression="step_1.summary.sample_count >= 5"
                    )
                )
            ]
        )

        template = pipeline_service.create_pipeline_template(db, create_req)
        print(f"✓ Created template: ID={template.id}, name={template.name}")
        print(f"  Steps: {len(template.steps)} steps")

        fetched = pipeline_service.get_pipeline_template(db, template.id)
        print(f"✓ Fetched template: {fetched.name}")

        templates = pipeline_service.list_pipeline_templates(db)
        print(f"✓ Listed {len(templates)} templates total")

        update_req = schemas.PipelineTemplateUpdate(
            description="Updated description"
        )
        updated = pipeline_service.update_pipeline_template(db, template, update_req)
        print(f"✓ Updated template description: {updated.description}")

        fetched2 = pipeline_service.get_pipeline_template(db, template.id)
        assert fetched2.description == "Updated description"
        print("✓ Update persisted correctly")

        pipeline_service.delete_pipeline_template(db, template)
        print(f"✓ Deleted template ID={template.id}")

        deleted = pipeline_service.get_pipeline_template(db, template.id)
        assert deleted is None
        print("✓ Template no longer exists after deletion")

    finally:
        db.close()


def test_duplicate_name_validation():
    separator("TEST 3: Duplicate Name Validation")
    db = SessionLocal()
    try:
        create_req = schemas.PipelineTemplateCreate(
            name="test_duplicate",
            description="First template",
            steps=[
                schemas.PipelineStepCreate(
                    name="step1",
                    step_type="spectrum",
                    input_params=[
                        schemas.PipelineStepInputParam(
                            name="sample_ids",
                            from_step="initial.sample_ids"
                        ),
                        schemas.PipelineStepInputParam(
                            name="reference_name",
                            value="test_ref"
                        )
                    ]
                )
            ]
        )

        template = pipeline_service.create_pipeline_template(db, create_req)
        print(f"✓ Created first template: {template.name}")

        try:
            template2 = pipeline_service.create_pipeline_template(db, create_req)
            print("✗ Should have failed for duplicate name")
        except Exception as e:
            db.rollback()
            print(f"✓ Correctly rejected duplicate name: {e}")

        pipeline_service.delete_pipeline_template(db, template)

    finally:
        db.close()


def test_invalid_param_reference():
    separator("TEST 4: Invalid Parameter Reference Validation")
    db = SessionLocal()
    try:
        invalid_req = schemas.PipelineTemplateCreate(
            name="test_invalid_ref",
            description="Test with invalid forward reference",
            steps=[
                schemas.PipelineStepCreate(
                    name="step1",
                    step_type="spectrum",
                    input_params=[
                        schemas.PipelineStepInputParam(
                            name="sample_ids",
                            from_step="step_2.output.sample_ids"
                        ),
                        schemas.PipelineStepInputParam(
                            name="reference_name",
                            value="test_ref"
                        )
                    ]
                ),
                schemas.PipelineStepCreate(
                    name="step2",
                    step_type="nmf_extract",
                    input_params=[
                        schemas.PipelineStepInputParam(
                            name="sample_ids",
                            from_step="initial.sample_ids"
                        ),
                        schemas.PipelineStepInputParam(
                            name="reference_name",
                            value="test_ref"
                        ),
                        schemas.PipelineStepInputParam(
                            name="k_value",
                            value=3
                        )
                    ]
                )
            ]
        )
        print("✗ Should have failed for forward reference")
    except ValueError as e:
        print(f"✓ Correctly rejected forward reference: {e}")
    finally:
        db.close()


def test_invalid_condition_reference():
    separator("TEST 5: Invalid Condition Reference Validation")
    db = SessionLocal()
    try:
        invalid_req = schemas.PipelineTemplateCreate(
            name="test_invalid_condition",
            description="Test with invalid condition reference",
            steps=[
                schemas.PipelineStepCreate(
                    name="step1",
                    step_type="spectrum",
                    input_params=[
                        schemas.PipelineStepInputParam(
                            name="sample_ids",
                            from_step="initial.sample_ids"
                        ),
                        schemas.PipelineStepInputParam(
                            name="reference_name",
                            value="test_ref"
                        )
                    ],
                    condition=schemas.PipelineStepCondition(
                        expression="step_5.output.count > 5"
                    )
                )
            ]
        )
        print("✗ Should have failed for invalid condition reference")
    except ValueError as e:
        print(f"✓ Correctly rejected invalid condition reference: {e}")
    finally:
        db.close()


def test_pipeline_execution():
    separator("TEST 6: Pipeline Execution")
    db = SessionLocal()
    try:
        samples = sample_service.list_samples(db, limit=10)
        sample_ids = [s.id for s in samples[:10]]
        print(f"Using {len(sample_ids)} samples for execution test")

        create_req = schemas.PipelineTemplateCreate(
            name="test_execution_pipeline",
            description="Test pipeline for execution",
            steps=[
                schemas.PipelineStepCreate(
                    name="build_spectrum",
                    step_type="spectrum",
                    input_params=[
                        schemas.PipelineStepInputParam(
                            name="sample_ids",
                            from_step="initial.sample_ids"
                        ),
                        schemas.PipelineStepInputParam(
                            name="reference_name",
                            from_step="initial.reference_name"
                        )
                    ]
                ),
                schemas.PipelineStepCreate(
                    name="ld_analysis",
                    step_type="ld_analysis",
                    input_params=[
                        schemas.PipelineStepInputParam(
                            name="reference_name",
                            from_step="step_1.output.reference_name"
                        )
                    ],
                    condition=schemas.PipelineStepCondition(
                        expression="step_1.summary.sample_count >= 5"
                    )
                )
            ]
        )

        template = pipeline_service.create_pipeline_template(db, create_req)
        print(f"✓ Created execution template: ID={template.id}")

        initial_params = {
            "sample_ids": sample_ids,
            "reference_name": SAMPLE_REFERENCE_NAME
        }

        print("Executing pipeline (this may take a moment)...")
        execution = pipeline_service.execute_pipeline(db, template, initial_params)
        print(f"✓ Pipeline executed: ID={execution.id}, status={execution.status}")
        print(f"  Total duration: {execution.total_duration_seconds:.2f}s")

        for step_exec in execution.step_executions:
            print(f"  Step {step_exec.step_index} ({step_exec.step_name}): "
                  f"{step_exec.status} - {step_exec.duration_seconds:.2f}s")
            if step_exec.output_summary:
                print(f"    Output summary: {step_exec.output_summary}")
            if step_exec.error_message:
                print(f"    Error: {step_exec.error_message}")

        assert execution.status == "completed"
        print("✓ Pipeline completed successfully")

        executions, total = pipeline_service.get_pipeline_executions(db, template.id)
        print(f"✓ Found {total} execution records for template")

        executions_filtered, total_filtered = pipeline_service.get_pipeline_executions(
            db, template.id, status="completed"
        )
        print(f"✓ Found {total_filtered} completed executions")

        executions_failed, total_failed = pipeline_service.get_pipeline_executions(
            db, template.id, status="failed"
        )
        print(f"✓ Found {total_failed} failed executions")

        fetched_exec = pipeline_service.get_pipeline_execution(db, execution.id)
        assert fetched_exec is not None
        assert fetched_exec.status == "completed"
        print("✓ Fetched execution details correctly")

        pipeline_service.delete_pipeline_template(db, template)

    except Exception as e:
        print(f"✗ Execution test failed: {e}")
        import traceback
        traceback.print_exc()
    finally:
        db.close()


def test_pipeline_skip_condition():
    separator("TEST 7: Pipeline Skip Condition")
    db = SessionLocal()
    try:
        samples = sample_service.list_samples(db, limit=3)
        sample_ids = [s.id for s in samples[:3]]
        print(f"Using {len(sample_ids)} samples for skip test")

        create_req = schemas.PipelineTemplateCreate(
            name="test_skip_pipeline",
            description="Test pipeline with skip condition",
            steps=[
                schemas.PipelineStepCreate(
                    name="build_spectrum",
                    step_type="spectrum",
                    input_params=[
                        schemas.PipelineStepInputParam(
                            name="sample_ids",
                            from_step="initial.sample_ids"
                        ),
                        schemas.PipelineStepInputParam(
                            name="reference_name",
                            from_step="initial.reference_name"
                        )
                    ]
                ),
                schemas.PipelineStepCreate(
                    name="ld_analysis",
                    step_type="ld_analysis",
                    input_params=[
                        schemas.PipelineStepInputParam(
                            name="reference_name",
                            from_step="step_1.output.reference_name"
                        )
                    ],
                    condition=schemas.PipelineStepCondition(
                        expression="step_1.summary.sample_count >= 100"
                    )
                )
            ]
        )

        template = pipeline_service.create_pipeline_template(db, create_req)
        print(f"✓ Created skip test template: ID={template.id}")

        initial_params = {
            "sample_ids": sample_ids,
            "reference_name": SAMPLE_REFERENCE_NAME
        }

        execution = pipeline_service.execute_pipeline(db, template, initial_params)
        print(f"✓ Pipeline executed: status={execution.status}")

        for step_exec in execution.step_executions:
            print(f"  Step {step_exec.step_index} ({step_exec.step_name}): "
                  f"{step_exec.status} - {step_exec.duration_seconds:.2f}s")
            if step_exec.output_summary:
                print(f"    Output summary: {step_exec.output_summary}")
            if step_exec.error_message:
                print(f"    Error: {step_exec.error_message}")

        assert len(execution.step_executions) == 2
        assert execution.step_executions[0].status == "completed"
        assert execution.step_executions[1].status == "skipped"
        print("✓ Step 2 correctly skipped due to condition")
        print(f"  Step 1 status: {execution.step_executions[0].status}")
        print(f"  Step 2 status: {execution.step_executions[1].status}")

        pipeline_service.delete_pipeline_template(db, template)

    except Exception as e:
        print(f"✗ Skip test failed: {e}")
        import traceback
        traceback.print_exc()
    finally:
        db.close()


def test_pipeline_failure():
    separator("TEST 8: Pipeline Failure Handling")
    db = SessionLocal()
    try:
        create_req = schemas.PipelineTemplateCreate(
            name="test_failure_pipeline",
            description="Test pipeline failure handling",
            steps=[
                schemas.PipelineStepCreate(
                    name="build_spectrum",
                    step_type="spectrum",
                    input_params=[
                        schemas.PipelineStepInputParam(
                            name="sample_ids",
                            value=[999999, 888888]
                        ),
                        schemas.PipelineStepInputParam(
                            name="reference_name",
                            value="nonexistent_ref"
                        )
                    ]
                ),
                schemas.PipelineStepCreate(
                    name="second_step",
                    step_type="ld_analysis",
                    input_params=[
                        schemas.PipelineStepInputParam(
                            name="reference_name",
                            value="test_ref"
                        )
                    ]
                )
            ]
        )

        template = pipeline_service.create_pipeline_template(db, create_req)
        print(f"✓ Created failure test template: ID={template.id}")

        initial_params = {"sample_ids": [1, 2], "reference_name": "test"}

        execution = pipeline_service.execute_pipeline(db, template, initial_params)
        print(f"✓ Pipeline executed: status={execution.status}")
        print(f"  Error message: {execution.error_message}")

        assert execution.status == "failed"
        assert len(execution.step_executions) == 2
        assert execution.step_executions[0].status == "failed"
        assert execution.step_executions[1].status == "pending"
        print("✓ Step 1 failed, Step 2 was not executed (correct behavior)")
        print(f"  Step 1 status: {execution.step_executions[0].status}")
        print(f"  Step 2 status: {execution.step_executions[1].status}")

        executions_filtered, total_filtered = pipeline_service.get_pipeline_executions(
            db, template.id, status="failed"
        )
        print(f"✓ Found {total_filtered} failed executions when filtering")

        pipeline_service.delete_pipeline_template(db, template)

    except Exception as e:
        print(f"✗ Failure test failed: {e}")
        import traceback
        traceback.print_exc()
    finally:
        db.close()


def test_list_step_types():
    separator("TEST 9: List Available Step Types")
    from app.services.pipeline_service import STEP_SERVICE_MAP

    print("Available step types:")
    for step_type, info in STEP_SERVICE_MAP.items():
        print(f"  - {step_type}:")
        print(f"    Required params: {info['required_params']}")
        print(f"    Optional params: {info['optional_params']}")

    assert len(STEP_SERVICE_MAP) == 5
    assert "spectrum" in STEP_SERVICE_MAP
    assert "nmf_extract" in STEP_SERVICE_MAP
    assert "scoring" in STEP_SERVICE_MAP
    assert "ld_analysis" in STEP_SERVICE_MAP
    assert "temporal_analysis" in STEP_SERVICE_MAP
    print("✓ All expected step types are available")


def test_retry_config_defaults():
    separator("TEST 10: Retry Config Defaults")
    db = SessionLocal()
    try:
        create_req = schemas.PipelineTemplateCreate(
            name="test_retry_defaults",
            description="Test default retry config",
            steps=[
                schemas.PipelineStepCreate(
                    name="step_no_retry",
                    step_type="spectrum",
                    input_params=[
                        schemas.PipelineStepInputParam(
                            name="sample_ids",
                            from_step="initial.sample_ids"
                        ),
                        schemas.PipelineStepInputParam(
                            name="reference_name",
                            from_step="initial.reference_name"
                        )
                    ]
                )
            ]
        )
        template = pipeline_service.create_pipeline_template(db, create_req)
        step = template.steps[0]
        assert step.get("retry_count", 0) == 0
        assert step.get("retry_interval_seconds", 0) == 0
        print("✓ Default retry_count=0 and retry_interval_seconds=0")

        pipeline_service.delete_pipeline_template(db, template)
    finally:
        db.close()


def test_retry_config_custom():
    separator("TEST 11: Custom Retry Config")
    db = SessionLocal()
    try:
        create_req = schemas.PipelineTemplateCreate(
            name="test_retry_custom",
            description="Test custom retry config",
            steps=[
                schemas.PipelineStepCreate(
                    name="step_with_retry",
                    step_type="spectrum",
                    input_params=[
                        schemas.PipelineStepInputParam(
                            name="sample_ids",
                            from_step="initial.sample_ids"
                        ),
                        schemas.PipelineStepInputParam(
                            name="reference_name",
                            from_step="initial.reference_name"
                        )
                    ],
                    retry_count=3,
                    retry_interval_seconds=2
                )
            ]
        )
        template = pipeline_service.create_pipeline_template(db, create_req)
        step = template.steps[0]
        assert step.get("retry_count") == 3
        assert step.get("retry_interval_seconds") == 2
        print("✓ Custom retry_count=3 and retry_interval_seconds=2 persisted correctly")

        pipeline_service.delete_pipeline_template(db, template)
    finally:
        db.close()


def test_retry_on_failure():
    separator("TEST 12: Retry On Failure (monkey-patched)")
    db = SessionLocal()
    try:
        call_count = {"n": 0}
        original_execute_step = pipeline_service._execute_step

        def flaky_execute_step(db, step_type, params):
            call_count["n"] += 1
            if call_count["n"] < 3:
                raise ValueError(f"Simulated transient error (attempt {call_count['n']})")
            return original_execute_step(db, step_type, params)

        samples = sample_service.list_samples(db, limit=5)
        sample_ids = [s.id for s in samples[:5]]

        create_req = schemas.PipelineTemplateCreate(
            name="test_retry_flaky",
            description="Test retry with flaky step",
            steps=[
                schemas.PipelineStepCreate(
                    name="build_spectrum",
                    step_type="spectrum",
                    input_params=[
                        schemas.PipelineStepInputParam(
                            name="sample_ids",
                            value=sample_ids
                        ),
                        schemas.PipelineStepInputParam(
                            name="reference_name",
                            value=SAMPLE_REFERENCE_NAME
                        )
                    ],
                    retry_count=3,
                    retry_interval_seconds=0
                )
            ]
        )

        template = pipeline_service.create_pipeline_template(db, create_req)
        print(f"✓ Created template with retry_count=3")

        import unittest.mock as mock
        with mock.patch.object(pipeline_service, "_execute_step", side_effect=flaky_execute_step):
            execution = pipeline_service.execute_pipeline(
                db, template, {"sample_ids": sample_ids, "reference_name": SAMPLE_REFERENCE_NAME}
            )

        print(f"  Execution status: {execution.status}")
        print(f"  _execute_step called: {call_count['n']} times")

        step_exec = execution.step_executions[0]
        print(f"  Step status: {step_exec.status}, retry_count: {step_exec.retry_count}")
        if step_exec.last_error:
            print(f"  Last error: {step_exec.last_error}")

        assert execution.status == "completed"
        assert call_count["n"] == 3
        assert step_exec.status == "completed"
        assert step_exec.retry_count == 2
        print("✓ Step succeeded after 2 retries (3 total calls)")
        print(f"✓ retry_count recorded correctly: {step_exec.retry_count}")

        pipeline_service.delete_pipeline_template(db, template)
    finally:
        db.close()


def test_retry_exhausted():
    separator("TEST 13: Retry Exhausted (all attempts fail)")
    db = SessionLocal()
    try:
        call_count = {"n": 0}

        def always_fail(db, step_type, params):
            call_count["n"] += 1
            raise ValueError(f"Permanent error (attempt {call_count['n']})")

        create_req = schemas.PipelineTemplateCreate(
            name="test_retry_exhausted",
            description="Test when all retries fail",
            steps=[
                schemas.PipelineStepCreate(
                    name="failing_step",
                    step_type="spectrum",
                    input_params=[
                        schemas.PipelineStepInputParam(
                            name="sample_ids",
                            value=[1, 2, 3]
                        ),
                        schemas.PipelineStepInputParam(
                            name="reference_name",
                            value=SAMPLE_REFERENCE_NAME
                        )
                    ],
                    retry_count=2,
                    retry_interval_seconds=0
                )
            ]
        )

        template = pipeline_service.create_pipeline_template(db, create_req)
        print(f"✓ Created template with retry_count=2 (3 total attempts)")

        import unittest.mock as mock
        with mock.patch.object(pipeline_service, "_execute_step", side_effect=always_fail):
            execution = pipeline_service.execute_pipeline(
                db, template, {"sample_ids": [1, 2, 3], "reference_name": SAMPLE_REFERENCE_NAME}
            )

        print(f"  Execution status: {execution.status}")
        print(f"  Error message: {execution.error_message}")
        print(f"  _execute_step called: {call_count['n']} times")

        step_exec = execution.step_executions[0]
        print(f"  Step status: {step_exec.status}, retry_count: {step_exec.retry_count}")
        print(f"  Step last_error: {step_exec.last_error}")

        assert execution.status == "failed"
        assert call_count["n"] == 3
        assert step_exec.status == "failed"
        assert step_exec.retry_count == 2
        assert "after 2 retries" in execution.error_message
        print("✓ Pipeline correctly failed after exhausting all 3 attempts")
        print("✓ 'after N retries' message included in error")

        pipeline_service.delete_pipeline_template(db, template)
    finally:
        db.close()


def test_resume_completed_rejected():
    separator("TEST 14: Resume Completed Execution Rejected")
    db = SessionLocal()
    try:
        samples = sample_service.list_samples(db, limit=5)
        sample_ids = [s.id for s in samples[:5]]

        create_req = schemas.PipelineTemplateCreate(
            name="test_resume_reject_completed",
            description="Test resume rejection for completed exec",
            steps=[
                schemas.PipelineStepCreate(
                    name="build_spectrum",
                    step_type="spectrum",
                    input_params=[
                        schemas.PipelineStepInputParam(
                            name="sample_ids",
                            value=sample_ids
                        ),
                        schemas.PipelineStepInputParam(
                            name="reference_name",
                            value=SAMPLE_REFERENCE_NAME
                        )
                    ]
                )
            ]
        )

        template = pipeline_service.create_pipeline_template(db, create_req)
        execution = pipeline_service.execute_pipeline(
            db, template, {"sample_ids": sample_ids, "reference_name": SAMPLE_REFERENCE_NAME}
        )
        assert execution.status == "completed"
        print(f"✓ Pipeline completed: ID={execution.id}")

        try:
            pipeline_service.resume_pipeline_execution(db, execution)
            print("✗ Should have rejected resuming completed execution")
            assert False
        except pipeline_service.PipelineResumeError as e:
            print(f"✓ Correctly rejected completed execution: {e}")

        pipeline_service.delete_pipeline_template(db, template)
    finally:
        db.close()


def test_resume_template_changed_rejected():
    separator("TEST 15: Resume Rejected When Template Modified")
    db = SessionLocal()
    try:
        create_req = schemas.PipelineTemplateCreate(
            name="test_resume_template_changed",
            description="Template will be modified after exec",
            steps=[
                schemas.PipelineStepCreate(
                    name="bad_step",
                    step_type="spectrum",
                    input_params=[
                        schemas.PipelineStepInputParam(
                            name="sample_ids",
                            value=[999999]
                        ),
                        schemas.PipelineStepInputParam(
                            name="reference_name",
                            value="nonexistent_ref"
                        )
                    ]
                ),
                schemas.PipelineStepCreate(
                    name="second_step",
                    step_type="ld_analysis",
                    input_params=[
                        schemas.PipelineStepInputParam(
                            name="reference_name",
                            value=SAMPLE_REFERENCE_NAME
                        )
                    ]
                )
            ]
        )

        template = pipeline_service.create_pipeline_template(db, create_req)
        execution = pipeline_service.execute_pipeline(
            db, template, {"sample_ids": [1], "reference_name": "x"}
        )
        assert execution.status == "failed"
        print(f"✓ Pipeline failed as expected: ID={execution.id}")

        update_req = schemas.PipelineTemplateUpdate(
            steps=[
                schemas.PipelineStepCreate(
                    name="renamed_step",
                    step_type="nmf_extract",
                    input_params=[
                        schemas.PipelineStepInputParam(
                            name="sample_ids",
                            value=[1]
                        ),
                        schemas.PipelineStepInputParam(
                            name="reference_name",
                            value=SAMPLE_REFERENCE_NAME
                        ),
                        schemas.PipelineStepInputParam(
                            name="k_value",
                            value=2
                        )
                    ]
                )
            ]
        )
        pipeline_service.update_pipeline_template(db, template, update_req)
        print("✓ Template modified (step count & type changed)")

        try:
            pipeline_service.resume_pipeline_execution(db, execution)
            print("✗ Should have rejected resume due to template change")
            assert False
        except pipeline_service.PipelineResumeError as e:
            print(f"✓ Correctly rejected: {e}")
            assert "Template has been modified" in str(e)

        pipeline_service.delete_pipeline_template(db, template)
    finally:
        db.close()


def test_resume_basic_workflow():
    separator("TEST 16: Resume Basic Workflow")
    db = SessionLocal()
    try:
        samples = sample_service.list_samples(db, limit=5)
        sample_ids = [s.id for s in samples[:5]]

        create_req = schemas.PipelineTemplateCreate(
            name="test_resume_basic",
            description="2-step pipeline: step1 succeeds, step2 fails then resumed",
            steps=[
                schemas.PipelineStepCreate(
                    name="build_spectrum",
                    step_type="spectrum",
                    input_params=[
                        schemas.PipelineStepInputParam(
                            name="sample_ids",
                            from_step="initial.sample_ids"
                        ),
                        schemas.PipelineStepInputParam(
                            name="reference_name",
                            from_step="initial.reference_name"
                        )
                    ]
                ),
                schemas.PipelineStepCreate(
                    name="ld_analysis",
                    step_type="ld_analysis",
                    input_params=[
                        schemas.PipelineStepInputParam(
                            name="reference_name",
                            from_step="step_1.output.reference_name"
                        )
                    ]
                )
            ]
        )

        template = pipeline_service.create_pipeline_template(db, create_req)

        original_execute = pipeline_service._execute_step
        call_count = {"n": 0}

        def flaky_execute(db, step_type, params):
            call_count["n"] += 1
            if step_type == "ld_analysis" and call_count["n"] < 3:
                raise ValueError("Simulated LD error on first run")
            return original_execute(db, step_type, params)

        import unittest.mock as mock
        with mock.patch.object(pipeline_service, "_execute_step", side_effect=flaky_execute):
            execution = pipeline_service.execute_pipeline(
                db, template,
                {"sample_ids": sample_ids, "reference_name": SAMPLE_REFERENCE_NAME}
            )

        assert execution.status == "failed", f"Expected failed, got {execution.status}"
        assert len(execution.step_executions) == 2
        assert execution.step_executions[0].status == "completed"
        assert execution.step_executions[1].status == "failed"
        print(f"✓ First run: step1 completed, step2 failed (as simulated)")
        print(f"  Step1 output summary: {execution.step_executions[0].output_summary}")
        orig_step1_duration = execution.step_executions[0].duration_seconds
        orig_step1_output = execution.step_executions[0].output_summary
        print(f"  Step1 duration: {orig_step1_duration:.2f}s")
        print(f"  _execute_step total calls in first run: {call_count['n']}")

        resumed_exec, skipped, resumed_from = pipeline_service.resume_pipeline_execution(
            db, execution
        )

        print(f"  Resumed from step: {resumed_from}")
        print(f"  Skipped completed steps: {skipped}")
        print(f"  Resume count: {resumed_exec.resume_count}")
        print(f"  Final status: {resumed_exec.status}")
        print(f"  Step1 status: {resumed_exec.step_executions[0].status}")
        print(f"  Step2 status: {resumed_exec.step_executions[1].status}")
        print(f"  Total _execute_step calls overall: {call_count['n']}")

        assert resumed_exec.status == "completed"
        assert resumed_from == 2
        assert skipped == 1
        assert resumed_exec.resume_count == 1
        assert resumed_exec.step_executions[0].status == "completed"
        assert resumed_exec.step_executions[0].output_summary == orig_step1_output
        assert resumed_exec.step_executions[0].duration_seconds == orig_step1_duration
        assert resumed_exec.step_executions[1].status == "completed"
        print("✓ Resume succeeded: step1 reused, step2 executed successfully")
        print("✓ Step1 output & duration unchanged (not re-executed)")

        pipeline_service.delete_pipeline_template(db, template)
    finally:
        db.close()


def test_resume_idempotent_all_completed():
    separator("TEST 17: Resume Idempotent - All Steps Already Done")
    db = SessionLocal()
    try:
        samples = sample_service.list_samples(db, limit=5)
        sample_ids = [s.id for s in samples[:5]]

        create_req = schemas.PipelineTemplateCreate(
            name="test_resume_idempotent",
            description="Single step that succeeds",
            steps=[
                schemas.PipelineStepCreate(
                    name="build_spectrum",
                    step_type="spectrum",
                    input_params=[
                        schemas.PipelineStepInputParam(
                            name="sample_ids",
                            value=sample_ids
                        ),
                        schemas.PipelineStepInputParam(
                            name="reference_name",
                            value=SAMPLE_REFERENCE_NAME
                        )
                    ]
                )
            ]
        )

        template = pipeline_service.create_pipeline_template(db, create_req)
        execution = pipeline_service.execute_pipeline(
            db, template, {"sample_ids": sample_ids, "reference_name": SAMPLE_REFERENCE_NAME}
        )
        assert execution.status == "completed"

        db.query(models.PipelineExecution).filter(
            models.PipelineExecution.id == execution.id
        ).update({"status": "failed"})
        db.commit()
        db.refresh(execution)
        print("✓ Manually set execution status to 'failed' for testing")

        resumed_exec, skipped, resumed_from = pipeline_service.resume_pipeline_execution(
            db, execution
        )

        assert resumed_exec.status == "completed"
        assert skipped == 1
        print(f"✓ Resume corrected status to 'completed' without re-running")
        print(f"  Skipped steps: {skipped}, resumed_from: {resumed_from}")

        pipeline_service.delete_pipeline_template(db, template)
    finally:
        db.close()


def cleanup_test_templates():
    separator("CLEANUP")
    db = SessionLocal()
    try:
        templates = db.query(models.PipelineTemplate).filter(
            models.PipelineTemplate.name.like("test_%")
        ).all()
        for t in templates:
            db.delete(t)
        db.commit()
        print(f"✓ Cleaned up {len(templates)} test templates")
    finally:
        db.close()


if __name__ == "__main__":
    print("\n" + "="*60)
    print("  PIPELINE MODULE TEST SUITE")
    print("="*60)

    init_database()

    try:
        test_step_type_validation()
        test_template_crud()
        test_duplicate_name_validation()
        test_invalid_param_reference()
        test_invalid_condition_reference()
        test_list_step_types()
        test_pipeline_execution()
        test_pipeline_skip_condition()
        test_pipeline_failure()
        test_retry_config_defaults()
        test_retry_config_custom()
        test_retry_on_failure()
        test_retry_exhausted()
        test_resume_completed_rejected()
        test_resume_template_changed_rejected()
        test_resume_basic_workflow()
        test_resume_idempotent_all_completed()

        print("\n" + "="*60)
        print("  ALL TESTS COMPLETED")
        print("="*60 + "\n")

    finally:
        cleanup_test_templates()
