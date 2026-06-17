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

        print("\n" + "="*60)
        print("  ALL TESTS COMPLETED")
        print("="*60 + "\n")

    finally:
        cleanup_test_templates()
