import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from datetime import datetime
from app.database import SessionLocal, Base, engine
from app import models
from app.services import genome_service, sample_service
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
    finally:
        db.close()


def test_sample_crud():
    separator("TEST 1: Sample CRUD Operations")
    db = SessionLocal()
    try:
        sample1_data = schemas.SampleCreate(
            name="Sample_Patient_001",
            species="Homo sapiens",
            collection_date=datetime(2025, 1, 15),
            notes="Blood sample from patient 001"
        )
        sample1 = sample_service.create_sample(db, sample1_data)
        print(f"✓ Created sample: ID={sample1.id}, name={sample1.name}, species={sample1.species}")

        sample2_data = schemas.SampleCreate(
            name="Sample_Patient_002",
            species="Homo sapiens",
            collection_date=datetime(2025, 2, 20),
            notes="Tissue sample from patient 002"
        )
        sample2 = sample_service.create_sample(db, sample2_data)
        print(f"✓ Created sample: ID={sample2.id}, name={sample2.name}, species={sample2.species}")

        sample3_data = schemas.SampleCreate(
            name="Sample_Control_001",
            species="Homo sapiens",
            notes="Healthy control sample"
        )
        sample3 = sample_service.create_sample(db, sample3_data)
        print(f"✓ Created sample: ID={sample3.id}, name={sample3.name}, species={sample3.species}")

        samples = sample_service.list_samples(db)
        print(f"✓ Listed {len(samples)} samples total")

        fetched = sample_service.get_sample_by_id(db, sample1.id)
        out = sample_service.sample_to_out(db, fetched, with_detail=True)
        print(f"✓ Sample detail: name={out.name}, alignments={out.alignment_count}, variants={out.variant_count}")

        update = schemas.SampleUpdate(notes="Updated: blood sample, confirmed diagnosis")
        updated = sample_service.update_sample(db, sample1, update)
        print(f"✓ Updated sample notes: {updated.notes}")

        return sample1.id, sample2.id, sample3.id
    finally:
        db.close()


def test_link_alignments(sample_ids):
    separator("TEST 2: Link Alignments to Samples & Build Variant Spectrum")
    db = SessionLocal()
    try:
        id1, id2, id3 = sample_ids
        sample1 = sample_service.get_sample_by_id(db, id1)
        sample2 = sample_service.get_sample_by_id(db, id2)
        sample3 = sample_service.get_sample_by_id(db, id3)

        alignments = db.query(models.AlignmentResult).all()
        alignment_ids = [a.id for a in alignments]
        print(f"Found {len(alignment_ids)} alignment results in database")

        half = len(alignment_ids) // 2
        ids_for_s1 = alignment_ids[:max(1, half)]
        ids_for_s2 = alignment_ids[max(1, half):max(2, len(alignment_ids)-1)]
        ids_for_s3 = alignment_ids[max(2, len(alignment_ids)-1):]

        added1, invalid1 = sample_service.link_alignments_to_sample(db, sample1, ids_for_s1)
        print(f"✓ Sample 1: linked {added1} alignments (invalid: {invalid1})")

        added2, invalid2 = sample_service.link_alignments_to_sample(db, sample2, ids_for_s2)
        print(f"✓ Sample 2: linked {added2} alignments (invalid: {invalid2})")

        added3, invalid3 = sample_service.link_alignments_to_sample(db, sample3, ids_for_s3 + ids_for_s1)
        print(f"✓ Sample 3: linked {added3} alignments (invalid: {invalid3})")

        sample1 = sample_service.get_sample_by_id(db, id1)
        sample2 = sample_service.get_sample_by_id(db, id2)
        sample3 = sample_service.get_sample_by_id(db, id3)

        out1 = sample_service.sample_to_out(db, sample1, with_detail=True)
        out2 = sample_service.sample_to_out(db, sample2, with_detail=True)
        out3 = sample_service.sample_to_out(db, sample3, with_detail=True)
        print(f"  Sample1: {out1.alignment_count} alignments, {out1.variant_count} variants")
        print(f"  Sample2: {out2.alignment_count} alignments, {out2.variant_count} variants")
        print(f"  Sample3: {out3.alignment_count} alignments, {out3.variant_count} variants")

    finally:
        db.close()


def test_sample_spectrum(sample_ids):
    separator("TEST 3: Query Sample Variant Spectrum with Filters")
    db = SessionLocal()
    try:
        id1, id2, id3 = sample_ids
        sample1 = sample_service.get_sample_by_id(db, id1)
        sample2 = sample_service.get_sample_by_id(db, id2)
        sample3 = sample_service.get_sample_by_id(db, id3)

        spectrum = sample_service.get_sample_spectrum(db, sample3)
        print(f"✓ Sample3 full spectrum: {spectrum.total_variants} variants")
        for v in spectrum.variants[:3]:
            print(f"  - pos={v.ref_pos}, type={v.variant_type}, {v.ref_base}→{v.alt_base}, gene={v.gene_name}, impact={v.impact}")
        if len(spectrum.variants) > 3:
            print(f"  ... and {len(spectrum.variants) - 3} more")

        spectrum_exon = sample_service.get_sample_spectrum(db, sample3, variant_type_filter="SNP")
        print(f"✓ Sample3 filtered (SNP only): {spectrum_exon.total_variants} variants")

        spectrum_gene = sample_service.get_sample_spectrum(db, sample3, gene_filter="SAMPLE_GENE")
        print(f"✓ Sample3 filtered (SAMPLE_GENE): {spectrum_gene.total_variants} variants")

    finally:
        db.close()


def test_population_frequency():
    separator("TEST 4: Population Frequency Calculation")
    db = SessionLocal()
    try:
        freq_all = sample_service.get_population_frequency(db, SAMPLE_REFERENCE_NAME)
        print(f"✓ Frequency analysis: {freq_all.total_variant_sites} variant sites across {freq_all.total_samples_analyzed} samples")
        for e in freq_all.entries[:5]:
            print(f"  - pos={e.ref_pos}, {e.ref_base}→{e.alt_base}, freq={e.frequency:.4f}, samples={e.sample_count}/{e.total_samples}, gene={e.gene_name}, impact={e.impact}")
        if len(freq_all.entries) > 5:
            print(f"  ... and {len(freq_all.entries) - 5} more")

        freq_common = sample_service.get_population_frequency(db, SAMPLE_REFERENCE_NAME, min_frequency=0.3)
        print(f"✓ Common variants (freq>=0.3): {freq_common.total_variant_sites} sites")

    finally:
        db.close()


def test_compare_samples(sample_ids):
    separator("TEST 5: Two-Sample Comparison")
    db = SessionLocal()
    try:
        id1, id2, id3 = sample_ids
        sample1 = sample_service.get_sample_by_id(db, id1)
        sample2 = sample_service.get_sample_by_id(db, id2)
        sample3 = sample_service.get_sample_by_id(db, id3)

        result = sample_service.compare_samples(db, sample1, sample3)
        print(f"✓ Compared {result.sample_a_name} vs {result.sample_b_name}")
        s = result.summary
        print(f"  Summary:")
        print(f"    - Only in A: {s.only_a_count}")
        print(f"    - Only in B: {s.only_b_count}")
        print(f"    - Shared: {s.shared_count}")
        print(f"    - Same site different alt: {s.same_site_diff_alt_count}")
        print(f"    - Jaccard similarity: {s.jaccard_similarity:.6f}")

        if result.shared:
            print(f"  Example shared variant: pos={result.shared[0].ref_pos}, {result.shared[0].ref_base}→{result.shared[0].alt_base}")
        if result.only_in_a:
            print(f"  Example only-in-A variant: pos={result.only_in_a[0].ref_pos}, {result.only_in_a[0].ref_base}→{result.only_in_a[0].alt_base}")
        if result.only_in_b:
            print(f"  Example only-in-B variant: pos={result.only_in_b[0].ref_pos}, {result.only_in_b[0].ref_base}→{result.only_in_b[0].alt_base}")
        if result.same_site_different_alt:
            print(f"  Same-site-different-alt example at pos={result.same_site_different_alt[0]['ref_pos']}")

    finally:
        db.close()


def test_hotspot_detection():
    separator("TEST 6: Variant Hotspot Detection")
    db = SessionLocal()
    try:
        hotspots = sample_service.find_hotspots(db, SAMPLE_REFERENCE_NAME, threshold_per_100bp=1.0)
        print(f"✓ Hotspot analysis with threshold 1.0/100bp: {hotspots.hotspot_count} hotspot regions")
        print(f"  Total samples analyzed: {hotspots.total_samples_analyzed}")
        for h in hotspots.hotspots[:3]:
            print(f"  - Gene={h.gene_name}, region={h.region_start}-{h.region_end} ({h.region_length}bp)")
            print(f"    Variants: {h.variant_count}, samples: {h.sample_count}, density: {h.density_per_100bp:.2f}/100bp")
            if h.top_variant:
                tv = h.top_variant
                print(f"    Top variant: pos={tv.ref_pos} {tv.ref_base}→{tv.alt_base}, freq={tv.frequency:.4f}")

        hotspots_strict = sample_service.find_hotspots(db, SAMPLE_REFERENCE_NAME, threshold_per_100bp=5.0)
        print(f"✓ Hotspot analysis with threshold 5.0/100bp: {hotspots_strict.hotspot_count} hotspot regions")

    finally:
        db.close()


def test_delete_sample(sample_ids):
    separator("TEST 7: Delete Sample (Alignment Results Preserved)")
    db = SessionLocal()
    try:
        id1, id2, id3 = sample_ids
        sample1 = sample_service.get_sample_by_id(db, id1)
        sample2 = sample_service.get_sample_by_id(db, id2)
        sample3 = sample_service.get_sample_by_id(db, id3)

        total_alignments_before = db.query(models.AlignmentResult).count()
        print(f"Alignments before deletion: {total_alignments_before}")

        sample2_name = sample2.name
        sample_service.delete_sample(db, sample2)
        print(f"✓ Deleted sample '{sample2_name}'")

        total_alignments_after = db.query(models.AlignmentResult).count()
        print(f"Alignments after deletion: {total_alignments_after}")
        assert total_alignments_before == total_alignments_after, "Alignment results should not be deleted!"
        print(f"✓ Confirmed: alignment results preserved ({total_alignments_after})")

        remaining = sample_service.list_samples(db)
        print(f"Remaining samples: {len(remaining)}")

    finally:
        db.close()


def cleanup():
    separator("CLEANUP")
    db = SessionLocal()
    try:
        db.query(models.SampleAlignmentLink).delete()
        db.query(models.SampleVariantSpectrum).delete()
        db.query(models.PopulationFrequencyCache).delete()
        db.query(models.PopulationFrequencyMeta).delete()
        db.query(models.Sample).delete()
        db.commit()
        print("✓ Cleaned up all test sample data")
    finally:
        db.close()


if __name__ == "__main__":
    print("="*60)
    print("SAMPLE VARIANT SPECTRUM & POPULATION ANALYSIS MODULE - TEST SUITE")
    print("="*60)

    init_database()
    cleanup()

    sample_ids = test_sample_crud()
    test_link_alignments(sample_ids)
    test_sample_spectrum(sample_ids)
    test_population_frequency()
    test_compare_samples(sample_ids)
    test_hotspot_detection()
    test_delete_sample(sample_ids)

    separator("ALL TESTS COMPLETED SUCCESSFULLY!")
