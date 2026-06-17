from sqlalchemy import Column, Integer, String, Text, ForeignKey, DateTime, Float, JSON, UniqueConstraint, Index, Boolean
from sqlalchemy.orm import relationship
from sqlalchemy.sql import func
from app.database import Base


class ReferenceSequence(Base):
    __tablename__ = "reference_sequences"

    id = Column(Integer, primary_key=True, index=True)
    name = Column(String, unique=True, index=True, nullable=False)
    sequence = Column(Text, nullable=False)
    length = Column(Integer, nullable=False)
    created_at = Column(DateTime(timezone=True), server_default=func.now())

    annotations = relationship("GeneAnnotation", back_populates="reference", cascade="all, delete-orphan")
    alignments = relationship("AlignmentResult", back_populates="reference", cascade="all, delete-orphan")


class GeneAnnotation(Base):
    __tablename__ = "gene_annotations"

    id = Column(Integer, primary_key=True, index=True)
    reference_id = Column(Integer, ForeignKey("reference_sequences.id"), nullable=False)
    start = Column(Integer, nullable=False)
    end = Column(Integer, nullable=False)
    feature_type = Column(String, nullable=False)
    gene_name = Column(String, nullable=False)
    strand = Column(String, default="+")

    reference = relationship("ReferenceSequence", back_populates="annotations")


class AlignmentResult(Base):
    __tablename__ = "alignment_results"

    id = Column(Integer, primary_key=True, index=True)
    reference_id = Column(Integer, ForeignKey("reference_sequences.id"), nullable=False)
    query_sequence = Column(Text, nullable=False)
    query_hash = Column(String, index=True, nullable=False)
    score = Column(Integer, nullable=False)
    ref_start = Column(Integer, nullable=False)
    ref_end = Column(Integer, nullable=False)
    query_start = Column(Integer, nullable=False)
    query_end = Column(Integer, nullable=False)
    cigar = Column(String, nullable=False)
    alignment_query = Column(Text, nullable=False)
    alignment_match = Column(Text, nullable=False)
    alignment_ref = Column(Text, nullable=False)
    created_at = Column(DateTime(timezone=True), server_default=func.now())

    reference = relationship("ReferenceSequence", back_populates="alignments")
    variants = relationship("Variant", back_populates="alignment", cascade="all, delete-orphan")
    batch_task_items = relationship("BatchTaskItem", back_populates="alignment")
    sample_links = relationship("SampleAlignmentLink", back_populates="alignment")


class Variant(Base):
    __tablename__ = "variants"

    id = Column(Integer, primary_key=True, index=True)
    alignment_id = Column(Integer, ForeignKey("alignment_results.id"), nullable=False)
    reference_id = Column(Integer, ForeignKey("reference_sequences.id"), nullable=False)
    variant_type = Column(String, nullable=False)
    ref_pos = Column(Integer, nullable=False)
    ref_base = Column(String, nullable=False)
    alt_base = Column(String, nullable=False)
    feature_type = Column(String)
    gene_name = Column(String)
    impact = Column(String)
    consequence = Column(String)
    codon_ref = Column(String)
    codon_alt = Column(String)
    aa_ref = Column(String)
    aa_alt = Column(String)

    alignment = relationship("AlignmentResult", back_populates="variants")

    __table_args__ = (
        Index("idx_variant_ref_pos", "reference_id", "ref_pos"),
        Index("idx_variant_five_tuple", "reference_id", "ref_pos", "variant_type", "ref_base", "alt_base"),
    )


class BatchTask(Base):
    __tablename__ = "batch_tasks"

    id = Column(Integer, primary_key=True, index=True)
    task_id = Column(String, unique=True, index=True, nullable=False)
    status = Column(String, default="pending")
    total = Column(Integer, default=0)
    completed = Column(Integer, default=0)
    created_at = Column(DateTime(timezone=True), server_default=func.now())
    completed_at = Column(DateTime(timezone=True))

    items = relationship("BatchTaskItem", back_populates="task", cascade="all, delete-orphan", order_by="BatchTaskItem.order_index")


class BatchTaskItem(Base):
    __tablename__ = "batch_task_items"

    id = Column(Integer, primary_key=True, index=True)
    task_id = Column(Integer, ForeignKey("batch_tasks.id"), nullable=False)
    order_index = Column(Integer, default=0)
    query_sequence = Column(Text, nullable=False)
    status = Column(String, default="pending")
    alignment_id = Column(Integer, ForeignKey("alignment_results.id"))

    task = relationship("BatchTask", back_populates="items")
    alignment = relationship("AlignmentResult", back_populates="batch_task_items")


class Sample(Base):
    __tablename__ = "samples"

    id = Column(Integer, primary_key=True, index=True)
    name = Column(String, unique=True, index=True, nullable=False)
    species = Column(String, nullable=False)
    collection_date = Column(DateTime(timezone=True), nullable=True)
    notes = Column(Text, default="")
    created_at = Column(DateTime(timezone=True), server_default=func.now())
    updated_at = Column(DateTime(timezone=True), server_default=func.now(), onupdate=func.now())

    alignment_links = relationship("SampleAlignmentLink", back_populates="sample", cascade="all, delete-orphan")
    spectrum_cache = relationship("SampleVariantSpectrum", back_populates="sample", cascade="all, delete-orphan")


class SampleAlignmentLink(Base):
    __tablename__ = "sample_alignment_links"

    id = Column(Integer, primary_key=True, index=True)
    sample_id = Column(Integer, ForeignKey("samples.id"), nullable=False)
    alignment_id = Column(Integer, ForeignKey("alignment_results.id"), nullable=False)
    linked_at = Column(DateTime(timezone=True), server_default=func.now())

    sample = relationship("Sample", back_populates="alignment_links")
    alignment = relationship("AlignmentResult", back_populates="sample_links")

    __table_args__ = (
        UniqueConstraint("sample_id", "alignment_id", name="uq_sample_alignment"),
    )


class SampleVariantSpectrum(Base):
    __tablename__ = "sample_variant_spectrum"

    id = Column(Integer, primary_key=True, index=True)
    sample_id = Column(Integer, ForeignKey("samples.id"), nullable=False)
    reference_id = Column(Integer, ForeignKey("reference_sequences.id"), nullable=False)
    ref_pos = Column(Integer, nullable=False)
    variant_type = Column(String, nullable=False)
    ref_base = Column(String, nullable=False)
    alt_base = Column(String, nullable=False)
    feature_type = Column(String)
    gene_name = Column(String)
    impact = Column(String)
    consequence = Column(String)
    aa_ref = Column(String)
    aa_alt = Column(String)
    source_alignment_ids = Column(JSON, default=list)

    sample = relationship("Sample", back_populates="spectrum_cache")

    __table_args__ = (
        UniqueConstraint(
            "sample_id", "reference_id", "ref_pos", "variant_type", "ref_base", "alt_base",
            name="uq_sample_variant_unique"
        ),
        Index("idx_spectrum_sample_ref", "sample_id", "reference_id"),
        Index("idx_spectrum_sample_gene", "sample_id", "gene_name"),
        Index("idx_spectrum_sample_impact", "sample_id", "impact"),
        Index("idx_spectrum_ref_pos", "reference_id", "ref_pos"),
    )


class PopulationFrequencyCache(Base):
    __tablename__ = "population_frequency_cache"

    id = Column(Integer, primary_key=True, index=True)
    reference_id = Column(Integer, ForeignKey("reference_sequences.id"), nullable=False)
    ref_pos = Column(Integer, nullable=False)
    variant_type = Column(String, nullable=False)
    ref_base = Column(String, nullable=False)
    alt_base = Column(String, nullable=False)
    sample_count = Column(Integer, nullable=False, default=0)
    total_samples = Column(Integer, nullable=False, default=0)
    frequency = Column(Float, nullable=False, default=0.0)
    impact = Column(String)
    gene_name = Column(String)
    last_updated = Column(DateTime(timezone=True), server_default=func.now())

    __table_args__ = (
        UniqueConstraint(
            "reference_id", "ref_pos", "variant_type", "ref_base", "alt_base",
            name="uq_population_frequency_unique"
        ),
        Index("idx_freq_ref_pos", "reference_id", "ref_pos"),
    )


class PopulationFrequencyMeta(Base):
    __tablename__ = "population_frequency_meta"

    id = Column(Integer, primary_key=True, index=True)
    reference_id = Column(Integer, ForeignKey("reference_sequences.id"), unique=True, nullable=False)
    total_samples_analyzed = Column(Integer, nullable=False, default=0)
    last_updated = Column(DateTime(timezone=True), server_default=func.now())


class PhyloTreeTask(Base):
    __tablename__ = "phylo_tree_tasks"

    id = Column(Integer, primary_key=True, index=True)
    task_id = Column(String, unique=True, index=True, nullable=False)
    status = Column(String, default="pending")
    sample_ids = Column(JSON, nullable=False)
    reference_name = Column(String, nullable=False)
    tree_method = Column(String, default="nj")
    bootstrap_replicates = Column(Integer, default=0)
    with_molecular_clock = Column(Boolean, default=False)
    total_steps = Column(Integer, default=100)
    completed_steps = Column(Integer, default=0)
    progress_message = Column(String, default="")
    result_id = Column(Integer, ForeignKey("phylo_tree_results.id"), nullable=True)
    error_message = Column(Text, nullable=True)
    created_at = Column(DateTime(timezone=True), server_default=func.now())
    completed_at = Column(DateTime(timezone=True), nullable=True)
    is_async = Column(Boolean, default=False)

    result = relationship("PhyloTreeResult", back_populates="task", uselist=False)


class PhyloTreeResult(Base):
    __tablename__ = "phylo_tree_results"

    id = Column(Integer, primary_key=True, index=True)
    cache_key = Column(String, unique=True, index=True, nullable=False)
    sample_ids = Column(JSON, nullable=False)
    reference_name = Column(String, nullable=False)
    tree_method = Column(String, nullable=False)
    bootstrap_replicates = Column(Integer, default=0)
    with_molecular_clock = Column(Boolean, default=False)
    newick_tree = Column(Text, nullable=False)
    tree_json = Column(JSON, nullable=False)
    distance_matrix = Column(JSON, nullable=False)
    sample_names = Column(JSON, nullable=False)
    distance_warnings = Column(JSON, nullable=True)
    molecular_clock_info = Column(JSON, nullable=True)
    data_hash = Column(String, nullable=False)
    is_stale = Column(Boolean, default=False)
    created_at = Column(DateTime(timezone=True), server_default=func.now())
    last_checked_at = Column(DateTime(timezone=True), server_default=func.now())

    task = relationship("PhyloTreeTask", back_populates="result", uselist=False)


class ScoringRule(Base):
    __tablename__ = "scoring_rules"

    id = Column(Integer, primary_key=True, index=True)
    name = Column(String, nullable=False)
    description = Column(Text, default="")
    condition_field = Column(String, nullable=False)
    condition_operator = Column(String, nullable=False)
    condition_value = Column(String, nullable=False)
    weight = Column(Float, nullable=False)
    enabled = Column(Boolean, default=True)
    created_at = Column(DateTime(timezone=True), server_default=func.now())
    updated_at = Column(DateTime(timezone=True), server_default=func.now(), onupdate=func.now())


class RuleSetVersion(Base):
    __tablename__ = "rule_set_versions"

    id = Column(Integer, primary_key=True, index=True)
    version = Column(Integer, nullable=False, default=1)
    updated_at = Column(DateTime(timezone=True), server_default=func.now())


class SamplePathogenicityScore(Base):
    __tablename__ = "sample_pathogenicity_scores"

    id = Column(Integer, primary_key=True, index=True)
    sample_id = Column(Integer, ForeignKey("samples.id"), nullable=False, unique=True)
    total_score = Column(Float, nullable=False)
    classification = Column(String, nullable=False)
    rule_version = Column(Integer, nullable=False)
    variant_scores = Column(JSON, nullable=False, default=list)
    scored_at = Column(DateTime(timezone=True), server_default=func.now())

    sample = relationship("Sample", backref="pathogenicity_score")


class TreeComparisonResult(Base):
    __tablename__ = "tree_comparison_results"

    id = Column(Integer, primary_key=True, index=True)
    cache_key = Column(String, unique=True, index=True, nullable=False)
    task_id_a = Column(String, nullable=False)
    task_id_b = Column(String, nullable=False)
    rf_distance = Column(Integer, nullable=False)
    normalized_rf_distance = Column(Float, nullable=False)
    matching_splits = Column(Integer, nullable=False)
    total_splits = Column(Integer, nullable=False)
    weight_dist = Column(Float, nullable=False)
    inconsistent_branches = Column(JSON, nullable=True)
    created_at = Column(DateTime(timezone=True), server_default=func.now())


class LDMatrixResult(Base):
    __tablename__ = "ld_matrix_results"

    id = Column(Integer, primary_key=True, index=True)
    cache_key = Column(String, unique=True, index=True, nullable=False)
    reference_name = Column(String, nullable=False)
    reference_id = Column(Integer, ForeignKey("reference_sequences.id"), nullable=False)
    total_samples = Column(Integer, nullable=False, default=0)
    snp_count = Column(Integer, nullable=False, default=0)
    snp_positions = Column(JSON, nullable=False, default=list)
    ld_matrix = Column(JSON, nullable=False, default=list)
    data_hash = Column(String, nullable=False)
    is_stale = Column(Boolean, default=False)
    created_at = Column(DateTime(timezone=True), server_default=func.now())
    last_checked_at = Column(DateTime(timezone=True), server_default=func.now())

    __table_args__ = (
        Index("idx_ld_matrix_ref", "reference_id"),
    )


class PrimerDesignCache(Base):
    __tablename__ = "primer_design_cache"

    id = Column(Integer, primary_key=True, index=True)
    cache_key = Column(String, unique=True, index=True, nullable=False)
    reference_name = Column(String, nullable=False)
    target_start = Column(Integer, nullable=False)
    target_end = Column(Integer, nullable=False)
    result_json = Column(Text, nullable=False)
    created_at = Column(DateTime(timezone=True), server_default=func.now())


class HaplotypeBlockResult(Base):
    __tablename__ = "haplotype_block_results"

    id = Column(Integer, primary_key=True, index=True)
    cache_key = Column(String, unique=True, index=True, nullable=False)
    reference_name = Column(String, nullable=False)
    reference_id = Column(Integer, ForeignKey("reference_sequences.id"), nullable=False)
    r2_threshold = Column(Float, nullable=False)
    total_samples = Column(Integer, nullable=False, default=0)
    snp_count = Column(Integer, nullable=False, default=0)
    block_count = Column(Integer, nullable=False, default=0)
    blocks = Column(JSON, nullable=False, default=list)
    data_hash = Column(String, nullable=False)
    is_stale = Column(Boolean, default=False)
    created_at = Column(DateTime(timezone=True), server_default=func.now())
    last_checked_at = Column(DateTime(timezone=True), server_default=func.now())

    __table_args__ = (
        Index("idx_haplotype_block_ref", "reference_id"),
    )


class ProteinDomain(Base):
    __tablename__ = "protein_domains"

    id = Column(Integer, primary_key=True, index=True)
    gene_name = Column(String, nullable=False, index=True)
    domain_name = Column(String, nullable=False)
    start_codon = Column(Integer, nullable=False)
    end_codon = Column(Integer, nullable=False)
    domain_type = Column(String, nullable=False)

    __table_args__ = (
        Index("idx_protein_domain_gene", "gene_name"),
    )


class SyntenyResult(Base):
    __tablename__ = "synteny_results"

    id = Column(Integer, primary_key=True, index=True)
    cache_key = Column(String, unique=True, index=True, nullable=False)
    ref_a_name = Column(String, nullable=False)
    ref_b_name = Column(String, nullable=False)
    ref_a_id = Column(Integer, ForeignKey("reference_sequences.id"), nullable=False)
    ref_b_id = Column(Integer, ForeignKey("reference_sequences.id"), nullable=False)
    seq_a_length = Column(Integer, nullable=False)
    seq_b_length = Column(Integer, nullable=False)
    anchor_length = Column(Integer, nullable=False)
    score_threshold_ratio = Column(Float, nullable=False)
    max_gap_ratio = Column(Float, nullable=False, default=3.0)
    total_anchors_aligned = Column(Integer, nullable=False, default=0)
    synteny_blocks = Column(JSON, nullable=False, default=list)
    rearrangements = Column(JSON, nullable=False, default=list)
    b_non_monotonic_transitions = Column(JSON, nullable=False, default=list)
    data_hash = Column(String, nullable=False)
    is_stale = Column(Boolean, default=False)
    created_at = Column(DateTime(timezone=True), server_default=func.now())
    last_checked_at = Column(DateTime(timezone=True), server_default=func.now())

    ref_a = relationship("ReferenceSequence", foreign_keys=[ref_a_id])
    ref_b = relationship("ReferenceSequence", foreign_keys=[ref_b_id])

    __table_args__ = (
        Index("idx_synteny_refs", "ref_a_id", "ref_b_id"),
    )


class QCRuleSet(Base):
    __tablename__ = "qc_rule_sets"

    id = Column(Integer, primary_key=True, index=True)
    name = Column(String, nullable=False)
    description = Column(Text, default="")
    is_active = Column(Boolean, default=False, index=True)
    created_at = Column(DateTime(timezone=True), server_default=func.now())
    updated_at = Column(DateTime(timezone=True), server_default=func.now(), onupdate=func.now())

    rules = relationship("QCRule", back_populates="rule_set", cascade="all, delete-orphan", order_by="QCRule.rule_index")

    __table_args__ = (
        Index("idx_qc_ruleset_active", "is_active"),
    )


class QCRule(Base):
    __tablename__ = "qc_rules"

    id = Column(Integer, primary_key=True, index=True)
    rule_set_id = Column(Integer, ForeignKey("qc_rule_sets.id"), nullable=False)
    rule_index = Column(Integer, nullable=False)
    rule_type = Column(String, nullable=False)
    min_alignment_score = Column(Integer, nullable=True)
    min_coverage_depth = Column(Integer, nullable=True)
    edge_margin = Column(Integer, nullable=True)
    created_at = Column(DateTime(timezone=True), server_default=func.now())

    rule_set = relationship("QCRuleSet", back_populates="rules")

    __table_args__ = (
        Index("idx_qc_rule_set", "rule_set_id"),
    )


class VariantQCResult(Base):
    __tablename__ = "variant_qc_results"

    id = Column(Integer, primary_key=True, index=True)
    sample_id = Column(Integer, ForeignKey("samples.id"), nullable=False)
    reference_id = Column(Integer, ForeignKey("reference_sequences.id"), nullable=False)
    ref_pos = Column(Integer, nullable=False)
    variant_type = Column(String, nullable=False)
    ref_base = Column(String, nullable=False)
    alt_base = Column(String, nullable=False)
    rule_set_id = Column(Integer, ForeignKey("qc_rule_sets.id"), nullable=False)
    qc_status = Column(String, nullable=False, index=True)
    failed_rules = Column(JSON, default=list)
    is_stale = Column(Boolean, default=False, index=True)
    created_at = Column(DateTime(timezone=True), server_default=func.now())

    __table_args__ = (
        UniqueConstraint(
            "sample_id", "reference_id", "ref_pos", "variant_type", "ref_base", "alt_base",
            name="uq_variant_qc_unique"
        ),
        Index("idx_qc_result_sample", "sample_id"),
        Index("idx_qc_result_sample_status", "sample_id", "qc_status"),
        Index("idx_qc_result_sample_stale", "sample_id", "is_stale"),
    )


class AuditLog(Base):
    __tablename__ = "audit_logs"

    id = Column(Integer, primary_key=True, index=True)
    timestamp = Column(DateTime(timezone=True), server_default=func.now(), index=True)
    operation_type = Column(String, index=True, nullable=False)
    resource_type = Column(String, index=True, nullable=False)
    resource_id = Column(String, index=True, nullable=True)
    summary = Column(String, nullable=False)
    method = Column(String, nullable=False)
    path = Column(String, nullable=False)
    status_code = Column(Integer, nullable=False)
    request_body = Column(Text, nullable=True)

    __table_args__ = (
        Index("idx_audit_resource", "resource_type", "resource_id"),
        Index("idx_audit_timestamp", "timestamp"),
    )


class MutationalSpectrumCache(Base):
    __tablename__ = "mutational_spectrum_cache"

    id = Column(Integer, primary_key=True, index=True)
    cache_key = Column(String, unique=True, index=True, nullable=False)
    sample_ids = Column(JSON, nullable=False)
    reference_name = Column(String, nullable=False)
    reference_id = Column(Integer, ForeignKey("reference_sequences.id"), nullable=False)
    sample_names = Column(JSON, nullable=False)
    mutation_types = Column(JSON, nullable=False)
    count_matrix = Column(JSON, nullable=False)
    data_hash = Column(String, nullable=False)
    is_stale = Column(Boolean, default=False)
    created_at = Column(DateTime(timezone=True), server_default=func.now())
    last_checked_at = Column(DateTime(timezone=True), server_default=func.now())

    reference = relationship("ReferenceSequence", foreign_keys=[reference_id])

    __table_args__ = (
        Index("idx_spectrum_cache_ref", "reference_id"),
    )


class NMFSignatureResult(Base):
    __tablename__ = "nmf_signature_results"

    id = Column(Integer, primary_key=True, index=True)
    cache_key = Column(String, unique=True, index=True, nullable=False)
    sample_ids = Column(JSON, nullable=False)
    reference_name = Column(String, nullable=False)
    k_value = Column(Integer, nullable=False)
    signature_matrix = Column(JSON, nullable=False)
    exposure_matrix = Column(JSON, nullable=False)
    mutation_types = Column(JSON, nullable=False)
    sample_names = Column(JSON, nullable=False)
    reconstruction_error = Column(Float, nullable=False)
    iterations = Column(Integer, nullable=False)
    data_hash = Column(String, nullable=False)
    is_stale = Column(Boolean, default=False)
    created_at = Column(DateTime(timezone=True), server_default=func.now())
    last_checked_at = Column(DateTime(timezone=True), server_default=func.now())


class OptimalKResult(Base):
    __tablename__ = "optimal_k_results"

    id = Column(Integer, primary_key=True, index=True)
    cache_key = Column(String, unique=True, index=True, nullable=False)
    sample_ids = Column(JSON, nullable=False)
    reference_name = Column(String, nullable=False)
    k_values = Column(JSON, nullable=False)
    reconstruction_errors = Column(JSON, nullable=False)
    cophenetic_correlations = Column(JSON, nullable=False)
    recommended_k = Column(Integer, nullable=False)
    data_hash = Column(String, nullable=False)
    is_stale = Column(Boolean, default=False)
    created_at = Column(DateTime(timezone=True), server_default=func.now())
    last_checked_at = Column(DateTime(timezone=True), server_default=func.now())


class ReferenceSignature(Base):
    __tablename__ = "reference_signatures"

    id = Column(Integer, primary_key=True, index=True)
    signature_id = Column(String, unique=True, index=True, nullable=False)
    name = Column(String, nullable=False)
    description = Column(Text, nullable=False)
    etiology = Column(String, nullable=False)
    mutation_types = Column(JSON, nullable=False)
    probabilities = Column(JSON, nullable=False)
    created_at = Column(DateTime(timezone=True), server_default=func.now())


VALID_PIPELINE_STEP_TYPES = [
    "spectrum", "nmf_extract", "scoring", "ld_analysis", "temporal_analysis"
]

PIPELINE_STATUSES = ["pending", "running", "completed", "skipped", "failed"]


class PipelineTemplate(Base):
    __tablename__ = "pipeline_templates"

    id = Column(Integer, primary_key=True, index=True)
    name = Column(String, unique=True, index=True, nullable=False)
    description = Column(Text, default="")
    steps = Column(JSON, nullable=False)
    created_at = Column(DateTime(timezone=True), server_default=func.now())
    updated_at = Column(DateTime(timezone=True), server_default=func.now(), onupdate=func.now())

    executions = relationship("PipelineExecution", back_populates="template", cascade="all, delete-orphan")


class PipelineExecution(Base):
    __tablename__ = "pipeline_executions"

    id = Column(Integer, primary_key=True, index=True)
    template_id = Column(Integer, ForeignKey("pipeline_templates.id"), nullable=False)
    initial_params = Column(JSON, nullable=False)
    status = Column(String, default="pending", index=True)
    total_duration_seconds = Column(Float, nullable=True)
    error_message = Column(Text, nullable=True)
    triggered_at = Column(DateTime(timezone=True), server_default=func.now())
    completed_at = Column(DateTime(timezone=True), nullable=True)
    template_snapshot = Column(JSON, nullable=True)
    resume_count = Column(Integer, default=0)

    template = relationship("PipelineTemplate", back_populates="executions")
    step_executions = relationship("PipelineStepExecution", back_populates="execution", cascade="all, delete-orphan", order_by="PipelineStepExecution.step_index")


class PipelineStepExecution(Base):
    __tablename__ = "pipeline_step_executions"

    id = Column(Integer, primary_key=True, index=True)
    execution_id = Column(Integer, ForeignKey("pipeline_executions.id"), nullable=False)
    step_index = Column(Integer, nullable=False)
    step_type = Column(String, nullable=False)
    step_name = Column(String, nullable=False)
    input_params = Column(JSON, nullable=False)
    status = Column(String, default="pending", index=True)
    output_summary = Column(JSON, nullable=True)
    output_data = Column(JSON, nullable=True)
    error_message = Column(Text, nullable=True)
    started_at = Column(DateTime(timezone=True), nullable=True)
    completed_at = Column(DateTime(timezone=True), nullable=True)
    duration_seconds = Column(Float, nullable=True)
    retry_count = Column(Integer, default=0)
    last_error = Column(Text, nullable=True)

    execution = relationship("PipelineExecution", back_populates="step_executions")
