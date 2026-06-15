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
