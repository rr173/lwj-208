from sqlalchemy import Column, Integer, String, Text, ForeignKey, DateTime, Float, JSON
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
