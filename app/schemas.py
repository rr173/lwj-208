from pydantic import BaseModel, Field
from typing import List, Optional
from datetime import datetime


class ReferenceSequenceBase(BaseModel):
    name: str
    sequence: str


class ReferenceSequenceCreate(ReferenceSequenceBase):
    pass


class ReferenceSequenceOut(BaseModel):
    id: int
    name: str
    length: int
    created_at: datetime

    class Config:
        from_attributes = True


class ReferenceSequenceDetail(ReferenceSequenceOut):
    sequence: str


class GeneAnnotationBase(BaseModel):
    seq_name: str
    start: int
    end: int
    feature_type: str
    gene_name: str


class GeneAnnotationCreate(GeneAnnotationBase):
    pass


class GeneAnnotationOut(BaseModel):
    id: int
    start: int
    end: int
    feature_type: str
    gene_name: str

    class Config:
        from_attributes = True


class AlignmentDetail(BaseModel):
    query: str
    match: str
    reference: str


class VariantOut(BaseModel):
    id: int
    variant_type: str
    ref_pos: int
    ref_base: str
    alt_base: str
    feature_type: Optional[str] = None
    gene_name: Optional[str] = None
    impact: Optional[str] = None
    consequence: Optional[str] = None
    aa_ref: Optional[str] = None
    aa_alt: Optional[str] = None

    class Config:
        from_attributes = True


class AlignmentResultOut(BaseModel):
    id: int
    reference_name: str
    score: int
    ref_start: int
    ref_end: int
    query_start: int
    query_end: int
    cigar: str
    alignment: AlignmentDetail
    variants: List[VariantOut] = []

    class Config:
        from_attributes = True


class AlignRequest(BaseModel):
    query_sequence: str = Field(..., max_length=1000)


class BatchAlignRequest(BaseModel):
    query_sequences: List[str] = Field(..., max_length=50)


class BatchTaskOut(BaseModel):
    task_id: str
    status: str
    total: int
    completed: int
    created_at: datetime
    completed_at: Optional[datetime] = None

    class Config:
        from_attributes = True


class BatchTaskResultOut(BatchTaskOut):
    results: List[AlignmentResultOut] = []


class StatsBucket(BaseModel):
    bucket_start: int
    bucket_end: int
    count: int


class VariantStatsOut(BaseModel):
    density_by_position: List[StatsBucket]
    count_by_impact: dict
    count_by_type: dict


class UploadFastaRequest(BaseModel):
    fasta_text: str


class UploadGffRequest(BaseModel):
    gff_text: str
