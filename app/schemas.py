from pydantic import BaseModel, Field
from typing import List, Optional, Dict, Any
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


class SampleCreate(BaseModel):
    name: str = Field(..., min_length=1, max_length=255)
    species: str = Field(..., min_length=1, max_length=255)
    collection_date: Optional[datetime] = None
    notes: Optional[str] = ""


class SampleUpdate(BaseModel):
    name: Optional[str] = None
    species: Optional[str] = None
    collection_date: Optional[datetime] = None
    notes: Optional[str] = None


class SampleOut(BaseModel):
    id: int
    name: str
    species: str
    collection_date: Optional[datetime] = None
    notes: str
    created_at: datetime
    updated_at: datetime
    alignment_count: int = 0
    variant_count: int = 0

    class Config:
        from_attributes = True


class SampleDetail(SampleOut):
    linked_alignment_ids: List[int] = []


class LinkAlignmentsRequest(BaseModel):
    alignment_ids: List[int] = Field(..., min_length=1)


class SpectrumVariantOut(BaseModel):
    reference_name: str
    ref_pos: int
    variant_type: str
    ref_base: str
    alt_base: str
    feature_type: Optional[str] = None
    gene_name: Optional[str] = None
    impact: Optional[str] = None
    consequence: Optional[str] = None
    aa_ref: Optional[str] = None
    aa_alt: Optional[str] = None
    source_alignment_ids: List[int] = []

    class Config:
        from_attributes = True


class SampleSpectrumOut(BaseModel):
    sample_id: int
    sample_name: str
    total_variants: int
    variants: List[SpectrumVariantOut] = []


class PopulationFrequencyEntry(BaseModel):
    ref_pos: int
    variant_type: str
    ref_base: str
    alt_base: str
    frequency: float
    sample_count: int
    total_samples: int
    impact: Optional[str] = None
    gene_name: Optional[str] = None


class PopulationFrequencyOut(BaseModel):
    reference_name: str
    total_samples_analyzed: int
    total_variant_sites: int
    entries: List[PopulationFrequencyEntry] = []


class CompareSummary(BaseModel):
    only_a_count: int
    only_b_count: int
    shared_count: int
    same_site_diff_alt_count: int
    jaccard_similarity: float


class CompareResultOut(BaseModel):
    sample_a_id: int
    sample_a_name: str
    sample_b_id: int
    sample_b_name: str
    summary: CompareSummary
    only_in_a: List[SpectrumVariantOut] = []
    only_in_b: List[SpectrumVariantOut] = []
    shared: List[SpectrumVariantOut] = []
    same_site_different_alt: List[Dict[str, Any]] = []


class HotspotEntry(BaseModel):
    region_start: int
    region_end: int
    region_length: int
    gene_name: str
    feature_type: Optional[str] = None
    variant_count: int
    sample_count: int
    density_per_100bp: float
    top_variant: Optional[PopulationFrequencyEntry] = None


class HotspotAnalysisOut(BaseModel):
    reference_name: str
    total_samples_analyzed: int
    threshold_per_100bp: float
    hotspot_count: int
    hotspots: List[HotspotEntry] = []
