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


class DistanceMatrixRequest(BaseModel):
    sample_ids: List[int] = Field(..., min_length=3, max_length=200)
    reference_name: str = Field(..., min_length=1)


DistanceWarning = Dict[str, Any]


class DistanceMatrixOut(BaseModel):
    sample_ids: List[int]
    sample_names: List[str]
    distance_matrix: List[List[float]]
    warnings: List[DistanceWarning] = []


class PhyloTreeRequest(BaseModel):
    sample_ids: List[int] = Field(..., min_length=3, max_length=200)
    reference_name: str = Field(..., min_length=1)
    tree_method: str = Field("nj", pattern="^(nj|upgma)$")
    bootstrap_replicates: int = Field(0, ge=0, le=1000)
    with_molecular_clock: bool = False


class PhyloTreeNode(BaseModel):
    name: Optional[str] = None
    children: List["PhyloTreeNode"] = []
    branch_length: Optional[float] = None
    bootstrap_support: Optional[float] = None
    divergence_time: Optional[float] = None
    is_leaf: bool = False


PhyloTreeNode.model_rebuild()


class MolecularClockInfo(BaseModel):
    slope: float
    intercept: float
    r_squared: float
    rate: float
    outlier_samples: List[str] = []
    warning: Optional[str] = None


class PhyloTreeOut(BaseModel):
    task_id: str
    status: str
    newick_tree: Optional[str] = None
    tree_json: Optional[PhyloTreeNode] = None
    sample_ids: List[int] = []
    sample_names: List[str] = []
    tree_method: str
    bootstrap_replicates: int
    distance_matrix: Optional[List[List[float]]] = []
    distance_warnings: List[DistanceWarning] = []
    molecular_clock: Optional[MolecularClockInfo] = None
    is_stale: bool = False
    created_at: Optional[datetime] = None
    progress: Optional[int] = None
    total: Optional[int] = None
    progress_message: Optional[str] = None
    error_message: Optional[str] = None


class PhyloTaskStatusOut(BaseModel):
    task_id: str
    status: str
    progress: int
    total: int
    progress_message: str
    is_async: bool
    created_at: datetime
    completed_at: Optional[datetime] = None
    error_message: Optional[str] = None


class TreeCompareRequest(BaseModel):
    task_id_a: str
    task_id_b: str


class InconsistentBranch(BaseModel):
    split_a: List[str]
    split_b: List[str]


class TreeCompareOut(BaseModel):
    rf_distance: int
    normalized_rf_distance: float
    matching_splits: int
    total_splits: int
    weighted_branch_length_distance: float
    inconsistent_branches: List[InconsistentBranch] = []
    sample_names_a: List[str] = []
    sample_names_b: List[str] = []


VALID_CONDITION_FIELDS = [
    "impact", "feature_type", "variant_type", "consequence", "population_frequency"
]
VALID_CONDITION_OPERATORS = ["eq", "ne", "gt", "gte", "lt", "lte"]
VALID_CLASSIFICATIONS = [
    "pathogenic", "likely_pathogenic", "uncertain_significance",
    "likely_benign", "benign",
]


class ScoringRuleCreate(BaseModel):
    name: str = Field(..., min_length=1, max_length=255)
    description: Optional[str] = ""
    condition_field: str = Field(..., pattern="|".join(VALID_CONDITION_FIELDS))
    condition_operator: str = Field(..., pattern="|".join(VALID_CONDITION_OPERATORS))
    condition_value: str = Field(..., min_length=1)
    weight: float
    enabled: bool = True


class ScoringRuleUpdate(BaseModel):
    name: Optional[str] = Field(None, min_length=1, max_length=255)
    description: Optional[str] = None
    condition_field: Optional[str] = Field(None, pattern="|".join(VALID_CONDITION_FIELDS))
    condition_operator: Optional[str] = Field(None, pattern="|".join(VALID_CONDITION_OPERATORS))
    condition_value: Optional[str] = Field(None, min_length=1)
    weight: Optional[float] = None
    enabled: Optional[bool] = None


class ScoringRuleOut(BaseModel):
    id: int
    name: str
    description: str
    condition_field: str
    condition_operator: str
    condition_value: str
    weight: float
    enabled: bool
    created_at: datetime
    updated_at: datetime

    class Config:
        from_attributes = True


class VariantScoreDetail(BaseModel):
    reference_id: int
    ref_pos: int
    variant_type: str
    ref_base: str
    alt_base: str
    feature_type: Optional[str] = None
    gene_name: Optional[str] = None
    impact: Optional[str] = None
    consequence: Optional[str] = None
    population_frequency: Optional[float] = None
    matched_rules: List[Dict[str, Any]] = []
    variant_score: float
    variant_classification: str


class SampleScoreOut(BaseModel):
    sample_id: int
    sample_name: str
    total_score: float
    classification: str
    rule_version: int
    current_rule_version: int
    score_may_be_outdated: bool
    total_variants: int
    variant_details: List[VariantScoreDetail] = []
    scored_at: datetime


class SampleScoreSummaryOut(BaseModel):
    sample_id: int
    sample_name: str
    total_score: float
    classification: str
    rule_version: int
    current_rule_version: int
    score_may_be_outdated: bool
    total_variants: int
    scored_at: datetime


class FilteredVariantOut(BaseModel):
    reference_id: int
    ref_pos: int
    variant_type: str
    ref_base: str
    alt_base: str
    feature_type: Optional[str] = None
    gene_name: Optional[str] = None
    impact: Optional[str] = None
    consequence: Optional[str] = None
    population_frequency: Optional[float] = None
    variant_score: float
    variant_classification: str
    matched_rules: List[Dict[str, Any]] = []


class ClassificationCount(BaseModel):
    classification: str
    count: int


class ScoreStatsOut(BaseModel):
    sample_id: int
    sample_name: str
    total_variants: int
    classification_distribution: List[ClassificationCount]
