from pydantic import BaseModel, Field, model_validator
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


class LDPair(BaseModel):
    pos_i: int
    pos_j: int
    d: float
    d_prime: float
    r_squared: float


class LDMatrixOut(BaseModel):
    reference_name: str
    total_samples: int
    snp_count: int
    original_snp_count: int
    filtered_snp_count: int
    snp_positions: List[int]
    ld_pairs: List[LDPair]
    start: Optional[int] = None
    end: Optional[int] = None
    min_maf: Optional[float] = None
    is_stale: bool
    created_at: Optional[datetime] = None


class HaplotypeBlock(BaseModel):
    block_index: int
    start_pos: int
    end_pos: int
    snp_count: int
    snp_positions: List[int]
    avg_r_squared: float


class HaplotypeBlockOut(BaseModel):
    reference_name: str
    r2_threshold: float
    total_samples: int
    snp_count: int
    original_snp_count: int
    filtered_snp_count: int
    block_count: int
    blocks: List[HaplotypeBlock]
    start: Optional[int] = None
    end: Optional[int] = None
    min_maf: Optional[float] = None
    is_stale: bool
    created_at: Optional[datetime] = None


class PrimerDesignRequest(BaseModel):
    reference_name: str = Field(..., min_length=1)
    target_start: int = Field(..., ge=0)
    target_end: int = Field(..., ge=0)


class BatchPrimerDesignItem(BaseModel):
    target_start: int = Field(..., ge=0)
    target_end: int = Field(..., ge=0)


class BatchPrimerDesignRequest(BaseModel):
    reference_name: str = Field(..., min_length=1)
    targets: List[BatchPrimerDesignItem] = Field(..., min_length=1, max_length=20)


class OffTargetSite(BaseModel):
    ref_start: int
    ref_end: int
    mismatches: int


class NonspecificAmplicon(BaseModel):
    fwd_site: Dict[str, Any]
    rev_site: Dict[str, Any]
    amplicon_start: int
    amplicon_end: int
    amplicon_length: int


class PrimerInfo(BaseModel):
    sequence: str
    length: int
    gc_percent: float
    tm: float
    bind_start: int
    bind_end: int
    off_target_sites: List[OffTargetSite] = []


class PrimerPairOut(BaseModel):
    rank: int
    forward_primer: PrimerInfo
    reverse_primer: PrimerInfo
    amplicon_length: int
    avg_tm_deviation: float
    tm_difference: float
    has_nonspecific_risk: bool
    nonspecific_amplicons: List[NonspecificAmplicon] = []


class PrimerDesignResult(BaseModel):
    reference_name: str
    target_start: int
    target_end: int
    flank_size: int
    primer_pairs: List[PrimerPairOut]
    cached: bool = False


class PrimerTmStats(BaseModel):
    all_fwd_tm: List[float]
    all_rev_tm: List[float]
    fwd_tm_std: float
    rev_tm_std: float
    combined_tm_std: float


class BatchPrimerDesignResultItem(BaseModel):
    target_index: int
    target_start: int
    target_end: int
    best_primer_pair: Optional[PrimerPairOut] = None
    success: bool
    error_message: Optional[str] = None
    cached: bool = False


class BatchPrimerDesignResult(BaseModel):
    reference_name: str
    total_targets: int
    success_count: int
    failure_count: int
    results: List[BatchPrimerDesignResultItem]
    tm_consistency: Optional[PrimerTmStats] = None


class PrimerDesignHistoryItem(BaseModel):
    id: int
    reference_name: str
    target_start: int
    target_end: int
    best_primer_pair: Optional[PrimerPairOut] = None
    created_at: datetime


class PrimerDesignHistoryList(BaseModel):
    reference_name: str
    total_records: int
    history: List[PrimerDesignHistoryItem]


class VariantTracingRequest(BaseModel):
    reference_name: str = Field(..., min_length=1)
    ref_pos: int = Field(..., ge=1)
    variant_type: str = Field(..., min_length=1)
    ref_base: str = Field(..., min_length=1)
    alt_base: str = Field(..., min_length=1)


class VariantTracingSampleEntry(BaseModel):
    sample_id: int
    sample_name: str
    collection_date: Optional[datetime] = None
    is_earliest: bool = False
    is_candidate_source: bool = False
    date_unknown: bool = False


class VariantTracingOut(BaseModel):
    reference_name: str
    ref_pos: int
    variant_type: str
    ref_base: str
    alt_base: str
    total_carrier_samples: int
    timeline: List[VariantTracingSampleEntry] = []


class TransmissionChainRequest(BaseModel):
    sample_ids: List[int] = Field(..., min_length=3, max_length=200)


class TransmissionNodeOut(BaseModel):
    sample_id: int
    sample_name: str
    collection_date: Optional[datetime] = None
    variant_count: int
    date_unknown: bool = False


class TransmissionEdgeOut(BaseModel):
    source_id: int
    source_name: str
    target_id: int
    target_name: str
    new_variants: List[SpectrumVariantOut] = []


class TransmissionChainOut(BaseModel):
    sample_ids: List[int]
    total_nodes: int
    total_edges: int
    nodes: List[TransmissionNodeOut] = []
    edges: List[TransmissionEdgeOut] = []
    excluded_no_date: List[int] = []


class ClusteringRequest(BaseModel):
    reference_name: str = Field(..., min_length=1)
    similarity_threshold: float = Field(0.7, ge=0.0, le=1.0)


class ClusterMemberOut(BaseModel):
    sample_id: int
    sample_name: str
    collection_date: Optional[datetime] = None
    variant_count: int


class ClusterOut(BaseModel):
    cluster_index: int
    member_count: int
    members: List[ClusterMemberOut] = []
    shared_variants: List[SpectrumVariantOut] = []
    avg_jaccard: float


class ClusteringOut(BaseModel):
    reference_name: str
    total_samples: int
    similarity_threshold: float
    cluster_count: int
    clusters: List[ClusterOut] = []


VALID_DOMAIN_TYPES = ["catalytic", "binding", "structural", "regulatory"]


class ProteinDomainItem(BaseModel):
    gene_name: str = Field(..., min_length=1)
    domain_name: str = Field(..., min_length=1)
    start_codon: int = Field(..., ge=1)
    end_codon: int = Field(..., ge=1)
    domain_type: str = Field(..., pattern="|".join(VALID_DOMAIN_TYPES))


class ProteinDomainBatchRequest(BaseModel):
    domains: List[ProteinDomainItem] = Field(..., min_length=1)


class ProteinDomainOut(BaseModel):
    id: int
    gene_name: str
    domain_name: str
    start_codon: int
    end_codon: int
    domain_type: str

    class Config:
        from_attributes = True


class HitDomainInfo(BaseModel):
    domain_name: str
    domain_type: str
    start_codon: int
    end_codon: int


class VariantDomainMappingOut(BaseModel):
    ref_pos: int
    gene_name: str
    codon_position: int
    hit_domains: List[HitDomainInfo] = []
    aa_ref: Optional[str] = None
    aa_alt: Optional[str] = None
    amino_acid_change: Optional[str] = None
    consequence: Optional[str] = None
    variant_type: str
    ref_base: str
    alt_base: str


class VariantImpactOut(BaseModel):
    ref_pos: int
    gene_name: str
    codon_position: int
    hit_domains: List[HitDomainInfo] = []
    aa_ref: Optional[str] = None
    aa_alt: Optional[str] = None
    amino_acid_change: Optional[str] = None
    consequence: Optional[str] = None
    variant_type: str
    ref_base: str
    alt_base: str
    risk_level: str
    risk_reason: Optional[str] = None


class SampleImpactAssessmentOut(BaseModel):
    sample_id: int
    sample_name: str
    total_exon_variants: int
    mapped_variants: int
    variants_with_domain_hit: int
    risk_summary: Dict[str, int] = {}
    variants: List[VariantImpactOut] = []


class DomainVariantCarrier(BaseModel):
    sample_id: int
    sample_name: str


class DomainVariantSummaryEntry(BaseModel):
    codon_position: int
    ref_pos: int
    ref_base: str
    alt_base: str
    variant_type: str
    consequence: Optional[str] = None
    aa_ref: Optional[str] = None
    aa_alt: Optional[str] = None
    amino_acid_change: Optional[str] = None
    risk_level: str
    risk_reason: Optional[str] = None
    carrier_samples: List[DomainVariantCarrier] = []


class DomainSummaryStats(BaseModel):
    total_variants: int
    unique_samples: int
    high_risk_count: int
    medium_risk_count: int
    low_risk_count: int
    unknown_risk_count: int


class DomainVariantSummaryOut(BaseModel):
    gene_name: str
    domain_name: str
    domain_type: str
    start_codon: int
    end_codon: int
    summary_stats: DomainSummaryStats
    variants: List[DomainVariantSummaryEntry] = []


class SyntenyBlockOut(BaseModel):
    a_start: int
    a_end: int
    b_start: int
    b_end: int
    anchor_count: int
    direction: str


class RearrangementEventOut(BaseModel):
    event_type: str
    a_start: int
    a_end: int
    b_start: int
    b_end: int
    anchor_count: int


class SyntenyAnalysisRequest(BaseModel):
    ref_a_name: str = Field(..., min_length=1)
    ref_b_name: str = Field(..., min_length=1)
    anchor_length: int = Field(50, ge=10, le=500)
    score_threshold_ratio: float = Field(1.5, ge=0.5, le=5.0)
    max_gap_ratio: float = Field(3.0, ge=1.0, le=10.0)


class SyntenyAnalysisOut(BaseModel):
    ref_a_name: str
    ref_b_name: str
    seq_a_length: int
    seq_b_length: int
    anchor_length: int
    score_threshold_ratio: float
    max_gap_ratio: float
    total_anchors_aligned: int
    synteny_blocks: List[SyntenyBlockOut] = []
    rearrangements: List[RearrangementEventOut] = []
    b_non_monotonic_transitions: List[int] = []
    is_stale: bool = False
    cached: bool = False
    created_at: Optional[datetime] = None


class SequenceInterval(BaseModel):
    sequence_name: str
    start: int
    end: int


class ConservedCoreRegion(BaseModel):
    region_index: int
    intervals: List[SequenceInterval]
    avg_anchor_score: float
    sequence_count: int
    total_length: int


class MultiSyntenyRequest(BaseModel):
    reference_names: List[str] = Field(..., min_length=3, max_length=10)
    anchor_length: int = Field(50, ge=10, le=500)
    score_threshold_ratio: float = Field(1.5, ge=0.5, le=5.0)
    max_gap_ratio: float = Field(3.0, ge=1.0, le=10.0)


class PairwiseSyntenySummary(BaseModel):
    ref_a_name: str
    ref_b_name: str
    block_count: int
    total_anchors: int
    cached: bool


class MultiSyntenyOut(BaseModel):
    reference_names: List[str]
    first_reference_name: str
    first_reference_length: int
    conserved_core_regions: List[ConservedCoreRegion] = []
    total_conserved_length: int
    conservation_ratio: float
    pairwise_summaries: List[PairwiseSyntenySummary] = []
    anchor_length: int
    score_threshold_ratio: float
    max_gap_ratio: float
    total_pairwise_comparisons: int
    cached_pairwise_count: int


class AuditLogOut(BaseModel):
    id: int
    timestamp: datetime
    operation_type: str
    resource_type: str
    resource_id: Optional[str] = None
    summary: str
    method: str
    path: str
    status_code: int

    class Config:
        from_attributes = True


class AuditLogListOut(BaseModel):
    total: int
    page: int
    page_size: int
    items: List[AuditLogOut] = []


class MutationalSpectrumRequest(BaseModel):
    sample_ids: List[int] = Field(..., min_length=5, max_length=200)
    reference_name: str = Field(..., min_length=1)


class MutationalSpectrumOut(BaseModel):
    reference_name: str
    sample_ids: List[int]
    sample_names: List[str]
    mutation_types: List[str]
    count_matrix: List[List[int]]
    cached: bool = False
    is_stale: bool = False
    created_at: Optional[datetime] = None


class NMFSignatureRequest(BaseModel):
    sample_ids: List[int] = Field(..., min_length=5, max_length=200)
    reference_name: str = Field(..., min_length=1)
    k_value: int = Field(..., ge=2, le=10)


class SignatureOut(BaseModel):
    signature_index: int
    mutation_types: List[str]
    probabilities: List[float]


class SampleExposureOut(BaseModel):
    sample_id: int
    sample_name: str
    signature_exposures: List[float]


class NMFSignatureOut(BaseModel):
    reference_name: str
    k_value: int
    sample_ids: List[int]
    sample_names: List[str]
    signatures: List[SignatureOut]
    exposures: List[SampleExposureOut]
    reconstruction_error: float
    iterations: int
    cached: bool = False
    is_stale: bool = False
    created_at: Optional[datetime] = None


class OptimalKRequest(BaseModel):
    sample_ids: List[int] = Field(..., min_length=5, max_length=200)
    reference_name: str = Field(..., min_length=1)


class KMetricEntry(BaseModel):
    k_value: int
    reconstruction_error: float
    cophenetic_correlation: float


class OptimalKOut(BaseModel):
    reference_name: str
    sample_ids: List[int]
    sample_names: List[str]
    metrics: List[KMetricEntry]
    recommended_k: int
    cached: bool = False
    is_stale: bool = False
    created_at: Optional[datetime] = None


class SignatureMatchRequest(BaseModel):
    signature: List[float] = Field(..., min_length=96, max_length=96)


class ReferenceSignatureOut(BaseModel):
    signature_id: str
    name: str
    description: str
    etiology: str


class SignatureMatchOut(BaseModel):
    reference_signature: ReferenceSignatureOut
    cosine_similarity: float
    is_reliable: bool


class SignatureComparisonOut(BaseModel):
    query_signature: List[float]
    top_matches: List[SignatureMatchOut]


class ReferenceSignatureListOut(BaseModel):
    total: int
    signatures: List[ReferenceSignatureOut] = []


class SignatureTemporalPoint(BaseModel):
    sample_id: int
    sample_name: str
    collection_date: datetime
    exposure: float


class LinearRegressionResult(BaseModel):
    slope: float
    intercept: float
    r_squared: float
    p_value: float
    trend: str


class SignatureTemporalAnalysis(BaseModel):
    signature_index: int
    temporal_points: List[SignatureTemporalPoint]
    regression: LinearRegressionResult


class SignatureTemporalAnalysisOut(BaseModel):
    reference_name: str
    k_value: int
    sample_ids: List[int]
    sample_names: List[str]
    samples_with_date: int
    samples_without_date: int
    excluded_sample_ids: List[int]
    signatures: List[SignatureTemporalAnalysis]
    total_signatures: int
    significant_rising: List[int] = []
    significant_falling: List[int] = []
    p_value_threshold: float = 0.05
    created_at: Optional[datetime] = None


class TemporalAnalysisRequest(BaseModel):
    sample_ids: List[int] = Field(..., min_length=5, max_length=200)
    reference_name: str = Field(..., min_length=1)
    k_value: Optional[int] = Field(None, ge=2, le=10)
    p_value_threshold: float = Field(0.05, gt=0, lt=1)


class SignatureExposureChange(BaseModel):
    signature_index: int
    start_exposure: float
    end_exposure: float
    absolute_change: float
    percent_change: Optional[float] = None
    percent_change_available: bool = True


class ExposureChangeRequest(BaseModel):
    sample_ids: List[int] = Field(..., min_length=5, max_length=200)
    reference_name: str = Field(..., min_length=1)
    start_date: datetime
    end_date: datetime
    k_value: Optional[int] = Field(None, ge=2, le=10)
    aggregation: str = Field("mean", pattern="^(mean|median)$")


class ExposureChangeOut(BaseModel):
    reference_name: str
    k_value: int
    start_date: datetime
    end_date: datetime
    aggregation: str
    start_sample_count: int
    end_sample_count: int
    changes: List[SignatureExposureChange]
    total_signatures: int


PIPELINE_STEP_TYPES = ["spectrum", "nmf_extract", "scoring", "ld_analysis", "temporal_analysis"]
PIPELINE_STATUSES = ["pending", "running", "completed", "skipped", "failed"]


class PipelineStepInputParam(BaseModel):
    name: str = Field(..., min_length=1)
    value: Optional[Any] = None
    from_step: Optional[str] = Field(None, description="Field path expression, e.g. 'step_1.output.sample_ids'")

    @model_validator(mode="after")
    def check_value_or_from_step(self):
        if self.value is None and self.from_step is None:
            raise ValueError("Either 'value' or 'from_step' must be provided")
        if self.value is not None and self.from_step is not None:
            raise ValueError("Only one of 'value' or 'from_step' can be provided")
        return self


class PipelineStepCondition(BaseModel):
    expression: str = Field(..., description="Condition expression, e.g. 'step_1.output.sample_count > 5'")


class PipelineStepCreate(BaseModel):
    name: str = Field(..., min_length=1, max_length=100)
    step_type: str = Field(..., pattern="|".join(PIPELINE_STEP_TYPES))
    input_params: List[PipelineStepInputParam] = []
    condition: Optional[PipelineStepCondition] = None
    retry_count: int = Field(0, ge=0, le=10, description="Number of automatic retries on failure (default 0)")
    retry_interval_seconds: int = Field(0, ge=0, le=3600, description="Seconds to wait between retries (default 0)")


class PipelineStepOut(BaseModel):
    name: str
    step_type: str
    input_params: List[PipelineStepInputParam]
    condition: Optional[PipelineStepCondition] = None
    retry_count: int = 0
    retry_interval_seconds: int = 0


class PipelineTemplateCreate(BaseModel):
    name: str = Field(..., min_length=1, max_length=255)
    description: Optional[str] = ""
    steps: List[PipelineStepCreate] = Field(..., min_length=1)

    @model_validator(mode="after")
    def check_step_order(self):
        step_names = set()
        for i, step in enumerate(self.steps):
            if step.name in step_names:
                raise ValueError(f"Duplicate step name: {step.name}")
            step_names.add(step.name)
            for param in step.input_params:
                if param.from_step:
                    if param.from_step.startswith("initial."):
                        continue
                    if not param.from_step.startswith("step_"):
                        raise ValueError(
                            f"Invalid from_step format: {param.from_step}. "
                            f"Must start with 'initial.' or 'step_'"
                        )
                    try:
                        step_ref = param.from_step.split(".")[0]
                        step_num = int(step_ref.replace("step_", ""))
                        if step_num < 1 or step_num > i:
                            raise ValueError(
                                f"Invalid from_step reference: {param.from_step}. "
                                f"Can only reference steps 1 to {i} (current step is {i+1})"
                            )
                    except (ValueError, IndexError):
                        raise ValueError(f"Invalid from_step format: {param.from_step}")
            if step.condition:
                expr = step.condition.expression
                valid_prefixes = [f"step_{s}.summary." for s in range(1, i + 1)]
                valid_prefixes += [f"step_{s}.output." for s in range(1, i + 1)]
                if not any(expr.startswith(prefix) for prefix in valid_prefixes):
                    raise ValueError(
                        f"Invalid condition expression: must reference previous steps "
                        f"using 'step_X.summary.field' or 'step_X.output.field' (step_1 to step_{i})"
                    )
        return self


class PipelineTemplateUpdate(BaseModel):
    name: Optional[str] = Field(None, min_length=1, max_length=255)
    description: Optional[str] = None
    steps: Optional[List[PipelineStepCreate]] = Field(None, min_length=1)


class PipelineTemplateOut(BaseModel):
    id: int
    name: str
    description: str
    steps: List[PipelineStepOut]
    created_at: datetime
    updated_at: datetime

    class Config:
        from_attributes = True


class PipelineExecutionTrigger(BaseModel):
    initial_params: Dict[str, Any] = Field(..., description="Initial parameters for the pipeline")


class PipelineStepExecutionOut(BaseModel):
    step_index: int
    step_name: str
    step_type: str
    status: str
    input_params: Dict[str, Any]
    output_summary: Optional[Dict[str, Any]] = None
    output_data: Optional[Dict[str, Any]] = None
    error_message: Optional[str] = None
    started_at: Optional[datetime] = None
    completed_at: Optional[datetime] = None
    duration_seconds: Optional[float] = None
    retry_count: int = 0
    last_error: Optional[str] = None

    class Config:
        from_attributes = True


class PipelineExecutionOut(BaseModel):
    id: int
    template_id: int
    template_name: str
    status: str
    initial_params: Dict[str, Any]
    total_duration_seconds: Optional[float] = None
    error_message: Optional[str] = None
    triggered_at: datetime
    completed_at: Optional[datetime] = None
    step_executions: List[PipelineStepExecutionOut] = []
    template_snapshot: Optional[List[Dict[str, Any]]] = None
    resume_count: int = 0

    class Config:
        from_attributes = True


class PipelineExecutionSummaryOut(BaseModel):
    id: int
    template_id: int
    template_name: str
    status: str
    initial_params: Dict[str, Any]
    total_duration_seconds: Optional[float] = None
    error_message: Optional[str] = None
    triggered_at: datetime
    completed_at: Optional[datetime] = None
    resume_count: int = 0

    class Config:
        from_attributes = True


class PipelineResumeResponse(BaseModel):
    id: int
    template_id: int
    template_name: str
    status: str
    initial_params: Dict[str, Any]
    total_duration_seconds: Optional[float] = None
    error_message: Optional[str] = None
    triggered_at: datetime
    completed_at: Optional[datetime] = None
    step_executions: List[PipelineStepExecutionOut] = []
    resume_count: int = 0
    resumed_from_step: Optional[int] = None
    skipped_completed_steps: int = 0


VALID_QC_RULE_TYPES = ["alignment_score", "coverage_depth", "edge_position"]


class QCRuleCreateItem(BaseModel):
    rule_type: str = Field(..., pattern="|".join(VALID_QC_RULE_TYPES))
    min_alignment_score: Optional[int] = Field(None, description="Minimum alignment score threshold (for alignment_score rule)")
    min_coverage_depth: Optional[int] = Field(None, description="Minimum independent detection count (for coverage_depth rule)")
    edge_margin: Optional[int] = Field(None, description="Exclude variants within N bases of sequence ends (for edge_position rule)")

    @model_validator(mode="after")
    def validate_rule_params(self):
        if self.rule_type == "alignment_score" and self.min_alignment_score is None:
            raise ValueError("min_alignment_score is required for alignment_score rule")
        if self.rule_type == "coverage_depth" and self.min_coverage_depth is None:
            raise ValueError("min_coverage_depth is required for coverage_depth rule")
        if self.rule_type == "edge_position" and self.edge_margin is None:
            raise ValueError("edge_margin is required for edge_position rule")
        return self


class QCRuleSetCreate(BaseModel):
    name: str = Field(..., min_length=1, max_length=255)
    description: Optional[str] = ""
    rules: List[QCRuleCreateItem] = Field(..., min_length=1)
    activate: bool = False


class QCRuleSetUpdate(BaseModel):
    name: Optional[str] = Field(None, min_length=1, max_length=255)
    description: Optional[str] = None
    rules: Optional[List[QCRuleCreateItem]] = Field(None, min_length=1)


class QCRuleOut(BaseModel):
    id: int
    rule_index: int
    rule_type: str
    min_alignment_score: Optional[int] = None
    min_coverage_depth: Optional[int] = None
    edge_margin: Optional[int] = None

    class Config:
        from_attributes = True


class QCRuleSetOut(BaseModel):
    id: int
    name: str
    description: str
    is_active: bool
    rules: List[QCRuleOut] = []
    created_at: datetime
    updated_at: datetime

    class Config:
        from_attributes = True


class QCRuleSetListOut(BaseModel):
    id: int
    name: str
    description: str
    is_active: bool
    rule_count: int
    created_at: datetime
    updated_at: datetime


class FailedRuleDetail(BaseModel):
    rule_id: int
    rule_type: str
    rule_description: str
    actual_value: Any
    threshold: Any


class VariantQCOut(BaseModel):
    reference_id: int
    ref_pos: int
    variant_type: str
    ref_base: str
    alt_base: str
    qc_status: str
    failed_rules: List[FailedRuleDetail] = []
    is_stale: bool = False


class SampleQCResultOut(BaseModel):
    sample_id: int
    sample_name: str
    rule_set_id: int
    rule_set_name: str
    total_variants: int
    pass_count: int
    fail_count: int
    stale_count: int
    qc_results: List[VariantQCOut] = []
    executed_at: datetime


class SampleQCSummaryOut(BaseModel):
    sample_id: int
    sample_name: str
    rule_set_id: Optional[int] = None
    rule_set_name: Optional[str] = None
    has_qc_results: bool
    is_stale: bool
    total_variants: int
    pass_count: int
    fail_count: int
    stale_count: int


class SpectrumVariantQCOut(BaseModel):
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
    qc_status: Optional[str] = None
    qc_failed_rules: List[Dict[str, Any]] = []
    qc_is_stale: bool = False

    class Config:
        from_attributes = True


class SampleSpectrumQCOut(BaseModel):
    sample_id: int
    sample_name: str
    total_variants: int
    pass_count: Optional[int] = None
    fail_count: Optional[int] = None
    stale_count: Optional[int] = None
    active_rule_set: Optional[str] = None
    variants: List[SpectrumVariantQCOut] = []
