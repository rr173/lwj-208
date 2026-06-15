from typing import List, Dict, Tuple, Optional
from datetime import datetime
import math
from .tree_building import TreeNode


def compute_root_to_tip_distances(root: TreeNode) -> Dict[str, float]:
    """
    Compute root-to-tip distances for all leaf nodes.
    
    Returns:
        {leaf_name: distance_from_root}
    """
    distances = {}
    
    def traverse(node: TreeNode, current_dist: float):
        if node.is_leaf:
            if node.name:
                distances[node.name] = current_dist
            return
        
        for child in node.children:
            bl = child.branch_length if child.branch_length is not None else 0.0
            traverse(child, current_dist + bl)
    
    traverse(root, 0.0)
    return distances


def linear_regression(
    x: List[float], 
    y: List[float]
) -> Tuple[float, float, float, List[float]]:
    """
    Perform simple linear regression: y = slope * x + intercept
    
    Returns:
        (slope, intercept, r_squared, residuals)
    """
    n = len(x)
    if n < 2:
        return 0.0, 0.0, 0.0, [0.0] * n
    
    mean_x = sum(x) / n
    mean_y = sum(y) / n
    
    ss_xy = sum((xi - mean_x) * (yi - mean_y) for xi, yi in zip(x, y))
    ss_xx = sum((xi - mean_x) ** 2 for xi in x)
    
    if ss_xx == 0:
        return 0.0, mean_y, 0.0, [yi - mean_y for yi in y]
    
    slope = ss_xy / ss_xx
    intercept = mean_y - slope * mean_x
    
    predicted = [slope * xi + intercept for xi in x]
    residuals = [yi - pred for yi, pred in zip(y, predicted)]
    
    ss_total = sum((yi - mean_y) ** 2 for yi in y)
    ss_residual = sum(r ** 2 for r in residuals)
    
    r_squared = 1 - (ss_residual / ss_total) if ss_total > 0 else 0.0
    
    return slope, intercept, r_squared, residuals


def date_to_year(date: datetime) -> float:
    """Convert datetime to decimal year."""
    year = date.year
    if hasattr(date, 'month'):
        year += (date.month - 1) / 12.0
    if hasattr(date, 'day'):
        year += (date.day - 1) / 365.25
    return year


def fit_molecular_clock(
    root: TreeNode,
    sample_dates: Dict[str, datetime],
    sample_names: List[str]
) -> Tuple[Optional[Dict], Optional[str]]:
    """
    Fit a molecular clock using root-to-tip distances vs collection dates.
    
    Args:
        root: The phylogenetic tree root
        sample_dates: {sample_name: collection_date}
        sample_names: List of all sample names in the tree
    
    Returns:
        (clock_info_dict, warning_message)
        clock_info contains: slope, intercept, r_squared, rate, outlier_samples
    """
    root_to_tip = compute_root_to_tip_distances(root)
    
    dated_samples = []
    for name in sample_names:
        if name in sample_dates and name in root_to_tip:
            dated_samples.append((name, sample_dates[name], root_to_tip[name]))
    
    if len(dated_samples) < 3:
        return None, f"Only {len(dated_samples)} samples have collection dates. At least 3 required for molecular clock."
    
    years = []
    distances = []
    names_in_regression = []
    
    for name, date, dist in dated_samples:
        year = date_to_year(date)
        years.append(year)
        distances.append(dist)
        names_in_regression.append(name)
    
    slope, intercept, r_squared, residuals = linear_regression(years, distances)
    
    outlier_indices = []
    if len(residuals) > 0:
        mean_residual = sum(abs(r) for r in residuals) / len(residuals)
        std_residual = math.sqrt(sum((r - mean_residual) ** 2 for r in residuals) / len(residuals))
        
        for i, r in enumerate(residuals):
            if abs(r) > 2 * (std_residual + mean_residual):
                outlier_indices.append(i)
    
    outlier_samples = [names_in_regression[i] for i in outlier_indices]
    
    warning = None
    if r_squared <= 0.5:
        warning = f"R² = {r_squared:.4f} <= 0.5, molecular clock fit is poor. Divergence time estimates may be unreliable."
    
    rate = abs(slope) if slope != 0 else 0.0
    
    clock_info = {
        "slope": slope,
        "intercept": intercept,
        "r_squared": r_squared,
        "rate": rate,
        "outlier_samples": outlier_samples,
        "warning": warning,
    }
    
    return clock_info, warning


def estimate_divergence_times(
    root: TreeNode,
    slope: float,
    intercept: float,
    sample_dates: Dict[str, datetime]
) -> TreeNode:
    """
    Estimate divergence times for all internal nodes using the molecular clock.
    
    Time is calculated as:
        For a node at distance d from root, the divergence time in years is:
            t = (d - intercept) / slope
    
    The tree is modified in-place with divergence_time set on each node.
    
    Args:
        root: The phylogenetic tree
        slope: Slope from root-to-tip regression (distance = slope * year + intercept)
        intercept: Intercept from regression
        sample_dates: {sample_name: collection_date} for known dates
    
    Returns:
        Tree with divergence_time set on all nodes
    """
    if slope == 0:
        return root
    
    def traverse(node: TreeNode, dist_from_root: float):
        if slope != 0:
            year = (dist_from_root - intercept) / slope
            node.divergence_time = year
        
        if node.is_leaf and node.name and node.name in sample_dates:
            known_year = date_to_year(sample_dates[node.name])
            node.divergence_time = known_year
        
        for child in node.children:
            bl = child.branch_length if child.branch_length is not None else 0.0
            traverse(child, dist_from_root + bl)
    
    traverse(root, 0.0)
    return root


def apply_molecular_clock(
    root: TreeNode,
    sample_dates: Dict[str, datetime],
    sample_names: List[str]
) -> Tuple[TreeNode, Optional[Dict], Optional[str]]:
    """
    Apply molecular clock analysis: fit clock and estimate divergence times.
    
    Args:
        root: Phylogenetic tree root
        sample_dates: {sample_name: collection_date}
        sample_names: List of all sample names
    
    Returns:
        (tree_with_divergence_times, clock_info, warning)
    """
    clock_info, warning = fit_molecular_clock(root, sample_dates, sample_names)
    
    if clock_info and clock_info["r_squared"] > 0.5:
        root = estimate_divergence_times(
            root,
            clock_info["slope"],
            clock_info["intercept"],
            sample_dates
        )
    
    return root, clock_info, warning
