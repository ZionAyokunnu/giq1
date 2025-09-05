import pandas as pd
import numpy as np
from collections import defaultdict
import math

def extract_gene_distribution(markov_profile, busco_id, pseudo_count=0.000001):
    """
    Extract probability distribution for a specific gene across all bins in the profile.
    
    Args:
        markov_profile: Profile from build_markov_profile()
        busco_id: Gene ID to extract distribution for
        pseudo_count: Small value to add to all bins
        
    Returns:
        dict: {bin_id: normalized_probability}
    """
    gene_distribution = {}
    
    # Collect all bins and their probabilities for this gene
    for bin_id, genes_data in markov_profile.items():
        if busco_id in genes_data:
            # Use average_percentage as the probability
            probability = genes_data[busco_id]['average_percentage'] / 100.0
        else:
            # Gene doesn't appear in this bin
            probability = 0.0
        
        gene_distribution[bin_id] = probability
    
    # Add pseudo-counts to all bins
    for bin_id in gene_distribution:
        gene_distribution[bin_id] += pseudo_count
    
    # Normalize so probabilities sum to 1.0
    total_prob = sum(gene_distribution.values())
    for bin_id in gene_distribution:
        gene_distribution[bin_id] /= total_prob
    
    return gene_distribution


def calculate_distribution_stats(gene_distribution):
    """
    Calculate statistics from a gene's probability distribution.
    
    Args:
        gene_distribution: {bin_id: normalized_probability}
        
    Returns:
        dict: Statistics including mean, std, 95% CI
    """
    # Convert bin_ids to numerical positions for calculation
    # Assuming bin_id format like "chr1_bin_5"
    positions = []
    probabilities = []
    
    for bin_id, prob in gene_distribution.items():
        # Extract bin number from bin_id
        try:
            bin_number = int(bin_id.split('_bin_')[1])
            positions.append(bin_number)
            probabilities.append(prob)
        except (IndexError, ValueError):
            # Skip bins that don't match expected format
            continue
    
    if len(positions) == 0:
        return None
    
    positions = np.array(positions)
    probabilities = np.array(probabilities)
    
    # Calculate weighted statistics
    mean_position = np.sum(positions * probabilities)
    variance = np.sum(probabilities * (positions - mean_position) ** 2)
    std_dev = np.sqrt(variance)
    
    # Calculate empirical 95% CI using percentiles
    # Create cumulative distribution
    sorted_indices = np.argsort(positions)
    sorted_positions = positions[sorted_indices]
    sorted_probs = probabilities[sorted_indices]
    cumulative_probs = np.cumsum(sorted_probs)
    
    # Find 2.5th and 97.5th percentiles
    ci_lower = np.interp(0.025, cumulative_probs, sorted_positions)
    ci_upper = np.interp(0.975, cumulative_probs, sorted_positions)
    
    return {
        'mean_position': mean_position,
        'std_deviation': std_dev,
        'variance': variance,
        'ci_95_lower': ci_lower,
        'ci_95_upper': ci_upper,
        'total_bins': len(positions)
    }


def calculate_bit_scores(gene_distribution, query_bin_id, total_genes):
    """
    Calculate bit scores for a gene using average approach.
    
    Args:
        gene_distribution: {bin_id: normalized_probability}
        query_bin_id: Bin where the query gene actually appears
        total_genes: Total number of genes (for E-value calculation)
        
    Returns:
        dict: Bit scores and E-value
    """
    # Position-specific bit score (for the bin where gene appears)
    position_specific_prob = gene_distribution.get(query_bin_id, 0.0)
    position_bit_score = -math.log2(position_specific_prob) if position_specific_prob > 0 else float('inf')
    
    # Overall match bit score (average across all bins)
    # Only include bins where gene actually appears (prob > pseudo_count threshold)
    meaningful_probs = [prob for prob in gene_distribution.values() if prob > 1e-5]  # Above pseudo-count
    
    if meaningful_probs:
        # Average (overall) bit score across bins where gene appears
        individual_bit_scores = [-math.log2(prob) for prob in meaningful_probs]
        overall_bit_score = sum(individual_bit_scores) / len(individual_bit_scores)  # Average
        overall_prob = 2**(-overall_bit_score)  # Convert back to probability
    else:
        overall_bit_score = float('inf')
        overall_prob = 0.0
    
    # E-value calculation (unchanged)
    e_value = total_genes * position_specific_prob if position_specific_prob > 0 else float('inf')
    
    return {
        'position_specific_bit_score': position_bit_score,
        'overall_match_bit_score': overall_bit_score,  # Now finite!
        'position_specific_probability': position_specific_prob,
        'overall_probability': overall_prob,
        'e_value': e_value
    }
    
    

def compare_query_genome_to_profile(query_bin_assignments, markov_profile):
    """
    Compare a query genome to the Markov profile.
    Now handles chromosome-grouped query data.
    
    Args:
        query_bin_assignments: {chromosome: {busco_id: [(bin_id, overlap_percentage), ...]}} for query genome
        markov_profile: Profile from build_markov_profile()
        
    Returns:
        dict: {chromosome: {busco_id: comparison_results}}
    """
    # Count total unique genes in profile for E-value calculation
    all_profile_genes = set()
    for genes_data in markov_profile.values():
        all_profile_genes.update(genes_data.keys())
    total_genes = len(all_profile_genes)
    
    comparison_results = {}
    
    # Process each chromosome separately
    for chromosome, gene_bin_assignments in query_bin_assignments.items():
        comparison_results[chromosome] = {}
        
        for busco_id, bin_overlaps in gene_bin_assignments.items():
            # For simplicity, use the bin with highest overlap for this gene
            if not bin_overlaps:
                continue
                
            primary_bin = max(bin_overlaps, key=lambda x: x[1])  # (bin_id, percentage)
            query_bin_id = primary_bin[0]
            
            # Extract gene distribution from profile
            gene_distribution = extract_gene_distribution(markov_profile, busco_id)
            
            # Calculate distribution statistics
            dist_stats = calculate_distribution_stats(gene_distribution)
            
            # Calculate deviation from expected position
            deviation = None
            std_deviations = None
            if dist_stats:
                # Extract bin number from query_bin_id
                try:
                    query_bin_number = int(query_bin_id.split('_bin_')[1])
                    deviation = abs(query_bin_number - dist_stats['mean_position'])
                    std_deviations = deviation / dist_stats['std_deviation'] if dist_stats['std_deviation'] > 0 else 0
                except (IndexError, ValueError):
                    pass
            
            # Calculate bit scores
            bit_scores = calculate_bit_scores(gene_distribution, query_bin_id, total_genes)
            
            # Compile results
            comparison_results[chromosome][busco_id] = {
                'query_bin': query_bin_id,
                'query_overlap_percentage': primary_bin[1],
                'expected_position': dist_stats['mean_position'] if dist_stats else None,
                'position_deviation': deviation,
                'standard_deviations': std_deviations,
                'distribution_stats': dist_stats,
                'bit_scores': bit_scores,
                'gene_distribution': gene_distribution
            }
    
    return comparison_results
