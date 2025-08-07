import pandas as pd
import numpy as np
from collections import defaultdict


from config.settings import CONFIG

def extract_gene_ranges(gene_distribution, probability_threshold=None):
    """
    Extract bins where gene has significant probability (above threshold).
    
    Args:
        gene_distribution: {bin_id: normalised_probability is commented out for no normalisation}
        probability_threshold: Minimum probability to consider (uses CONFIG if None)
        
    Returns:
        list: [(bin_number, probability)] sorted by probability (descending)
    """
    if probability_threshold is None:
        probability_threshold = CONFIG['probability_threshold_for_target']
    
    significant_bins = []
    
    for bin_id, probability in gene_distribution.items():
        if probability >= probability_threshold:
            try:
                chromosome, bin_part = bin_id.split('_bin_')
                bin_number = int(bin_part)
                significant_bins.append((bin_id, bin_number, probability))
            except (IndexError, ValueError):
                continue
    
    significant_bins.sort(key=lambda x: (-x[1], x[0]))
    
    return significant_bins



def extract_gene_ranges(gene_distribution, probability_threshold=None):
    """
    Extract bins where gene has significant probability (above threshold).
    
    Args:
        gene_distribution: {bin_id: normalised_probability}
        probability_threshold: Minimum probability to consider (uses CONFIG if None)
        
    Returns:
        list: [(bin_id, bin_number, probability)] sorted by probability (descending)
    """
    if probability_threshold is None:
        probability_threshold = CONFIG['probability_threshold_for_target']
    
    significant_bins = []
    
    for bin_id, probability in gene_distribution.items():
        if probability >= probability_threshold:
            try:
                chromosome, bin_part = bin_id.split('_bin_')
                bin_number = int(bin_part)
                significant_bins.append((bin_id, bin_number, probability))
            except (IndexError, ValueError):
                continue
    
    # Sort by probability (descending), then by bin_number
    significant_bins.sort(key=lambda x: (-x[2], x[1]))
    
    return significant_bins


def extract_current_ranges(bin_overlaps):
    """
    Extract current bin ranges from query genome bin assignments.
    
    Args:
        bin_overlaps: [(bin_id, overlap_percentage), ...]
        
    Returns:
        list: [(bin_id, bin_number, overlap_percentage)] sorted by percentage (descending)
    """
    current_bins = []
    
    for bin_id, overlap_percentage in bin_overlaps:
        try:
            chromosome, bin_part = bin_id.split('_bin_')
            bin_number = int(bin_part)
            current_bins.append((bin_id, bin_number, overlap_percentage))
        except (IndexError, ValueError):
            continue
    
    # Sort by overlap percentage (descending), then by bin_number
    current_bins.sort(key=lambda x: (-x[2], x[1]))
    
    return current_bins


def filter_same_chromosome(current_ranges, target_ranges):
    """
    Filter target ranges to only include bins from same chromosomes as current ranges.
    
    Args:
        current_ranges: [(bin_id, bin_number, overlap_percentage)]
        target_ranges: [(bin_id, bin_number, probability)]
        
    Returns:
        tuple: (filtered_target_ranges, interchromosomal_targets)
    """
    # Extract chromosomes from current ranges
    current_chromosomes = set()
    for bin_id, _, _ in current_ranges:
        chromosome = bin_id.split('_bin_')[0]
        current_chromosomes.add(chromosome)
    
    # Filter target ranges
    filtered_targets = []
    interchromosomal_targets = []
    
    for bin_id, bin_number, probability in target_ranges:
        chromosome = bin_id.split('_bin_')[0]
        if chromosome in current_chromosomes:
            filtered_targets.append((bin_id, bin_number, probability))
        else:
            interchromosomal_targets.append((bin_id, bin_number, probability))
    
    return filtered_targets, interchromosomal_targets


def calculate_gene_movement(current_ranges, target_ranges):
    """
    Calculate movement for a gene with chromosome awareness.
    
    Args:
        current_ranges: [(bin_id, bin_number, overlap_percentage)] from query genome
        target_ranges: [(bin_id, bin_number, probability)] from profile
        
    Returns:
        dict: Movement analysis results with structural variation flags
    """
    if not current_ranges or not target_ranges:
        return {
            'movements': [],
            'mean_movement': 0.0,
            'pairings': [],
            'unpaired_current': current_ranges.copy(),
            'unpaired_target': target_ranges.copy(),
            'total_pairs': 0,
            'structural_variation': None,
            'interchromosomal_targets': []
        }
    
    # Filter for same-chromosome targets
    filtered_targets, interchromosomal_targets = filter_same_chromosome(current_ranges, target_ranges)
    
    # Check for structural variations
    structural_variation = None
    if not filtered_targets and interchromosomal_targets:
        structural_variation = "translocation"
    elif filtered_targets and interchromosomal_targets:
        structural_variation = "mixed_synteny"
    
    # If no same-chromosome targets, flag as translocation
    if not filtered_targets:
        return {
            'movements': [],
            'mean_movement': 0.0,
            'pairings': [],
            'unpaired_current': current_ranges,
            'unpaired_target': target_ranges,
            'total_pairs': 0,
            'structural_variation': structural_variation,
            'interchromosomal_targets': interchromosomal_targets
        }
    
    # Proceed with movement calculation using filtered targets
    pairings = []
    movements = []
    remaining_current = current_ranges.copy()
    remaining_target = filtered_targets.copy()
    
    # First pairing: Highest probability current ↔ Highest probability target
    if remaining_current and remaining_target:
        current_bin_id, current_bin_num, current_prob = remaining_current.pop(0)
        target_bin_id, target_bin_num, target_prob = remaining_target.pop(0)
        
        movement = target_bin_num - current_bin_num
        pairings.append({
            'current_bin_id': current_bin_id,
            'current_bin_num': current_bin_num,
            'target_bin_id': target_bin_id,
            'target_bin_num': target_bin_num,
            'current_prob': current_prob,
            'target_prob': target_prob,
            'movement': movement,
            'pairing_type': 'highest_probability'
        })
        movements.append(movement)
    
    # Subsequent pairings: Minimize distance with adjacency constraint
    while remaining_current and remaining_target:
        current_bin_id, current_bin_num, current_prob = remaining_current.pop(0)
        
        # Find target bin that minimizes distance and maintains adjacency
        paired_targets = [p['target_bin_num'] for p in pairings]
        
        best_target = None
        best_distance = float('inf')
        best_index = -1
        
        for i, (target_bin_id, target_bin_num, target_prob) in enumerate(remaining_target):
            # Check adjacency to already paired targets
            is_adjacent = any(abs(target_bin_num - paired_target) == 1 for paired_target in paired_targets)
            
            # If no adjacency, check gene wholeness constraint
            if not is_adjacent and len(paired_targets) > 0:
                min_paired = min(paired_targets)
                max_paired = max(paired_targets)
                
                # Allow some flexibility for gene wholeness (±2 range)
                if not (min_paired - 2 <= target_bin_num <= max_paired + 2):
                    continue
            
            distance = abs(target_bin_num - current_bin_num)
            if distance < best_distance:
                best_distance = distance
                best_target = (target_bin_id, target_bin_num, target_prob)
                best_index = i
        
        # If no adjacency-constrained target found, find closest remaining target
        if best_target is None and remaining_target:
            for i, (target_bin_id, target_bin_num, target_prob) in enumerate(remaining_target):
                distance = abs(target_bin_num - current_bin_num)
                if distance < best_distance:
                    best_distance = distance
                    best_target = (target_bin_id, target_bin_num, target_prob)
                    best_index = i
        
        if best_target is not None:
            target_bin_id, target_bin_num, target_prob = remaining_target.pop(best_index)
            movement = target_bin_num - current_bin_num
            
            pairings.append({
                'current_bin_id': current_bin_id,
                'current_bin_num': current_bin_num,
                'target_bin_id': target_bin_id,
                'target_bin_num': target_bin_num,
                'current_prob': current_prob,
                'target_prob': target_prob,
                'movement': movement,
                'pairing_type': 'minimize_distance'
            })
            movements.append(movement)
    
    # Calculate mean movement
    mean_movement = np.mean(movements) if movements else 0.0
    
    return {
        'movements': movements,
        'mean_movement': mean_movement,
        'pairings': pairings,
        'unpaired_current': remaining_current,
        'unpaired_target': remaining_target,
        'total_pairs': len(pairings),
        'structural_variation': structural_variation,
        'interchromosomal_targets': interchromosomal_targets
    }


def analyse_query_movements(query_bin_assignments, markov_profile):
    """        
    Analyse movements per chromosome with structural variation detection.
    
    Args:
        query_bin_assignments: {chromosome: {busco_id: [(bin_id, overlap_percentage)]}}
        markov_profile: Profile from build_markov_profile()
        
    Returns:
        dict: {chromosome: {busco_id: movement_analysis}} with structural variation flags
    """
    movement_results = {}
    structural_variations = {
        'translocations': [],
        'mixed_synteny': [],
        'total_genes_analyzed': 0,
        'intrachromosomal_genes': 0
    }
    
    # Process each chromosome separately
    for chromosome, gene_bin_assignments in query_bin_assignments.items():
        movement_results[chromosome] = {}
        
        for busco_id, bin_overlaps in gene_bin_assignments.items():
            structural_variations['total_genes_analyzed'] += 1
            
            # Extract current ranges (with chromosome info)
            current_ranges = extract_current_ranges(bin_overlaps)
            
            # Extract gene distribution from profile
            gene_distribution = extract_gene_distribution(markov_profile, busco_id)
            
            # Extract target ranges (with chromosome info)
            target_ranges = extract_gene_ranges(gene_distribution)
            
            # Calculate movement with chromosome awareness
            movement_analysis = calculate_gene_movement(current_ranges, target_ranges)
            
            # Track structural variations
            sv_type = movement_analysis['structural_variation']
            if sv_type == 'translocation':
                structural_variations['translocations'].append({
                    'chromosome': chromosome,
                    'busco_id': busco_id,
                    'current_ranges': current_ranges,
                    'interchromosomal_targets': movement_analysis['interchromosomal_targets']
                })
            elif sv_type == 'mixed_synteny':
                structural_variations['mixed_synteny'].append({
                    'chromosome': chromosome,
                    'busco_id': busco_id,
                    'current_ranges': current_ranges,
                    'interchromosomal_targets': movement_analysis['interchromosomal_targets']
                })
            else:
                structural_variations['intrachromosomal_genes'] += 1
            
            movement_results[chromosome][busco_id] = {
                'current_ranges': current_ranges,
                'target_ranges': target_ranges,
                'movement_analysis': movement_analysis,
                'gene_distribution': gene_distribution
            }
    
    return movement_results, structural_variations





def extract_gene_distribution(markov_profile, busco_id, pseudo_count=0.000001):
    """
    Extract and normalise gene distribution (commented out) for movement analysis.
    """
    gene_distribution = {}
    
    for bin_id, genes_data in markov_profile.items():
        if busco_id in genes_data:
            probability = genes_data[busco_id]['average_percentage'] / 100.0
        else:
            probability = 0.0
        gene_distribution[bin_id] = probability
    
    for bin_id in gene_distribution:
        gene_distribution[bin_id] += pseudo_count
    
    # total_prob = sum(gene_distribution.values())
    # for bin_id in gene_distribution:
    #     gene_distribution[bin_id] /= total_prob
    
    return gene_distribution


def get_movement_summary(movement_results):
    """
    Get summary statistics of movement analysis.
    Now handles chromosome structure.
    """
    total_genes = 0
    genes_with_movement = 0
    all_movements = []
    
    # Handle chromosome structure: {chromosome: {busco_id: result}}
    for chromosome, gene_results in movement_results.items():
        for busco_id, result in gene_results.items():
            total_genes += 1
            
            if result['movement_analysis']['mean_movement'] != 0:
                genes_with_movement += 1
            
            if result['movement_analysis']['movements']:
                all_movements.extend(result['movement_analysis']['movements'])
    
    summary = {
        'total_genes': total_genes,
        'genes_with_movement': genes_with_movement,
        'genes_no_movement': total_genes - genes_with_movement,
        'total_movements': len(all_movements),
        'mean_absolute_movement': np.mean([abs(m) for m in all_movements]) if all_movements else 0,
        'movement_distribution': {
            'left_movements': sum(1 for m in all_movements if m < 0),
            'right_movements': sum(1 for m in all_movements if m > 0),
            'no_movement': sum(1 for m in all_movements if m == 0)
        }
    }
    
    return summary
    



if __name__ == "__main__":
    print("working")