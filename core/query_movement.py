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
                bin_number = int(bin_id.split('_bin_')[1])
                significant_bins.append((bin_number, probability))
            except (IndexError, ValueError):
                continue
    

    significant_bins.sort(key=lambda x: (-x[1], x[0]))
    
    return significant_bins







# def extract_gene_ranges(markov_profile, busco_id, probability_threshold=0.7):
#     """Extract bins using RAW profile percentages, not normalized."""
#     significant_bins = []
    
#     for bin_id, genes_data in markov_profile.items():
#         if busco_id in genes_data:
#             # Use RAW percentage from profile
#             raw_percentage = genes_data[busco_id]['average_percentage']
#             if raw_percentage >= probability_threshold * 100:
#                 bin_number = int(bin_id.split('_bin_')[1])
#                 significant_bins.append((bin_number, raw_percentage/100.0))
    
#     return sorted(significant_bins, key=lambda x: (-x[1], x[0]))


def extract_current_ranges(bin_overlaps):
    """
    Extract current bin ranges from query genome bin assignments.
    
    Args:
        bin_overlaps: [(bin_id, overlap_percentage), ...]
        
    Returns:
        list: [(bin_number, overlap_percentage)] sorted by percentage (descending)
    """
    current_bins = []
    
    for bin_id, overlap_percentage in bin_overlaps:
        try:
            bin_number = int(bin_id.split('_bin_')[1])
            current_bins.append((bin_number, overlap_percentage))
        except (IndexError, ValueError):
            continue
    
    current_bins.sort(key=lambda x: (-x[1], x[0]))
    
    return current_bins


def calculate_gene_movement(current_ranges, target_ranges):
    """
    Calculate movement for a gene using the pairing algorithm.
    
    Args:
        current_ranges: [(bin_number, overlap_percentage)] from query genome
        target_ranges: [(bin_number, probability)] from profile
        
    Returns:
        dict: Movement analysis results
    """
    if not current_ranges or not target_ranges:
        return {
            'movements': [],
            'mean_movement': 0.0,
            'pairings': [],
            'unpaired_current': current_ranges.copy(),
            'unpaired_target': target_ranges.copy(),
            'total_pairs': 0
        }
    
    pairings = []
    movements = []
    remaining_current = current_ranges.copy()
    remaining_target = target_ranges.copy()
    
    # First pairing: Highest probability current ---- Highest probability target
    if remaining_current and remaining_target:
        current_bin, current_prob = remaining_current.pop(0)
        target_bin, target_prob = remaining_target.pop(0)
        
        movement = target_bin - current_bin
        pairings.append({
            'current_bin': current_bin,
            'target_bin': target_bin,
            'current_prob': current_prob,
            'target_prob': target_prob,
            'movement': movement,
            'pairing_type': 'highest_probability'
        })
        movements.append(movement)
    
    # Subsequent pairings
    while remaining_current and remaining_target:
        current_bin, current_prob = remaining_current.pop(0)  # Next highest current
        
        # Find target bin that minimizes distance and is adjacent to already paired targets
        paired_targets = [p['target_bin'] for p in pairings]
        
        best_target = None
        best_distance = float('inf')
        best_index = -1
        
        for i, (target_bin, target_prob) in enumerate(remaining_target):
  
            is_adjacent = any(abs(target_bin - paired_target) == 1 for paired_target in paired_targets)
            
            # If no adjacent constraint can be satisfied, allow any pairing ---????
            if not is_adjacent and len(paired_targets) > 0:
                #gene wholeness is important.
                min_paired = min(paired_targets)
                max_paired = max(paired_targets)
                
                if not (min_paired - 2 <= target_bin <= max_paired + 2):
                    continue
            
            distance = abs(target_bin - current_bin)
            if distance < best_distance:
                best_distance = distance
                best_target = (target_bin, target_prob)
                best_index = i
        

        if best_target is None and remaining_target:
            for i, (target_bin, target_prob) in enumerate(remaining_target):
                distance = abs(target_bin - current_bin)
                if distance < best_distance:
                    best_distance = distance
                    best_target = (target_bin, target_prob)
                    best_index = i
        
        if best_target is not None:
            target_bin, target_prob = remaining_target.pop(best_index)
            movement = target_bin - current_bin
            
            pairings.append({
                'current_bin': current_bin,
                'target_bin': target_bin,
                'current_prob': current_prob,
                'target_prob': target_prob,
                'movement': movement,
                'pairing_type': 'minimize_distance'
            })
            movements.append(movement)
    
    #Mean movement
    mean_movement = np.mean(movements) if movements else 0.0
    
    return {
        'movements': movements,
        'mean_movement': mean_movement,
        'pairings': pairings,
        'unpaired_current': remaining_current,
        'unpaired_target': remaining_target,
        'total_pairs': len(pairings)
    }


def analyse_query_movements(query_bin_assignments, markov_profile):
    """        
    Analyse movements and returns:
        dict: {busco_id: movement_analysis}
    """
    movement_results = {}
    
    for busco_id, bin_overlaps in query_bin_assignments.items():

        current_ranges = extract_current_ranges(bin_overlaps)
        
        gene_distribution = extract_gene_distribution(markov_profile, busco_id)
        
        target_ranges = extract_gene_ranges(gene_distribution)
        
        movement_analysis = calculate_gene_movement(current_ranges, target_ranges)
        
        movement_results[busco_id] = {
            'current_ranges': current_ranges,
            'target_ranges': target_ranges,
            'movement_analysis': movement_analysis,
            'gene_distribution': gene_distribution
        }
    
    return movement_results


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

    total_genes = len(movement_results)
    genes_with_movement = sum(1 for r in movement_results.values() 
                             if r['movement_analysis']['mean_movement'] != 0)
    
    all_movements = []
    for result in movement_results.values():
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