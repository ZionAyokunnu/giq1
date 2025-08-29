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
    Extract bins where gene has significant probability (updated for hybrid support).
    """
    if probability_threshold is None:
        probability_threshold = CONFIG['probability_threshold_for_target']
    
    significant_bins = []
    
    for bin_id, probability in gene_distribution.items():
        if probability >= probability_threshold:
            try:
                if '_rank_' in bin_id:
                    # Ordinal data: extract rank number
                    chromosome, rank_part = bin_id.split('_rank_')
                    rank_number = int(rank_part)
                    significant_bins.append((bin_id, rank_number, probability))
                else:
                    # Positional data: extract bin number
                    chromosome, bin_part = bin_id.split('_bin_')
                    bin_number = int(bin_part)
                    significant_bins.append((bin_id, bin_number, probability))
            except (IndexError, ValueError):
                continue
    
    # Sort by probability (descending), then by number
    significant_bins.sort(key=lambda x: (-x[2], x[1]))
    
    return significant_bins


def extract_current_ranges(bin_overlaps):
    """
    Extract current ranges from hybrid bin assignments (positional + ordinal).
    """
    current_bins = []
    
    # Debug: Check the structure of bin_overlaps
    print(f"DEBUG - bin_overlaps type: {type(bin_overlaps)}")
    print(f"DEBUG - bin_overlaps: {bin_overlaps}")
    
    for item in bin_overlaps:
        print(f"DEBUG - Processing item: {item}, type: {type(item)}")
        
        # Check if item is a tuple/list with exactly 2 elements
        if isinstance(item, (tuple, list)) and len(item) == 2:
            bin_id, overlap_percentage = item
        else:
            print(f"DEBUG - Skipping malformed item: {item}")
            continue
            
        try:
            if '_rank_' in bin_id:
                # Ordinal data
                chromosome, rank_part = bin_id.split('_rank_')
                rank_number = int(rank_part)
                current_bins.append((bin_id, rank_number, overlap_percentage))
            else:
                # Positional data
                chromosome, bin_part = bin_id.split('_bin_')
                bin_number = int(bin_part)
                current_bins.append((bin_id, bin_number, overlap_percentage))
        except (IndexError, ValueError) as e:
            print(f"DEBUG - Error processing {bin_id}: {e}")
            continue
    
    # Sort by overlap percentage (descending), then by number
    current_bins.sort(key=lambda x: (-x[2], x[1]))
    
    return current_bins


def filter_same_chromosome(current_ranges, target_ranges):
    """
    Filter target ranges to include bins from same chromosomes (updated for hybrid).
    """
    # Extract chromosomes from current ranges
    current_chromosomes = set()
    for bin_id, _, _ in current_ranges:
        if '_rank_' in bin_id:
            chromosome = bin_id.split('_rank_')[0]
        else:
            chromosome = bin_id.split('_bin_')[0]
        current_chromosomes.add(chromosome)
    
    # Filter target ranges
    filtered_targets = []
    interchromosomal_targets = []
    
    for bin_id, bin_number, probability in target_ranges:
        if '_rank_' in bin_id:
            chromosome = bin_id.split('_rank_')[0]
        else:
            chromosome = bin_id.split('_bin_')[0]
            
        if chromosome in current_chromosomes:
            filtered_targets.append((bin_id, bin_number, probability))
        else:
            interchromosomal_targets.append((bin_id, bin_number, probability))

    print(f"Current chromosomes: {current_chromosomes}")
    print(f"Target chromosome samples: {[bin_id.split('_rank_' if '_rank_' in bin_id else '_bin_')[0] for bin_id, _, _ in target_ranges[:5]]}")
    print(f"Filtered targets count: {len(filtered_targets)}")

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


def analyse_query_movements(query_bin_assignments, markov_profile, input_format='hybrid'):
    """        
    Analyse movements per chromosome with structural variation detection.
    
    Args:
        query_bin_assignments: {chromosome: {busco_id: [(bin_id, overlap_percentage)]}}
        markov_profile: Profile from build_markov_profile()
        
    Returns:
        dict: {chromosome: {busco_id: movement_analysis}} with structural variation flags
    """
    
    # Convert hybrid format to legacy format if needed
    if input_format == 'hybrid':
        converted_assignments = convert_hybrid_to_legacy_format(query_bin_assignments)
    else:
        converted_assignments = query_bin_assignments
    
    movement_results = {}
    structural_variations = {
        'translocations': [],
        'mixed_synteny': [],
        'total_genes_analyzed': 0,
        'intrachromosomal_genes': 0
    }
    
    # Process each chromosome separately
    for chromosome, gene_bin_assignments in converted_assignments.items():
        movement_results[chromosome] = {}
        
        for busco_id, bin_overlaps in gene_bin_assignments.items():
            structural_variations['total_genes_analyzed'] += 1
            
            # Extract current ranges (with chromosome info)
            current_ranges = extract_current_ranges(bin_overlaps)
            
            # Debug
            print(f"DEBUG - current_ranges for {busco_id}:")
            for bin_id, bin_num, overlap in current_ranges:
                print(f"  {bin_id}: bin_{bin_num}, overlap={overlap}")
    
            # Extract gene distribution from profile
            gene_distribution = extract_gene_distribution_hybrid_aware(markov_profile, busco_id)
            
            # Debug
            print(f"DEBUG - gene_distribution for {busco_id}:")
            for bin_id, prob in gene_distribution.items():
                if prob > 0.01:  # Only show significant probabilities
                    print(f"  {bin_id}: {prob}")
        
            # Extract target ranges (with chromosome info)
            target_ranges = extract_gene_ranges(gene_distribution)
            
            # Debug
            print(f"DEBUG - target_ranges for {busco_id}:")
            for bin_id, bin_num, prob in target_ranges:
                if '_rank_' in bin_id:
                    target_chromosome = bin_id.split('_rank_')[0]  # ← Use different variable name
                else:
                    target_chromosome = bin_id.split('_bin_')[0]   # ← Use different variable name
                print(f"  {bin_id}: {prob}")
    
            # Calculate movement with chromosome awareness
            try:
                movement_analysis = calculate_gene_movement(current_ranges, target_ranges)
            except Exception as e:
                print(f"ERROR in movement calculation: {e}")
                print(f"current_ranges: {current_ranges}")
                print(f"target_ranges: {target_ranges}")
                raise
            
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
            
            # # Ensure the chromosome key exists in movement_results
            # if chromosome not in movement_results:
            #     movement_results[chromosome] = {}
            
            movement_results[chromosome][busco_id] = {
                'current_ranges': current_ranges,
                'target_ranges': target_ranges,
                'movement_analysis': movement_analysis,
                'gene_distribution': gene_distribution
            }
    
    return movement_results, structural_variations
    
def convert_hybrid_to_legacy_format(hybrid_assignments):
    """
    Convert hybrid format to legacy format, combining positional and ordinal data.
    
    Args:
        hybrid_assignments: {genome_id: {chromosome: {busco_id: hybrid_data}}}
        
    Returns:
        dict: {chromosome: {busco_id: [(bin_id, overlap_percentage)]}}
    """
    converted = {}
    
    # Extract the single genome's data (query genome)
    genome_data = list(hybrid_assignments.values())[0]
    
    for chromosome, gene_assignments in genome_data.items():
        converted[chromosome] = {}
        
        for busco_id, hybrid_data in gene_assignments.items():
            # Start with positional bins
            combined_bins = hybrid_data['positional_bins'].copy()
            
            # Add ordinal data as a "bin" if available
            if hybrid_data['ordinal_window']:
                # Convert ordinal window to a pseudo-bin entry
                # Use 100% overlap since ordinal represents discrete rank position
                combined_bins.append((hybrid_data['ordinal_window'], 100.0))
            
            # If no positional bins but has ordinal, ensure we have data
            if not combined_bins and hybrid_data['ordinal_window']:
                combined_bins = [(hybrid_data['ordinal_window'], 100.0)]
            
            converted[chromosome][busco_id] = combined_bins
    
    return converted


def extract_gene_distribution(markov_profile, busco_id, pseudo_count=0.000001):
    """
    Extract gene distribution - handles both summary and direct gene entry formats.
    """
    gene_distribution = {}
    
    # Debug: Check if busco_id exists in profile
    print(f"DEBUG - Looking for busco_id: {busco_id}")
    print(f"DEBUG - Profile has {len(markov_profile)} bins")
    
    for bin_id, genes_data in markov_profile.items():
        probability = 0.0
        
        if busco_id in genes_data:
            print(f"DEBUG - Found {busco_id} in {bin_id}")
            print(f"DEBUG - genes_data[{busco_id}] keys: {list(genes_data[busco_id].keys())}")
            
            gene_entry = genes_data[busco_id]
            
            # Check different probability field names
            if 'average_percentage' in gene_entry:
                probability = gene_entry['average_percentage'] / 100.0
            elif 'frequency_percentage' in gene_entry:
                probability = gene_entry['frequency_percentage'] / 100.0
            elif 'overlap_percentage' in gene_entry:
                probability = gene_entry['overlap_percentage'] / 100.0
            else:
                # If gene exists but no percentage field, use summary data
                if 'position_summary' in genes_data:
                    probability = genes_data['position_summary']['average_percentage'] / 100.0
                elif 'rank_summary' in genes_data:
                    probability = genes_data['rank_summary']['average_percentage'] / 100.0
        
        gene_distribution[bin_id] = probability
    
    # Add pseudo-counts
    for bin_id in gene_distribution:
        gene_distribution[bin_id] += pseudo_count
    
    return gene_distribution


def extract_gene_distribution_hybrid_aware(profile, busco_id, pseudo_count=0.000001):
    """
    Extract gene distribution from hybrid profile, combining positional and ordinal data.
    """
    if 'positional_profile' in profile and 'ordinal_profile' in profile:
        # Hybrid format - combine both profiles
        combined_distribution = {}
        
        # Add positional data
        positional_dist = extract_gene_distribution(profile['positional_profile'], busco_id, pseudo_count)
        combined_distribution.update(positional_dist)
        
        # Add ordinal data
        ordinal_dist = extract_gene_distribution(profile['ordinal_profile'], busco_id, pseudo_count)
        combined_distribution.update(ordinal_dist)
        
        return combined_distribution
        
    elif 'positional_profile' in profile:
        # Hybrid format - positional only
        return extract_gene_distribution(profile['positional_profile'], busco_id, pseudo_count)
        
    else:
        # Legacy format
        return extract_gene_distribution(profile, busco_id, pseudo_count)

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