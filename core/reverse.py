import pandas as pd
import numpy as np
from collections import defaultdict

from .query_movement import extract_gene_distribution
from config.settings import CONFIG

def extract_movement_sequence(movement_results):
    """
    Extract ordered movement sequence from movement analysis results.
    
    Args:
        movement_results: Output from analyse_query_movements()
        
    Returns:
        list: [(gene_id, current_position, mean_movement)] ordered by current position
    """
    movement_sequence = []
    
    for gene_id, result in movement_results.items():
        if result['current_ranges']:
            current_position = result['current_ranges'][0][0]  # (bin_number, probability)
            mean_movement = result['movement_analysis']['mean_movement']
            movement_sequence.append((gene_id, current_position, mean_movement))
    
    movement_sequence.sort(key=lambda x: x[1])
    
    return movement_sequence


def detect_adjacency_inversions(movement_sequence):
    """
    Detect +/- adjacencies in the movement sequence.
    
    Args:
        movement_sequence: [(gene_id, position, movement)]
        
    Returns:
        list: [(index1, index2)] of adjacent pairs with opposite signs
    """
    adjacency_inversions = []
    
    for i in range(len(movement_sequence) - 1):
        _, _, movement1 = movement_sequence[i]
        _, _, movement2 = movement_sequence[i + 1]
        
        if ((movement1 > 0 and movement2 < 0) or 
            (movement1 < 0 and movement2 > 0)):
            adjacency_inversions.append((i, i + 1))
    
    return adjacency_inversions


# def detect_extended(movement_sequence):

#     flip_patterns = []
    
#     for i in range(len(movement_sequence) - 2):
#         pattern = []
        
#         for j in range(i, len(movement_sequence)):
#             _, _, movement = movement_sequence[j]
#             pattern.append(movement)
            
#             # +n, +(n-1), ..., +1, 0, -1, ..., -(n-1), -n
#             # Or +n, 0, -n
#             if len(pattern) >= 3:
#                 flip_size = detect_flip_in_pattern(pattern)
#                 if flip_size > 0:
#                     flip_patterns.append((i, j, flip_size))
#                     break
    
#     return flip_patterns



def detect_extended(movement_sequence):
    """
    Detect extended flip patterns - check all possible segments.
    """
    flip_patterns = []
    
    # Check all possible starting positions and lengths
    for start_i in range(len(movement_sequence)):
        for end_j in range(start_i + 1, len(movement_sequence)):  # Minimum length 2
            # Extract movement pattern for this segment
            pattern = []
            for k in range(start_i, end_j + 1):
                _, _, movement = movement_sequence[k]
                pattern.append(movement)
            
            # Test if this pattern is flippable
            flip_indicator = detect_flip_in_pattern(pattern)
            if flip_indicator > 0:
                flip_patterns.append((start_i, end_j, flip_indicator))
    
    # Sort by segment length (descending) to prioritize larger flips
    flip_patterns.sort(key=lambda x: (x[1] - x[0]), reverse=True)
    
    return flip_patterns


#OLD
# def detect_flip_in_pattern(pattern):

#     #Detect if a pattern represents a flip and return its size.

#     if len(pattern) < 3:
#         return 0
    
#     # +n, 0, -n
#     if len(pattern) == 3:
#         if (pattern[0] > 0 and pattern[1] == 0 and pattern[2] < 0 and
#             abs(pattern[0]) == abs(pattern[2])):
#             return abs(pattern[0])
    
#     # +n, +(n-1), ..., +1, 0, -1, ..., -(n-1), -n
#     mid_index = None
#     for i, val in enumerate(pattern):
#         if val == 0:
#             mid_index = i
#             break
    
#     if mid_index is None:
#         return 0
    
#     left_part = pattern[:mid_index]
#     right_part = pattern[mid_index + 1:]
    
#     if not all(val > 0 for val in left_part):
#         return 0
    
#     if not all(val < 0 for val in right_part):
#         return 0
    
#     if len(left_part) != len(right_part):
#         return 0
    
#     expected_left = list(range(len(left_part), 0, -1))
#     expected_right = list(range(-1, -len(right_part) - 1, -1))
    
#     if left_part == expected_left and right_part == expected_right:
#         return len(left_part)
    
#     return 0


def detect_flip_in_pattern(pattern):
    """
    Detect flip patterns based on corrected conditions with NO silent assumptions.
    
    Arching Rules:
    - Incrementalism: Values decrease toward center 
    - Symmetry: Equal counts on both sides of fulcrum
    
    Type 1: Simple Adjacent Patterns
    - Direct adjacent opposites: +n, -m (any magnitudes)
    - Single/multiple zeros between: +n, 0, -m or +n, 0, 0, 0, -m
    
    Type 2: Incremental Patterns  
    - Symmetric perfect decreasing: +3, +2, +1, 0, -1, -2, -3
    - Symmetric non-perfect decreasing: +7, +5, +2, 0, -1, -4, -6
    
    Type 4: Multiple Zeros Pattern
    - Zeros within positive/negative sections: +4, 0, +3, 0, -2, -4, -7
    """
    if len(pattern) < 2:
        return 0
    
    # Type 1: Simple Adjacent Patterns (length 2)
    if len(pattern) == 2:
        # Direct adjacent opposites: +n, -m (any magnitudes)
        if (pattern[0] > 0 and pattern[1] < 0) or (pattern[0] < 0 and pattern[1] > 0):
            return 1  # Return 1 to indicate simple flip (not magnitude-based)
        return 0
    
    # Type 1: Simple patterns with zeros between (length 3+)
    if len(pattern) >= 3:
        first_val = pattern[0]
        last_val = pattern[-1]
        
        # Must have opposite signs and non-zero
        if first_val != 0 and last_val != 0:
            if (first_val > 0 and last_val < 0) or (first_val < 0 and last_val > 0):
                # Check if all middle elements are zeros
                middle = pattern[1:-1]
                if all(val == 0 for val in middle):
                    return 1  # Simple flip with zeros between
    
    # Type 2 & 4: Incremental Patterns (with or without zeros in sections)
    # Find the fulcrum (first zero encountered)
    fulcrum_index = None
    for i, val in enumerate(pattern):
        if val == 0:
            fulcrum_index = i
            break
    
    # If no zero found, cannot be an incremental pattern
    if fulcrum_index is None:
        return 0
    
    # Split into left and right sections around fulcrum
    left_section = pattern[:fulcrum_index]
    right_section = pattern[fulcrum_index + 1:]
    
    # SYMMETRY CHECK: Must have equal counts on both sides
    if len(left_section) != len(right_section):
        return 0
    
    # If either section is empty, not a valid incremental pattern
    if len(left_section) == 0:
        return 0
    
    # Separate zeros from non-zeros in each section
    left_nonzero = [val for val in left_section if val != 0]
    right_nonzero = [val for val in right_section if val != 0]
    
    # Must have at least one non-zero value on each side for incremental
    if not left_nonzero or not right_nonzero:
        return 0
    
    # SIGN CHECK: Left section non-zeros must be positive, right must be negative
    if not all(val > 0 for val in left_nonzero):
        return 0
    if not all(val < 0 for val in right_nonzero):
        return 0
    
    # INCREMENTALISM CHECK: Values should decrease toward center
    # Left side: should be non-increasing (allowing equal values)
    for i in range(len(left_nonzero) - 1):
        if left_nonzero[i] < left_nonzero[i + 1]:  # Must be decreasing or equal
            return 0
    
    # Right side: should be non-decreasing in absolute value (toward center)
    right_abs = [abs(val) for val in right_nonzero]
    for i in range(len(right_abs) - 1):
        if right_abs[i] < right_abs[i + 1]:  # Must be decreasing or equal in magnitude
            return 0
    
    # Check for Type 2.1: Perfect incremental patterns
    if is_perfect_incremental(left_nonzero, right_nonzero):
        return len(left_section)  # Use section length as flip indicator
    
    # Type 2.2 & 4: Non-perfect incremental or zeros within sections
    return len(left_section)  # Use section length as flip indicator


def is_perfect_incremental(left_nonzero, right_nonzero):
    """
    Check if pattern is perfect incremental: +n, +(n-1), ..., +1, 0, -1, ..., -(n-1), -n
    """
    if len(left_nonzero) != len(right_nonzero):
        return False
    
    n = len(left_nonzero)
    expected_left = list(range(n, 0, -1))  # [n, n-1, ..., 1]
    expected_right = list(range(-1, -n - 1, -1))  # [-1, -2, ..., -n]
    
    return left_nonzero == expected_left and right_nonzero == expected_right



#OLD
# def apply_adjacency_inversion(movement_sequence, index1, index2):


#     updated_sequence = movement_sequence.copy()
    
#     gene1_id, pos1, move1 = updated_sequence[index1]
#     gene2_id, pos2, move2 = updated_sequence[index2]
    
#     ##debug print
#     print(f"  BEFORE: {gene1_id}({move1}), {gene2_id}({move2})")
    
    
#     updated_sequence[index1] = (gene2_id, pos1, move2)
#     updated_sequence[index2] = (gene1_id, pos2, move1)
    
#     new_move1 = move1 - 1 if move1 > 0 else move1 + 1 if move1 < 0 else 0
#     new_move2 = move2 - 1 if move2 > 0 else move2 + 1 if move2 < 0 else 0
    
#     updated_sequence[index1] = (gene2_id, pos1, new_move2)
#     updated_sequence[index2] = (gene1_id, pos2, new_move1)
    
#     ##debug print
#     print(f"  AFTER:  {gene2_id}({new_move2}), {gene1_id}({new_move1})")
    
#     inversion_record = {
#         'type': 'adjacency',
#         'positions': [pos1, pos2],
#         'genes': [gene1_id, gene2_id],
#         'gene_inversions': 2
#     }
    
#     return updated_sequence, inversion_record


def apply_adjacency_inversion(movement_sequence, index1, index2):
    """
    Apply adjacency inversion with correct movement updates.
    """
    updated_sequence = movement_sequence.copy()
    
    gene1_id, pos1, move1 = updated_sequence[index1]
    gene2_id, pos2, move2 = updated_sequence[index2]
    
    print(f"  ADJACENCY BEFORE: {gene1_id}({move1}) <-> {gene2_id}({move2})")
    
    # Step 1: Swap genes (keep original positions)
    updated_sequence[index1] = (gene2_id, pos1, move2)
    updated_sequence[index2] = (gene1_id, pos2, move1)
    
    # Step 2: Update movements (each gene moved 1 position)
    new_move1 = move1 - 1 if move1 > 0 else move1 + 1 if move1 < 0 else 0
    new_move2 = move2 - 1 if move2 > 0 else move2 + 1 if move2 < 0 else 0
    
    updated_sequence[index1] = (gene2_id, pos1, new_move2)
    updated_sequence[index2] = (gene1_id, pos2, new_move1)
    
    print(f"  ADJACENCY AFTER:  {gene2_id}({new_move2}) <-> {gene1_id}({new_move1})")
    
    inversion_record = {
        'type': 'adjacency',
        'positions': [pos1, pos2],
        'genes': [gene1_id, gene2_id],
        'gene_inversions': 2
    }
    
    return updated_sequence, inversion_record


#OLD
# def apply_flip_inversion(movement_sequence, start_index, end_index, flip_size):

#     updated_sequence = movement_sequence.copy()
    
#     segment = updated_sequence[start_index:end_index + 1]
    
#     reversed_genes = [(gene_id, pos, move) for gene_id, _, move in reversed([s for s in segment])]
    
#     positions = [pos for _, pos, _ in segment]
#     for i, (gene_id, _, move) in enumerate(reversed_genes):
#         reversed_genes[i] = (gene_id, positions[i], move)
    
#     for i in range(len(reversed_genes)):
#         gene_id, pos, move = reversed_genes[i]
#         new_move = move - flip_size if move > 0 else move + flip_size if move < 0 else 0
#         reversed_genes[i] = (gene_id, pos, new_move)
    
#     updated_sequence[start_index:end_index + 1] = reversed_genes
    
#     positions = [pos for _, pos, _ in segment]
#     genes = [gene_id for gene_id, _, _ in segment]
    
#     inversion_record = {
#         'type': 'flip',
#         'flip_size': flip_size,
#         'positions': positions,
#         'genes': genes,
#         'gene_inversions': len(segment)
#     }
    
#     return updated_sequence, inversion_record


def apply_flip_inversion(movement_sequence, start_index, end_index, flip_indicator):
    """
    Apply flip inversion with CORRECT positional distance calculations.
    
    Process:
    1. Reverse the gene order in the segment
    2. Calculate how far each gene moved from its original position
    3. Update each gene's movement by the distance it traveled
    """
    updated_sequence = movement_sequence.copy()
    
    # Extract the segment to flip
    segment = updated_sequence[start_index:end_index + 1]
    segment_length = len(segment)
    
    print(f"  FLIP BEFORE: {[(gene_id, move) for gene_id, _, move in segment]}")
    
    # Step 1: Create the reversed segment
    reversed_segment = []
    original_positions = [pos for _, pos, _ in segment]
    
    # Reverse gene order but keep original positions
    for i, (gene_id, _, move) in enumerate(reversed(segment)):
        # Gene goes to the position at index i in the segment
        new_position = original_positions[i]
        reversed_segment.append((gene_id, new_position, move))
    
    # Step 2: Calculate movement updates based on positional distances
    final_segment = []
    for i, (gene_id, pos, old_move) in enumerate(reversed_segment):
        # Calculate how far this gene moved during the flip
        original_index = (segment_length - 1) - i  # Where it was in the original segment
        new_index = i  # Where it is now in the reversed segment
        distance_moved = abs(new_index - original_index)
        
        # Update movement: reduce by the distance moved
        if old_move > 0:
            new_move = old_move - distance_moved
        elif old_move < 0:
            new_move = old_move + distance_moved  
        else:
            new_move = 0  # Zero stays zero
        
        final_segment.append((gene_id, pos, new_move))
        
        print(f"    {gene_id}: moved {distance_moved} positions, {old_move} → {new_move}")
    
    # Step 3: Update the sequence
    updated_sequence[start_index:end_index + 1] = final_segment
    
    print(f"  FLIP AFTER:  {[(gene_id, move) for gene_id, _, move in final_segment]}")
    
    # Create inversion record
    positions = [pos for _, pos, _ in segment]
    genes = [gene_id for gene_id, _, _ in segment]
    
    inversion_record = {
        'type': 'flip',
        'flip_indicator': flip_indicator,
        'segment_length': segment_length,
        'positions': positions,
        'genes': genes,
        'gene_inversions': len(segment)
    }
    
    return updated_sequence, inversion_record


# def iterative_detection(movement_sequence, max_iterations=100):

#     current_sequence = movement_sequence.copy()
#     inversion_events = []
#     iteration = 0
    
#     while iteration < max_iterations:
#         iteration += 1
        
        
#         #debug
#         print(f"Iteration {iteration}: Starting...")
#         ##
        
        
#         applied_inversion = False
        
#         ##debug
#         movements = [move for _, _, move in current_sequence if move != 0]
#         print(f"  Non-zero movements: {len(movements)} (sample: {movements[:5]})")
        
#         ##
        
#         flip_patterns = detect_extended(current_sequence)
#         if flip_patterns:
            
#                 ##debig print
#                 print(f"  Found {len(flip_patterns)} flip patterns")
                
#                 ###
                
                
#                 start_idx, end_idx, flip_size = flip_patterns[0]
#                 current_sequence, inversion_record = apply_flip_inversion(
#                     current_sequence, start_idx, end_idx, flip_size
#                 )
#                 inversion_record['iteration'] = iteration
#                 inversion_events.append(inversion_record)
#                 applied_inversion = True
        
#         elif not applied_inversion:    
            
#             ##debug print
#             print(f"  Checking adjacency patterns...")
            
#             adjacency_inversions = detect_adjacency_inversions(current_sequence)
#             if adjacency_inversions:
                
#                     ##debug print
#                     print(f"  Found {len(adjacency_inversions)} adjacency patterns")

#                     index1, index2 = adjacency_inversions[0]
#                     current_sequence, inversion_record = apply_adjacency_inversion(
#                         current_sequence, index1, index2
#                     )
#                     inversion_record['iteration'] = iteration
#                     inversion_events.append(inversion_record)
#                     applied_inversion = True
#         if not applied_inversion:
            
#             ##debug print
#             print(f"  No inversions found - terminating")

#             break
        
#         ##debug print
#         print(f"Iteration {iteration}: Applied inversion")
#     print(f"Iterative detection completed after {iteration} iterations")
    
#     total_events = len(inversion_events)
#     total_gene_inversions = sum(event['gene_inversions'] for event in inversion_events)
#     adjacency_events = sum(1 for event in inversion_events if event['type'] == 'adjacency')
#     flip_events = sum(1 for event in inversion_events if event['type'] == 'flip')
    
#     return {
#         'final_sequence': current_sequence,
#         'inversion_events': inversion_events,
#         'iterations': iteration,
#         'total_events': total_events,
#         'total_gene_inversions': total_gene_inversions,
#         'adjacency_events': adjacency_events,
#         'flip_events': flip_events,
#         'converged': not applied_inversion
#     }









##Detect in Batches



def find_non_overlapping_adjacencies(adjacency_inversions):
    """
    Find adjacencies that don't share any indices to avoid conflicts.
    
    Args:
        adjacency_inversions: [(index1, index2)] list of adjacent pairs
        
    Returns:
        list: [(index1, index2)] non-overlapping pairs that can be applied safely
    """
    non_overlapping = []
    used_indices = set()
    
    # Sort by first index to process left-to-right
    sorted_adjacencies = sorted(adjacency_inversions, key=lambda x: x[0])
    
    for index1, index2 in sorted_adjacencies:
        # Check if either index is already used
        if index1 not in used_indices and index2 not in used_indices:
            non_overlapping.append((index1, index2))
            used_indices.add(index1)
            used_indices.add(index2)
    
    return non_overlapping


def iterative_detection(movement_sequence, max_iterations=1000):
    """
    Optimized iterative detection with batch processing of non-overlapping inversions.
    
    Args:
        movement_sequence: [(gene_id, position, movement)]
        max_iterations: Maximum number of iterations
        
    Returns:
        dict: Complete inversion analysis results
    """
    print(f"Starting optimized iterative detection with {len(movement_sequence)} genes")
    
    current_sequence = movement_sequence.copy()
    inversion_events = []
    iteration = 0
    
    while iteration < max_iterations:
        iteration += 1
        print(f"Iteration {iteration}: Starting...")
        applied_inversion = False
        batch_count = 0
        
        # Debug: Check what movements exist
        movements = [move for _, _, move in current_sequence if move != 0]
        print(f"  Non-zero movements: {len(movements)} (sample: {movements[:5] if movements else []})")
        
        # Check for flip patterns first (more efficient)
        print(f"  Checking flip patterns...")
        flip_patterns = detect_extended(current_sequence)
        if flip_patterns:
            print(f"  Found {len(flip_patterns)} flip patterns")
            
            # Apply first flip pattern (flips are usually more complex, handle one at a time)
            start_idx, end_idx, flip_size = flip_patterns[0]
            current_sequence, inversion_record = apply_flip_inversion(
                current_sequence, start_idx, end_idx, flip_size
            )
            inversion_record['iteration'] = iteration
            inversion_record['batch_position'] = 0
            inversion_events.append(inversion_record)
            applied_inversion = True
            print(f"  Applied flip inversion: {inversion_record['genes']}")
        
        # If no flip patterns, check for adjacency inversions with batch processing
        elif not applied_inversion:
            print(f"  Checking adjacency patterns...")
            adjacency_inversions = detect_adjacency_inversions(current_sequence)
            
            if adjacency_inversions:
                print(f"  Found {len(adjacency_inversions)} adjacency patterns")
                
                # Find non-overlapping adjacencies for batch processing
                non_overlapping = find_non_overlapping_adjacencies(adjacency_inversions)
                print(f"  Non-overlapping adjacencies: {len(non_overlapping)}")
                
                # Apply all non-overlapping adjacencies in this iteration
                for batch_idx, (index1, index2) in enumerate(non_overlapping):
                    # Debug: Show what's being swapped
                    gene1_id, pos1, move1 = current_sequence[index1]
                    gene2_id, pos2, move2 = current_sequence[index2]
                    print(f"    Batch {batch_idx}: {gene1_id}({move1}) <-> {gene2_id}({move2})")
                    
                    # Apply adjacency inversion
                    current_sequence, inversion_record = apply_adjacency_inversion(
                        current_sequence, index1, index2
                    )
                    
                    # Record the event
                    inversion_record['iteration'] = iteration
                    inversion_record['batch_position'] = batch_idx
                    inversion_record['total_in_batch'] = len(non_overlapping)
                    inversion_events.append(inversion_record)
                    batch_count += 1
                
                applied_inversion = True
                print(f"  Applied {batch_count} adjacency inversions in batch")
        
        # Termination condition
        if not applied_inversion:
            print(f"  No inversions found - terminating")
            break
        
        # Progress check: stop if very few large movements remain
        if iteration % 50 == 0:
            large_movements = sum(1 for _, _, move in current_sequence if abs(move) > 2.0)
            print(f"  Progress check: {large_movements} genes with |movement| > 2.0")
            if large_movements < 50:  # Most work done
                print(f"  Early termination: most large movements resolved")
                break
        
        print(f"Iteration {iteration}: Applied {batch_count if batch_count > 0 else 1} inversions")
    
    # Calculate summary statistics
    total_events = len(inversion_events)
    total_gene_inversions = sum(event['gene_inversions'] for event in inversion_events)
    adjacency_events = sum(1 for event in inversion_events if event['type'] == 'adjacency')
    flip_events = sum(1 for event in inversion_events if event['type'] == 'flip')
    
    print(f"Optimized iterative detection completed after {iteration} iterations")
    print(f"  Total events: {total_events}")
    print(f"  Total gene inversions: {total_gene_inversions}")
    print(f"  Adjacency events: {adjacency_events}, Flip events: {flip_events}")
    
    return {
        'final_sequence': current_sequence,
        'inversion_events': inversion_events,
        'iterations': iteration,
        'total_events': total_events,
        'total_gene_inversions': total_gene_inversions,
        'adjacency_events': adjacency_events,
        'flip_events': flip_events,
        'converged': not applied_inversion
    }





















def get_permutable_positions(gene_distribution, min_probability=None, max_positions=None):
    """
    Get permutable positions for a gene based on probability thresholds.
    
    Args:
        gene_distribution: {bin_id: normalized_probability}
        min_probability: Minimum probability threshold (uses CONFIG if None)
        max_positions: Maximum positions to return (uses CONFIG if None)
        
    Returns:
        list: [(bin_number, probability)] sorted by probability (descending)
    """
    if min_probability is None:
        min_probability = CONFIG['permutable_positions_threshold']
    if max_positions is None:
        max_positions = CONFIG['max_permutable_positions']
    
    permutable_positions = []
    
    for bin_id, probability in gene_distribution.items():
        if probability >= min_probability:
            try:
                bin_number = int(bin_id.split('_bin_')[1])
                permutable_positions.append((bin_number, probability))
            except (IndexError, ValueError):
                continue
    
    permutable_positions.sort(key=lambda x: (-x[1], x[0]))
    
    return permutable_positions[:max_positions]


def calculate_position_probability(gene_id, target_position, markov_profile):
    """
    Calculate probability of a gene being at a specific position.
    
    Args:
        gene_id: Gene identifier
        target_position: Target bin number
        markov_profile: Markov profile
        
    Returns:
        float: Probability of gene being at target position
    """
    gene_distribution = extract_gene_distribution(markov_profile, gene_id)
    
    for bin_id, probability in gene_distribution.items():
        try:
            bin_number = int(bin_id.split('_bin_')[1])
            if bin_number == target_position:
                return probability
        except (IndexError, ValueError):
            continue
    
    return 0.0 


def evaluate_inversion_step_probability(inversion_step, markov_profile, threshold=None):
    """
    Evaluate if an inversion step meets probability thresholds.
    
    Args:
        inversion_step: Inversion event record with genes and positions
        markov_profile: Markov profile
        threshold: Probability threshold (uses CONFIG if None)
        
    Returns:
        dict: Evaluation results
    """
    if threshold is None:
        threshold = CONFIG['permutable_positions_threshold']
    
    gene_probabilities = {}
    passed_genes = []
    failed_genes = []
    
    genes = inversion_step['genes']
    positions = inversion_step['positions']
    
    for i, gene_id in enumerate(genes):
        final_position = positions[i]
        probability = calculate_position_probability(gene_id, final_position, markov_profile)
        
        gene_probabilities[gene_id] = {
            'position': final_position,
            'probability': probability,
            'passed': probability >= threshold
        }
        
        if probability >= threshold:
            passed_genes.append(gene_id)
        else:
            failed_genes.append(gene_id)
    
    overall_probability = 1.0
    for prob_data in gene_probabilities.values():
        overall_probability *= prob_data['probability']
    
    bit_score = -math.log2(overall_probability) if overall_probability > 0 else float('inf')
    
    return {
        'gene_probabilities': gene_probabilities,
        'passed_genes': passed_genes,
        'failed_genes': failed_genes,
        'overall_probability': overall_probability,
        'bit_score': bit_score,
        'threshold_passed': len(failed_genes) == 0
    }


def generate_gene_specific_pathways(target_gene, current_position, markov_profile, other_genes_positions, initial_sequence):
    """
    Generate alternative pathways for moving a specific gene to its permutable positions.
    
    Args:
        target_gene: Gene that needs alternative pathway
        current_position: Current position of the target gene
        markov_profile: Markov profile
        other_genes_positions: {gene_id: current_position} for other genes
        initial_sequence: Initial movement sequence [(gene_id, position, movement)]
        
    Returns:
        list: Alternative pathways ranked by combined bit score
    """
    # Get permutable positions for target gene
    gene_distribution = extract_gene_distribution(markov_profile, target_gene)
    permutable_positions = get_permutable_positions(gene_distribution)
    
    alternative_pathways = []
    
    for target_position, target_probability in permutable_positions:
        if target_position == current_position:
            continue  # Already at this position
        
        # Generate inversion steps to move target_gene from current to target position
        pathway_steps = generate_inversion_steps_for_gene(
            target_gene, current_position, target_position, other_genes_positions
        )
        
        # Evaluate pathway including target gene and all affected genes
        pathway_evaluation = evaluate_pathway_steps(pathway_steps, markov_profile, initial_sequence)
        
        alternative_pathways.append({
            'target_gene': target_gene,
            'target_position': target_position,
            'target_probability': target_probability,
            'pathway_steps': pathway_steps,
            'evaluation': pathway_evaluation
        })
    
    # Sort by combined bit score (lower is better)
    alternative_pathways.sort(key=lambda x: x['evaluation']['combined_bit_score'])
    
    return alternative_pathways


def generate_inversion_steps_for_gene(gene_id, current_pos, target_pos, other_genes_positions):
    """
    Generate inversion steps to move a specific gene from current to target position.
    Uses the full iterative detection algorithm with gene-specific constraint.
    
    Args:
        gene_id: Gene to move
        current_pos: Current position
        target_pos: Target position
        other_genes_positions: {gene_id: position} for other genes
        
    Returns:
        list: Inversion steps that involve the target gene
    """
    movement_sequence = []
    
    target_movement = target_pos - current_pos
    movement_sequence.append((gene_id, current_pos, target_movement))
    
    for other_gene, pos in other_genes_positions.items():
        movement_sequence.append((other_gene, pos, 0.0)) 
    
    movement_sequence.sort(key=lambda x: x[1])
    
    gene_specific_events = iterative_detection_gene_specific(movement_sequence, gene_id, target_pos)
    
    return gene_specific_events


def iterative_detection_gene_specific(movement_sequence, target_gene, target_position, max_iterations=100):
    """
    Run iterative detection but only apply inversions that involve the target gene.
    
    Args:
        movement_sequence: [(gene_id, position, movement)]
        target_gene: Gene that must be involved in all inversions
        target_position: Target position for the gene
        max_iterations: Maximum iterations
        
    Returns:
        list: Inversion events that involve target gene
    """
    current_sequence = movement_sequence.copy()
    gene_specific_events = []
    iteration = 0
    
    while iteration < max_iterations:
        iteration += 1
        applied_inversion = False
        
        flip_patterns = detect_extended(current_sequence)
        if flip_patterns:
  
            for start_idx, end_idx, flip_size in flip_patterns:
                segment_genes = [current_sequence[i][0] for i in range(start_idx, end_idx + 1)]
                
                if target_gene in segment_genes:
                    current_sequence, inversion_record = apply_flip_inversion(
                        current_sequence, start_idx, end_idx, flip_size
                    )
                    inversion_record['iteration'] = iteration
                    inversion_record['gene_specific_target'] = target_gene
                    gene_specific_events.append(inversion_record)
                    applied_inversion = True
                    break
        

        if not applied_inversion:
            adjacency_inversions = detect_adjacency_inversions(current_sequence)
            if adjacency_inversions:

                for index1, index2 in adjacency_inversions:
                    gene1 = current_sequence[index1][0]
                    gene2 = current_sequence[index2][0]
                    
                    if target_gene == gene1 or target_gene == gene2:
                        current_sequence, inversion_record = apply_adjacency_inversion(
                            current_sequence, index1, index2
                        )
                        inversion_record['iteration'] = iteration
                        inversion_record['gene_specific_target'] = target_gene
                        gene_specific_events.append(inversion_record)
                        applied_inversion = True
                        break
        

        target_gene_current_pos = None
        for gene_name, pos, movement in current_sequence:
            if gene_name == target_gene:
                target_gene_current_pos = pos
                break
        
        if target_gene_current_pos == target_position:
            break 
        

        if not applied_inversion:
            break
    
    return gene_specific_events


def evaluate_pathway_steps(pathway_steps, markov_profile, initial_sequence):
    """
    Evaluate all steps in a pathway, including target gene and all affected genes.
    
    Args:
        pathway_steps: List of inversion steps
        markov_profile: Markov profile
        initial_sequence: Initial movement sequence [(gene_id, position, movement)]
        
    Returns:
        dict: Pathway evaluation
    """
    step_evaluations = []
    total_bit_score = 0.0
    all_affected_genes = set()
    
    for step in pathway_steps:
        step_eval = evaluate_inversion_step_probability(step, markov_profile)
        step_evaluations.append(step_eval)
        total_bit_score += step_eval['bit_score']
        
        all_affected_genes.update(step['genes'])

    final_positions = track_all_gene_final_positions(pathway_steps, initial_sequence)
    

    combined_probability = 1.0
    affected_gene_probabilities = {}
    
    for gene_id in all_affected_genes:
        final_position = final_positions.get(gene_id)
        if final_position is not None:
            gene_probability = calculate_position_probability(gene_id, final_position, markov_profile)
            
            affected_gene_probabilities[gene_id] = {
                'final_position': final_position,
                'probability': gene_probability
            }
            
            combined_probability *= gene_probability
    
    combined_bit_score = -math.log2(combined_probability) if combined_probability > 0 else float('inf')
    
    return {
        'step_evaluations': step_evaluations,
        'total_bit_score': total_bit_score,
        'combined_bit_score': combined_bit_score,
        'affected_genes': list(all_affected_genes),
        'affected_gene_probabilities': affected_gene_probabilities,
        'final_positions': final_positions,
        'pathway_valid': combined_probability > 0
    }


def get_gene_final_position(gene_id, pathway_steps, initial_sequence):
    """
    Get the final position of a gene after all pathway steps by applying each step sequentially.
    
    Args:
        gene_id: Gene identifier
        pathway_steps: List of inversion steps
        initial_sequence: Initial movement sequence [(gene_id, position, movement)]
        
    Returns:
        int: Final position of the gene
    """
    current_sequence = initial_sequence.copy()
    
    for step in pathway_steps:
        if step['type'] == 'flip':
            # Find the positions in current sequence
            start_pos = min(step['positions'])
            end_pos = max(step['positions'])
            

            start_idx = None
            end_idx = None
            for i, (gene, pos, movement) in enumerate(current_sequence):
                if pos == start_pos and start_idx is None:
                    start_idx = i
                if pos == end_pos:
                    end_idx = i
            
            if start_idx is not None and end_idx is not None:
                current_sequence, _ = apply_flip_inversion(
                    current_sequence, start_idx, end_idx, step.get('flip_size', 1)
                )
        
        elif step['type'] == 'adjacency':

            genes_in_step = step['genes']
            if len(genes_in_step) >= 2:
                gene1, gene2 = genes_in_step[0], genes_in_step[1]
                
                index1 = None
                index2 = None
                for i, (gene, pos, movement) in enumerate(current_sequence):
                    if gene == gene1:
                        index1 = i
                    elif gene == gene2:
                        index2 = i
                
                if index1 is not None and index2 is not None:
                    current_sequence, _ = apply_adjacency_inversion(
                        current_sequence, index1, index2
                    )
    
    for gene, pos, movement in current_sequence:
        if gene == gene_id:
            return pos
    
    return None 


def track_all_gene_final_positions(pathway_steps, initial_sequence):
    """
    Track final positions of all genes after pathway steps.
    
    Args:
        pathway_steps: List of inversion steps
        initial_sequence: Initial movement sequence [(gene_id, position, movement)]
        
    Returns:
        dict: {gene_id: final_position}
    """
    current_sequence = initial_sequence.copy()
    
    for step in pathway_steps:
        if step['type'] == 'flip':
            start_pos = min(step['positions'])
            end_pos = max(step['positions'])
            
            start_idx = None
            end_idx = None
            for i, (gene, pos, movement) in enumerate(current_sequence):
                if pos == start_pos and start_idx is None:
                    start_idx = i
                if pos == end_pos:
                    end_idx = i
            
            if start_idx is not None and end_idx is not None:
                current_sequence, _ = apply_flip_inversion(
                    current_sequence, start_idx, end_idx, step.get('flip_size', 1)
                )
        
        elif step['type'] == 'adjacency':
            genes_in_step = step['genes']
            if len(genes_in_step) >= 2:
                gene1, gene2 = genes_in_step[0], genes_in_step[1]
                
                index1 = None
                index2 = None
                for i, (gene, pos, movement) in enumerate(current_sequence):
                    if gene == gene1:
                        index1 = i
                    elif gene == gene2:
                        index2 = i
                
                if index1 is not None and index2 is not None:
                    current_sequence, _ = apply_adjacency_inversion(
                        current_sequence, index1, index2
                    )
    
    final_positions = {}
    for gene, pos, movement in current_sequence:
        final_positions[gene] = pos
    
    return final_positions

def probability_weighted_inversion_analysis(movement_results, markov_profile):
    """
    Perform probability-weighted inversion analysis with alternative pathways.
    
    Args:
        movement_results: Output from movement analysis
        markov_profile: Markov profile
        
    Returns:
        dict: Complete probability-weighted analysis
    """
    standard_analysis = check_events_iteration(movement_results)
    
    # Evaluate each inversion step for probability
    evaluated_steps = []
    problematic_genes = []
    
    for event in standard_analysis['inversion_events']:
        step_eval = evaluate_inversion_step_probability(event, markov_profile)
        
        evaluated_steps.append({
            'original_event': event,
            'probability_evaluation': step_eval
        })
        
        if not step_eval['threshold_passed']:
            problematic_genes.extend(step_eval['failed_genes'])
    
    alternative_pathways = {}
    
    movement_sequence = extract_movement_sequence(movement_results)
    current_positions = {gene_id: pos for gene_id, pos, _ in movement_sequence}
    
    for gene_id in set(problematic_genes):
        if gene_id in current_positions:
            other_positions = {g: p for g, p in current_positions.items() if g != gene_id}
            
            alternatives = generate_gene_specific_pathways(
                gene_id, current_positions[gene_id], markov_profile, other_positions, movement_sequence
            )
            
            alternative_pathways[gene_id] = alternatives
            
    total_standard_bit_score = sum(step['probability_evaluation']['bit_score'] 
                                  for step in evaluated_steps)
    
    return {
        'standard_analysis': standard_analysis,
        'evaluated_steps': evaluated_steps,
        'total_standard_bit_score': total_standard_bit_score,
        'problematic_genes': list(set(problematic_genes)),
        'alternative_pathways': alternative_pathways,
        'has_alternatives': len(alternative_pathways) > 0
    }


def check_events_iteration(movement_results):

    movement_sequence = extract_movement_sequence(movement_results)
    
    inversion_analysis = iterative_detection(movement_sequence)
    
    inversion_analysis['original_sequence'] = movement_sequence
    
    return inversion_analysis


if __name__ == "__main__":
    print("done")
