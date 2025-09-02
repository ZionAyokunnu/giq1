import pandas as pd
import numpy as np
from collections import defaultdict
import math
from typing import Dict

from .query_movement import extract_gene_distribution
from .transposition import detect_and_apply_translocations
from config.settings import CONFIG


def extract_movement_sequence(movement_results):
    """
    Extract ordered movement sequence from movement analysis results.
    
      Returns:
        list: [(gene_id, current_position, mean_movement)] ordered by current position

    """
    movement_sequence = []
    
    for chromosome, gene_results in movement_results.items():
        for gene_id, result in gene_results.items():
            if result['current_ranges']:
                current_position = result['current_ranges'][0][0] 
                mean_movement = result['movement_analysis']['mean_movement']
                movement_sequence.append((gene_id, current_position, mean_movement))
    
    
    movement_sequence.sort(key=lambda x: x[1])
    
    return movement_sequence




def detect_adjacency_inversions(movement_sequence):
    """
    Detect +/- adjacencies in the movement sequence.
    
    Args:
        movement_sequence: [(gene_id, position, movement, target_position)]
        
    Returns:
        list: [(index1, index2)] of adjacent pairs with opposite signs
    """
    adjacency_inversions = []
    
    for i in range(len(movement_sequence) - 1):
        _, _, movement1, _ = movement_sequence[i]
        _, _, movement2, _ = movement_sequence[i + 1]
        
        if ((movement1 > 0 and movement2 < 0) or 
            (movement1 < 0 and movement2 > 0)):
            adjacency_inversions.append((i, i + 1))
    
    return adjacency_inversions


def calculate_total_movement(sequence):
    """Calculate total absolute movement for a sequence."""
    return sum(abs(move) for _, _, move, _ in sequence)



def detect_extended(movement_sequence):
    """
    Detect extended flip patterns including subsequences.
    """
    flip_patterns = []
    
    # Check all possible segments (including subsequences)
    for start_i in range(len(movement_sequence)):
        for end_j in range(start_i + 1, len(movement_sequence)):  # Minimum length 2
            # Extract movement pattern for this segment
            pattern = []
            for k in range(start_i, end_j + 1):
                _, _, movement, _ = movement_sequence[k]
                pattern.append(movement)
            
            # Test if this pattern is flippable
            flip_indicator = detect_flip_in_pattern(pattern)
            if flip_indicator > 0:
                flip_patterns.append((start_i, end_j, flip_indicator))
    
    # Sort by flip indicator (larger patterns first), then by length
    flip_patterns.sort(key=lambda x: (x[2], x[1] - x[0]), reverse=True)
    
    return flip_patterns

def detect_flip_in_pattern(pattern):
    """
    Detect flip patterns based on comprehensive classification.
    
    Returns flip_indicator > 0 if pattern is valid, 0 if not valid.
    """
    if len(pattern) < 2:
        return 0
    
    # Type 1: Simple Adjacent Patterns (length 2 only)
    if len(pattern) == 2:
        if (pattern[0] > 0 and pattern[1] < 0) or (pattern[0] < 0 and pattern[1] > 0):
            return 1
        return 0
    
    # Type 2A: Odd length incremental with natural fulcrum
    if len(pattern) % 2 == 1:
        return detect_odd_length_incremental(pattern)
    
    # Type 2B: Even length incremental with synthetic fulcrum  
    else:
        return detect_even_length_incremental(pattern)


def detect_odd_length_incremental(pattern):
    """Handle odd-length patterns with natural fulcrum."""
    fulcrum_index = len(pattern) // 2
    
    # Fulcrum must be zero
    if pattern[fulcrum_index] != 0:
        return 0
    
    # Split into left and right sections
    left_section = pattern[:fulcrum_index]
    right_section = pattern[fulcrum_index + 1:]
    
    # Must have equal non-empty sections
    if len(left_section) != len(right_section) or len(left_section) == 0:
        return 0
    
    return validate_incremental_sections(left_section, right_section)


def detect_even_length_incremental(pattern):
    """Handle even-length patterns with synthetic fulcrum."""
    mid_point = len(pattern) // 2
    
    # Split into left and right sections
    left_section = pattern[:mid_point]
    right_section = pattern[mid_point:]
    
    # Must have equal non-empty sections
    if len(left_section) != len(right_section) or len(left_section) == 0:
        return 0
    
    return validate_incremental_sections(left_section, right_section)


def validate_incremental_sections(left_section, right_section):
    """Validate that sections meet incremental pattern requirements."""
    # Get non-zero values from each section
    left_nonzero = [val for val in left_section if val != 0]
    right_nonzero = [val for val in right_section if val != 0]
    
    # Must have at least one non-zero value on each side
    if not left_nonzero or not right_nonzero:
        return 0
    
    # Endpoints (first left, last right) must be non-zero and opposite signs
    if left_section[0] == 0 or right_section[-1] == 0:
        return 0
    
    if not ((left_section[0] > 0 and right_section[-1] < 0) or 
            (left_section[0] < 0 and right_section[-1] > 0)):
        return 0
    
    # Sign separation: left section non-zeros same sign, right section opposite
    if left_section[0] > 0:
        # Left should be positive, right should be negative
        if not all(val >= 0 for val in left_section):
            return 0
        if not all(val <= 0 for val in right_section):
            return 0
    else:
        # Left should be negative, right should be positive  
        if not all(val <= 0 for val in left_section):
            return 0
        if not all(val >= 0 for val in right_section):
            return 0
    
    # Incrementalism: magnitudes should decrease toward center
    if not check_incrementalism(left_nonzero, right_nonzero):
            return 0
    
    return len(left_section)  # Return section length as flip indicator


def check_incrementalism(left_nonzero, right_nonzero):
    """Check if magnitudes decrease toward center."""
    # Left side: should decrease in magnitude from outside to center
    for i in range(len(left_nonzero) - 1):
        if abs(left_nonzero[i]) < abs(left_nonzero[i + 1]):
            return False
    
    # Right side: should decrease in magnitude from center to outside
    # So when read from outside to center, should also decrease
    right_reversed = list(reversed(right_nonzero))
    for i in range(len(right_reversed) - 1):
        if abs(right_reversed[i]) < abs(right_reversed[i + 1]):
            return False
    
    return True


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


def apply_adjacency_inversion(movement_sequence, index1, index2):
    """
    Apply adjacency inversion with correct movement updates.
    """
    updated_sequence = movement_sequence.copy()
    
    gene1_id, pos1, move1, target1 = updated_sequence[index1]
    gene2_id, pos2, move2, target2 = updated_sequence[index2]
    
    print(f"  ADJACENCY BEFORE: {gene1_id}({move1}) <-> {gene2_id}({move2})")
    
    # Step 1: Swap genes (keep original positions)
    updated_sequence[index1] = (gene2_id, pos1, move2, target2)
    updated_sequence[index2] = (gene1_id, pos2, move1, target1)
    
    # Step 2: Update movements (each gene moved 1 position)
    new_move1 = move1 - 1 if move1 > 0 else move1 + 1 if move1 < 0 else 0
    new_move2 = move2 - 1 if move2 > 0 else move2 + 1 if move2 < 0 else 0
    
    updated_sequence[index1] = (gene2_id, pos1, new_move2, target2)
    updated_sequence[index2] = (gene1_id, pos2, new_move1, target1)
    
    print(f"  ADJACENCY AFTER:  {gene2_id}({new_move2}) <-> {gene1_id}({new_move1})")
    
    inversion_record = {
        'type': 'adjacency',
        'positions': [pos1, pos2],
        'genes': [gene1_id, gene2_id],
        'gene_inversions': 2
    }
    
    return updated_sequence, inversion_record


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
    
    print(f"  FLIP BEFORE: {[(gene_id, move) for gene_id, _, move, _ in segment]}")
    
    # Step 1: Create the reversed segment
    reversed_segment = []
    original_positions = [pos for _, pos, _, _ in segment]
    
    # Reverse gene order but keep original positions
    for i, (gene_id, _, move, target) in enumerate(reversed(segment)):
        # Gene goes to the position at index i in the segment
        new_position = original_positions[i]
        reversed_segment.append((gene_id, new_position, move, target))
    
    # Step 2: Calculate movement updates based on positional distances
    final_segment = []
    for i, (gene_id, pos, old_move, target) in enumerate(reversed_segment):
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
        
        final_segment.append((gene_id, pos, new_move, target))
        
        print(f"    {gene_id}: moved {distance_moved} positions, {old_move} → {new_move}")
    
    # Step 3: Update the sequence
    updated_sequence[start_index:end_index + 1] = final_segment
    
    print(f"  FLIP AFTER:  {[(gene_id, move) for gene_id, _, move, _ in final_segment]}")
    
    # Create inversion record
    positions = [pos for _, pos, _, _ in segment]
    genes = [gene_id for gene_id, _, _, _ in segment]
    
    inversion_record = {
        'type': 'flip',
        'flip_indicator': flip_indicator,
        'segment_length': segment_length,
        'positions': positions,
        'genes': genes,
        'gene_inversions': len(segment)
    }
    
    return updated_sequence, inversion_record


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


def find_non_overlapping_flips(flip_patterns):
    """
    Find flip patterns that don't overlap in their index ranges to avoid conflicts.
    
    Args:
        flip_patterns: [(start_i, end_j, flip_indicator)] flip pattern tuples
        
    Returns:
        list: [(start_i, end_j, flip_indicator)] non-overlapping flips that can be applied safely
    """
    if not flip_patterns:
        return []
    
    # Sort by flip size (descending) to prioritize larger flips
    sorted_flips = sorted(flip_patterns, key=lambda x: x[2], reverse=True)
    
    non_overlapping = []
    used_ranges = set()
    
    for start_i, end_j, flip_indicator in sorted_flips:
        # Check if this flip range overlaps with any already used ranges
        overlap_found = False
        for used_start, used_end in used_ranges:
            # Check if ranges overlap: (start_i <= used_end) and (end_j >= used_start)
            if start_i <= used_end and end_j >= used_start:
                overlap_found = True
                break
        
        if not overlap_found:
            non_overlapping.append((start_i, end_j, flip_indicator))
            used_ranges.add((start_i, end_j))
    
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
    # Define recalculation function
    def recalculate_movements(current_sequence, target_positions):
        """Recalculate movements based on current vs target positions"""
        updated_sequence = []
        for gene_id, current_rank, old_movement, target_pos in current_sequence:
            target_rank = target_positions[gene_id]  # linearis rank
            new_movement = target_rank - current_rank  # target - current (correct direction)
            updated_sequence.append((gene_id, current_rank, new_movement, target_pos))
        return updated_sequence
    
    # Get target positions (linearis ranks) from original movement sequence
    target_positions = {gene_id: target_pos for gene_id, _, _, target_pos in movement_sequence}
    
    # DEBUG: Check initial state
    print("=== ALGORITHM DEBUG ===")
    print(f"Initial sequence length: {len(movement_sequence)}")
    total_movement = sum(abs(move) for _, _, move, _ in movement_sequence)
    print(f"Initial total movement: {total_movement}")
    non_zero = sum(1 for _, _, move, _ in movement_sequence if move != 0)
    print(f"Non-zero movements: {non_zero}")
    large_movements = [move for _, _, move, _ in movement_sequence if abs(move) > 10]
    print(f"Large movements (>10): {len(large_movements)} - sample: {large_movements[:10]}")
    
    # DEBUG: Check sample movements
    print("Sample movements:")
    for i in range(min(20, len(movement_sequence))):
        gene_id, pos, movement, target_pos = movement_sequence[i]
        print(f"  {gene_id}: source={pos}, target={target_pos}, movement={movement}")
    
    # DEBUG: Check if movements are meaningful
    if total_movement == 0:
        print("CRITICAL ERROR: No initial movement detected!")
        print("This means the algorithm has no work to do.")
        return {
            'final_sequence': movement_sequence,
            'inversion_events': [],
            'total_events': 0,
            'total_gene_inversions': 0,
            'iterations': 0,
            'converged': True,
            'final_total_movement': 0,
            'final_non_zero_movements': 0,
            'final_large_movements': 0
        }
    
    print("========================")
    
    print(f"Starting optimized iterative detection with {len(movement_sequence)} genes")
    
    current_sequence = movement_sequence.copy()
    inversion_events = []
    iteration = 0
    
    # Track previous states to detect cycles
    previous_states = []
    max_cycle_length = 10  # Maximum cycle length to detect
    
    while iteration < max_iterations:
        iteration += 1
        print(f"Iteration {iteration}: Starting...")
        applied_inversion = False
        batch_count = 0
        
        # Enhanced convergence monitoring
        current_total_movement = calculate_total_movement(current_sequence)
        non_zero_movements = [move for _, _, move, _ in current_sequence if move != 0]
        large_movements = sum(1 for _, _, move, _ in current_sequence if abs(move) > 2.0)
        
        # Calculate sum of positive and negative movements separately
        positive_movements = [move for _, _, move, _ in current_sequence if move > 0]
        negative_movements = [move for _, _, move, _ in current_sequence if move < 0]
        sum_positive = sum(positive_movements)
        sum_negative = abs(sum(negative_movements))
        
        print(f"  CONVERGENCE STATUS:")
        print(f"    Total movement: {current_total_movement:.2f}")
        print(f"    Non-zero movements: {len(non_zero_movements)} (sample: {non_zero_movements[:5] if non_zero_movements else []})")
        print(f"    Large movements (|>2|): {large_movements}")
        print(f"    Sum of + movements: {sum_positive:.2f} ({len(positive_movements)} genes)")
        print(f"    Sum of - movements: {sum_negative:.2f} ({len(negative_movements)} genes)")
        print(f"    Movement balance: {sum_positive - sum_negative:.2f}")
        
        # Track movement reduction
        if iteration > 1:
            movement_reduction = previous_total_movement - current_total_movement
            print(f"    Movement reduction: {movement_reduction:.2f}")
            if movement_reduction < 0.1:  # Very small reduction
                print(f"    WARNING: Minimal movement reduction - may be converging")
        
        previous_total_movement = current_total_movement
        
        # Check for flip patterns first (more efficient)
        print(f"  Checking flip patterns...")
        flip_patterns = detect_extended(current_sequence)
        print(f"  Found {len(flip_patterns)} flip patterns")
            
        # DEBUG: Check adjacency patterns too
        adjacency_patterns = detect_adjacency_inversions(current_sequence)
        print(f"  Found {len(adjacency_patterns)} adjacency patterns")
        
        if flip_patterns:
            
            # Find non-overlapping flips for batch processing
            non_overlapping_flips = find_non_overlapping_flips(flip_patterns)
            print(f"  Non-overlapping flips: {len(non_overlapping_flips)}")
            
            # Calculate current total movement before applying inversions
            current_total_movement = calculate_total_movement(current_sequence)
            
            # Test each flip and only apply those that reduce total movement
            valid_flips = []
            print(f"  DEBUG: Testing {len(non_overlapping_flips)} flips for biological validity...")
            for start_idx, end_idx, flip_size in non_overlapping_flips:
                # Create a test sequence with this flip applied
                test_sequence = current_sequence.copy()
                test_sequence, _ = apply_flip_inversion(test_sequence, start_idx, end_idx, flip_size)
                
                # Check if this flip reduces total movement
                new_total_movement = calculate_total_movement(test_sequence)
                
                # Strict acceptance by default, flexible if configured
                movement_change = new_total_movement - current_total_movement
                flexible_threshold = CONFIG.get('flexible_threshold', 0.0)  # Default to strict
                
                if movement_change <= flexible_threshold:
                    valid_flips.append((start_idx, end_idx, flip_size))
                    if movement_change <= 0:
                        print(f"    Flip accepted: {current_total_movement} -> {new_total_movement}")
                    else:
                        print(f"    Flip accepted (flexible): {current_total_movement} -> {new_total_movement} (+{movement_change:.2f})")
                else:
                    # Get sample genes from the flip to show what's being rejected
                    sample_genes = [current_sequence[i][0] for i in range(start_idx, min(start_idx + 3, end_idx + 1))]
                    sample_movements = [current_sequence[i][2] for i in range(start_idx, min(start_idx + 3, end_idx + 1))]
                    print(f"    Flip rejected: {current_total_movement} -> {new_total_movement} (increase: +{movement_change:.2f})")
                    print(f"      Sample genes: {sample_genes} with movements: {sample_movements}")
            
            # Apply all valid flips
            if valid_flips:
                for batch_idx, (start_idx, end_idx, flip_size) in enumerate(valid_flips):
                    genes_in_flip = [current_sequence[i][0] for i in range(start_idx, end_idx + 1)]
                    print(f"    Batch {batch_idx}: Flip {genes_in_flip[0]}...{genes_in_flip[-1]} (size: {flip_size})")
                    
            current_sequence, inversion_record = apply_flip_inversion(
                current_sequence, start_idx, end_idx, flip_size
            )
            inversion_record['iteration'] = iteration
            inversion_record['batch_position'] = batch_idx
            inversion_record['total_in_batch'] = len(valid_flips)
            inversion_events.append(inversion_record)
            batch_count += 1
            applied_inversion = True
            print(f"  Applied {batch_count} biologically valid flip inversions")
        
        # Process adjacencies only if no flips were applied
        if not applied_inversion:
            print(f"  Checking adjacency patterns...")
            adjacency_inversions = detect_adjacency_inversions(current_sequence)
            
            if adjacency_inversions:
                print(f"  Found {len(adjacency_inversions)} adjacency patterns")
                
                # Find non-overlapping adjacencies for batch processing
                non_overlapping = find_non_overlapping_adjacencies(adjacency_inversions)
                print(f"  Non-overlapping adjacencies: {len(non_overlapping)}")
                
                # Calculate current total movement
                current_total_movement = calculate_total_movement(current_sequence)
                
                # Test each adjacency and only apply those that reduce total movement
                valid_adjacencies = []
                for index1, index2 in non_overlapping:
                    # Create a test sequence with this adjacency applied
                    test_sequence = current_sequence.copy()
                    test_sequence, _ = apply_adjacency_inversion(test_sequence, index1, index2)
                    
                    # Check if this adjacency reduces total movement
                    new_total_movement = calculate_total_movement(test_sequence)
                    
                    # Strict acceptance for adjacencies by default, flexible if configured
                    movement_change = new_total_movement - current_total_movement
                    flexible_threshold = CONFIG.get('flexible_threshold', 0.0)  # Default to strict
                    
                    if movement_change <= flexible_threshold:
                        valid_adjacencies.append((index1, index2))
                        if movement_change <= 0:
                            print(f"    Adjacency accepted: {current_total_movement} -> {new_total_movement}")
                        else:
                            print(f"    Adjacency accepted (flexible): {current_total_movement} -> {new_total_movement} (+{movement_change:.2f})")
                    else:
                        gene1_id, _, move1 = current_sequence[index1]
                        gene2_id, _, move2 = current_sequence[index2]
                        print(f"    Adjacency rejected: {gene1_id}({move1})<->{gene2_id}({move2}) increase: +{movement_change:.2f}")
                        print(f"      Pattern: [{move1}, {move2}] - should be opposite signs for adjacency inversion")
                
                # Apply all valid adjacencies
                if valid_adjacencies:
                    for batch_idx, (index1, index2) in enumerate(valid_adjacencies):
                        gene1_id, pos1, move1 = current_sequence[index1]
                        gene2_id, pos2, move2 = current_sequence[index2]
                        print(f"    Batch {batch_idx}: {gene1_id}({move1}) <-> {gene2_id}({move2})")
                        
                    current_sequence, inversion_record = apply_adjacency_inversion(
                        current_sequence, index1, index2
                    )
                    
                    inversion_record['iteration'] = iteration
                    inversion_record['batch_position'] = batch_idx
                    inversion_record['total_in_batch'] = len(valid_adjacencies)
                    inversion_events.append(inversion_record)
                    batch_count += 1
                
                applied_inversion = True
                print(f"  Applied {batch_count} biologically valid adjacency inversions")
        
        # Termination condition
        if not applied_inversion:
            print(f"  No inversions found - terminating")
            
            #TRANSPOSITION
            # Check for translocation patterns before terminating
            # print(f"  Checking for intrachromosomal translocations...")
            # translocation_events = detect_and_apply_translocations(current_sequence, iteration)
            
            # if translocation_events:
            #     inversion_events.extend(translocation_events)
            #     applied_inversion = True
            #     print(f"  Applied {len(translocation_events)} translocations")
            # else:
            #     print(f"  No beneficial translocations possible - terminating")
            break
        
        # Progress check: stop if very few large movements remain
        if iteration % 50 == 0:
            large_movements = sum(1 for _, _, move, _ in current_sequence if abs(move) > 2.0)
            print(f"  Progress check: {large_movements} genes with |movement| > 2.0")
            if large_movements < 50:  # Most work done
                print(f"  Early termination: most large movements resolved")
                break
        
        # Cycle detection: check if we're stuck in a loop
        current_state_hash = hash(str(current_sequence))
        if current_state_hash in previous_states:
            cycle_start = previous_states.index(current_state_hash)
            cycle_length = len(previous_states) - cycle_start
            print(f"  CYCLE DETECTED! State repeated after {cycle_length} iterations.")
            print(f"  Stopping to prevent infinite loop. Total inversions: {len(inversion_events)}")
            break
            
        # Store current state (keep only last max_cycle_length states)
        previous_states.append(current_state_hash)
        if len(previous_states) > max_cycle_length:
            previous_states.pop(0)
        
        # Determine how many inversions were applied
        if batch_count > 0:
            inversions_applied = batch_count
        else:
            inversions_applied = 1
            
        print(f"Iteration {iteration}: Applied {inversions_applied} inversions")
        
        # Recalculate movements after each iteration to maintain accuracy
        current_sequence = recalculate_movements(current_sequence, target_positions)
        print(f"  Movement recalculation completed")
    

    if iteration % 10 == 0:
        remaining_movements = sum(1 for _, _, move, _ in current_sequence if abs(move) > 0.1)
        print(f"  Progress: {remaining_movements} genes still need movement")
    
    # Calculate summary statistics
    total_events = len(inversion_events)
    total_gene_inversions = sum(event['gene_inversions'] for event in inversion_events)
    adjacency_events = sum(1 for event in inversion_events if event['type'] == 'adjacency')
    flip_events = sum(1 for event in inversion_events if event['type'] == 'flip')
    
    # Final convergence analysis
    final_total_movement = calculate_total_movement(current_sequence)
    final_non_zero = sum(1 for _, _, move, _ in current_sequence if move != 0)
    final_large_movements = sum(1 for _, _, move, _ in current_sequence if abs(move) > 2.0)
    
    # Calculate final sum of positive and negative movements separately
    final_positive_movements = [move for _, _, move, _ in current_sequence if move > 0]
    final_negative_movements = [move for _, _, move, _ in current_sequence if move < 0]
    final_sum_positive = sum(final_positive_movements)
    final_sum_negative = abs(sum(final_negative_movements))
    final_movement_balance = final_sum_positive - final_sum_negative
    
    print(f"Optimized iterative detection completed after {iteration} iterations")
    print(f"  Total events: {total_events}")
    print(f"  Total gene inversions: {total_gene_inversions}")
    print(f"  Adjacency events: {adjacency_events}, Flip events: {flip_events}")
    
    # Convergence summary
    print(f"  CONVERGENCE SUMMARY:")
    print(f"    Final total movement: {final_total_movement:.2f}")
    print(f"    Final non-zero movements: {final_non_zero}")
    print(f"    Final large movements (|>2|): {final_large_movements}")
    print(f"    Final sum of + movements: {final_sum_positive:.2f} ({len(final_positive_movements)} genes)")
    print(f"    Final sum of - movements: {final_sum_negative:.2f} ({len(final_negative_movements)} genes)")
    print(f"    Final movement balance: {final_movement_balance:.2f}")
    
    # Determine convergence type
    if iteration >= max_iterations:
        print(f"    CONVERGENCE TYPE: Max iterations reached ({max_iterations})")
        converged = False
    elif final_non_zero == 0:
        print(f"    CONVERGENCE TYPE: Perfect convergence (all movements = 0)")
        converged = True
    elif final_large_movements < 10:
        print(f"    CONVERGENCE TYPE: Good convergence (few large movements remaining)")
        converged = True
    else:
        print(f"    CONVERGENCE TYPE: Partial convergence (many movements remaining)")
        converged = False
    
    return {
        'final_sequence': current_sequence,
        'inversion_events': inversion_events,
        'iterations': iteration,
        'total_events': total_events,
        'total_gene_inversions': total_gene_inversions,
        'adjacency_events': adjacency_events,
        'flip_events': flip_events,
        'converged': converged,
        'final_total_movement': final_total_movement,
        'final_non_zero_movements': final_non_zero,
        'final_large_movements': final_large_movements
    }


def create_pairwise_movement_sequence_per_chromosome(genome1_df, genome2_df, config):
    """
    Create separate movement sequences per chromosome pair.
    
    Returns:
        dict: {chromosome_pair_name: movement_sequence}
    """
    use_positions = config.get('use_genomic_positions', False)
    common_genes = set(genome1_df['busco_id']) & set(genome2_df['busco_id'])
    
    if len(common_genes) == 0:
        raise ValueError("No common genes found between the two genomes")
    
    genome1_grouped = genome1_df.groupby('sequence')
    genome2_grouped = genome2_df.groupby('sequence')
    
    chromosome_sequences = {}
    
    for chr1, chr1_genes in genome1_grouped:
        chr1_busco_ids = set(chr1_genes['busco_id'])
        
        # Find best matching chromosome
        best_match_chr = None
        best_overlap = 0
        
        for chr2, chr2_genes in genome2_grouped:
            chr2_busco_ids = set(chr2_genes['busco_id'])
            overlap = len(chr1_busco_ids & chr2_busco_ids)
            
            if overlap > best_overlap:
                best_overlap = overlap
                best_match_chr = chr2
        
        print(f"DEBUG: Best match for {chr1}: {best_match_chr} with {best_overlap} genes")  # ADD THIS
        print(f"DEBUG: Threshold check: {best_overlap} > {config.get('shared_genes_threshold', 50)}? {best_overlap > config.get('shared_genes_threshold', 50)}")  # ADD THIS
    
        
        if best_match_chr and best_overlap > config.get('shared_genes_threshold', 50):  # Configurable threshold for chromosome pairing
            chr1_data = chr1_genes.sort_values('gene_start')
            chr2_data = genome2_grouped.get_group(best_match_chr).sort_values('gene_start')
            chr_common_genes = chr1_busco_ids & set(chr2_data['busco_id'])
            
            print(f"DEBUG: use_positions = {use_positions}")
            print(f"DEBUG: {chr1} has {len(chr1_data)} genes, {best_match_chr} has {len(chr2_data)} genes")
            print(f"DEBUG: Common genes: {len(chr_common_genes)}")
            
            #Filter data to common genes only
            chr1_common_data = chr1_data[chr1_data['busco_id'].isin(chr_common_genes)].sort_values('gene_start')
            chr2_common_data = chr2_data[chr2_data['busco_id'].isin(chr_common_genes)].sort_values('gene_start')
            
            print(f"DEBUG: After filtering - chr1: {len(chr1_common_data)} genes, chr2: {len(chr2_common_data)} genes")
            # Create sequential positions from filtered data
            
            if use_positions:
                # Use genomic coordinates (in kb)
                chr1_positions = {row['busco_id']: row['gene_start'] // 1000 
                                for _, row in chr1_data.iterrows()}
                chr2_positions = {row['busco_id']: row['gene_start'] // 1000 
                                for _, row in chr2_data.iterrows()}
            else:
                # Use ordinal ranks (current behaviour)
                # Create sequential positions for common genes only
                chr1_positions = {row['busco_id']: idx for idx, (_, row) in enumerate(chr1_common_data.iterrows())}
                chr2_positions = {row['busco_id']: idx for idx, (_, row) in enumerate(chr2_common_data.iterrows())}
                
            # Calculate movements for this chromosome pair
            chromosome_movement_sequence = []
            
       

            
            # DEBUG: Print ranking information
            print(f"DEBUG: {chr1} has {len(chr1_data)} genes, {best_match_chr} has {len(chr2_data)} genes")
            print(f"DEBUG: Common genes: {len(chr_common_genes)}")
            print(f"DEBUG: First 5 chr1_sequential_positions: {list(chr1_positions.items())[:5]}")
            print(f"DEBUG: First 5 chr2_sequential_positions: {list(chr2_positions.items())[:5]}")
            
            for gene_id in chr_common_genes:
                pos1 = chr1_positions[gene_id]  # source position
                pos2 = chr2_positions[gene_id]  # target position  
                movement = pos2 - pos1  # target - source (positive = move right)
                chromosome_movement_sequence.append((gene_id, pos1, movement, pos2))  # 4-tuple
            
            # Sort by position
            chromosome_movement_sequence.sort(key=lambda x: x[1])
            
            chromosome_pair_name = f"{chr1}_vs_{best_match_chr}"
            chromosome_sequences[chromosome_pair_name] = chromosome_movement_sequence
            
            print(f"Chromosome pair {chromosome_pair_name}: {len(chromosome_movement_sequence)} genes")
    
    return chromosome_sequences


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
    
    integrated_analysis = integrate_best_alternatives(
        standard_analysis, 
        alternative_pathways, 
        evaluated_steps
    )
        
    return {
        'standard_analysis': standard_analysis,
        'integrated_analysis': integrated_analysis,
        'evaluated_steps': evaluated_steps,
        'total_standard_bit_score': total_standard_bit_score,
        'problematic_genes': list(set(problematic_genes)),
        'alternative_pathways': alternative_pathways,
        'has_alternatives': len(alternative_pathways) > 0
    }


def integrate_best_alternatives(standard_analysis, alternative_pathways, evaluated_steps):
    """Replace failed inversions with best alternatives"""
    
    final_events = standard_analysis['inversion_events'].copy()
    
    # For each failed event, replace with best alternative
    for i, step in enumerate(evaluated_steps):
        if not step['probability_evaluation']['threshold_passed']:

            failed_genes = step['probability_evaluation']['failed_genes']
            
            # Replace with alternatives for these genes
            for gene_id in failed_genes:
                if gene_id in alternative_pathways and alternative_pathways[gene_id]:
                    best_alternative = alternative_pathways[gene_id][0] 
  
                    final_events[i] = best_alternative['pathway_steps'][0]  # Use best pathway
    
    return {
        'final_sequence': final_events, 
        'inversion_events': final_events,
        'total_events': len(final_events),
        'converged': True 
    }


def check_events_iteration(movement_results):
    """
    Run inversion detection separately for each chromosome
    """
    all_chromosome_results = {}
    combined_events = []
    
    for chromosome, gene_results in movement_results.items():
        chromosome_sequence = extract_chromosome_movement_sequence(gene_results)
        
        chromosome_analysis = iterative_detection(chromosome_sequence)
        
        all_chromosome_results[chromosome] = chromosome_analysis
        combined_events.extend(chromosome_analysis['inversion_events'])
    
    return {
        'inversion_events': combined_events,
        'chromosome_results': all_chromosome_results,
        'total_events': len(combined_events),
        'total_gene_inversions': sum(event['gene_inversions'] for event in combined_events),
        'adjacency_events': sum(1 for event in combined_events if event['type'] == 'adjacency'),
        'flip_events': sum(1 for event in combined_events if event['type'] == 'flip'),
        'converged': all(r['converged'] for r in all_chromosome_results.values()),
        'iterations': max(r.get('iterations', 0) for r in all_chromosome_results.values())
    }

def extract_chromosome_movement_sequence(gene_results):
    """Extract movement sequence for a single chromosome"""
    movement_sequence = []
    
    for gene_id, result in gene_results.items():
        if result['current_ranges']:
            current_position = result['current_ranges'][0][0]
            mean_movement = result['movement_analysis']['mean_movement']
            movement_sequence.append((gene_id, current_position, mean_movement))
    
    movement_sequence.sort(key=lambda x: x[1])
    return movement_sequence



if __name__ == "__main__":
    print("done")
