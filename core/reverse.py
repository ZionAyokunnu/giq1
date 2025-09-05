import pandas as pd
import numpy as np
from collections import defaultdict
import math
from typing import Dict

from .query_movement import extract_gene_distribution
from .transposition import detect_and_apply_translocations
from config.settings import CONFIG
from .support import apply_batch_with_sequential_fallback, find_non_overlapping_flips, apply_concurrent_batch_flips, apply_concurrent_batch_adjacencies, validate_segment_independence
from .rules import detect_odd_length_incremental, detect_even_length_incremental, detect_extended, detect_flip_in_pattern

from logging import getLogger

logger = getLogger(__name__)

FOCUS_GENES = [
    'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J'
]  # Genes to trace through the algorithm

def debug_focus_gene(message, **kwargs):
    """Centralized debug function for focus gene tracking"""
    print(f" FOCUS_GENE_DEBUG: {message}")
    for key, value in kwargs.items():
        print(f"   {key}: {value}")

def debug_focus_gene_in_sequence(sequence, context=""):
    """Debug focus genes in any sequence"""
    print(f" FOCUS_GENE_DEBUG: Checking sequence {context}")
    found_genes = []
    for i, (gene_id, pos, move, target) in enumerate(sequence):
        if gene_id in FOCUS_GENES:
            print(f"   Found focus gene {gene_id} at index {i}: position={pos}, movement={move}, target={target}")
            found_genes.append((gene_id, i, pos, move, target))
    if not found_genes:
        print(f"   Focus genes NOT FOUND in sequence {context}")
    return found_genes

def debug_focus_gene_operation(operation_type, gene_id, old_pos, new_pos, old_move, new_move, target, context=""):
    """Debug focus gene operations (flips, adjacencies, etc.)"""
    if gene_id in FOCUS_GENES:
        print(f"🔍 FOCUS GENE OPERATION: {operation_type}")
        print(f"   Gene: {gene_id}")
        print(f"   Position: {old_pos} → {new_pos}")
        print(f"   Movement: {old_move} → {new_move}")
        print(f"   Target: {target}")
        print(f"   Context: {context}")
        print(f"   Progress: {abs(new_move)}/{abs(target - new_pos)} units remaining")


def extract_movement_sequence(movement_results):
    """
    Extract ordered movement sequence from movement analysis results.
    
      Returns:
        list: [(gene_id, current_position, mean_movement)] ordered by current position

    """
    movement_sequence = []
    
    debug_focus_gene("Starting extract_movement_sequence", movement_results_keys=list(movement_results.keys()))
    
    for chromosome, gene_results in movement_results.items():
        debug_focus_gene(f"Processing chromosome {chromosome}", gene_count=len(gene_results))
        
        for gene_id, result in gene_results.items():
            if result['current_ranges']:
                current_position = result['current_ranges'][0][0] 
                mean_movement = result['movement_analysis']['mean_movement']
                movement_sequence.append((gene_id, current_position, mean_movement))
                
                if gene_id in FOCUS_GENES:
                    debug_focus_gene(f"Found focus gene {gene_id} in {chromosome}", 
                                   current_position=current_position, 
                                   mean_movement=mean_movement,
                                   result=result)
    
    
    movement_sequence.sort(key=lambda x: x[1])
    
    debug_focus_gene_in_sequence(movement_sequence, "final sorted")
    
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
    
    debug_focus_gene("Starting detect_adjacency_inversions", sequence_length=len(movement_sequence))
    debug_focus_gene_in_sequence(movement_sequence, "input")
    
    for i in range(len(movement_sequence) - 1):
        gene1_id, _, movement1, _ = movement_sequence[i]
        gene2_id, _, movement2, _ = movement_sequence[i + 1]
        
        # Check if focus gene is involved in this adjacency
        if gene1_id in FOCUS_GENES or gene2_id in FOCUS_GENES:
            debug_focus_gene(f"Focus gene adjacency check at index {i}", 
                           gene1=gene1_id, movement1=movement1,
                           gene2=gene2_id, movement2=movement2,
                           has_opposite_signs=((movement1 > 0 and movement2 < 0) or (movement1 < 0 and movement2 > 0)))
            
            # Check if this adjacency would help the focus gene
            if ((movement1 > 0 and movement2 < 0) or (movement1 < 0 and movement2 > 0)):
                print(f"🔍 FOCUS GENE ADJACENCY POTENTIAL: {gene1_id}({movement1}) <-> {gene2_id}({movement2})")
                if gene1_id in FOCUS_GENES:
                    print(f"   {gene1_id} would move from pos {i} to pos {i+1}")
                if gene2_id in FOCUS_GENES:
                    print(f"   {gene2_id} would move from pos {i+1} to pos {i}")
        
        if ((movement1 > 0 and movement2 < 0) or 
            (movement1 < 0 and movement2 > 0)):
            adjacency_inversions.append((i, i + 1))
            
            if gene1_id in FOCUS_GENES or gene2_id in FOCUS_GENES:
                debug_focus_gene(f"Focus gene adjacency FOUND at index {i}", 
                               gene1=gene1_id, movement1=movement1,
                               gene2=gene2_id, movement2=movement2)
                print(f"🔍 FOCUS GENE ADJACENCY SELECTED: {gene1_id}({movement1}) <-> {gene2_id}({movement2}) at indices {i}-{i+1}")
    
    debug_focus_gene("Finished detect_adjacency_inversions", 
                    total_adjacencies=len(adjacency_inversions))
    
    return adjacency_inversions


def calculate_total_movement(sequence):
    """Calculate total absolute movement for a sequence."""
    return sum(abs(move) for _, _, move, _ in sequence)







def apply_adjacency_inversion(movement_sequence, index1, index2):
    """Apply adjacency inversion using working script logic."""
    updated_sequence = movement_sequence.copy()
    
    # Extract the 2-gene segment (always length 2 for adjacency)
    segment = [updated_sequence[index1], updated_sequence[index2]]
    start_index = index1
    segment_length = 2
    
    # Apply the working script's core logic
    for i in range(segment_length // 2):  # range(1) = [0]
        left_idx = i  # 0
        right_idx = segment_length - 1 - i  # 1
        
        # Get genes from segment
        left_gene_id, left_pos, left_move, left_target = segment[left_idx]
        right_gene_id, right_pos, right_move, right_target = segment[right_idx]
        
        # Calculate new positions (CRUCIAL: use start_index + idx pattern)
        new_left_pos = start_index + left_idx  # index1 + 0 = index1
        new_right_pos = start_index + right_idx  # index1 + 1 = index2
        
        # CRITICAL SWAP: right gene goes to left position, left gene goes to right position
        segment[left_idx] = (right_gene_id, new_left_pos, right_move, right_target)
        segment[right_idx] = (left_gene_id, new_right_pos, left_move, left_target)
    
    # Put segment back into sequence
    updated_sequence[index1] = segment[0]
    updated_sequence[index2] = segment[1]
    
    # MISSING PIECE: Immediate movement recalculation within adjacency
    print(f"    🔧 RECALCULATING movements within adjacency:")
    for i, (gene_id, new_pos, old_move, target_pos) in enumerate(segment):
        new_movement = target_pos - new_pos
        segment[i] = (gene_id, new_pos, new_movement, target_pos)
        print(f"      {gene_id}: {old_move} → {new_movement}")
        
        # Debug focus genes
        if gene_id in FOCUS_GENES:
            debug_focus_gene_operation("ADJACENCY", gene_id, new_pos, new_pos, old_move, new_movement, target_pos, f"adjacency[{index1}-{index2}]")
    
    # Put updated segment back into sequence
    updated_sequence[index1] = segment[0]
    updated_sequence[index2] = segment[1]
    
    # Create record
    inversion_record = {
        'type': 'adjacency',
        'positions': [left_pos, right_pos],
        'genes': [left_gene_id, right_gene_id],
        'gene_inversions': 2
    }
    
    return updated_sequence, inversion_record


def apply_flip_inversion(movement_sequence, start_index, end_index, flip_indicator):
    """Apply flip inversion using working script logic."""
    updated_sequence = movement_sequence.copy()
    
    # Extract the segment to flip
    segment = updated_sequence[start_index:end_index + 1]
    segment_length = len(segment)
    
    print(f"  FLIP BEFORE: {[(gene_id, move) for gene_id, _, move, _ in segment]}")
    
    # DEBUG: Show the actual sequence positions and genes
    print(f"  🔍 SEQUENCE DEBUG - Segment [{start_index}-{end_index}]:")
    print(f"     Segment length: {segment_length}")
    print(f"     All genes in segment:")
    for i in range(start_index, end_index + 1):
        if i < len(movement_sequence):
            gene_id, pos, move, target = movement_sequence[i]
            focus_marker = " 🔍" if gene_id in FOCUS_GENES else ""
            print(f"       [{i}]: {gene_id} at pos={pos}, move={move}, target={target}{focus_marker}")
        else:
            print(f"       [{i}]: INDEX OUT OF RANGE!")
    
    # DEBUG: Check for focus genes in this segment
    focus_genes_in_segment = [g for g, _, _, _ in segment if g in FOCUS_GENES]
    if focus_genes_in_segment:
        print(f"     🔍 FOCUS GENES in this segment: {focus_genes_in_segment}")
    
    # DEBUG: Calculate magnitude before flip
    magnitude_before = sum(abs(move) for _, _, move, _ in segment)
    print(f"  🔍 MAGNITUDE ANALYSIS - Before flip:")
    print(f"     Segment [{start_index}-{end_index}]: {segment_length} genes")
    print(f"     Magnitude before: {magnitude_before}")
    print(f"     Individual movements: {[(gene_id, move) for gene_id, _, move, _ in segment]}")
    
    # Apply the working script's EXACT core logic
    for i in range(segment_length // 2):
        left_idx = i
        right_idx = segment_length - 1 - i
        
        # Get genes from segment  
        left_gene_id, left_pos, left_move, left_target = segment[left_idx]
        right_gene_id, right_pos, right_move, right_target = segment[right_idx]
        
        # Calculate new positions (CRUCIAL: use start_index + idx pattern)
        new_left_pos = start_index + left_idx
        new_right_pos = start_index + right_idx
        
        # CRITICAL SWAP: Apply the working script's exact swapping
        segment[left_idx] = (right_gene_id, new_left_pos, right_move, right_target)
        segment[right_idx] = (left_gene_id, new_right_pos, left_move, left_target)
    
    # IMPORTANT: Handle fulcrum for odd-length segments (from EDCBA example)
    # The middle gene stays in place - this is automatic in the loop logic
    
    # Put segment back into sequence
    updated_sequence[start_index:end_index + 1] = segment
    
    # MISSING PIECE: Immediate movement recalculation within flip
    print(f"    🔧 RECALCULATING movements within flip:")
    
    # DEBUG: Calculate magnitude after flip
    magnitude_after = 0
    print(f"  🔍 MAGNITUDE ANALYSIS - After flip:")
    print(f"     Individual movements after flip:")
    
    for i, (gene_id, new_pos, old_move, target_pos) in enumerate(segment):
        new_movement = target_pos - new_pos
        segment[i] = (gene_id, new_pos, new_movement, target_pos)
        magnitude_after += abs(new_movement)
        print(f"      {gene_id}: {old_move} → {new_movement} (|{new_movement}| = {abs(new_movement)})")
        
        # Debug focus genes
        if gene_id in FOCUS_GENES:
            debug_focus_gene_operation("FLIP", gene_id, new_pos, new_pos, old_move, new_movement, target_pos, f"flip[{start_index}-{end_index}]")
    
    print(f"     Magnitude after: {magnitude_after}")
    print(f"     Magnitude change: {magnitude_after - magnitude_before}")
    
    # Put updated segment back into sequence
    updated_sequence[start_index:end_index + 1] = segment
    
    print(f"  FLIP AFTER:  {[(gene_id, move) for gene_id, _, move, _ in updated_sequence[start_index:end_index + 1]]}")
    
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



def find_non_overlapping_adjacencies(adjacency_inversions, current_sequence):
    """
    Find adjacencies that are truly adjacent AND don't overlap.
    
    Args:
        adjacency_inversions: [(index1, index2)] list of adjacent pairs
        current_sequence: Current movement sequence to validate against
        
    Returns:
        list: [(index1, index2)] non-overlapping pairs that are truly adjacent
    """
    non_overlapping = []
    used_indices = set()
    
    debug_focus_gene("Starting find_non_overlapping_adjacencies", 
                    total_adjacencies=len(adjacency_inversions),
                    sequence_length=len(current_sequence))
    
    # DEBUG: Print all input adjacencies
    logger.debug(f"  DEBUG: Input adjacencies: {adjacency_inversions}")
    
    # Sort by first index to process left-to-right
    sorted_adjacencies = sorted(adjacency_inversions, key=lambda x: x[0])
    logger.debug(f"  DEBUG: Sorted adjacencies: {sorted_adjacencies}")
    
    for index1, index2 in sorted_adjacencies:
        logger.debug(f"  DEBUG: Processing adjacency {index1}-{index2}")
        
        # CRITICAL: Verify they're still adjacent (consecutive indices)
        if abs(index2 - index1) != 1:
            logger.debug(f"  DEBUG: Adjacency {index1}-{index2} REJECTED - not consecutive indices (diff={abs(index2 - index1)})")
            debug_focus_gene(f"Adjacency {index1}-{index2} REJECTED - not consecutive indices")
            continue  # Skip non-adjacent pairs
            
        # Check if focus gene is involved in this adjacency
        focus_gene_involved = False
        try:
            gene1_id = current_sequence[index1][0]
            gene2_id = current_sequence[index2][0]
            focus_gene_involved = gene1_id in FOCUS_GENES or gene2_id in FOCUS_GENES
            logger.debug(f"  DEBUG: Genes at {index1}-{index2}: {gene1_id} <-> {gene2_id}")
            print(f"  DEBUG: Focus gene involved: {focus_gene_involved}")
        except IndexError:
            logger.debug(f"  DEBUG: Adjacency {index1}-{index2} REJECTED - index out of bounds")
            debug_focus_gene(f"Adjacency {index1}-{index2} REJECTED - index out of bounds")
            continue
        
        if focus_gene_involved:
            debug_focus_gene(f"Focus gene adjacency in non_overlapping check", 
                           index1=index1, index2=index2,
                           gene1=gene1_id, gene2=gene2_id,
                           index1_used=index1 in used_indices,
                           index2_used=index2 in used_indices,
                           will_include=index1 not in used_indices and index2 not in used_indices)
        
        # Check if either index is already used
        index1_used = index1 in used_indices
        index2_used = index2 in used_indices
        logger.debug(f"  DEBUG: Index usage check - index1({index1}): {index1_used}, index2({index2}): {index2_used}")
        
        if index1 not in used_indices and index2 not in used_indices:
            logger.debug(f"  DEBUG: Both indices available, proceeding with validation")
            
            # CRITICAL: Double-check the genes are consecutive in current sequence
            gene1_id = current_sequence[index1][0]
            gene2_id = current_sequence[index2][0]
            
            # Verify this is still a valid adjacency pattern
            gene1_movement = current_sequence[index1][2]
            gene2_movement = current_sequence[index2][2]
            
            logger.debug(f"  DEBUG: Movement validation - {gene1_id}({gene1_movement}) <-> {gene2_id}({gene2_movement})")
            
            # Check if movements have opposite signs (required for adjacency inversion)
            has_opposite_signs = ((gene1_movement > 0 and gene2_movement < 0) or 
                                 (gene1_movement < 0 and gene2_movement > 0))
            
            logger.debug(f"  DEBUG: Opposite signs check: {has_opposite_signs}")
            
            if not has_opposite_signs:
                print(f"  DEBUG: Adjacency {index1}-{index2} REJECTED - no opposite signs")
                debug_focus_gene(f"Adjacency {index1}-{index2} REJECTED - no opposite signs: {gene1_id}({gene1_movement}) <-> {gene2_id}({gene2_movement})")
                continue
            
            non_overlapping.append((index1, index2))
            used_indices.add(index1)
            used_indices.add(index2)
            
            logger.debug(f"  DEBUG: Adjacency {index1}-{index2} ACCEPTED - added to non_overlapping")
            print(f"  DEBUG: Updated used_indices: {used_indices}")
            
            if focus_gene_involved:
                debug_focus_gene(f"Focus gene adjacency INCLUDED in non_overlapping", 
                               index1=index1, index2=index2,
                               gene1=gene1_id, gene2=gene2_id,
                               gene1_movement=gene1_movement,
                               gene2_movement=gene2_movement)
        else:
            logger.debug(f"  DEBUG: Adjacency {index1}-{index2} REJECTED - index already used")
            if focus_gene_involved:
                debug_focus_gene(f"Focus gene adjacency EXCLUDED from non_overlapping", 
                               index1=index1, index2=index2,
                               reason="index already used")
    
    print(f"  DEBUG: Final non_overlapping result: {non_overlapping}")
    debug_focus_gene("Finished find_non_overlapping_adjacencies", 
                    non_overlapping_count=len(non_overlapping))
    
    return non_overlapping




def recalculate_movements(current_sequence, target_positions):
    """Recalculate movements based on current vs target positions"""
    print(f"  DEBUG: recalculate_movements called")
    
    debug_focus_gene("Starting recalculate_movements", 
                    sequence_length=len(current_sequence),
                    target_positions_keys=list(target_positions.keys())[:5])
    
    updated_sequence = []
    for gene_id, current_rank, old_movement, target_pos in current_sequence:
        target_rank = target_positions[gene_id]  # linearis rank
        new_movement = target_rank - current_rank  # target - current (correct direction)
        updated_sequence.append((gene_id, current_rank, new_movement, target_pos))
        
        # DEBUG: Track focus genes in recalculate_movements
        if gene_id in FOCUS_GENES:
            print(f"  DEBUG {gene_id} - recalculate_movements:")
            print(f"    current_rank: {current_rank}")
            print(f"    target_rank: {target_rank}")
            print(f"    old_movement: {old_movement}")
            print(f"    new_movement: {new_movement}")
            print(f"    target_positions[{gene_id}]: {target_positions.get(gene_id, 'NOT_FOUND')}")
            
            debug_focus_gene(f"Focus gene {gene_id} movement recalculation", 
                           current_rank=current_rank,
                           target_rank=target_rank,
                           old_movement=old_movement,
                           new_movement=new_movement,
                           movement_change=new_movement - old_movement)
    
    debug_focus_gene_in_sequence(updated_sequence, "after recalculation")
    
    return updated_sequence
    

def iterative_detection(movement_sequence, max_iterations=1000):
    """
    Optimized iterative detection with batch processing of non-overlapping inversions.
    
    Args:
        movement_sequence: [(gene_id, position, movement)]
        max_iterations: Maximum number of iterations
        
    Returns:
        dict: Complete inversion analysis results
    """
    # DEBUG: Show the complete movement sequence at the start
    print(f"🔍 DEBUG: iterative_detection called with {len(movement_sequence)} genes")
    print(f"🔍 DEBUG: First 10 genes in movement_sequence:")
    for i, (gene_id, pos, move, target) in enumerate(movement_sequence[:10]):
        print(f"  [{i}]: {gene_id} at pos={pos}, move={move}, target={target}")
    
    # Check for focus genes in the input
    focus_genes_found = [g for g, _, _, _ in movement_sequence if g in FOCUS_GENES]
    if focus_genes_found:
        print(f"🔍 DEBUG: Focus genes found in input: {focus_genes_found}")
        for gene_id, pos, move, target in movement_sequence:
            if gene_id in FOCUS_GENES:
                print(f"  {gene_id}: pos={pos}, move={move}, target={target}")
    else:
        print(f"🔍 DEBUG: NO focus genes found in input movement_sequence")
    
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
    for i in range(min(3, len(movement_sequence))):
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
        # DEBUG: Track focus genes through each iteration
        for gene_id, pos, movement, target_pos in current_sequence:
            if gene_id in FOCUS_GENES:
                print(f"  DEBUG {gene_id} - Iteration {iteration}:")
                print(f"    Current position: {pos}")
                print(f"    Current movement: {movement}")
                print(f"    Target position: {target_pos}")
                print(f"    Distance to target: {target_pos - pos}")
                print(f"    Target_positions[{gene_id}]: {target_positions.get(gene_id, 'NOT_FOUND')}")
                
                debug_focus_gene(f"Iteration {iteration} start", 
                               position=pos, movement=movement, target=target_pos,
                               distance_to_target=target_pos - pos)
                
                # Track progress towards convergence
                if movement == 0:
                    print(f"    ✅ CONVERGED!")
                elif abs(movement) < abs(target_pos - pos):
                    print(f"    📈 IMPROVING (movement reduced)")
                else:
                    print(f"    📉 DEGRADING (movement increased)")
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
            # non_overlapping_flips = find_non_overlapping_flips_with_constraint_sorting(
            #     flip_patterns, current_sequence, target_positions
            # )
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
                movement_change = new_total_movement - current_total_movement
                
                # Strict acceptance by default, flexible if configured
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
                # 🔍 VALIDATE INDEPENDENCE
                validate_segment_independence(current_sequence, valid_flips, 'flip')
                
                # Apply with sequential fallback
                current_sequence, flip_records = apply_batch_with_sequential_fallback(
                    current_sequence, valid_flips, 'flip', target_positions
                )
                
                # Add iteration info to records
                for record in flip_records:
                    record['iteration'] = iteration
                    inversion_events.append(record)
                
                applied_inversion = True
                print(f"  Applied {len(flip_records)} flip inversions")
        
        # Process adjacencies only if no flips were applied
        if not applied_inversion:
            print(f"  Checking adjacency patterns...")
            adjacency_inversions = detect_adjacency_inversions(current_sequence)
            
            if adjacency_inversions:
                print(f"  Found {len(adjacency_inversions)} adjacency patterns")
                
                # Find non-overlapping adjacencies for batch processing
                non_overlapping = find_non_overlapping_adjacencies(adjacency_inversions, current_sequence)
                print(f"  Non-overlapping adjacencies: {len(non_overlapping)}")
                
                # Calculate current total movement
                current_total_movement = calculate_total_movement(current_sequence)
                
                # Test each adjacency and only apply those that reduce total movement
                valid_adjacencies = []
                debug_focus_gene(f"Iteration {iteration} - Starting adjacency evaluation", 
                               total_adjacencies=len(non_overlapping))
                
                for index1, index2 in non_overlapping:
                    # Create a test sequence with this adjacency applied
                    test_sequence = current_sequence.copy()
                    test_sequence, _ = apply_adjacency_inversion(test_sequence, index1, index2)
                    
                    # Check if this adjacency reduces total movement
                    new_total_movement = calculate_total_movement(test_sequence)
                    
                    # Strict acceptance for adjacencies by default, flexible if configured
                    movement_change = new_total_movement - current_total_movement
                    flexible_threshold = CONFIG.get('flexible_threshold', 0.0)  # Default to strict
                    
                    # Check if focus gene is involved in this adjacency
                    gene1_id, _, _, _ = current_sequence[index1]
                    gene2_id, _, _, _ = current_sequence[index2]
                    focus_gene_involved = gene1_id in FOCUS_GENES or gene2_id in FOCUS_GENES
                    
                    # Debug ALL adjacencies being evaluated
                    debug_focus_gene(f"Iteration {iteration} - Adjacency evaluation loop", 
                                   index1=index1, index2=index2,
                                   gene1=gene1_id, gene2=gene2_id,
                                   focus_gene_involved=focus_gene_involved,
                                   current_movement=current_total_movement,
                                   new_movement=new_total_movement,
                                   movement_change=movement_change)
                    
                    print(f" DEBUG: Evaluating adjacency {index1}-{index2}: {gene1_id}({current_sequence[index1][2]}) <-> {gene2_id}({current_sequence[index2][2]})")
                    
                    if focus_gene_involved:
                        debug_focus_gene(f"Iteration {iteration} - Focus gene adjacency evaluation", 
                                       gene1=gene1_id, gene2=gene2_id,
                                       index1=index1, index2=index2,
                                       current_movement=current_total_movement,
                                       new_movement=new_total_movement,
                                       movement_change=movement_change,
                                       threshold=flexible_threshold,
                                       will_accept=movement_change <= flexible_threshold)
                    
                    if focus_gene_involved:
                        debug_focus_gene(f"Iteration {iteration} - Focus gene adjacency evaluation", 
                                       gene1=gene1_id, gene2=gene2_id,
                                       index1=index1, index2=index2,
                                       current_movement=current_total_movement,
                                       new_movement=new_total_movement,
                                       movement_change=movement_change,
                                       threshold=flexible_threshold,
                                       will_accept=movement_change <= flexible_threshold)
                    
                    if movement_change <= flexible_threshold:
                        valid_adjacencies.append((index1, index2))
                        if movement_change <= 0:
                            print(f"    Adjacency accepted: {current_total_movement} -> {new_total_movement}")
                        else:
                            print(f"    Adjacency accepted (flexible): {current_total_movement} -> {new_total_movement} (+{movement_change:.2f})")
                    else:
                        gene1_id, _, move1, _ = current_sequence[index1]
                        gene2_id, _, move2, _ = current_sequence[index2]
                        print(f"    Adjacency rejected: {gene1_id}({move1})<->{gene2_id}({move2}) increase: +{movement_change:.2f}")
                        print(f"      Pattern: [{move1}, {move2}] - should be opposite signs for adjacency inversion")
                        
                        if focus_gene_involved:
                            debug_focus_gene(f"Iteration {iteration} - Focus gene adjacency REJECTED", 
                                           gene1=gene1_id, gene2=gene2_id,
                                           move1=move1, move2=move2,
                                           movement_change=movement_change,
                                           threshold=flexible_threshold)
                
                debug_focus_gene(f"Iteration {iteration} - Finished adjacency evaluation", 
                               valid_adjacencies=len(valid_adjacencies))
                
                # Apply all valid adjacencies
                if valid_adjacencies:
                    # 🔍 VALIDATE INDEPENDENCE
                    validate_segment_independence(current_sequence, valid_adjacencies, 'adjacency')
                    
                    # Apply with sequential fallback
                    current_sequence, adjacency_records = apply_batch_with_sequential_fallback(
                        current_sequence, valid_adjacencies, 'adjacency', target_positions
                    )
                    
                    # Add iteration info to records
                    for record in adjacency_records:
                        record['iteration'] = iteration
                        inversion_events.append(record)
                    
                    applied_inversion = True
                    print(f"  Applied {len(adjacency_records)} adjacency inversions")
        
        # Termination condition
        if not applied_inversion:
            print(f"  No inversions found - terminating")
    
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
            
            debug_focus_gene("Cycle detected - checking focus gene state", 
                           cycle_length=cycle_length,
                           total_inversions=len(inversion_events))
            debug_focus_gene_in_sequence(current_sequence, "at cycle detection")
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
        # DEBUG: Capture movement values before recalculation (distance tracking results)
        # Use a separate variable to avoid state interference
        debug_distance_tracking = None
        if batch_count > 0:
            debug_distance_tracking = {gene_id: move for gene_id, _, move, _ in current_sequence}
        
        # NOTE: Movement recalculation is already handled within apply functions
        # The apply_flip_inversion and apply_adjacency_inversion functions
        # already recalculate movements correctly after each operation
        # No additional recalculation needed here
        
        # DEBUG: Compare distance tracking vs recalculation for genes that were inverted
        if batch_count > 0 and debug_distance_tracking:  # Only if inversions were applied
            print(f"  DEBUG: Comparing distance tracking vs recalculation for {batch_count} inversions")
            for event in inversion_events[-batch_count:]:  # Get the events from this iteration
                if event['type'] == 'flip':
                    genes_in_flip = event['genes']
                    print(f"    Flip event genes: {genes_in_flip}")
                    # Compare the two approaches using separate lookups
                    for gene_id in genes_in_flip:
                        distance_tracked = debug_distance_tracking.get(gene_id, 'N/A')
                        # Find recalculated value without modifying current_sequence
                        recalculated = None
                        for gid, _, move, _ in current_sequence:
                            if gid == gene_id:
                                recalculated = move
                                break
                        if recalculated is None:
                            recalculated = 'N/A'
                        
                        if distance_tracked != 'N/A' and recalculated != 'N/A':
                            if abs(distance_tracked - recalculated) < 0.001:
                                print(f"      {gene_id}: EQUAL - Distance tracked: {distance_tracked}, Recalculated: {recalculated}")
                            else:
                                print(f"      {gene_id}: NOT EQUAL - Distance tracked: {distance_tracked}, Recalculated: {recalculated}")
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
    
    # Determine convergence type first
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
    
    # Final focus gene analysis
    debug_focus_gene("Algorithm completion - final state", 
                    total_iterations=iteration,
                    total_events=total_events,
                    final_total_movement=final_total_movement,
                    converged=converged)
    debug_focus_gene_in_sequence(current_sequence, "final state")
    
    # DEBUG: Track focus genes in final state
    print(f"\n🔍 DEBUG: Focus genes final state after iteration {iteration}:")
    focus_genes_final = [g for g, _, _, _ in current_sequence if g in FOCUS_GENES]
    if focus_genes_final:
        for gene_id, pos, movement, target_pos in current_sequence:
            if gene_id in FOCUS_GENES:
                print(f"  {gene_id}: pos={pos}, movement={movement}, target={target_pos}")
                if movement == 0:
                    print(f"    ✅ CONVERGED!")
                else:
                    print(f"    ❌ NOT CONVERGED - movement={movement}")
    else:
        print(f"  ❌ NO focus genes found in final state")
    
    # Convergence summary
    print(f"  CONVERGENCE SUMMARY:")
    print(f"    Final total movement: {final_total_movement:.2f}")
    print(f"    Final non-zero movements: {final_non_zero}")
    print(f"    Final large movements (|>2|): {final_large_movements}")
    print(f"    Final sum of + movements: {final_sum_positive:.2f} ({len(final_positive_movements)} genes)")
    print(f"    Final sum of - movements: {final_sum_negative:.2f} ({len(final_negative_movements)} genes)")
    print(f"    Final movement balance: {final_movement_balance:.2f}")
    
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
    
    # DEBUG: Check for focus genes in the input data
    focus_genes_in_input = [g for g in common_genes if g in FOCUS_GENES]
    if focus_genes_in_input:
        print(f"🔍 DEBUG: Focus genes found in input data: {focus_genes_in_input}")
        for gene_id in focus_genes_in_input:
            genome1_gene = genome1_df[genome1_df['busco_id'] == gene_id]
            genome2_gene = genome2_df[genome2_df['busco_id'] == gene_id]
            if not genome1_gene.empty and not genome2_gene.empty:
                print(f"  {gene_id}: genome1_chr={genome1_gene.iloc[0]['sequence']}, genome2_chr={genome2_gene.iloc[0]['sequence']}")
    else:
        print(f"🔍 DEBUG: NO focus genes found in input data")
    
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
            
            # DEBUG: Check complete similarities at each stage
            # print(f"DEBUG: chr1_common_data genes: {sorted(chr1_common_data['busco_id'].tolist())}")
            # print(f"DEBUG: chr2_common_data genes: {sorted(chr2_common_data['busco_id'].tolist())}")
            # print(f"DEBUG: chr_common_genes: {sorted(chr_common_genes)}")
            print(f"DEBUG: Length check - chr1_common: {len(chr1_common_data)}, chr2_common: {len(chr2_common_data)}, common_genes: {len(chr_common_genes)}")
            
            # Check if gene orders are identical after sorting
            chr1_gene_order = chr1_common_data['busco_id'].tolist()
            chr2_gene_order = chr2_common_data['busco_id'].tolist()
            print(f"DEBUG: Gene order identical? {chr1_gene_order == chr2_gene_order}")
            if chr1_gene_order != chr2_gene_order:
                print(f"DEBUG: First 5 chr1 order: {chr1_gene_order[:5]}")
                print(f"DEBUG: First 5 chr2 order: {chr2_gene_order[:5]}")
            
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
                
            # DEBUG: Check position mappings
            print(f"DEBUG: chr1_positions length: {len(chr1_positions)}")
            print(f"DEBUG: chr2_positions length: {len(chr2_positions)}")
            # print(f"DEBUG: chr1_positions keys: {sorted(chr1_positions.keys())}")
            # print(f"DEBUG: chr2_positions keys: {sorted(chr2_positions.keys())}")
            
            # Check for any genes with same position but different movement
            for gene_id in chr_common_genes:
                pos1 = chr1_positions[gene_id]
                pos2 = chr2_positions[gene_id]
                # if pos1 == pos2:
                #     print(f"DEBUG: WARNING - {gene_id}: same position {pos1} in both genomes")
                
            # DEBUG: Check for focus genes in the movement sequence being created
            focus_genes_in_chr = [g for g in chr_common_genes if g in FOCUS_GENES]
            if focus_genes_in_chr:
                print(f"🔍 DEBUG: Focus genes found in {chr1}_vs_{best_match_chr}: {focus_genes_in_chr}")
                for gene_id in focus_genes_in_chr:
                    pos1 = chr1_positions[gene_id]
                    pos2 = chr2_positions[gene_id]
                    movement = pos2 - pos1
                    print(f"  {gene_id}: chr1_pos={pos1}, chr2_pos={pos2}, movement={movement}")
            else:
                print(f"🔍 DEBUG: NO focus genes found in {chr1}_vs_{best_match_chr}")
            
            # Calculate movements for this chromosome pair
            chromosome_movement_sequence = []
            
       

            
            # DEBUG: Print ranking information
            print(f"DEBUG: {chr1} has {len(chr1_data)} genes, {best_match_chr} has {len(chr2_data)} genes")
            print(f"DEBUG: Common genes: {len(chr_common_genes)}")
            # print(f"DEBUG: First 5 chr1_sequential_positions: {list(chr1_positions.items())[:5]}")
            # print(f"DEBUG: First 5 chr2_sequential_positions: {list(chr2_positions.items())[:5]}")
            
            for gene_id in chr_common_genes:
                pos1 = chr1_positions[gene_id]  # source position
                pos2 = chr2_positions[gene_id]  # target position  
                movement = pos2 - pos1  # target - source (positive = move right)
                # movement = pos1 - pos2  # source - target (positive = move left)
                chromosome_movement_sequence.append((gene_id, pos1, movement, pos2))  # 4-tuple
                
                # DEBUG: Add comprehensive tracing for 4164at7147
                if gene_id == '4164at7147':
                    print(f"FOUND 4164at7147 in movement sequence creation!")
                    print(f"=== TRACING 4164at7147 THROUGH PIPELINE ===")
                    print(f"Step 1 - Chromosome Pairing:")
                    print(f"  chr1: {chr1}, chr2: {best_match_chr}")
                    print(f"  overlap: {len(chr_common_genes)}")
                    print(f"  best_overlap: {best_overlap}")
                    print(f"  threshold: {config.get('shared_genes_threshold', 50)}")
                    
                    print(f"Step 2 - Gene Filtering:")
                    print(f"  chr1_common_data length: {len(chr1_common_data)}")
                    print(f"  chr2_common_data length: {len(chr2_common_data)}")
                    print(f"  chr_common_genes length: {len(chr_common_genes)}")
                    
                    # Find 4164at7147 in both datasets
                    chr1_4164 = chr1_common_data[chr1_common_data['busco_id'] == '4164at7147']
                    chr2_4164 = chr2_common_data[chr2_common_data['busco_id'] == '4164at7147']
                    
                    if not chr1_4164.empty:
                        print(f"  4164at7147 in chr1: gene_start={chr1_4164.iloc[0]['gene_start']}, gene_end={chr1_4164.iloc[0]['gene_end']}")
                    if not chr2_4164.empty:
                        print(f"  4164at7147 in chr2: gene_start={chr2_4164.iloc[0]['gene_start']}, gene_end={chr2_4164.iloc[0]['gene_end']}")
                    
                    print(f"Step 3 - Position Assignment:")
                    print(f"  chr1_positions keys: {list(chr1_positions.keys())[:10]}...")
                    print(f"  chr2_positions keys: {list(chr2_positions.keys())[:10]}...")
                    
                    if '4164at7147' in chr1_positions and '4164at7147' in chr2_positions:
                        pos1 = chr1_positions['4164at7147']
                        pos2 = chr2_positions['4164at7147']
                        movement = pos2 - pos1
                        print(f"  4164at7147: chr1_pos={pos1}, chr2_pos={pos2}, movement={movement}")
                        
                        # Show surrounding genes for context
                        chr1_genes = list(chr1_positions.keys())
                        chr2_genes = list(chr2_positions.keys())
                        
                        start_idx = max(0, pos1-2)
                        end_idx = min(len(chr1_genes), pos1+3)
                        print(f"  chr1 context around 4164at7147: {chr1_genes[start_idx:end_idx]}")
                        
                        start_idx = max(0, pos2-2)
                        end_idx = min(len(chr2_genes), pos2+3)
                        print(f"  chr2 context around 4164at7147: {chr2_genes[start_idx:end_idx]}")
                    
                    print(f"Step 4 - Movement Sequence Creation:")
                    print(f"  Final movement sequence for 4164at7147: {movement}")
                    print(f"  Target position for 4164at7147: {pos2}")
                    print("=" * 50)
            
            # Sort by position
            chromosome_movement_sequence.sort(key=lambda x: x[1])
            
            chromosome_pair_name = f"{chr1}_vs_{best_match_chr}"
            chromosome_sequences[chromosome_pair_name] = chromosome_movement_sequence
            
            print(f"Chromosome pair {chromosome_pair_name}: {len(chromosome_movement_sequence)} genes")
    
    return chromosome_sequences

