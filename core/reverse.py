import pandas as pd
import numpy as np
from collections import defaultdict
import math
from typing import Dict, List

from config.settings import CONFIG
from core.rules import (
    detect_adjacency_inversions,
    detect_extended,
)
from core.support import (
    validate_segment_independence,
    apply_batch_with_sequential_fallback,
    find_non_overlapping_adjacencies,
    find_non_overlapping_flips
)


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
    
    for i, (gene_id, new_pos, old_move, target_pos) in enumerate(segment):
        new_movement = target_pos - new_pos
        segment[i] = (gene_id, new_pos, new_movement, target_pos)

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
    
    for i in range(start_index, end_index + 1):
        if i < len(movement_sequence):
            gene_id, pos, move, target = movement_sequence[i]
   
        else:
            print(f"       [{i}]: INDEX OUT OF RANGE!")

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

    # The middle gene stays in place - this is automatic in the loop logic
    
    # Put segment back into sequence
    updated_sequence[start_index:end_index + 1] = segment

    for i, (gene_id, new_pos, old_move, target_pos) in enumerate(segment):
        new_movement = target_pos - new_pos
        segment[i] = (gene_id, new_pos, new_movement, target_pos)

    # Put updated segment back into sequence
    updated_sequence[start_index:end_index + 1] = segment
   
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
        
        if best_match_chr and best_overlap > config.get('shared_genes_threshold', 50):  # Configurable threshold for chromosome pairing
            chr1_data = chr1_genes.sort_values('gene_start')
            chr2_data = genome2_grouped.get_group(best_match_chr).sort_values('gene_start')
            chr_common_genes = chr1_busco_ids & set(chr2_data['busco_id'])

            #Filter data to common genes only
            chr1_common_data = chr1_data[chr1_data['busco_id'].isin(chr_common_genes)].sort_values('gene_start')
            chr2_common_data = chr2_data[chr2_data['busco_id'].isin(chr_common_genes)].sort_values('gene_start')
           
            # Check if gene orders are identical after sorting
            chr1_gene_order = chr1_common_data['busco_id'].tolist()
            chr2_gene_order = chr2_common_data['busco_id'].tolist()
        
            
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
                
           
            # Check for any genes with same position but different movement
            for gene_id in chr_common_genes:
                pos1 = chr1_positions[gene_id]
                pos2 = chr2_positions[gene_id]
               
            # Calculate movements for this chromosome pair
            chromosome_movement_sequence = []
            
       
            for gene_id in chr_common_genes:
                pos1 = chr1_positions[gene_id]  # source position
                pos2 = chr2_positions[gene_id]  # target position  
                movement = pos2 - pos1  # target - source (positive = move right)
                # movement = pos1 - pos2  # source - target (positive = move left)
                chromosome_movement_sequence.append((gene_id, pos1, movement, pos2))  # 4-tuple
                
            # Sort by position
            chromosome_movement_sequence.sort(key=lambda x: x[1])
            
            chromosome_pair_name = f"{chr1}_vs_{best_match_chr}"
            chromosome_sequences[chromosome_pair_name] = chromosome_movement_sequence
        
    return chromosome_sequences



def recalculate_movements(current_sequence, target_positions):
    """Recalculate movements based on current vs target positions"""

    updated_sequence = []
    for gene_id, current_rank, old_movement, target_pos in current_sequence:
        target_rank = target_positions[gene_id]  # linearis rank
        new_movement = target_rank - current_rank  # target - current (correct direction)
        updated_sequence.append((gene_id, current_rank, new_movement, target_pos))
        
    return updated_sequence
    

def iterative_detection(movement_sequence, max_iterations=1000):
    
    target_positions = {gene_id: target_pos for gene_id, _, _, target_pos in movement_sequence}
    

    print("========================")
    
    print(f"Starting iterative detection with {len(movement_sequence)} genes")
    
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

        positive_movements = [move for _, _, move, _ in current_sequence if move > 0]
        negative_movements = [move for _, _, move, _ in current_sequence if move < 0]
        sum_positive = sum(positive_movements)
        sum_negative = abs(sum(negative_movements))
        
        print(f"  CONVERGENCE STATUS:")
        print(f"    Total movement: {current_total_movement:.2f}")
        print(f"    Sum of + movements: {sum_positive:.2f} ({len(positive_movements)} genes)")
        print(f"    Sum of - movements: {sum_negative:.2f} ({len(negative_movements)} genes)")
        
        previous_total_movement = current_total_movement
        
        # Check for flip patterns first (more efficient)

        flip_patterns = detect_extended(current_sequence)

            
 
        if flip_patterns:
            
            non_overlapping_flips = find_non_overlapping_flips(flip_patterns)
            
            # Calculate current total movement before applying inversions
            current_total_movement = calculate_total_movement(current_sequence)
            
            # Test each flip and only apply those that reduce total movement
            valid_flips = []
           
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
                else:
                    # Get sample genes from the flip to show what's being rejected
                    sample_genes = [current_sequence[i][0] for i in range(start_idx, min(start_idx + 3, end_idx + 1))]
                    sample_movements = [current_sequence[i][2] for i in range(start_idx, min(start_idx + 3, end_idx + 1))]
             
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
              
        
        # Process adjacencies only if no flips were applied
        if not applied_inversion:
            adjacency_inversions = detect_adjacency_inversions(current_sequence)
            
            if adjacency_inversions:
                
                # Find non-overlapping adjacencies for batch processing
                non_overlapping = find_non_overlapping_adjacencies(adjacency_inversions, current_sequence)
            
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
                    
                    # Check if focus gene is involved in this adjacency
                    gene1_id, _, _, _ = current_sequence[index1]
                    gene2_id, _, _, _ = current_sequence[index2]
                    
                    if movement_change <= flexible_threshold:
                        valid_adjacencies.append((index1, index2))
                    else:
                        gene1_id, _, move1, _ = current_sequence[index1]
                        gene2_id, _, move2, _ = current_sequence[index2]
                
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
        
        # Termination condition
        if not applied_inversion:
            
            #TRANSPOSITION
            # Check for translocation patterns before terminating
            # translocation_events = detect_and_apply_translocations(current_sequence, iteration)
            
            # if translocation_events:
            #     inversion_events.extend(translocation_events)
            #     applied_inversion = True
            # else:
            break
        
        # Progress check: stop if very few large movements remain
        if iteration % 50 == 0:
            large_movements = sum(1 for _, _, move, _ in current_sequence if abs(move) > 2.0)
            if large_movements < 50:  # Most work done
            
                break
        
        # Cycle detection: check if we're stuck in a loop
        current_state_hash = hash(str(current_sequence))
        if current_state_hash in previous_states:
            cycle_start = previous_states.index(current_state_hash)
            cycle_length = len(previous_states) - cycle_start
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
            


    if iteration % 10 == 0:
        remaining_movements = sum(1 for _, _, move, _ in current_sequence if abs(move) > 0.1)
    
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
    
    # Determine convergence type first
    if iteration >= max_iterations:
        converged = False
    elif final_non_zero == 0:
        converged = True
    elif final_large_movements < 10:
        converged = True
    else:
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



