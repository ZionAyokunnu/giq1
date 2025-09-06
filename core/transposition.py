"""
Transposition Detection and Application Module
Implements sophisticated transposition pattern detection with contiguity resolution
"""

import logging
from typing import List, Dict, Tuple, Optional
from pathlib import Path

logger = logging.getLogger(__name__)

def detect_transposition_patterns(movement_sequence):
    """
    Detect transposition patterns following detect_extended() pattern.
    
    Requirements:
    1. Symmetrical magnitude blocks: +6+6-2-2-2-2-2-2 (12|12)
    2. Opposite signs: positive block and negative block
    3. Contiguous same values within each block
    4. Perfect: already contiguous, Almost-perfect: needs contiguity resolution
        
    Returns:
        list: [TranspositionPattern] sorted by magnitude (largest first)
    """
    transposition_patterns = []
    
    print(f"    DEBUG: Scanning for transposition patterns in sequence of length {len(movement_sequence)}")
    
    # Scan all possible segments (following detect_extended pattern)
    for start_i in range(len(movement_sequence)):
        for end_j in range(start_i + 1, len(movement_sequence)):  # Minimum length 2
            
            # Extract movement pattern for this segment
            segment_movements = []
            segment_genes = []
            for k in range(start_i, end_j + 1):
                gene_id, _, movement, _ = movement_sequence[k]
                segment_movements.append(movement)
                segment_genes.append(gene_id)
            
            # DEBUG: Show all segments being tested
            if len(segment_movements) >= 3:  # Only show segments of length 3+
                print(f"    DEBUG: Testing segment [{start_i}-{end_j}]: {segment_genes} with movements {segment_movements}")
            
            # Test if this segment is a valid transposition pattern
            pattern_result = detect_transposition_in_segment(
                segment_movements, start_i, end_j, segment_genes
            )
            
            if pattern_result:
                pattern_result['start_idx'] = start_i
                pattern_result['end_idx'] = end_j
                pattern_result['genes'] = segment_genes
                transposition_patterns.append(pattern_result)
                print(f"    DEBUG: Found transposition pattern [{start_i}-{end_j}]: {segment_genes} with movements {segment_movements}")
                print(f"    DEBUG:   Magnitude: {pattern_result['total_magnitude']}, Perfect: {pattern_result['is_perfect']}, Almost-perfect: {pattern_result['is_almost_perfect']}")
    
    # Sort by total magnitude (largest first)
    transposition_patterns.sort(key=lambda x: x['total_magnitude'], reverse=True)
    
    print(f"    DEBUG: Found {len(transposition_patterns)} total transposition patterns")
    for i, pattern in enumerate(transposition_patterns):
        print(f"    DEBUG: Pattern {i+1}: [{pattern['start_idx']}-{pattern['end_idx']}] magnitude={pattern['total_magnitude']}, perfect={pattern['is_perfect']}")
    
    # Filter out overlapping patterns - keep largest non-overlapping ones
    filtered_patterns = []
    used_positions = set()

    for pattern in transposition_patterns:
        pattern_positions = set(range(pattern['start_idx'], pattern['end_idx'] + 1))
        
        # Check if this pattern overlaps with any already selected pattern
        if not pattern_positions.intersection(used_positions):
            filtered_patterns.append(pattern)
            used_positions.update(pattern_positions)
            print(f"    DEBUG: Selected pattern [{pattern['start_idx']}-{pattern['end_idx']}] magnitude={pattern['total_magnitude']}")
        else:
            print(f"    DEBUG: Skipped overlapping pattern [{pattern['start_idx']}-{pattern['end_idx']}] magnitude={pattern['total_magnitude']}")

    print(f"    DEBUG: Filtered from {len(transposition_patterns)} to {len(filtered_patterns)} non-overlapping patterns")

    # Replace the original list with filtered results
    transposition_patterns = filtered_patterns
    
    return transposition_patterns

def detect_transposition_in_segment(movements, start_idx, end_idx, genes):
    """
    Test if a movement segment is a valid transposition pattern.
    
    Pattern: +6+6-2-2-2-2-2-2 where left magnitude = right magnitude
    
    Returns:
        dict: Pattern info if valid, None if invalid
    """
    if len(movements) < 2:
        return None
    
    # Split into positive and negative sections
    positive_section = []
    negative_section = []
    positive_indices = []
    negative_indices = []
    
    for i, move in enumerate(movements):
        if move > 0:
            positive_section.append(move)
            positive_indices.append(i)
        elif move < 0:
            negative_section.append(move)
            negative_indices.append(i)
        # Skip zeros
    
    if not positive_section or not negative_section:
        return None  # Must have both positive and negative
    
    # Check segment count symmetry for translocation vs inversion
    positive_count = len(positive_section)
    negative_count = len(negative_section)
    segment_count_symmetrical = (positive_count == negative_count)
    
    # For translocation/transposition: segments must NOT be symmetrical
    if segment_count_symmetrical:
        print(f"    DEBUG: Rejecting pattern - symmetrical segment count ({positive_count}+{negative_count}), likely inversion")
        return None  # This is likely an inversion, not a transposition
    
    print(f"    DEBUG: ✅ Non-symmetrical segment count ({positive_count}+{negative_count}) - checking further...")
    
    # Check magnitude symmetry
    positive_magnitude = sum(positive_section)
    negative_magnitude = abs(sum(negative_section))
    
    # If magnitude is not symmetrical, it's almost-perfect (needs contiguity resolution)
    magnitude_symmetrical = (positive_magnitude == negative_magnitude)
    
    # Check contiguity within each section
    positive_contiguous = check_contiguity(positive_indices)
    negative_contiguous = check_contiguity(negative_indices)
    
    # Check for same values within each block
    positive_same_values = len(set(positive_section)) == 1  # All same value
    negative_same_values = len(set([abs(x) for x in negative_section])) == 1  # All same magnitude
    
    if not (positive_same_values and negative_same_values):
        print(f"    DEBUG: Rejecting pattern - not same values: positive_same={positive_same_values}, negative_same={negative_same_values}")
        print(f"    DEBUG:   Positive section: {positive_section}, Negative section: {negative_section}")
        return None  # Must have same values in each block
    
    # Determine if perfect or almost-perfect
    # Perfect: magnitude symmetrical AND both blocks contiguous
    is_perfect = magnitude_symmetrical and positive_contiguous and negative_contiguous
    # Almost-perfect: either magnitude not symmetrical OR only one block contiguous
    is_almost_perfect = not magnitude_symmetrical or (not is_perfect and (positive_contiguous or negative_contiguous))
    
    if not (is_perfect or is_almost_perfect):
        return None  # Must be at least almost-perfect
    
    return {
        'total_magnitude': positive_magnitude,
        'positive_block': positive_section,
        'negative_block': negative_section,
        'positive_indices': positive_indices,
        'negative_indices': negative_indices,
        'is_perfect': is_perfect,
        'is_almost_perfect': is_almost_perfect,
        'positive_contiguous': positive_contiguous,
        'negative_contiguous': negative_contiguous,
        'magnitude_symmetrical': magnitude_symmetrical
    }

def check_contiguity(indices):
    """Check if indices are contiguous (consecutive)"""
    if len(indices) <= 1:
        return True
    
    for i in range(len(indices) - 1):
        if indices[i + 1] - indices[i] != 1:
            return False
    return True

# Your transposition detection and application functions will go here
def detect_and_apply_translocations(movement_sequence):
    """
    Main function to detect and apply translocations in a single call.
    This is a wrapper that combines detection and application.
    
    Args:
        movement_sequence: [(gene_id, position, movement, target_position)]
        
    Returns:
        list: List of translocation events applied (empty if none)
    """
    print(f"    DEBUG: detect_and_apply_translocations called with {len(movement_sequence)} genes")
    
    # Use the existing translocation detection logic
    from .transposition import detect_translocation_patterns, group_into_translocation_blocks, apply_translocation
    
    current_total_movement = sum(abs(move) for _, _, move, _ in movement_sequence)
    
    # Detect translocation candidates using existing logic
    candidates = detect_translocation_patterns(movement_sequence)
    
    if not candidates:
        print(f"    DEBUG: No translocation candidates found")
        return []
    
    print(f"    DEBUG: Found {len(candidates)} translocation candidates")
    
    # Group into blocks using existing logic
    blocks = group_into_translocation_blocks(candidates, movement_sequence)
    
    if not blocks:
        print(f"    DEBUG: No translocation blocks formed")
        return []
    
    print(f"    DEBUG: Grouped into {len(blocks)} potential translocation blocks")
    
    # Test and apply beneficial translocations
    applied_translocations = []
    
    for block in blocks[:3]:  # Limit to top 3 blocks per call
        # Find optimal insertion point using existing logic
        insertion_pos, invert_orientation = find_optimal_insertion_point(block, movement_sequence)
        
        # Test translocation
        test_sequence, translocation_record = apply_translocation(
            movement_sequence, block, insertion_pos, invert_orientation
        )
        
        new_total_movement = sum(abs(move) for _, _, move, _ in test_sequence)
        movement_reduction = current_total_movement - new_total_movement
        
        if movement_reduction > 0.1:  # Only apply if significant improvement
            print(f"    DEBUG: Translocation accepted: {current_total_movement:.1f} -> {new_total_movement:.1f}")
            print(f"    DEBUG: Block: {len(block['genes'])} genes, avg_movement: {block['avg_movement']:.1f}")
            
            # Update record
            translocation_record['movement_reduction'] = movement_reduction
            applied_translocations.append(translocation_record)
            
            # Update sequence for next potential translocation
            movement_sequence = test_sequence
            current_total_movement = new_total_movement
        else:
            print(f"    DEBUG: Translocation rejected: would change movement {current_total_movement:.1f} -> {new_total_movement:.1f}")
    
    print(f"    DEBUG: Applied {len(applied_translocations)} translocations")
    return applied_translocations

def find_contiguity_resolution_operation(trans_pattern, current_sequence, max_contiguity_iterations, 
                                       use_whole_sequence=False):
    """
    Find flip/adjacency operation that can resolve contiguity for almost-perfect transposition.
    
    Args:
        trans_pattern: Almost-perfect transposition pattern
        current_sequence: Current movement sequence
        max_contiguity_iterations: Max attempts to find resolution
        use_whole_sequence: If True, search entire sequence (slower but comprehensive)
        
    Returns:
        dict: Operation that resolves contiguity, or None
    """
    print(f"    DEBUG: Finding contiguity resolution for pattern with magnitude {trans_pattern['total_magnitude']}")
    print(f"    DEBUG: Pattern contiguity - positive_contiguous: {trans_pattern.get('positive_contiguous', 'N/A')}, negative_contiguous: {trans_pattern.get('negative_contiguous', 'N/A')}")
    print(f"    DEBUG: Pattern indices - positive: {trans_pattern.get('positive_indices', 'N/A')}, negative: {trans_pattern.get('negative_indices', 'N/A')}")
    
    # Identify which block needs contiguity resolution
    # For almost-perfect patterns, we need to resolve magnitude asymmetry
    if not trans_pattern.get('magnitude_symmetrical', True):
        print(f"    DEBUG: Pattern has magnitude asymmetry - needs contiguity resolution")
        # Identify which block has lower magnitude
        pos_mag = sum(trans_pattern['positive_block'])
        neg_mag = abs(sum(trans_pattern['negative_block']))
        
        if pos_mag < neg_mag:
            target_block_type = 'positive'
            target_indices = trans_pattern['positive_indices']
        else:
            target_block_type = 'negative'  
            target_indices = trans_pattern['negative_indices']
        
        print(f"    DEBUG: {target_block_type} block has lower magnitude, needs expansion")
        # Continue with gap resolution logic to expand the sparse block
    elif not trans_pattern['positive_contiguous']:
        target_indices = trans_pattern['positive_indices']
        target_block_type = 'positive'
    elif not trans_pattern['negative_contiguous']:
        target_indices = trans_pattern['negative_indices']
        target_block_type = 'negative'
    else:
        print(f"    DEBUG: Pattern is already perfect - no contiguity resolution needed")
        return None
    
    print(f"    DEBUG: {target_block_type} block needs contiguity at indices {target_indices}")
    
    # Determine search range based on the target block
    if use_whole_sequence:
        search_start = 0
        search_end = len(current_sequence) - 1
        print(f"    DEBUG: Using WHOLE SEQUENCE search [0-{search_end}] (slower but comprehensive)")
    else:
        # Create search range based on flexibility FROM THE BLOCK
        from config.settings import CONFIG
        tolerance = CONFIG.get('contiguity_gap_tolerance', 10)
        block_start = min(target_indices)
        block_end = max(target_indices)

        search_start = max(0, block_start - tolerance)
        search_end = block_end + tolerance  # Don't bound by sequence length - we want the full flexibility

        print(f"    DEBUG: Block-based search range [{search_start}-{search_end}] from block [{block_start}-{block_end}] with tolerance {tolerance}")
    
    # Find exact expansion positions (previous and next indices)
    expansion_positions = []

    # Get the target movement value that new genes must have
    target_movement_value = None
    if target_block_type == 'positive':
        target_movement_value = trans_pattern['positive_block'][0]
    else:
        target_movement_value = trans_pattern['negative_block'][0]

    print(f"    DEBUG: Target block {target_indices} needs genes with movement value {target_movement_value}")

    # Check exact next position after the block
    max_target_pos = max(target_indices)
    next_pos = max_target_pos + 1
    if next_pos < len(current_sequence):
        expansion_positions.append(next_pos)
        print(f"    DEBUG: Next expansion position: {next_pos}")

    # Check exact previous position before the block  
    min_target_pos = min(target_indices)
    prev_pos = min_target_pos - 1
    if prev_pos >= 0:
        expansion_positions.append(prev_pos)
        print(f"    DEBUG: Previous expansion position: {prev_pos}")

    print(f"    DEBUG: Expansion positions to check: {expansion_positions}")
    
    # For each expansion position, check if flip/adjacency can make it have the target movement value
    for attempt in range(max_contiguity_iterations):
        for expansion_pos in expansion_positions:
            
            # Try flip resolution first
            flip_op = check_expansion_flip_resolution(expansion_pos, target_movement_value, current_sequence, 
                                                    search_start, search_end)
            if flip_op:
                print(f"    DEBUG: Found flip resolution for expansion: {flip_op}")
                return flip_op
            
            # Try adjacency second
            adjacency_op = check_expansion_adjacency_resolution(expansion_pos, target_movement_value, current_sequence, 
                                                              search_start, search_end)
            if adjacency_op:
                print(f"    DEBUG: Found adjacency resolution for expansion: {adjacency_op}")
                return adjacency_op
    
    # NO AUTOMATIC FALLBACK - just return None if nothing found
    print(f"    DEBUG: No contiguity resolution found after {max_contiguity_iterations} attempts")
    return None


def check_adjacency_resolution(gap_start, gap_end, current_sequence, trans_pattern, 
                             search_start, search_end):
    """Check if adjacency can resolve contiguity gap within search range"""
    # For adjacency, gap must be exactly 1 position
    if gap_end - gap_start != 0:
        return None
    
    gap_pos = gap_start
    if gap_pos >= len(current_sequence) or gap_pos < search_start or gap_pos > search_end:
        return None
    
    # Check if gene at gap position belongs to the transposition pattern
    gap_gene_id = current_sequence[gap_pos][0]
    pattern_genes = [current_sequence[i][0] for i in trans_pattern['positive_indices'] + trans_pattern['negative_indices']]
    
    if gap_gene_id not in pattern_genes:
        return None
    
    # Check for adjacent pairs within search range
    for adj_idx in [gap_pos - 1, gap_pos + 1]:
        if search_start <= adj_idx <= search_end and 0 <= adj_idx < len(current_sequence):
            # Check if this forms a valid adjacency pattern
            gene1_id, _, move1, _ = current_sequence[min(gap_pos, adj_idx)]
            gene2_id, _, move2, _ = current_sequence[max(gap_pos, adj_idx)]
            
            # Must have opposite signs for adjacency
            if (move1 > 0 and move2 < 0) or (move1 < 0 and move2 > 0):
                return {
                    'type': 'adjacency',
                    'index1': min(gap_pos, adj_idx),
                    'index2': max(gap_pos, adj_idx),
                    'resolves_gap': (gap_start, gap_end),
                    'search_method': 'whole_sequence' if search_start == 0 and search_end == len(current_sequence)-1 else 'extended_pattern'
                }
    
    return None


def check_expansion_flip_resolution(expansion_pos, target_movement_value, current_sequence, 
                                  search_start, search_end):
    """Check if flip within search range can make gene at expansion_pos have target_movement_value"""
    if expansion_pos >= len(current_sequence):
        return None
    
    print(f"    DEBUG: Checking if flip can make position {expansion_pos} have movement {target_movement_value}")
    print(f"    DEBUG: Search range [{search_start}-{search_end}], expansion position: {expansion_pos}")
    
    # Find ALL flip patterns in the search range
    for flip_start in range(search_start, search_end):
        for flip_end in range(flip_start + 1, search_end + 1):
            
            if flip_end - flip_start < 1:
                continue
            
            # Extract pattern for this potential flip (only if within sequence bounds)
            if flip_start >= len(current_sequence) or flip_end >= len(current_sequence):
                continue
                
            flip_movements = [current_sequence[k][2] for k in range(flip_start, flip_end + 1)]
            
            from core.rules import detect_flip_in_pattern
            flip_indicator = detect_flip_in_pattern(flip_movements)
            
            if flip_indicator > 0:
                # SIMULATE the flip to verify it achieves target value
                test_sequence = current_sequence.copy()
                from core.reverse import apply_flip_inversion
                
                try:
                    simulated_sequence, _ = apply_flip_inversion(test_sequence, flip_start, flip_end, flip_indicator)
                    
                    # Check if expansion_pos now has target_movement_value (only if within bounds)
                    if expansion_pos >= len(simulated_sequence):
                        continue
                        
                    new_movement = simulated_sequence[expansion_pos][2]
                    
                    print(f"    DEBUG: Simulated flip [{flip_start}-{flip_end}]: position {expansion_pos} movement {current_sequence[expansion_pos][2]} → {new_movement}")
                    
                    if new_movement == target_movement_value:
                        print(f"    DEBUG: ✅ Flip simulation successful! Position {expansion_pos} now has target movement {target_movement_value}")
                        return {
                            'type': 'flip',
                            'start_idx': flip_start,
                            'end_idx': flip_end,
                            'flip_size': flip_indicator,
                            'expansion_pos': expansion_pos,
                            'target_value': target_movement_value,
                            'verified': True
                        }
                    else:
                        print(f"    DEBUG: ❌ Flip simulation failed: got {new_movement}, needed {target_movement_value}")
                        
                except Exception as e:
                    print(f"    DEBUG: Flip simulation error: {e}")
                    continue
    
    print(f"    DEBUG: No flip in search range [{search_start}-{search_end}] achieves target movement {target_movement_value} at position {expansion_pos}")
    return None


def check_expansion_adjacency_resolution(expansion_pos, target_movement_value, current_sequence,
                                       search_start, search_end):
    """Check if adjacency within search range can make gene at expansion_pos have target_movement_value"""
    if expansion_pos >= len(current_sequence):
        return None
    
    print(f"    DEBUG: Checking if adjacency can make position {expansion_pos} have movement {target_movement_value}")
    
    # Check adjacent positions within search range
    for adj_pos in [expansion_pos - 1, expansion_pos + 1]:
        if search_start <= adj_pos <= search_end and 0 <= adj_pos < len(current_sequence):
            
            gene1_move = current_sequence[min(expansion_pos, adj_pos)][2]
            gene2_move = current_sequence[max(expansion_pos, adj_pos)][2]
            
            # Must have opposite signs for adjacency
            if (gene1_move > 0 and gene2_move < 0) or (gene1_move < 0 and gene2_move > 0):
                
                # SIMULATE the adjacency to verify it achieves target value
                test_sequence = current_sequence.copy()
                from core.reverse import apply_adjacency_inversion
                
                try:
                    index1 = min(expansion_pos, adj_pos)
                    index2 = max(expansion_pos, adj_pos)
                    
                    simulated_sequence, _ = apply_adjacency_inversion(test_sequence, index1, index2)
                    
                    # Check if expansion_pos now has target_movement_value
                    new_movement = simulated_sequence[expansion_pos][2]
                    
                    print(f"    DEBUG: Simulated adjacency [{index1}-{index2}]: position {expansion_pos} movement {current_sequence[expansion_pos][2]} → {new_movement}")
                    
                    if new_movement == target_movement_value:
                        print(f"    DEBUG: ✅ Adjacency simulation successful! Position {expansion_pos} now has target movement {target_movement_value}")
                        return {
                            'type': 'adjacency',
                            'index1': index1,
                            'index2': index2,
                            'expansion_pos': expansion_pos,
                            'target_value': target_movement_value,
                            'verified': True
                        }
                    else:
                        print(f"    DEBUG: ❌ Adjacency simulation failed: got {new_movement}, needed {target_movement_value}")
                        
                except Exception as e:
                    print(f"    DEBUG: Adjacency simulation error: {e}")
                    continue
    
    return None


def check_flip_resolution(gap_start, gap_end, current_sequence, trans_pattern, 
                        search_start, search_end):
    """Check if flip can resolve contiguity gap within search range"""
    # Try different flip ranges that include the gap, within search bounds
    for flip_start in range(max(search_start, gap_start - 3), min(search_end, gap_start + 2)):
        for flip_end in range(max(search_start, gap_end - 1), min(search_end + 1, gap_end + 4)):
            
            if flip_end - flip_start < 1 or flip_start < 0 or flip_end >= len(current_sequence):
                continue
            
            # Extract pattern for this potential flip
            flip_movements = []
            for k in range(flip_start, flip_end + 1):
                _, _, movement, _ = current_sequence[k]
                flip_movements.append(movement)
            
            # Check if this is a valid flip pattern (reuse existing logic)
            from core.rules import detect_flip_in_pattern
            flip_indicator = detect_flip_in_pattern(flip_movements)
            
            if flip_indicator > 0:
                return {
                    'type': 'flip',
                    'start_idx': flip_start,
                    'end_idx': flip_end,
                    'flip_size': flip_indicator,
                    'resolves_gap': (gap_start, gap_end),
                    'search_method': 'whole_sequence' if search_start == 0 and search_end == len(current_sequence)-1 else 'extended_pattern'
            }
    
    return None

def apply_transposition_inversion(movement_sequence, start_index, end_index, trans_pattern):
    """
    Apply transposition by moving blocks according to head1-tail2 → head2-tail1 rule.
    
    Args:
        movement_sequence: Current sequence
        start_index: Start of transposition segment  
        end_index: End of transposition segment
        trans_pattern: Transposition pattern info
        
    Returns:
        tuple: (updated_sequence, inversion_record)
    """
    print(f"    DEBUG: Applying transposition on segment [{start_index}-{end_index}]")
    print(f"    DEBUG: Pattern magnitude: {trans_pattern['total_magnitude']}")
    
    # Calculate total movement before transposition
    total_movement_before = sum(abs(move) for _, _, move, _ in movement_sequence)
    print(f"    DEBUG: Total movement before transposition: {total_movement_before}")
    
    updated_sequence = movement_sequence.copy()
    
    # Extract the segment
    segment = updated_sequence[start_index:end_index + 1]
    
    print(f"    DEBUG: Original segment [{start_index}-{end_index}]:")
    for i, (gene_id, pos, move, target) in enumerate(segment):
        print(f"      [{start_index + i}] {gene_id}: pos={pos}, move={move}, target={target}")
    
    # Identify positive and negative blocks within segment
    positive_block_genes = []
    negative_block_genes = []
    
    for i, (gene_id, pos, move, target) in enumerate(segment):
        if move > 0:
            positive_block_genes.append((gene_id, pos, move, target, start_index + i))
        elif move < 0:
            negative_block_genes.append((gene_id, pos, move, target, start_index + i))
    
    print(f"    DEBUG: Positive block: {len(positive_block_genes)} genes")
    for gene_id, pos, move, target, orig_idx in positive_block_genes:
        print(f"      {gene_id}: pos={pos}, move={move}, target={target}")
    
    print(f"    DEBUG: Negative block: {len(negative_block_genes)} genes")
    for gene_id, pos, move, target, orig_idx in negative_block_genes:
        print(f"      {gene_id}: pos={pos}, move={move}, target={target}")
    
    # Apply head1-tail2 → head2-tail1 transformation
    # Positive block (head1) moves to where negative block (tail2) target positions are
    # Negative block (head2) moves to where positive block (tail1) target positions are
    
    new_segment = [None] * len(segment)
    
    # Place positive block at negative block's target positions
    negative_targets = sorted([target for _, _, _, target, _ in negative_block_genes])
    print(f"    DEBUG: Negative targets (sorted): {negative_targets}")
    
    for i, (gene_id, old_pos, old_move, target, orig_idx) in enumerate(positive_block_genes):
        if i < len(negative_targets):
            new_pos = start_index + (orig_idx - start_index)  # Keep relative position for now
            new_movement = negative_targets[i] - new_pos  # Recalculate movement
            new_segment[orig_idx - start_index] = (gene_id, new_pos, new_movement, negative_targets[i])
            print(f"    DEBUG: {gene_id} positive → negative target: pos={new_pos}, move={new_movement}, target={negative_targets[i]}")
    
    # Place negative block at positive block's target positions  
    positive_targets = sorted([target for _, _, _, target, _ in positive_block_genes])
    print(f"    DEBUG: Positive targets (sorted): {positive_targets}")
    
    for i, (gene_id, old_pos, old_move, target, orig_idx) in enumerate(negative_block_genes):
        if i < len(positive_targets):
            new_pos = start_index + (orig_idx - start_index)  # Keep relative position for now
            new_movement = positive_targets[i] - new_pos  # Recalculate movement
            new_segment[orig_idx - start_index] = (gene_id, new_pos, new_movement, positive_targets[i])
            print(f"    DEBUG: {gene_id} negative → positive target: pos={new_pos}, move={new_movement}, target={positive_targets[i]}")
    
    # Handle any remaining genes (zeros, etc.)
    for i, item in enumerate(segment):
        if new_segment[i] is None:
            gene_id, pos, move, target = item
            new_pos = start_index + i
            new_movement = target - new_pos
            new_segment[i] = (gene_id, new_pos, new_movement, target)
            print(f"    DEBUG: {gene_id} unchanged: pos={new_pos}, move={new_movement}, target={target}")
    
    # Put updated segment back into sequence
    updated_sequence[start_index:end_index + 1] = new_segment
    
    # Calculate total movement after transposition
    total_movement_after = sum(abs(move) for _, _, move, _ in updated_sequence)
    movement_change = total_movement_after - total_movement_before
    
    print(f"    DEBUG: Total movement after transposition: {total_movement_after}")
    print(f"    DEBUG: Movement change: {movement_change:+.1f}")
    
    print(f"    DEBUG: Updated segment [{start_index}-{end_index}]:")
    for i, (gene_id, pos, move, target) in enumerate(new_segment):
        print(f"      [{start_index + i}] {gene_id}: pos={pos}, move={move}, target={target}")
    
    # Create record
    genes_involved = [gene_id for gene_id, _, _, _, _ in positive_block_genes + negative_block_genes]
    positions_involved = list(range(start_index, end_index + 1))
    
    transposition_record = {
        'type': 'transposition',
        'start_idx': start_index,
        'end_idx': end_index,
        'total_magnitude': trans_pattern['total_magnitude'],
        'positive_block_size': len(positive_block_genes),
        'negative_block_size': len(negative_block_genes),
        'positions': positions_involved,
        'genes': genes_involved,
        'gene_inversions': len(genes_involved),
        'movement_before': total_movement_before,
        'movement_after': total_movement_after,
        'movement_change': movement_change
    }
    
    print(f"    DEBUG: Transposition applied - {len(genes_involved)} genes rearranged")
    print(f"    DEBUG: Movement change: {total_movement_before:.1f} → {total_movement_after:.1f} ({movement_change:+.1f})")
    
    return updated_sequence, transposition_record
