def calculate_total_movement(sequence):
    """Calculate total absolute movement for a sequence."""
    return sum(abs(move) for _, _, move, _ in sequence)



def detect_translocation_patterns(sequence, min_movement_threshold=10):
    """
    Detect genes that need large movements that can't be resolved by inversions.
    
    Args:
        sequence: [(gene_id, position, movement)]
        min_movement_threshold: Minimum movement distance to consider for translocation
        
    Returns:
        list: [(gene_id, position, movement)] of genes needing translocation
    """
    translocation_candidates = []
    
    for gene_id, position, movement in sequence:
        if abs(movement) >= min_movement_threshold:
            translocation_candidates.append((gene_id, position, movement))
    
    # Sort by movement magnitude (descending)
    translocation_candidates.sort(key=lambda x: abs(x[2]), reverse=True)
    
    return translocation_candidates


def group_into_translocation_blocks(candidates, sequence, max_gap=3):
    """
    Group genes into contiguous blocks for translocation.
    
    Args:
        candidates: [(gene_id, position, movement)] genes needing translocation
        sequence: Full sequence for position lookup
        max_gap: Maximum gap between genes to include in same block
        
    Returns:
        list: [{'genes': [...], 'start_pos': int, 'end_pos': int, 'target_direction': int}]
    """
    if not candidates:
        return []
    
    # Create position mapping
    pos_to_gene = {pos: (gene_id, movement) for gene_id, pos, movement, _ in sequence}
    
    blocks = []
    current_block = []
    
    for gene_id, position, movement in candidates:
        if not current_block:
            current_block = [{'gene_id': gene_id, 'position': position, 'movement': movement}]
        else:
            last_pos = current_block[-1]['position']
            
            # Check if this gene is contiguous with the current block
            if abs(position - last_pos) <= max_gap + 1:
                # Check if movement direction is similar
                last_movement = current_block[-1]['movement']
                if (movement > 0 and last_movement > 0) or (movement < 0 and last_movement < 0):
                    current_block.append({'gene_id': gene_id, 'position': position, 'movement': movement})
                    continue
            
            # Finalize current block and start new one
            if len(current_block) >= 1:  # Minimum block size
                avg_movement = sum(g['movement'] for g in current_block) / len(current_block)
                blocks.append({
                    'genes': current_block,
                    'start_pos': min(g['position'] for g in current_block),
                    'end_pos': max(g['position'] for g in current_block),
                    'target_direction': 1 if avg_movement > 0 else -1,
                    'avg_movement': avg_movement
                })
            
            current_block = [{'gene_id': gene_id, 'position': position, 'movement': movement}]
    
    # Don't forget the last block
    if len(current_block) >= 1:
        avg_movement = sum(g['movement'] for g in current_block) / len(current_block)
        blocks.append({
            'genes': current_block,
            'start_pos': min(g['position'] for g in current_block),
            'end_pos': max(g['position'] for g in current_block),
            'target_direction': 1 if avg_movement > 0 else -1,
            'avg_movement': avg_movement
        })
    
    # Sort blocks by movement magnitude (largest first)
    blocks.sort(key=lambda x: abs(x['avg_movement']), reverse=True)
    
    return blocks


def find_optimal_insertion_point(block, sequence):
    """
    Find the optimal insertion point for a translocation block.
    
    Args:
        block: Block dictionary with genes and movement info
        sequence: Current sequence
        
    Returns:
        tuple: (insertion_position, orientation_change)
    """
    avg_movement = block['avg_movement']
    block_size = len(block['genes'])
    current_start = block['start_pos']
    
    # Calculate target position based on average movement
    target_center = current_start + avg_movement
    
    # Find best insertion point near target
    best_position = max(0, min(len(sequence) - block_size, int(target_center)))
    
    # Determine if orientation should be changed
    # For simplicity, keep original orientation unless movement is very large
    orientation_change = abs(avg_movement) > 50  # Invert if very large movement
    
    return best_position, orientation_change


def apply_translocation(sequence, block, insertion_pos, invert_orientation=False):
    """
    Apply translocation by moving a block of genes to a new position.
    
    Args:
        sequence: [(gene_id, position, movement)]
        block: Block dictionary
        insertion_pos: Where to insert the block
        invert_orientation: Whether to invert the block orientation
        
    Returns:
        tuple: (updated_sequence, translocation_record)
    """
    updated_sequence = sequence.copy()
    
    # Extract genes in the block
    block_genes = [(g['gene_id'], g['position'], g['movement']) for g in block['genes']]
    
    # Remove block genes from their current positions
    remaining_genes = []
    for gene_id, pos, movement, _ in updated_sequence:
        if gene_id not in [g[0] for g in block_genes]:
            remaining_genes.append((gene_id, pos, movement))
    
    # Invert block order if requested
    if invert_orientation:
        block_genes = list(reversed(block_genes))
    
    # Insert block at new position
    final_sequence = []
    block_inserted = False
    
    for i, (gene_id, pos, movement) in enumerate(remaining_genes):
        # Insert block when we reach the target position
        if i == insertion_pos and not block_inserted:
            for block_gene_id, _, block_movement in block_genes:
                # Calculate new movement after translocation
                new_position = len(final_sequence)
                # Movement should be reduced since we moved toward target
                new_movement = block_movement * 0.3  # Significant reduction
                final_sequence.append((block_gene_id, new_position, new_movement))
            block_inserted = True
        
        # Add the regular gene with updated position
        new_position = len(final_sequence)
        final_sequence.append((gene_id, new_position, movement))
    
    # If block wasn't inserted yet, add it at the end
    if not block_inserted:
        for block_gene_id, _, block_movement in block_genes:
            new_position = len(final_sequence)
            new_movement = block_movement * 0.3
            final_sequence.append((block_gene_id, new_position, new_movement))
    
    # Create translocation record
    translocation_record = {
        'type': 'translocation',
        'genes': [g[0] for g in block_genes],
        'original_positions': [g[1] for g in block_genes],
        'insertion_position': insertion_pos,
        'orientation_inverted': invert_orientation,
        'gene_inversions': len(block_genes),
        'block_size': len(block_genes)
    }
    
    return final_sequence, translocation_record


def detect_and_apply_translocations(sequence, iteration):
    """
    Main function to detect and apply beneficial translocations.
    
    Args:
        sequence: [(gene_id, position, movement)]
        iteration: Current iteration number
        
    Returns:
        list: List of translocation events applied
    """
    current_total_movement = calculate_total_movement(sequence)
    
    # Detect translocation candidates
    candidates = detect_translocation_patterns(sequence)
    
    if not candidates:
        return []
    
    print(f"    Found {len(candidates)} translocation candidates")
    
    # Group into blocks
    blocks = group_into_translocation_blocks(candidates, sequence)
    
    if not blocks:
        return []
    
    print(f"    Grouped into {len(blocks)} potential translocation blocks")
    
    # Test each block for beneficial translocation
    applied_translocations = []
    
    for block in blocks[:3]:  # Limit to top 3 blocks per iteration
        # Find optimal insertion point
        insertion_pos, invert_orientation = find_optimal_insertion_point(block, sequence)
        
        # Test translocation
        test_sequence, translocation_record = apply_translocation(
            sequence, block, insertion_pos, invert_orientation
        )
        
        new_total_movement = calculate_total_movement(test_sequence)
        movement_reduction = current_total_movement - new_total_movement
        
        if movement_reduction > 0.1:  # Only apply if significant improvement
            print(f"    Translocation accepted: {current_total_movement:.1f} -> {new_total_movement:.1f} "
                  f"(reduction: {movement_reduction:.1f})")
            print(f"      Block: {len(block['genes'])} genes, avg_movement: {block['avg_movement']:.1f}")
            
            # Apply this translocation
            sequence = test_sequence
            current_total_movement = new_total_movement
            
            # Update record with iteration info
            translocation_record['iteration'] = iteration
            translocation_record['movement_reduction'] = movement_reduction
            
            applied_translocations.append(translocation_record)
        else:
            print(f"    Translocation rejected: would change movement {current_total_movement:.1f} -> {new_total_movement:.1f}")
    
    return applied_translocations






def analyze_transposition_patterns(non_converged_df):
    """
    Analyze patterns in non-converged genes to identify potential transposition candidates.
    
    Args:
        non_converged_df: DataFrame of non-converged genes
        
    Returns:
        dict: Analysis results with transposition patterns
    """
    patterns = {
        'chromosome_transpositions': {},
        'position_clusters': {},
        'gene_groups': []
    }
    
    # Analyze chromosome differences
    chr_diff = non_converged_df[non_converged_df['convergence_status'] == 'CHROMOSOME_DIFFERENT']
    if len(chr_diff) > 0:
        chr_mapping = chr_diff.groupby(['linearis_chromosome', 'converged_rufipes_chromosome']).size()
        patterns['chromosome_transpositions'] = chr_mapping.to_dict()
    
    # Analyze position differences
    pos_diff = non_converged_df[non_converged_df['convergence_status'] == 'POSITION_DIFFERENT']
    if len(pos_diff) > 0:
        # Group by chromosome and analyze position differences
        for chr_name in pos_diff['linearis_chromosome'].unique():
            chr_data = pos_diff[pos_diff['linearis_chromosome'] == chr_name]
            patterns['position_clusters'][chr_name] = {
                'total_genes': len(chr_data),
                'avg_position_diff': chr_data['position_difference'].mean(),
                'max_position_diff': chr_data['position_difference'].max(),
                'min_position_diff': chr_data['position_difference'].min()
            }
    
    return patterns
