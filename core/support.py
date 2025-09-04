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
    

    # Sort by first index to process left-to-right
    sorted_adjacencies = sorted(adjacency_inversions, key=lambda x: x[0])
    
    for index1, index2 in sorted_adjacencies:
       
        # CRITICAL: Verify they're still adjacent (consecutive indices)
        if abs(index2 - index1) != 1:
          
            continue  # Skip non-adjacent pairs
            
        # Check if focus gene is involved in this adjacency
        focus_gene_involved = False
        try:
            gene1_id = current_sequence[index1][0]
            gene2_id = current_sequence[index2][0]
     
        except IndexError:
          
            continue
        
        if index1 not in used_indices and index2 not in used_indices:  
            
            # CRITICAL: Double-check the genes are consecutive in current sequence
            gene1_id = current_sequence[index1][0]
            gene2_id = current_sequence[index2][0]
            
            # Verify this is still a valid adjacency pattern
            gene1_movement = current_sequence[index1][2]
            gene2_movement = current_sequence[index2][2]
            
            # Check if movements have opposite signs (required for adjacency inversion)
            has_opposite_signs = ((gene1_movement > 0 and gene2_movement < 0) or 
                                 (gene1_movement < 0 and gene2_movement > 0))
            
            if not has_opposite_signs:
                  continue
            
            non_overlapping.append((index1, index2))
            used_indices.add(index1)
            used_indices.add(index2)
            
        else:
            pass  # Index already used
    
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
        
        else:
            pass  # Overlap detected
    
    return non_overlapping


def apply_batch_with_sequential_fallback(current_sequence, batch_operations, operation_type='flip', target_positions=None):
    """
    Apply batch operations with sequential fallback for interdependent operations.
    
    1. First try concurrent batch (if all independent)
    2. If conflicts found, apply sequentially in optimal order
    """
    # Lazy import to avoid circular dependency
    from core.reverse import apply_flip_inversion, apply_adjacency_inversion
    
    # Step 1: Check for segment independence
    is_independent = validate_segment_independence(current_sequence, batch_operations, operation_type)
    
    if is_independent:
      
        if operation_type == 'flip':
            return apply_concurrent_batch_flips(current_sequence, batch_operations)
        else:
            return apply_concurrent_batch_adjacencies(current_sequence, batch_operations)
    
    else:
      
        updated_sequence = current_sequence.copy()
        all_records = []
        
        # Apply operations one by one
        for i, operation in enumerate(batch_operations):
            
            if operation_type == 'flip':
                start_idx, end_idx, flip_size = operation
                updated_sequence, record = apply_flip_inversion(updated_sequence, start_idx, end_idx, flip_size)
            else:
                index1, index2 = operation
                updated_sequence, record = apply_adjacency_inversion(updated_sequence, index1, index2)
            
            record['sequential_step'] = i + 1
            record['total_sequential_steps'] = len(batch_operations)
            all_records.append(record)
            
            pass
        
        return updated_sequence, all_records


def apply_concurrent_batch_flips(current_sequence, batch_operations):
    """
    Apply multiple flip operations concurrently (when they are independent).
    """
    # Lazy import to avoid circular dependency
    from core.reverse import apply_flip_inversion
    
    updated_sequence = current_sequence.copy()
    all_records = []
    
    for i, (start_idx, end_idx, flip_size) in enumerate(batch_operations):
        updated_sequence, record = apply_flip_inversion(updated_sequence, start_idx, end_idx, flip_size)
        record['batch_position'] = i
        record['total_in_batch'] = len(batch_operations)
        all_records.append(record)
    
    return updated_sequence, all_records


def apply_concurrent_batch_adjacencies(current_sequence, batch_operations):
    """
    Apply multiple adjacency operations concurrently (when they are independent).
    """
    # Lazy import to avoid circular dependency
    from core.reverse import apply_adjacency_inversion
    
    updated_sequence = current_sequence.copy()
    all_records = []
    
    for i, (index1, index2) in enumerate(batch_operations):
        updated_sequence, record = apply_adjacency_inversion(updated_sequence, index1, index2)
        record['batch_position'] = i
        record['total_in_batch'] = len(batch_operations)
        all_records.append(record)
    
    return updated_sequence, all_records


def validate_segment_independence(current_sequence, batch_operations, operation_type='flip'):
    """
    Validate that batch operations have independent movement segments.
    Each gene's current→target path must not overlap with any other gene's path.
    """
    
    # Extract all genes affected by batch operations
    affected_genes = []
    
    for operation in batch_operations:
        if operation_type == 'flip':
            start_idx, end_idx, flip_size = operation
            for i in range(start_idx, end_idx + 1):
                gene_id, current_pos, movement, target_pos = current_sequence[i]
                affected_genes.append({
                    'gene_id': gene_id,
                    'current': current_pos,
                    'target': target_pos,
                    'operation': f'flip[{start_idx}:{end_idx}]'
                })
        elif operation_type == 'adjacency':
            index1, index2 = operation
            for idx in [index1, index2]:
                gene_id, current_pos, movement, target_pos = current_sequence[idx]
                affected_genes.append({
                    'gene_id': gene_id,
                    'current': current_pos,
                    'target': target_pos,
                    'operation': f'adjacency[{index1},{index2}]'
                })
    
    # Calculate movement segments for each gene
    segments = []
    for gene in affected_genes:
        segment_start = min(gene['current'], gene['target'])
        segment_end = max(gene['current'], gene['target'])
        
        segments.append({
            'gene_id': gene['gene_id'],
            'operation': gene['operation'],
            'current': gene['current'],
            'target': gene['target'],
            'segment_start': segment_start,
            'segment_end': segment_end,
            'segment_span': segment_end - segment_start
        })
    
    # Check for segment overlaps
    overlaps = []
    for i in range(len(segments)):
        for j in range(i + 1, len(segments)):
            seg1, seg2 = segments[i], segments[j]
            
            # Check if segments overlap
            has_overlap = (seg1['segment_start'] <= seg2['segment_end']) and \
                         (seg2['segment_start'] <= seg1['segment_end'])
            
            if has_overlap:
                overlap_start = max(seg1['segment_start'], seg2['segment_start'])
                overlap_end = min(seg1['segment_end'], seg2['segment_end'])
                
                overlaps.append({
                    'gene1': seg1['gene_id'],
                    'gene2': seg2['gene_id'],
                    'operation1': seg1['operation'],
                    'operation2': seg2['operation'],
                    'overlap_start': overlap_start,
                    'overlap_end': overlap_end,
                    'overlap_size': overlap_end - overlap_start
                })
    
