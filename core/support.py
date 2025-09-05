import logging

logger = logging.getLogger(__name__)





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
    
    print(f"  DEBUG: find_non_overlapping_flips - Input patterns: {flip_patterns}")
    
    # Sort by flip size (descending) to prioritize larger flips
    sorted_flips = sorted(flip_patterns, key=lambda x: x[2], reverse=True)
    print(f"  DEBUG: find_non_overlapping_flips - Sorted patterns: {sorted_flips}")
    
    non_overlapping = []
    used_ranges = set()
    
    for start_i, end_j, flip_indicator in sorted_flips:
        print(f"  DEBUG: Processing flip pattern {start_i}-{end_j} (size: {flip_indicator})")
        
        # Check if this flip range overlaps with any already used ranges
        overlap_found = False
        for used_start, used_end in used_ranges:
            # Check if ranges overlap: (start_i <= used_end) and (end_j >= used_start)
            if start_i <= used_end and end_j >= used_start:
                print(f"  DEBUG: Flip {start_i}-{end_j} overlaps with {used_start}-{used_end}")
                overlap_found = True
                break
        
        if not overlap_found:
            non_overlapping.append((start_i, end_j, flip_indicator))
            used_ranges.add((start_i, end_j))
            print(f"  DEBUG: Flip {start_i}-{end_j} ACCEPTED - added to non_overlapping")
            print(f"  DEBUG: Updated used_ranges: {used_ranges}")
        else:
            print(f"  DEBUG: Flip {start_i}-{end_j} REJECTED - overlap detected")
    
    print(f"  DEBUG: find_non_overlapping_flips - Final result: {non_overlapping}")
    return non_overlapping



def apply_batch_with_sequential_fallback(current_sequence, batch_operations, operation_type='flip', target_positions=None):
    """
    Apply batch operations with sequential fallback for interdependent operations.
    
    1. First try concurrent batch (if all independent)
    2. If conflicts found, apply sequentially in optimal order
    """
    from core.reverse import apply_adjacency_inversion, apply_flip_inversion
    
    # Step 1: Check for segment independence
    is_independent = validate_segment_independence(current_sequence, batch_operations, operation_type)
    
    if is_independent:
        # All operations are independent - apply concurrently
        print(f"  🎯 CONCURRENT: Applying {len(batch_operations)} independent {operation_type}s")
        if operation_type == 'flip':
            return apply_concurrent_batch_flips(current_sequence, batch_operations)
        else:
            return apply_concurrent_batch_adjacencies(current_sequence, batch_operations)
    
    else:
        # Operations have conflicts - apply sequentially
        print(f"  🔄 SEQUENTIAL FALLBACK: Applying {len(batch_operations)} interdependent {operation_type}s one by one")
        
        updated_sequence = current_sequence.copy()
        all_records = []
        
        # Apply operations one by one
        for i, operation in enumerate(batch_operations):
            print(f"    Step {i+1}/{len(batch_operations)}: Applying {operation}")
            
            if operation_type == 'flip':
                start_idx, end_idx, flip_size = operation
                updated_sequence, record = apply_flip_inversion(updated_sequence, start_idx, end_idx, flip_size)
            else:
                index1, index2 = operation
                updated_sequence, record = apply_adjacency_inversion(updated_sequence, index1, index2)
            
            record['sequential_step'] = i + 1
            record['total_sequential_steps'] = len(batch_operations)
            all_records.append(record)
            
            # NOTE: Movement recalculation is already handled within apply functions
            # No need for additional recalculation here
            pass
        
        return updated_sequence, all_records


def apply_concurrent_batch_flips(current_sequence, batch_operations):
    """
    Apply multiple flip operations concurrently (when they are independent).
    """
    
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
    print(f"🔍 SEGMENT INDEPENDENCE: Validating {len(batch_operations)} {operation_type} operations")
    
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
    
    # Report results
    if len(overlaps) == 0:
        print(f"  ✅ SEGMENT INDEPENDENCE CONFIRMED")
        print(f"     • {len(segments)} genes with independent movement paths")
        print(f"     • No segment overlaps found")
        
        # Show segments for debugging
        print("     • Movement segments:")
        for seg in segments:
            print(f"       - {seg['gene_id']}: [{seg['segment_start']}-{seg['segment_end']}] ({seg['current']}→{seg['target']})")
            
        return True
    else:
        print(f"  ❌ SEGMENT INDEPENDENCE VIOLATED")
        print(f"     • {len(overlaps)} segment overlaps found:")
        for overlap in overlaps:
            print(f"       - {overlap['gene1']} vs {overlap['gene2']}: overlap at [{overlap['overlap_start']}-{overlap['overlap_end']}]")
            print(f"         ({overlap['operation1']} conflicts with {overlap['operation2']})")
        
        return False

