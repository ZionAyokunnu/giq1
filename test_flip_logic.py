def apply_flip_inversion_ACTUAL(movement_sequence, start_index, end_index, flip_indicator):
    """
    ACTUAL FLIPPING: Physically reverse gene order by correct index reassignment.
    This is the proper implementation using your logic.
    """
    updated_sequence = movement_sequence.copy()
    segment = updated_sequence[start_index:end_index + 1]
    
    print(f"🔄 ACTUAL FLIP: segment[{start_index}:{end_index+1}] ({len(segment)} genes)")
    print(f"   BEFORE: {[(gene_id, pos, move) for gene_id, pos, move, _ in segment]}")
    
    # Step 1: Extract original positions (these stay the same - just reassigned)
    original_positions = [pos for _, pos, _, _ in segment]
    print(f"   Positions to reassign: {original_positions}")
    
    # Step 2: ACTUAL FLIPPING using your correct logic
    final_segment = []
    
    for i in range(len(segment)):
        # YOUR LOGIC: Get gene from opposite end (actual physical reversal)
        reversed_gene_index = len(segment) - 1 - i
        gene_id, original_pos, old_move, target = segment[reversed_gene_index]
        
        # YOUR LOGIC: Assign this gene to position at index i
        new_position = original_positions[i]
        
        # Calculate movement update based on actual distance moved
        actual_distance_moved = new_position - original_pos
        new_move = old_move - actual_distance_moved
        
        print(f"   Gene {gene_id}: index {reversed_gene_index}→{i}, pos {original_pos}→{new_position}, move {old_move}→{new_move}")
        
        # Verify: new_move should equal target - new_position
        expected = target - new_position
        assert abs(new_move - expected) < 0.001, f"Math error: {gene_id} calculated {new_move} but expected {expected}"
        
        final_segment.append((gene_id, new_position, new_move, target))
    
    # Update the main sequence
    updated_sequence[start_index:end_index + 1] = final_segment
    
    print(f"   AFTER:  {[(gene_id, pos, move) for gene_id, pos, move, _ in final_segment]}")
    
    # Verify gene order is actually physically reversed
    original_gene_order = [gene_id for gene_id, _, _, _ in segment]
    final_gene_order = [gene_id for gene_id, _, _, _ in final_segment]
    expected_reversed = list(reversed(original_gene_order))
    
    print(f"   Gene order: {original_gene_order} → {final_gene_order}")
    print(f"   Correctly reversed: {final_gene_order == expected_reversed} ✓")
    
    return updated_sequence, {
        'type': 'flip',
        'flip_indicator': flip_indicator,
        'genes': original_gene_order,
        'positions': original_positions,
        'gene_inversions': len(segment)
    }


def test_all_pattern_types_with_actual_flipping():
    """Test all pattern types from codebase with actual physical flipping"""
    
    print("🧪 TESTING ALL PATTERN TYPES WITH ACTUAL FLIPPING")
    print("=" * 80)
    
    # Pattern Type 1: Simple Adjacent (+/-)
    print("\n📍 TYPE 1: Simple Adjacent Pattern [+5, -3]")
    sequence1 = [('A', 10, 5, 15), ('B', 11, -3, 8)]
    result1, record1 = apply_flip_inversion_ACTUAL(sequence1, 0, 1, 1)
    
    # Pattern Type 2A: Odd Length Incremental [4, 2, 0, -1, -3]
    print("\n📍 TYPE 2A: Odd Length Incremental [4, 2, 0, -1, -3]")
    sequence2 = [('C', 20, 4, 24), ('D', 21, 2, 23), ('E', 22, 0, 22), 
                ('F', 23, -1, 22), ('G', 24, -3, 21)]
    result2, record2 = apply_flip_inversion_ACTUAL(sequence2, 0, 4, 2)
    
    # Pattern Type 2B: Even Length Incremental [6, 3, -2, -5]
    print("\n📍 TYPE 2B: Even Length Incremental [6, 3, -2, -5]")
    sequence3 = [('H', 30, 6, 36), ('I', 31, 3, 34), ('J', 32, -2, 30), ('K', 33, -5, 28)]
    result3, record3 = apply_flip_inversion_ACTUAL(sequence3, 0, 3, 2)
    
    # Complex Pattern: Large Odd Length [8, 5, 2, 0, -1, -4, -7]
    print("\n📍 COMPLEX: Large Odd Length [8, 5, 2, 0, -1, -4, -7]")
    sequence4 = [('L', 40, 8, 48), ('M', 41, 5, 46), ('N', 42, 2, 44), ('O', 43, 0, 43),
                ('P', 44, -1, 43), ('Q', 45, -4, 41), ('R', 46, -7, 39)]
    result4, record4 = apply_flip_inversion_ACTUAL(sequence4, 0, 6, 3)
    
    print("\n✅ ALL PATTERN TYPES TESTED WITH ACTUAL PHYSICAL FLIPPING")


def verify_your_index_mapping_logic():
    """Verify your index mapping logic is mathematically correct"""
    print("\n🔍 VERIFYING YOUR INDEX MAPPING LOGIC")
    print("=" * 50)
    
    segment_lengths = [2, 3, 4, 5, 6, 7]
    
    for n in segment_lengths:
        print(f"\nSegment length {n}:")
        print("Index mapping (i → reversed_index):")
        for i in range(n):
            reversed_index = n - 1 - i
            print(f"  i={i} → reversed_index={reversed_index}")
        
        # Verify this creates proper reversal
        original_indices = list(range(n))
        mapped_indices = [n - 1 - i for i in range(n)]
        print(f"  Original: {original_indices}")
        print(f"  Mapped:   {mapped_indices}")
        print(f"  Correctly reversed: {mapped_indices == list(reversed(original_indices))}")


def verify_flip_mathematics():
    """
    Test function to verify the flip mathematics with a concrete example.
    """
    print("🧮 VERIFYING FLIP MATHEMATICS")
    print("=" * 50)
    
    # Example: 5 genes in positions [10, 11, 12, 13, 14] with targets [15, 16, 17, 18, 19]
    # This means movements are [5, 5, 5, 5, 5]
    
    test_sequence = [
        ('GeneA', 10, 5, 15),  # at pos 10, needs to go to 15, movement = +5
        ('GeneB', 11, 5, 16),  # at pos 11, needs to go to 16, movement = +5  
        ('GeneC', 12, 5, 17),  # at pos 12, needs to go to 17, movement = +5
        ('GeneD', 13, 5, 18),  # at pos 13, needs to go to 18, movement = +5
        ('GeneE', 14, 5, 19),  # at pos 14, needs to go to 19, movement = +5
    ]
    
    print("BEFORE FLIP:")
    for gene_id, pos, move, target in test_sequence:
        print(f"  {gene_id}: pos={pos}, target={target}, movement={move}")
    
    print(f"\nTotal movement before: {sum(abs(move) for _, _, move, _ in test_sequence)}")
    
    # Apply flip to entire sequence
    flipped_sequence, record = apply_flip_inversion_ACTUAL(test_sequence, 0, 4, 2)
    
    print("\nAFTER FLIP:")
    for gene_id, pos, move, target in flipped_sequence:
        print(f"  {gene_id}: pos={pos}, target={target}, movement={move}")
        # Verify movement is correct
        expected = target - pos
        print(f"    Expected movement: {expected}, Actual: {move}, Match: {abs(expected - move) < 0.001}")
    
    print(f"\nTotal movement after: {sum(abs(move) for _, _, move, _ in flipped_sequence)}")
    
    # Check if positions are correctly reversed
    original_positions = [10, 11, 12, 13, 14]
    final_positions = [pos for _, pos, _, _ in flipped_sequence]
    print(f"\nPosition mapping verification:")
    print(f"  Original: {original_positions}")
    print(f"  Final:    {final_positions}")
    print(f"  Expected: {original_positions}  (positions stay the same)")
    
    # Check if gene order is correctly reversed
    original_genes = [gene for gene, _, _, _ in test_sequence]
    final_genes = [gene for gene, _, _, _ in flipped_sequence]
    expected_genes = list(reversed(original_genes))
    print(f"\nGene order verification:")
    print(f"  Original: {original_genes}")
    print(f"  Final:    {final_genes}")
    print(f"  Expected: {expected_genes}")
    print(f"  Correct:  {final_genes == expected_genes}")


if __name__ == "__main__":
    test_all_pattern_types_with_actual_flipping()
    verify_your_index_mapping_logic()
    verify_flip_mathematics()





"""
Target order: ABCDEFG (positions 0,1,2,3,4,5,6) * Current order: FABCDG

Target order: ABCDEFG (positions 0,1,2,3,4,5,6)
Current order: FABCDGE (positions 0,1,2,3,4,5,6)

Movement pattern: [+5, -1, -1, -1, -1, +1, -2]

F@pos0 → target@pos5 (move +5)
A@pos1 → target@pos0 (move -1)
B@pos2 → target@pos1 (move -1)
C@pos3 → target@pos2 (move -1)
D@pos4 → target@pos3 (move -1)
G@pos5 → target@pos6 (move +1)
E@pos6 → target@pos4 (move -2)

"""

#For each valid adjacent pattern, calculate the movement for the segment

all_segment = [('F', 0, +5, 5), ('A', 1, -1, 0), ('B', 2, -1, 1), 
           ('C', 3, -1, 2), ('D', 4, -1, 3), ('G', 5, +1, 6), 
           ('E', 6, -2, 4)]

#for the first adjacent pattern:
segment = [('F', 0, +5, 5), ('A', 1, -1, 0)]

start_index = 0
end_index = 0
segment_length = 2

for i in range(segment_length//2):
    
    #range(2//2) = range(1) = [0], so i = 0
    left_idx = i
    #left_idx = 0
    
    right_idx = segment_length - 1 - i
    #right_idx = 2 - 1 - 0 = 1
    
    left_gene_id, left_pos, left_move, left_target = segment[left_idx]
    #segment[0] = ('F', 0, +5, 5)
    #left_gene_id = 'F', left_pos = 0, left_move = +5, left_target = 5
    
    right_gene_id, right_pos, right_move, right_target = segment[right_idx]
    #right_gene_id = 'A', right_pos = 1, right_move = -1, right_target = 0
    
    new_left_pos = start_index + left_idx
    #new_left_pos = 0 + 0 = 0
    
    new_right_pos = start_index + right_idx
    #new_right_pos = 0 + 1 = 1
    
    segment[left_idx] = (right_gene_id, new_left_pos, right_move, right_target)
    #segment(0) = ('A',0, -1, 0)
    
    segment[right_idx] = (left_gene_id, new_right_pos, left_move, left_target)
    #segment(1) = ('F', 1, +5, 5)
    
    """
    Recalculate movement
    
    segement(0) = ('A', 0, 0, 0)
    segment(1) = ('F', 1, 4, 5)
    
    """
    
    """
    Total movement for validation 
    
    Before = |+5|+|-1| adding magnitude = 6
    After = |0|+|4| adding magnitude. = 4
    
    Hence = valid move and correct position for adjacency inversions.
    
    """
    
    

#### continue for other valid adjacency segments in batch processing

#update sequence

all_segment = [('A', 0, 0, 0), ('F', 1, 4, 5), ('B', 2, -1, 1), 
           ('C', 3, -1, 2), ('D', 4, -1, 3), ('G', 5, +1, 6), 
           ('E', 6, -2, 4)]