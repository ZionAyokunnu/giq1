from core.reverse import iterative_detection

def create_movement_sequence(target_sequence, current_sequence):

    movement_sequence = []
    
    for i, gene in enumerate(current_sequence):
        # Find position in target sequence
        if gene in target_sequence:
            target_pos = target_sequence.index(gene)
            movement = target_pos - i  # target - current
            movement_sequence.append((gene, i, movement, target_pos))
        else:
            print(f"Warning: Gene {gene} NOT FOUND in target sequence")
    
    return movement_sequence


def run_algorithm_test(movement_sequence=None, max_iterations=10):
    """
    Run the algorithm test with any movement sequence.
    
    Args:
        movement_sequence (list, optional): Pre-created movement sequence
        max_iterations (int): Maximum iterations for the algorithm
    """
    
    if movement_sequence is None:
        raise ValueError("movement_sequence must be provided - no default data in formats.py")
    
    print("\n" + "="*50)
    print("SEQUENCE EVOLUTION ANALYSIS")
    print("="*50)
    
    # Display initial state
    print("\n=== INITIAL STATE ===")
    print(f"Sequence: {[gene for gene, _, _, _ in movement_sequence]}")
    print(f"Movements: {format_movement_values(movement_sequence)}")
    
    # Run the algorithm with detailed tracking
    results = iterative_detection(movement_sequence, max_iterations=max_iterations)
    
    print("\n" + "="*50)
    print("RESULTS")
    print("="*50)
    
    final_sequence = results['final_sequence']
    
    # Display final state
    print("\n=== FINAL STATE ===")
    print(f"Sequence: {[gene for gene, _, _, _ in final_sequence]}")
    print(f"Movements: {format_movement_values(final_sequence)}")
    print(f"Formatted: {format_sequence_with_movements(final_sequence)}")
    
    print(f"\nTotal inversions applied: {len(results['inversion_events'])}")
    print(f"Converged: {results['converged']}")
    print(f"Final total movement: {results['final_total_movement']}")
    
    # Analyze inversion events to show sequence evolution
    if results['inversion_events']:
        print("\n" + "="*50)
        print("SEQUENCE EVOLUTION THROUGH INVERSIONS")
        print("="*50)
        
        current_seq = movement_sequence.copy()
        
        for i, event in enumerate(results['inversion_events']):
            print(f"\n--- INVERSION {i+1} ---")
            print(f"Type: {event['type']}")
            
            if event['type'] == 'flip':
                print(f"Flip indicator: {event['flip_indicator']}")
                print(f"Segment length: {event['segment_length']}")
                print(f"Positions: {event['positions']}")
                print(f"Genes: {event['genes']}")
                print(f"Gene inversions: {event['gene_inversions']}")
                
                # Apply the inversion to get the new sequence
                start, end = event['positions'][0], event['positions'][-1]
                segment_genes = [gene for gene, _, _, _ in current_seq[start:end+1]]
                reversed_segment = segment_genes[::-1]
                
                # Create new sequence with reversed segment
                new_seq = []
                for j, (gene, pos, move, target) in enumerate(current_seq):
                    if start <= j <= end:
                        # Gene is in the flipped segment
                        segment_idx = j - start
                        new_gene = reversed_segment[segment_idx]
                        # Find new position and recalculate movement
                        new_pos = j
                        new_target = next(t for g, _, _, t in current_seq if g == new_gene)
                        new_move = new_target - new_pos
                        new_seq.append((new_gene, new_pos, new_move, new_target))
                    else:
                        new_seq.append((gene, pos, move, target))
                
                current_seq = new_seq
            elif event['type'] == 'adjacency':
                print(f"Adjacency positions: {event['positions']}")
                print(f"Genes: {event['genes']}")
                
                # For adjacency inversions, swap two adjacent genes
                pos1, pos2 = event['positions'][0], event['positions'][1]
                new_seq = current_seq.copy()
                new_seq[pos1], new_seq[pos2] = new_seq[pos2], new_seq[pos1]
                
                # Recalculate movements for swapped genes
                gene1, pos1_old, move1_old, target1 = new_seq[pos1]
                gene2, pos2_old, move2_old, target2 = new_seq[pos2]
                
                new_move1 = target1 - pos1
                new_move2 = target2 - pos2
                
                new_seq[pos1] = (gene1, pos1, new_move1, target1)
                new_seq[pos2] = (gene2, pos2, new_move2, target2)
                
                current_seq = new_seq
            
            print(f"Sequence after inversion: {[gene for gene, _, _, _ in current_seq]}")
            print(f"Movements: {format_movement_values(current_seq)}")
            print(f"Formatted: {format_sequence_with_movements(current_seq)}")
            
            # Try to identify pattern type
            movements = [move for _, _, move, _ in current_seq]
            non_zero_movements = [m for m in movements if m != 0]
            
            if len(non_zero_movements) >= 2:
                # Look for patterns
                if len(non_zero_movements) == 2 and non_zero_movements[0] == -non_zero_movements[1]:
                    print(f"Pattern: TYPE 1 - Perfect Adjacency")
                elif len(non_zero_movements) >= 3:
                    # Check for flip patterns
                    positive_moves = [m for m in non_zero_movements if m > 0]
                    negative_moves = [m for m in non_zero_movements if m < 0]
                    
                    if len(positive_moves) > 0 and len(negative_moves) > 0:
                        print(f"Pattern: TYPE 2 - Flip Pattern (Pos: {positive_moves}, Neg: {negative_moves})")
                    else:
                        print(f"Pattern: TYPE 3 - Unidirectional")
                else:
                    print(f"Pattern: TYPE 4 - Mixed")
            else:
                print(f"Pattern: Converged or Minimal")
    
    return results
                
def format_sequence_with_movements(sequence):
    """Format sequence showing genes with their movement values"""
    formatted = []
    for gene, pos, move, target in sequence:
        if move >= 0:
            formatted.append(f"{gene}(+{move})")
        else:
            formatted.append(f"{gene}({move})")
    return formatted

def format_movement_values(sequence):
    """Extract just the movement values"""
    movements = []
    for gene, pos, move, target in sequence:
        if move >= 0:
            movements.append(f"+{move}")
        else:
            movements.append(f"{move}")
    return movements