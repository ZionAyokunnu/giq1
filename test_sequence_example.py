#!/usr/bin/env python3

import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from core.formats import run_algorithm_test, create_movement_sequence, terminal_sequence_tracker_line

def custom_sequence():
    
    # Target sequence
    target_sequence = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J']
    
    # Current sequence 
    current_sequence = ['H', 'I', 'J', 'A', 'B', 'F', 'E', 'D', 'C', 'G']
    
    print("=" * 50)
    print("TEST ANALYSIS ")
    print(f"Target: {' '.join(target_sequence)}")
    print(f"Current:  {' '.join(current_sequence)}")
    print()
    
    # Create movement sequence
    movement_sequence = create_movement_sequence(target_sequence, current_sequence)
    
    # Run algorithm test with custom sequence
    results = run_algorithm_test(movement_sequence, max_iterations=20)
    
    # Create sequence states for visualization
    sequence_states = []
    
    # Add initial state
    sequence_states.append(movement_sequence)
    
    # Reconstruct intermediate states from inversion events
    current_seq = movement_sequence.copy()
    
    if 'inversion_events' in results:
        for event in results['inversion_events']:
            if event['type'] == 'flip':
                # Apply flip inversion
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
                        new_target = next(t for g, _, _, t in movement_sequence if g == new_gene)
                        new_move = new_target - new_pos
                        new_seq.append((new_gene, new_pos, new_move, new_target))
                    else:
                        new_seq.append((gene, pos, move, target))
                
                current_seq = new_seq
                sequence_states.append(current_seq.copy())
                
            elif event['type'] == 'adjacency':
                # Apply adjacency inversion
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
                sequence_states.append(current_seq.copy())
    
    # Add final state if different from last intermediate
    if sequence_states and sequence_states[-1] != results.get('final_sequence', []):
        sequence_states.append(results.get('final_sequence', []))
    
    # Show complete sequence tracker only
    terminal_sequence_tracker_line(sequence_states)
    
    return results

if __name__ == "__main__":
    # Run main custom sequence
    print("Running main custom sequence...")
    custom_sequence()





