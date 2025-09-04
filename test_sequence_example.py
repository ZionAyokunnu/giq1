#!/usr/bin/env python3

import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from core.reverse import iterative_detection, detect_adjacency_inversions, detect_flip_in_pattern, find_non_overlapping_adjacencies, find_non_overlapping_flips

def create_test_sequences():
    
    # Target sequence (linearis)
    target_sequence = ['A', 'B', 'C', 'D', 'E', 'F', 'D', 'E', 'F', 'G', 'H', 'I', 'J']
    
    # Current sequence (rufipes) 
    current_sequence = ['H', 'I', 'J', 'A', 'B', 'F', 'E', 'D', 'C', 'G']
    
    print("TEST ANALYSIS ")
    print(f"Target: {' '.join(target_sequence)}")
    print(f"Current:  {' '.join(current_sequence)}")
    print()
    
    # Create movement sequence
    movement_sequence = []
    
    for i, gene in enumerate(current_sequence):
        # Find position in target sequence
        if gene in target_sequence:
            target_pos = target_sequence.index(gene)
            movement = target_pos - i  # target - current
            movement_sequence.append((gene, i, movement, target_pos))
            print(f"Gene {gene}: current_pos={i}, target_pos={target_pos}, movement={movement}")
        else:
            print(f"Gene {gene}: NOT FOUND in target sequence")
    
    print()
    print()
    print()
    print()
    print()
    print()
    print()
    print()
    print("Movement sequence:")
    
    #Space for my eyes
    
    
    
    for gene, pos, move, target in movement_sequence:
        
        print(f"  {gene}: pos={pos}, movement={move}, target={target}")
    
    return movement_sequence

def run_algorithm_test():
    
    movement_sequence = create_test_sequences()
    
    # Run the algorithm
    results = iterative_detection(movement_sequence, max_iterations=10)
    
    print("\n" + "="*50)
    print("RESULTS")
    print("="*50)
    print("="*50)
    print("="*50)
    print("="*50)
    print("="*50)
    
    final_sequence = results['final_sequence']
    print("Final sequence:")
    for gene, pos, move, target in final_sequence:
        print(f"  {gene}: pos={pos}, movement={move}, target={target}")
    
    print(f"\nTotal inversions applied: {len(results['inversion_events'])}")
    print(f"Converged: {results['converged']}")
    print(f"Final total movement: {results['final_total_movement']}")

if __name__ == "__main__":
    run_algorithm_test()
