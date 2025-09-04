#!/usr/bin/env python3

import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from core.formats import run_algorithm_test, create_movement_sequence

def custom_sequence():
    
    # Target sequence (linearis)
    target_sequence = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J']
    
    # Current sequence (rufipes) 
    current_sequence = ['H', 'I', 'J', 'A', 'B', 'F', 'E', 'D', 'C', 'G']
    
    print("TEST ANALYSIS ")
    print(f"Target: {' '.join(target_sequence)}")
    print(f"Current:  {' '.join(current_sequence)}")
    print()
    
    # Create movement sequence
    movement_sequence = create_movement_sequence(target_sequence, current_sequence)
    
    print("Initial movement sequence:")
    for gene, pos, move, target in movement_sequence:
        print(f"  {gene}: pos={pos}, movement={move}, target={target}")
    
    # Run algorithm test with custom sequence
    results = run_algorithm_test(movement_sequence, max_iterations=10)
    
    return results

if __name__ == "__main__":
    # Run main custom sequence
    print("Running main custom sequence...")
    custom_sequence()





