#!/usr/bin/env python3

import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from core.formats import run_algorithm_test, create_movement_sequence

def print_gene_tracker_with_diagonals(sequence_states, tracked_gene='A'):
    
    # Track the gene's position through each iteration
    gene_positions = []
    for state in sequence_states:
        genes = [g for g, _, _, _ in state]
        if tracked_gene in genes:
            pos = genes.index(tracked_gene)
            gene_positions.append(pos)
        else:
            gene_positions.append(-1)  # Gene not found
    
    sequence_length = len(sequence_states[0])
    iterations = len(sequence_states)
    
    # Create a grid for the diagram
    grid_width = sequence_length * 4  # 4 chars per position
    grid_height = iterations * 3      # 3 rows per iteration
    
    # Initialize grid with spaces
    grid = [[' ' for _ in range(grid_width)] for _ in range(grid_height)]
    
    # Place gene positions and draw diagonal connections
    for i, pos in enumerate(gene_positions):
        if pos >= 0:  # Valid position
            row = i * 3
            col = pos * 4 + 1
            
            # Place the gene letter
            if row < grid_height and col < grid_width:
                grid[row][col] = tracked_gene
            
            # Draw diagonal line to next position
            if i < len(gene_positions) - 1 and gene_positions[i + 1] >= 0:
                next_pos = gene_positions[i + 1]
                next_col = next_pos * 4 + 1
                
                # Use Unicode diagonal characters
                if next_col > col:  # Moving right
                    # Use ╲ (U+2572) for down-right diagonal
                    if row + 1 < grid_height and col + 1 < grid_width:
                        grid[row + 1][col + 1] = '╲'
                    if row + 2 < grid_height and col + 2 < grid_width:
                        grid[row + 2][col + 2] = '╲'
                elif next_col < col:  # Moving left
                    # Use ╱ (U+2571) for down-left diagonal
                    if row + 1 < grid_height and col - 1 >= 0:
                        grid[row + 1][col - 1] = '╱'
                    if row + 2 < grid_height and col - 2 >= 0:
                        grid[row + 2][col - 2] = '╱'
                else:  # Same position
                    # Use │ (vertical line) for staying in place
                    if row + 1 < grid_height:
                        grid[row + 1][col] = '│'
                    if row + 2 < grid_height:
                        grid[row + 2][col] = '│'
    
    # Print the grid
    for i, row in enumerate(grid):
        if i % 3 == 0:  # Gene position row
            iteration = i // 3
            line = ''.join(row).rstrip()
            print(f"Iter {iteration}: {line}")
        elif i % 3 == 1 or i % 3 == 2:  # Connection rows
            line = ''.join(row).rstrip()
            if line.strip():
                print(f"        {line}")

def draw_line_dda(x0, y0, x1, y1):
    points = []
    
    # Calculate differences
    dx = x1 - x0
    dy = y1 - y0
    
    # Calculate steps needed
    steps = max(abs(dx), abs(dy))
    
    if steps == 0:
        return [(x0, y0)]
    
    # Calculate increments
    x_increment = dx / steps
    y_increment = dy / steps
    
    # Generate line points
    x, y = x0, y0
    for i in range(steps + 1):
        points.append((round(x), round(y)))
        x += x_increment
        y += y_increment
    
    return points

def print_complete_sequence_tracker(sequence_states):

    print("COMPLETE SEQUENCE TRACKER")
    print("=" * 60)
    
    iterations = len(sequence_states)
    sequence_length = len(sequence_states[0])
    
    # Track all genes through iterations
    all_genes = [g for g, _, _, _ in sequence_states[0]]
    
    # Calculate movement distances for all genes
    gene_positions = {}
    movement_distances = {}
    
    for gene in all_genes:
        positions = []
        distances = []
        for state in sequence_states:
            genes = [g for g, _, _, _ in state]
            if gene in genes:
                pos = genes.index(gene)
                positions.append(pos)
            else:
                positions.append(-1)
        
        gene_positions[gene] = positions
        
        # Calculate distances between consecutive positions
        for i in range(len(positions) - 1):
            if positions[i] >= 0 and positions[i + 1] >= 0:
                distance = abs(positions[i + 1] - positions[i])
                distances.append(distance)
            else:
                distances.append(0)
        movement_distances[gene] = distances
    
    # Calculate row allocations for each gene
    row_allocations = {}
    for gene in all_genes:
        allocations = []
        for distance in movement_distances[gene]:
            if distance == 0:
                rows = 1  # Minimal spacing for no movement
            elif distance <= 2:
                rows = 3  # Short movements get 3 rows
            elif distance <= 5:
                rows = distance + 1  # Medium movements
            else:
                rows = min(distance, 8)  # Cap at 8 rows for very long movements
            allocations.append(rows)
        row_allocations[gene] = allocations
    
    # Calculate total grid height (use maximum allocation for each iteration)
    total_rows = iterations  # Gene position rows
    for i in range(iterations - 1):
        max_allocation = 0
        for gene in all_genes:
            if i < len(row_allocations[gene]):
                max_allocation = max(max_allocation, row_allocations[gene][i])
        total_rows += max_allocation
    
    # Create dynamic grid
    grid_width = sequence_length * 4
    grid = [[' ' for _ in range(grid_width)] for _ in range(total_rows)]
    
    # Place all genes and draw connecting lines
    current_row = 0
    for i in range(iterations):
        # Place all genes at current iteration
        for gene in all_genes:
            pos = gene_positions[gene][i]
            if pos >= 0:
                col = pos * 4 + 1
                if current_row < total_rows and col < grid_width:
                    grid[current_row][col] = gene
        
        # Draw connecting lines to next iteration
        if i < iterations - 1:
            max_allocation = 0
            for gene in all_genes:
                if i < len(row_allocations[gene]):
                    max_allocation = max(max_allocation, row_allocations[gene][i])
            
            next_row = current_row + max_allocation + 1
            
            # Draw lines for each gene
            for gene in all_genes:
                current_pos = gene_positions[gene][i]
                next_pos = gene_positions[gene][i + 1]
                
                if current_pos >= 0 and next_pos >= 0:
                    current_col = current_pos * 4 + 1
                    next_col = next_pos * 4 + 1
                    
                    # Use DDA algorithm for proper line drawing
                    line_points = draw_line_dda(current_col, current_row, next_col, next_row)
                    
                    # Draw the line with appropriate characters
                    for point_col, point_row in line_points[1:-1]:  # Skip endpoints
                        if 0 <= point_row < total_rows and 0 <= point_col < grid_width:
                            # Choose character based on line direction
                            dx = next_col - current_col
                            dy = next_row - current_row
                            
                            if abs(dx) > abs(dy):  # More horizontal
                                grid[point_row][point_col] = '─'
                            elif abs(dy) > abs(dx):  # More vertical
                                grid[point_row][point_col] = '│'
                            elif dx > 0 and dy > 0:  # Down-right diagonal
                                grid[point_row][point_col] = '╲'
                            elif dx < 0 and dy > 0:  # Down-left diagonal
                                grid[point_row][point_col] = '╱'
                            else:
                                grid[point_row][point_col] = '+'
        
        # Move to next iteration row
        if i < iterations - 1:
            max_allocation = 0
            for gene in all_genes:
                if i < len(row_allocations[gene]):
                    max_allocation = max(max_allocation, row_allocations[gene][i])
            current_row += max_allocation + 1
        else:
            current_row += 1
    
    # Print the grid with iteration labels
    current_row = 0
    for i in range(iterations):
        line = ''.join(grid[current_row]).rstrip()
        print(f"Iter {i}: {line}")
        
        # Print connection lines
        if i < iterations - 1:
            max_allocation = 0
            for gene in all_genes:
                if i < len(row_allocations[gene]):
                    max_allocation = max(max_allocation, row_allocations[gene][i])
            
            for j in range(1, max_allocation + 1):
                if current_row + j < total_rows:
                    line = ''.join(grid[current_row + j]).rstrip()
                    if line.strip():
                        print(f"        {line}")
            current_row += max_allocation + 1
        else:
            current_row += 1

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
    
    # Run algorithm test with custom sequence
    results = run_algorithm_test(movement_sequence, max_iterations=5)
    
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
    
    # Show diagonal tracking for key genes
    print("\n" + "=" * 60)
    
    # Show complete sequence tracker only
    print_complete_sequence_tracker(sequence_states)
    
    return results

if __name__ == "__main__":
    # Run main custom sequence
    print("Running main custom sequence...")
    custom_sequence()





