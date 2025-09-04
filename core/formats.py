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
    
    # Run the algorithm with detailed tracking
    results = iterative_detection(movement_sequence, max_iterations=max_iterations)
    
    print(f"\nTotal inversions applied: {len(results['inversion_events'])}")
    print(f"Converged: {results['converged']}")
    
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

def terminal_sequence_tracker_line(sequence_states):

    print("COMPLETE SEQUENCE TRACKER")
    print("=" * 50)
    
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