

def detect_odd_length_incremental(pattern):
    """Handle odd-length patterns with natural fulcrum."""
    fulcrum_index = len(pattern) // 2
    
    # Fulcrum must be zero
    if pattern[fulcrum_index] != 0:
        return 0
    
    # Split into left and right sections
    left_section = pattern[:fulcrum_index]
    right_section = pattern[fulcrum_index + 1:]
    
    # Must have equal non-empty sections
    if len(left_section) != len(right_section) or len(left_section) == 0:
        return 0
    
    return validate_incremental_sections(left_section, right_section)


def detect_even_length_incremental(pattern):
    """Handle even-length patterns with synthetic fulcrum."""
    mid_point = len(pattern) // 2
    
    # Split into left and right sections
    left_section = pattern[:mid_point]
    right_section = pattern[mid_point:]
    
    # Must have equal non-empty sections
    if len(left_section) != len(right_section) or len(left_section) == 0:
        return 0
    
    return validate_incremental_sections(left_section, right_section)


def validate_incremental_sections(left_section, right_section):
    """Validate that sections meet incremental pattern requirements."""
    # Get non-zero values from each section
    left_nonzero = [val for val in left_section if val != 0]
    right_nonzero = [val for val in right_section if val != 0]
    
    # Must have at least one non-zero value on each side
    if not left_nonzero or not right_nonzero:
        return 0
    
    # Endpoints (first left, last right) must be non-zero and opposite signs
    if left_section[0] == 0 or right_section[-1] == 0:
        return 0
    
    if not ((left_section[0] > 0 and right_section[-1] < 0) or 
            (left_section[0] < 0 and right_section[-1] > 0)):
        return 0
    
    # Sign separation: left section non-zeros same sign, right section opposite
    if left_section[0] > 0:
        # Left should be positive, right should be negative
        if not all(val >= 0 for val in left_section):
            return 0
        if not all(val <= 0 for val in right_section):
            return 0
    else:
        # Left should be negative, right should be positive  
        if not all(val <= 0 for val in left_section):
            return 0
        if not all(val >= 0 for val in right_section):
            return 0
    
    # Incrementalism: magnitudes should decrease toward center
    if not check_incrementalism(left_nonzero, right_nonzero):
            return 0
    
    return len(left_section)  # Return section length as flip indicator


def check_incrementalism(left_nonzero, right_nonzero):
    """Check if magnitudes decrease toward center."""
    # Left side: should decrease in magnitude from outside to center
    for i in range(len(left_nonzero) - 1):
        if abs(left_nonzero[i]) < abs(left_nonzero[i + 1]):
            return False
    
    # Right side: should decrease in magnitude from center to outside
    # So when read from outside to center, should also decrease
    right_reversed = list(reversed(right_nonzero))
    for i in range(len(right_reversed) - 1):
        if abs(right_reversed[i]) < abs(right_reversed[i + 1]):
            return False
    
    return True


def is_perfect_incremental(left_nonzero, right_nonzero):
    """
    Check if pattern is perfect incremental: +n, +(n-1), ..., +1, 0, -1, ..., -(n-1), -n
    """
    if len(left_nonzero) != len(right_nonzero):
        return False
    
    n = len(left_nonzero)
    expected_left = list(range(n, 0, -1))  # [n, n-1, ..., 1]
    expected_right = list(range(-1, -n - 1, -1))  # [-1, -2, ..., -n]
    
    return left_nonzero == expected_left and right_nonzero == expected_right


def detect_extended(movement_sequence):
    """
    Detect extended flip patterns including subsequences.
    """
    flip_patterns = []
    
    # Check all possible segments (including subsequences)
    for start_i in range(len(movement_sequence)):
        for end_j in range(start_i + 1, len(movement_sequence)):  # Minimum length 2
            # Extract movement pattern for this segment
            pattern = []
            for k in range(start_i, end_j + 1):
                _, _, movement, _ = movement_sequence[k]
                pattern.append(movement)
            
            # Test if this pattern is flippable
            flip_indicator = detect_flip_in_pattern(pattern)
            if flip_indicator > 0:
                flip_patterns.append((start_i, end_j, flip_indicator))
    
    # Sort by flip indicator (larger patterns first), then by length
    flip_patterns.sort(key=lambda x: (x[2], x[1] - x[0]), reverse=True)
    
    return flip_patterns

def detect_flip_in_pattern(pattern):
    """
    Detect flip patterns based on comprehensive classification.
    
    Returns flip_indicator > 0 if pattern is valid, 0 if not valid.
    """
    if len(pattern) < 2:
        return 0
    
    # DEBUG: Show pattern being analyzed
    print(f"  🔍 FLIP PATTERN ANALYSIS:")
    print(f"     Pattern: {pattern}")
    print(f"     Length: {len(pattern)}")
    
    # Type 1: Simple Adjacent Patterns (length 2 only)
    if len(pattern) == 2:
        if (pattern[0] > 0 and pattern[1] < 0) or (pattern[0] < 0 and pattern[1] > 0):
            return 1
        return 0
    
    # Type 2A: Odd length incremental with natural fulcrum
    if len(pattern) % 2 == 1:
        return detect_odd_length_incremental(pattern)
    
    # Type 2B: Even length incremental with synthetic fulcrum  
    else:
        return detect_even_length_incremental(pattern)
