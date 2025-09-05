# GIQ - Gene Inversion Quantifier

A Python package for analyzing chromosomal rearrangements through gene inversion algorithms.

## Features

- **Gene Inversion Detection**: Identifies flip and adjacency inversions in gene sequences
- **Iterative Optimization**: Converges gene sequences to target arrangements
- **Movement Analysis**: Tracks gene movements and calculates total displacement
- **Batch Processing**: Handles multiple chromosome comparisons
- **Visualization**: ASCII-based sequence evolution tracking for very short test letter sequence

## Installation

### From GitHub
```bash
pip install git+https://github.com/ZionAyokunnu/giq2.git
```

### From PyPI (when published)
```bash
pip install giq2
```

## Quick Start

### Basic Usage
```python
from giq2.core.formats import create_movement_sequence, run_algorithm_test

# Define target and current sequences
target_sequence = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J']
current_sequence = ['H', 'I', 'J', 'A', 'B', 'F', 'E', 'D', 'C', 'G']

# Create movement sequence
movement_sequence = create_movement_sequence(target_sequence, current_sequence)

# Run algorithm
results = run_algorithm_test(movement_sequence, max_iterations=20)
print(f"Converged: {results['converged']}")
print(f"Total inversions: {len(results['inversion_events'])}")
```

```
Running main custom sequence...
==================================================
TEST ANALYSIS 
Target: A B C D E F G H I J
Current:  H I J A B F E D C G

========================

Total inversions applied: 10
Converged: True
COMPLETE SEQUENCE TRACKER
==================================================
Iter 0:  H   I─  J   A  ─B   F   E   D   C   G
         │     ── │ │ ──     │   │   │   │   │
         │       ─────       │   │   │   │   │
         │       ──│──       │   │   │   │   │
         │     ── │ │ ──     │   │   │   │   │
Iter 1:  H   B─  A   J─ ─I   F   E   D  ─C   G
         │   │   │     ─────  │ │  ─────     │
         │   │   │        ────│ │────        │
         │   │   │           ─────           │
         │   │   │           ──│──           │
         │   │   │        ────│ │────        │
         │   │   │     ─────  │ │  ─────     │
Iter 2:  H─  B  ─A   C─  D   E   F   I  ─J   G
           ─ │ ─     │   │   │   │   │   │   │
            ─│─      │   │   │   │   │   │   │
           ─ │ ─     │   │   │   │   │   │   │
Iter 3:  A─  B  ─H   C   D   E   F   I   J   G
         │   │    ╲ ╱    │   │   │   │   │   │
         │   │     ╱     │   │   │   │   │   │
         │   │    ╱ ╲    │   │   │   │   │   │
Iter 4:  A   B   C   H   D   E   F   I   J   G
         │   │   │    ╲ ╱    │   │   │   │   │
         │   │   │     ╱     │   │   │   │   │
         │   │   │    ╱ ╲    │   │   │   │   │
Iter 5:  A   B   C   D   H   E   F   I   J   G
         │   │   │   │    ╲ ╱    │   │   │   │
         │   │   │   │     ╱     │   │   │   │
         │   │   │   │    ╱ ╲    │   │   │   │
Iter 6:  A   B   C   D   E   H   F   I   J   G
         │   │   │   │   │   │   │   │    ╲ ╱
         │   │   │   │   │   │   │   │     ╱
         │   │   │   │   │   │   │   │    ╱ ╲
Iter 7:  A   B   C   D   E   H   F   I   G   J
         │   │   │   │   │    ╲ ╱    │   │   │
         │   │   │   │   │     ╱     │   │   │
         │   │   │   │   │    ╱ ╲    │   │   │
Iter 8:  A   B   C   D   E   F   H   I   G   J
         │   │   │   │   │   │   │    ╲ ╱    │
         │   │   │   │   │   │   │     ╱     │
         │   │   │   │   │   │   │    ╱ ╲    │
Iter 9:  A   B   C   D   E   F   H   G   I   J
         │   │   │   │   │   │    ╲ ╱    │   │
         │   │   │   │   │   │     ╱     │   │
         │   │   │   │   │   │    ╱ ╲    │   │
Iter 10:  A   B   C   D   E   F   G   H   I   J
```

### Command Line Interface
```bash
# Pairwise comparison
giq2 pairwise-comparison file1.tsv file2.tsv -o output_directory

# Test sequence
python -m giq2.test_sequence_example
```

## Documentation

For detailed documentation, its coming right away!

## License

This project is licensed under construction 🥳

## Contributing

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Add tests
5. Submit a pull request

## Citation

```
Supervisors: Kamil, Sasha, Arif, Sam
```

If you use this software in your research, please cite:

```
Ayokunnu, Z. (2024). GIQ: Gene Inversion Quantifier. 
https://github.com/ZionAyokunnu/giq2
```
