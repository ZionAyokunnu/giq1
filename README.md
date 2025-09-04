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
pip install git+https://github.com/ZionAyokunnu/giq22.git
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

If you use this software in your research, please cite:

```
Ayokunnu, Z. (2024). GIQ: Gene Inversion Quantifier. 
https://github.com/ZionAyokunnu/giq2
```
