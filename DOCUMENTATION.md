# GIQ - Gene Inversion Quantifier: Comprehensive Documentation

## Overview

GIQ (Gene Inversion Quantifier) is a Python package for analyzing chromosomal rearrangements through gene inversion algorithms. It identifies flip and adjacency inversions in gene sequences, converges gene sequences to target arrangements, and tracks gene movements to calculate total displacement.

## Architecture Overview

The GIQ package follows a modular architecture with the following main components:

```
giq/
├── main.py                 # Main CLI entry point
├── __main__.py            # Module entry point
├── __init__.py            # Package initialization
├── core/                  # Core analysis modules
│   ├── busco_processor.py    # BUSCO data processing
│   ├── convergence_analysis.py # Convergence analysis
│   ├── formats.py            # Format handling and testing
│   ├── outputs.py           # Output generation and metrics
│   ├── reverse.py           # Reverse/inversion algorithms
│   ├── rules.py             # Pattern detection rules
│   └── support.py           # Support functions
├── config/                # Configuration management
│   └── settings.py          # Configuration settings
├── contextual/            # Contextual analysis
│   └── metrics.py           # BUSCO contextual metrics
├── utils/                 # Utility functions
│   └── file_utils.py        # File operations
├── simple_dotplotter.py   # Visualization tool
└── test_sequence_example.py # Test examples
```

## Main Pipeline Flow

### 1. Entry Points

#### `main.py` - CLI Interface
**Purpose**: Main command-line interface for genome comparison
**Key Functions**:
- `main()`: CLI argument parsing and command dispatch
- `pairwise_comparison_command()`: Direct pairwise genome comparison

**Input**: Command-line arguments (genome files, output directory, parameters)
**Output**: Analysis results, reports, and visualizations

#### `__main__.py` - Module Entry
**Purpose**: Allows running the package as a module
**Function**: `main()` - Delegates to main.py

### 2. Core Analysis Pipeline

#### BUSCO Data Processing (`core/busco_processor.py`)

**`parse_busco_table(busco_path, config)`**
- **Input**: 
  - `busco_path`: String path to BUSCO TSV file
  - `config`: Dictionary with configuration settings
  - **Example Input**:
    ```python
    busco_path = "/path/to/full_table.tsv"
    config = {
        'first_species_name': 'Dioctria_linearis',
        'enable_duplicate_handling': False
    }
    ```
- **Output**: pandas DataFrame with parsed BUSCO data
  - **Example Output Structure**:
    ```python
    DataFrame({
        'busco_id': ['OZ002740.1_001', 'OZ002740.1_002', ...],
        'status': ['Complete', 'Complete', ...],
        'sequence': ['OZ002740.1', 'OZ002740.1', ...],
        'gene_start': [1000, 5000, ...],
        'gene_end': [2000, 6000, ...],
        'strand': ['+', '-', ...],
        'score': [95.2, 98.7, ...],
        'length': [1000, 1000, ...],
        'line_number': [1, 2, ...]
    })
    ```
- **Function**: Parses BUSCO TSV files, handles N/A coordinates, filters valid entries

**`filter_busco_genes(busco_df, config, quality_info=None, species_name=None)`**
- **Input**: 
  - `busco_df`: pandas DataFrame from parse_busco_table()
  - `config`: Dictionary with filtering settings
  - `quality_info`: Optional quality metrics dictionary
  - `species_name`: Optional string species name
  - **Example Input**:
    ```python
    busco_df = DataFrame({
        'busco_id': ['gene1', 'gene2', 'gene3'],
        'status': ['Complete', 'Fragmented', 'Complete'],
        'sequence': ['chr1', 'chr1', 'chr2'],
        'gene_start': [1000, 2000, 3000],
        'gene_end': [1500, 2500, 3500],
        'strand': ['+', '-', '+']
    })
    config = {'busco_status_filter': ['Complete']}
    ```
- **Output**: Filtered DataFrame containing only Complete genes
  - **Example Output Structure**:
    ```python
    DataFrame({
        'busco_id': ['gene1', 'gene3'],
        'status': ['Complete', 'Complete'],
        'sequence': ['chr1', 'chr2'],
        'gene_start': [1000, 3000],
        'gene_end': [1500, 3500],
        'strand': ['+', '+']
    })
    ```
- **Function**: Applies BUSCO status filtering (default: Complete genes only)

**`detect_flips(df1, df2, config)`**
- **Input**: 
  - `df1`: pandas DataFrame (genome 1 BUSCO data)
  - `df2`: pandas DataFrame (genome 2 BUSCO data)
  - `config`: Dictionary with analysis settings
  - **Example Input**:
    ```python
    df1 = DataFrame({
        'busco_id': ['gene1', 'gene2', 'gene3'],
        'sequence': ['chr1', 'chr1', 'chr1'],
        'gene_start': [1000, 2000, 3000],
        'gene_end': [1500, 2500, 3500],
        'strand': ['+', '+', '+']
    })
    df2 = DataFrame({
        'busco_id': ['gene1', 'gene2', 'gene3'],
        'sequence': ['chr1', 'chr1', 'chr1'],
        'gene_start': [1000, 2000, 3000],
        'gene_end': [1500, 2500, 3500],
        'strand': ['+', '-', '+']  # gene2 is flipped
    })
    ```
- **Output**: Tuple of (inversion_results_df, joined_df)
  - **Example Output Structure**:
    ```python
    # inversion_results_df
    DataFrame({
        'chr1': ['chr1'],
        'chr2': ['chr1'],
        'total_genes': [3],
        'flipped_genes': [1],
        'flip_rate': [0.33],
        'correlation': [0.95],
        'p_value': [0.01],
        'strand_consistency': [0.67],
        'inversion_type': ['Minor strand differences (<50%)'],
        'single_gene_inversions': [1],
        'multi_gene_inversions': [0],
        'largest_inversion': [1],
        'total_inversion_events': [1]
    })
    
    # joined_df
    DataFrame({
        'busco_id': ['gene1', 'gene2', 'gene3'],
        'chr1': ['chr1', 'chr1', 'chr1'],
        'chr2': ['chr1', 'chr1', 'chr1'],
        'start1': [1000, 2000, 3000],
        'end1': [1500, 2500, 3500],
        'start2': [1000, 2000, 3000],
        'end2': [1500, 2500, 3500],
        'strand1': ['+', '+', '+'],
        'strand2': ['+', '-', '+'],
        'is_flipped': [False, True, False]
    })
    ```
- **Function**: Detects strand flips between genomes, calculates correlations, identifies inversion events

**`correct_chromosome_orientation(df1, df2, joined_df, config)`**
- **Input**: Two BUSCO DataFrames, joined DataFrame, configuration
- **Output**: Tuple of corrected DataFrames
- **Function**: Corrects chromosome orientation based on gene position correlations

**`track_strand_changes_per_iteration(initial_genome_df, final_converged_df, inversion_events)`**
- **Input**: Initial genome DataFrame, converged genome DataFrame, inversion events list
- **Output**: Dictionary with strand tracking and validation results
- **Function**: Tracks strand changes through iterations and validates against biological expectations

**`generate_strand_debug_tsv(initial_genome_df, inversion_events, output_path, genome1_name, genome2_name)`**
- **Input**: Initial genome DataFrame, inversion events, output path, genome names
- **Output**: Tuple of (debug_path, debug_df)
- **Function**: Generates detailed TSV showing strand changes for each gene through iterations

#### Movement Sequence Creation (`core/reverse.py`)

**`create_pairwise_movement_sequence_per_chromosome(genome1_df, genome2_df, config)`**
- **Input**: 
  - `genome1_df`: pandas DataFrame (genome 1 BUSCO data)
  - `genome2_df`: pandas DataFrame (genome 2 BUSCO data)
  - `config`: Dictionary with analysis settings
  - **Example Input**:
    ```python
    genome1_df = DataFrame({
        'busco_id': ['gene1', 'gene2', 'gene3', 'gene4'],
        'sequence': ['chr1', 'chr1', 'chr2', 'chr2'],
        'gene_start': [1000, 2000, 1000, 2000],
        'gene_end': [1500, 2500, 1500, 2500]
    })
    genome2_df = DataFrame({
        'busco_id': ['gene1', 'gene2', 'gene3', 'gene4'],
        'sequence': ['chr1', 'chr1', 'chr2', 'chr2'],
        'gene_start': [2000, 1000, 2000, 1000],  # Different order
        'gene_end': [2500, 1500, 2500, 1500]
    })
    config = {
        'use_genomic_positions': False,
        'shared_genes_threshold': 2
    }
    ```
- **Output**: Dictionary mapping chromosome pairs to movement sequences
  - **Example Output Structure**:
    ```python
    {
        'chr1_vs_chr1': [
            ('gene1', 0, 1, 1),   # (gene_id, source_pos, movement, target_pos)
            ('gene2', 1, -1, 0)
        ],
        'chr2_vs_chr2': [
            ('gene3', 0, 1, 1),
            ('gene4', 1, -1, 0)
        ]
    }
    ```
- **Function**: Creates movement sequences per chromosome pair, handles chromosome pairing based on gene overlap

**`extract_movement_sequence(movement_results)`**
- **Input**: Movement analysis results dictionary
- **Output**: List of (gene_id, current_position, mean_movement) tuples
- **Function**: Extracts ordered movement sequence from analysis results

**`calculate_total_movement(sequence)`**
- **Input**: Movement sequence list
- **Output**: Total absolute movement value
- **Function**: Calculates total absolute movement for a sequence

#### Inversion Detection and Application (`core/reverse.py`)

**`iterative_detection(movement_sequence, max_iterations=1000)`**
- **Input**: 
  - `movement_sequence`: List of (gene_id, position, movement, target_position) tuples
  - `max_iterations`: Integer maximum iterations (default: 1000)
  - **Example Input**:
    ```python
    movement_sequence = [
        ('gene1', 0, 2, 2),    # gene1 needs to move +2 positions
        ('gene2', 1, -1, 0),   # gene2 needs to move -1 position
        ('gene3', 2, 0, 2),   # gene3 is already in correct position
        ('gene4', 3, -2, 1)   # gene4 needs to move -2 positions
    ]
    max_iterations = 100
    ```
- **Output**: Dictionary with comprehensive analysis results
  - **Example Output Structure**:
    ```python
    {
        'final_sequence': [
            ('gene1', 2, 0, 2),   # Final converged positions
            ('gene2', 0, 0, 0),
            ('gene3', 2, 0, 2),
            ('gene4', 1, 0, 1)
        ],
        'inversion_events': [
            {
                'type': 'flip',
                'iteration': 1,
                'flip_indicator': 2,
                'segment_length': 4,
                'positions': [0, 1, 2, 3],
                'genes': ['gene1', 'gene2', 'gene3', 'gene4'],
                'gene_inversions': 4
            }
        ],
        'iterations': 1,
        'total_events': 1,
        'total_gene_inversions': 4,
        'adjacency_events': 0,
        'flip_events': 1,
        'converged': True,
        'final_total_movement': 0,
        'final_non_zero_movements': 0,
        'final_large_movements': 0
    }
    ```
- **Function**: Main iterative algorithm that detects and applies inversions until convergence

**`apply_flip_inversion(movement_sequence, start_index, end_index, flip_indicator)`**
- **Input**: Movement sequence, start/end indices, flip indicator
- **Output**: Tuple of (updated_sequence, inversion_record)
- **Function**: Applies flip inversion to a segment of genes

**`apply_adjacency_inversion(movement_sequence, index1, index2)`**
- **Input**: Movement sequence, two adjacent indices
- **Output**: Tuple of (updated_sequence, inversion_record)
- **Function**: Swaps two adjacent genes with opposite movement signs

#### Pattern Detection Rules (`core/rules.py`)

**`detect_adjacency_inversions(movement_sequence)`**
- **Input**: 
  - `movement_sequence`: List of (gene_id, position, movement, target_position) tuples
  - **Example Input**:
    ```python
    movement_sequence = [
        ('gene1', 0, 1, 1),    # Positive movement
        ('gene2', 1, -1, 0),   # Negative movement (adjacent opposite)
        ('gene3', 2, 2, 4),    # Positive movement
        ('gene4', 3, 1, 4)     # Positive movement
    ]
    ```
- **Output**: List of (index1, index2) tuples for adjacent pairs with opposite signs
  - **Example Output Structure**:
    ```python
    [(0, 1)]  # gene1 and gene2 are adjacent with opposite movement signs
    ```
- **Function**: Detects +/- adjacencies in the movement sequence

**`detect_extended(movement_sequence)`**
- **Input**: 
  - `movement_sequence`: List of (gene_id, position, movement, target_position) tuples
  - **Example Input**:
    ```python
    movement_sequence = [
        ('gene1', 0, 3, 3),    # Pattern: +3, +2, +1, 0, -1, -2, -3
        ('gene2', 1, 2, 3),
        ('gene3', 2, 1, 3),
        ('gene4', 3, 0, 3),
        ('gene5', 4, -1, 3),
        ('gene6', 5, -2, 3),
        ('gene7', 6, -3, 3)
    ]
    ```
- **Output**: List of (start_i, end_j, flip_indicator) tuples
  - **Example Output Structure**:
    ```python
    [(0, 6, 3)]  # Flip pattern from index 0 to 6 with flip indicator 3
    ```
- **Function**: Detects extended flip patterns including subsequences

**`detect_flip_in_pattern(pattern)`**
- **Input**: 
  - `pattern`: List of movement values
  - **Example Input**:
    ```python
    # Odd-length incremental pattern
    pattern = [3, 2, 1, 0, -1, -2, -3]  # Perfect incremental with natural fulcrum
    
    # Even-length incremental pattern  
    pattern = [2, 1, -1, -2]  # Perfect incremental with synthetic fulcrum
    
    # Simple adjacent pattern
    pattern = [1, -1]  # Simple swap
    ```
- **Output**: Flip indicator (>0 if valid, 0 if not)
  - **Example Output Structure**:
    ```python
    3  # Valid flip pattern with flip indicator 3
    2  # Valid flip pattern with flip indicator 2  
    1  # Valid flip pattern with flip indicator 1
    0  # Invalid pattern
    ```
- **Function**: Detects flip patterns based on comprehensive classification

**`detect_odd_length_incremental(pattern)`**
- **Input**: Movement pattern
- **Output**: Flip indicator
- **Function**: Handles odd-length patterns with natural fulcrum

**`detect_even_length_incremental(pattern)`**
- **Input**: Movement pattern
- **Output**: Flip indicator
- **Function**: Handles even-length patterns with synthetic fulcrum

**`validate_incremental_sections(left_section, right_section)`**
- **Input**: Left and right sections of pattern
- **Output**: Flip indicator
- **Function**: Validates that sections meet incremental pattern requirements

**`check_incrementalism(left_nonzero, right_nonzero)`**
- **Input**: Non-zero values from left and right sections
- **Output**: Boolean
- **Function**: Checks if magnitudes decrease toward center

**`is_perfect_incremental(left_nonzero, right_nonzero)`**
- **Input**: Non-zero values from left and right sections
- **Output**: Boolean
- **Function**: Checks if pattern is perfect incremental: +n, +(n-1), ..., +1, 0, -1, ..., -(n-1), -n

#### Support Functions (`core/support.py`)

**`find_non_overlapping_adjacencies(adjacency_inversions, current_sequence)`**
- **Input**: Adjacency inversions list, current sequence
- **Output**: List of non-overlapping adjacent pairs
- **Function**: Finds adjacencies that are truly adjacent AND don't overlap

**`find_non_overlapping_flips(flip_patterns)`**
- **Input**: Flip patterns list
- **Output**: List of non-overlapping flips
- **Function**: Finds flip patterns that don't overlap in their index ranges

**`apply_batch_with_sequential_fallback(current_sequence, batch_operations, operation_type, target_positions)`**
- **Input**: Current sequence, batch operations, operation type, target positions
- **Output**: Tuple of (updated_sequence, all_records)
- **Function**: Applies batch operations with sequential fallback for interdependent operations

**`validate_segment_independence(current_sequence, batch_operations, operation_type)`**
- **Input**: Current sequence, batch operations, operation type
- **Output**: Boolean indicating independence
- **Function**: Validates that batch operations have independent movement segments

**`apply_concurrent_batch_flips(current_sequence, batch_operations)`**
- **Input**: Current sequence, batch operations
- **Output**: Tuple of (updated_sequence, all_records)
- **Function**: Applies multiple flip operations concurrently when independent

**`apply_concurrent_batch_adjacencies(current_sequence, batch_operations)`**
- **Input**: Current sequence, batch operations
- **Output**: Tuple of (updated_sequence, all_records)
- **Function**: Applies multiple adjacency operations concurrently when independent

### 3. Output Generation (`core/outputs.py`)

**`calculate_comprehensive_metrics(results, start_time=None, end_time=None)`**
- **Input**: 
  - `results`: Dictionary with analysis results from iterative_detection()
  - `start_time`: Optional float timestamp (default: None)
  - `end_time`: Optional float timestamp (default: None)
  - **Example Input**:
    ```python
    results = {
        'distance_metrics': {
            'total_events': 5,
            'total_gene_inversions': 12,
            'iterations': 3,
            'chromosome_results': {
                'chr1_vs_chr1': {
                    'total_events': 3,
                    'total_gene_inversions': 8,
                    'iterations': 2,
                    'converged': True,
                    'final_total_movement': 0,
                    'final_non_zero_movements': 0
                }
            }
        }
    }
    start_time = 1234567890.0
    end_time = 1234567895.0
    ```
- **Output**: Comprehensive metrics dictionary
  - **Example Output Structure**:
    ```python
    {
        'summary': {
            'total_inversion_events': 5,
            'total_genes_rearranged': 12,
            'total_genes_processed': 50,
            'algorithm_iterations': 3,
            'chromosomes_processed': 1,
            'convergence_rate': '100.0%',
            'overall_movement_reduction': '95.2%'
        },
        'event_breakdown': {
            'adjacency_inversions': 2,
            'segment_flips': 3,
            'avg_genes_per_event': '2.4',
            'events_per_iteration': '1.7'
        },
        'efficiency': {
            'avg_events_per_chromosome': '5.0',
            'rearrangement_density': '24.0%',
            'convergence_success_rate': '100.0%',
            'movement_reduction_efficiency': '19.0 units/event'
        },
        'biological_significance': {
            'perfectly_converged_chromosomes': 1,
            'synteny_preserved_chromosomes': 0,
            'highly_rearranged_chromosomes': 1,
            'total_movement_eliminated': 95,
            'remaining_movement': 5
        },
        'performance': {
            'execution_time_seconds': 5.0,
            'genes_per_second': 10.0,
            'events_per_second': 1.0,
            'chromosomes_per_second': 0.2
        },
        'chromosome_details': {
            'chr1_vs_chr1': {
                'genes_processed': 50,
                'events_applied': 3,
                'genes_rearranged': 8,
                'iterations_needed': 2,
                'convergence_achieved': True,
                'initial_movement': 100,
                'final_movement': 5,
                'movement_reduction': 95,
                'reduction_percentage': 95.0,
                'remaining_unresolved': 2
            }
        }
    }
    ```
- **Function**: Calculates detailed metrics for genomic rearrangement analysis

**`print_detailed_metrics(results, start_time=None, end_time=None)`**
- **Input**: Analysis results, optional timing information
- **Output**: None (prints to console)
- **Function**: Prints comprehensive analysis metrics in formatted output

**`generate_iteration_report(results, output_path, genome1_name, genome2_name)`**
- **Input**: Analysis results, output path, genome names
- **Output**: Report file path
- **Function**: Generates detailed Markdown report of all iterations and flip patterns

**`save_pattern_analysis(results, output_path, genome1_name, genome2_name)`**
- **Input**: Analysis results, output path, genome names
- **Output**: None (saves CSV file)
- **Function**: Saves detailed pattern analysis as CSV

**`save_inversion_events(events, output_path, comparison_name)`**
- **Input**: Inversion events list, output path, comparison name
- **Output**: None (saves CSV file)
- **Function**: Saves inversion events to CSV with iteration details

**`save_converged_genome(final_sequence, genome1_df, genome2_df, chromosome_pair, output_path)`**
- **Input**: Final sequence, genome DataFrames, chromosome pair, output path
- **Output**: None (saves TSV file)
- **Function**: Saves the final converged genome sequence in original format

**`get_converged_genome_data(final_sequence, genome1_df, genome2_df, chromosome_pair)`**
- **Input**: Final sequence, genome DataFrames, chromosome pair
- **Output**: List of converged data rows
- **Function**: Gets converged genome data for a single chromosome without saving

**`save_combined_converged_genome(all_converged_data, output_path, genome1_name, genome2_name)`**
- **Input**: All converged data, output path, genome names
- **Output**: None (saves TSV file)
- **Function**: Saves all converged genome data together in a single TSV file

**`save_stage_data(data, stage_name, output_dir, description)`**
- **Input**: Data (DataFrame or dict), stage name, output directory, description
- **Output**: None (saves file)
- **Function**: Saves data from each stage as CSV or JSON with description

### 4. Convergence Analysis (`core/convergence_analysis.py`)

**`create_single_convergence_tsv(genome1_df, converged_df, movement_sequences, all_results, output_path, genome1_name, genome2_name)`**
- **Input**: Genome DataFrames, movement sequences, results, output path, genome names
- **Output**: Analysis DataFrame
- **Function**: Creates single TSV with convergence analysis format
- **Output Format**: gene_id, {genome1_name}_convergence_positions, {genome2_name}_target_position, convergence_status, chr_concatenated, current_movement_value

**`create_convergence_analysis_tsv(converged_df, linearis_df, output_path, genome1_name, genome2_name)`**
- **Input**: Converged DataFrame, linearis DataFrame, output path, genome names
- **Output**: Analysis DataFrame
- **Function**: Creates detailed convergence analysis TSV comparing converged vs original genomes

**`analyze_transposition_patterns(non_converged_df)`**
- **Input**: Non-converged genes DataFrame
- **Output**: Dictionary with transposition patterns
- **Function**: Analyzes patterns in non-converged genes to identify potential transposition candidates

### 5. Format Handling (`core/formats.py`)

**`create_movement_sequence(target_sequence, current_sequence)`**
- **Input**: 
  - `target_sequence`: List of gene IDs in target order
  - `current_sequence`: List of gene IDs in current order
  - **Example Input**:
    ```python
    target_sequence = ['A', 'B', 'C', 'D', 'E']
    current_sequence = ['C', 'A', 'E', 'B', 'D']
    ```
- **Output**: Movement sequence list
  - **Example Output Structure**:
    ```python
    [
        ('C', 0, 2, 2),  # C is at position 0, needs to move +2 to position 2
        ('A', 1, -1, 0), # A is at position 1, needs to move -1 to position 0
        ('E', 2, 2, 4),  # E is at position 2, needs to move +2 to position 4
        ('B', 3, -2, 1), # B is at position 3, needs to move -2 to position 1
        ('D', 4, -2, 3)  # D is at position 4, needs to move -2 to position 3
    ]
    ```
- **Function**: Creates movement sequence from target and current sequences

**`run_algorithm_test(movement_sequence=None, max_iterations=10)`**
- **Input**: Movement sequence, maximum iterations
- **Output**: Analysis results dictionary
- **Function**: Runs the algorithm test with any movement sequence

**`terminal_sequence_tracker_line(sequence_states)`**
- **Input**: List of sequence states
- **Output**: None (prints to console)
- **Function**: Creates ASCII-based sequence evolution tracking visualization

**`format_sequence_with_movements(sequence)`**
- **Input**: Movement sequence
- **Output**: List of formatted strings
- **Function**: Formats sequence showing genes with their movement values

**`format_movement_values(sequence)`**
- **Input**: Movement sequence
- **Output**: List of movement value strings
- **Function**: Extracts just the movement values

**`draw_line_dda(x0, y0, x1, y1)`**
- **Input**: Start and end coordinates
- **Output**: List of line points
- **Function**: Digital Differential Analyzer algorithm for line drawing

### 6. Configuration Management (`config/settings.py`)

**CONFIG Dictionary**
- **Purpose**: Central configuration for all GIQ operations
- **Example Structure**:
  ```python
  CONFIG = {
      # File paths
      'first_busco_path': 'template/full_table.tsv',
      'second_busco_path': 'template/full_table.tsv',
      'first_species_name': 'Dioctria_linearis',
      'second_species_name': 'Dioctria_rufipes',
      
      # BUSCO filtering
      'busco_status_filter': ['Complete'],
      'enable_duplicate_handling': False,
      
      # Output settings
      'base_output_dir': 'results',
      'inversion_summary_csv': 'inversion_summary.csv',
      
      # Algorithm parameters
      'max_iterations': 1000,
      'flexible_threshold': 0.0,
      'shared_genes_threshold': 300,
      'use_genomic_positions': False,
      
      # Visualization
      'plot_width': 15,
      'plot_height': 10,
      'dpi': 300,
      'generate_dotplots': True,
      
      # Analysis parameters
      'position_bin_size_kb': 1.5,
      'probability_threshold_for_target': 0.3
  }
  ```

### 7. Contextual Analysis (`contextual/metrics.py`)

**`compute_inversion_rate_per_mb_busco(inversion_df, joined_df)`**
- **Input**: Inversion DataFrame, joined DataFrame
- **Output**: Dictionary with rate metrics
- **Function**: Computes inversion rate per megabase using BUSCO data
- **Metrics**: Total inversions, flipped genes, genome coverage, rates by type and chromosome

**`_analyse_gene_density_correlation(inversion_df, joined_df, config)`**
- **Input**: Inversion DataFrame, joined DataFrame, configuration
- **Output**: Dictionary with correlation analysis
- **Function**: Analyzes correlation between gene density and inversion density

### 8. Utility Functions (`utils/file_utils.py`)

**`create_output_directory(config)`**
- **Input**: Configuration dictionary
- **Output**: Base directory path
- **Function**: Creates output directory structure with subdirectories (plots, data, reports, debug, cache)

### 9. Visualization (`simple_dotplotter.py`)

**`load_busco_file(file_path, name)`**
- **Input**: File path, genome name
- **Output**: BUSCO DataFrame
- **Function**: Loads BUSCO file handling both formats (with/without comments)

**`create_simple_busco_dotplot(genome1_df, genome2_df, genome1_name, genome2_name, output_path, chr_order1=None, chr_order2=None)`**
- **Input**: Genome DataFrames, names, output path, optional chromosome orders
- **Output**: Correlation coefficient
- **Function**: Creates simple diagonal dotplot using linearized genomic coordinates

**`linearize_real_genome_coordinates(normal_df, custom_order=None)`**
- **Input**: Normal DataFrame, optional custom chromosome order
- **Output**: Linearized DataFrame
- **Function**: Linearizes using real genomic coordinates only

**`create_comparison_summary(genome1_df, genome2_df, genome1_name, genome2_name, correlation, output_dir)`**
- **Input**: Genome DataFrames, names, correlation, output directory
- **Output**: Summary file path
- **Function**: Creates text summary of the comparison

### 10. Testing (`test_sequence_example.py`)

**`custom_sequence()`**
- **Input**: None
- **Output**: Analysis results
- **Function**: Runs test analysis with predefined sequences
- **Test Data**: 
  - Target: ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J']
  - Current: ['H', 'I', 'J', 'A', 'B', 'F', 'E', 'D', 'C', 'G']

## Data Flow Pipeline

### 1. Input Processing
1. **BUSCO File Parsing**: `parse_busco_table()` reads and parses BUSCO TSV files
2. **Gene Filtering**: `filter_busco_genes()` filters to Complete genes only
3. **Chromosome Pairing**: `create_pairwise_movement_sequence_per_chromosome()` pairs chromosomes based on gene overlap

### 2. Movement Analysis
1. **Movement Sequence Creation**: Creates (gene_id, source_pos, movement, target_pos) tuples
2. **Pattern Detection**: `detect_extended()` and `detect_adjacency_inversions()` identify inversion patterns
3. **Iterative Application**: `iterative_detection()` applies inversions until convergence

### 3. Convergence Analysis
1. **Convergence Tracking**: Monitors total movement reduction and convergence status
2. **Strand Validation**: `track_strand_changes_per_iteration()` validates biological expectations
3. **Pattern Analysis**: Analyzes inversion patterns and types

### 4. Output Generation
1. **Metrics Calculation**: `calculate_comprehensive_metrics()` computes detailed statistics
2. **Report Generation**: Creates Markdown reports and CSV analyses
3. **Visualization**: Generates dotplots and sequence tracking visualizations

## Key Algorithm Components

### Inversion Types
1. **Adjacency Inversions**: Simple swaps of adjacent genes with opposite movement signs
2. **Flip Inversions**: Multi-gene segment reversals with specific patterns:
   - Simple Adjacent (length 2)
   - Odd-length Incremental (natural fulcrum)
   - Even-length Incremental (synthetic fulcrum)

### Convergence Criteria
1. **Perfect Convergence**: All genes have zero movement
2. **Partial Convergence**: Few large movements remain (<10 genes with |movement| > 2)
3. **Iteration Limit**: Maximum iterations reached
4. **Cycle Detection**: Algorithm stuck in repeating patterns

### Pattern Detection Rules
1. **Incremental Patterns**: Magnitudes decrease toward center
2. **Sign Separation**: Left and right sections have opposite signs
3. **Fulcrum Requirements**: Natural fulcrum must be zero, synthetic fulcrum splits evenly

## Usage Examples

### Command Line Usage
```bash
# Pairwise comparison
giq2 pairwise-comparison file1.tsv file2.tsv -o output_directory

# With custom parameters
giq2 pairwise-comparison file1.tsv file2.tsv -o output_directory --iterations 5 --mode position --shared-genes 200
```

### Python API Usage
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

### Test Sequence Example
```bash
# Run test sequence
python -m giq2.test_sequence_example
```

## Output Files

### Main Outputs
1. **Stages Directory**: Intermediate analysis files
   - `1_parsed_{genome_name}.csv`: Parsed BUSCO data
   - `2_movement_sequence_{chromosome_pair}.csv`: Movement sequences
   - `3_distance_metrics.json`: Summary metrics
   - `4_inversion_events_{chromosome_pair}.csv`: Inversion events
   - `5_converged_genome_{genome1}_vs_{genome2}.tsv`: Converged genome

2. **Analysis Reports**:
   - `iteration_report_{genome1}_vs_{genome2}.md`: Detailed iteration report
   - `pattern_analysis_{genome1}_vs_{genome2}.csv`: Pattern analysis
   - `convergence_analysis_{genome1}_vs_{genome2}.tsv`: Convergence analysis
   - `strand_{genome1}_vs_{genome2}.tsv`: Strand change tracking

3. **Results Summary**:
   - `pairwise_results_{genome1}_vs_{genome2}.json`: Complete results

## Dependencies

### Required Packages
- `pandas>=1.3.0`: Data manipulation and analysis
- `numpy>=1.20.0`: Numerical computations
- `matplotlib`: Visualization (for dotplotter)
- `scipy`: Statistical functions (Pearson correlation)

### Optional Packages
- `pytest>=6.0`: Testing
- `black>=21.0`: Code formatting
- `flake8>=3.8`: Linting

## Configuration Options

### Algorithm Parameters
- `max_iterations`: Maximum algorithm iterations (default: 1000)
- `flexible_threshold`: Allow movement increases up to this value (default: 0.0 = strict)
- `shared_genes_threshold`: Minimum shared genes for chromosome pairing (default: 300)
- `use_genomic_positions`: Use gene ranks or genomic positions (default: rank)

### Output Options
- `base_output_dir`: Base output directory (default: 'results')
- `generate_dotplots`: Generate visualization plots (default: True)
- `plot_width`, `plot_height`: Plot dimensions
- `dpi`: Plot resolution

### BUSCO Filtering
- `busco_status_filter`: BUSCO statuses to include (default: ['Complete'])
- `enable_duplicate_handling`: Handle duplicate entries (default: False)

## Performance Considerations

### Optimization Features
1. **Batch Processing**: Concurrent application of independent operations
2. **Sequential Fallback**: Sequential processing when operations conflict
3. **Early Termination**: Stop when convergence criteria met
4. **Cycle Detection**: Prevent infinite loops

### Memory Management
1. **Lazy Loading**: Load data only when needed
2. **Efficient Data Structures**: Use pandas DataFrames for large datasets
3. **Garbage Collection**: Clean up intermediate results

## Error Handling

### Common Issues
1. **No Common Genes**: Check BUSCO file formats and gene IDs
2. **Parsing Errors**: Verify TSV format and column structure
3. **Memory Issues**: Reduce dataset size or increase system memory
4. **Convergence Issues**: Adjust parameters or check data quality

### Debugging Features
1. **Verbose Logging**: Detailed progress information
2. **Intermediate Outputs**: Save data at each stage
3. **Strand Validation**: Biological consistency checks
4. **Pattern Analysis**: Detailed inversion pattern reporting

## Future Enhancements

### Planned Features
1. **Translocation Detection**: Handle inter-chromosomal rearrangements
2. **Multi-Genome Analysis**: Compare multiple genomes simultaneously
3. **Advanced Visualization**: Interactive plots and animations
4. **Machine Learning**: Pattern recognition and prediction

### Extensibility
1. **Plugin Architecture**: Custom pattern detectors
2. **API Extensions**: Additional analysis modules
3. **Format Support**: Additional input/output formats
4. **Integration**: Bioinformatics pipeline integration

---

*This documentation provides a comprehensive overview of the GIQ package. For specific implementation details, refer to the source code and inline documentation.*
