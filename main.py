#!/usr/bin/env python3
"""
GIQ2 - Genome Inversion Quantifier
Main entry point.

Author: Zion Ayokunnu
Supervisors: Kamil, Sasha, Arif, Sam
Version: 1.0.0


"""

import sys
import argparse
import logging
import os
from pathlib import Path
import pandas as pd
import json
from typing import Dict, List
import logging

logger = logging.getLogger(__name__)

from config.settings import CONFIG
from utils.file_utils import create_output_directory
from core import (
    parse_busco_table,
    filter_busco_genes,
    create_pairwise_movement_sequence_per_chromosome,
    iterative_detection,
    save_stage_data,
    print_detailed_metrics,
    generate_iteration_report,
    save_pattern_analysis,
    track_strand_changes_per_iteration,
    generate_strand_debug_tsv,
    create_single_convergence_tsv,
    save_inversion_events,
    save_combined_converged_genome,
    get_converged_genome_data
    
)



def pairwise_comparison_command(busco_file1: str, busco_file2: str, output_dir: str, config_overrides: Dict = None):
    """Direct pairwise genome comparison using iterative movement analysis."""
    import time
    
    # Start timing
    start_time = time.time()
    
    logger.info("=" * 60)
    logger.info("PAIRWISE GENOME COMPARISON")
    logger.info("=" * 60)
    
    config = CONFIG.copy()
    if config_overrides:
        config.update(config_overrides)
    
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # Extract genome names
    genome1_name = Path(busco_file1).stem
    genome2_name = Path(busco_file2).stem
    
    logger.info(f"Genome 1: {genome1_name}")
    logger.info(f"Genome 2: {genome2_name}")
    
    # Parse and process both genomes
    genome1_df = parse_busco_table(busco_file1, config)
    genome1_filtered = filter_busco_genes(genome1_df, config)
    
    genome2_df = parse_busco_table(busco_file2, config)
    genome2_filtered = filter_busco_genes(genome2_df, config)

    
    save_stage_data(genome1_filtered, f'1_parsed_{genome1_name}', output_path, f"Parsed genome 1: {genome1_name}")
    save_stage_data(genome2_filtered, f'1_parsed_{genome2_name}', output_path, f"Parsed genome 2: {genome2_name}")
    
    # Create movement sequence directly from position differences
    chromosome_sequences = create_pairwise_movement_sequence_per_chromosome(genome1_filtered, genome2_filtered, config)
    

        # Run iterative detection on each chromosome pair
    all_results = {}
    total_events = 0
    total_gene_inversions = 0
    all_converged_data = []  # Collect all converged data
    
    
    for chromosome_pair, movement_sequence in chromosome_sequences.items():
        
        # Save chromosome-specific movement sequence
        save_stage_data(pd.DataFrame(movement_sequence, columns=['gene_id', 'source_pos', 'movement', 'target_pos']), 
                       f'2_movement_sequence_{chromosome_pair}', output_path, 
                       f"Movement sequence for {chromosome_pair}: {len(movement_sequence)} genes")

            

        # Run iterative detection on this chromosome
        chr_analysis = iterative_detection(movement_sequence, config.get('max_iterations', 1000))
        all_results[chromosome_pair] = chr_analysis
        
        
        # Accumulate totals
        total_events += chr_analysis['total_events']
        total_gene_inversions += chr_analysis['total_gene_inversions']
        
        # Save chromosome-specific events
        save_inversion_events(chr_analysis['inversion_events'], output_path, chromosome_pair)
        
        # Collect converged genome data (don't save individually)
        chr_converged_data = get_converged_genome_data(chr_analysis['final_sequence'], genome1_filtered, genome2_filtered, 
                                                      chromosome_pair)
        all_converged_data.extend(chr_converged_data)
    
    # Create combined distance metrics from per-chromosome results
    distance_metrics = {
        'total_events': total_events,
        'total_gene_inversions': total_gene_inversions,
        'chromosome_results': all_results,
        'iterations': max([result['iterations'] for result in all_results.values()]) if all_results else 0,
        'converged': all([result['converged'] for result in all_results.values()]) if all_results else True
    }
    
    # Save combined results summary
    save_stage_data(distance_metrics, '3_distance_metrics', output_path, "Distance metrics summary")
    
        # Save combined converged genome
    converged_df = None  # Initialize converged_df
    if all_converged_data:
        save_combined_converged_genome(all_converged_data, output_path, genome1_name, genome2_name)
        
        
        # Create detailed convergence analysis
        converged_df = pd.DataFrame(all_converged_data)
        create_single_convergence_tsv(genome1_df, converged_df, chromosome_sequences, all_results, output_path, genome1_name, genome2_name)
    
    # Initialize results (moved outside the conditional block)
    results = {
        'genome1': genome1_name,
        'genome2': genome2_name,
        'distance_metrics': distance_metrics,
        'chromosome_results': all_results,
        'timing': {
            'start_time': start_time,
            'end_time': time.time()
        }
    }
    
    # Perform strand tracking validation
    try:
        # Collect all inversion events from all chromosomes
        all_inversion_events = []
        for chr_name, chr_data in all_results.items():
            if 'inversion_events' in chr_data:
                all_inversion_events.extend(chr_data['inversion_events'])
        
        # Perform strand tracking validation only if converged_df exists
        if converged_df is not None:
            strand_validation = track_strand_changes_per_iteration(genome1_df, converged_df, all_inversion_events)
            
            # Add strand validation results to the main results
            results['strand_validation'] = strand_validation
            
            # Print strand validation summary
            validation = strand_validation['validation_results']

            # Generate detailed strand debug TSV
            strand_path, debug_df = generate_strand_debug_tsv(
                genome2_filtered,  # Converged genome original strand i.e baseline strand
                all_inversion_events,
                output_path,
                genome1_name,
                genome2_name
            )
            print(f"Strand TSV saved: {strand_path}")
            
    except Exception as e:
        print(f"Warning: Strand tracking validation failed: {e}")
        results['strand_validation'] = None
    
    results_file = output_path / f'pairwise_results_{genome1_name}_vs_{genome2_name}.json'
    with open(results_file, 'w') as f:
        json.dump(results, f, indent=2, default=str)
    
    logger.info(f"Results saved to: {results_file}")
    return results





def main():
    """Main CLI entry point"""
    parser = argparse.ArgumentParser(description='Genome Inversion Quantification Tools')
    subparsers = parser.add_subparsers(dest='command', help='Available commands')
    
    # Pairwise comparison command
    pairwise_parser = subparsers.add_parser('pairwise-comparison', help='Direct pairwise genome comparison')
    pairwise_parser.add_argument('genome1_busco', help='First BUSCO table file')
    pairwise_parser.add_argument('genome2_busco', help='Second BUSCO table file (Target)') 
    pairwise_parser.add_argument('-o', '--output', required=True, help='Output directory')
    pairwise_parser.add_argument('--iterations', type=int, default=3, help='Max iterations (default: 3)')
    pairwise_parser.add_argument('--mode', choices=['rank', 'position'], default='rank', help='Use gene ranks or genomic positions (default: rank)')
    pairwise_parser.add_argument('--flexible', type=float, default=0.0, help='Allow movement increases up to this value (default: 0.0 = strict)')
    pairwise_parser.add_argument('--shared-genes', type=int, default=300, help='Minimum shared genes required for chromosome pairing (default: 300)')
    pairwise_parser.add_argument('--bin-size', type=int, default=200, help='IGNORE: Bin size in kb (default: 200)')
    
    args = parser.parse_args()
    
    if args.command == 'pairwise-comparison':

        try:
            config_overrides = {
                'position_bin_size_kb': args.bin_size,
                'max_iterations': args.iterations,
                'use_genomic_positions': args.mode == 'position',
                'flexible_threshold': args.flexible,
                'shared_genes_threshold': args.shared_genes
            }
            results = pairwise_comparison_command(
                args.genome1_busco, 
                args.genome2_busco, 
                args.output, 
                config_overrides
            )
           
           
            print_detailed_metrics(results, results.get('timing', {}).get('start_time'), results.get('timing', {}).get('end_time'))
           
            # Extract genome names from file paths
            genome1_name = Path(args.genome1_busco).stem
            genome2_name = Path(args.genome2_busco).stem

            
            # Generate detailed reports
            try:
                generate_iteration_report(results, args.output, genome1_name, genome2_name)
            except Exception as e:
                print(f"Error generating iteration report: {e}")
            
            try:
                save_pattern_analysis(results, args.output, genome1_name, genome2_name)
                print(f"Pattern analysis saved successfully")
            except Exception as e:
                print(f"Error saving pattern analysis: {e}")
           
           
            
            print(f"\nPairwise comparison completed!")
            print(f"Distance: {results['distance_metrics']['total_events']} events")
            print(f"Gene inversions: {results['distance_metrics']['total_gene_inversions']}")
            print(f"Iterations: {results['distance_metrics']['iterations']}")
            return 0
        except Exception as e:
            print(f"Pairwise comparison failed: {e}")
            return 1
    else:
        parser.print_help()
        return 1




if __name__ == "__main__":
   sys.exit(main())
