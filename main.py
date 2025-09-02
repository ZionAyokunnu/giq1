#!/usr/bin/env python3
"""
GIQ2 - Genome Inversion Quantifier
Main entry point.

Author: Zion Ayokunnu
Supervisors: Kamil, Sasha, Arif, Sam
Version: 1.0.0

1. align-chr: Align chromosomes across genomes using RagTag
2. build-profile: Build Markov profile from multiple training genomes
3. analyze-query: Analyze query genome against existing profile
"""

import sys
import argparse
import logging
import os
from pathlib import Path
import pandas as pd
import json
from typing import Dict, List

from config.settings import CONFIG
from utils.file_utils import create_output_directory
from core import (
    parse_busco_table,
    filter_busco_genes,
    get_genome_name_from_fasta,
    create_pairwise_movement_sequence_per_chromosome,
    iterative_detection,
    build_profile_command_hybrid,
    analyze_query_command_hybrid,
    create_convergence_analysis_tsv,
    save_stage_data,
    multi_reference_alignment_pipeline
    
)

from contextual.metrics import (
    compute_inversion_rate_per_mb_busco,
    _analyse_gene_density_correlation
)

from visualisation.plots import (
    create_linearised_dotplot,
    create_busco_dotplot
)

from core.convergence_analysis import create_single_convergence_tsv

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)




def pairwise_comparison_command(busco_file1: str, busco_file2: str, output_dir: str, config_overrides: Dict = None):
    """Direct pairwise genome comparison using iterative movement analysis."""
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
    
    # DEBUG: Print filtered data immediately after filtering
    print(f"DEBUG: Genome1 filtered - {len(genome1_filtered)} genes")
    print(f"DEBUG: Genome1 chromosomes: {genome1_filtered['sequence'].unique()}")
    print(f"DEBUG: Genome2 filtered - {len(genome2_filtered)} genes")
    print(f"DEBUG: Genome2 chromosomes: {genome2_filtered['sequence'].unique()}")
    
    # DEBUG: Check specific chromosome data
    chr1_740 = genome1_filtered[genome1_filtered['sequence'] == 'OZ005740.1']
    chr2_740 = genome2_filtered[genome2_filtered['sequence'] == 'OZ002740.1']
    print(f"DEBUG: OZ005740.1 genes: {len(chr1_740)}")
    print(f"DEBUG: OZ002740.1 genes: {len(chr2_740)}")
    print(f"DEBUG: First 5 genes in OZ005740.1:")
    print(chr1_740[['busco_id', 'gene_start']].head())
    
    save_stage_data(genome1_filtered, f'1_parsed_{genome1_name}', output_path, f"Parsed genome 1: {genome1_name}")
    save_stage_data(genome2_filtered, f'1_parsed_{genome2_name}', output_path, f"Parsed genome 2: {genome2_name}")
    
    # Create movement sequence directly from position differences
    chromosome_sequences = create_pairwise_movement_sequence_per_chromosome(genome1_filtered, genome2_filtered, config)
    
    # save_stage_data(pd.DataFrame(chromosome_sequences, columns=['gene_id', 'rank', 'movement']), 
    #                '2_movement_sequence', output_path, f"Movement sequence: {len(chromosome_sequences)} genes")

    
    
    # Note: We only do per-chromosome analysis, not combined analysis
    
    
        # Run iterative detection on each chromosome pair
    all_results = {}
    total_events = 0
    total_gene_inversions = 0
    all_converged_data = []  # Collect all converged data
    
    for chromosome_pair, movement_sequence in chromosome_sequences.items():
        print(f"Processing {chromosome_pair}...")
        
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
    if all_converged_data:
        save_combined_converged_genome(all_converged_data, output_path, genome1_name, genome2_name)
        
        # Create detailed convergence analysis
        converged_df = pd.DataFrame(all_converged_data)
        create_single_convergence_tsv(genome1_df, converged_df, chromosome_sequences, all_results, output_path, genome1_name, genome2_name)
    
    results = {
        'genome1': genome1_name,
        'genome2': genome2_name,
        'distance_metrics': distance_metrics,
        'chromosome_results': all_results
    }
    
    results_file = output_path / f'pairwise_results_{genome1_name}_vs_{genome2_name}.json'
    with open(results_file, 'w') as f:
        json.dump(results, f, indent=2, default=str)
    
    logger.info(f"Results saved to: {results_file}")
    return results



def save_inversion_events(events: List[Dict], output_path: Path, comparison_name: str):
    """Save inversion events to CSV with iteration details."""
    if not events:
        logger.info("No inversion events to save")
        return
    
    events_records = []
    for event in events:
        events_records.append({
            'iteration': event['iteration'],
            'type': event['type'],
            'genes_involved': len(event['genes']),
            'gene_list': ';'.join(event['genes']),
            'positions': ';'.join(map(str, event['positions'])),
            'gene_inversions': event['gene_inversions']
        })
    
    events_df = pd.DataFrame(events_records)
    save_stage_data(events_df, f'4_inversion_events_{comparison_name}', output_path, 
                   f"Inversion events: {len(events)} events")


def save_converged_genome(final_sequence, genome1_df, genome2_df, chromosome_pair, output_path):
    """
    Save the final converged genome sequence in the same format as the original genome.
    
    Args:
        final_sequence: List of (gene_id, rank, movement) tuples from iterative detection
        genome1_df: Original genome 1 DataFrame
        genome2_df: Original genome 2 DataFrame  
        chromosome_pair: Name of the chromosome pair being processed
        output_path: Output directory path
    """
    # Extract chromosome names from pair name
    chr1_name, chr2_name = chromosome_pair.split('_vs_')
    
    # Get the target genome (genome2) data for this chromosome
    target_chr_data = genome2_df[genome2_df['sequence'] == chr2_name].copy()
    
    # Create a mapping from gene_id to final rank
    gene_to_final_rank = {}
    for gene_id, final_rank, _, _ in final_sequence:
        gene_to_final_rank[gene_id] = final_rank
    
    # Update the target chromosome data with final ranks
    converged_data = []
    
    for _, row in target_chr_data.iterrows():
        gene_id = row['busco_id']
        if gene_id in gene_to_final_rank:
            # Create new row with updated gene_start based on final rank
            new_row = row.copy()
            final_rank = gene_to_final_rank[gene_id]
            
            # Use the final rank as the gene_start (linearized position)
            # Multiply by a size factor to create proper genomic coordinates
            size_factor = 288000  # 1Mb spacing between genes
            new_row['gene_start'] = final_rank * size_factor
            # Set gene_end to be gene_start + size_factor to maintain spacing
            new_row['gene_end'] = (final_rank + 1) * size_factor
            
            converged_data.append(new_row)
    
    if converged_data:
        converged_df = pd.DataFrame(converged_data)
        
        # Sort by final gene_start to maintain order
        converged_df = converged_df.sort_values('gene_start')
        
        # Save as TSV in the same format as original
        output_file = output_path / f'5_converged_genome_{chromosome_pair}.tsv'
        converged_df.to_csv(output_file, sep='\t', index=False)
        
        logger.info(f"Converged genome saved: {output_file} ({len(converged_df)} genes)")
        
        # Also save as CSV for easier plotting
        csv_file = output_path / f'5_converged_genome_{chromosome_pair}.csv'
        converged_df.to_csv(csv_file, index=False)
        
        logger.info(f"Converged genome CSV saved: {csv_file}")
    else:
        logger.warning(f"No converged data for {chromosome_pair}")


def get_converged_genome_data(final_sequence, genome1_df, genome2_df, chromosome_pair):
    """
    Get converged genome data for a single chromosome without saving.
    
    Args:
        final_sequence: List of (gene_id, rank, movement) tuples from iterative detection
        genome1_df: Original genome 1 DataFrame
        genome2_df: Original genome 2 DataFrame  
        chromosome_pair: Name of the chromosome pair being processed
        
    Returns:
        List of converged data rows
    """
    # Extract chromosome names from pair name
    chr1_name, chr2_name = chromosome_pair.split('_vs_')
    
    # Get the target genome (genome2) data for this chromosome - we're converging rufipes
    target_chr_data = genome2_df[genome2_df['sequence'] == chr2_name].copy()
    
    # Create a mapping from gene_id to final rank
    gene_to_final_rank = {}
    for gene_id, final_rank, _, _ in final_sequence:
        gene_to_final_rank[gene_id] = final_rank
    
    # Update the target chromosome data with final ranks
    converged_data = []
    
    for _, row in target_chr_data.iterrows():
        gene_id = row['busco_id']
        if gene_id in gene_to_final_rank:
            # Create new row with updated gene_start based on final rank
            new_row = row.copy()
            final_rank = gene_to_final_rank[gene_id]
            
            # Use the final rank as the gene_start (linearized position)
            # Multiply by a size factor to create proper genomic coordinates
            size_factor = 288000  # 1Mb spacing between genes
            new_row['gene_start'] = final_rank * size_factor
            # Set gene_end to be gene_start + size_factor to maintain spacing
            new_row['gene_end'] = (final_rank + 1) * size_factor
            
            converged_data.append(new_row)
    
    return converged_data


def save_combined_converged_genome(all_converged_data, output_path, genome1_name, genome2_name):
    """
    Save all converged genome data together in a single TSV file.
    
    Args:
        all_converged_data: List of all converged data rows from all chromosomes
        output_path: Output directory path
        genome1_name: Name of genome 1
        genome2_name: Name of genome 2
    """
    if not all_converged_data:
        logger.warning("No converged data to save")
        return
    
    # Create DataFrame from all converged data
    converged_df = pd.DataFrame(all_converged_data)
    
    # Sort by sequence (chromosome) then by gene_start to maintain order
    converged_df = converged_df.sort_values(['sequence', 'gene_start'])
    
    # Save as TSV in the same format as original
    output_file = output_path / f'5_converged_genome_{genome1_name}_vs_{genome2_name}.tsv'
    converged_df.to_csv(output_file, sep='\t', index=False)
    
    logger.info(f"Combined converged genome saved: {output_file} ({len(converged_df)} genes from {len(converged_df['sequence'].unique())} chromosomes)")
    
    # Show chromosome distribution
    chr_counts = converged_df['sequence'].value_counts().sort_index()
    for chr_name, count in chr_counts.items():
        logger.info(f"  {chr_name}: {count} genes")


















#Giq1 stuff


# def save_results(inversion_df: pd.DataFrame,
#                 joined_df: pd.DataFrame,
#                 contextual_metrics: Dict,
#                 output_dir: Path,
#                 config: Dict):
#     """Save all analysis results to files"""
  
#     data_dir = output_dir / 'data'
#     data_dir.mkdir(exist_ok=True)
  
#     inversion_df.to_csv(data_dir / 'inversion_analysis.csv', index=False)
#     joined_df.to_csv(data_dir / 'joined_gene_data.csv', index=False)
  
#     import json
#     with open(data_dir / 'contextual_metrics.json', 'w') as f:
#         json.dump(contextual_metrics, f, indent=2, default=str)
      
      
# def run_busco_inversion_analysis(config: Dict = None) -> Dict:
#     """   
#     Args:
#         config: Configuration dictionary (uses default CONFIG if None)
      
#     """
  
#     if config is None:
#         config = CONFIG
  
#     output_dir = create_output_directory(config)
#     logger.info(f"Output directory: {output_dir}")
  
#     species1_name = config.get('first_species_name', 'Species1')
#     species2_name = config.get('second_species_name', 'Species2')
  
#     logger.info("\n" + "-" * 40)
#     logger.info(f"Analyzing inversions between {species1_name} and {species2_name}")
#     logger.info("-" * 40)
  
#     try:


#         busco_df1 = parse_busco_table(config['first_busco_path'], config)
#         busco_df2 = parse_busco_table(config['second_busco_path'], config)
  
      
#         filtered_df1 = filter_busco_genes(busco_df1, config)
#         filtered_df2 = filter_busco_genes(busco_df2, config)
      
#         logger.info("\n" + "-" * 40)
 #         inversion_df, joined_df = detect_flips(filtered_df1, filtered_df2, config)
      
#         logger.info(f"Detected inversions across {len(inversion_df)} chromosome pairs")
#         logger.info(f"Total genes compared: {len(joined_df)}")
#         logger.info(f"Total flipped genes: {joined_df['is_flipped'].sum()}")
      
      
#         contextual_metrics = compute_contextual_metrics(
#             inversion_df, joined_df, config
#         )
      
#         logger.info("\n" + "-" * 40)


      
#         visualisation_results = create_analysis_visualisations(
#             inversion_df, joined_df, output_dir, config
#         )
      
#         save_results(
#             inversion_df, joined_df, contextual_metrics, output_dir, config
#         )
      
#         report_path = generate_analysis_report(
#             inversion_df, joined_df, contextual_metrics,
#             output_dir, config
#         )
      
#         results = {
#             'inversion_df': inversion_df,
#             'joined_df': joined_df,
#             'contextual_metrics': contextual_metrics,
#             'visualisation_results': visualisation_results,
#             'output_dir': output_dir,
#             'report_path': report_path,
#             'config': config,
#             'species_names': [species1_name, species2_name]
#         }
      
#         logger.info(f"Report available at: {report_path}")
#         logger.info("=" * 80)
      
#         return results
      
#     except Exception as e:
#         logger.error(f"Analysis failed: {str(e)}")
#         if config.get('enable_debug_output', False):
#             import traceback
#             traceback.print_exc()
#         raise      
      
      
# def compute_contextual_metrics(inversion_df: pd.DataFrame,
#                              joined_df: pd.DataFrame,
#                              config: Dict) -> Dict:
  
#     contextual_metrics = {}
  
#     rate_metrics = compute_inversion_rate_per_mb_busco(inversion_df, joined_df)
#     contextual_metrics['inversion_rates'] = rate_metrics
  
#     gene_density_corr = _analyse_gene_density_correlation(inversion_df, joined_df, config)
#     contextual_metrics['gene_density_correlation'] = gene_density_corr


  
#     return contextual_metrics




# def create_analysis_visualisations(inversion_df: pd.DataFrame,
#                                  joined_df: pd.DataFrame,
#                                  output_dir: Path,
#                                  config: Dict) -> Dict:
  
#     plots_dir = output_dir / 'plots'
#     plots_dir.mkdir(exist_ok=True)
  
#     visualisation_results = {}
  
#     species1_name = config.get('first_species_name', 'Species1')
#     species2_name = config.get('second_species_name', 'Species2')
  
#     import matplotlib.pyplot as plt
#     plt.style.use('default')
#     plt.rcParams['figure.dpi'] = config.get('dpi', 300)
#     plt.rcParams['font.size'] = config.get('font_size', 12)
  
#     try:
#         try:
#             create_busco_dotplot(joined_df, plots_dir, config)
#             visualisation_results['dotplot'] = plots_dir / 'busco_dotplot.png'
#         except Exception as e:
#             logger.warning(f"create_busco_dotplot failed: {e}")
          
#         try:
#             fig, ax = create_linearised_dotplot(
#                 joined_df=joined_df,
#                 plots_dir=plots_dir,
#                 config=config['dotplot_config']  # Pass the nested dotplot config
#             )
#             visualisation_results['linearized_dotplot'] = plots_dir / 'linearized_busco_dotplot.png'
          
#         except Exception as e:
#             logger.warning(f"Linearized dotplot failed: {e}")
          
#     except Exception as e:
#         logger.warning(f"No Vis worked! Ouch!")








# def generate_analysis_report(inversion_df: pd.DataFrame,
#                            joined_df: pd.DataFrame,
#                            contextual_metrics: Dict,
#                            output_dir: Path,
#                            config: Dict) -> Path:
  
#     reports_dir = output_dir / 'reports'
#     reports_dir.mkdir(exist_ok=True)
  
#     report_path = reports_dir / 'analysis_report.md'
  
#     species1_name = config.get('first_species_name', 'Species1')
#     species2_name = config.get('second_species_name', 'Species2')
  
#     with open(report_path, 'w') as f:
#         f.write("# GIQ1 - Genome Inversion Analysis Report\n\n")
#         f.write(f"**Generated:** {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
#         f.write(f"**Species Comparison:** {species1_name} vs {species2_name}\n\n")
      
#         f.write("## Analysis Summary\n\n")
#         f.write(f"- **Total chromosome pairs analyzed:** {len(inversion_df)}\n")
#         f.write(f"- **Pairs with inversions:** {len(inversion_df[inversion_df['flipped_genes'] > 0])}\n")
#         f.write(f"- **Total genes compared:** {len(joined_df)}\n")
#         f.write(f"- **Total flipped genes:** {joined_df['is_flipped'].sum()}\n")
#         f.write(f"- **Overall flip rate:** {joined_df['is_flipped'].mean():.1%}\n\n")
      
#         f.write("## Inversion Types Distribution\n\n")
#         inversion_types = inversion_df['inversion_type'].value_counts()
#         for inv_type, count in inversion_types.items():
#             f.write(f"- **{inv_type}:** {count} pairs\n")
#         f.write("\n")
      
#         f.write("## Contextual Metrics\n\n")
#         if 'inversion_rates' in contextual_metrics:
#             rates = contextual_metrics['inversion_rates']
#             f.write(f"- **Genome coverage (Mb):** {rates.get('genome_coverage_mb', 0):.2f}\n")
#             f.write(f"- **Inversions per Mb:** {rates.get('inversion_events_per_mb', 0):.3f}\n")
#             f.write(f"- **Flipped genes per Mb:** {rates.get('flipped_genes_per_mb', 0):.3f}\n")
#         f.write("\n")
      
#         f.write("## Chromosome Pair Details\n\n")
#         f.write("| Chr1 | Chr2 | Total Genes | Flipped | Rate | Correlation | Type |\n")
#         f.write("|------|------|-------------|---------|------|-------------|------|\n")
      
#         for _, row in inversion_df.iterrows():
#             f.write(f"| {row['chr1']} | {row['chr2']} | {row['total_genes']} | "
#                    f"{row['flipped_genes']} | {row['flip_rate']:.2f} | "
#                    f"{row.get('correlation', 'N/A')} | {row['inversion_type']} |\n")
      
#         f.write(f"\n---\n*Report generated by GIQ1 v1.0.0*\n")
  
#     logger.info(f"Analysis report saved to {report_path}")
#     return report_path






















































def main():
   """Main CLI entry point"""
   parser = argparse.ArgumentParser(description='Genome Inversion Quantification Tools')
   subparsers = parser.add_subparsers(dest='command', help='Available commands')
  
   # Align chromosomes command
   align_parser = subparsers.add_parser('align-chr', help='Align chromosomes across genomes using RagTag')
   align_parser.add_argument('genome_fastas', nargs='+', help='Query genome FASTA files')
   align_parser.add_argument('-o', '--output', required=True, help='Output directory for alignment results')
   align_parser.add_argument('--busco-files', nargs='+', help='BUSCO TSV files (same order as FASTA files)') 
   
   # Build profile command
   profile_parser = subparsers.add_parser('build-profile', help='Build Markov profile from training genomes')
   profile_parser.add_argument('busco_files', nargs='+', help='BUSCO table files for training genomes')
   profile_parser.add_argument('-o', '--output', required=True, help='Output directory for profile')
   profile_parser.add_argument('--bin-size', type=int, default=2, help='Bin size in kb (default: 2)')
   profile_parser.add_argument('--method', choices=['average', 'range'], default='range', help='Profile calculation method')
   profile_parser.add_argument('--chr-map', help='Path to chromosome mappings JSON file from align-chr command')

  
   # Analyse query command
   query_parser = subparsers.add_parser('analyze-query', help='Analyze query genome against profile')
   query_parser.add_argument('query_busco', help='BUSCO table file for query genome')
   query_parser.add_argument('profile', help='Saved Markov profile JSON file')
   query_parser.add_argument('-o', '--output', required=True, help='Output directory for analysis')
   query_parser.add_argument('--threshold', type=float, default=0.7, help='Probability threshold (default: 0.7)')
   query_parser.add_argument('--bin-size', type=int, help='Bin size in kb (overrides profile default)') 
    
  
   # Pairwise comparison command
   pairwise_parser = subparsers.add_parser('pairwise-comparison', help='Direct pairwise genome comparison')
   pairwise_parser.add_argument('genome1_busco', help='First BUSCO table file')
   pairwise_parser.add_argument('genome2_busco', help='Second BUSCO table file') 
   pairwise_parser.add_argument('-o', '--output', required=True, help='Output directory')
   pairwise_parser.add_argument('--iterations', type=int, default=1000, help='Max iterations (default: 1000)')
   pairwise_parser.add_argument('--bin-size', type=int, default=200, help='Bin size in kb (default: 200)')
   pairwise_parser.add_argument('--mode', choices=['rank', 'position'], default='rank', help='Use gene ranks or genomic positions (default: rank)')
   pairwise_parser.add_argument('--flexible', type=float, default=0.0, help='Allow movement increases up to this value (default: 0.0 = strict)')
   pairwise_parser.add_argument('--shared-genes', type=int, default=50, help='Minimum shared genes required for chromosome pairing (default: 50)')
  
   args = parser.parse_args()
  
   if args.command == 'align-chr':
        try:
           
            busco_mapping = None
            if args.busco_files:
                if len(args.busco_files) != len(args.genome_fastas):
                    print("Error: Number of BUSCO files must match number of FASTA files")
                    return 1
                
                busco_mapping = {}
                for fasta_path, busco_path in zip(args.genome_fastas, args.busco_files):
                    genome_name = get_genome_name_from_fasta(fasta_path)
                    busco_mapping[genome_name] = busco_path
            
            mappings_file = multi_reference_alignment_pipeline(
                args.genome_fastas,
                args.output,
                busco_mapping
            )


            print(f"\nChromosome alignment completed successfully!")
            print(f"Mappings file: {mappings_file}")
            print(f"Use this mappings file with build-profile: --chr-map {mappings_file}")
            return 0
        except Exception as e:
           print(f"Chromosome alignment failed: {e}")
           return 1
  
   elif args.command == 'build-profile':
       config_overrides = {
           'position_bin_size_kb': args.bin_size,
           'profile_calculation_method': args.method
       }
      
       try:
           profile_data = build_profile_command_hybrid(
               args.busco_files,
               args.output,
               config_overrides,
               args.chr_map
           )
           print(f"\nProfile building completed successfully!")
           print(f"Profile saved to: {args.output}/markov_profile.json")
           if args.chr_map:
               print(f"Chromosome standardization: Applied")
           else:
               print(f"Chromosome standardization: Not applied")
           return 0
       except Exception as e:
           print(f"Profile building failed: {e}")
           return 1
  
   elif args.command == 'analyze-query':
       config_overrides = {
           'permutable_positions_threshold': args.threshold
       }
       
       if args.bin_size is not None:
            config_overrides['position_bin_size_kb'] = args.bin_size
        
       try:
           results = analyze_query_command_hybrid(args.query_busco, args.profile, args.output, config_overrides)
           print(f"\nQuery analysis completed successfully!")
           print(f"Results saved to: {args.output}")
           return 0
       except Exception as e:
           print(f"Query analysis failed: {e}")
           return 1
  
        
   elif args.command == 'pairwise-comparison':
       try:
           config_overrides = {
               'position_bin_size_kb': args.bin_size,
               'max_iterations': args.iterations,
               'use_genomic_positions': args.mode == 'position',
               'flexible_threshold': args.flexible,
               'shared_genes_threshold': args.shared_genes
           }
           print(f"DEBUG: args.mode = '{args.mode}', use_genomic_positions = {config_overrides['use_genomic_positions']}")
           results = pairwise_comparison_command(
               args.genome1_busco, 
               args.genome2_busco, 
               args.output, 
               config_overrides
           )
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
