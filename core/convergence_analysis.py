"""
Convergence Analysis Module
Handles detailed analysis of genome convergence results
"""

import pandas as pd
from pathlib import Path
import logging

logger = logging.getLogger(__name__)


def create_single_convergence_tsv(genome1_df, converged_df, movement_sequences, converged_sequences, output_path, genome1_name, genome2_name):
    """
    Create a single TSV with the exact format requested:
    gene_id, {genome1_name}_positions, {genome2_name}_convergence_position, convergence_status, chr_concatenated, current_movement_value
    
    Args:
        genome1_df: DataFrame of original genome 1 (linearis)
        converged_df: DataFrame of converged genome 2 (rufipes)
        movement_sequences: Dictionary of movement sequences per chromosome
        output_path: Output directory path
        genome1_name: Name of genome 1 (linearis)
        genome2_name: Name of genome 2 (rufipes)
    """
    # Get shared genes between both genomes
    genome1_genes = set(genome1_df['busco_id'])
    converged_genes = set(converged_df['busco_id'])
    shared_genes = genome1_genes & converged_genes
    
    logger.info(f"Creating single convergence TSV with {len(shared_genes)} shared genes")
    
    # Create mappings for quick lookup
    converged_gene_to_chr = {row['busco_id']: row['sequence'] for _, row in converged_df.iterrows()}
    genome1_gene_to_chr = {row['busco_id']: row['sequence'] for _, row in genome1_df.iterrows()}
    
    # Create rank and movement value mappings from movement sequences (initial positions)
    genome1_ranks = {}
    converged_ranks = {}
    movement_values = {}
    
    # First, get initial positions from movement sequences
    for chr_pair, sequence in movement_sequences.items():
        for gene_id, rank, movement, target_pos in sequence:
            genome1_ranks[gene_id] = rank
    
    # Get final converged ranks and movement values from converged sequences (algorithm's final state)
    for chr_pair, result in converged_sequences.items():
        if isinstance(result, dict):
            # Extract the final sequence from the result dictionary
            if 'final_sequence' in result:
                final_sequence = result['final_sequence']
                if isinstance(final_sequence, list):
                    for gene_id, rank, movement, target_pos in final_sequence:
                        converged_ranks[gene_id] = target_pos  # Use target_pos, not rank!
                        # Use the algorithm's final movement values (this is the key fix!)
                        movement_values[gene_id] = movement
        elif isinstance(result, tuple) and len(result) >= 2:
            # Extract the final sequence from the result tuple
            final_sequence = result[1]  # The converged sequence
            if isinstance(final_sequence, list):
                for gene_id, rank, movement, target_pos in final_sequence:
                    converged_ranks[gene_id] = target_pos  # Use target_pos, not rank!
                    # Use the algorithm's final movement values (this is the key fix!)
                    movement_values[gene_id] = movement
    
    # Prepare data for the TSV
    analysis_data = []
    
    for gene_id in sorted(shared_genes):
        # Get ranks and chromosomes
        genome1_rank = genome1_ranks.get(gene_id, 'N/A')
        genome1_chr = genome1_gene_to_chr.get(gene_id, 'N/A')
        converged_rank = converged_ranks.get(gene_id, 'N/A')
        converged_chr = converged_gene_to_chr.get(gene_id, 'N/A')
        
        # Use the algorithm's final movement values (from iterative_detection)
        current_movement = movement_values.get(gene_id, 0)
        
        # DEBUG: Add detailed tracing for problematic genes
        if gene_id in ['5858at7147', '1072at7147', '4511at7147', '5845at7147', '4606at7147', '4164at7147', '5263at7147']:
            print(f"DEBUG {gene_id}:")
            print(f"  genome1_rank: {genome1_rank}")
            print(f"  converged_rank: {converged_rank}")
            print(f"  algorithm_final_movement: {current_movement}")
            print(f"  actual_movement_calculation: {converged_rank - genome1_rank if isinstance(converged_rank, (int, float)) and isinstance(genome1_rank, (int, float)) else 'N/A'}")
            print(f"  movement_values dict contains: {gene_id in movement_values}")
            if gene_id in movement_values:
                print(f"  movement_values[{gene_id}]: {movement_values[gene_id]}")
        
        # Determine convergence status
        if genome1_chr == converged_chr:
            if isinstance(genome1_rank, (int, float)) and isinstance(converged_rank, (int, float)):
                pos_diff = abs(converged_rank - genome1_rank)
                if pos_diff <= 1:  # Within 1 rank
                    convergence_status = "CONVERGED"
                else:
                    convergence_status = "POSITION_DIFFERENT"
            else:
                convergence_status = "POSITION_UNKNOWN"
        else:
            convergence_status = "CHROMOSOME_DIFFERENT"
        
        # Create concatenated chromosome string
        chr_concatenated = f"{genome1_chr}-{converged_chr}"
        
        analysis_data.append({
            'gene_id': gene_id,
            f'{genome1_name}_positions': genome1_rank,
            f'{genome2_name}_convergence_position': converged_rank,
            'convergence_status': convergence_status,
            'chr_concatenated': chr_concatenated,
            'current_movement_value': current_movement
        })
    
    # Create DataFrame
    analysis_df = pd.DataFrame(analysis_data)
    
    # Sort by chromosome, then by converged genome positions (final order)
    analysis_df = analysis_df.sort_values(['chr_concatenated', f'{genome2_name}_convergence_position'])
    
    # Save the single TSV
    output_file = output_path / f'convergence_analysis_{genome1_name}_vs_{genome2_name}.tsv'
    analysis_df.to_csv(output_file, sep='\t', index=False)
    
    logger.info(f"Single convergence TSV saved: {output_file}")
    logger.info(f"Total shared genes analyzed: {len(analysis_df)}")
    
    # Print summary
    status_counts = analysis_df['convergence_status'].value_counts()
    logger.info("Convergence summary:")
    for status, count in status_counts.items():
        percentage = (count / len(analysis_df)) * 100
        logger.info(f"  {status}: {count} genes ({percentage:.1f}%)")
    
    return analysis_df


def create_convergence_analysis_tsv(converged_df, linearis_df, output_path, genome1_name, genome2_name):
    """
    Create detailed convergence analysis TSV comparing converged rufipes vs linearis.
    
    Args:
        converged_df: DataFrame of converged rufipes genome
        linearis_df: DataFrame of original linearis genome
        output_path: Output directory path
        genome1_name: Name of genome 1 (linearis)
        genome2_name: Name of genome 2 (rufipes)
    """
    # Create a comprehensive comparison
    analysis_data = []
    
    # Get all unique genes from both genomes
    all_genes = set(converged_df['busco_id']) | set(linearis_df['busco_id'])
    
    # Create mappings for quick lookup
    converged_gene_to_pos = {row['busco_id']: row['gene_start'] for _, row in converged_df.iterrows()}
    converged_gene_to_chr = {row['busco_id']: row['sequence'] for _, row in converged_df.iterrows()}
    linearis_gene_to_pos = {row['busco_id']: row['gene_start'] for _, row in linearis_df.iterrows()}
    linearis_gene_to_chr = {row['busco_id']: row['sequence'] for _, row in linearis_df.iterrows()}
    
    for gene_id in sorted(all_genes):
        # Check if gene exists in both genomes
        in_converged = gene_id in converged_gene_to_pos
        in_linearis = gene_id in linearis_gene_to_pos
        
        # Get positions and chromosomes
        converged_pos = converged_gene_to_pos.get(gene_id, 'N/A')
        converged_chr = converged_gene_to_chr.get(gene_id, 'N/A')
        linearis_pos = linearis_gene_to_pos.get(gene_id, 'N/A')
        linearis_chr = linearis_gene_to_chr.get(gene_id, 'N/A')
        
        # Determine convergence status
        if in_converged and in_linearis:
            # Both genomes have this gene
            if converged_chr == linearis_chr:
                # Same chromosome - check if positions are similar
                if isinstance(converged_pos, (int, float)) and isinstance(linearis_pos, (int, float)):
                    pos_diff = abs(converged_pos - linearis_pos)
                    if pos_diff <= 1000:  # Within 1kb
                        convergence_status = "CONVERGED"
                    else:
                        convergence_status = "POSITION_DIFFERENT"
                else:
                    convergence_status = "POSITION_UNKNOWN"
            else:
                convergence_status = "CHROMOSOME_DIFFERENT"
        elif in_converged:
            convergence_status = "ONLY_IN_CONVERGED"
        elif in_linearis:
            convergence_status = "ONLY_IN_LINEARIS"
        else:
            convergence_status = "UNKNOWN"
        
        # Calculate position difference for analysis
        pos_diff = None
        if isinstance(converged_pos, (int, float)) and isinstance(linearis_pos, (int, float)):
            pos_diff = converged_pos - linearis_pos
        
        analysis_data.append({
            'gene_id': gene_id,
            'converged_rufipes_position': converged_pos,
            'converged_rufipes_chromosome': converged_chr,
            'linearis_position': linearis_pos,
            'linearis_chromosome': linearis_chr,
            'position_difference': pos_diff,
            'convergence_status': convergence_status,
            'in_converged': in_converged,
            'in_linearis': in_linearis
        })
    
    # Create DataFrame and save
    analysis_df = pd.DataFrame(analysis_data)
    
    # Sort by convergence status and position difference
    analysis_df = analysis_df.sort_values(['convergence_status', 'position_difference'])
    
    # Save detailed analysis
    output_file = output_path / f'6_convergence_analysis_{genome1_name}_vs_{genome2_name}.tsv'
    analysis_df.to_csv(output_file, sep='\t', index=False)
    
    # Create summary
    status_counts = analysis_df['convergence_status'].value_counts()
    total_genes = len(analysis_df)
    
    logger.info(f"Convergence analysis saved: {output_file}")
    logger.info(f"Total genes analyzed: {total_genes}")
    logger.info("Convergence summary:")
    for status, count in status_counts.items():
        percentage = (count / total_genes) * 100
        logger.info(f"  {status}: {count} genes ({percentage:.1f}%)")
    
    # Save non-converged genes separately for transposition analysis
    non_converged = analysis_df[analysis_df['convergence_status'].isin(['POSITION_DIFFERENT', 'CHROMOSOME_DIFFERENT'])]
    if len(non_converged) > 0:
        non_converged_file = output_path / f'7_non_converged_genes_{genome1_name}_vs_{genome2_name}.tsv'
        non_converged.to_csv(non_converged_file, sep='\t', index=False)
        logger.info(f"Non-converged genes saved: {non_converged_file} ({len(non_converged)} genes)")
    
    return analysis_df


def analyze_transposition_patterns(non_converged_df):
    """
    Analyze patterns in non-converged genes to identify potential transposition candidates.
    
    Args:
        non_converged_df: DataFrame of non-converged genes
        
    Returns:
        dict: Analysis results with transposition patterns
    """
    patterns = {
        'chromosome_transpositions': {},
        'position_clusters': {},
        'gene_groups': []
    }
    
    # Analyze chromosome differences
    chr_diff = non_converged_df[non_converged_df['convergence_status'] == 'CHROMOSOME_DIFFERENT']
    if len(chr_diff) > 0:
        chr_mapping = chr_diff.groupby(['linearis_chromosome', 'converged_rufipes_chromosome']).size()
        patterns['chromosome_transpositions'] = chr_mapping.to_dict()
    
    # Analyze position differences
    pos_diff = non_converged_df[non_converged_df['convergence_status'] == 'POSITION_DIFFERENT']
    if len(pos_diff) > 0:
        # Group by chromosome and analyze position differences
        for chr_name in pos_diff['linearis_chromosome'].unique():
            chr_data = pos_diff[pos_diff['linearis_chromosome'] == chr_name]
            patterns['position_clusters'][chr_name] = {
                'total_genes': len(chr_data),
                'avg_position_diff': chr_data['position_difference'].mean(),
                'max_position_diff': chr_data['position_difference'].max(),
                'min_position_diff': chr_data['position_difference'].min()
            }
    
    return patterns
