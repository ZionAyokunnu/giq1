import numpy as np
import logging
import time
from datetime import datetime
import pandas as pd
import json
from pathlib import Path
from typing import Dict, List


logger = logging.getLogger(__name__)       
            
            
            
def generate_iteration_report(results, output_path, genome1_name, genome2_name):
    """Generate detailed Markdown report of all iterations and flip patterns"""
    
    # Ensure output_path is a Path object
    if isinstance(output_path, str):
        output_path = Path(output_path)
    
    report_lines = []
    
    # Header
    report_lines.extend([
        f"# Genomic Rearrangement Analysis Report",
        f"**Genomes:** {genome1_name} → {genome2_name}",
        f"**Analysis Date:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
        "",
        "---",
        ""
    ])
    
    # Summary - ensure they are numbers
    total_events = int(results['distance_metrics']['total_events'])
    total_genes = int(results['distance_metrics']['total_gene_inversions'])
    iterations = int(results['distance_metrics']['iterations'])
    
    report_lines.extend([
        "## Summary",
        f"- **Total Events:** {total_events}",
        f"- **Genes Rearranged:** {total_genes}",
        f"- **Iterations:** {iterations}",
        f"- **Chromosomes Processed:** {len(results['distance_metrics']['chromosome_results'])}",
        "",
        "---",
        ""
    ])
    
    # Per-chromosome analysis
    for chr_name, chr_data in results['distance_metrics']['chromosome_results'].items():
        report_lines.extend([
            f"## {chr_name}",
            f"- **Genes:** {len(chr_data.get('final_sequence', []))}",
            f"- **Events:** {chr_data.get('total_events', 0)}",
            f"- **Iterations:** {chr_data.get('iterations', 0)}",
            f"- **Convergence:** {'✅ Complete' if chr_data.get('converged', False) else '⚠️ Partial'}",
            ""
        ])
        
        # Iteration details
        if 'inversion_events' in chr_data:
            current_iteration = 0
            
            for event_idx, event in enumerate(chr_data['inversion_events']):
                iteration = event.get('iteration', 1)
                
                # New iteration header
                if iteration != current_iteration:
                    current_iteration = iteration
                    report_lines.extend([
                        f"### Iteration {iteration}",
                        ""
                    ])
                
                # Event details
                event_type = event.get('type', 'unknown')
                genes = event.get('genes', [])
                positions = event.get('positions', [])
                
                report_lines.extend([
                    f"#### Event {event_idx + 1}: {event_type.upper()} Inversion",
                    ""
                ])
                
                # Pattern analysis
                if event_type == 'flip':
                    flip_size = event.get('flip_indicator', 1)
                    segment_length = len(genes)
                    
                    # Determine pattern type
                    if segment_length == 2:
                        pattern_type = "Simple Adjacent Swap"
                    elif segment_length % 2 == 1:
                        pattern_type = f"Odd-length Incremental (Natural Fulcrum, Size {flip_size})"
                    else:
                        pattern_type = f"Even-length Incremental (Synthetic Fulcrum, Size {flip_size})"
                    
                    report_lines.extend([
                        f"**Pattern Type:** {pattern_type}",
                        f"**Segment Length:** {segment_length} genes",
                        f"**Flip Size:** {flip_size}",
                        ""
                    ])
                
                elif event_type == 'adjacency':
                    report_lines.extend([
                        f"**Pattern Type:** Adjacent Gene Swap",
                        f"**Genes Involved:** 2",
                        ""
                    ])
                
                # Gene list
                if len(genes) <= 20:
                    # Show all genes for small events
                    report_lines.extend([
                        "**Genes Involved:**",
                        "```"
                    ])
                    for i, gene in enumerate(genes):
                        pos = positions[i] if i < len(positions) else "N/A"
                        report_lines.append(f"{i+1:2d}. {gene} (position {pos})")
                    report_lines.extend(["```", ""])
                else:
                    # Show summary for large events
                    report_lines.extend([
                        f"**Genes Involved:** {len(genes)} genes",
                        "",
                        "**First 10 genes:**",
                        "```"
                    ])
                    for i in range(min(10, len(genes))):
                        gene = genes[i]
                        pos = positions[i] if i < len(positions) else "N/A"
                        report_lines.append(f"{i+1:2d}. {gene} (position {pos})")
                    
                    if len(genes) > 10:
                        report_lines.append(f"... and {len(genes) - 10} more genes")
                    
                    report_lines.extend(["```", ""])
                
                # Movement pattern (if available from debug info)
                # Ensure positions are numbers
                numeric_positions = [int(p) for p in positions if p is not None and str(p).isdigit()]
                if numeric_positions:
                    pos_range = f"{min(numeric_positions)}-{max(numeric_positions)}"
                else:
                    pos_range = "N/A"
                
                report_lines.extend([
                    "**Biological Effect:**",
                    f"- Reversed gene order in segment spanning positions {pos_range}",
                    f"- Reduced total movement by bringing genes closer to their target positions",
                    ""
                ])
        
        report_lines.extend(["---", ""])
    
    # Write report
    report_path = output_path / f"iteration_report_{genome1_name}_vs_{genome2_name}.md"
    
    with open(report_path, 'w') as f:
        f.write('\n'.join(report_lines))
    
    print(f"Detailed iteration report saved: {report_path}")
    return report_path

def save_pattern_analysis(results, output_path, genome1_name, genome2_name):
    """Save detailed pattern analysis as CSV"""
    
    # Ensure output_path is a Path object
    if isinstance(output_path, str):
        output_path = Path(output_path)
    
    pattern_data = []
    
    for chr_name, chr_data in results['distance_metrics']['chromosome_results'].items():
        if 'inversion_events' in chr_data:
            for event_idx, event in enumerate(chr_data['inversion_events']):
                
                event_type = event.get('type', 'unknown')
                genes = event.get('genes', [])
                positions = event.get('positions', [])
                iteration = event.get('iteration', 1)
                
                # Determine pattern specifics
                if event_type == 'flip':
                    flip_size = event.get('flip_indicator', 1)
                    segment_length = len(genes)
                    
                    if segment_length == 2:
                        pattern_type = "simple_adjacent"
                    elif segment_length % 2 == 1:
                        pattern_type = f"odd_incremental_size_{flip_size}"
                    else:
                        pattern_type = f"even_incremental_size_{flip_size}"
                else:
                    pattern_type = "adjacency_swap"
                    flip_size = 1
                    segment_length = 2
                
                # Ensure positions are numbers for calculations
                numeric_positions = [int(p) for p in positions if p is not None and str(p).isdigit()]
                
                pattern_data.append({
                    'chromosome': chr_name,
                    'iteration': iteration,
                    'event_number': event_idx + 1,
                    'pattern_type': pattern_type,
                    'event_type': event_type,
                    'segment_length': segment_length,
                    'flip_size': flip_size,
                    'genes_involved': len(genes),
                    'gene_list': ';'.join(genes),
                    'position_range': f"{min(numeric_positions)}-{max(numeric_positions)}" if numeric_positions else "N/A",
                    'start_position': min(numeric_positions) if numeric_positions else None,
                    'end_position': max(numeric_positions) if numeric_positions else None
                })
    
    # Save as CSV
    if pattern_data:
        pattern_df = pd.DataFrame(pattern_data)
        pattern_path = output_path / f"pattern_analysis_{genome1_name}_vs_{genome2_name}.csv"
        pattern_df.to_csv(pattern_path, index=False)
        print(f"Pattern analysis CSV saved: {pattern_path}")
        
        
        
        
        
        
        

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



def save_stage_data(data, stage_name: str, output_dir: Path, description: str = ""):
    """Save data from each stage as CSV with description"""
    stage_dir = output_dir / 'stages'
    stage_dir.mkdir(exist_ok=True)
    
    if isinstance(data, pd.DataFrame):
        filepath = stage_dir / f'{stage_name}.csv'
        data.to_csv(filepath, index=False)
        logger.info(f"Stage {stage_name}: Saved {len(data)} records to {filepath}")
        if description:
            logger.info(f"  Description: {description}")
    elif isinstance(data, dict):
        filepath = stage_dir / f'{stage_name}.json'
        with open(filepath, 'w') as f:
            json.dump(data, f, indent=2, default=str)
        logger.info(f"Stage {stage_name}: Saved dictionary to {filepath}")
        if description:
            logger.info(f"  Description: {description}")
    else:
        logger.warning(f"Stage {stage_name}: Unknown data type {type(data)}")
        
        
        


def save_hybrid_profile_data(hybrid_profile, output_path, description=""):
    """Save hybrid profile data as separate CSV files for cross-checking"""
    
    if 'positional_profile' in hybrid_profile:
        # Save positional profile
        positional_records = []
        for bin_id, genes_data in hybrid_profile['positional_profile'].items():
            for busco_id, profile_data in genes_data.items():
                # Skip summary entries - only process actual gene entries
                if busco_id in ['position_summary', 'rank_summary']:
                    continue
                
                # Get summary data for this bin
                summary_data = genes_data.get('position_summary', {})
                
                positional_records.append({
                    'bin_id': bin_id,
                    'busco_id': busco_id,
                    'average_percentage': summary_data.get('average_percentage', 0.0),
                    'percentage_range_min': summary_data.get('percentage_range', (0, 0))[0],
                    'percentage_range_max': summary_data.get('percentage_range', (0, 0))[1],
                    'genome_frequency': summary_data.get('genome_frequency', '0/0'),
                    'genome_count': summary_data.get('genome_count', 0),
                    'total_genomes': summary_data.get('total_genomes', 0),
                    'profile_type': 'positional'
                })
        
        positional_df = pd.DataFrame(positional_records)
        save_stage_data(positional_df, '6a_positional_profile', output_path, 
                       f"Positional profile: {len(positional_records)} bin-gene entries")
        
        # Save ordinal profile (similar fix)
        ordinal_records = []
        for rank_id, genes_data in hybrid_profile['ordinal_profile'].items():
            for busco_id, profile_data in genes_data.items():
                # Skip summary entries - only process actual gene entries
                if busco_id in ['position_summary', 'rank_summary']:
                    continue
                
                # Get summary data for this rank
                summary_data = genes_data.get('rank_summary', {})
                
                ordinal_records.append({
                    'rank_id': rank_id,
                    'busco_id': busco_id,
                    'frequency_percentage': summary_data.get('average_percentage', 0.0),
                    'genome_frequency': summary_data.get('genome_frequency', '0/0'),
                    'genome_count': summary_data.get('genome_count', 0),
                    'total_genomes': summary_data.get('total_genomes', 0),
                    'profile_type': 'ordinal'
                })
        
        ordinal_df = pd.DataFrame(ordinal_records)
        save_stage_data(ordinal_df, '6b_ordinal_profile', output_path,
                       f"Ordinal profile: {len(ordinal_records)} rank-gene entries")
        
        # Save hybrid summary
        save_stage_data(hybrid_profile['hybrid_summary'], '6c_hybrid_summary', output_path,
                       "Hybrid profile summary statistics")
    
    else:
        # Legacy profile format (unchanged)
        profile_records = []
        for bin_id, genes_data in hybrid_profile.items():
            for busco_id, profile_data in genes_data.items():
                profile_records.append({
                    'bin_id': bin_id,
                    'busco_id': busco_id,
                    'average_percentage': profile_data['average_percentage'],
                    'percentage_range_min': profile_data['percentage_range'][0],
                    'percentage_range_max': profile_data['percentage_range'][1],
                    'genome_frequency': profile_data['genome_frequency'],
                    'genome_count': profile_data['genome_count'],
                    'total_genomes': profile_data['total_genomes'],
                    'profile_type': 'legacy'
                })
        
        profile_df = pd.DataFrame(profile_records)
        save_stage_data(profile_df, '6_legacy_profile', output_path,
                       f"Legacy profile: {len(profile_records)} bin-gene entries")


def save_hybrid_assignments(hybrid_assignments, output_path, description=""):
    """Save hybrid bin assignments as CSV for cross-checking"""
    
    for genome_id, chromosomes in hybrid_assignments.items():
        all_assignment_records = []
        
        for chromosome, gene_assignments in chromosomes.items():
            for busco_id, hybrid_data in gene_assignments.items():
                # Save positional bins
                for bin_id, overlap_percentage in hybrid_data['positional_bins']:
                    all_assignment_records.append({
                        'busco_id': busco_id,
                        'chromosome': chromosome,
                        'assignment_type': 'positional',
                        'bin_or_rank_id': bin_id,
                        'overlap_or_frequency': overlap_percentage,
                        'gene_position': hybrid_data['gene_position'],
                        'relative_position': hybrid_data['relative_position']
                    })
                
                # Save ordinal window
                if hybrid_data['ordinal_window']:
                    all_assignment_records.append({
                        'busco_id': busco_id,
                        'chromosome': chromosome,
                        'assignment_type': 'ordinal',
                        'bin_or_rank_id': hybrid_data['ordinal_window'],
                        'overlap_or_frequency': 100.0,  # 100% presence in rank
                        'gene_position': hybrid_data['gene_position'],
                        'relative_position': hybrid_data['relative_position']
                    })
        
        assignment_df = pd.DataFrame(all_assignment_records)
        save_stage_data(
            assignment_df,
            f'5_hybrid_assignments_{genome_id}',
            output_path,
            f"Hybrid assignments for {genome_id}: {len(all_assignment_records)} assignments"
        )
