"""
BUSCO Contextual Metrics
"""

import pandas as pd
import numpy as np
import logging
from pathlib import Path
from typing import Dict, List, Any, Optional, Union
from collections import defaultdict

from config import (
    CONFIG
)

logger = logging.getLogger(__name__)


def compute_inversion_rate_per_mb_busco(inversion_df: pd.DataFrame, #Using Species A as the reference mb span
                                       joined_df: pd.DataFrame) -> Dict:

        
    total_inversions = len(inversion_df[inversion_df['flipped_genes'] > 0])
    total_flipped_genes = inversion_df['flipped_genes'].sum()
    total_genes = len(joined_df)
    
    genome_coverage_estimate = 0
    for _, row in inversion_df.iterrows():
        if row['total_genes'] > 1:
            chr_genes = joined_df[
                (joined_df['chr1'] == row['chr1']) & 
                (joined_df['chr2'] == row['chr2'])
            ]
            if len(chr_genes) > 0:
                chr1_span = chr_genes['end1'].max() - chr_genes['start1'].min()
                genome_coverage_estimate += chr1_span
    
    genome_size_mb = genome_coverage_estimate / 1_000_000 if genome_coverage_estimate > 0 else 0            #Print genome span for species 1 (reference Genome)
    
    if genome_size_mb == 0:
        logger.warning("Cannot estimate genome size from BUSCO data")
        return {}
    

    inversions_per_mb = total_inversions / genome_size_mb                           #Print
    flipped_genes_per_mb = total_flipped_genes / genome_size_mb                    #Print
    

    rates_by_type = {}
    if 'inversion_type' in inversion_df.columns:
        for inv_type in inversion_df['inversion_type'].unique():
            type_data = inversion_df[inversion_df['inversion_type'] == inv_type]
            type_inversions = len(type_data[type_data['flipped_genes'] > 0])
            rates_by_type[inv_type] = type_inversions / genome_size_mb              #Print
    

    rates_by_chr = {}
    for chrom in inversion_df['chr1'].unique():
        chr_data = inversion_df[inversion_df['chr1'] == chrom]
        chr_inversions = len(chr_data[chr_data['flipped_genes'] > 0])
        rates_by_chr[chrom] = chr_inversions / genome_size_mb if genome_size_mb > 0 else 0      #Print
    
    rate_metrics = {
        'total_chromosome_pairs_with_inversions': total_inversions,
        'total_flipped_genes': total_flipped_genes,
        'total_analysed_genes': total_genes,
        'genome_coverage_mb': genome_size_mb,
        'inversion_events_per_mb': inversions_per_mb,
        'flipped_genes_per_mb': flipped_genes_per_mb,
        'rates_by_type': rates_by_type,
        'rates_by_chromosome': rates_by_chr,
        'global_flip_rate': total_flipped_genes / total_genes if total_genes > 0 else 0
    }
    
    logger.info(f"Inversion rate: {inversions_per_mb:.3f} events per Mb, "
               f"{flipped_genes_per_mb:.3f} flipped genes per Mb")
    return rate_metrics




def _analyse_gene_density_correlation(inversion_df: pd.DataFrame, 
                                    joined_df: pd.DataFrame, config) -> Dict:
    """
    analyse correlation between gene density (per mb) and inversion density
    """
    
    species1_name = config.get('first_species_name', 'Species A')
    species2_name = config.get('second_species_name', 'Species B')

    species_a_gene_density = {}  
    species_b_gene_density = {} 

    for chrom in joined_df['chr1'].unique():
        chr_genes = joined_df[joined_df['chr1'] == chrom]
        
        if len(chr_genes) > 1:
            chr_span = chr_genes['end1'].max() - chr_genes['start1'].min()
            density = len(chr_genes) / (chr_span / 1000000) if chr_span > 0 else 0
            species_a_gene_density[chrom] = density             #Print
    
    for chrom in joined_df['chr2'].unique():
        chr_genes = joined_df[joined_df['chr2'] == chrom]
        
        if len(chr_genes) > 1:
            chr_span = chr_genes['end2'].max() - chr_genes['start2'].min()
            density = len(chr_genes) / (chr_span / 1000000) if chr_span > 0 else 0
            species_b_gene_density[chrom] = density
    
    species_a_inversion_density = {}  # Inversions of A relative to B
    species_b_inversion_density = {}  # Inversions of B relative to A
    
    for chrom in species_a_gene_density.keys():
        chr_inv_pairs = inversion_df[inversion_df['chr1'] == chrom]
        total_flipped = chr_inv_pairs['flipped_genes'].sum()
        
        chr_genes = joined_df[joined_df['chr1'] == chrom]
        chr_span = chr_genes['end1'].max() - chr_genes['start1'].min()
        span_mb = chr_span / 1000000 if chr_span > 0 else 1
        
        species_a_inversion_density[chrom] = total_flipped / span_mb
    

    for chrom in species_b_gene_density.keys():
        chr_inv_pairs = inversion_df[inversion_df['chr2'] == chrom]
        total_flipped = chr_inv_pairs['flipped_genes'].sum()
        
        chr_genes = joined_df[joined_df['chr2'] == chrom]
        chr_span = chr_genes['end2'].max() - chr_genes['start2'].min()
        span_mb = chr_span / 1000000 if chr_span > 0 else 1
        
        species_b_inversion_density[chrom] = total_flipped / span_mb
    
    correlations = {}
    
    if len(species_a_gene_density) > 2:
        common_chroms_a = set(species_a_gene_density.keys()) & set(species_a_inversion_density.keys())
        if len(common_chroms_a) > 2:
            gene_vals_a = [species_a_gene_density[chrom] for chrom in common_chroms_a]
            inv_vals_a = [species_a_inversion_density[chrom] for chrom in common_chroms_a]
            
            if np.std(gene_vals_a) > 0 and np.std(inv_vals_a) > 0:
                corr_a = np.corrcoef(gene_vals_a, inv_vals_a)[0, 1]
                correlations['species_a_within_genome'] = {
                    'correlation': corr_a,
                    'chromosomes': list(common_chroms_a),
                    'description': '{species1_name} gene density vs {species1_name} inversion density'
                }
    
    if len(species_b_gene_density) > 2:
        common_chroms_b = set(species_b_gene_density.keys()) & set(species_b_inversion_density.keys())
        if len(common_chroms_b) > 2:
            gene_vals_b = [species_b_gene_density[chrom] for chrom in common_chroms_b]
            inv_vals_b = [species_b_inversion_density[chrom] for chrom in common_chroms_b]
            
            if np.std(gene_vals_b) > 0 and np.std(inv_vals_b) > 0:
                corr_b = np.corrcoef(gene_vals_b, inv_vals_b)[0, 1]
                correlations['species_b_within_genome'] = {
                    'correlation': corr_b,
                    'chromosomes': list(common_chroms_b),
                    'description': '{species2_name} gene density vs {species2_name} inversion density'
                }

    
    result = {
        'gene_density_species_a': species_a_gene_density,
        'gene_density_species_b': species_b_gene_density,
        'inversion_density_species_a': species_a_inversion_density,
        'inversion_density_species_b': species_b_inversion_density,
    }
    
    return result



def calculate_comprehensive_metrics(results, start_time=None, end_time=None):
    """Calculate detailed metrics for genomic rearrangement analysis"""
    execution_time = (end_time - start_time) if start_time and end_time else None
    
    # Basic counts - ensure they are numbers
    total_events = int(results['distance_metrics']['total_events'])
    total_genes = int(results['distance_metrics']['total_gene_inversions'])
    iterations = int(results['distance_metrics']['iterations'])
    
    # Calculate per-chromosome metrics and aggregate event types
    chr_results = results['distance_metrics']['chromosome_results']
    chr_metrics = {}
    
    total_initial_movement = 0
    total_final_movement = 0
    
    # CALCULATE EVENT TYPES FROM CHROMOSOME DATA INSTEAD OF EXPECTING PRE-AGGREGATED
    total_adjacency_events = 0
    total_flip_events = 0
    
    for chr_name, chr_data in chr_results.items():
        # Safely get event type counts from each chromosome
        chr_adjacency = int(chr_data.get('adjacency_events', 0))
        chr_flip = int(chr_data.get('flip_events', 0))
        
        total_adjacency_events += chr_adjacency
        total_flip_events += chr_flip
        
        # Handle movement calculations with missing data detection
        initial_movement = float(chr_data.get('initial_total_movement', 0))
        final_movement = float(chr_data.get('final_total_movement', 0))
        events_applied = int(chr_data.get('total_events', 0))
        genes_rearranged = int(chr_data.get('total_gene_inversions', 0))
        iterations_needed = int(chr_data.get('iterations', 0))
        remaining_unresolved = int(chr_data.get('final_non_zero_movements', 0))
        
        # Fix missing initial movement data
        if initial_movement == 0 and (events_applied > 0 or final_movement > 0):
            if final_movement > 0:
                estimated_initial = final_movement + (genes_rearranged * 2)
                print(f"WARNING: Missing initial movement for {chr_name}, estimating {estimated_initial}")
                initial_movement = float(estimated_initial)
            elif events_applied > 0:
                estimated_initial = float(genes_rearranged * 3)
                print(f"WARNING: Missing initial movement for {chr_name}, estimating {estimated_initial} from {genes_rearranged} genes rearranged")
                initial_movement = estimated_initial
        
        movement_reduction = initial_movement - final_movement
        reduction_percentage = (movement_reduction / initial_movement * 100) if initial_movement > 0 else (100.0 if final_movement == 0 else 0.0)
        
        chr_metrics[chr_name] = {
            'genes_processed': len(chr_data.get('final_sequence', [])),
            'events_applied': events_applied,
            'genes_rearranged': genes_rearranged,
            'iterations_needed': iterations_needed,
            'convergence_achieved': chr_data.get('converged', False),
            'initial_movement': initial_movement,
            'final_movement': final_movement,
            'movement_reduction': movement_reduction,
            'reduction_percentage': reduction_percentage,
            'remaining_unresolved': remaining_unresolved,
            'adjacency_events': chr_adjacency,  # Include for transparency
            'flip_events': chr_flip
        }
        
        total_initial_movement += initial_movement
        total_final_movement += final_movement
    
    # Now we have correct totals calculated from actual data
    adjacency_events = total_adjacency_events
    flip_events = total_flip_events
    
    # Verify our calculation matches the expected total
    calculated_total = adjacency_events + flip_events
    if calculated_total != total_events:
        print(f"WARNING: Event count mismatch - calculated {calculated_total}, expected {total_events}")
    
    # Rest of calculations...
    total_genes_processed = sum(cm['genes_processed'] for cm in chr_metrics.values())
    convergence_rate = sum(1 for cm in chr_metrics.values() if cm['convergence_achieved']) / len(chr_metrics) if chr_metrics else 0
    avg_events_per_chromosome = total_events / len(chr_metrics) if chr_metrics else 0
    
    genes_per_event = total_genes / total_events if total_events > 0 else 0
    events_per_iteration = total_events / iterations if iterations > 0 else 0
    total_movement_reduction = total_initial_movement - total_final_movement
    overall_reduction_percentage = (total_movement_reduction / total_initial_movement * 100) if total_initial_movement > 0 else 100
    
    perfectly_converged = sum(1 for cm in chr_metrics.values() if cm['final_movement'] == 0)
    conserved_chromosomes = sum(1 for cm in chr_metrics.values() if cm['events_applied'] == 0)
    highly_rearranged = sum(1 for cm in chr_metrics.values() if cm['events_applied'] >= 5)
    
    # Performance metrics
    performance_metrics = {}
    if execution_time and execution_time > 0.001:
        exec_time = float(execution_time)
        performance_metrics = {
            'execution_time_seconds': exec_time,
            'genes_per_second': total_genes_processed / exec_time,
            'events_per_second': total_events / exec_time,
            'chromosomes_per_second': len(chr_metrics) / exec_time
        }
    
    metrics = {
        'summary': {
            'total_inversion_events': total_events,
            'total_genes_rearranged': total_genes,
            'total_genes_processed': total_genes_processed,
            'algorithm_iterations': iterations,
            'chromosomes_processed': len(chr_metrics),
            'convergence_rate': f"{convergence_rate:.1%}",
            'overall_movement_reduction': f"{overall_reduction_percentage:.1f}%"
        },
        'event_breakdown': {
            'adjacency_inversions': adjacency_events,  # Now calculated correctly
            'segment_flips': flip_events,              # Now calculated correctly
            'avg_genes_per_event': f"{genes_per_event:.1f}",
            'events_per_iteration': f"{events_per_iteration:.1f}"
        },
        'efficiency': {
            'avg_events_per_chromosome': f"{avg_events_per_chromosome:.1f}",
            'rearrangement_density': f"{(total_genes/total_genes_processed*100):.1f}%" if total_genes_processed > 0 else "0%",
            'convergence_success_rate': f"{convergence_rate:.1%}",
            'movement_reduction_efficiency': f"{(total_movement_reduction/max(total_events,1)):.1f} units/event"
        },
        'biological_significance': {
            'perfectly_converged_chromosomes': perfectly_converged,
            'synteny_preserved_chromosomes': conserved_chromosomes,
            'highly_rearranged_chromosomes': highly_rearranged,
            'total_movement_eliminated': int(total_movement_reduction),
            'remaining_movement': int(total_final_movement)
        },
        'performance': performance_metrics,
        'chromosome_details': chr_metrics
    }
    
    return metrics     
            
        
