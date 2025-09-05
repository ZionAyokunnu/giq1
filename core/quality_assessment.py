import numpy as np
import logging
import time
from datetime import datetime
import pandas as pd
from pathlib import Path
from contextual.metrics import compute_inversion_rate_per_mb_busco, calculate_comprehensive_metrics

logger = logging.getLogger(__name__)


def assess_assembly_quality(fasta_path, busco_df, config):
    
    quality_metrics = {}
    # BUSCO completeness
    if len(busco_df) > 0:
        complete_buscos = len(busco_df[busco_df['status'] == 'Complete'])
        fragmented_buscos = len(busco_df[busco_df['status'] == 'Fragmented'])
        missing_buscos = len(busco_df[busco_df['status'] == 'Missing'])
        duplicated_buscos = len(busco_df[busco_df['status'] == 'Duplicated'])
        total_buscos = len(busco_df)
        
        completeness = complete_buscos / total_buscos if total_buscos > 0 else 0
        fragmentation = fragmented_buscos / total_buscos if total_buscos > 0 else 0
        duplication = duplicated_buscos / total_buscos if total_buscos > 0 else 0
        
        quality_metrics = {
            'busco_completeness': completeness,
            'busco_fragmentation': fragmentation,
            'busco_duplication': duplication,
            'busco_missing': missing_buscos / total_buscos if total_buscos > 0 else 0
        }

    
    
    logger.info(f"  Assembly quality: {quality_metrics}")
    
    return {
        'metrics': quality_metrics,
    }



def print_detailed_metrics(results, start_time=None, end_time=None):
    """Print comprehensive analysis metrics"""
    metrics = calculate_comprehensive_metrics(results, start_time, end_time)
    
    print(f"\nGENOMIC REARRANGEMENT ANALYSIS COMPLETE")
    print(f"=" * 50)
    
    # Summary stats
    summary = metrics['summary']
    print(f"SUMMARY:")
    print(f"  • Inversion events applied: {summary['total_inversion_events']}")
    print(f"  • Genes rearranged: {summary['total_genes_rearranged']}")
    print(f"  • Total genes processed: {summary['total_genes_processed']}")
    print(f"  • Algorithm iterations: {summary['algorithm_iterations']}")
    print(f"  • Chromosomes processed: {summary['chromosomes_processed']}")
    print(f"  • Convergence achieved: {summary['convergence_rate']}")
    print(f"  • Movement reduction: {summary['overall_movement_reduction']}")
    
    # Event breakdown
    events = metrics['event_breakdown']
    print(f"\nREARRANGEMENT OPERATIONS:")
    print(f"  • Adjacent gene swaps: {events['adjacency_inversions']}")
    print(f"  • Multi-gene segment flips: {events['segment_flips']}")
    print(f"  • Average genes per operation: {events['avg_genes_per_event']}")
    print(f"  • Operations per iteration: {events['events_per_iteration']}")
    
    # Efficiency metrics  
    efficiency = metrics['efficiency']
    print(f"\nALGORITHM EFFICIENCY:")
    print(f"  • Events per chromosome: {efficiency['avg_events_per_chromosome']}")
    print(f"  • Genome rearrangement density: {efficiency['rearrangement_density']}")
    print(f"  • Successful convergence rate: {efficiency['convergence_success_rate']}")
    print(f"  • Movement reduction per event: {efficiency['movement_reduction_efficiency']}")
    
    # Biological significance
    bio = metrics['biological_significance']
    print(f"\nBIOLOGICAL SIGNIFICANCE:")
    print(f"  • Perfectly converged chromosomes: {bio['perfectly_converged_chromosomes']}")
    print(f"  • Synteny preserved chromosomes: {bio['synteny_preserved_chromosomes']}")
    print(f"  • Highly rearranged chromosomes: {bio['highly_rearranged_chromosomes']}")
    print(f"  • Total movement eliminated: {bio['total_movement_eliminated']} units")
    print(f"  • Remaining movement: {bio['remaining_movement']} units")
    
    # Performance metrics
    if metrics['performance']:
        perf = metrics['performance']
        print(f"\nPERFORMANCE METRICS:")
        print(f"  • Execution time: {perf['execution_time_seconds']:.2f} seconds")
        print(f"  • Processing rate: {perf['genes_per_second']:.1f} genes/sec")
        print(f"  • Event rate: {perf['events_per_second']:.2f} events/sec")
        print(f"  • Chromosome rate: {perf['chromosomes_per_second']:.2f} chr/sec")
    
    # Per-chromosome breakdown
    print(f"\nPER-CHROMOSOME RESULTS:")
    for chr_name, chr_data in metrics['chromosome_details'].items():
        status = "CONVERGED" if chr_data['convergence_achieved'] else "PARTIAL"
        print(f"  • {chr_name}: {chr_data['genes_processed']} genes, {chr_data['events_applied']} events, {chr_data['iterations_needed']} iterations ({status})")
        print(f"    └─ Movement: {chr_data['initial_movement']} → {chr_data['final_movement']} ({chr_data['reduction_percentage']:.1f}% reduction)")
        if chr_data['remaining_unresolved'] > 0:
            print(f"    └─ {chr_data['remaining_unresolved']} genes still need rearrangement")
            
            
            
            
            
