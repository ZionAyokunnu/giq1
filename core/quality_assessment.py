"""
Quality assessment module for the Genome Inversion analyser1

"""

import numpy as np
import logging
import time
from datetime import datetime
import pandas as pd
from pathlib import Path


logger = logging.getLogger(__name__)


def assess_assembly_quality(fasta_path, busco_df, config):
    """Comprehensive assembly quality assessment"""
    
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
            
            
            
            
            
            
            
            
            
            
            #MARKDOWN
            
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
        
        
        