#!/usr/bin/env python3
"""
GIQ2 - Genome Inversion Quantifier
Main entry point.

Author: Zion Ayokunnu
Supervisors: Kamil Jaron, Sasha, Arif
Version: 1.0.0


Separate commands for genome inversion analysis:
1. build-profile: Build Markov profile from multiple training genomes
2. analyze-query: Analyze query genome against existing profile
"""

import sys
import argparse
import logging
from pathlib import Path
import pandas as pd
import json
from typing import Dict, List

from config.settings import CONFIG
from utils.file_utils import create_output_directory
from core import (
    parse_busco_table,
    filter_busco_genes,
    correct_strand_orientation,
    process_genomes_binning,
    compare_query_genome_to_profile,
    analyse_query_movements,
    get_movement_summary,
    check_events_iteration,
    probability_weighted_inversion_analysis,
    detect_flips
)

from contextual.metrics import (
    compute_inversion_rate_per_mb_busco,
    _analyse_gene_density_correlation
)


from visualisation.plots import (
    create_linearised_dotplot,
    create_busco_dotplot
)

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


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


def build_profile_command(busco_files: List[str], output_dir: str, config_overrides: Dict = None):

    config = CONFIG.copy()
    if config_overrides:
        config.update(config_overrides)
    
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    save_stage_data(config, '0_profile_config', output_path, "Profile building configuration")

    parsed_genomes = {}
    for i, busco_file in enumerate(busco_files):
        busco_path = Path(busco_file)
        genome_id = busco_path.stem  # Use filename as genome ID

        busco_df = parse_busco_table(str(busco_path), config)
        parsed_genomes[genome_id] = busco_df
        
        save_stage_data(
            busco_df, 
            f'1_parsed_{genome_id}', 
            output_path,
            f"Parsed BUSCO table for {genome_id}: {len(busco_df)} total genes"
        )
    filtered_genomes = {}
    for genome_id, busco_df in parsed_genomes.items():
        filtered_df = filter_busco_genes(busco_df, config)
        filtered_genomes[genome_id] = filtered_df
        
        save_stage_data(
            filtered_df,
            f'2_filtered_{genome_id}',
            output_path,
            f"Filtered BUSCO genes for {genome_id}: {len(filtered_df)} complete genes"
        )
    
    corrected_genomes = {}
    for genome_id, filtered_df in filtered_genomes.items():
        corrected_df = correct_strand_orientation(filtered_df)
        corrected_genomes[genome_id] = corrected_df
        
        save_stage_data(
            corrected_df,
            f'3_corrected_{genome_id}',
            output_path,
            f"Strand-corrected genes for {genome_id}"
        )

    all_bin_assignments = process_genomes_binning(corrected_genomes, config.get('position_bin_size_kb', 100))
    
    for genome_id, bin_assignments in all_bin_assignments.items():

        bin_records = []
        for busco_id, bin_overlaps in bin_assignments.items():
            for bin_id, overlap_percentage in bin_overlaps:
                bin_records.append({
                    'busco_id': busco_id,
                    'bin_id': bin_id,
                    'overlap_percentage': overlap_percentage
                })
        
        bin_df = pd.DataFrame(bin_records)
        save_stage_data(
            bin_df,
            f'4_bins_{genome_id}',
            output_path,
            f"Bin assignments for {genome_id}: {len(bin_records)} gene-bin mappings"
        )
    

    from core.profile import build_markov_profile
    
    markov_profile = build_markov_profile(all_bin_assignments, config.get('profile_calculation_method', 'average'))
    
    profile_records = []
    for bin_id, genes_data in markov_profile.items():
        for busco_id, profile_data in genes_data.items():
            profile_records.append({
                'bin_id': bin_id,
                'busco_id': busco_id,
                'average_percentage': profile_data['average_percentage'],
                'percentage_range_min': profile_data['percentage_range'][0],
                'percentage_range_max': profile_data['percentage_range'][1],
                'genome_frequency': profile_data['genome_frequency'],
                'genome_count': profile_data['genome_count'],
                'total_genomes': profile_data['total_genomes']
            })
    
    profile_df = pd.DataFrame(profile_records)
    save_stage_data(
        profile_df,
        '5_markov_profile',
        output_path,
        f"Markov profile: {len(profile_records)} gene-bin probability entries"
    )
    
    profile_data = {
        'markov_profile': markov_profile,
        'config': config,
        'training_genomes': list(parsed_genomes.keys()),
        'total_bins': len(markov_profile),
        'total_genes': len(set(r['busco_id'] for r in profile_records))
    }
    
    profile_file = output_path / 'markov_profile.json'
    with open(profile_file, 'w') as f:
        json.dump(profile_data, f, indent=2, default=str)
    
    logger.info("=" * 60)

    logger.info(f"Profile saved to: {profile_file}")
    logger.info(f"Training genomes: {', '.join(parsed_genomes.keys())}")
    logger.info(f"Total bins: {len(markov_profile)}")
    logger.info(f"Unique genes: {len(set(r['busco_id'] for r in profile_records))}")
    logger.info("=" * 60)
    
    return profile_data


def analyze_query_command(query_busco_file: str, profile_file: str, output_dir: str, config_overrides: Dict = None):

    logger.info("=" * 60)

    logger.info(f"Loading profile from: {profile_file}")
    with open(profile_file, 'r') as f:
        profile_data = json.load(f)
    
    markov_profile = profile_data['markov_profile']
    config = profile_data['config'].copy()
    
    if config_overrides:
        config.update(config_overrides)
    
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    query_path = Path(query_busco_file)
    query_id = query_path.stem

    logger.info(f"Profile training genomes: {', '.join(profile_data['training_genomes'])}")
    
    query_config = config.copy()
    query_config['query_genome'] = query_id
    query_config['profile_info'] = {
        'training_genomes': profile_data['training_genomes'],
        'total_bins': profile_data['total_bins'],
        'total_genes': profile_data['total_genes']
    }
    save_stage_data(query_config, '0_query_config', output_path, f"Query analysis configuration for {query_id}")
    
    query_busco_df = parse_busco_table(query_busco_file, config)
    save_stage_data(query_busco_df, f'1_parsed_{query_id}', output_path, f"Parsed query genome {query_id}")
    
    query_filtered_df = filter_busco_genes(query_busco_df, config)
    save_stage_data(query_filtered_df, f'2_filtered_{query_id}', output_path, f"Filtered query genome {query_id}")
    
    query_corrected_df = correct_strand_orientation(query_filtered_df)
    save_stage_data(query_corrected_df, f'3_corrected_{query_id}', output_path, f"Strand-corrected query genome {query_id}")
    
    query_bin_assignments = process_genomes_binning({query_id: query_corrected_df}, config.get('position_bin_size_kb', 100))
    
    query_bin_records = []
    for busco_id, bin_overlaps in query_bin_assignments[query_id].items():
        for bin_id, overlap_percentage in bin_overlaps:
            query_bin_records.append({
                'busco_id': busco_id,
                'bin_id': bin_id,
                'overlap_percentage': overlap_percentage
            })
    
    query_bin_df = pd.DataFrame(query_bin_records)
    save_stage_data(query_bin_df, f'4_bins_{query_id}', output_path, f"Query bin assignments for {query_id}")
    
    comparison_results = compare_query_genome_to_profile(query_bin_assignments[query_id], markov_profile)
    
    comparison_records = []
    for busco_id, result in comparison_results.items():
        comparison_records.append({
            'busco_id': busco_id,
            'query_bin': result['query_bin'],
            'query_overlap_percentage': result['query_overlap_percentage'],
            'expected_position': result['expected_position'],
            'position_deviation': result['position_deviation'],
            'standard_deviations': result['standard_deviations'],
            'position_specific_bit_score': result['bit_scores']['position_specific_bit_score'],
            'overall_match_bit_score': result['bit_scores']['overall_match_bit_score'],
            'e_value': result['bit_scores']['e_value']
        })
    
    comparison_df = pd.DataFrame(comparison_records)
    save_stage_data(comparison_df, f'5_comparison_{query_id}', output_path, f"Profile comparison for {query_id}")
    
    movement_results = analyse_query_movements(query_bin_assignments[query_id], markov_profile)

    movement_records = []
    for busco_id, result in movement_results.items():
        movement_records.append({
            'busco_id': busco_id,
            'current_ranges': str(result['current_ranges']),
            'target_ranges': str(result['target_ranges']),
            'mean_movement': result['movement_analysis']['mean_movement'],
            'total_pairs': result['movement_analysis']['total_pairs']
        })
    
    movement_df = pd.DataFrame(movement_records)
    save_stage_data(movement_df, f'6_movement_{query_id}', output_path, f"Movement analysis for {query_id}")
    
    movement_summary = get_movement_summary(movement_results)
    save_stage_data(movement_summary, f'7_movement_summary_{query_id}', output_path, f"Movement summary for {query_id}")
    
    inversion_analysis = check_events_iteration(movement_results)
    
    if inversion_analysis['inversion_events']:
        inversion_records = []
        for event in inversion_analysis['inversion_events']:
            inversion_records.append({
                'iteration': event['iteration'],
                'type': event['type'],
                'genes': str(event['genes']),
                'positions': str(event['positions']),
                'gene_inversions': event['gene_inversions']
            })
        
        inversion_events_df = pd.DataFrame(inversion_records)
        save_stage_data(inversion_events_df, f'8_inversion_events_{query_id}', output_path, f"Inversion events for {query_id}")
    
    inversion_summary = {
        'total_events': inversion_analysis['total_events'],
        'total_gene_inversions': inversion_analysis['total_gene_inversions'],
        'adjacency_events': inversion_analysis['adjacency_events'],
        'flip_events': inversion_analysis['flip_events'],
        'converged': inversion_analysis['converged'],
        'iterations': inversion_analysis['iterations']
    }
    save_stage_data(inversion_summary, f'9_inversion_summary_{query_id}', output_path, f"Inversion summary for {query_id}")
    
    probability_analysis = probability_weighted_inversion_analysis(movement_results, markov_profile)
    
    if probability_analysis['evaluated_steps']:
        prob_records = []
        for i, step in enumerate(probability_analysis['evaluated_steps']):
            prob_records.append({
                'step_index': i,
                'event_type': step['original_event']['type'],
                'genes': str(step['original_event']['genes']),
                'threshold_passed': step['probability_evaluation']['threshold_passed'],
                'overall_probability': step['probability_evaluation']['overall_probability'],
                'bit_score': step['probability_evaluation']['bit_score'],
                'failed_genes': str(step['probability_evaluation']['failed_genes'])
            })
        
        prob_analysis_df = pd.DataFrame(prob_records)
        save_stage_data(prob_analysis_df, f'10_probability_analysis_{query_id}', output_path, f"Probability analysis for {query_id}")
    
    if probability_analysis['alternative_pathways']:
        alt_records = []
        for gene_id, pathways in probability_analysis['alternative_pathways'].items():
            for i, pathway in enumerate(pathways):
                alt_records.append({
                    'gene_id': gene_id,
                    'pathway_index': i,
                    'target_position': pathway['target_position'],
                    'target_probability': pathway['target_probability'],
                    'combined_bit_score': pathway['evaluation']['combined_bit_score'],
                    'pathway_valid': pathway['evaluation']['pathway_valid']
                })
        
        alt_pathways_df = pd.DataFrame(alt_records)
        save_stage_data(alt_pathways_df, f'11_alternative_pathways_{query_id}', output_path, f"Alternative pathways for {query_id}")
    
    final_results = {
        'query_genome': query_id,
        'profile_info': profile_data,
        'movement_summary': movement_summary,
        'inversion_summary': inversion_summary,
        'probability_analysis_summary': {
            'total_standard_bit_score': probability_analysis['total_standard_bit_score'],
            'problematic_genes': probability_analysis['problematic_genes'],
            'has_alternatives': probability_analysis['has_alternatives']
        }
    }
    
    results_file = output_path / f'analysis_results_{query_id}.json'
    with open(results_file, 'w') as f:
        json.dump(final_results, f, indent=2, default=str)
    
    logger.info("=" * 60)
    logger.info("QUERY ANALYSIS COMPLETE")
    logger.info("=" * 60)
    logger.info(f"Query genome: {query_id}")
    logger.info(f"Total inversion events: {inversion_summary['total_events']}")
    logger.info(f"Total gene inversions: {inversion_summary['total_gene_inversions']}")
    logger.info(f"Converged: {inversion_summary['converged']}")
    logger.info(f"Problematic genes: {len(probability_analysis['problematic_genes'])}")
    logger.info(f"Results saved to: {results_file}")
    logger.info(f"All stage data saved to: {output_path / 'stages'}")
    
    return final_results














def save_results(inversion_df: pd.DataFrame,
                joined_df: pd.DataFrame,
                contextual_metrics: Dict,
                output_dir: Path,
                config: Dict):
    """Save all analysis results to files"""
    
    data_dir = output_dir / 'data'
    data_dir.mkdir(exist_ok=True)
    
    inversion_df.to_csv(data_dir / 'inversion_analysis.csv', index=False)
    joined_df.to_csv(data_dir / 'joined_gene_data.csv', index=False)
    
    import json
    with open(data_dir / 'contextual_metrics.json', 'w') as f:
        json.dump(contextual_metrics, f, indent=2, default=str)



def run_busco_inversion_analysis(config: Dict = None) -> Dict:
    """    
    Args:
        config: Configuration dictionary (uses default CONFIG if None)
        
    """
    
    if config is None:
        config = CONFIG
    
    output_dir = create_output_directory(config)
    logger.info(f"Output directory: {output_dir}")
    
    species1_name = config.get('first_species_name', 'Species1')
    species2_name = config.get('second_species_name', 'Species2')
    
    logger.info("\n" + "-" * 40)
    logger.info(f"Analyzing inversions between {species1_name} and {species2_name}")
    logger.info("-" * 40)
    
    try:

        busco_df1 = parse_busco_table(config['first_busco_path'], config)
        busco_df2 = parse_busco_table(config['second_busco_path'], config)
    
        
        filtered_df1 = filter_busco_genes(busco_df1, config)
        filtered_df2 = filter_busco_genes(busco_df2, config)
        
        logger.info("\n" + "-" * 40)
  
        inversion_df, joined_df = detect_flips(filtered_df1, filtered_df2, config)
        
        logger.info(f"Detected inversions across {len(inversion_df)} chromosome pairs")
        logger.info(f"Total genes compared: {len(joined_df)}")
        logger.info(f"Total flipped genes: {joined_df['is_flipped'].sum()}")
        
        
        contextual_metrics = compute_contextual_metrics(
            inversion_df, joined_df, config
        )
        
        logger.info("\n" + "-" * 40)

        
        visualisation_results = create_analysis_visualisations(
            inversion_df, joined_df, output_dir, config
        )
        
        save_results(
            inversion_df, joined_df, contextual_metrics, output_dir, config
        )
        
        report_path = generate_analysis_report(
            inversion_df, joined_df, contextual_metrics, 
            output_dir, config
        )
        
        results = {
            'inversion_df': inversion_df,
            'joined_df': joined_df,
            'contextual_metrics': contextual_metrics,
            'visualisation_results': visualisation_results,
            'output_dir': output_dir,
            'report_path': report_path,
            'config': config,
            'species_names': [species1_name, species2_name]
        }
        
        logger.info(f"Report available at: {report_path}")
        logger.info("=" * 80)
        
        return results
        
    except Exception as e:
        logger.error(f"Analysis failed: {str(e)}")
        if config.get('enable_debug_output', False):
            import traceback
            traceback.print_exc()
        raise


def compute_contextual_metrics(inversion_df: pd.DataFrame, 
                             joined_df: pd.DataFrame, 
                             config: Dict) -> Dict:
    
    contextual_metrics = {}
    
    rate_metrics = compute_inversion_rate_per_mb_busco(inversion_df, joined_df)
    contextual_metrics['inversion_rates'] = rate_metrics
    
    gene_density_corr = _analyse_gene_density_correlation(inversion_df, joined_df, config)
    contextual_metrics['gene_density_correlation'] = gene_density_corr

    
    return contextual_metrics


def create_analysis_visualisations(inversion_df: pd.DataFrame,
                                 joined_df: pd.DataFrame,
                                 output_dir: Path,
                                 config: Dict) -> Dict:
    
    plots_dir = output_dir / 'plots'
    plots_dir.mkdir(exist_ok=True)
    
    visualisation_results = {}
    
    species1_name = config.get('first_species_name', 'Species1')
    species2_name = config.get('second_species_name', 'Species2')
    
    import matplotlib.pyplot as plt
    plt.style.use('default')
    plt.rcParams['figure.dpi'] = config.get('dpi', 300)
    plt.rcParams['font.size'] = config.get('font_size', 12)
    
    try:
        try:
            create_busco_dotplot(joined_df, plots_dir, config)
            visualisation_results['dotplot'] = plots_dir / 'busco_dotplot.png'
        except Exception as e:
            logger.warning(f"create_busco_dotplot failed: {e}")
            
        try:
            fig, ax = create_linearised_dotplot(
                joined_df=joined_df, 
                plots_dir=plots_dir,
                config=config['dotplot_config']  # Pass the nested dotplot config
            )
            visualisation_results['linearized_dotplot'] = plots_dir / 'linearized_busco_dotplot.png'
            
        except Exception as e:
            logger.warning(f"Linearized dotplot failed: {e}")
            
    except Exception as e:
        logger.warning(f"No Vis worked! Ouch!")




def generate_analysis_report(inversion_df: pd.DataFrame,
                           joined_df: pd.DataFrame,
                           contextual_metrics: Dict,
                           output_dir: Path,
                           config: Dict) -> Path:
    
    reports_dir = output_dir / 'reports'
    reports_dir.mkdir(exist_ok=True)
    
    report_path = reports_dir / 'analysis_report.md'
    
    species1_name = config.get('first_species_name', 'Species1')
    species2_name = config.get('second_species_name', 'Species2')
    
    with open(report_path, 'w') as f:
        f.write("# GIQ1 - Genome Inversion Analysis Report\n\n")
        f.write(f"**Generated:** {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        f.write(f"**Species Comparison:** {species1_name} vs {species2_name}\n\n")
        
        f.write("## Analysis Summary\n\n")
        f.write(f"- **Total chromosome pairs analyzed:** {len(inversion_df)}\n")
        f.write(f"- **Pairs with inversions:** {len(inversion_df[inversion_df['flipped_genes'] > 0])}\n")
        f.write(f"- **Total genes compared:** {len(joined_df)}\n")
        f.write(f"- **Total flipped genes:** {joined_df['is_flipped'].sum()}\n")
        f.write(f"- **Overall flip rate:** {joined_df['is_flipped'].mean():.1%}\n\n")
        
        f.write("## Inversion Types Distribution\n\n")
        inversion_types = inversion_df['inversion_type'].value_counts()
        for inv_type, count in inversion_types.items():
            f.write(f"- **{inv_type}:** {count} pairs\n")
        f.write("\n")
        
        f.write("## Contextual Metrics\n\n")
        if 'inversion_rates' in contextual_metrics:
            rates = contextual_metrics['inversion_rates']
            f.write(f"- **Genome coverage (Mb):** {rates.get('genome_coverage_mb', 0):.2f}\n")
            f.write(f"- **Inversions per Mb:** {rates.get('inversion_events_per_mb', 0):.3f}\n")
            f.write(f"- **Flipped genes per Mb:** {rates.get('flipped_genes_per_mb', 0):.3f}\n")
        f.write("\n")
        
        f.write("## Chromosome Pair Details\n\n")
        f.write("| Chr1 | Chr2 | Total Genes | Flipped | Rate | Correlation | Type |\n")
        f.write("|------|------|-------------|---------|------|-------------|------|\n")
        
        for _, row in inversion_df.iterrows():
            f.write(f"| {row['chr1']} | {row['chr2']} | {row['total_genes']} | "
                   f"{row['flipped_genes']} | {row['flip_rate']:.2f} | "
                   f"{row.get('correlation', 'N/A')} | {row['inversion_type']} |\n")
        
        f.write(f"\n---\n*Report generated by GIQ1 v1.0.0*\n")
    
    logger.info(f"Analysis report saved to {report_path}")
    return report_path

























            

def main():
    """Main CLI entry point"""
    parser = argparse.ArgumentParser(description='Genome Inversion Analysis Tools')
    subparsers = parser.add_subparsers(dest='command', help='Available commands')
    
    # Build profile command
    profile_parser = subparsers.add_parser('build-profile', help='Build Markov profile from training genomes')
    profile_parser.add_argument('busco_files', nargs='+', help='BUSCO table files for training genomes')
    profile_parser.add_argument('-o', '--output', required=True, help='Output directory for profile')
    profile_parser.add_argument('--bin-size', type=int, default=100, help='Bin size in kb (default: 100)')
    profile_parser.add_argument('--method', choices=['average', 'range'], default='average', help='Profile calculation method')
    
    # Analyze query command
    query_parser = subparsers.add_parser('analyze-query', help='Analyze query genome against profile')
    query_parser.add_argument('query_busco', help='BUSCO table file for query genome')
    query_parser.add_argument('profile', help='Saved Markov profile JSON file')
    query_parser.add_argument('-o', '--output', required=True, help='Output directory for analysis')
    query_parser.add_argument('--threshold', type=float, default=0.5, help='Probability threshold (default: 0.5)')
    
    args = parser.parse_args()
    
    if args.command == 'build-profile':
        config_overrides = {
            'position_bin_size_kb': args.bin_size,
            'profile_calculation_method': args.method
        }
        
        try:
            profile_data = build_profile_command(args.busco_files, args.output, config_overrides)
            print(f"Profile building completed successfully!")
            print(f"Profile saved to: {args.output}/markov_profile.json")
            return 0
        except Exception as e:
            print(f"Profile building failed: {e}")
            return 1
    
    elif args.command == 'analyze-query':
        config_overrides = {
            'permutable_positions_threshold': args.threshold
        }
        
        try:
            results = analyze_query_command(args.query_busco, args.profile, args.output, config_overrides)
            print(f"Query analysis completed successfully!")
            print(f"Results saved to: {args.output}")
            return 0
        except Exception as e:
            print(f"Query analysis failed: {e}")
            return 1
    
    else:
        parser.print_help()
        return 1


if __name__ == "__main__":
    sys.exit(main())
