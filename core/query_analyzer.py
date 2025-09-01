"""
Query Analyzer Module
Handles query genome analysis against Markov profiles
"""

import pandas as pd
from pathlib import Path
import logging
import json
from typing import Dict

from config.settings import CONFIG
from core import (
    parse_busco_table,
    filter_busco_genes,
    compare_query_genome_to_profile,
    analyse_query_movements,
    get_movement_summary,
    load_chromosome_mappings,
    group_genomes_by_chromosome,
    process_genomes_binning_hybrid,
    check_events_iteration,
    probability_weighted_inversion_analysis,
    load_chromosome_mappings,
    group_genomes_by_chromosome,
    process_genomes_binning_hybrid,  
    save_stage_data,
    standardize_chromosome_names_in_binning,
)

logger = logging.getLogger(__name__)

def analyze_query_command_hybrid(query_busco_file: str, profile_file: str, output_dir: str, config_overrides: Dict = None):
    """
    Analyze query genome against hybrid profile (mostly unchanged logic)
    """
    logger.info("=" * 60)
    logger.info("ANALYZING QUERY GENOME AGAINST HYBRID PROFILE")
    logger.info("=" * 60)

    # Load profile (now supports hybrid format)
    logger.info(f"Loading profile from: {profile_file}")
    with open(profile_file, 'r') as f:
        profile_data = json.load(f)
    
    # Extract the appropriate profile format
    if 'hybrid_profile' in profile_data:
        profile = profile_data['hybrid_profile']
        profile_type = profile_data.get('profile_type', 'unknown')
        logger.info(f"Loaded {profile_type} profile")
        
     
        markov_profile = profile
        logger.info(f"Using full {profile_type} profile for hybrid movement analysis")
    else:
        # Legacy format
        markov_profile = profile_data['markov_profile']
        profile_type = 'legacy'
        logger.info("Loaded legacy profile")
    
    config = profile_data['config'].copy()
    if config_overrides:
        config.update(config_overrides)
    
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    query_path = Path(query_busco_file)
    query_id = query_path.stem

    logger.info(f"Query genome: {query_id}")
    logger.info(f"Profile training genomes: {', '.join(profile_data['training_genomes'])}")
    
    # Continue with existing query analysis pipeline...
    # (The rest of the analysis remains the same since we're using positional profile)
    
    query_config = config.copy()
    query_config['query_genome'] = query_id
    query_config['profile_info'] = {
        'training_genomes': profile_data['training_genomes'],
        'profile_type': profile_type
    }
    save_stage_data(query_config, '0_query_config', output_path, f"Query analysis configuration for {query_id}")
    
    # Parse and process query genome (standard pipeline)
    query_busco_df = parse_busco_table(query_busco_file, config)
    save_stage_data(query_busco_df, f'1_parsed_{query_id}', output_path, f"Parsed query genome {query_id}")
    
    query_filtered_df = filter_busco_genes(query_busco_df, config)
    save_stage_data(query_filtered_df, f'2_filtered_{query_id}', output_path, f"Filtered query genome {query_id}")
    
    save_stage_data(query_filtered_df, f'3_corrected_{query_id}', output_path, f"Strand-corrected query genome {query_id}")
    
    query_grouped = group_genomes_by_chromosome({query_id: query_filtered_df})
    query_chromosomes = query_grouped[query_id]
    
    # Apply chromosome standardization if needed
    if profile_data.get('chromosome_mappings_used', False):
        chr_map_file = config.get('chromosome_mappings_file')
        if chr_map_file and Path(chr_map_file).exists():
            logger.info("Applying chromosome standardization to query genome")
            chromosome_mappings = load_chromosome_mappings(chr_map_file)
            query_standardized = standardize_chromosome_names_in_binning(
                {query_id: query_chromosomes}, chromosome_mappings
            )
            query_chromosomes = query_standardized[query_id]
    
    # Debug: Check chromosome names after standardization
    logger.info("Chromosome names after standardization:")
    for chr_name in query_chromosomes.keys():
        logger.info(f"  {chr_name}")
    
    # For now, use standard binning (future: hybrid query binning)

    query_bin_assignments = process_genomes_binning_hybrid(
        {query_id: query_chromosomes}, 
        config.get('position_bin_size_kb', 100)
    )
    
    query_bin_records = []
    for chromosome, bin_assignments in query_bin_assignments[query_id].items():
        for busco_id, hybrid_data in bin_assignments.items():
            # Process positional bins
            for bin_id, overlap_percentage in hybrid_data['positional_bins']:
                query_bin_records.append({
                    'busco_id': busco_id,
                    'bin_id': bin_id,
                    'overlap_percentage': overlap_percentage
                })
            
            # Process ordinal data if available
            if hybrid_data['ordinal_window']:
                query_bin_records.append({
                    'busco_id': busco_id,
                    'bin_id': hybrid_data['ordinal_window'],
                    'overlap_percentage': 100.0  # Ordinal represents discrete rank position
                })

    query_bin_df = pd.DataFrame(query_bin_records)
    save_stage_data(query_bin_df, f'4_bins_{query_id}', output_path, f"Query bin assignments for {query_id}")
    
    comparison_results = compare_query_genome_to_profile(query_bin_assignments[query_id], markov_profile)

    # debug:
    print("DEBUG - comparison_results structure:")
    print(f"Type: {type(comparison_results)}")
    for key, value in list(comparison_results.items())[:2]:
        print(f"Key: {key}, Value type: {type(value)}")
        if isinstance(value, dict):
            for subkey in list(value.keys())[:2]:
                print(f"  Subkey: {subkey}")
    
    comparison_records = []
    for chromosome, gene_results in comparison_results.items():
        for busco_id, result in gene_results.items():
            comparison_records.append({
                'busco_id': busco_id,
                'chromosome': chromosome,
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
    
    # Extract the query genome data from the nested structure
    query_genome_data = query_bin_assignments[query_id]
    
    # Debug: Check the structure
    print("DEBUG - query_genome_data structure:")
    print(f"Type: {type(query_genome_data)}")
    if isinstance(query_genome_data, dict):
        for key, value in list(query_genome_data.items())[:2]:
            print(f"Key: {key}, Value type: {type(value)}")
            if isinstance(value, dict):
                for subkey in list(value.keys())[:2]:
                    print(f"  Subkey: {subkey}")
    
    # Convert the hybrid data structure to the format expected by movement analysis
    # The data is already in the right structure, just need to convert hybrid format to legacy
    converted_assignments = {}
    for chromosome, gene_assignments in query_genome_data.items():
        converted_assignments[chromosome] = {}
        for busco_id, hybrid_data in gene_assignments.items():
            # Start with positional bins
            combined_bins = hybrid_data['positional_bins'].copy()
            
            # Add ordinal data as a "bin" if available
            if hybrid_data['ordinal_window']:
                # Convert ordinal window to a pseudo-bin entry
                # Use 100% overlap since ordinal represents discrete rank position
                combined_bins.append((hybrid_data['ordinal_window'], 100.0))
            
            # If no positional bins but has ordinal, ensure we have data
            if not combined_bins and hybrid_data['ordinal_window']:
                combined_bins = [(hybrid_data['ordinal_window'], 100.0)]
            
            converted_assignments[chromosome][busco_id] = combined_bins
    
    print("DEBUG - About to call analyse_query_movements")
    print(f"DEBUG - converted_assignments type: {type(converted_assignments)}")
    print(f"DEBUG - converted_assignments keys: {list(converted_assignments.keys())[:3]}")
    
    # Debug: Check the profile structure
    print("DEBUG - markov_profile structure:")
    print(f"Type: {type(markov_profile)}")
    if isinstance(markov_profile, dict):
        print(f"Keys: {list(markov_profile.keys())[:5]}")
        if 'positional_profile' in markov_profile:
            sample_bin = next(iter(markov_profile['positional_profile'].keys()))
            sample_genes = markov_profile['positional_profile'][sample_bin]
            print(f"Sample bin {sample_bin} structure:")
            if sample_genes:
                sample_gene = next(iter(sample_genes.keys()))
                print(f"  Sample gene {sample_gene}: {sample_genes[sample_gene]}")
    
    try:
        movement_results, structural_variations = analyse_query_movements(converted_assignments, markov_profile, input_format='legacy')
        print("DEBUG - analyse_query_movements completed successfully")
    except Exception as e:
        print(f"DEBUG - Error in analyse_query_movements: {e}")
        import traceback
        traceback.print_exc()
        raise
    
    movement_records = []
    for chromosome, gene_results in movement_results.items():
        for busco_id, result in gene_results.items():
            movement_records.append({
                'busco_id': busco_id,
                'chromosome': chromosome,  # Add chromosome info
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
    
    # logger.info("Query analysis complete - using positional profile")
    # logger.info("TODO: Add ordinal movement analysis in future update")
    
    # return {"status": "completed", "profile_type": profile_type}

