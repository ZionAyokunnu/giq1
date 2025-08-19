#!/usr/bin/env python3
"""
GIQ2 - Genome Inversion Quantifier
Main entry point.

Author: Zion Ayokunnu
Supervisors: Kamil Jaron, Sasha, Arif
Version: 1.0.0

Separate commands for genome inversion analysis:
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
    compare_query_genome_to_profile,
    analyse_query_movements,
    get_movement_summary,
    check_events_iteration,
    probability_weighted_inversion_analysis,
    detect_flips,
    load_chromosome_mappings,
    standardize_chromosome_name_unified,
    correct_chromosome_orientation,
    group_genomes_by_chromosome,
    process_genomes_binning_hybrid,
    build_markov_profile,
    get_genome_name_from_fasta,
    run_all_pairwise_alignments,
    create_unified_mappings_multi_reference,
    save_unified_mappings
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


def save_hybrid_profile_data(hybrid_profile, output_path, description=""):
    """Save hybrid profile data as separate CSV files for cross-checking"""
    
    if 'positional_profile' in hybrid_profile:
        # Save positional profile
        positional_records = []
        for bin_id, genes_data in hybrid_profile['positional_profile'].items():
            for busco_id, profile_data in genes_data.items():
                positional_records.append({
                    'bin_id': bin_id,
                    'busco_id': busco_id,
                    'average_percentage': profile_data['average_percentage'],
                    'percentage_range_min': profile_data['percentage_range'][0],
                    'percentage_range_max': profile_data['percentage_range'][1],
                    'genome_frequency': profile_data['genome_frequency'],
                    'genome_count': profile_data['genome_count'],
                    'total_genomes': profile_data['total_genomes'],
                    'profile_type': 'positional'
                })
        
        positional_df = pd.DataFrame(positional_records)
        save_stage_data(positional_df, '6a_positional_profile', output_path, 
                       f"Positional profile: {len(positional_records)} bin-gene entries")
        
        # Save ordinal profile
        ordinal_records = []
        for rank_id, genes_data in hybrid_profile['ordinal_profile'].items():
            for busco_id, profile_data in genes_data.items():
                ordinal_records.append({
                    'rank_id': rank_id,
                    'busco_id': busco_id,
                    'frequency_percentage': profile_data['frequency_percentage'],
                    'genome_frequency': profile_data['genome_frequency'],
                    'genome_count': profile_data['genome_count'],
                    'total_genomes': profile_data['total_genomes'],
                    'profile_type': 'ordinal'
                })
        
        ordinal_df = pd.DataFrame(ordinal_records)
        save_stage_data(ordinal_df, '6b_ordinal_profile', output_path,
                       f"Ordinal profile: {len(ordinal_records)} rank-gene entries")
        
        # Save hybrid summary
        save_stage_data(hybrid_profile['hybrid_summary'], '6c_hybrid_summary', output_path,
                       "Hybrid profile summary statistics")
    
    else:
        # Legacy profile format
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


def standardize_chromosome_names_in_binning(grouped_genomes, chromosome_mappings):
    """Apply chromosome name standardization to grouped genomes before binning"""
    standardized_genomes = {}
    
    for genome_id, chromosomes in grouped_genomes.items():
        standardized_chromosomes = {}
        
        for original_chr_name, genes_df in chromosomes.items():
            standard_chr_name = standardize_chromosome_name_unified(
                chromosome_mappings, genome_id, original_chr_name
            )
            logger.info(f"  {genome_id}: {original_chr_name} → {standard_chr_name}")
            
            standardized_df = genes_df.copy()
            standardized_df['sequence'] = standard_chr_name
            
            if standard_chr_name in standardized_chromosomes:
                standardized_chromosomes[standard_chr_name] = pd.concat([
                    standardized_chromosomes[standard_chr_name], 
                    standardized_df
                ], ignore_index=True)
            else:
                standardized_chromosomes[standard_chr_name] = standardized_df
        
        standardized_genomes[genome_id] = standardized_chromosomes
    
    return standardized_genomes


def build_profile_command_hybrid(busco_files: List[str], output_dir: str, config_overrides: Dict = None, chr_map_file: str = None):
    """
    Build hybrid Markov profile from multiple training genomes
    """
    logger.info("=" * 60)
    logger.info("BUILDING HYBRID MARKOV PROFILE")
    logger.info("=" * 60)
    
    config = CONFIG.copy()
    if config_overrides:
        config.update(config_overrides)
    
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # Load chromosome mappings if provided
    chromosome_mappings = None
    if chr_map_file:
        logger.info(f"Loading chromosome mappings from: {chr_map_file}")
        chromosome_mappings = load_chromosome_mappings(chr_map_file)
        logger.info(f"Loaded mappings for {len(chromosome_mappings['genome_mappings'])} genomes")
        config['chromosome_mappings_file'] = chr_map_file

    save_stage_data(config, '0_profile_config', output_path, "Hybrid profile building configuration")

    # Stage 1: Parse BUSCO files
    logger.info("Step 1: Parsing BUSCO files")
    parsed_genomes = {}
    for busco_file in busco_files:
        busco_path = Path(busco_file)
        genome_id = busco_path.stem
        
        busco_df = parse_busco_table(str(busco_path), config)
        parsed_genomes[genome_id] = busco_df
        
        save_stage_data(
            busco_df, 
            f'1_parsed_{genome_id}', 
            output_path,
            f"Parsed BUSCO table for {genome_id}: {len(busco_df)} total genes"
        )

    # Stage 2: Filter genes
    logger.info("Step 2: Filtering BUSCO genes")
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
    
    # Stage 3: Correct strand orientation
    logger.info("Step 3: Correcting strand orientations")
    corrected_genomes = {}
    for genome_id, filtered_df in filtered_genomes.items():
        corrected_df = correct_chromosome_orientation(filtered_df)
        corrected_genomes[genome_id] = corrected_df
        
        save_stage_data(
            corrected_df,
            f'3_corrected_{genome_id}',
            output_path,
            f"Strand-corrected genes for {genome_id}"
        )
    
    # Stage 4: Group by chromosome
    logger.info("Step 4: Grouping by chromosomes")
    grouped_genomes = group_genomes_by_chromosome(corrected_genomes)
    
    # Stage 4b: Apply chromosome standardization if mappings provided
    if chromosome_mappings:
        logger.info("Step 4b: Applying chromosome name standardization")
        grouped_genomes = standardize_chromosome_names_in_binning(grouped_genomes, chromosome_mappings)
        
        chr_info = {}
        for genome_id, chromosomes in grouped_genomes.items():
            chr_info[genome_id] = {
                'standardized_chromosomes': list(chromosomes.keys()),
                'gene_counts_per_chr': {chr_name: len(genes_df) for chr_name, genes_df in chromosomes.items()}
            }
        
        save_stage_data(
            chr_info,
            '4b_standardized_chromosomes',
            output_path,
            "Standardized chromosome names and gene distributions"
        )

    # Stage 5: Process hybrid binning
    logger.info("Step 5: Processing hybrid binning (positional + ordinal)")
    hybrid_assignments = process_genomes_binning_hybrid(
        grouped_genomes, 
        config.get('position_bin_size_kb', 100),
        chromosome_mappings
    )
    
    # Save hybrid assignments as CSV for cross-checking
    save_hybrid_assignments(hybrid_assignments, output_path, "Hybrid bin and rank assignments")
    
    # Stage 6: Build hybrid Markov profile
    logger.info("Step 6: Building hybrid Markov profile")
    hybrid_profile = build_markov_profile(hybrid_assignments, config.get('profile_calculation_method', 'average'))
    
    # Save hybrid profile as separate CSV files for cross-checking
    save_hybrid_profile_data(hybrid_profile, output_path, "Hybrid Markov profile data")
    
    # Save final profile data
    profile_data = {
        'hybrid_profile': hybrid_profile,
        'config': config,
        'training_genomes': list(parsed_genomes.keys()),
        'chromosome_mappings_used': chr_map_file is not None,
        'profile_type': 'hybrid' if 'positional_profile' in hybrid_profile else 'legacy'
    }
    
    if 'positional_profile' in hybrid_profile:
        profile_data.update({
            'total_positional_bins': len(hybrid_profile['positional_profile']),
            'total_ordinal_ranks': len(hybrid_profile['ordinal_profile']),
            'total_genes': hybrid_profile['hybrid_summary']['overlap_genes']
        })
    else:
        profile_data.update({
            'total_bins': len(hybrid_profile),
            'total_genes': len(set(gene for genes_data in hybrid_profile.values() for gene in genes_data.keys()))
        })
    
    profile_file = output_path / 'hybrid_profile.json'
    with open(profile_file, 'w') as f:
        json.dump(profile_data, f, indent=2, default=str)
    
    logger.info("=" * 60)
    logger.info("HYBRID PROFILE BUILDING COMPLETE")
    logger.info("=" * 60)
    logger.info(f"Profile saved to: {profile_file}")
    logger.info(f"Training genomes: {', '.join(parsed_genomes.keys())}")
    
    if 'positional_profile' in hybrid_profile:
        summary = hybrid_profile['hybrid_summary']
        logger.info(f"Profile type: HYBRID")
        logger.info(f"Positional bins: {summary['positional_stats']['total_bins']}")
        logger.info(f"Ordinal ranks: {summary['ordinal_stats']['total_ranks']}")
        logger.info(f"Overlapping genes: {summary['overlap_genes']}")
    else:
        logger.info(f"Profile type: LEGACY")
        logger.info(f"Total bins: {profile_data['total_bins']}")
        logger.info(f"Total genes: {profile_data['total_genes']}")
    
    if chromosome_mappings:
        logger.info(f"Chromosome standardization: Applied")
        logger.info(f"Standard chromosomes: {len(chromosome_mappings['standard_chromosomes'])}")
    else:
        logger.info(f"Chromosome standardization: Not applied")
        
    logger.info(f"All stage data saved to: {output_path / 'stages'}")
    logger.info("=" * 60)
    
    return profile_data




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
        
        # For analysis, we'll use the positional profile (ordinal analysis comes later)
        if 'positional_profile' in profile:
            markov_profile = profile['positional_profile']
            logger.info("Using positional profile for current analysis pipeline")
        else:
            markov_profile = profile
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
    
    query_corrected_df = correct_chromosome_orientation(query_filtered_df)
    save_stage_data(query_corrected_df, f'3_corrected_{query_id}', output_path, f"Strand-corrected query genome {query_id}")
    
    query_grouped = group_genomes_by_chromosome({query_id: query_corrected_df})
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
    
    # For now, use standard binning (future: hybrid query binning)
    from core.profile import process_genomes_binning
    query_bin_assignments = process_genomes_binning(
        {query_id: query_chromosomes}, 
        config.get('position_bin_size_kb', 100)
    )
    
    # Continue with existing analysis pipeline...
    # (Rest of the function remains largely unchanged)
    
    logger.info("Query analysis complete - using positional profile")
    logger.info("TODO: Add ordinal movement analysis in future update")
    
    return {"status": "completed", "profile_type": profile_type}







def multi_reference_alignment_pipeline(genome_fastas, output_dir):
    """Complete multi-directional alignment pipeline"""
    
    print("=" * 80)
    print("MULTI-DIRECTIONAL CHROMOSOME ALIGNMENT PIPELINE")
    print("=" * 80)
    print(f"Input genomes: {len(genome_fastas)}")
    for fasta in genome_fastas:
        print(f"  - {get_genome_name_from_fasta(fasta)}")
    print(f"Output directory: {output_dir}")
    
    # Create output directory structure
    os.makedirs(output_dir, exist_ok=True)
    stages_dir = os.path.join(output_dir, 'stages')
    os.makedirs(stages_dir, exist_ok=True)
    
    # Step 1: Run all pairwise alignments
    print("\nSTEP 1: Running pairwise RagTag alignments...")
    all_alignments = run_all_pairwise_alignments(genome_fastas, output_dir)
    
    # Save stage 1 data
    stage1_file = os.path.join(stages_dir, '1_pairwise_alignments.json')
    with open(stage1_file, 'w') as f:
        # Serialize alignment data (excluding file paths for cleaner JSON)
        serializable_alignments = {}
        for key, data in all_alignments.items():
            serializable_alignments[key] = {
                'query_genome': data['query_genome'],
                'reference_genome': data['reference_genome'],
                'mappings': data['mappings'],
                'agp_file': data['agp_file']
            }
        json.dump(serializable_alignments, f, indent=2)
    print(f"Saved pairwise alignments: {stage1_file}")
    
    # Step 2: Create unified mappings
    print("\nSTEP 2: Creating unified multi-reference mappings...")
    unified_mappings = create_unified_mappings_multi_reference(all_alignments, genome_fastas)
    
    # Step 3: Save results
    print("\nSTEP 3: Saving mappings...")
    mappings_file = save_unified_mappings(unified_mappings, output_dir)
    
    # Save stage 2 data  
    stage2_file = os.path.join(stages_dir, '2_unified_mappings.json')
    with open(stage2_file, 'w') as f:
        json.dump({
            'genome_mappings': dict(unified_mappings['genome_mappings']),
            'standard_chromosomes': sorted(list(unified_mappings['standard_chromosomes'])),
            'homology_levels': unified_mappings['homology_levels'],
            'statistics': unified_mappings['statistics']
        }, f, indent=2)
    print(f"Saved unified mappings: {stage2_file}")
    
    # Save summary
    summary = {
        'pipeline_type': 'multi_directional',
        'input_genomes': [get_genome_name_from_fasta(f) for f in genome_fastas],
        'total_genomes': len(genome_fastas),
        'pairwise_alignments': len(all_alignments),
        'standard_chromosomes': len(unified_mappings['standard_chromosomes']),
        'homology_levels': {
            level: len(chrs) for level, chrs in unified_mappings['homology_levels'].items()
        },
        'mappings_file': mappings_file,
        'output_directory': output_dir
    }
    
    summary_file = os.path.join(output_dir, 'alignment_summary.json')
    with open(summary_file, 'w') as f:
        json.dump(summary, f, indent=2)
    
    print("\n" + "=" * 80)
    print("MULTI-DIRECTIONAL ALIGNMENT PIPELINE COMPLETE")
    print("=" * 80)
    print(f"Total genomes: {len(genome_fastas)}")
    print(f"Pairwise alignments: {len(all_alignments)}")
    print(f"Standard chromosomes: {len(unified_mappings['standard_chromosomes'])}")
    print("\nHierarchy breakdown:")
    for level, chrs in unified_mappings['homology_levels'].items():
        print(f"  {level}: {len(chrs)} chromosomes")
    print(f"\nMappings file: {mappings_file}")
    print(f"Summary file: {summary_file}")
    print("=" * 80)
    
    return mappings_file














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
  
   args = parser.parse_args()
  
   if args.command == 'align-chr':
       try:
           mappings_file = multi_reference_alignment_pipeline(
               args.genome_fastas,
               args.output
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
  
   else:
       parser.print_help()
       return 1




if __name__ == "__main__":
   sys.exit(main())
