import sys
import argparse
import logging
import os
from pathlib import Path
import pandas as pd
import json
from typing import Dict, List

from config.settings import CONFIG
from .busco_processor import parse_busco_table, filter_busco_genes
from .profile import process_genomes_binning_hybrid, build_markov_profile_hybrid, filter_chromosomes_for_genome, group_genomes_by_chromosome, build_markov_profile
from .chr_aligner import load_chromosome_mappings, standardize_chromosome_name_unified, run_all_pairwise_alignments, create_unified_mappings_multi_reference, save_unified_mappings, get_genome_name_from_fasta
from .outputs import save_stage_data, save_hybrid_profile_data, save_hybrid_assignments



def filter_chromosomes_by_gene_count_before_grouping(filtered_df, genome_id, min_genes=100):
    """Filter out chromosomes with too few genes before grouping"""
    
    # Count genes per chromosome
    gene_counts = filtered_df['sequence'].value_counts()
    
    print(f"  Chromosome gene counts for {genome_id}:")
    for chr_name, count in gene_counts.items():
        print(f"    {chr_name}: {count} genes")
    
    # Filter chromosomes
    valid_chromosomes = []
    for chr_name, count in gene_counts.items():
        # Filter small contigs aggressively
        if 'CATVHY010000' in chr_name:
            if count < min_genes:
                print(f"    ✗ Filtered contig {chr_name}: {count} < {min_genes} genes")
                continue
        
        # Keep chromosomes with sufficient genes
        if count >= min_genes:
            valid_chromosomes.append(chr_name)
            print(f"    ✓ Kept {chr_name}: {count} genes")
        else:
            print(f"    ✗ Filtered {chr_name}: {count} < {min_genes} genes")
    
    # Return filtered DataFrame
    return filtered_df[filtered_df['sequence'].isin(valid_chromosomes)]




logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


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
        
        # Filter chromosomes by gene count (NEW)
        filtered_df = filter_chromosomes_by_gene_count_before_grouping(
            filtered_df, genome_id, min_genes=100
        )
        
        filtered_genomes[genome_id] = filtered_df
        
        save_stage_data(
            filtered_df,
            f'2_filtered_{genome_id}',
            output_path,
            f"Filtered BUSCO genes for {genome_id}: {len(filtered_df)} complete genes"
        )
    
    
    # Stage 4: Group by chromosome
    logger.info("Step 4: Grouping by chromosomes")
    grouped_genomes = group_genomes_by_chromosome(filtered_genomes)
    
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

    # CREATE BUSCO FILE MAPPING
    busco_file_mapping = {}
    parsed_genomes = {}
    
    for busco_file in busco_files:
        busco_path = Path(busco_file)
        genome_id = busco_path.stem  # Get genome name from filename
        
        # Add to mapping
        busco_file_mapping[genome_id] = str(busco_path)  # ← CREATE MAPPING
        
        # Parse BUSCO files (existing logic)
        busco_df = parse_busco_table(str(busco_path), config)
        parsed_genomes[genome_id] = busco_df
    
    print(f"Created BUSCO file mapping for: {list(busco_file_mapping.keys())}")
    
    # Stage 5: Process hybrid binning
    logger.info("Step 5: Processing hybrid binning (positional + ordinal)")
    hybrid_assignments = process_genomes_binning_hybrid(
        grouped_genomes, 
        config.get('position_bin_size_kb', 100),
        chromosome_mappings,
        busco_file_mapping
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
            'total_genes': len(set(
                gene for genes_data in hybrid_profile['positional_profile'].values() 
                for gene in genes_data.keys() 
                if gene not in ['position_summary', 'rank_summary']
            ))
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




def multi_reference_alignment_pipeline(genome_fastas, output_dir, busco_file_mapping=None):
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
    unified_mappings = create_unified_mappings_multi_reference(all_alignments, genome_fastas, busco_file_mapping) 
    
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

