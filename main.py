#!/usr/bin/env python3
"""
GIQ1 - Genome Inversion Quantifier
Main entry point.

Author: Zion Ayokunnu
Supervisors: Kamil Jaron, Sasha, Arif
Version: 1.0.0
"""

import sys
import logging
from pathlib import Path
import pandas as pd
import numpy as np
from typing import Dict, Optional
from config import CONFIG
from utils import create_output_directory
from core import parse_busco_table, filter_busco_genes, detect_inversions, run_pairwise_ragtag
from contextual.metrics import (
    compute_inversion_rate_per_mb_busco,
    _analyze_gene_density_correlation
)

from visualisation import (
    create_busco_dotplot,
    create_linearised_dotplot
)

logging.basicConfig(
    level=logging.INFO, 
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def run_busco_inversion_analysis(config: Dict = None) -> Dict:
    if config is None:
        config = CONFIG
    
    output_dir = create_output_directory(config)
    logger.info(f"Output directory: {output_dir}")
    
    species1_name = config.get('first_species_name', 'Species1')
    species2_name = config.get('second_species_name', 'Species2')
    
    try:

        busco_df1 = parse_busco_table(config['first_busco_path'], config)
        busco_df2 = parse_busco_table(config['second_busco_path'], config)
    
        
        filtered_df1 = filter_busco_genes(busco_df1, config)
        filtered_df2 = filter_busco_genes(busco_df2, config)
        
        
        chr_names_1 = set(filtered_df1['sequence'])
        chr_names_2 = set(filtered_df2['sequence'])
        common_chrs = chr_names_1 & chr_names_2
  
  
        logger.info(f"Genome 1 chromosomes: {len(chr_names_1)} - {sorted(list(chr_names_1))[:5]}")
        logger.info(f"Genome 2 chromosomes: {len(chr_names_2)} - {sorted(list(chr_names_2))[:5]}")
        logger.info(f"Common chromosomes: {len(common_chrs)} - {sorted(list(common_chrs))[:5]}")
        
        compatibility_threshold = 0.7 * min(len(chr_names_1), len(chr_names_2))
        
        
        logger.info(f"Compatibility threshold: {compatibility_threshold}")
        logger.info(f"Common chromosomes found: {len(common_chrs)}")
        
        
        if len(common_chrs) < compatibility_threshold:
            logger.info("Chromosome names incompatible - running RagTag alignment...")
            
            if 'first_fasta_path' not in config or 'second_fasta_path' not in config:
                logger.warning("FASTA file paths not provided - skipping RagTag alignment")
                chromosome_pairs = [(chr, chr) for chr in common_chrs]  # Fix: ensure it's always set
            else:
                chromosome_mapping = run_pairwise_ragtag(config)
                chromosome_pairs = list(chromosome_mapping.items())
                
                logger.info(f"Sample RagTag mappings:")
                for i, (query_chr, ref_chr) in enumerate(chromosome_pairs[:5]):
                    logger.info(f"  {query_chr} → {ref_chr}")

                # DEBUG: Check if mapped chromosomes exist in BUSCO data
                mapped_query_chrs = set([pair[0] for pair in chromosome_pairs])
                mapped_ref_chrs = set([pair[1] for pair in chromosome_pairs])

                logger.info(f"BUSCO genome1 chromosomes: {sorted(list(chr_names_1))[:5]}")
                logger.info(f"BUSCO genome2 chromosomes: {sorted(list(chr_names_2))[:5]}")
                logger.info(f"RagTag mapped query chrs: {sorted(list(mapped_query_chrs))[:5]}")
                logger.info(f"RagTag mapped ref chrs: {sorted(list(mapped_ref_chrs))[:5]}")

                overlap1 = chr_names_1 & mapped_ref_chrs
                overlap2 = chr_names_2 & mapped_query_chrs
                logger.info(f"Overlap genome1 BUSCO ∩ RagTag ref: {len(overlap1)}")
                logger.info(f"Overlap genome2 BUSCO ∩ RagTag query: {len(overlap2)}")

                # Convert mapping to pairs (quick fix)
                chromosome_pairs = list(chromosome_mapping.items())  # [(genome2_chr, genome1_chr), ...]
                
                logger.info(f"After RagTag mapping: {len(chromosome_pairs)} chromosome pairs")
                
                valid_pairs = []
                for query_chr, ref_chr in chromosome_mapping.items():
                    # Only keep pairs where both chromosomes exist in BUSCO data
                    if query_chr in chr_names_2 and ref_chr in chr_names_1:
                        valid_pairs.append((ref_chr, query_chr))  # (genome1_chr, genome2_chr)

                chromosome_pairs = valid_pairs  # ← Use filtered pairs, not all 491!

                logger.info(f"Valid chromosome pairs after filtering: {len(chromosome_pairs)}")
                for pair in chromosome_pairs[:5]:
                    logger.info(f"  {pair[0]} ↔ {pair[1]}")
        else:
            chromosome_pairs = [(chr, chr) for chr in common_chrs]  # Fix: same format
            logger.info("Using direct chromosome name matching - chromosomes are compatible")
        
        logger.info("Applying per-chromosome synteny filtering...")
        logger.info(f"Starting with {len(filtered_df1)} genes in genome1, {len(filtered_df2)} genes in genome2")
        
        syntenic_df1_parts = []
        syntenic_df2_parts = []
        
        total_syntenic_genes = 0
        
        for chr1_name, chr2_name in chromosome_pairs:
            chr1_genes = filtered_df1[filtered_df1['sequence'] == chr1_name]
            chr2_genes = filtered_df2[filtered_df2['sequence'] == chr2_name]
            
            # Find syntenic genes (present on both corresponding chromosomes)
            common_genes = set(chr1_genes['busco_id']) & set(chr2_genes['busco_id'])
            
            if len(common_genes) > 0:
                logger.info(f"  {chr1_name} ↔ {chr2_name}: {len(chr1_genes)} vs {len(chr2_genes)} → {len(common_genes)} syntenic genes")
                syntenic_chr1 = chr1_genes[chr1_genes['busco_id'].isin(common_genes)]
                syntenic_chr2 = chr2_genes[chr2_genes['busco_id'].isin(common_genes)]
                
                syntenic_df1_parts.append(syntenic_chr1)
                syntenic_df2_parts.append(syntenic_chr2)
                
                total_syntenic_genes += len(common_genes)
        
        if not syntenic_df1_parts:
            raise ValueError("No syntenic genes found between genomes")
        
        syntenic_df1 = pd.concat(syntenic_df1_parts, ignore_index=True)
        syntenic_df2 = pd.concat(syntenic_df2_parts, ignore_index=True)
        
        
        logger.info(f"Synteny filtering results:")
        logger.info(f"  Before: {len(filtered_df1)} + {len(filtered_df2)} = {len(filtered_df1) + len(filtered_df2)} total genes")
        logger.info(f"  After: {len(syntenic_df1)} + {len(syntenic_df2)} = {len(syntenic_df1) + len(syntenic_df2)} total genes")
        logger.info(f"  Filtered out: {(len(filtered_df1) + len(filtered_df2)) - (len(syntenic_df1) + len(syntenic_df2))} genes")
        
        inversion_df, joined_df = detect_inversions(syntenic_df1, syntenic_df2, config)

        
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
    
    gene_density_corr = _analyze_gene_density_correlation(inversion_df, joined_df, config)
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
            logger.info(f"Data going to visualizations:")
            logger.info(f"  joined_df size: {len(joined_df)}")
            logger.info(f"  Unique chromosomes chr1: {joined_df['chr1'].nunique()}")
            logger.info(f"  Unique chromosomes chr2: {joined_df['chr2'].nunique()}")
            logger.info(f"  Sample chromosome pairs: {list(zip(joined_df['chr1'].unique()[:3], joined_df['chr2'].unique()[:3]))}")
            logger.info(f"  All BUSCO genes in complete status: {joined_df.shape}")
            fig, ax = create_linearised_dotplot(
                joined_df=joined_df, 
                plots_dir=plots_dir,
                config=config['dotplot_config']  # Pass the nested dotplot config
            )
            visualisation_results['linearised_dotplot'] = plots_dir / 'linearised_busco_dotplot.png'
            
        except Exception as e:
            logger.warning(f"linearised dotplot failed: {e}")
        
    except Exception as e:
        logger.error(f"visualisation creation failed: {e}")
    
    return visualisation_results

def save_results(inversion_df: pd.DataFrame,
                joined_df: pd.DataFrame,
                contextual_metrics: Dict,
                output_dir: Path,
                config: Dict):
    
    data_dir = output_dir / 'data'
    data_dir.mkdir(exist_ok=True)
    
    inversion_df.to_csv(data_dir / 'inversion_analysis.csv', index=False)
    joined_df.to_csv(data_dir / 'joined_gene_data.csv', index=False)
    
    import json
    with open(data_dir / 'contextual_metrics.json', 'w') as f:
        json.dump(contextual_metrics, f, indent=2, default=str)
    
    logger.info(f"Results saved to {data_dir}")


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
    """Main entry point when run as script"""
    
    CONFIG.update({
    'first_busco_path': '/Users/zionayokunnu/Documents/Bibionidae/busco-tables/Bibio_marci.tsv',
    'second_busco_path': '/Users/zionayokunnu/Documents/Bibionidae/busco-tables/Dilophus_febrilis.tsv',
    'first_species_name': 'Bibio_marci',
    'second_species_name': 'Dilophus_febrilis',
    'first_fasta_path' : '/Users/zionayokunnu/Documents/Bibionidae/fasta/Bibio_marci.fna',
    'second_fasta_path' : '/Users/zionayokunnu/Documents/Bibionidae/fasta/Dilophus_febrilis.fna'
    })
    
    import random
    random.seed(42)
    np.random.seed(42)
    
    try:
        results = run_busco_inversion_analysis()
        
        print(f"Species analyzed: {' vs '.join(results['species_names'])}")
        print(f"Total inversions detected: {results['inversion_df']['flipped_genes'].sum()}")
        print(f"Output directory: {results['output_dir']}")
        print(f"Report: {results['report_path']}")
        print("=" * 80)
        
        return 0
        
    except Exception as e:
        logger.error(f"Analysis failed: {str(e)}")
        print(f"\nERROR: {str(e)}")
        
        if CONFIG.get('enable_debug_output', False):
            import traceback
            traceback.print_exc()
        
        return 1


if __name__ == "__main__":
    sys.exit(main())
    
