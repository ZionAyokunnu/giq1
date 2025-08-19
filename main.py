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
from core import parse_busco_table, filter_busco_genes, detect_flips, run_pairwise_ragtag, correct_chromosome_orientation
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
   
   data_dir = output_dir / 'data'
   data_dir.mkdir(parents=True, exist_ok=True)
   
   species1_name = config.get('first_species_name', 'Species1')
   species2_name = config.get('second_species_name', 'Species2')
   
   try:
       busco_df1 = parse_busco_table(config['first_busco_path'], config)
       busco_df2 = parse_busco_table(config['second_busco_path'], config)
       
       filtered_df1 = filter_busco_genes(busco_df1, config, species_name=species1_name)
       filtered_df2 = filter_busco_genes(busco_df2, config, species_name=species2_name)
       
       inversion_df, joined_df = detect_flips(filtered_df1, filtered_df2, config)
       
       corrected_df1, corrected_df2 = correct_chromosome_orientation(filtered_df1, filtered_df2, joined_df, config)
       
       chr_names_1 = set(corrected_df1['sequence'])
       chr_names_2 = set(corrected_df2['sequence'])
       common_chrs = chr_names_1 & chr_names_2
       
       chr_analysis_data = [{
           'species': species1_name,
           'chromosome_count': len(chr_names_1),
           'chromosomes': sorted(list(chr_names_1))
       }, {
           'species': species2_name,
           'chromosome_count': len(chr_names_2),
           'chromosomes': sorted(list(chr_names_2))
       }]
       
       compatibility_threshold = 0.7 * min(len(chr_names_1), len(chr_names_2))
       
       if len(common_chrs) < compatibility_threshold:
           logger.info("Chromosome names incompatible - running RagTag alignment...")
           
           if 'first_fasta_path' not in config or 'second_fasta_path' not in config:
               chromosome_pairs = [(chr, chr) for chr in common_chrs]
           else:
               chromosome_mapping = run_pairwise_ragtag(config)
               
               valid_pairs = []
               for query_chr, ref_chr in chromosome_mapping.items():
                   if query_chr in chr_names_2 and ref_chr in chr_names_1:
                       chr1_genes = corrected_df1[corrected_df1['sequence'] == ref_chr]
                       chr2_genes = corrected_df2[corrected_df2['sequence'] == query_chr]
                       common_genes = set(chr1_genes['busco_id']) & set(chr2_genes['busco_id'])
                       
                       if len(common_genes) > 10:
                           valid_pairs.append((ref_chr, query_chr, len(common_genes)))
               
               valid_pairs_data = []
               for chr1, chr2, count in valid_pairs:
                   valid_pairs_data.append({
                       'chromosome1': chr1,
                       'chromosome2': chr2,
                       'common_gene_count': count
                   })
               
               valid_pairs.sort(key=lambda x: x[2], reverse=True) 
               
               used_chr1 = set()
               used_chr2 = set()
               final_pairs = []
               
               for chr1, chr2, count in valid_pairs:
                   if chr1 not in used_chr1 and chr2 not in used_chr2:
                       final_pairs.append((chr1, chr2))
                       used_chr1.add(chr1)
                       used_chr2.add(chr2)
                       logger.info(f"  Selected pair: {chr1} ↔ {chr2} ({count} genes)")
               
               chromosome_pairs = final_pairs
               
               final_pairs_data = []
               for chr1, chr2 in final_pairs:
                   final_pairs_data.append({
                       'chromosome1': chr1,
                       'chromosome2': chr2
                   })
               
       else:
           chromosome_pairs = [(chr, chr) for chr in common_chrs]
           logger.info("Using direct chromosome name matching - chromosomes are compatible")
           
       syntenic_df1_parts = []
       syntenic_df2_parts = []
       synteny_summary_data = []
       
       for chr1_name, chr2_name in chromosome_pairs:
           chr1_genes = corrected_df1[corrected_df1['sequence'] == chr1_name]
           chr2_genes = corrected_df2[corrected_df2['sequence'] == chr2_name]
           
           common_genes = set(chr1_genes['busco_id']) & set(chr2_genes['busco_id'])
           
           if len(common_genes) > 0:
               syntenic_chr1 = chr1_genes[chr1_genes['busco_id'].isin(common_genes)]
               syntenic_chr2 = chr2_genes[chr2_genes['busco_id'].isin(common_genes)]
               
               syntenic_df1_parts.append(syntenic_chr1)
               syntenic_df2_parts.append(syntenic_chr2)
       
       if not syntenic_df1_parts:
           raise ValueError("No syntenic genes found between genomes")
       
       syntenic_df1 = pd.concat(syntenic_df1_parts, ignore_index=True)
       syntenic_df2 = pd.concat(syntenic_df2_parts, ignore_index=True)
       
       final_inversion_df, final_joined_df = detect_flips(syntenic_df1, syntenic_df2, config)
       
       contextual_metrics = compute_contextual_metrics(
           final_inversion_df, final_joined_df, config
       )
       
       contextual_metrics_data = []
       for metric_category, metric_data in contextual_metrics.items():
           if isinstance(metric_data, dict):
               for metric_name, metric_value in metric_data.items():
                   contextual_metrics_data.append({
                       'metric_category': metric_category,
                       'metric_name': metric_name,
                       'metric_value': metric_value
                   })
           else:
               contextual_metrics_data.append({
                   'metric_category': metric_category,
                   'metric_name': 'value',
                   'metric_value': metric_data
               })
       
       visualisation_results = create_analysis_visualisations(
           final_inversion_df, final_joined_df, output_dir, config
       )
       
       save_results(
           final_inversion_df, final_joined_df, contextual_metrics, output_dir, config
       )
       
       report_path = generate_analysis_report(
           final_inversion_df, final_joined_df, contextual_metrics, 
           output_dir, config
       )
       
       results = {
           'inversion_df': final_inversion_df,
           'joined_df': final_joined_df,
           'contextual_metrics': contextual_metrics,
           'visualisation_results': visualisation_results,
           'output_dir': output_dir,
           'report_path': report_path,
           'config': config,
           'species_names': [species1_name, species2_name]
       }
       
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
    

    data_dir = output_dir / 'data'
    data_dir.mkdir(parents=True, exist_ok=True)
    
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
    

    
    import json
    with open(data_dir / 'contextual_metrics.json', 'w') as f:
        json.dump(contextual_metrics, f, indent=2, default=str)
    



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
        f.write("# GIQ1 - Genome Inversion Quantification Report\n\n")
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