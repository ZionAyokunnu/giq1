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
from core import parse_busco_table, filter_busco_genes, detect_flips, complete_chromosome_analysis, correct_chromosome_orientation
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
        

        syntenic_df1, syntenic_df2, chromosome_pairs = complete_chromosome_analysis(
            filtered_df1, filtered_df2, config
        )
        

        inversion_df, joined_df = detect_flips(syntenic_df1, syntenic_df2, config)
        
        corrected_df1, corrected_df2 = correct_chromosome_orientation(
            syntenic_df1, syntenic_df2, joined_df, config
        )
        
        final_inversion_df, final_joined_df = detect_flips(corrected_df1, corrected_df2, config)
        

        contextual_metrics = compute_contextual_metrics(
            final_inversion_df, final_joined_df, config
        )
        
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
            'chromosome_pairs': chromosome_pairs,
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
  all_inversion_events = config.get('all_inversion_events', [])
  
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
      
      f.write("## Inversion Event Analysis\n\n")
      single_gene_events = [e for e in all_inversion_events if e['type'] == 'single-gene']
      multi_gene_events = [e for e in all_inversion_events if e['type'] == 'multi-gene']
      largest_event = max([e['gene_count'] for e in all_inversion_events]) if all_inversion_events else 0
      
      f.write(f"- **Total inversion events detected:** {len(all_inversion_events)}\n")
      f.write(f"- **Single-gene inversions:** {len(single_gene_events)}\n")
      f.write(f"- **Multi-gene inversions:** {len(multi_gene_events)}\n")
      f.write(f"- **Largest inversion:** {largest_event} genes\n\n")
      
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
      f.write("| Chr1 | Chr2 | Total Genes | Flipped | Rate | Single | Multi | Largest | Type |\n")
      f.write("|------|------|-------------|---------|------|--------|-------|---------|------|\n")
      
      for _, row in inversion_df.iterrows():
          f.write(f"| {row['chr1']} | {row['chr2']} | {row['total_genes']} | "
                 f"{row['flipped_genes']} | {row['flip_rate']:.2f} | "
                 f"{row['single_gene_inversions']} | {row['multi_gene_inversions']} | "
                 f"{row['largest_inversion']} | {row['inversion_type']} |\n")
      
      if len(all_inversion_events) > 0:
          f.write("\n## Detailed Inversion Events\n\n")
          f.write(f"| Chr Pair | Type | Genes | {species1_name} Coordinates | {species2_name} Coordinates | Span (bp) |\n")
          f.write("|----------|------|-------|----------------------|----------------------|----------|\n")
          
          for event in sorted(all_inversion_events, key=lambda x: x['gene_count'], reverse=True)[:50]:
              coord1 = f"{event['start_pos1']:,} - {event['end_pos1']:,}"
              coord2 = f"{event['start_pos2']:,} - {event['end_pos2']:,}"
              f.write(f"| {event['chr1']}-{event['chr2']} | {event['type']} | {event['gene_count']} | "
                     f"{coord1} | {coord2} | {event['span_bp']:,} |\n")
          
          if len(all_inversion_events) > 50:
              f.write(f"\n*Showing top 50 events. Total: {len(all_inversion_events)} events.*\n")
      
      f.write(f"\n---\n*Report generated by GIQ1 v1.0.0*\n")
  
  logger.info(f"Analysis report saved to {report_path}")
  
  return report_path



def main():
    
    CONFIG.update({
    'first_busco_path': '/Users/zionayokunnu/Documents/Bibionidae/busco-tables/Dioctria_linearis.tsv',
    'second_busco_path': '/Users/zionayokunnu/Documents/Bibionidae/busco-tables/Dioctria_rufipes.tsv',
    'first_species_name': 'Dioctria_linearis',
    'second_species_name': 'Dioctria_rufipes',
    'first_fasta_path' : '/Users/zionayokunnu/Documents/Bibionidae/fasta/Dioctria_linearis.fna',
    'second_fasta_path' : '/Users/zionayokunnu/Documents/Bibionidae/fasta/Dioctria_rufipes.fna'
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