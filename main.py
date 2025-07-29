#!/usr/bin/env python3
"""
 main entry point for inversion analysis
"""

import sys
import random
import logging
from pathlib import Path
import pandas as pd
import numpy as np

from giq1.config import (
    CONFIG,
)

from giq1.utils import (
    create_output_directory
)

from giq1.core import (
    parse_busco_table,
    filter_busco_genes,
    detect_inversions,
    
)

from giq1.visualisation import (

)


logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger(__name__)


def detect_inversions(config=None):
    config = CONFIG
    
    logger.info("-" * 80)
    logger.info('Starting')
    logger.info("-" * 80)
    
    output_dir = create_output_directory(config)
    
    species1_name = config.get('first_species_name', 'Species1')
    species2_name = config.get('second_species_name', 'Species2')
    
    try:
        first_busco_raw = parse_busco_table(config['first_busco_path'], config)
        second_busco_raw = parse_busco_table(config['second_busco_path'], config)
        
        registry.register_file(
            'first_busco_raw',
            first_busco_raw,
            'csv',
            f"Raw BUSCO results for {species1_name} genome",
            metadata={'genome': {species1_name}, 'source': config['first_busco_path']}
        )
        
        registry.register_file(
            'second_busco_raw',
            second_busco_raw,
            'csv',
            f"Raw BUSCO results for {species2_name} genome",
            metadata={'genome': {species2_name}, 'source': config['second_busco_path']}
        )
        
        logger.info("\n" + "=" * 80)
        logger.info("ANALYSIS")
        logger.info("=" * 80)
        
        inversion_df = detect_inversions(busco1, busco2, config):
            busco1 = ''
            busco2 = '' 
            return inversion_df
        
        save_enhanced_results(output_dir, results_dict, config)
        
        exporter.export_analysis_summary(results_dict)
        
        create_enhanced_visualisations(output_dir, results_dict, config)
        plots_dir = output_dir / 'plots'
        if plots_dir.exists():
            for plot_file in plots_dir.glob('*.png'):
                registry.register_file(
                    f'plot_{plot_file.stem}',
                    {'path': str(plot_file)},
                    'json',
                    f'Visualisation: {plot_file.stem}',
                    dependencies=['inversions', 'synteny_blocks', 'inversion_events'],
                    metadata={'plot_type': 'png', 'visualisation': True}
                )
                
        def add_contextual_metrics_to_busco_analysis(results_df: pd.DataFrame, 
                                                joined_df: pd.DataFrame,
                                                fasta_file: str = None,
                                                repeat_file: str = None) -> Dict:
            """
            Add contextual metrics to your existing BUSCO inversion analysis
            
            Args:
                results_df: Output from your detect_inversions function
                joined_df: Joined gene data from your detect_inversions function
                fasta_file: Optional path to genome FASTA for GC analysis
                repeat_file: Optional path to repeat elements file
                
            Returns:
                Dictionary with all contextual metrics
            """
            
            contextual_metrics = {}
            
            # 1. Inversion rate per Mb
            print("Computing inversion rates...")
            rate_metrics = compute_inversion_rate_per_mb_busco(results_df, joined_df)
            contextual_metrics['inversion_rates'] = rate_metrics
            
            # 2. Repeat correlation (if repeat file provided)
            if repeat_file:
                print("Analyzing repeat correlation...")
                repeat_correlation = analyze_repeat_correlation_busco(results_df, joined_df, repeat_file)
                contextual_metrics['repeat_correlation'] = repeat_correlation
            else:
                print("Analyzing gene density correlation (repeat proxy)...")
                gene_density_corr = _analyze_gene_density_correlation(results_df, joined_df)
                contextual_metrics['gene_density_correlation'] = gene_density_corr
            
            # 3. GC correlation (if FASTA file provided)
            if fasta_file:
                print("Analyzing GC content correlation...")
                gc_correlation = analyze_gc_correlation_busco(results_df, joined_df, fasta_file)
                contextual_metrics['gc_correlation'] = gc_correlation
            
            return contextual_metrics
        


        logger.info("\n" + ">" * 80)
        logger.info("ANALYSIS SUCCESSFUL,THANK YOU FOR USING THE GIA1 SYSTEM!")
        logger.info("<" * 80)
        

        return {
            'ortholog_df': ortholog_df,
            'paralog_df': paralog_df,
            'inversion_df': inversion_df,
            'output_dir': output_dir,
            'config': config,
            'registry': registry,
            'file_manifest': str(manifest_path),
            'integrity_report': integrity_report
        }
        
    except Exception as e:
        logger.error(f"Analysis failed: {str(e)}")
        if config.get('enable_debug_output', False):
            import traceback
            traceback.print_exc()
        raise


def save_enhanced_results(output_dir, results, config):
    """Legacy data saving"""
    data_dir = output_dir / 'data'
    
    results['inversion_df'].to_csv(data_dir / Path(config['inversion_summary_csv']).name, index=False)



def generate_comprehensive_report(output_dir, results, registry=None):
    """Analyses report summary and registry"""
    reports_dir = output_dir / 'reports'
    reports_dir.mkdir(exist_ok=True)
    
    if registry:
        report_path = reports_dir / 'analysis_report_with_registry.md'
        
        with open(report_path, 'w') as f:
            f.write("# Genome Inversion Analysis Report\n\n")
            f.write(f"Generated: {pd.Timestamp.now()}\n\n")
            f.write(f"This is between the genomes of {results['config']['first_species_name']} and {results['config']['second_species_name']}.\n\n")

            f.write(f"- Total registered files: {len(registry.list_files())}\n")
            f.write(f"- Registry location: {registry.registry_file}\n")
            
            for file_type, count in sorted(file_types.items()):
                f.write(f"- {file_type}: {count} files\n")
            
            f.write("Analysis Summary\n\n")
            f.write(f"{len(df1)} genes in {species1_name}, and {len(df2)} in {species2_name}")
            f.write(f"{len(joined_genes)} common genes with {total_flipped} flipped whole genes")
            f.write(f"{len(chromosome_dict)} chromosome pairs across {species1_name} and {species2_name}")
            f.write(f"   {chr1} vs {chr2}: {total_genes} genes, {flipped_genes} flipped (rate: {flip_rate:.2f})")
            f.write(f"   Strand consistency is {strand_consistency:.2f} across the {species1_name} and {species2_name} genome")
            f.write(f"   Correlation between start coords per chromosomes of the two genomes: {correlation:.3f} (p={p_value:.3f})")
            f.write(f"   - Total chromosome pairs: {len(results_df)}")
            f.write(f"   - Pairs with inversions: {len(results_df[results_df['flipped_genes'] > 0])}")
            f.write(f"   - Total flipped genes: {results_df['flipped_genes'].sum()}")
            f.write(inversion_results[['chr1', 'chr2', 'total_genes', 'flipped_genes', 
                           'correlation', 'strand_consistency', 'inversion_type']])
            f.write(f"Genome size per mb (using {species1_name} as reference) is {genome_size_mb}")
            f.write(f"Inversion (chromosome pair) per mb = {inversion_per_mb}")
            f.write(f"Inversion (each flipped gene) per mb = {flipped_gene_per_mb}")
            f.write(f"Inversion rate by type = {rate_by_type[inv_type]}")
            f.write(f"Inversion rate by chromosome = {rates_by_chr[chrom]}")
            f.write(f"{species1_name} density is {gene_density_species_a}")
            f.write(f"{species2_name} density is {gene_density_species_b}")
        

def enhanced_visualizations(output_dir, results_dict, config):
    """
    Create comprehensive visualizations including dotplots, chromosome map, summary plots
    """    
    plots_dir = os.path.join(output_dir, 'plots')
    os.makedirs(plots_dir, exist_ok=True)
    
    plt.style.use('default')
    plt.rcParams['figure.dpi'] = config.get('dpi', 300)
    plt.rcParams['font.size'] = config.get('font_size', 12)
    
    create_busco_dotplot(results_dict['?'], plots_dir, config)
    create_ortholog_quality_plots(results_dict['?'], plots_dir, config)
    create_synteny_block_plots(results_dict['?'], plots_dir, config)
    create_inversion_landscape_plot(results_dict['?'], results_dict['?'], plots_dir, config)
    create_statistics_summary_plot(results_dict, ax)
    create_quality_summary_plot(results_dict, ax):
    create_chromosome_mapping_overview(results_dict, ax)
    create_synteny_summary_plot(results_dict, ax)
    create_inversion_summary_plot(results_dict, ax)
    _create_synteny_plots(self, all_results: Dict, output_dir: Path)
    _create_fallback_synteny_plots(self, all_results: Dict, output_dir: Path)
    create_annotated_phylogeny(self, all_results: Dict, species_stats: Dict, 
                                output_dir: Path)
    _create_simple_tree(self, species_list: List[str])
    _create_tree_plot(self, tree: 'Tree', output_file: Path, inversion_data: Dict)
    _create_matplotlib_tree_plot(self, inversion_data: Dict, output_file: Path)
    _create_tree_heatmap(self, inversion_data: Dict, output_file: Path)
    create_busco_phylogenetic_tree(self, all_results: Dict, species_stats: Dict, output_dir: Path)
    create_annotated_tree_plot(self, tree, output_file, node_metrics)
    create_circular_synteny_plot(self,
                            ortholog_df: pd.DataFrame,
                            inversion_df: pd.DataFrame,
                            species1: str,
                            species2: str)
    create_chromosome_comparison_plot(self,
                                    ortholog_df: pd.DataFrame,
                                    inversion_df: pd.DataFrame,
                                    species1: str,
                                    species2: str)



if __name__ == "__main__":
    random.seed(42)
    np.random.seed(42)
    config = CONFIG
    species1_name = CONFIG['first_species_name']
    species2_name = CONFIG['second_species_name']
    
    try:
        results = detect_inversions(config)
        
        print("\n" + "=" * 80)
        print(f"Analysis between {config['first_species_name']} and {config['second_species_name']}")
        
    except Exception as e:
        logger.error(f"Analysis failed: {str(e)}")
        print(f"\nError: {str(e)}")
        if config.get('enable_debug_output', False):
            import traceback
            traceback.print_exc()
    
    print("\n" + "=" * 80)