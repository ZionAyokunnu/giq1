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

# Import GIQ1 modules
from giq1.config.settings import CONFIG
from giq1.utils.file_utils import create_output_directory
from giq1.core.busco_processor import parse_busco_table, filter_busco_genes, detect_inversions
from giq1.core.quality_assessment import assess_assembly_quality
from giq1.contextual.metrics import (
    compute_inversion_rate_per_mb_busco,
    _analyze_gene_density_correlation
)

# Import working visualisation functions
from giq1.visualisation import (
    create_busco_dotplot,
    # create_ortholog_quality_plots,
    # create_inversion_landscape_plot,
    # create_statistics_summary_plot,
    # create_quality_summary_plot,
    # create_chromosome_mapping_overview,
    # create_synteny_summary_plot,
    # create_inversion_summary_plot,
    # create_synteny_block_plots,
    # create_synteny_plots,
    # create_fallback_synteny_plots,
    # create_annotated_phylogeny,
    # create_simple_tree,
    # create_tree_plot,
    # create_matplotlib_tree_plot,
    # create_tree_heatmap,
    # create_busco_phylogenetic_tree,
    # create_annotated_tree_plot,
    # create_circular_synteny_plot,
    # create_chromosome_comparison_plot,
)

logging.basicConfig(
    level=logging.INFO, 
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


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
  
        inversion_df, joined_df = detect_inversions(filtered_df1, filtered_df2, config)
        
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
            quality_metrics, output_dir, config
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
        


        # try:
        #     create_ortholog_quality_plots(joined_df, plots_dir, config)
        #     visualisation_results['quality_plots'] = plots_dir / 'ortholog_quality_assessment.png'
        # except Exception as e:
        #     logger.warning(f"create_ortholog_quality_plots failed: {e}")
        
        # # 3. Inversion landscape - uses inversion_df, joined_df (ortholog_df), plots_dir, config
        # logger.info("Creating inversion landscape plot...")
        # try:
        #     create_inversion_landscape_plot(inversion_df, joined_df, plots_dir, config)
        #     visualisation_results['landscape'] = plots_dir / 'inversion_landscape.png'
        # except Exception as e:
        #     logger.warning(f"create_inversion_landscape_plot failed: {e}")
        
        # # 4. Synteny block plots - uses joined_df (ortholog_df), plots_dir, config
        # logger.info("Creating synteny block plots...")
        # try:
        #     create_synteny_block_plots(joined_df, plots_dir, config)
        #     visualisation_results['synteny_blocks'] = plots_dir / 'synteny_block_analysis.png'
        # except Exception as e:
        #     logger.warning(f"create_synteny_block_plots failed: {e}")
        
        # # 5. Create summary plots with results dictionary
        # logger.info("Creating summary plots...")
        # try:
        #     results_dict = {
        #         'ortholog_df': joined_df,  # For functions that need the actual dataframe
        #         'synteny_df': pd.DataFrame(),  # Empty since you don't have synteny blocks
        #         'inversion_df': inversion_df,
        #         'total_genes': len(joined_df),
        #         'flipped_genes': joined_df['is_flipped'].sum(),
        #         'chromosome_pairs': len(inversion_df),
        #         'pairs_with_inversions': len(inversion_df[inversion_df['flipped_genes'] > 0]),
        #         'flip_rate': joined_df['is_flipped'].mean(),
        #         'species1_name': species1_name,
        #         'species2_name': species2_name
        #     }
            
        #     fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        #     create_statistics_summary_plot(results_dict, axes[0, 0])
        #     create_quality_summary_plot(results_dict, axes[0, 1])
        #     create_chromosome_mapping_overview(results_dict, axes[1, 0])
        #     create_synteny_summary_plot(results_dict, axes[1, 1])
            
        #     plt.tight_layout()
        #     summary_plot_path = plots_dir / 'analysis_summary.png'
        #     plt.savefig(summary_plot_path, dpi=300, bbox_inches='tight')
        #     plt.close()
        #     visualisation_results['summary'] = summary_plot_path
        # except Exception as e:
        #     logger.warning(f"create_summary_plots failed: {e}")
        
        # # 6. Inversion summary plot - uses results_dict, ax (but we'll create our own plot)
        # logger.info("Creating inversion summary plot...")
        # try:
        #     fig, ax = plt.subplots(figsize=(8, 6))
        #     create_inversion_summary_plot(results_dict, ax)
        #     plt.tight_layout()
        #     inv_summary_path = plots_dir / 'inversion_summary.png'
        #     plt.savefig(inv_summary_path, dpi=300, bbox_inches='tight')
        #     plt.close()
        #     visualisation_results['inversion_summary'] = inv_summary_path
        # except Exception as e:
        #     logger.warning(f"create_inversion_summary_plot failed: {e}")
        
        # # 7. Synteny plots - uses all_results dict, output_dir, config
        # logger.info("Creating synteny plots...")
        # try:
        #     all_results = {
        #         f"{species1_name}_vs_{species2_name}": {
        #             'full_results': True,
        #             'species_pair': (species1_name, species2_name),
        #             'ortholog_df': joined_df,
        #             'inversion_df': inversion_df
        #         }
        #     }
        #     synteny_result = create_synteny_plots(all_results, plots_dir, config)
        #     visualisation_results.update(synteny_result)
        # except Exception as e:
        #     logger.warning(f"create_synteny_plots failed: {e}")
        
        # # 8. Fallback synteny plots - uses all_results dict, output_dir, config
        # logger.info("Creating fallback synteny plots...")
        # try:
        #     fallback_result = create_fallback_synteny_plots(all_results, plots_dir, config)
        #     visualisation_results.update(fallback_result)
        # except Exception as e:
        #     logger.warning(f"create_fallback_synteny_plots failed: {e}")
        
        # # 9. Simple tree - uses species_list, config
        # logger.info("Creating simple tree...")
        # try:
        #     species_list = [species1_name, species2_name]
        #     tree_result = create_simple_tree(species_list, config)
        #     if tree_result:
        #         visualisation_results['simple_tree'] = tree_result
        # except Exception as e:
        #     logger.warning(f"create_simple_tree failed: {e}")
        
        # # 10. Circular synteny plot - FIXED: Create dummy class instance
        # logger.info("Creating circular synteny plot...")
        # try:
        #     class DummyVisualizer:
        #         def __init__(self, output_dir):
        #             self.output_dir = output_dir
            
        #     viz_instance = DummyVisualizer(plots_dir)
        #     circular_result = create_circular_synteny_plot(
        #         viz_instance,
        #         joined_df,  # ortholog_df
        #         inversion_df,  # inversion_df  
        #         species1_name,
        #         species2_name
        #     )
        #     visualisation_results['circular_synteny'] = circular_result
        # except Exception as e:
        #     logger.warning(f"create_circular_synteny_plot failed: {e}")
        
        # # 11. Chromosome comparison plot - FIXED: Create dummy class instance
        # logger.info("Creating chromosome comparison plot...")
        # try:
        #     viz_instance = DummyVisualizer(plots_dir)
        #     comparison_result = create_chromosome_comparison_plot(
        #         viz_instance,
        #         joined_df,  # ortholog_df
        #         inversion_df,  # inversion_df
        #         species1_name,
        #         species2_name
        #     )
        #     visualisation_results['chromosome_comparison'] = comparison_result
        # except Exception as e:
        #     logger.warning(f"create_chromosome_comparison_plot failed: {e}")
        
        # # ===== FUNCTIONS THAT NEED ADDITIONAL DATA (PLACEHOLDERS) =====
        
        # # 12. Annotated phylogeny - NEEDS: all_results + species_stats
        # logger.info("Creating annotated phylogeny... [PLACEHOLDER - NEEDS SPECIES_STATS]")
        # try:
        #     # PLACEHOLDER: You need species statistics
        #     species_stats = {
        #         species1_name: {
        #             'quality': {'quality_score': 0.8, 'metrics': {'total_length': 100000000}},
        #             'genome_size': 100000000
        #         },
        #         species2_name: {
        #             'quality': {'quality_score': 0.8, 'metrics': {'total_length': 100000000}},
        #             'genome_size': 100000000
        #         }
        #     }
            
        #     annotated_result = create_annotated_phylogeny(all_results, species_stats, plots_dir, config)
        #     visualisation_results.update(annotated_result)
        # except Exception as e:
        #     logger.warning(f"create_annotated_phylogeny failed: {e}")
        
        # # 13. Tree plot - NEEDS: tree object + inversion_annotations
        # logger.info("Creating tree plot... [PLACEHOLDER - NEEDS TREE OBJECT]")
        # try:
        #     # PLACEHOLDER: You need a tree object
        #     tree = None  # TODO: Provide tree object (from create_simple_tree or phylogenetic analysis)
        #     inversion_annotations = {
        #         species1_name: {'total_inversions': inversion_df['flipped_genes'].sum(), 'rate_per_mb': 0.1},
        #         species2_name: {'total_inversions': inversion_df['flipped_genes'].sum(), 'rate_per_mb': 0.1}
        #     }
        #     plot_file = plots_dir / 'tree_plot.png'
            
        #     if tree:
        #         create_tree_plot(tree, plot_file, inversion_annotations, config)
        #         visualisation_results['tree_plot'] = plot_file
        #     else:
        #         logger.warning("Skipping tree_plot - no tree object")
        # except Exception as e:
        #     logger.warning(f"create_tree_plot failed: {e}")
        
        # # 14. Matplotlib tree plot - NEEDS: inversion_data dict
        # logger.info("Creating matplotlib tree plot... [USING AVAILABLE DATA]")
        # try:
        #     inversion_data_for_tree = {
        #         species1_name: {
        #             'rate_per_mb': inversion_df['flipped_genes'].sum() / 100,  # Estimate
        #             'total_inversions': inversion_df['flipped_genes'].sum(),
        #             'genome_size': 100000000,
        #             'normalized_score': 0.5
        #         },
        #         species2_name: {
        #             'rate_per_mb': inversion_df['flipped_genes'].sum() / 100,  # Estimate
        #             'total_inversions': inversion_df['flipped_genes'].sum(),
        #             'genome_size': 100000000,
        #             'normalized_score': 0.5
        #         }
        #     }
        #     output_file = plots_dir / 'matplotlib_tree.png'
            
        #     create_matplotlib_tree_plot(inversion_data_for_tree, output_file, config)
        #     visualisation_results['matplotlib_tree'] = output_file
        # except Exception as e:
        #     logger.warning(f"create_matplotlib_tree_plot failed: {e}")
        
        # # 15. Tree heatmap - NEEDS: inversion_annotations dict
        # logger.info("Creating tree heatmap... [USING AVAILABLE DATA]")
        # try:
        #     inversion_annotations = {
        #         species1_name: {
        #             'total_inversions': inversion_df['flipped_genes'].sum(),
        #             'rate_per_mb': inversion_df['flipped_genes'].sum() / 100,
        #             'normalized_score': 0.5,
        #             'genome_size': 100000000
        #         },
        #         species2_name: {
        #             'total_inversions': inversion_df['flipped_genes'].sum(),
        #             'rate_per_mb': inversion_df['flipped_genes'].sum() / 100,
        #             'normalized_score': 0.5,
        #             'genome_size': 100000000
        #         }
        #     }
        #     heatmap_file = plots_dir / 'tree_heatmap.png'
            
        #     create_tree_heatmap(inversion_annotations, heatmap_file, config)
        #     visualisation_results['tree_heatmap'] = heatmap_file
        # except Exception as e:
        #     logger.warning(f"create_tree_heatmap failed: {e}")
        
        # # 16. BUSCO phylogenetic tree - NEEDS: all_results + species_stats
        # logger.info("Creating BUSCO phylogenetic tree... [USING AVAILABLE DATA]")
        # try:
        #     busco_tree_result = create_busco_phylogenetic_tree(all_results, species_stats, plots_dir, config)
        #     visualisation_results.update(busco_tree_result)
        # except Exception as e:
        #     logger.warning(f"create_busco_phylogenetic_tree failed: {e}")
        
        # # 17. Annotated tree plot - NEEDS: tree + node_metrics
        # logger.info("Creating annotated tree plot... [PLACEHOLDER - NEEDS TREE]")
        # try:
        #     # Use simple tree if available
        #     tree = create_simple_tree([species1_name, species2_name], config)
        #     node_metrics = {
        #         species1_name: {
        #             'inversion_count': inversion_df['flipped_genes'].sum(),
        #             'inversion_rate_per_mb': inversion_df['flipped_genes'].sum() / 100
        #         },
        #         species2_name: {
        #             'inversion_count': inversion_df['flipped_genes'].sum(),
        #             'inversion_rate_per_mb': inversion_df['flipped_genes'].sum() / 100
        #         }
        #     }
        #     plot_file = plots_dir / 'annotated_tree_plot.png'
            
        #     if tree:
        #         create_annotated_tree_plot(tree, plot_file, node_metrics, config)
        #         visualisation_results['annotated_tree_plot'] = plot_file
        #     else:
        #         logger.warning("Skipping annotated_tree_plot - no tree")
        # except Exception as e:
        #     logger.warning(f"create_annotated_tree_plot failed: {e}")
        
        # logger.info(f"Created {len([v for v in visualisation_results.values() if v is not None])} visualisation files")
        
    except Exception as e:
        logger.error(f"visualisation creation failed: {e}")
    
    return visualisation_results

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
    
    logger.info(f"Results saved to {data_dir}")


def generate_analysis_report(inversion_df: pd.DataFrame,
                           joined_df: pd.DataFrame,
                           contextual_metrics: Dict,
                           quality_metrics: Dict,
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
    'first_busco_path': '/Users/zionayokunnu/Documents/Bibionidae/busco-tables/Dioctria_linearis.tsv',
    'second_busco_path': '/Users/zionayokunnu/Documents/Bibionidae/busco-tables/Dioctria_rufipes.tsv',
    'first_species_name': 'Dioctria_linearis',
    'second_species_name': 'Dioctria_rufipes',
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
    
