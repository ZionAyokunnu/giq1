#!/usr/bin/env python3
"""
GIQ1 - Genome Inversion Quantifier
Main entry point for BUSCO-based chromosomal inversion analysis

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

# Import working visualization functions
from giq1.visualisation import (
    create_busco_dotplot,
    create_ortholog_quality_plots,
    create_inversion_landscape_plot,
    create_statistics_summary_plot,
    create_quality_summary_plot,
    create_chromosome_mapping_overview,
    create_synteny_summary_plot,
    create_inversion_summary_plot,
    create_synteny_block_plots,
    create_synteny_plots,
    create_fallback_synteny_plots,
    create_annotated_phylogeny,
    create_simple_tree,
    create_tree_plot,
    create_matplotlib_tree_plot,
    create_tree_heatmap,
    create_busco_phylogenetic_tree,
    create_annotated_tree_plot,
    create_circular_synteny_plot,
    create_chromosome_comparison_plot,
)

# Set up logging
logging.basicConfig(
    level=logging.INFO, 
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def run_busco_inversion_analysis(config: Dict = None) -> Dict:
    """
    Main function to run complete BUSCO inversion analysis pipeline
    
    Args:
        config: Configuration dictionary (uses default CONFIG if None)
        
    Returns:
        Dictionary with analysis results and file paths
    """
    
    if config is None:
        config = CONFIG
    
    logger.info("=" * 80)
    logger.info("GIQ1 - GENOME INVERSION QUANTIFIER")
    logger.info("BUSCO-based Chromosomal Inversion Analysis")
    logger.info("=" * 80)
    
    # Create output directory
    output_dir = create_output_directory(config)
    logger.info(f"Output directory: {output_dir}")
    
    # Get species names
    species1_name = config.get('first_species_name', 'Species1')
    species2_name = config.get('second_species_name', 'Species2')
    
    logger.info(f"Analyzing inversions between {species1_name} and {species2_name}")
    
    try:
        # Step 1: Parse BUSCO tables
        logger.info("\n" + "-" * 40)
        logger.info("STEP 1: PARSING BUSCO TABLES")
        logger.info("-" * 40)
        
        busco_df1 = parse_busco_table(config['first_busco_path'], config)
        busco_df2 = parse_busco_table(config['second_busco_path'], config)
        
        logger.info(f"Parsed {len(busco_df1)} genes from {species1_name}")
        logger.info(f"Parsed {len(busco_df2)} genes from {species2_name}")
        
        # Step 2: Filter BUSCO genes
        logger.info("\n" + "-" * 40)
        logger.info("STEP 2: FILTERING BUSCO GENES")
        logger.info("-" * 40)
        
        filtered_df1 = filter_busco_genes(busco_df1, config)
        filtered_df2 = filter_busco_genes(busco_df2, config)
        
        logger.info(f"Filtered to {len(filtered_df1)} complete genes from {species1_name}")
        logger.info(f"Filtered to {len(filtered_df2)} complete genes from {species2_name}")
        
        # Step 3: Quality assessment (optional)
        quality_metrics = {}
        if config.get('assess_quality', True):
            logger.info("\n" + "-" * 40)
            logger.info("STEP 3: QUALITY ASSESSMENT")
            logger.info("-" * 40)
            
            fasta1_path = config.get('first_fasta_path')
            fasta2_path = config.get('second_fasta_path')
            
            if fasta1_path:
                quality1 = assess_assembly_quality(fasta1_path, filtered_df1, config)
                quality_metrics[species1_name] = quality1
                
            if fasta2_path:
                quality2 = assess_assembly_quality(fasta2_path, filtered_df2, config)
                quality_metrics[species2_name] = quality2
        
        # Step 4: Detect inversions
        logger.info("\n" + "-" * 40)
        logger.info("STEP 4: DETECTING INVERSIONS")
        logger.info("-" * 40)
        
        inversion_df, joined_df = detect_inversions(filtered_df1, filtered_df2, config)
        
        logger.info(f"Detected inversions across {len(inversion_df)} chromosome pairs")
        logger.info(f"Total genes compared: {len(joined_df)}")
        logger.info(f"Total flipped genes: {joined_df['is_flipped'].sum()}")
        
        # Step 5: Contextual metrics
        logger.info("\n" + "-" * 40)
        logger.info("STEP 5: CONTEXTUAL METRICS")
        logger.info("-" * 40)
        
        contextual_metrics = compute_contextual_metrics(
            inversion_df, joined_df, config
        )
        
        # Step 6: Create visualizations
        logger.info("\n" + "-" * 40)
        logger.info("STEP 6: CREATING VISUALIZATIONS")
        logger.info("-" * 40)
        
        visualization_results = create_analysis_visualizations(
            inversion_df, joined_df, output_dir, config
        )
        
        # Step 7: Save results
        logger.info("\n" + "-" * 40)
        logger.info("STEP 7: SAVING RESULTS")
        logger.info("-" * 40)
        
        save_results(
            inversion_df, joined_df, contextual_metrics, 
            quality_metrics, output_dir, config
        )
        
        # Step 8: Generate summary report
        logger.info("\n" + "-" * 40)
        logger.info("STEP 8: GENERATING REPORT")
        logger.info("-" * 40)
        
        report_path = generate_analysis_report(
            inversion_df, joined_df, contextual_metrics, 
            quality_metrics, output_dir, config
        )
        
        # Compile final results
        results = {
            'inversion_df': inversion_df,
            'joined_df': joined_df,
            'contextual_metrics': contextual_metrics,
            'quality_metrics': quality_metrics,
            'visualization_results': visualization_results,
            'output_dir': output_dir,
            'report_path': report_path,
            'config': config,
            'species_names': [species1_name, species2_name]
        }
        
        logger.info("\n" + "=" * 80)
        logger.info("ANALYSIS COMPLETED SUCCESSFULLY!")
        logger.info(f"Results saved to: {output_dir}")
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
    """Compute contextual metrics for the analysis"""
    
    contextual_metrics = {}
    
    # 1. Inversion rate per Mb
    logger.info("Computing inversion rates...")
    rate_metrics = compute_inversion_rate_per_mb_busco(inversion_df, joined_df)
    contextual_metrics['inversion_rates'] = rate_metrics
    
    # 2. Gene density correlation
    logger.info("Analyzing gene density correlations...")
    gene_density_corr = _analyze_gene_density_correlation(inversion_df, joined_df, config)
    contextual_metrics['gene_density_correlation'] = gene_density_corr
    
    # 3. Additional metrics can be added here
    # TODO: Add repeat correlation and GC correlation when files are available
    
    return contextual_metrics


def create_analysis_visualizations(inversion_df: pd.DataFrame,
                                 joined_df: pd.DataFrame,
                                 output_dir: Path,
                                 config: Dict) -> Dict:
    """Create all analysis visualizations"""
    
    plots_dir = output_dir / 'plots'
    plots_dir.mkdir(exist_ok=True)
    
    visualization_results = {}
    
    # Configure matplotlib
    import matplotlib.pyplot as plt
    plt.style.use('default')
    plt.rcParams['figure.dpi'] = config.get('dpi', 300)
    plt.rcParams['font.size'] = config.get('font_size', 12)
    
    try:
        # 1. BUSCO dotplot
        logger.info("Creating BUSCO synteny dotplot...")
        create_busco_dotplot(joined_df, plots_dir, config)
        visualization_results['dotplot'] = plots_dir / 'busco_synteny_dotplot.png'
        
        # 2. Ortholog quality plots
        logger.info("Creating ortholog quality plots...")
        create_ortholog_quality_plots(joined_df, plots_dir, config)
        visualization_results['quality_plots'] = plots_dir / 'ortholog_quality_plots.png'
        
        # 3. Inversion landscape
        logger.info("Creating inversion landscape plot...")
        create_inversion_landscape_plot(inversion_df, joined_df, plots_dir, config)
        visualization_results['landscape'] = plots_dir / 'inversion_landscape.png'
        
        # 4. Summary plots
        logger.info("Creating summary plots...")
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        
        # Prepare results dictionary for summary plots
        results_dict = {
            'total_genes': len(joined_df),
            'flipped_genes': joined_df['is_flipped'].sum(),
            'chromosome_pairs': len(inversion_df),
            'pairs_with_inversions': len(inversion_df[inversion_df['flipped_genes'] > 0]),
            'flip_rate': joined_df['is_flipped'].mean(),
            'species1_name': config.get('first_species_name', 'Species1'),
            'species2_name': config.get('second_species_name', 'Species2')
        }
        
        create_statistics_summary_plot(results_dict, axes[0, 0])
        create_quality_summary_plot(results_dict, axes[0, 1])
        create_chromosome_mapping_overview(results_dict, axes[1, 0])
        create_synteny_summary_plot(results_dict, axes[1, 1])
        
        plt.tight_layout()
        summary_plot_path = plots_dir / 'analysis_summary.png'
        plt.savefig(summary_plot_path, dpi=300, bbox_inches='tight')
        plt.close()
        
        visualization_results['summary'] = summary_plot_path
        
        logger.info(f"Created {len(visualization_results)} visualization files")
        
    except Exception as e:
        logger.error(f"Visualization creation failed: {e}")
        # Continue with analysis even if visualizations fail
    
    return visualization_results


def save_results(inversion_df: pd.DataFrame,
                joined_df: pd.DataFrame,
                contextual_metrics: Dict,
                quality_metrics: Dict,
                output_dir: Path,
                config: Dict):
    """Save all analysis results to files"""
    
    data_dir = output_dir / 'data'
    data_dir.mkdir(exist_ok=True)
    
    # Save main dataframes
    inversion_df.to_csv(data_dir / 'inversion_analysis.csv', index=False)
    joined_df.to_csv(data_dir / 'joined_gene_data.csv', index=False)
    
    # Save contextual metrics as JSON
    import json
    with open(data_dir / 'contextual_metrics.json', 'w') as f:
        json.dump(contextual_metrics, f, indent=2, default=str)
    
    # Save quality metrics if available
    if quality_metrics:
        with open(data_dir / 'quality_metrics.json', 'w') as f:
            json.dump(quality_metrics, f, indent=2, default=str)
    
    logger.info(f"Results saved to {data_dir}")


def generate_analysis_report(inversion_df: pd.DataFrame,
                           joined_df: pd.DataFrame,
                           contextual_metrics: Dict,
                           quality_metrics: Dict,
                           output_dir: Path,
                           config: Dict) -> Path:
    """Generate comprehensive analysis report"""
    
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
    
    # Set random seeds for reproducibility
    import random
    random.seed(42)
    np.random.seed(42)
    
    try:
        # Run the analysis
        results = run_busco_inversion_analysis()
        
        # Print final summary
        print("\n" + "=" * 80)
        print("ANALYSIS COMPLETED SUCCESSFULLY!")
        print("=" * 80)
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