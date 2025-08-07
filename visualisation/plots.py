"""
visualisation module for the Genome Inversion Analyzer
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import logging
import seaborn as sns
from pathlib import Path
import subprocess
import tempfile
import shutil
import json
from typing import Dict, List, Optional, Tuple, Any, Union
import matplotlib.cm as cm
from matplotlib.colors import Normalize
from config.settings import CONFIG

logger = logging.getLogger(__name__)

try:
    from ete3 import Tree
    ETE3_AVAILABLE = True
except ImportError:
    ETE3_AVAILABLE = False
    logger.warning("ete3 not available - tree manipulation will be limited")

try:
    from ete3 import Tree, TreeStyle, NodeStyle, faces, AttrFace, CircleFace
    ETE3_AVAILABLE = True
except ImportError as e:
    logger.warning(f"ete3 import issue: {e}")
    try:
        from ete3 import Tree
        ETE3_AVAILABLE = True
    except ImportError:
        ETE3_AVAILABLE = False
        logger.warning("ete3 not available - tree manipulation will be limited")


def create_linearised_dotplot(joined_df, plots_dir, config=None):
    
    plot_config = {**CONFIG, **(config or {})}
    
    # Convert plots_dir to Path and ensure it exists
    plots_dir = Path(plots_dir)
    plots_dir.mkdir(parents=True, exist_ok=True)
    
    # Get chromosome sizes and linearise coordinates
    genome_data = get_chromosome_sizes_from_sheet()
    joined_df_linear, species1_data, species2_data = linearise_coordinates(joined_df, genome_data)
    
    # Create the plot
    fig, ax = plt.subplots(figsize=plot_config['figure_size'])
    
    # Split data by inversion status
    syntenic_data = joined_df_linear[joined_df_linear['is_flipped'] == False]
    inverted_data = joined_df_linear[joined_df_linear['is_flipped'] == True]
    
    # Plot syntenic points (blue) - using linearised coordinates
    if len(syntenic_data) > 0:
        ax.scatter(
            syntenic_data['linear_start1'], 
            syntenic_data['linear_start2'],
            c=plot_config['synteny_color'], 
            s=plot_config['point_size'],
            alpha=plot_config['point_alpha'],
            label='Syntenic',
            edgecolors='none'
        )
    
    # Plot inverted points (red) - using linearised coordinates
    if len(inverted_data) > 0:
        ax.scatter(
            inverted_data['linear_start1'], 
            inverted_data['linear_start2'],
            c=plot_config['inversion_color'], 
            s=plot_config['point_size'],
            alpha=plot_config['point_alpha'],
            label='Inverted',
            edgecolors='none'
        )
    
    # Add chromosome boundary lines and labels
    add_chromosome_boundaries(ax, species1_data, species2_data)
    
    add_chromosome_labels(ax, species1_data, species2_data, joined_df_linear)
        
    # Format axes to show positions in Mb
    ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f'{x/1e6:.0f}'))
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda y, p: f'{y/1e6:.0f}'))
    
    # Set axis labels
    ax.set_xlabel(f'{CONFIG["first_species_name"]} (Mb)', 
                 fontsize=plot_config['font_size_title'])
    ax.set_ylabel(f'{CONFIG["second_species_name"]} (Mb)', 
                 fontsize=plot_config['font_size_title'])
    
    # Format tick labels and add legend
    ax.tick_params(axis='both', labelsize=plot_config['font_size_labels'])
    ax.legend()

    # Apply styling
    ax.set_facecolor('white')
    ax.grid(False)
    for spine in ax.spines.values():
        spine.set_linewidth(0.8)
        spine.set_color('black')
    
    # Save plots
    png_path = plots_dir / "linearised_busco_dotplot.png"
    fig.savefig(png_path, dpi=plot_config['dpi'], bbox_inches='tight', facecolor='white')
    
    plt.tight_layout()
    return fig, ax



def add_chromosome_labels(ax, species1_data, species2_data, joined_df_linear):
    """Add chromosome labels on secondary axes"""
    
    species1_data = species1_data.sort_values('chromosome')
    chr1_positions = []
    chr1_labels = []
    cumulative_offset = 0
    
    for _, row in species1_data.iterrows():
        chr_size = row['chromsome_size_b']
        midpoint = cumulative_offset + (chr_size / 2)
        chr1_positions.append(midpoint)
        chr1_labels.append(row['chromosome'])
        cumulative_offset += chr_size
    
    species2_data = species2_data.sort_values('chromosome')
    chr2_positions = []
    chr2_labels = []
    cumulative_offset = 0
    
    for _, row in species2_data.iterrows():
        chr_size = row['chromsome_size_b']
        midpoint = cumulative_offset + (chr_size / 2)
        chr2_positions.append(midpoint)
        chr2_labels.append(row['chromosome'])
        cumulative_offset += chr_size
    
    # Create secondary axes
    ax2 = ax.secondary_xaxis('top')
    ax3 = ax.secondary_yaxis('right')
    
    # Set chromosome labels on top axis (species 1)
    ax2.set_xticks(chr1_positions)
    ax2.set_xticklabels(chr1_labels, 
                       rotation=45, 
                       ha='left',
                       fontsize=12)
    
    # Set chromosome labels on right axis (species 2)
    ax3.set_yticks(chr2_positions)
    ax3.set_yticklabels(chr2_labels, 
                       fontsize=12)
    
    # Remove tick marks (keep only labels)
    ax2.tick_params(length=0)
    ax3.tick_params(length=0)


def linearise_coordinates(joined_df, genome_data):
    """Convert per-chromosome coordinates to linearised genome coordinates"""
    
    # Create chromosome offset mapping for each species
    species1_offsets = {}
    species2_offsets = {}
    
    # Calculate cumulative offsets for species 1
    species1_data = genome_data[genome_data['species'] == CONFIG['first_species_name']]
    species1_data = species1_data.sort_values('chromosome')
    cumulative_offset = 0
    
    for _, row in species1_data.iterrows():
        species1_offsets[row['chromosome']] = cumulative_offset
        cumulative_offset += row['chromsome_size_b']
    
    # Calculate cumulative offsets for species 2  
    species2_data = genome_data[genome_data['species'] == CONFIG['second_species_name']]
    species2_data = species2_data.sort_values('chromosome')
    cumulative_offset = 0
    
    for _, row in species2_data.iterrows():
        species2_offsets[row['chromosome']] = cumulative_offset
        cumulative_offset += row['chromsome_size_b']
    
    # Add linearised coordinates to joined_df
    joined_df_linear = joined_df.copy()
    joined_df_linear['linear_start1'] = joined_df_linear.apply(
        lambda row: species1_offsets.get(row['chr1'], 0) + row['start1'], axis=1
    )
    joined_df_linear['linear_start2'] = joined_df_linear.apply(
        lambda row: species2_offsets.get(row['chr2'], 0) + row['start2'], axis=1
    )
    
    return joined_df_linear, species1_data, species2_data



def get_chromosome_sizes_from_sheet():
    """Get chromosome sizes from the Google Sheet"""
    sheet_url = "https://docs.google.com/spreadsheets/d/1K01wVWkMW-m6yT9zDX8gDekp-OECubE-9HcmD8RnmkM/edit?gid=1940964825#gid=1940964825"
    # Convert to CSV export URL
    csv_url = sheet_url.replace('/edit?gid=', '/export?format=csv&gid=')
    
    genome_data = pd.read_csv(csv_url)
    
    # Filter for your species
    target_species = ["Dioctria_linearis", "Dioctria_rufipes"]
    species_data = genome_data[genome_data['species'].isin(target_species)]
    
    return species_data[['species', 'chromosome', 'chromsome_size_b']]





def add_chromosome_boundaries(ax, species1_data, species2_data):
    """Add vertical and horizontal lines at chromosome boundaries"""
    
    # Calculate cumulative ends for vertical lines (species 1)
    species1_data = species1_data.sort_values('chromosome')
    cumulative_end = 0
    
    for i, (_, row) in enumerate(species1_data.iterrows()):
        cumulative_end += row['chromsome_size_b']
        # Add vertical line at end of each chromosome (except the last one)
        if i < len(species1_data) - 1:
            ax.axvline(x=cumulative_end, color='grey', linewidth=0.8, alpha=0.7)
    
    # Calculate cumulative ends for horizontal lines (species 2)  
    species2_data = species2_data.sort_values('chromosome')
    cumulative_end = 0
    
    for i, (_, row) in enumerate(species2_data.iterrows()):
        cumulative_end += row['chromsome_size_b']
        # Add horizontal line at end of each chromosome (except the last one)
        if i < len(species2_data) - 1:
            ax.axhline(y=cumulative_end, color='grey', linewidth=0.8, alpha=0.7)
            
            
            
            
            
            
            


def add_synteny_block_lines(df, ax, config): #Helper for create_busco_dotplot
    """Add lines connecting blocks on dotplot"""

    for (first_chr, second_chr), group in df.groupby(['chr1', 'chr2']):
        if len(group) >= 3:
            syntenic_group = group[group['strand1'] == group['strand2']]
            if len(syntenic_group) >= 2:

                sorted_group = syntenic_group.sort_values('first_order')
                ax.plot(sorted_group['first_order'], sorted_group['second_order'], 
                       'gray', alpha=0.5, linewidth=0.5, zorder=0)


def create_busco_dotplot(ortholog_df, plots_dir, config):
    """Create BUSCO-based dotplot showing gene order relationships"""
    if len(ortholog_df) == 0:
        return
    
    plt.figure(figsize=config.get('dotplot_size', (12, 10)))
    
    df = ortholog_df.copy()
    
    first_order = {}
    second_order = {}
    
    for chr_name in df['chr1'].unique():
        chr_genes = df[df['chr1'] == chr_name].sort_values('start1')
        first_order.update({row['busco_id']: i for i, (_, row) in enumerate(chr_genes.iterrows())})
    
    for chr_name in df['chr2'].unique():
        chr_genes = df[df['chr2'] == chr_name].sort_values('start2')
        second_order.update({row['busco_id']: i for i, (_, row) in enumerate(chr_genes.iterrows())})
    
    df['first_order'] = df['busco_id'].map(first_order)
    df['second_order'] = df['busco_id'].map(second_order)
    
    # Remove genes without order (shouldn't happen but safety check)
    # df = df.dropna(subset=['first_order', 'second_order'])
    
    if len(df) == 0:
        return
    
    syntenic_mask = df['strand1'] == df['strand2']
    
    syntenic_data = df[syntenic_mask]
    if len(syntenic_data) > 0:
        alpha_vals = syntenic_data['confidence'].values if config.get('confidence_alpha', True) and 'confidence' in df.columns else 0.7
        plt.scatter(syntenic_data['first_order'], syntenic_data['second_order'], 
                   c=config.get('synteny_color', '#1f77b4'), alpha=alpha_vals, s=20, label='Syntenic')
    
    inverted_data = df[~syntenic_mask]
    if len(inverted_data) > 0:
        alpha_vals = inverted_data['confidence'].values if config.get('confidence_alpha', True) and 'confidence' in df.columns else 0.7
        plt.scatter(inverted_data['first_order'], inverted_data['second_order'], 
                   c=config.get('inversion_color', '#d62728'), alpha=alpha_vals, s=20, label='Inverted')
    
    if config.get('show_synteny_blocks', True):
        add_synteny_block_lines(df, plt.gca(), config)
    

    plt.xlabel(f"Gene Order in {CONFIG['first_species_name']}")
    plt.ylabel(f"Gene Order in {CONFIG['second_species_name']}")
    plt.title('BUSCO Synteny Dotplot\n(Blue=Syntenic, Red=Inverted)')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(plots_dir, 'busco_dotplot.png'), dpi=300, bbox_inches='tight')
    plt.close()


# def create_ortholog_quality_plots(ortholog_df, plots_dir, config):
#     """Create plots showing ortholog quality metrics"""
#     fig, axes = plt.subplots(2, 2, figsize=(12, 10))
#     fig.suptitle('Ortholog Quality Assessment', fontsize=14, fontweight='bold')
    
#     # 1. Similarity distribution
#     if 'similarity' in ortholog_df.columns:
#         axes[0, 0].hist(ortholog_df['similarity'], bins=30, alpha=0.7, edgecolor='black', color='skyblue')
#         axes[0, 0].axvline(ortholog_df['similarity'].mean(), color='red', linestyle='--', label=f'Mean: {ortholog_df["similarity"].mean():.3f}')
#         axes[0, 0].set_xlabel('Similarity Score')
#         axes[0, 0].set_ylabel('Number of Orthologs')
#         axes[0, 0].set_title('Similarity Score Distribution')
#         axes[0, 0].legend()
#         axes[0, 0].grid(True, alpha=0.3)
    
#     # 2. Confidence distribution
#     if 'confidence' in ortholog_df.columns:
#         axes[0, 1].hist(ortholog_df['confidence'], bins=30, alpha=0.7, edgecolor='black', color='lightgreen')
#         axes[0, 1].axvline(ortholog_df['confidence'].mean(), color='red', linestyle='--', label=f'Mean: {ortholog_df["confidence"].mean():.3f}')
#         axes[0, 1].set_xlabel('Confidence Score')
#         axes[0, 1].set_ylabel('Number of Orthologs')
#         axes[0, 1].set_title('Confidence Score Distribution')
#         axes[0, 1].legend()
#         axes[0, 1].grid(True, alpha=0.3)
    
#     # 3. Alignment method usage
#     if 'alignment_method' in ortholog_df.columns:
#         method_counts = ortholog_df['alignment_method'].value_counts()
#         # Use Set3 colormap to generate distinct colors
#         cmap = cm.get_cmap('Set3')
#         colors = [cmap(i) for i in range(len(method_counts))]
#         wedges, texts, autotexts = axes[1, 0].pie(method_counts.values, labels=method_counts.index, 
#                                                   autopct='%1.1f%%', colors=colors)
#         axes[1, 0].set_title('Alignment Methods Used')
    
#     # 4. Gene length vs similarity scatter
#     if 'similarity' in ortholog_df.columns and 'first_length' in ortholog_df.columns:
#         scatter = axes[1, 1].scatter(ortholog_df['first_length'], ortholog_df['similarity'], 
#                                    alpha=0.6, s=10, c=ortholog_df.get('confidence', 'blue'))
#         axes[1, 1].set_xlabel('Gene Length (bp)')
#         axes[1, 1].set_ylabel('Similarity Score')
#         axes[1, 1].set_title('Gene Length vs Similarity')
#         axes[1, 1].grid(True, alpha=0.3)
        
#         if 'confidence' in ortholog_df.columns:
#             cbar = plt.colorbar(scatter, ax=axes[1, 1])
#             cbar.set_label('Confidence')
    
#     plt.tight_layout()
#     plt.savefig(os.path.join(plots_dir, 'ortholog_quality_assessment.png'), dpi=300, bbox_inches='tight')
#     plt.close()
    
#     logger.info("    ✓ Ortholog quality plots created")


# def create_inversion_landscape_plot(inversion_df, ortholog_df, plots_dir, config):
#     """Create inversion landscape showing inversions across chromosomes"""
#     if len(inversion_df) == 0:
#         return
    
#     fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 10))
#     fig.suptitle('Inversion Landscape Analysis', fontsize=14, fontweight='bold')
    
#     # 1. Inversion size distribution
#     if 'size_genes' in inversion_df.columns:
#         ax1.hist(inversion_df['size_genes'], bins=20, alpha=0.7, edgecolor='black', color='crimson')
#         ax1.set_xlabel('Inversion Size (genes)')
#         ax1.set_ylabel('Number of Inversions')
#         ax1.set_title('Inversion Size Distribution')
#         ax1.grid(True, alpha=0.3)
    
#     # 2. Inversion types
#     if 'inversion_type' in inversion_df.columns:
#         type_counts = inversion_df['inversion_type'].value_counts()
#         cmap = cm.get_cmap('Reds')
#         bars = ax2.bar(range(len(type_counts)), type_counts.values, color=cmap(0.7))
#         ax2.set_xticks(range(len(type_counts)))
#         ax2.set_xticklabels(type_counts.index, rotation=45, ha='right')
#         ax2.set_ylabel('Number of Inversions')
#         ax2.set_title('Inversion Types')
#         ax2.grid(True, alpha=0.3)
        
#         # Add value labels on bars
#         for bar, value in zip(bars, type_counts.values):
#             height = bar.get_height()
#             ax2.text(bar.get_x() + bar.get_width()/2., height + 0.1,
#                     f'{value}', ha='center', va='bottom')
    
#     plt.tight_layout()
#     plt.savefig(os.path.join(plots_dir, 'inversion_landscape.png'), dpi=300, bbox_inches='tight')
#     plt.close()
    
#     logger.info("    ✓ Inversion landscape plot created")


# def create_synteny_block_plots(ortholog_df, plots_dir, config):
#     """Create synteny block analysis plots from BUSCO ortholog data"""
#     import matplotlib.pyplot as plt
#     import os
#     import logging
#     import numpy as np
    
#     logger = logging.getLogger(__name__)
    
#     # FIXED: Create synteny block metrics from BUSCO ortholog data instead of expecting pre-calculated blocks
    
#     # Calculate block_size as number of genes per chromosome pair
#     chr_pair_counts = ortholog_df.groupby(['chr1', 'chr2']).size().reset_index(name='block_size')
    
#     # Calculate strand_consistency per chromosome pair
#     strand_consistency = ortholog_df.groupby(['chr1', 'chr2']).apply(
#         lambda x: (x['strand1'] == x['strand2']).mean()
#     ).reset_index(name='strand_consistency')
    
#     # Calculate position correlation per chromosome pair (if enough genes)
#     position_correlations = []
#     synteny_types = []
    
#     for (chr1, chr2), group in ortholog_df.groupby(['chr1', 'chr2']):
#         if len(group) >= 3:
#             corr = np.corrcoef(group['start1'], group['start2'])[0, 1]
#             position_correlations.append(corr if not np.isnan(corr) else 0)
#         else:
#             position_correlations.append(0)
            
#         # Determine synteny type based on strand consistency and correlation
#         strand_cons = (group['strand1'] == group['strand2']).mean()
#         if strand_cons > 0.8:
#             synteny_types.append('Conserved')
#         elif strand_cons < 0.2:
#             synteny_types.append('Inverted')
#         else:
#             synteny_types.append('Mixed')
    
#     # Create synthetic synteny_df from BUSCO data
#     synteny_df = chr_pair_counts.copy()
#     synteny_df['strand_consistency'] = strand_consistency['strand_consistency'].values
#     synteny_df['position_correlation'] = position_correlations
#     synteny_df['synteny_type'] = synteny_types
    
#     # Rest of the function stays the same
#     fig, axes = plt.subplots(2, 2, figsize=(12, 10))
#     fig.suptitle('Synteny Block Analysis', fontsize=14, fontweight='bold')
    
#     # 1. Block size distribution
#     axes[0, 0].hist(synteny_df['block_size'], bins=20, alpha=0.7, edgecolor='black', color='orange')
#     axes[0, 0].set_xlabel('Block Size (genes)')
#     axes[0, 0].set_ylabel('Number of Blocks')
#     axes[0, 0].set_title('Synteny Block Size Distribution')
#     axes[0, 0].grid(True, alpha=0.3)
    
#     # 2. Synteny types
#     if 'synteny_type' in synteny_df.columns:
#         type_counts = synteny_df['synteny_type'].value_counts()
#         cmap = cm.get_cmap('Pastel1')
#         colors = [cmap(i) for i in range(len(type_counts))]
#         axes[0, 1].pie(type_counts.values, labels=type_counts.index, autopct='%1.1f%%', colors=colors)
#         axes[0, 1].set_title('Synteny Types')
    
#     # 3. Position correlation distribution
#     if 'position_correlation' in synteny_df.columns:
#         axes[1, 0].hist(synteny_df['position_correlation'], bins=20, alpha=0.7, edgecolor='black', color='lightcoral')
#         axes[1, 0].set_xlabel('Position Correlation')
#         axes[1, 0].set_ylabel('Number of Blocks')
#         axes[1, 0].set_title('Position Correlation Distribution')
#         axes[1, 0].grid(True, alpha=0.3)
    
#     # 4. Strand consistency distribution
#     if 'strand_consistency' in synteny_df.columns:
#         axes[1, 1].hist(synteny_df['strand_consistency'], bins=20, alpha=0.7, edgecolor='black', color='mediumpurple')
#         axes[1, 1].set_xlabel('Strand Consistency')
#         axes[1, 1].set_ylabel('Number of Blocks')
#         axes[1, 1].set_title('Strand Consistency Distribution')
#         axes[1, 1].grid(True, alpha=0.3)
    
#     plt.tight_layout()
#     plt.savefig(os.path.join(plots_dir, 'synteny_block_analysis.png'), dpi=300, bbox_inches='tight')
#     plt.close()
    
#     logger.info("    ✓ Synteny block plots created")


# def create_statistics_summary_plot(results_dict, ax):
#     """Create overall statistics summary"""
#     stats = []
#     labels = []
    
#     if 'ortholog_df' in results_dict and not results_dict['ortholog_df'].empty:
#         stats.extend([
#             len(results_dict['ortholog_df']),
#             results_dict['ortholog_df']['similarity'].mean() if 'similarity' in results_dict['ortholog_df'].columns else 0,
#             results_dict['ortholog_df']['confidence'].mean() if 'confidence' in results_dict['ortholog_df'].columns else 0
#         ])
#         labels.extend(['Orthologs', 'Avg Similarity', 'Avg Confidence'])
    
#     if 'synteny_df' in results_dict:
#         stats.append(len(results_dict['synteny_df']))
#         labels.append('Synteny Blocks')
    
#     if 'inversion_df' in results_dict:
#         stats.append(len(results_dict['inversion_df']))
#         labels.append('Inversions')
    
#     if 'rearrangement_df' in results_dict:
#         stats.append(len(results_dict['rearrangement_df']))
#         labels.append('Rearrangements')
    
#     if stats:
#         # normalise stats for display (except counts)
#         display_stats = []
#         for i, (stat, label) in enumerate(zip(stats, labels)):
#             if 'Avg' in label:
#                 display_stats.append(stat * 100)  # Convert to percentage
#             else:
#                 display_stats.append(stat)
        
#         cmap = cm.get_cmap('Set3')
#         colors = [cmap(i) for i in range(len(stats))]
#         bars = ax.bar(labels, display_stats, color=colors, alpha=0.8)
        
#         # Add value labels on bars
#         for bar, stat, label in zip(bars, stats, labels):
#             height = bar.get_height()
#             if 'Avg' in label:
#                 ax.text(bar.get_x() + bar.get_width()/2., height + 0.5,
#                        f'{stat:.3f}', ha='center', va='bottom', fontsize=10)
#             else:
#                 ax.text(bar.get_x() + bar.get_width()/2., height + 0.5,
#                        f'{int(stat)}', ha='center', va='bottom', fontsize=10)
        
#         ax.set_title('Analysis Summary Statistics')
#         ax.set_ylabel('Count / Score')
#         plt.setp(ax.get_xticklabels(), rotation=45, ha='right')


# def create_quality_summary_plot(results_dict, ax):
#     """Create quality summary plot"""
#     if 'first_quality' in results_dict and 'second_quality' in results_dict:
#         genomes = ['First Genome', 'Second Genome']
#         scores = [results_dict['first_quality']['quality_score'], 
#                  results_dict['second_quality']['quality_score']]
#         classes = [results_dict['first_quality']['quality_class'], 
#                   results_dict['second_quality']['quality_class']]
        
#         # Color by quality class
#         color_map = {'high': 'green', 'medium': 'orange', 'low': 'red', 'fragmented': 'darkred'}
#         colors = [color_map.get(cls, 'gray') for cls in classes]
        
#         bars = ax.bar(genomes, scores, color=colors, alpha=0.7)
        
#         # Add score labels and quality class
#         for bar, score, cls in zip(bars, scores, classes):
#             height = bar.get_height()
#             ax.text(bar.get_x() + bar.get_width()/2., height + 0.02,
#                    f'{score:.3f}\n({cls})', ha='center', va='bottom', fontsize=10)
        
#         ax.set_title('Assembly Quality Scores')
#         ax.set_ylabel('Quality Score')
#         ax.set_ylim(0, 1.1)
#         ax.grid(True, alpha=0.3)


# def create_chromosome_mapping_overview(results_dict, ax):
#     """Create chromosome mapping overview"""
#     if 'ortholog_df' not in results_dict or results_dict['ortholog_df'].empty:
#         ax.text(0.5, 0.5, 'No ortholog data available', ha='center', va='center', transform=ax.transAxes)
#         ax.set_title('Chromosome Mapping Overview')
#         return
    
#     df = results_dict['ortholog_df']
    
#     # Create chromosome mapping matrix
#     first_chrs = sorted(df['chr1'].unique())
#     second_chrs = sorted(df['chr2'].unique())
    
#     # Create mapping matrix
#     mapping_matrix = np.zeros((len(first_chrs), len(second_chrs)))
    
#     for i, first_chr in enumerate(first_chrs):
#         for j, second_chr in enumerate(second_chrs):
#             count = len(df[(df['chr1'] == first_chr) & (df['chr2'] == second_chr)])
#             mapping_matrix[i, j] = count
    
#     # Create heatmap
#     im = ax.imshow(mapping_matrix, cmap='YlOrRd', aspect='auto')
    
#     # Set ticks and labels
#     ax.set_xticks(range(len(second_chrs)))
#     ax.set_yticks(range(len(first_chrs)))
#     ax.set_xticklabels([chr[:10] for chr in second_chrs], rotation=45, ha='right')
#     ax.set_yticklabels([chr[:10] for chr in first_chrs])
    
#     # Add colorbar
#     cbar = plt.colorbar(im, ax=ax, fraction=0.02)
#     cbar.set_label('Number of Orthologs')
    
#     ax.set_title('Chromosome Mapping Matrix')
#     ax.set_xlabel('Second Genome Chromosomes')
#     ax.set_ylabel('First Genome Chromosomes')


# def create_synteny_summary_plot(results_dict, ax):
#     """Create synteny summary plot"""
#     if 'synteny_df' not in results_dict or results_dict['synteny_df'].empty:
#         ax.text(0.5, 0.5, 'No synteny\ndata', ha='center', va='center', transform=ax.transAxes)
#         ax.set_title('Synteny Summary')
#         return
    
#     synteny_df = results_dict['synteny_df']
#     if 'synteny_type' in synteny_df.columns:
#         type_counts = synteny_df['synteny_type'].value_counts()
#         cmap = cm.get_cmap('Pastel1')
#         colors = [cmap(i) for i in range(len(type_counts))]
#         ax.pie(type_counts.values, labels=type_counts.index, autopct='%1.0f%%', 
#                colors=colors, textprops={'fontsize': 8})
#     ax.set_title('Synteny Types', fontsize=10)


# def create_inversion_summary_plot(results_dict, ax):
#     """Create inversion summary plot"""
#     if 'inversion_df' not in results_dict or results_dict['inversion_df'].empty:
#         ax.text(0.5, 0.5, 'No inversion\ndata', ha='center', va='center', transform=ax.transAxes)
#         ax.set_title('Inversion Summary')
#         return
    
#     inversion_df = results_dict['inversion_df']
#     if 'size_genes' in inversion_df.columns:
#         ax.hist(inversion_df['size_genes'], bins=10, alpha=0.7, color='crimson', edgecolor='black')
#         ax.set_xlabel('Size (genes)', fontsize=8)
#         ax.set_ylabel('Count', fontsize=8)
#     ax.set_title('Inversion Sizes', fontsize=10)

    
# def create_synteny_plots(all_results: Dict, output_dir: Path, config: Dict) -> Dict:
#     """Create publication-quality synteny plots using synteny_plotter"""
    
#     # FIXED: Converted from class method to standalone function - removed 'self' parameter and added 'config' parameter
#     import logging
#     from pathlib import Path
#     import pandas as pd
#     import subprocess
#     import shutil
    
#     logger = logging.getLogger(__name__)
#     logger.info("  📊 Creating curved synteny plots...")
    
#     synteny_dir = output_dir / 'synteny_curves'
#     synteny_dir.mkdir(exist_ok=True)
    
#     # FIXED: Use config parameter instead of self.config
#     synteny_plotter_dir = config.get('external_tools', {}).get('synteny_plotter')
    
#     if synteny_plotter_dir:
#         synteny_plotter_path = Path(synteny_plotter_dir).resolve() / 'scripts' / 'generate_synteny_plot.R'
#         logger.info(f"DEBUG: Looking for R script at: {synteny_plotter_path.absolute()}")
#     else:
#         synteny_plotter_path = None
    
#     if not synteny_plotter_path or not synteny_plotter_path.exists():
#         logger.warning(f"Synteny plotter R script not found: {synteny_plotter_path}")
#         logger.info("  Nope 1 Failed: Creating fallback matplotlib synteny plots...")
#         # FIXED: Call standalone function instead of self.method
#         return create_fallback_synteny_plots.nope(all_results, synteny_dir, config)
    
#     results = {}
    
#     for pair_name, pair_data in all_results.items():
#         if 'full_results' not in pair_data:
#             continue
            
#         try:
#             species1, species2 = pair_data['species_pair']
#             logger.info(f"    • Creating synteny plot: {pair_name}")
            
#             # Convert data to synteny_plotter format
#             # FIXED: Call standalone helper function instead of self.method
#             temp_files = _prepare_synteny_plotter_input(
#                 pair_data['ortholog_df'], 
#                 pair_data['inversion_df'],
#                 species1, species2,
#                 synteny_dir,
#                 config  # FIXED: Pass config parameter
#             )
            
#             # Call synteny_plotter
#             output_plot = synteny_dir / f'curved_synteny_{pair_name}.png'
            
#             # FIXED: Call standalone helper function instead of self.method
#             success = _run_synteny_plotter(temp_files, output_plot, species1, species2, config)
            
#             if success:
#                 results[pair_name] = output_plot
#                 logger.info(f"      Created: {output_plot}")
#             else:
#                 logger.warning(f"      ❌ Nope 2 Failed: {pair_name}, creating fallback")
#                 # FIXED: Call standalone function instead of self.method  
#                 fallback_plot = create_single_fallback_plot.nope(
#                     pair_data['ortholog_df'], pair_data['inversion_df'], 
#                     species1, species2, synteny_dir, config
#                 )
#                 if fallback_plot:
#                     results[pair_name] = fallback_plot
            
#         except Exception as e:
#             logger.error(f"      ❌ Error creating synteny plot for {pair_name}: {e}")
    
#     logger.info(f"  ✅ Created {len(results)} synteny plots")
#     return results


# def _prepare_synteny_plotter_input(ortholog_df: pd.DataFrame, inversion_df: pd.DataFrame,
#                                 species1: str, species2: str, output_dir: Path, config: Dict) -> Dict:
#     """Convert ortholog data to synteny_plotter BUSCO TSV format - EXACT match to working test files"""
    
#     synteny_plotter_dir = Path(config['external_tools']['synteny_plotter'])
    
#     # Create temp directory inside synteny_plotter
#     temp_dir = synteny_plotter_dir / f'temp_{species1}_vs_{species2}'
#     temp_dir.mkdir(exist_ok=True)
    
#     # The R script expects BUSCO files that EXACTLY match the working test format
#     busco1_data = []
#     busco2_data = []
    
#     for _, row in ortholog_df.iterrows():
#         # Reference species - match EXACT format of working Melitaea_cinxia.tsv
#         busco1_data.append([
#             row['busco_id'],      # Column 1: Busco id
#             'Complete',           # Column 2: Status
#             row['chr1'],         # Column 3: Sequence (chromosome)
#             int(row['start1']),  # Column 4: Gene Start
#             int(row['end1']),    # Column 5: Gene End
#             row['strand1'],      # Column 6: Strand
#             500.0,               # Column 7: Score (dummy value)
#             abs(int(row['end1']) - int(row['start1'])),  # Column 8: Length
#             f"https://www.orthodb.org/v10?query={row['busco_id']}",  # Column 9: OrthoDB url
#             row['busco_id']      # Column 10: Description
#         ])
        
#         # Query species - same format
#         busco2_data.append([
#             row['busco_id'],      # Column 1: Busco id
#             'Complete',           # Column 2: Status
#             row['chr2'],         # Column 3: Sequence (chromosome)
#             int(row['start2']),  # Column 4: Gene Start
#             int(row['end2']),    # Column 5: Gene End
#             row['strand2'],      # Column 6: Strand
#             500.0,               # Column 7: Score (dummy value)
#             abs(int(row['end2']) - int(row['start2'])),  # Column 8: Length
#             f"https://www.orthodb.org/v10?query={row['busco_id']}",  # Column 9: OrthoDB url
#             row['busco_id']      # Column 10: Description
#         ])
    
#     # Write BUSCO TSV files in EXACT format as working test files
#     busco1_file = temp_dir / f'{species1}.tsv'
#     busco2_file = temp_dir / f'{species2}.tsv'
    
#     # Write with EXACT BUSCO headers as in working test file
#     with open(busco1_file, 'w') as f:
#         f.write("# BUSCO version is: 5.0.0 \n")
#         f.write("# The lineage dataset is: custom (Creation date: 2025-07-25, number of BUSCOs: custom)\n")
#         f.write("# Busco id\tStatus\tSequence\tGene Start\tGene End\tStrand\tScore\tLength\tOrthoDB url\tDescription\n")
#         for row in busco1_data:
#             f.write('\t'.join(map(str, row)) + '\n')

#     with open(busco2_file, 'w') as f:
#         f.write("# BUSCO version is: 5.0.0 \n")
#         f.write("# The lineage dataset is: custom (Creation date: 2025-07-25, number of BUSCOs: custom)\n")
#         f.write("# Busco id\tStatus\tSequence\tGene Start\tGene End\tStrand\tScore\tLength\tOrthoDB url\tDescription\n")
#         for row in busco2_data:
#             f.write('\t'.join(map(str, row)) + '\n')
    
#     # Create chromosome info files - EXACT format as working test (only 3 columns!)
#     # Working format: chr, order, invert (NO length, NO direction columns)
    
#     chr1_info = {}
#     chr1_names = sorted(ortholog_df['chr1'].unique())
#     for i, chr_name in enumerate(chr1_names):
#         chr1_info[chr_name] = {
#             'chr': chr_name,
#             'order': i + 1,
#             'invert': 'F'  # Use 'F' not 'FALSE' to match working format
#         }
    
#     # For query species - same 3-column format
#     chr2_info = {}
#     chr2_names = sorted(ortholog_df['chr2'].unique())
#     for i, chr_name in enumerate(chr2_names):
#         chr2_info[chr_name] = {
#             'chr': chr_name,
#             'order': i + 1,
#             'invert': 'F'  # Use 'F' not 'FALSE'
#         }
    
#     # Write chromosome info files with EXACT 3-column format
#     chrom1_file = temp_dir / f'{species1}_info.tsv'
#     chrom2_file = temp_dir / f'{species2}_info.tsv'
    
#     # Convert to DataFrame with exact column order
#     chr1_df = pd.DataFrame(list(chr1_info.values()))
#     chr2_df = pd.DataFrame(list(chr2_info.values()))
    
#     # Ensure EXACT column order as working test: chr, order, invert
#     column_order = ['chr', 'order', 'invert']
#     chr1_df = chr1_df[column_order]
#     chr2_df = chr2_df[column_order]
    
#     chr1_df.to_csv(chrom1_file, sep='\t', index=False, header=True)
#     chr2_df.to_csv(chrom2_file, sep='\t', index=False, header=True)
    
#     logger.info(f"Created BUSCO files matching working test format:")
#     logger.info(f"  Reference: {busco1_file} ({len(busco1_data)} genes)")
#     logger.info(f"  Query: {busco2_file} ({len(busco2_data)} genes)")
#     logger.info(f"  Ref chromosomes: {chrom1_file} ({len(chr1_info)} chromosomes)")
#     logger.info(f"  Query chromosomes: {chrom2_file} ({len(chr2_info)} chromosomes)")
    
#     # Debug: Show format matches
#     logger.info("Verifying file format matches working test:")
#     with open(busco1_file, 'r') as f:
#         lines = f.readlines()
#         logger.info(f"  Header line 1: {lines[0].strip()}")
#         logger.info(f"  Header line 2: {lines[1].strip()}")
#         logger.info(f"  Header line 3: {lines[2].strip()}")
#         if len(lines) > 3:
#             logger.info(f"  Data line 1: {lines[3].strip()}")
    
#     return {
#         'busco1': busco1_file,
#         'busco2': busco2_file, 
#         'chrom1': chrom1_file,
#         'chrom2': chrom2_file,
#         'temp_dir': temp_dir
#     }



# def _run_synteny_plotter(temp_files: Dict, output_plot: Path, 
#                     species1: str, species2: str, config: Dict) -> bool:
#     """Run the R synteny plotter with correct argument format"""
    
#     synteny_plotter_dir = Path(config['external_tools']['synteny_plotter'])
#     r_script_path = synteny_plotter_dir / 'scripts' / 'generate_synteny_plot.R'
    
#     if not r_script_path.exists():
#         logger.error(f"R script not found: {r_script_path}")
#         return False
    
#     # Build R command with CORRECT arguments based on the R script:
#     # The R script expects: -busco1, -chrom1, -busco_list, -chrom_list, -o, -f, -g, -alpha
    
#     cmd = [
#         'Rscript', 'scripts/generate_synteny_plot.R',
#         '-busco1', str(temp_files['busco1'].relative_to(synteny_plotter_dir)),      # Reference BUSCO file
#         '-chrom1', str(temp_files['chrom1'].relative_to(synteny_plotter_dir)),      # Reference chromosome info
#         '-busco_list', str(temp_files['busco2'].relative_to(synteny_plotter_dir)),  # Query BUSCO file(s) - list
#         '-chrom_list', str(temp_files['chrom2'].relative_to(synteny_plotter_dir)),  # Query chromosome info(s) - list
#         '-o', str(output_plot.with_suffix('').name)  # Output prefix (no extension)
#     ]
    
#     # Add optional parameters from config
#     synteny_config = config.get('synteny_visualisation', {})
#     options = synteny_config.get('options', {})
    
#     # Add filter threshold (minimum BUSCOs per chromosome)
#     filter_threshold = options.get('filter_threshold', 3)  # Default 3, R script default is 5
#     cmd.extend(['-f', str(filter_threshold)])
    
#     # Add gap parameter (distance between chromosomes)
#     gap = options.get('gap', 6)  # Default 6, matches R script default
#     cmd.extend(['-g', str(gap)])
    
#     # Add alpha parameter (transparency)
#     alpha = options.get('alpha', 0)  # Default 0, matches R script default
#     cmd.extend(['-alpha', str(alpha)])
    
#     try:
#         logger.info(f"Running synteny_plotter command:")
#         logger.info(f"  {' '.join(cmd)}")
#         logger.info(f"Working directory: {synteny_plotter_dir}")
        
#         result = subprocess.run(
#             cmd, 
#             capture_output=True, 
#             text=True, 
#             check=True,
#             cwd=str(synteny_plotter_dir)  # CRITICAL: Run from synteny_plotter directory
#         )
        
#         logger.info("R script completed successfully")
#         if result.stdout:
#             logger.info(f"R stdout: {result.stdout}")
        
#         # R script creates PDF output, look for it
#         possible_outputs = [
#             synteny_plotter_dir / f"{output_plot.with_suffix('').name}.pdf",
#             synteny_plotter_dir / f"{output_plot.with_suffix('').name}.png",
#             synteny_plotter_dir / f"{output_plot.with_suffix('').name}.svg"
#         ]
        
#         output_created = False
#         for possible_output in possible_outputs:
#             if possible_output.exists():
#                 # Move to expected location with PNG extension (or keep original)
#                 final_output = output_plot.with_suffix('.png' if possible_output.suffix == '.pdf' else possible_output.suffix)
#                 if possible_output.suffix == '.pdf':
#                     # Convert PDF to PNG if needed (or just copy)
#                     shutil.copy(str(possible_output), str(final_output.with_suffix('.pdf')))
#                     logger.info(f"PDF output saved: {final_output.with_suffix('.pdf')}")
#                 else:
#                     shutil.move(str(possible_output), str(final_output))
#                     logger.info(f"Moved output: {possible_output.name} -> {final_output}")
#                 output_created = True
#                 break
        
#         # Clean up temp directory
#         # if temp_files['temp_dir'].exists():
#         #     shutil.rmtree(temp_files['temp_dir'])
#         #     logger.info("Cleaned up temp files")
        
#         if not output_created:
#             logger.warning("R script completed but no output file found")
#             # List files in directory for debugging
#             created_files = list(synteny_plotter_dir.glob(f"{output_plot.with_suffix('').name}*"))
#             if created_files:
#                 logger.info(f"Files found: {[f.name for f in created_files]}")
        
#         return output_created
        
#     except subprocess.CalledProcessError as e:
#         logger.error(f"R synteny plotter failed with error:")
#         logger.error(f"Return code: {e.returncode}")
#         logger.error(f"STDERR: {e.stderr}")
#         logger.error(f"STDOUT: {e.stdout}")
#         logger.error(f"Command: {' '.join(cmd)}")
        
#         # Clean up temp directory even on failure
#         # if temp_files['temp_dir'].exists():
#         #     shutil.rmtree(temp_files['temp_dir'])
        
#         return False
#     except Exception as e:
#         logger.error(f"Error running synteny plotter: {e}")
        
#         # Clean up temp directory
#         # if temp_files['temp_dir'].exists():
#         #     shutil.rmtree(temp_files['temp_dir'])
        
#         return False

# def create_fallback_synteny_plots(all_results: Dict, output_dir: Path, config: Dict) -> Dict:
#     """Create matplotlib-based synteny plots as fallback"""
    
#     # FIXED: Converted from class method to standalone function - removed 'self' parameter and added 'config' parameter
#     import logging
#     from pathlib import Path
#     from typing import Dict
    
#     logger = logging.getLogger(__name__)
#     logger.info("    Creating matplotlib fallback synteny plots...")
    
#     results = {}
    
#     for pair_name, pair_data in all_results.items():
#         if 'full_results' not in pair_data:
#             continue
            
#         try:
#             species1, species2 = pair_data['species_pair']
            
#             # FIXED: Call standalone function instead of self.method
#             plot_file = create_single_fallback_plot.nope(
#                 pair_data['ortholog_df'], 
#                 pair_data['inversion_df'],
#                 species1, species2, 
#                 output_dir,
#                 config  # FIXED: Pass config parameter
#             )
            
#             if plot_file:
#                 results[pair_name] = plot_file
#                 logger.info(f"    ✅ Fallback plot: {pair_name}")
                
#         except Exception as e:
#             logger.error(f"    ❌ Fallback plot failed for {pair_name}: {e}")
    
#     return results

# def create_single_fallback_plot(ortholog_df: pd.DataFrame, inversion_df: pd.DataFrame,
#                             species1: str, species2: str, output_dir: Path, config: Dict) -> Optional[Path]:
#     """Create a single matplotlib synteny plot with curved connections"""
    
#     # FIXED: Converted from class method to standalone function - removed 'self' parameter and added 'config' parameter
#     import matplotlib.pyplot as plt
#     import matplotlib.cm as cm
#     from matplotlib.colors import normalize
#     import numpy as np
#     import pandas as pd
#     from pathlib import Path
#     from typing import Optional, Dict
    
#     if ortholog_df.empty:
#         print(f"    No ortholog data for {species1} vs {species2}")
#         return None
    
#     try:
#         print(f"    Creating matplotlib plot for {species1} vs {species2}")
#         fig, ax = plt.subplots(figsize=(14, 10))
        
#         # Use correct column names directly
#         chr1_list = sorted(ortholog_df['chr1'].unique())
#         chr2_list = sorted(ortholog_df['chr2'].unique())
        
#         print(f"    Chromosomes: {len(chr1_list)} vs {len(chr2_list)}")
        
#         # Create chromosome positions
#         chr1_positions = {chr_name: i for i, chr_name in enumerate(chr1_list)}
#         chr2_positions = {chr_name: i for i, chr_name in enumerate(chr2_list)}
        
#         y1_base = 0.7
#         y2_base = 0.3
        
#         # Plot chromosome backbones
#         for i, chr_name in enumerate(chr1_list):
#             ax.plot([0, 1], [y1_base + i*0.05, y1_base + i*0.05], 'k-', linewidth=3, alpha=0.7)
#             display_name = str(chr_name)[:10] + "..." if len(str(chr_name)) > 10 else str(chr_name)
#             ax.text(-0.05, y1_base + i*0.05, display_name, ha='right', va='center', fontsize=8)
        
#         for i, chr_name in enumerate(chr2_list):
#             ax.plot([0, 1], [y2_base - i*0.05, y2_base - i*0.05], 'k-', linewidth=3, alpha=0.7)
#             display_name = str(chr_name)[:10] + "..." if len(str(chr_name)) > 10 else str(chr_name)
#             ax.text(-0.05, y2_base - i*0.05, display_name, ha='right', va='center', fontsize=8)
        
#         # Plot connections (limit to avoid overcrowding)
#         connections_plotted = 0
#         max_connections = min(200, len(ortholog_df))  # Limit for readability
        
#         for _, row in ortholog_df.head(max_connections).iterrows():
#             try:
#                 chr1_idx = chr1_positions[row['chr1']]
#                 chr2_idx = chr2_positions[row['chr2']]
                
#                 y1 = y1_base + chr1_idx * 0.05
#                 y2 = y2_base - chr2_idx * 0.05
                
#                 # Simple position mapping
#                 x1 = 0.1 + 0.8 * (connections_plotted / max_connections)
#                 x2 = 0.1 + 0.8 * (connections_plotted / max_connections)
                
#                 # FIXED: Calculate similarity from BUSCO data instead of expecting 'similarity' column
#                 if 'similarity' in row:
#                     similarity = row['similarity']
#                 else:
#                     # Calculate similarity based on strand consistency and gene length ratio
#                     strand_match = 1.0 if row['strand1'] == row['strand2'] else 0.5
#                     length1 = row['end1'] - row['start1']
#                     length2 = row['end2'] - row['start2']
#                     length_ratio = min(length1, length2) / max(length1, length2) if max(length1, length2) > 0 else 1.0
#                     similarity = (strand_match + length_ratio) / 2.0
                
#                 # Color by similarity
#                 color = cm.get_cmap('viridis')(similarity)
#                 alpha = 0.6
                
#                 # Create curved connection
#                 x_curve = np.linspace(x1, x2, 20)
#                 y_curve = []
                
#                 for x in x_curve:
#                     t = (x - x1) / (x2 - x1) if x2 != x1 else 0
#                     y = y1 * (1 - t) + y2 * t + 0.1 * np.sin(np.pi * t)
#                     y_curve.append(y)
                
#                 ax.plot(x_curve, y_curve, color=color, alpha=alpha, linewidth=1)
#                 connections_plotted += 1
                
#             except (KeyError, ValueError) as e:
#                 continue
        
#         print(f"    Plotted {connections_plotted} connections")
        
#         # Styling
#         ax.set_xlim(-0.2, 1.1)
#         ax.set_ylim(min(y2_base - len(chr2_list)*0.05, y1_base - 0.1), 
#                 max(y1_base + len(chr1_list)*0.05, y2_base + 0.1))
        
#         ax.set_title(f'Synteny Plot: {species1} vs {species2}', fontsize=16, fontweight='bold')
#         ax.text(0.5, y1_base + len(chr1_list)*0.05 + 0.05, species1, 
#             ha='center', va='bottom', fontsize=14, fontweight='bold')
#         ax.text(0.5, y2_base - len(chr2_list)*0.05 - 0.05, species2, 
#             ha='center', va='top', fontsize=14, fontweight='bold')
        
#         ax.axis('off')
        
#         # Add colorbar
#         sm = cm.ScalarMappable(cmap=cm.get_cmap('viridis'), norm=normalize(vmin=0, vmax=1))
#         sm.set_array([])
#         cbar = plt.colorbar(sm, ax=ax, shrink=0.6, pad=0.02)
#         cbar.set_label('Similarity', rotation=270, labelpad=15)
        
#         plt.tight_layout()
        
#         # Save plot
#         plot_file = output_dir / f'matplotlib_synteny_{species1}_vs_{species2}.png'
#         plt.savefig(plot_file, dpi=300, bbox_inches='tight')
#         plt.close()
        
#         print(f"    ✅ Saved: {plot_file}")
#         return plot_file
        
#     except Exception as e:
#         print(f"    ❌ Plot creation failed: {e}")
#         import traceback
#         traceback.print_exc()
#         if 'fig' in locals():
#             plt.close()
#         return None

# def create_annotated_phylogeny(all_results: Dict, species_stats: Dict, 
#                                 output_dir: Path, config: Dict) -> Dict:
#     """Create annotated phylogenetic tree from Diptera tree"""
    
#     # FIXED: Converted from class method to standalone function - removed 'self' parameter and added 'config' parameter
#     import logging
#     from pathlib import Path
#     from typing import Dict
    
#     logger = logging.getLogger(__name__)
#     logger.info("  🌳 Creating annotated phylogenetic tree...")
    
#     tree_dir = output_dir / 'annotated_trees'
#     tree_dir.mkdir(exist_ok=True)
    
#     # FIXED: Use config parameter instead of self.config
#     tree_config = config.get('tree_annotation', {})
#     source_tree_path = tree_config.get('source_tree_path')
    
#     if not source_tree_path or not Path(source_tree_path).exists():
#         logger.error(f"Source tree not found: {source_tree_path}")
#         return {}
    
#     try:
#         # Load and prune tree
#         target_species = list(species_stats.keys())
        
#         # FIXED: Call standalone helper function instead of self.method
#         pruned_tree = _prune_diptera_tree(source_tree_path, target_species, config)
        
#         if not pruned_tree:
#             logger.error("Failed to prune tree")
#             return {}
        
#         # Calculate inversion annotations
#         # FIXED: Call standalone helper function instead of self.method
#         inversion_annotations = _calculate_inversion_annotations(all_results, species_stats, config)
        
#         # Annotate tree with inversion data
#         # FIXED: Call standalone helper function instead of self.method
#         annotated_tree = _annotate_tree_with_inversions(pruned_tree, inversion_annotations, config)
        
#         # Create visualisations
#         results = {}
        
#         # Save annotated tree
#         tree_file = tree_dir / 'annotated_bibionidae_tree.newick'
        
#         # FIXED: Check ETE3_AVAILABLE as standalone variable or import check
#         try:
#             from ete3 import Tree
#             ETE3_AVAILABLE = True
#         except ImportError:
#             ETE3_AVAILABLE = False
            
#         if ETE3_AVAILABLE and annotated_tree:
#             annotated_tree.write(outfile=str(tree_file))
#             results['newick'] = tree_file
        
#         # Create tree plot
#         plot_file = tree_dir / 'annotated_tree_plot.png'
#         # FIXED: Call standalone function instead of self.method
#         create_tree_plot(annotated_tree, plot_file, inversion_annotations, config)
#         results['plot'] = plot_file
        
#         # Create tree heatmap
#         heatmap_file = tree_dir / 'tree_inversion_heatmap.png'
#         # FIXED: Call standalone function instead of self.method
#         create_tree_heatmap(inversion_annotations, heatmap_file, config)
#         results['heatmap'] = heatmap_file
        
#         logger.info(f"  ✅ Created annotated tree: {len(results)} outputs")
#         return results
        
#     except Exception as e:
#         logger.error(f"Tree annotation failed: {e}")
#         import traceback
#         traceback.print_exc()
#         return {}

# def _prune_diptera_tree(tree_path: str, target_species: List[str], config: Dict) -> Optional['Tree']:
#     """Prune large Diptera tree to target species"""
    
#     # FIXED: Converted from class method to standalone function - removed 'self' parameter and added 'config' parameter
#     import logging
#     from typing import List, Optional
    
#     logger = logging.getLogger(__name__)
    
#     # FIXED: Check ETE3 availability as import check instead of class variable
#     try:
#         from ete3 import Tree
#         ETE3_AVAILABLE = True
#     except ImportError:
#         ETE3_AVAILABLE = False
    
#     if not ETE3_AVAILABLE:
#         logger.warning("ete3 not available - creating simple tree structure")
#         # FIXED: Call standalone function instead of self.method
#         return create_simple_tree.nope(target_species, config)
    
#     try:
#         # Load tree
#         tree = Tree(tree_path)
        
#         # Find target species in tree
#         found_species = {}
#         for species in target_species:
#             # Try different name formats
#             possible_names = [
#                 species,
#                 species.replace('_', ' '),
#                 ' '.join(species.split('_')[:2])  # Genus species format
#             ]
            
#             for name in possible_names:
#                 nodes = tree.search_nodes(name=name)
#                 if nodes:
#                     found_species[species] = nodes[0]
#                     break
        
#         if len(found_species) < 2:
#             logger.warning(f"Not enough species found in tree: {found_species.keys()}")
#             logger.info("Creating simple tree structure")
#             # FIXED: Call standalone function instead of self.method
#             return create_simple_tree.nope(target_species, config)
        
#         logger.info(f"Found {len(found_species)} species in tree: {list(found_species.keys())}")
        
#         # Get common ancestor and prune
#         common_ancestor = tree.get_common_ancestor(list(found_species.values()))
#         pruned_tree = common_ancestor.copy()
        
#         # Keep only target species
#         target_nodes = [found_species[sp] for sp in found_species.keys()]
#         pruned_tree.prune(target_nodes, preserve_branch_length=True)
        
#         return pruned_tree
        
#     except Exception as e:
#         logger.error(f"Tree pruning failed: {e}")
#         # FIXED: Call standalone function instead of self.method
#         return create_simple_tree.nope(target_species, config)

# def create_simple_tree(species_list: List[str], config: Dict) -> Optional['Tree']:
#     """Create simple tree structure if ete3 not available or pruning fails"""
    
#     logger = logging.getLogger(__name__)
    
#     # FIXED: Check ETE3 availability as import check instead of class variable
#     try:
#         from ete3 import Tree
#         ETE3_AVAILABLE = True
#     except ImportError:
#         ETE3_AVAILABLE = False
    
#     if not ETE3_AVAILABLE:
#         logger.info("Creating matplotlib-based tree visualisation")
#         return None
    
#     try:
#         # Create simple balanced tree
#         if len(species_list) == 4:
#             tree_str = f"(({species_list[0]}:0.1,{species_list[1]}:0.1):0.1,({species_list[2]}:0.1,{species_list[3]}:0.1):0.1);"
#         elif len(species_list) == 3:
#             tree_str = f"(({species_list[0]}:0.1,{species_list[1]}:0.1):0.1,{species_list[2]}:0.1);"
#         elif len(species_list) == 2:
#             tree_str = f"({species_list[0]}:0.1,{species_list[1]}:0.1);"
#         else:
#             # For more species, create a simple star topology
#             species_with_branch = [f"{sp}:0.1" for sp in species_list]
#             tree_str = f"({','.join(species_with_branch)});"
        
#         tree = Tree(tree_str)
#         logger.info(f"Created simple tree for {len(species_list)} species")
#         return tree
        
#     except Exception as e:
#         logger.error(f"Simple tree creation failed: {e}")
#         return None

# def _calculate_inversion_annotations(all_results: Dict, species_stats: Dict, config: Dict) -> Dict:
#     """Calculate inversion metrics for tree annotation"""
    
#     # FIXED: Converted from class method to standalone function - removed 'self' parameter and added 'config' parameter
#     from typing import Dict
    
#     inversion_data = {}
    
#     for species_name, stats in species_stats.items():
#         # FIXED: Handle different possible data structures for genome_size
#         try:
#             genome_size = stats['quality']['metrics'].get('total_length', 1000000)
#         except (KeyError, TypeError):
#             # Fallback if stats structure is different
#             genome_size = stats.get('genome_size', 1000000)
        
#         # Collect all inversions involving this species
#         species_inversions = []
#         for pair_name, pair_data in all_results.items():
#             if 'full_results' in pair_data and species_name in pair_data['species_pair']:
#                 inv_df = pair_data['inversion_df']
#                 if not inv_df.empty:
#                     species_inversions.extend(inv_df.to_dict('records'))
        
#         # Calculate metrics
#         total_inversions = len(species_inversions)
#         rate_per_mb = total_inversions / (genome_size / 1_000_000) if genome_size > 0 else 0
        
#         # FIXED: Handle different possible data structures for quality_score
#         try:
#             quality_score = stats['quality']['quality_score']
#         except (KeyError, TypeError):
#             quality_score = stats.get('quality_score', 1.0)
        
#         inversion_data[species_name] = {
#             'total_inversions': total_inversions,
#             'rate_per_mb': rate_per_mb,
#             'genome_size': genome_size,
#             'quality_score': quality_score
#         }
    
#     # Calculate normalised scores
#     if inversion_data:
#         max_rate = max(data['rate_per_mb'] for data in inversion_data.values())
#         for species_name in inversion_data:
#             if max_rate > 0:
#                 inversion_data[species_name]['normalised_score'] = inversion_data[species_name]['rate_per_mb'] / max_rate
#             else:
#                 inversion_data[species_name]['normalised_score'] = 0
    
#     return inversion_data

# def _annotate_tree_with_inversions(tree: 'Tree', inversion_data: Dict, config: Dict) -> 'Tree':
#     """Annotate tree nodes with inversion data"""
    
#     # FIXED: Converted from class method to standalone function - removed 'self' parameter and added 'config' parameter
#     import logging
#     import numpy as np
#     from typing import Dict
    
#     logger = logging.getLogger(__name__)
    
#     # FIXED: Check ETE3 availability as import check
#     try:
#         from ete3 import Tree
#         ETE3_AVAILABLE = True
#     except ImportError:
#         ETE3_AVAILABLE = False
    
#     if not ETE3_AVAILABLE or not tree:
#         return tree
    
#     try:
#         # Annotate leaf nodes
#         for leaf in tree.get_leaves():
#             species_name = leaf.name
#             if species_name in inversion_data:
#                 data = inversion_data[species_name]
#                 leaf.add_features(
#                     inversion_count=data['total_inversions'],
#                     inversion_rate=data['rate_per_mb'],
#                     normalised_score=data['normalised_score'],
#                     genome_size=data['genome_size']
#                 )
        
#         # Annotate internal nodes (ancestral state reconstruction)
#         for node in tree.traverse():
#             if not node.is_leaf():
#                 # Calculate average metrics for this clade
#                 leaves = node.get_leaves()
#                 if leaves:
#                     avg_rate = np.mean([getattr(leaf, 'inversion_rate', 0) for leaf in leaves])
#                     total_count = sum([getattr(leaf, 'inversion_count', 0) for leaf in leaves])
#                     avg_normalised = np.mean([getattr(leaf, 'normalised_score', 0) for leaf in leaves])
                    
#                     node.add_features(
#                         inversion_count=total_count,
#                         inversion_rate=avg_rate,
#                         normalised_score=avg_normalised
#                     )
        
#         return tree
        
#     except Exception as e:
#         logger.error(f"Tree annotation failed: {e}")
#         return tree

# def create_tree_plot(tree: 'Tree', output_file: Path, inversion_data: Dict, config: Dict):
#     """Create publication-quality tree plot with inversion annotations"""
    
#     # FIXED: Converted from class method to standalone function - removed 'self' parameter and added 'config' parameter
#     import logging
#     from pathlib import Path
#     from typing import Dict
    
#     logger = logging.getLogger(__name__)
    
#     try:
#         # FIXED: Check ETE3 availability
#         try:
#             from ete3 import Tree
#             ETE3_AVAILABLE = True
#         except ImportError:
#             ETE3_AVAILABLE = False
            
#         if ETE3_AVAILABLE and tree:
#             # Use ete3 for tree plotting
#             # FIXED: Call standalone function instead of self.method
#             create_ete3_tree_plot(tree, output_file, inversion_data, config)
#         else:
#             # Fallback to matplotlib
#             # FIXED: Call standalone function instead of self.method
#             create_matplotlib_tree_plot.nope(inversion_data, output_file, config)
            
#     except Exception as e:
#         logger.error(f"Tree plot creation failed: {e}")
#         # Create simple fallback
#         # FIXED: Call standalone function instead of self.method
#         create_matplotlib_tree_plot.nope(inversion_data, output_file, config)

# def create_ete3_tree_plot(tree: 'Tree', output_file: Path, inversion_data: Dict, config: Dict):
#     """Create tree plot using ete3 with better error handling"""
    
#     # FIXED: Converted from class method to standalone function - removed 'self' parameter and added 'config' parameter
#     import logging
#     from pathlib import Path
#     from typing import Dict
    
#     logger = logging.getLogger(__name__)
    
#     try:
#         from ete3 import TreeStyle, NodeStyle, AttrFace
#     except ImportError as e:
#         logger.error(f"Cannot import ete3 components: {e}")
#         return False

#     try:
#         # Create tree style
#         ts = TreeStyle()
#         ts.show_leaf_name = True
#         ts.show_branch_length = True
#         ts.show_branch_support = True
#         ts.mode = "r"  # rectangular
#         ts.scale = 200
        
#         # Configure tree style
#         ts.branch_vertical_margin = 10
#         ts.scale = 120
        
#         # Color nodes based on inversion rate
#         for node in tree.traverse():
#             ns = NodeStyle()
            
#             if node.is_leaf():
#                 # Color leaf nodes by inversion rate
#                 if hasattr(node, 'normalised_score'):
#                     score = node.normalised_score
#                     color_intensity = int(255 * (1 - score))  # Higher score = redder
#                     ns.bgcolor = f"rgb(255,{color_intensity},{color_intensity})"
#                     ns.size = 10
#                 else:
#                     ns.bgcolor = "lightblue"
#                     ns.size = 8
                
#                 # Add inversion count as text
#                 if hasattr(node, 'inversion_count'):
#                     count_face = AttrFace("inversion_count", fsize=8)
#                     node.add_face(count_face, column=1, position="branch-right")
#             else:
#                 # Internal nodes
#                 ns.size = 5
#                 ns.shape = "square"
#                 if hasattr(node, 'normalised_score'):
#                     score = node.normalised_score
#                     color_intensity = int(255 * (1 - score))
#                     ns.bgcolor = f"rgb(255,{color_intensity},{color_intensity})"
#                 else:
#                     ns.bgcolor = "gray"
            
#             node.set_style(ns)
        
#         # Render tree
#         tree.render(str(output_file), tree_style=ts, dpi=300)
#         logger.info(f"    ✅ ETE3 tree plot created: {output_file}")
#         return True
        
#     except Exception as e:
#         logger.error(f"ETE3 tree plot creation failed: {e}")
#         return False
    
# def create_matplotlib_tree_plot(inversion_data: Dict, output_file: Path, config: Dict):
#     """Create simple tree plot using matplotlib with proper linkage"""
    
#     # FIXED: Converted from class method to standalone function - removed 'self' parameter and added 'config' parameter
#     import matplotlib.pyplot as plt
#     import matplotlib.cm as cm
#     from matplotlib.colors import normalise
#     import numpy as np
#     import logging
#     from pathlib import Path
#     from typing import Dict
    
#     logger = logging.getLogger(__name__)
    
#     try:
#         from scipy.cluster.hierarchy import dendrogram, linkage
#         from scipy.spatial.distance import pdist
        
#         fig, ax = plt.subplots(figsize=(12, 8))
        
#         species_names = list(inversion_data.keys())
#         n_species = len(species_names)
        
#         if n_species < 2:
#             ax.text(0.5, 0.5, "Not enough species for tree", ha='center', va='center')
#             plt.savefig(output_file, dpi=300, bbox_inches='tight')
#             plt.close()
#             return
        
#         # Create distance matrix based on inversion rates
#         rates = [inversion_data[sp]['rate_per_mb'] for sp in species_names]
        
#         if n_species == 2:
#             # Special case for 2 species
#             distance = abs(rates[0] - rates[1]) if rates[0] != rates[1] else 0.1
#             linkage_matrix = np.array([[0, 1, distance, 2]])
#         else:
#             # For 3+ species, create proper distance matrix
#             rate_matrix = np.array(rates).reshape(-1, 1)
#             distances = pdist(rate_matrix, metric='euclidean')
#             linkage_matrix = linkage(distances, method='average')  # Use 'average' instead of 'ward'
        
#         # Create dendrogram
#         dendro = dendrogram(linkage_matrix, labels=species_names, ax=ax, orientation='top')
        
#         # Color-code by inversion rates
#         max_rate = max(rates) if rates else 1
#         for i, species in enumerate(species_names):
#             rate = inversion_data[species]['rate_per_mb']
#             normalised_rate = rate / max_rate if max_rate > 0 else 0
#             color = cm.get_cmap('Reds')(normalised_rate)
            
#             # Find the corresponding x position in the dendrogram
#             if i < len(dendro['icoord']):
#                 x_pos = (dendro['icoord'][i][1] + dendro['icoord'][i][2]) / 2
#                 ax.text(x_pos, -0.1, f"{rate:.2f}", 
#                     rotation=45, ha='right', va='top', color=color, fontweight='bold')
        
#         ax.set_title('Bibionidae Phylogeny with Inversion Rates', fontsize=14, fontweight='bold')
#         ax.set_ylabel('Distance')
#         ax.set_xlabel('Species')
        
#         # Add colorbar
#         sm = cm.ScalarMappable(cmap='Reds', norm=normalise(vmin=0, vmax=max_rate))
#         sm.set_array([])
#         cbar = plt.colorbar(sm, ax=ax, shrink=0.6)
#         cbar.set_label('Inversion Rate (per Mb)', rotation=270, labelpad=15)
        
#         plt.tight_layout()
#         plt.savefig(output_file, dpi=300, bbox_inches='tight')
#         plt.close()
        
#         logger.info(f"    ✅ Matplotlib tree plot created: {output_file}")
#         return True
        
#     except Exception as e:
#         logger.error(f"Matplotlib tree plot failed: {e}")
#         # Create simple bar plot fallback
#         fig, ax = plt.subplots(figsize=(10, 6))
#         species_names = list(inversion_data.keys())
#         rates = [inversion_data[sp]['rate_per_mb'] for sp in species_names]
        
#         max_rate = max(rates) if rates else 1
#         colors = [cm.get_cmap('Reds')(r/max_rate) for r in rates]
#         bars = ax.bar(species_names, rates, color=colors)
        
#         ax.set_title('Inversion Rates by Species', fontsize=14, fontweight='bold')
#         ax.set_ylabel('Inversions per Mb')
#         ax.set_xlabel('Species')
#         plt.xticks(rotation=45)
#         plt.tight_layout()
#         plt.savefig(output_file, dpi=300, bbox_inches='tight')
#         plt.close()
#         return True

# def create_tree_heatmap(inversion_data: Dict, output_file: Path, config: Dict):
#     """Create heatmap showing inversion metrics across species"""
    
#     # FIXED: Converted from class method to standalone function - removed 'self' parameter and added 'config' parameter
#     import matplotlib.pyplot as plt
#     import pandas as pd
#     import numpy as np
#     import logging
#     from pathlib import Path
#     from typing import Dict
    
#     logger = logging.getLogger(__name__)
    
#     try:
#         # Try to import seaborn, fallback to matplotlib if not available
#         try:
#             import seaborn as sns
#             SEABORN_AVAILABLE = True
#         except ImportError:
#             SEABORN_AVAILABLE = False
#             logger.warning("Seaborn not available, using matplotlib heatmap")
        
#         # Prepare data for heatmap
#         species_names = list(inversion_data.keys())
#         metrics = ['total_inversions', 'rate_per_mb', 'normalised_score', 'genome_size']
        
#         heatmap_data = []
#         for metric in metrics:
#             row = [inversion_data[sp].get(metric, 0) for sp in species_names]
#             # normalise each metric to 0-1 scale
#             max_val = max(row) if max(row) > 0 else 1
#             normalised_row = [val / max_val for val in row]
#             heatmap_data.append(normalised_row)
        
#         heatmap_df = pd.DataFrame(heatmap_data, 
#                                 index=['Total Inversions', 'Rate per Mb', 'normalised Score', 'Genome Size'],
#                                 columns=species_names)
        
#         # Create heatmap
#         fig, ax = plt.subplots(figsize=(10, 6))
        
#         if SEABORN_AVAILABLE:
#             # Use seaborn if available
#             sns.heatmap(heatmap_df, annot=True, fmt='.2f', cmap='RdYlBu_r', 
#                         cbar_kws={'label': 'normalised Value'}, ax=ax)
#         else:
#             # FIXED: Fallback matplotlib heatmap if seaborn not available
#             im = ax.imshow(heatmap_data, cmap='RdYlBu_r', aspect='auto')
            
#             # Add colorbar
#             cbar = plt.colorbar(im, ax=ax)
#             cbar.set_label('normalised Value')
            
#             # Add annotations
#             for i in range(len(metrics)):
#                 for j in range(len(species_names)):
#                     text = ax.text(j, i, f'{heatmap_data[i][j]:.2f}',
#                                  ha="center", va="center", color="black")
            
#             # Set ticks and labels
#             ax.set_xticks(range(len(species_names)))
#             ax.set_xticklabels(species_names)
#             ax.set_yticks(range(len(metrics)))
#             ax.set_yticklabels(['Total Inversions', 'Rate per Mb', 'normalised Score', 'Genome Size'])
        
#         ax.set_title('Inversion Metrics Heatmap', fontsize=14, fontweight='bold')
#         ax.set_xlabel('Species')
#         ax.set_ylabel('Metrics')
        
#         plt.tight_layout()
#         plt.savefig(output_file, dpi=300, bbox_inches='tight')
#         plt.close()
        
#         logger.info(f"    ✅ Tree heatmap created: {output_file}")
        
#     except Exception as e:
#         logger.error(f"Tree heatmap creation failed: {e}")
    
# def create_busco_phylogenetic_tree(all_results: Dict, species_stats: Dict, output_dir: Path, config: Dict) -> Dict:
#     """Create phylogenetic tree from existing Newick file with configurable node annotations"""
    
#     # FIXED: Converted from class method to standalone function - removed 'self' parameter and added 'config' parameter
#     import logging
#     from pathlib import Path
#     from typing import Dict
    
#     logger = logging.getLogger(__name__)
#     logger.info("  🌳 Creating BUSCO phylogenetic tree from existing file...")
    
#     # FIXED: Use config parameter instead of self.config
#     tree_config = config.get('tree_annotation', {})
#     source_tree = tree_config.get('source_tree_path', 'diptera_clean_20species.newick')
    
#     if not Path(source_tree).exists():
#         logger.warning(f"Source tree file not found: {source_tree}")
#         return {}
    
#     try:
#         from ete3 import Tree
        
#         # Load the tree
#         tree = Tree(source_tree)
#         logger.info(f"    • Loaded tree with {len(tree)} nodes")
        
#         # Get species names
#         species_names = list(species_stats.keys())
#         logger.info(f"    • Target species: {species_names}")
        
#         # Prune tree to target species (if enabled)
#         if tree_config.get('prune_to_target_species', True):
#             # FIXED: Call standalone function instead of self.method
#             tree = prune_tree_to_species(tree, species_names, config)
#             if tree is None:
#                 logger.warning("Tree pruning failed")
#                 return {}
        
#         # Calculate metrics for node annotation
#         # FIXED: Call standalone function instead of self.method
#         node_metrics = calculate_tree_node_metrics(all_results, species_stats, config)
        
#         # Annotate tree with metrics
#         # FIXED: Call standalone function instead of self.method
#         annotated_tree = annotate_tree_nodes(tree, node_metrics, config)
        
#         # Create tree plots
#         tree_dir = output_dir / 'busco_phylogenetic_trees'
#         tree_dir.mkdir(exist_ok=True)
        
#         results = {}
        
#         # Save annotated tree
#         tree_file = tree_dir / 'busco_phylogenetic_tree.newick'
#         annotated_tree.write(format=1, outfile=str(tree_file))
#         results['newick'] = tree_file
        
#         # Create tree plot
#         plot_file = tree_dir / 'busco_phylogenetic_tree.png'
#         # FIXED: Call standalone function instead of self.method
#         create_annotated_tree_plot(annotated_tree, plot_file, node_metrics, config)
#         results['plot'] = plot_file
        
#         logger.info(f"    ✅ BUSCO phylogenetic tree created: {len(results)} outputs")
#         return results
        
#     except ImportError:
#         logger.error("ete3 required for phylogenetic tree manipulation")
#         return {}
#     except Exception as e:
#         logger.error(f"BUSCO phylogenetic tree creation failed: {e}")
#         return {}

# def calculate_tree_node_metrics(all_results: Dict, species_stats: Dict, config: Dict) -> Dict:
#     """Calculate metrics for tree node annotation from config"""
    
#     # FIXED: Converted from class method to standalone function - removed 'self' parameter and added 'config' parameter
#     from typing import Dict
    
#     # FIXED: Use config parameter instead of self.config
#     tree_config = config.get('tree_annotation', {})
#     metrics_config = tree_config.get('annotation_metrics', {})
    
#     node_metrics = {}
    
#     for species_name, stats in species_stats.items():
#         metrics = {}
        
#         # Default metric: inversion rate per MB
#         if metrics_config.get('inversion_rate_per_mb', True):
#             # FIXED: Handle different possible data structures for genome_size
#             try:
#                 genome_size = stats['quality']['metrics'].get('total_length', 1000000)
#             except (KeyError, TypeError):
#                 genome_size = stats.get('genome_size', 1000000)
                
#             species_inversions = sum(1 for pair_results in all_results.values() 
#                                 if 'full_results' in pair_results and 
#                                 species_name in pair_results['species_pair'] and
#                                 len(pair_results['inversion_df']) > 0)
#             metrics['inversion_rate_per_mb'] = species_inversions / (genome_size / 1_000_000) if genome_size > 0 else 0
        
#         # Inversion count
#         if metrics_config.get('inversion_count', True):
#             species_inversions = []
#             for pair_results in all_results.values():
#                 if 'full_results' in pair_results and species_name in pair_results['species_pair']:
#                     species_inversions.extend(pair_results['inversion_df'].to_dict('records'))
#             metrics['inversion_count'] = len(species_inversions)
        
#         # normalised inversion score
#         if metrics_config.get('normalised_inversion_score', True):
#             if 'inversion_rate_per_mb' in metrics:
#                 all_rates = [node_metrics.get(sp, {}).get('inversion_rate_per_mb', 0) 
#                         for sp in species_stats.keys()]
#                 max_rate = max(all_rates + [metrics['inversion_rate_per_mb']])
#                 metrics['normalised_score'] = metrics['inversion_rate_per_mb'] / max_rate if max_rate > 0 else 0
        
#         # Assembly quality
#         if metrics_config.get('assembly_quality', False):
#             # FIXED: Handle different possible data structures for quality_score
#             try:
#                 metrics['quality_score'] = stats['quality'].get('quality_score', 0.5)
#             except (KeyError, TypeError):
#                 metrics['quality_score'] = stats.get('quality_score', 0.5)
        
#         # Genome size
#         if metrics_config.get('genome_size', False):
#             # FIXED: Handle different possible data structures for genome size
#             try:
#                 metrics['genome_size_mb'] = stats['quality']['metrics'].get('total_length', 0) / 1_000_000
#             except (KeyError, TypeError):
#                 metrics['genome_size_mb'] = stats.get('genome_size', 0) / 1_000_000
        
#         node_metrics[species_name] = metrics
    
#     return node_metrics


# def prune_tree_to_species(tree: 'Tree', target_species: List[str], config: Dict) -> Optional['Tree']:
#     """Prune tree to target species"""
    
#     # FIXED: Converted from class method to standalone function - removed 'self' parameter and added 'config' parameter
#     import logging
#     from typing import List, Optional
    
#     logger = logging.getLogger(__name__)
    
#     try:
#         # Get all leaf names
#         leaf_names = [leaf.name for leaf in tree.get_leaves()]
        
#         # Find matches (allowing partial matches)
#         species_in_tree = []
#         for species in target_species:
#             matches = [leaf for leaf in leaf_names if species in leaf or leaf in species]
#             if matches:
#                 species_in_tree.append(matches[0])  # Take first match
#             else:
#                 logger.warning(f"Species {species} not found in tree")
        
#         if len(species_in_tree) < 2:
#             logger.error(f"Too few species found in tree: {species_in_tree}")
#             return None
        
#         # Prune tree
#         tree.prune(species_in_tree)
#         logger.info(f"    • Pruned tree to {len(species_in_tree)} species: {species_in_tree}")
        
#         return tree
        
#     except Exception as e:
#         logger.error(f"Tree pruning failed: {e}")
#         return None
    

# def annotate_tree_nodes(tree: 'Tree', node_metrics: Dict, config: Dict) -> 'Tree':
#     """Annotate tree nodes with calculated metrics"""
    
#     # FIXED: Converted from class method to standalone function - removed 'self' parameter and added 'config' parameter
#     from typing import Dict
    
#     for leaf in tree.get_leaves():
#         leaf_name = leaf.name
        
#         # Find matching species (allowing partial matches)
#         matching_species = None
#         for species in node_metrics.keys():
#             if species in leaf_name or leaf_name in species:
#                 matching_species = species
#                 break
        
#         if matching_species and matching_species in node_metrics:
#             metrics = node_metrics[matching_species]
            
#             # Add metrics as node features
#             for metric_name, metric_value in metrics.items():
#                 setattr(leaf, metric_name, metric_value)
    
#     return tree

# def create_annotated_tree_plot(tree, output_file, node_metrics, config: Dict):
#     """Create tree plot using ete3 render without GUI components"""
    
#     # FIXED: Converted from class method to standalone function - removed 'self' parameter and added 'config' parameter
#     import logging
#     from pathlib import Path
#     from typing import Dict
    
#     logger = logging.getLogger(__name__)
    
#     try:
#         # Use basic ete3 render - NO GUI components
#         tree.render(str(output_file), w=800, h=600, dpi=300)
#         logger.info(f" • BUSCO tree rendered: {output_file}")
        
#     except Exception as e:
#         logger.error(f"ete3 render failed: {e}, falling back to matplotlib")
#         # Fallback to matplotlib
#         import matplotlib.pyplot as plt
#         import numpy as np
        
#         fig, ax = plt.subplots(1, 1, figsize=(12, 8))
        
#         leaves = tree.get_leaves()
#         species_names = [leaf.name.replace('_', ' ') for leaf in leaves]
        
#         for i, name in enumerate(species_names):
#             ax.text(1.0, i, name, ha='left', va='center', fontsize=11, weight='bold')
#             ax.plot([0.8, 1.0], [i, i], 'k-', linewidth=2)
            
#             # Add inversion count if available
#             inversions = getattr(leaves[i], 'inversions', 0)
#             ax.text(1.5, i, f"({inversions} inv)", ha='left', va='center', 
#                 fontsize=9, color='red')
        
#         ax.set_xlim(0, 2.0)
#         ax.set_ylim(-0.5, len(leaves) - 0.5)
#         ax.set_title('BUSCO Phylogenetic Tree', fontsize=14, weight='bold')
#         ax.axis('off')
        
#         plt.tight_layout()
#         plt.savefig(output_file, dpi=300, bbox_inches='tight')
#         plt.close()
        
#         logger.info(f" • Matplotlib fallback tree saved: {output_file}")
        
                
# def create_circular_synteny_plot(self,
#                             ortholog_df: pd.DataFrame,
#                             inversion_df: pd.DataFrame,
#                             species1: str,
#                             species2: str) -> Path:
#     """
#     Create circular synteny plot (Circos-style)
    
#     Args:
#         ortholog_df: Ortholog pairs DataFrame
#         inversion_df: Inversion events DataFrame
#         species1: First species name
#         species2: Second species name
        
#     Returns:
#         Path to generated circular plot
#     """
#     logger.info(f"Creating circular synteny plot for {species1} vs {species2}")
    
#     fig, ax = plt.subplots(figsize=(12, 12), subplot_kw=dict(projection='polar'))
    
#     # Get chromosome lists
#     chr1_list = sorted(ortholog_df['chr1'].unique())
#     chr2_list = sorted(ortholog_df['chr2'].unique())
    
#     # Create angular positions for chromosomes
#     n_chr1 = len(chr1_list)
#     n_chr2 = len(chr2_list)
    
#     # First species on top half, second species on bottom half
#     chr1_angles = np.linspace(0, np.pi, n_chr1)
#     chr2_angles = np.linspace(np.pi, 2*np.pi, n_chr2)
    
#     # Plot chromosome arcs
#     radius = 1.0
#     for i, chr_name in enumerate(chr1_list):
#         ax.plot([chr1_angles[i], chr1_angles[i]], [0.8, radius], 'k-', linewidth=3)
#         ax.text(chr1_angles[i], radius + 0.1, chr_name, ha='center', va='center')
    
#     for i, chr_name in enumerate(chr2_list):
#         ax.plot([chr2_angles[i], chr2_angles[i]], [0.8, radius], 'k-', linewidth=3)
#         ax.text(chr2_angles[i], radius + 0.1, chr_name, ha='center', va='center')
    
#     # Plot synteny connections
#     for _, ortholog in ortholog_df.iterrows():
#         chr1_idx = chr1_list.index(ortholog['chr1'])
#         chr2_idx = chr2_list.index(ortholog['chr2'])
        
#         angle1 = chr1_angles[chr1_idx]
#         angle2 = chr2_angles[chr2_idx]
        
#         # Determine color based on strand consistency
#         same_strand = ortholog['strand1'] == ortholog['strand2']
#         color = CONFIG['syri_colors']['syn'] if same_strand else CONFIG['syri_colors']['inv']
#         alpha = 0.1
        
#         # Draw connection arc
#         angles = np.linspace(angle1, angle2, 100)
#         radii = 0.7 * np.sin(np.pi * (angles - angle1) / (angle2 - angle1))
        
#         ax.plot(angles, radii, color=color, alpha=alpha, linewidth=0.5)
    
#     # Highlight inversions
#     for _, inversion in inversion_df.iterrows():
#         try:
#             chr1_idx = chr1_list.index(inversion['chr1'])
#             chr2_idx = chr2_list.index(inversion['chr2'])
            
#             angle1 = chr1_angles[chr1_idx]
#             angle2 = chr2_angles[chr2_idx]
            
#             # Draw thick inversion connection
#             angles = np.linspace(angle1, angle2, 100)
#             radii = 0.6 * np.sin(np.pi * (angles - angle1) / (angle2 - angle1))
            
#             ax.plot(angles, radii, color=CONFIG['syri_colors']['inv'], linewidth=2, alpha=0.8)
#         except (ValueError, KeyError):
#             continue
    
#     # Customize circular plot
#     ax.set_ylim(0, 1.2)
#     ax.set_title(f'Circular Synteny: {species1} vs {species2}', fontsize=14, fontweight='bold', pad=20)
#     ax.grid(False)
#     ax.set_rticks([])
#     ax.set_thetagrids([])
    
#     # Save plot
#     plot_file = self.output_dir / f'circular_synteny_{species1}_vs_{species2}.png'
#     plt.tight_layout()
#     plt.savefig(plot_file, dpi=300, bbox_inches='tight')
#     plt.close()
    
#     logger.info(f"Created circular synteny plot: {plot_file}")
#     return plot_file
    
# def create_chromosome_comparison_plot(self,
#                                     ortholog_df: pd.DataFrame,
#                                     inversion_df: pd.DataFrame,
#                                     species1: str,
#                                     species2: str) -> Path:
#     """
#     Create detailed chromosome-by-chromosome comparison plot
    
#     Args:
#         ortholog_df: Ortholog pairs DataFrame
#         inversion_df: Inversion events DataFrame
#         species1: First species name
#         species2: Second species name
        
#     Returns:
#         Path to generated plot
#     """
#     logger.info(f"Creating chromosome comparison plot for {species1} vs {species2}")
    
#     # Get unique chromosome pairs
#     chr_pairs = ortholog_df.groupby(['chr1', 'chr2']).size().reset_index(name='count')
#     chr_pairs = chr_pairs.sort_values('count', ascending=False)
    
#     n_pairs = min(len(chr_pairs), 9)  # Show top 9 pairs
    
#     fig, axes = plt.subplots(3, 3, figsize=(15, 15))
#     axes = axes.flatten()
    
#     for i, (_, pair) in enumerate(chr_pairs.head(n_pairs).iterrows()):
#         ax = axes[i]
        
#         chr1 = pair['chr1']
#         chr2 = pair['chr2']
        
#         # Get orthologs for this chromosome pair
#         pair_orthologs = ortholog_df[
#             (ortholog_df['chr1'] == chr1) & 
#             (ortholog_df['chr2'] == chr2)
#         ]
        
#         # Get inversions for this pair
#         pair_inversions = inversion_df[
#             (inversion_df['chr1'] == chr1) & 
#             (inversion_df['chr2'] == chr2)
#         ]
        
#         # Create dot plot for this chromosome pair
#         self.create_chromosome_dotplot(ax, pair_orthologs, pair_inversions, chr1, chr2)
        
#         ax.set_title(f'{chr1} vs {chr2}\n({len(pair_orthologs)} orthologs)', fontsize=10)
    
#     # Hide unused subplots
#     for i in range(n_pairs, len(axes)):
#         axes[i].set_visible(False)
    
#     plt.suptitle(f'Chromosome Comparisons: {species1} vs {species2}', fontsize=14, fontweight='bold')
#     plt.tight_layout()
    
#     # Save plot
#     plot_file = self.output_dir / f'chromosome_comparison_{species1}_vs_{species2}.png'
#     plt.savefig(plot_file, dpi=300, bbox_inches='tight')
#     plt.close()
    
#     logger.info(f"Created chromosome comparison plot: {plot_file}")
#     return plot_file

# # Helper methods
# def create_chromosome_positions(chr_list: List[str], 
#                                chromosome_lengths: Dict[str, int] = None,
#                                genome: str = 'first') -> Dict[str, float]:
#     """Create cumulative positions for chromosomes"""
#     positions = {}
#     current_pos = 0
    
#     for chr_name in chr_list:
#         positions[chr_name] = current_pos
        
#         # Use actual length if available, otherwise estimate
#         if chromosome_lengths and chr_name in chromosome_lengths:
#             chr_length = chromosome_lengths[chr_name] / 1_000_000  # Convert to Mb
#         else:
#             chr_length = 10  # Default 10 Mb spacing
        
#         current_pos += chr_length + 2  # Add 2 Mb gap between chromosomes
    
#     return positions

# def plot_chromosome_backbones(ax, chr1_positions: Dict, chr2_positions: Dict,
#                              species1: str, species2: str):
#     """Plot chromosome backbone representations"""
#     y1_base = len(chr1_positions) + 2
#     y2_base = 1
    
#     # Plot species 1 chromosomes (top)
#     for chr_name, x_pos in chr1_positions.items():
#         ax.plot([x_pos, x_pos + 8], [y1_base, y1_base], 'k-', linewidth=4)
#         ax.text(x_pos + 4, y1_base + 0.3, chr_name, ha='center', va='bottom', fontsize=8)
    
#     # Plot species 2 chromosomes (bottom)
#     for chr_name, x_pos in chr2_positions.items():
#         ax.plot([x_pos, x_pos + 8], [y2_base, y2_base], 'k-', linewidth=4)
#         ax.text(x_pos + 4, y2_base - 0.3, chr_name, ha='center', va='top', fontsize=8)
    
#     # Add species labels
#     ax.text(-2, y1_base, species1, ha='right', va='center', fontsize=12, fontweight='bold')
#     ax.text(-2, y2_base, species2, ha='right', va='center', fontsize=12, fontweight='bold')

# def plot_synteny_blocks(ax, synteny_df: pd.DataFrame, ortholog_df: pd.DataFrame,
#                        chr1_positions: Dict, chr2_positions: Dict):
#     """Plot synteny blocks as colored regions"""
#     y1_base = len(chr1_positions) + 2
#     y2_base = 1
    
#     for _, synteny in synteny_df.iterrows():
#         chr1 = synteny.get('chr1', '')
#         chr2 = synteny.get('chr2', '')
        
#         if chr1 in chr1_positions and chr2 in chr2_positions:
#             # Get block boundaries from orthologs
#             block_orthologs = ortholog_df[
#                 (ortholog_df['chr1'] == chr1) & 
#                 (ortholog_df['chr2'] == chr2)
#             ]
            
#             if len(block_orthologs) > 0:
#                 x1_start = chr1_positions[chr1] + block_orthologs['start1'].min() / 1_000_000
#                 x1_end = chr1_positions[chr1] + block_orthologs['end1'].max() / 1_000_000
#                 x2_start = chr2_positions[chr2] + block_orthologs['start2'].min() / 1_000_000
#                 x2_end = chr2_positions[chr2] + block_orthologs['end2'].max() / 1_000_000
                
#                 # Determine block color based on synteny type
#                 synteny_type = synteny.get('synteny_type', 'colinear')
#                 color = CONFIG['syri_colors']['inv'] if synteny_type == 'inverted' else CONFIG['syri_colors']['syn']
                
#                 # Draw synteny block regions
#                 rect1 = patches.Rectangle((x1_start, y1_base - 0.1), x1_end - x1_start, 0.2, 
#                                         facecolor=color, alpha=0.6, edgecolor='none')
#                 rect2 = patches.Rectangle((x2_start, y2_base - 0.1), x2_end - x2_start, 0.2,
#                                         facecolor=color, alpha=0.6, edgecolor='none')
                
#                 ax.add_patch(rect1)
#                 ax.add_patch(rect2)

# def plot_inversions(ax, inversion_df: pd.DataFrame, 
#                    chr1_positions: Dict, chr2_positions: Dict):
#     """Plot specific inversion events"""
#     y1_base = len(chr1_positions) + 2
#     y2_base = 1
    
#     for _, inversion in inversion_df.iterrows():
#         chr1 = inversion.get('chr1', '')
#         chr2 = inversion.get('chr2', '')
        
#         if chr1 in chr1_positions and chr2 in chr2_positions:
#             x1 = chr1_positions[chr1] + inversion.get('start1', 0) / 1_000_000
#             x2 = chr2_positions[chr2] + inversion.get('start2', 0) / 1_000_000
            
#             # Draw inversion connection with distinctive style
#             ax.plot([x1, x2], [y1_base, y2_base], 
#                     color=CONFIG['syri_colors']['inv'], linewidth=2, alpha=0.8, linestyle='--')
            
#             # Add inversion markers
#             ax.scatter([x1], [y1_base], color=CONFIG['syri_colors']['inv'], s=50, marker='v', zorder=5)
#             ax.scatter([x2], [y2_base], color=CONFIG['syri_colors']['inv'], s=50, marker='^', zorder=5)

# def plot_ortholog_connections(ax, ortholog_df: pd.DataFrame,
#                              chr1_positions: Dict, chr2_positions: Dict):
#     """Plot individual ortholog connections"""
#     y1_base = len(chr1_positions) + 2
#     y2_base = 1
    
#     # Sample orthologs to avoid overcrowding
#     if len(ortholog_df) > 1000:
#         sample_orthologs = ortholog_df.sample(n=1000, random_state=42)
#     else:
#         sample_orthologs = ortholog_df
    
#     lines = []
#     colors = []
    
#     for _, ortholog in sample_orthologs.iterrows():
#         chr1 = ortholog['chr1']
#         chr2 = ortholog['chr2']
        
#         if chr1 in chr1_positions and chr2 in chr2_positions:
#             x1 = chr1_positions[chr1] + ortholog['start1'] / 1_000_000
#             x2 = chr2_positions[chr2] + ortholog['start2'] / 1_000_000
            
#             lines.append([(x1, y1_base), (x2, y2_base)])
            
#             # Color by strand consistency
#             same_strand = ortholog['strand1'] == ortholog['strand2']
#             color = CONFIG['syri_colors']['syn'] if same_strand else CONFIG['syri_colors']['inv']
#             colors.append(color)
    
#     # Plot all connections at once for efficiency
#     lc = LineCollection(lines, colors=colors, alpha=0.3, linewidths=0.5)
#     ax.add_collection(lc)

# def add_syri_legend(ax):
#     """Add SyRI-style legend to plot"""
#     legend_elements = [
#         plt.Line2D([0], [0], color=CONFIG['syri_colors']['syn'], lw=2, label='Syntenic'),
#         plt.Line2D([0], [0], color=CONFIG['syri_colors']['inv'], lw=2, label='Inverted'),
#         plt.Line2D([0], [0], color=CONFIG['syri_colors']['inv'], lw=2, linestyle='--', label='Inversion Event')
#     ]
    
#     ax.legend(handles=legend_elements, loc='upper right', frameon=True, 
#               fancybox=True, shadow=True)

# def create_chromosome_dotplot(ax, orthologs: pd.DataFrame, inversions: pd.DataFrame,
#                              chr1: str, chr2: str):
#     """Create dot plot for chromosome pair comparison"""
#     if len(orthologs) == 0:
#         ax.text(0.5, 0.5, 'No orthologs', ha='center', va='center', transform=ax.transAxes)
#         return
    
#     # Plot orthologs as dots
#     same_strand = orthologs['strand1'] == orthologs['strand2']
    
#     # Syntenic orthologs
#     syn_orthologs = orthologs[same_strand]
#     if len(syn_orthologs) > 0:
#         ax.scatter(syn_orthologs['start1'], syn_orthologs['start2'],
#                    c=CONFIG['syri_colors']['syn'], alpha=0.6, s=10, label='Syntenic')
    
#     # Inverted orthologs
#     inv_orthologs = orthologs[~same_strand]
#     if len(inv_orthologs) > 0:
#         ax.scatter(inv_orthologs['start1'], inv_orthologs['start2'],
#                    c=CONFIG['syri_colors']['inv'], alpha=0.6, s=10, label='Inverted')
    
#     # Highlight specific inversions
#     if len(inversions) > 0:
#         ax.scatter(inversions['start1'], inversions['start2'],
#                    c=CONFIG['syri_colors']['inv'], s=50, marker='D', 
#                    edgecolors='black', linewidth=0.5, label='Inversion Event')
    
#     ax.set_xlabel(f'{chr1} position')
#     ax.set_ylabel(f'{chr2} position')
#     ax.legend()
    
#     # Add diagonal line for perfect synteny
#     if len(orthologs) > 0:
#         min_pos = min(orthologs['start1'].min(), orthologs['start2'].min())
#         max_pos = max(orthologs['end1'].max(), orthologs['end2'].max())
#         ax.plot([min_pos, max_pos], [min_pos, max_pos], 'k--', alpha=0.3, linewidth=1)

# def integrate_with_syri(output_dir: Union[str, Path],
#                        pairwise_results: Dict[str, Any],
#                        species_pairs: List[Tuple[str, str]] = None) -> Dict[str, Any]:
#     """
#     Complete SyRI integration for multi-species analysis
    
#     Args:
#         output_dir: Output directory for SyRI files
#         pairwise_results: Dictionary of pairwise analysis results
#         species_pairs: Optional list of specific pairs to process
        
#     Returns:
#         Dictionary with SyRI integration results
#     """
#     logger.info("Starting complete SyRI integration")
    
#     # Create output directory
#     syri_dir = Path(output_dir) / 'syri_integration'
#     syri_dir.mkdir(parents=True, exist_ok=True)
    
#     # Process all species pairs or specified ones
#     if species_pairs is None:
#         species_pairs = []
#         for pair_key in pairwise_results.keys():
#             if '_vs_' in pair_key and 'error' not in pairwise_results[pair_key]:
#                 species1, species2 = pair_key.split('_vs_')
#                 species_pairs.append((species1, species2))
    
#     syri_results = {}
    
#     for species1, species2 in species_pairs:
#         try:
#             logger.info(f"Processing SyRI integration: {species1} vs {species2}")
            
#             pair_key = f"{species1}_vs_{species2}"
#             pair_data = pairwise_results.get(pair_key, {})
            
#             if 'error' in pair_data:
#                 logger.warning(f"Skipping {pair_key} due to analysis error")
#                 continue
            
#             # Extract data
#             ortholog_df = pair_data.get('ortholog_df', pd.DataFrame())
#             synteny_df = pair_data.get('synteny_df', pd.DataFrame())
#             inversion_df = pair_data.get('inversion_df', pd.DataFrame())
            
#             # Create visualisations using the standalone functions
#             fig, ax = plt.subplots(figsize=(15, 8))
            
#             # Get chromosome lists
#             chr1_list = sorted(ortholog_df['chr1'].unique()) if len(ortholog_df) > 0 else []
#             chr2_list = sorted(ortholog_df['chr2'].unique()) if len(ortholog_df) > 0 else []
            
#             if chr1_list and chr2_list:
#                 # Create chromosome positions
#                 chr1_positions = create_chromosome_positions(chr1_list)
#                 chr2_positions = create_chromosome_positions(chr2_list)
                
#                 # Plot components
#                 plot_chromosome_backbones(ax, chr1_positions, chr2_positions, species1, species2)
#                 plot_synteny_blocks(ax, synteny_df, ortholog_df, chr1_positions, chr2_positions)
#                 plot_inversions(ax, inversion_df, chr1_positions, chr2_positions)
#                 plot_ortholog_connections(ax, ortholog_df, chr1_positions, chr2_positions)
#                 add_syri_legend(ax)
                
#                 ax.set_title(f'Synteny Analysis: {species1} vs {species2}')
#                 ax.set_xlim(-5, max(max(chr1_positions.values(), default=0), 
#                                    max(chr2_positions.values(), default=0)) + 15)
#                 ax.set_ylim(0, len(chr1_positions) + 4)
#                 ax.axis('off')
            
#             # Save plot
#             syri_plot_path = syri_dir / f"{species1}_vs_{species2}_syri_style.png"
#             plt.savefig(syri_plot_path, dpi=300, bbox_inches='tight')
#             plt.close()
            
#             syri_results[pair_key] = {
#                 'syri_style_plot': str(syri_plot_path),
#                 'species_pair': (species1, species2)
#             }
            
#             logger.info(f"Completed SyRI integration: {pair_key}")
            
#         except Exception as e:
#             logger.error(f"SyRI integration failed for {species1} vs {species2}: {e}")
#             syri_results[f"{species1}_vs_{species2}"] = {'error': str(e)}
    
#     logger.info(f"SyRI integration completed for {len(species_pairs)} species pairs")
#     return syri_results