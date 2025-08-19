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
    
    # Create data directory for TSV outputs
    data_dir = plots_dir.parent / 'data'
    data_dir.mkdir(parents=True, exist_ok=True)
    
    # Get chromosome sizes and linearise coordinates
    genome_data = get_chromosome_sizes_from_sheet(CONFIG)
    
    # STAGE 11: Save genome data from sheet
    genome_data_file = data_dir / 'stage_11_genome_data_from_sheet.tsv'
    genome_data.to_csv(genome_data_file, sep='\t', index=False)
    logger.info(f"Saved genome data from sheet: {genome_data_file}")
    
    joined_df_linear, species1_data, species2_data = linearise_coordinates(joined_df, genome_data, plot_config)
    
    # Create the plot
    fig, ax = plt.subplots(figsize=plot_config['figure_size'])
    
    # Split data by inversion status
    syntenic_data = joined_df_linear[joined_df_linear['is_flipped'] == False]
    inverted_data = joined_df_linear[joined_df_linear['is_flipped'] == True]
    
    # STAGE 12: Save linearized plot data separated by inversion status
    syntenic_file = data_dir / 'stage_12_linearized_syntenic_data.tsv'
    syntenic_data.to_csv(syntenic_file, sep='\t', index=False)
    logger.info(f"Saved linearized syntenic data: {syntenic_file}")
    
    inverted_file = data_dir / 'stage_12_linearized_inverted_data.tsv'
    inverted_data.to_csv(inverted_file, sep='\t', index=False)
    logger.info(f"Saved linearized inverted data: {inverted_file}")
    
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
    add_chromosome_boundaries(ax, species1_data, species2_data, plot_config)
    
    add_chromosome_labels(ax, species1_data, species2_data, joined_df_linear, plot_config)
        
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



def add_chromosome_labels(ax, species1_data, species2_data, joined_df_linear, config=None):
    """Add chromosome labels on secondary axes"""
    
    # Create data directory for TSV outputs
    if config:
        output_dir = config.get('output_dir', '.')
        data_dir = Path(output_dir) / 'data'
        data_dir.mkdir(parents=True, exist_ok=True)
    
    species1_data = species1_data.sort_values('chromosome')
    chr1_positions = []
    chr1_labels = []
    cumulative_offset = 0
    
    chr1_label_data = []
    for _, row in species1_data.iterrows():
        chr_size = row['chromsome_size_b']
        midpoint = cumulative_offset + (chr_size / 2)
        chr1_positions.append(midpoint)
        chr1_labels.append(row['chromosome'])
        
        chr1_label_data.append({
            'chromosome': row['chromosome'],
            'chromosome_size': chr_size,
            'cumulative_start': cumulative_offset,
            'cumulative_end': cumulative_offset + chr_size,
            'label_position': midpoint
        })
        
        cumulative_offset += chr_size
    
    species2_data = species2_data.sort_values('chromosome')
    chr2_positions = []
    chr2_labels = []
    cumulative_offset = 0
    
    chr2_label_data = []
    for _, row in species2_data.iterrows():
        chr_size = row['chromsome_size_b']
        midpoint = cumulative_offset + (chr_size / 2)
        chr2_positions.append(midpoint)
        chr2_labels.append(row['chromosome'])
        
        chr2_label_data.append({
            'chromosome': row['chromosome'],
            'chromosome_size': chr_size,
            'cumulative_start': cumulative_offset,
            'cumulative_end': cumulative_offset + chr_size,
            'label_position': midpoint
        })
        
        cumulative_offset += chr_size
    
    # STAGE 13: Save chromosome label positioning data
    if config:
        chr1_labels_df = pd.DataFrame(chr1_label_data)
        chr1_labels_file = data_dir / 'stage_13_chromosome1_label_positions.tsv'
        chr1_labels_df.to_csv(chr1_labels_file, sep='\t', index=False)
        logger.info(f"Saved chromosome1 label positions: {chr1_labels_file}")
        
        chr2_labels_df = pd.DataFrame(chr2_label_data)
        chr2_labels_file = data_dir / 'stage_13_chromosome2_label_positions.tsv'
        chr2_labels_df.to_csv(chr2_labels_file, sep='\t', index=False)
        logger.info(f"Saved chromosome2 label positions: {chr2_labels_file}")
    
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


def linearise_coordinates(joined_df, genome_data, config=None):
    """Convert per-chromosome coordinates to linearised genome coordinates"""
    
    # Create data directory for TSV outputs
    if config:
        output_dir = config.get('output_dir', '.')
        data_dir = Path(output_dir) / 'data'
        data_dir.mkdir(parents=True, exist_ok=True)
    
    # Create chromosome offset mapping for each species
    species1_offsets = {}
    species2_offsets = {}
    
    # Calculate cumulative offsets for species 1
    species1_data = genome_data[genome_data['species'] == CONFIG['first_species_name']]
    species1_data = species1_data.sort_values('chromosome')
    cumulative_offset = 0
    
    species1_offset_data = []
    for _, row in species1_data.iterrows():
        species1_offsets[row['chromosome']] = cumulative_offset
        species1_offset_data.append({
            'chromosome': row['chromosome'],
            'chromosome_size': row['chromsome_size_b'],
            'cumulative_offset': cumulative_offset,
            'cumulative_end': cumulative_offset + row['chromsome_size_b']
        })
        cumulative_offset += row['chromsome_size_b']
    
    # Calculate cumulative offsets for species 2  
    species2_data = genome_data[genome_data['species'] == CONFIG['second_species_name']]
    species2_data = species2_data.sort_values('chromosome')
    cumulative_offset = 0
    
    species2_offset_data = []
    for _, row in species2_data.iterrows():
        species2_offsets[row['chromosome']] = cumulative_offset
        species2_offset_data.append({
            'chromosome': row['chromosome'],
            'chromosome_size': row['chromsome_size_b'],
            'cumulative_offset': cumulative_offset,
            'cumulative_end': cumulative_offset + row['chromsome_size_b']
        })
        cumulative_offset += row['chromsome_size_b']
    
    # STAGE 14: Save chromosome offset data
    if config:
        species1_offsets_df = pd.DataFrame(species1_offset_data)
        species1_offsets_file = data_dir / 'stage_14_species1_chromosome_offsets.tsv'
        species1_offsets_df.to_csv(species1_offsets_file, sep='\t', index=False)
        logger.info(f"Saved species1 chromosome offsets: {species1_offsets_file}")
        
        species2_offsets_df = pd.DataFrame(species2_offset_data)
        species2_offsets_file = data_dir / 'stage_14_species2_chromosome_offsets.tsv'
        species2_offsets_df.to_csv(species2_offsets_file, sep='\t', index=False)
        logger.info(f"Saved species2 chromosome offsets: {species2_offsets_file}")
    
    # Add linearised coordinates to joined_df
    joined_df_linear = joined_df.copy()
    joined_df_linear['linear_start1'] = joined_df_linear.apply(
        lambda row: species1_offsets.get(row['chr1'], 0) + row['start1'], axis=1
    )
    joined_df_linear['linear_start2'] = joined_df_linear.apply(
        lambda row: species2_offsets.get(row['chr2'], 0) + row['start2'], axis=1
    )
    
    # STAGE 15: Save complete linearized coordinates
    if config:
        linearized_file = data_dir / 'stage_15_complete_linearized_coordinates.tsv'
        joined_df_linear.to_csv(linearized_file, sep='\t', index=False)
        logger.info(f"Saved complete linearized coordinates: {linearized_file}")
    
    return joined_df_linear, species1_data, species2_data



def get_chromosome_sizes_from_sheet(CONFIG):
    """Get chromosome sizes from the Google Sheet"""
    sheet_url = "https://docs.google.com/spreadsheets/d/1K01wVWkMW-m6yT9zDX8gDekp-OECubE-9HcmD8RnmkM/edit?gid=1940964825#gid=1940964825"
    # Convert to CSV export URL
    csv_url = sheet_url.replace('/edit?gid=', '/export?format=csv&gid=')
    
    genome_data = pd.read_csv(csv_url)
    
    # Filter for your species
    target_species = {CONFIG["first_species_name"], CONFIG["second_species_name"]}
    species_data = genome_data[genome_data['species'].isin(target_species)]
    
    return species_data[['species', 'chromosome', 'chromsome_size_b']]





def add_chromosome_boundaries(ax, species1_data, species2_data, config=None):
    """Add vertical and horizontal lines at chromosome boundaries"""
    
    # Create data directory for TSV outputs
    if config:
        output_dir = config.get('output_dir', '.')
        data_dir = Path(output_dir) / 'data'
        data_dir.mkdir(parents=True, exist_ok=True)
    
    # Calculate cumulative ends for vertical lines (species 1)
    species1_data = species1_data.sort_values('chromosome')
    cumulative_end = 0
    
    species1_boundaries = []
    for i, (_, row) in enumerate(species1_data.iterrows()):
        cumulative_end += row['chromsome_size_b']
        species1_boundaries.append({
            'chromosome': row['chromosome'],
            'boundary_position': cumulative_end,
            'is_last': i == len(species1_data) - 1
        })
        # Add vertical line at end of each chromosome (except the last one)
        if i < len(species1_data) - 1:
            ax.axvline(x=cumulative_end, color='grey', linewidth=0.8, alpha=0.7)
    
    # Calculate cumulative ends for horizontal lines (species 2)  
    species2_data = species2_data.sort_values('chromosome')
    cumulative_end = 0
    
    species2_boundaries = []
    for i, (_, row) in enumerate(species2_data.iterrows()):
        cumulative_end += row['chromsome_size_b']
        species2_boundaries.append({
            'chromosome': row['chromosome'],
            'boundary_position': cumulative_end,
            'is_last': i == len(species2_data) - 1
        })
        # Add horizontal line at end of each chromosome (except the last one)
        if i < len(species2_data) - 1:
            ax.axhline(y=cumulative_end, color='grey', linewidth=0.8, alpha=0.7)
    
    # STAGE 16: Save chromosome boundary data
    if config:
        species1_boundaries_df = pd.DataFrame(species1_boundaries)
        species1_boundaries_file = data_dir / 'stage_16_species1_chromosome_boundaries.tsv'
        species1_boundaries_df.to_csv(species1_boundaries_file, sep='\t', index=False)
        logger.info(f"Saved species1 chromosome boundaries: {species1_boundaries_file}")
        
        species2_boundaries_df = pd.DataFrame(species2_boundaries)
        species2_boundaries_file = data_dir / 'stage_16_species2_chromosome_boundaries.tsv'
        species2_boundaries_df.to_csv(species2_boundaries_file, sep='\t', index=False)
        logger.info(f"Saved species2 chromosome boundaries: {species2_boundaries_file}")
            
            
            
            
            
            
            


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
    
    # Create data directory for TSV outputs
    output_dir = config.get('output_dir', '.')
    data_dir = Path(output_dir) / 'data'
    data_dir.mkdir(parents=True, exist_ok=True)
    
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
    
    # STAGE 17: Save gene order data for dotplot
    dotplot_data_file = data_dir / 'stage_17_dotplot_gene_order_data.tsv'
    df.to_csv(dotplot_data_file, sep='\t', index=False)
    logger.info(f"Saved dotplot gene order data: {dotplot_data_file}")
    
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
    
    # STAGE 18: Save syntenic vs inverted data for dotplot
    dotplot_syntenic_file = data_dir / 'stage_18_dotplot_syntenic_data.tsv'
    syntenic_data.to_csv(dotplot_syntenic_file, sep='\t', index=False)
    logger.info(f"Saved dotplot syntenic data: {dotplot_syntenic_file}")
    
    dotplot_inverted_file = data_dir / 'stage_18_dotplot_inverted_data.tsv'
    inverted_data.to_csv(dotplot_inverted_file, sep='\t', index=False)
    logger.info(f"Saved dotplot inverted data: {dotplot_inverted_file}")
    
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