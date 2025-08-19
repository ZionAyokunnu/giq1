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
    

    plots_dir = Path(plots_dir)
    plots_dir.mkdir(parents=True, exist_ok=True)
    

    data_dir = plots_dir.parent / 'data'
    data_dir.mkdir(parents=True, exist_ok=True)
    

    genome_data = get_chromosome_sizes_from_sheet(CONFIG)
    

    
    joined_df_linear, species1_data, species2_data = linearise_coordinates(joined_df, genome_data, plot_config)
    

    fig, ax = plt.subplots(figsize=plot_config['figure_size'])
    

    syntenic_data = joined_df_linear[joined_df_linear['is_flipped'] == False]
    inverted_data = joined_df_linear[joined_df_linear['is_flipped'] == True]
    
    
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
    

    add_chromosome_boundaries(ax, species1_data, species2_data, plot_config)
    
    add_chromosome_labels(ax, species1_data, species2_data, joined_df_linear, plot_config)

    ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f'{x/1e6:.0f}'))
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda y, p: f'{y/1e6:.0f}'))
    

    ax.set_xlabel(f'{CONFIG["first_species_name"]} (Mb)', 
                 fontsize=plot_config['font_size_title'])
    ax.set_ylabel(f'{CONFIG["second_species_name"]} (Mb)', 
                 fontsize=plot_config['font_size_title'])
    

    ax.tick_params(axis='both', labelsize=plot_config['font_size_labels'])
    ax.legend()


    ax.set_facecolor('white')
    ax.grid(False)
    for spine in ax.spines.values():
        spine.set_linewidth(0.8)
        spine.set_color('black')
    

    png_path = plots_dir / "busco_dotplot.png"
    fig.savefig(png_path, dpi=plot_config['dpi'], bbox_inches='tight', facecolor='white')
    
    plt.tight_layout()
    return fig, ax



def add_chromosome_labels(ax, species1_data, species2_data, joined_df_linear, config=None):

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
    

    ax2 = ax.secondary_xaxis('top')
    ax3 = ax.secondary_yaxis('right')
    

    ax2.set_xticks(chr1_positions)
    ax2.set_xticklabels(chr1_labels, 
                       rotation=45, 
                       ha='left',
                       fontsize=12)
    

    ax3.set_yticks(chr2_positions)
    ax3.set_yticklabels(chr2_labels, 
                       fontsize=12)
    

    ax2.tick_params(length=0)
    ax3.tick_params(length=0)


def linearise_coordinates(joined_df, genome_data, config=None):

    if config:
        output_dir = config.get('output_dir', '.')
        data_dir = Path(output_dir) / 'data'
        data_dir.mkdir(parents=True, exist_ok=True)
    

    species1_offsets = {}
    species2_offsets = {}
    

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
    


    joined_df_linear = joined_df.copy()
    joined_df_linear['linear_start1'] = joined_df_linear.apply(
        lambda row: species1_offsets.get(row['chr1'], 0) + row['start1'], axis=1
    )
    joined_df_linear['linear_start2'] = joined_df_linear.apply(
        lambda row: species2_offsets.get(row['chr2'], 0) + row['start2'], axis=1
    )
    
    return joined_df_linear, species1_data, species2_data



def get_chromosome_sizes_from_sheet(CONFIG):

    sheet_url = "https://docs.google.com/spreadsheets/d/1K01wVWkMW-m6yT9zDX8gDekp-OECubE-9HcmD8RnmkM/edit?gid=1940964825#gid=1940964825"

    csv_url = sheet_url.replace('/edit?gid=', '/export?format=csv&gid=')
    
    genome_data = pd.read_csv(csv_url)
    

    target_species = {CONFIG["first_species_name"], CONFIG["second_species_name"]}
    species_data = genome_data[genome_data['species'].isin(target_species)]
    
    return species_data[['species', 'chromosome', 'chromsome_size_b']]





def add_chromosome_boundaries(ax, species1_data, species2_data, config=None):


    if config:
        output_dir = config.get('output_dir', '.')
        data_dir = Path(output_dir) / 'data'
        data_dir.mkdir(parents=True, exist_ok=True)
    

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

        if i < len(species1_data) - 1:
            ax.axvline(x=cumulative_end, color='grey', linewidth=0.8, alpha=0.7)
    
 
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
    if len(ortholog_df) == 0:
        return
    
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
    plt.savefig(os.path.join(plots_dir, 'crunched_dotplot.png'), dpi=300, bbox_inches='tight')
    plt.close()