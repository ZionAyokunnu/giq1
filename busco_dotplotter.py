#!/usr/bin/env python3
"""
Independent linear BUSCO dotplotter for comparing gene orders between two BUSCO TSV files.
"""


"""
script:

python3 busco_dotplotter.py \
  agora_ancestral_genome.tsv \
  giq_ancestral_genome.tsv \
  agora_vs_giq_comparison.png \
  --name1 "AGORA Ancestral" \
  --name2 "GIQ Ancestral"


"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse
from pathlib import Path
import seaborn as sns

def load_busco_tsv(file_path: str, name: str):
    """Load and process BUSCO TSV file"""
    print(f"Loading {name} from: {file_path}")
    
    try:
        df = pd.read_csv(file_path, sep='\t')
        print(f"  Loaded {len(df)} genes")
        
        # Ensure required columns exist
        required_cols = ['busco_id', 'sequence', 'gene_start', 'gene_end']
        missing_cols = [col for col in required_cols if col not in df.columns]
        
        if missing_cols:
            raise ValueError(f"Missing required columns: {missing_cols}")
        
        # Add 'start' column for compatibility with linearised functions
        df['start'] = df['gene_start']
        
        # Sort by position within each chromosome
        df = df.sort_values(['sequence', 'gene_start'])
        
        # Add cumulative position for linear plotting
        df['linear_position'] = 0
        cumulative_pos = 0
        
        for chrom in sorted(df['sequence'].unique()):
            chrom_mask = df['sequence'] == chrom
            chrom_genes = df[chrom_mask]
            
            # Assign linear positions
            positions = np.arange(cumulative_pos, cumulative_pos + len(chrom_genes))
            df.loc[chrom_mask, 'linear_position'] = positions
            cumulative_pos += len(chrom_genes)
            
            print(f"  {chrom}: {len(chrom_genes)} genes")
        
        return df
        
    except Exception as e:
        print(f"Error loading {name}: {e}")
        raise


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


def linearise_genome_coordinates(genome_df, chromosome_sizes_df=None):
    """Convert per-chromosome coordinates to linearised genome coordinates"""
    
    genome_linear = genome_df.copy()
    
    if chromosome_sizes_df is not None:
        # Use provided chromosome sizes
        chr_offsets = {}
        cumulative_offset = 0
        
        # Sort chromosomes for consistent ordering
        sorted_chrs = sorted(chromosome_sizes_df['chromosome'].unique())
        
        for chr_name in sorted_chrs:
            chr_offsets[chr_name] = cumulative_offset
            chr_size = chromosome_sizes_df[chromosome_sizes_df['chromosome'] == chr_name]['chromosome_size_bp'].iloc[0]
            cumulative_offset += chr_size
    else:
        # Estimate from gene positions
        chr_offsets = {}
        cumulative_offset = 0
        
        # Group by chromosome and find max position
        chr_groups = genome_df.groupby('sequence')['start'].agg(['min', 'max'])
        sorted_chrs = sorted(chr_groups.index)
        
        for chr_name in sorted_chrs:
            chr_offsets[chr_name] = cumulative_offset
            chr_size = chr_groups.loc[chr_name, 'max'] - chr_groups.loc[chr_name, 'min'] + 10000  # Add padding
            cumulative_offset += chr_size
    
    # Add linearised coordinates
    genome_linear['linear_position'] = genome_linear.apply(
        lambda row: chr_offsets.get(row['sequence'], 0) + row['start'], axis=1
    )
    
    return genome_linear, chr_offsets


def create_standard_dotplot(g1_common, g2_common, genome1_name, genome2_name, ax):
    """Create the standard dotplot"""
    
    # Create mapping between genomes
    g1_positions = dict(zip(g1_common['busco_id'], g1_common['start']))
    g2_positions = dict(zip(g2_common['busco_id'], g2_common['start']))
    
    # Prepare data for plotting
    plot_data = []
    for gene in set(g1_common['busco_id']) & set(g2_common['busco_id']):
        if gene in g1_positions and gene in g2_positions:
            plot_data.append({
                'genome1_pos': g1_positions[gene],
                'genome2_pos': g2_positions[gene]
            })
    
    plot_df = pd.DataFrame(plot_data)
    
    # Create scatter plot
    ax.scatter(plot_df['genome1_pos'], plot_df['genome2_pos'], 
               alpha=0.7, s=40, c='blue', edgecolors='darkblue', linewidth=0.5)
    
    ax.set_xlabel(f'{genome1_name} Position')
    ax.set_ylabel(f'{genome2_name} Position')
    ax.set_title('Standard Gene Position Comparison')
    ax.grid(True, alpha=0.3)
    
    # Add correlation
    if len(plot_df) > 1:
        correlation = plot_df['genome1_pos'].corr(plot_df['genome2_pos'])
        ax.text(0.05, 0.95, f'R = {correlation:.3f}', 
                transform=ax.transAxes, bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))


def add_chromosome_boundaries(ax, g1_offsets, g2_offsets):
    """Add lines at chromosome boundaries"""
    
    # Add vertical lines for genome1 chromosome boundaries
    for chr_name, offset in g1_offsets.items():
        if offset > 0:  # Don't add line at position 0
            ax.axvline(x=offset, color='grey', linewidth=0.8, alpha=0.5, linestyle='--')
    
    # Add horizontal lines for genome2 chromosome boundaries  
    for chr_name, offset in g2_offsets.items():
        if offset > 0:  # Don't add line at position 0
            ax.axhline(y=offset, color='grey', linewidth=0.8, alpha=0.5, linestyle='--')


def add_chromosome_labels_simple(ax, g1_offsets, g2_offsets):
    """Add simple chromosome labels"""
    
    # Get current axis limits
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    
    # Add chromosome names as text annotations
    sorted_g1_chrs = sorted(g1_offsets.items(), key=lambda x: x[1])
    for i, (chr_name, offset) in enumerate(sorted_g1_chrs):
        if offset < xlim[1]:
            ax.text(offset + (xlim[1] * 0.02), ylim[1] * 0.98, chr_name, 
                   rotation=90, ha='left', va='top', fontsize=8, alpha=0.7)
    
    sorted_g2_chrs = sorted(g2_offsets.items(), key=lambda x: x[1])
    for i, (chr_name, offset) in enumerate(sorted_g2_chrs):
        if offset < ylim[1]:
            ax.text(xlim[1] * 0.98, offset + (ylim[1] * 0.02), chr_name, 
                   rotation=0, ha='right', va='bottom', fontsize=8, alpha=0.7)


def create_linearised_dotplot(g1_common, g2_common, genome1_name, genome2_name, ax, chromosome_sizes_df=None):
    """Create the linearised genome dotplot"""
    
    # Linearise coordinates for both genomes
    g1_linear, g1_offsets = linearise_genome_coordinates(g1_common, chromosome_sizes_df)
    g2_linear, g2_offsets = linearise_genome_coordinates(g2_common, chromosome_sizes_df)
    
    # Create mapping between genomes using linear coordinates
    g1_positions = dict(zip(g1_linear['busco_id'], g1_linear['linear_position']))
    g2_positions = dict(zip(g2_linear['busco_id'], g2_linear['linear_position']))
    
    # Prepare data for plotting
    plot_data = []
    for gene in set(g1_linear['busco_id']) & set(g2_linear['busco_id']):
        if gene in g1_positions and gene in g2_positions:
            plot_data.append({
                'genome1_linear': g1_positions[gene],
                'genome2_linear': g2_positions[gene]
            })
    
    plot_df = pd.DataFrame(plot_data)
    
    # Create scatter plot
    ax.scatter(plot_df['genome1_linear'], plot_df['genome2_linear'], 
               alpha=0.6, s=30, c='red', edgecolors='darkred', linewidth=0.3)
    
    # Add chromosome boundaries
    add_chromosome_boundaries(ax, g1_offsets, g2_offsets)
    
    # Format axes to show positions in Mb
    ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f'{x/1e6:.1f}'))
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda y, p: f'{y/1e6:.1f}'))
    
    ax.set_xlabel(f'{genome1_name} Linearised Position (Mb)')
    ax.set_ylabel(f'{genome2_name} Linearised Position (Mb)')
    ax.set_title('Linearised Genome Comparison\n(Chromosomes end-to-end)')
    ax.grid(True, alpha=0.2)
    
    # Add chromosome labels
    add_chromosome_labels_simple(ax, g1_offsets, g2_offsets)
    
    # Calculate correlation
    correlation = 0.0
    if len(plot_df) > 1:
        correlation = plot_df['genome1_linear'].corr(plot_df['genome2_linear'])
        ax.text(0.05, 0.95, f'Linear R = {correlation:.3f}', 
                transform=ax.transAxes, bbox=dict(boxstyle="round,pad=0.3", facecolor="yellow", alpha=0.8))
    
    return {
        'linear_correlation': correlation,
        'total_linearised_genes': len(plot_df),
        'g1_offsets': g1_offsets,
        'g2_offsets': g2_offsets
    }


def create_stats_panel(genome1_df, genome2_df, g1_common, g2_common, genome1_name, genome2_name, ax, linear_stats):
    """Create statistics panel"""
    
    ax.axis('off')
    
    stats_text = f"""
GENOME COMPARISON SUMMARY

{genome1_name}: {len(genome1_df):,} total genes
{genome2_name}: {len(genome2_df):,} total genes

Common BUSCO genes: {len(g1_common):,}

CHROMOSOMES:
{genome1_name}: {len(g1_common['sequence'].unique())} chromosomes
{genome2_name}: {len(g2_common['sequence'].unique())} chromosomes

LINEARISED ANALYSIS:
Linear correlation: {linear_stats['linear_correlation']:.3f}

INTERPRETATION:
• R > 0.8: High synteny conservation
• R 0.5-0.8: Moderate rearrangement
• R < 0.5: Extensive rearrangement

The linearised plot shows gene order
across entire genomes arranged 
end-to-end by chromosome.
"""
    
    ax.text(0.05, 0.95, stats_text, transform=ax.transAxes, 
           fontsize=10, verticalalignment='top', fontfamily='monospace')


def create_chromosome_distributions(g1_common, g2_common, genome1_name, genome2_name, ax):
    """Create chromosome distribution comparison"""
    
    chr_counts1 = g1_common['sequence'].value_counts()
    chr_counts2 = g2_common['sequence'].value_counts()
    
    # Create side-by-side bar plot
    x = np.arange(len(chr_counts1))
    width = 0.35
    
    ax.bar(x - width/2, chr_counts1.values, width, label=genome1_name, alpha=0.7, color='lightblue')
    
    # Align genome2 chromosomes if they match
    chr_counts2_aligned = []
    for chr_name in chr_counts1.index:
        if chr_name in chr_counts2.index:
            chr_counts2_aligned.append(chr_counts2[chr_name])
        else:
            chr_counts2_aligned.append(0)
    
    ax.bar(x + width/2, chr_counts2_aligned, width, label=genome2_name, alpha=0.7, color='lightgreen')
    
    ax.set_xlabel('Chromosome')
    ax.set_ylabel('Common Gene Count')
    ax.set_title('Chromosome-wise Gene Distribution')
    ax.set_xticks(x)
    ax.set_xticklabels(chr_counts1.index, rotation=45, ha='right')
    ax.legend()


def create_chromosome_info_panel(g1_common, g2_common, genome1_name, genome2_name, ax):
    """Create detailed chromosome information panel"""
    
    ax.axis('off')
    
    info_text = "CHROMOSOME OVERLAP ANALYSIS:\n\n"
    
    # Calculate chromosome-to-chromosome gene overlaps
    for chr1 in sorted(g1_common['sequence'].unique()):
        chr1_genes = set(g1_common[g1_common['sequence'] == chr1]['busco_id'])
        
        best_match = ""
        best_overlap = 0
        
        for chr2 in sorted(g2_common['sequence'].unique()):
            chr2_genes = set(g2_common[g2_common['sequence'] == chr2]['busco_id'])
            overlap = len(chr1_genes & chr2_genes)
            
            if overlap > best_overlap:
                best_overlap = overlap
                best_match = chr2
        
        if best_overlap > 0:
            info_text += f"{chr1} → {best_match}: {best_overlap} genes\n"
    
    ax.text(0.05, 0.95, info_text, transform=ax.transAxes, 
           fontsize=9, verticalalignment='top', fontfamily='monospace')


def create_error_plot(genome1_df, genome2_df, genome1_name, genome2_name, output_path):
    """Create error plot when no common genes found"""
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    
    # Empty main plot
    ax_main = axes[0, 0]
    ax_main.text(0.5, 0.5, 'NO COMMON GENES FOUND\nCheck BUSCO ID formats', 
                ha='center', va='center', transform=ax_main.transAxes, 
                fontsize=16, color='red', weight='bold')
    ax_main.set_xlabel(f'{genome1_name} Linear Position')
    ax_main.set_ylabel(f'{genome2_name} Linear Position')
    ax_main.set_title(f'Linear Gene Order Comparison\n{genome1_name} vs {genome2_name}')
    
    # Show some sample BUSCO IDs from each genome
    ax_chr1 = axes[0, 1]
    ax_chr1.axis('off')
    sample1 = list(genome1_df['busco_id'].head(10))
    ax_chr1.text(0.05, 0.95, f'{genome1_name} Sample IDs:\n' + '\n'.join(sample1), 
                transform=ax_chr1.transAxes, fontsize=10, verticalalignment='top')
    
    ax_chr2 = axes[1, 0]
    ax_chr2.axis('off')
    sample2 = list(genome2_df['busco_id'].head(10))
    ax_chr2.text(0.05, 0.95, f'{genome2_name} Sample IDs:\n' + '\n'.join(sample2), 
                transform=ax_chr2.transAxes, fontsize=10, verticalalignment='top')
    
    # Error summary
    ax_stats = axes[1, 1]
    ax_stats.axis('off')
    stats_text = f"""
ERROR: NO COMMON GENES FOUND

{genome1_name}: {len(genome1_df)} genes
{genome2_name}: {len(genome2_df)} genes

Possible causes:
• Different BUSCO ID formats
• Different gene sets
• Parsing errors

Check that both files contain 
valid BUSCO IDs like:
100019at7147, 77167at7147, etc.
"""
    ax_stats.text(0.05, 0.95, stats_text, transform=ax_stats.transAxes, fontsize=9, verticalalignment='top', fontfamily='monospace')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Dotplot saved to: {output_path}")
    
    # Return summary statistics
    return {
        'total_genes_genome1': len(genome1_df),
        'total_genes_genome2': len(genome2_df),
        'common_genes': 0,
        'correlation': 0.0,
        'chromosomes_genome1': len(genome1_df['sequence'].unique()),
        'chromosomes_genome2': len(genome2_df['sequence'].unique())
    }


def create_linear_dotplot_with_linearised(genome1_df, genome2_df, genome1_name, genome2_name, output_path, chromosome_sizes_df=None):
    """Create both standard and linearised dotplots comparing two genomes"""
    
    # Find common genes
    common_genes = set(genome1_df['busco_id']) & set(genome2_df['busco_id'])
    print(f"Common genes found: {len(common_genes)}")
    
    if len(common_genes) == 0:
        print("Warning: No common genes found between the two genomes!")
        # Keep the original error handling code here
        return create_error_plot(genome1_df, genome2_df, genome1_name, genome2_name, output_path)
    
    # Filter to common genes only
    g1_common = genome1_df[genome1_df['busco_id'].isin(common_genes)].copy()
    g2_common = genome2_df[genome2_df['busco_id'].isin(common_genes)].copy()
    
    # Create figure with 3 subplots: standard dotplot, linearised dotplot, and stats
    fig = plt.figure(figsize=(20, 12))
    
    # Define grid layout
    gs = fig.add_gridspec(2, 3, width_ratios=[1, 1.2, 0.8], height_ratios=[1, 1])
    
    # Standard dotplot (top left)
    ax_standard = fig.add_subplot(gs[0, 0])
    create_standard_dotplot(g1_common, g2_common, genome1_name, genome2_name, ax_standard)
    
    # Linearised dotplot (top right and bottom right - spans both)
    ax_linear = fig.add_subplot(gs[:, 1])
    linear_stats = create_linearised_dotplot(g1_common, g2_common, genome1_name, genome2_name, ax_linear, chromosome_sizes_df)
    
    # Statistics and chromosome info (top right)
    ax_stats = fig.add_subplot(gs[0, 2])
    create_stats_panel(genome1_df, genome2_df, g1_common, g2_common, genome1_name, genome2_name, ax_stats, linear_stats)
    
    # Chromosome distributions (bottom left)
    ax_chr = fig.add_subplot(gs[1, 0])
    create_chromosome_distributions(g1_common, g2_common, genome1_name, genome2_name, ax_chr)
    
    # Bottom right space for additional chromosome info
    ax_chr_info = fig.add_subplot(gs[1, 2])
    create_chromosome_info_panel(g1_common, g2_common, genome1_name, genome2_name, ax_chr_info)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"Enhanced dotplot saved to: {output_path}")
    
    return linear_stats


def create_linear_dotplot(genome1_df, genome2_df, genome1_name, genome2_name, output_path):
    """Create linear dotplot comparing two genomes - Fixed to handle empty common genes"""
    
    # Find common genes
    common_genes = set(genome1_df['busco_id']) & set(genome2_df['busco_id'])
    print(f"Common genes found: {len(common_genes)}")
    
    if len(common_genes) == 0:
        print("Warning: No common genes found between the two genomes!")
        
        # Still create a basic comparison plot
        fig, axes = plt.subplots(2, 2, figsize=(16, 12))
        
        # Empty main plot
        ax_main = axes[0, 0]
        ax_main.text(0.5, 0.5, 'NO COMMON GENES FOUND\nCheck BUSCO ID formats', 
                    ha='center', va='center', transform=ax_main.transAxes, 
                    fontsize=16, color='red', weight='bold')
        ax_main.set_xlabel(f'{genome1_name} Linear Position')
        ax_main.set_ylabel(f'{genome2_name} Linear Position')
        ax_main.set_title(f'Linear Gene Order Comparison\n{genome1_name} vs {genome2_name}')
        
        # Show some sample BUSCO IDs from each genome
        ax_chr1 = axes[0, 1]
        ax_chr1.axis('off')
        sample1 = list(genome1_df['busco_id'].head(10))
        ax_chr1.text(0.05, 0.95, f'{genome1_name} Sample IDs:\n' + '\n'.join(sample1), 
                    transform=ax_chr1.transAxes, fontsize=10, verticalalignment='top')
        
        ax_chr2 = axes[1, 0]
        ax_chr2.axis('off')
        sample2 = list(genome2_df['busco_id'].head(10))
        ax_chr2.text(0.05, 0.95, f'{genome2_name} Sample IDs:\n' + '\n'.join(sample2), 
                    transform=ax_chr2.transAxes, fontsize=10, verticalalignment='top')
        
        # Error summary
        ax_stats = axes[1, 1]
        ax_stats.axis('off')
        stats_text = f"""
ERROR: NO COMMON GENES FOUND

{genome1_name}: {len(genome1_df)} genes
{genome2_name}: {len(genome2_df)} genes

Possible causes:
• Different BUSCO ID formats
• Different gene sets
• Parsing errors

Check that both files contain 
valid BUSCO IDs like:
100019at7147, 77167at7147, etc.
"""
        ax_stats.text(0.05, 0.95, stats_text, transform=ax_stats.transAxes, 
                     fontsize=12, verticalalignment='top', color='red')
        
        plt.tight_layout()
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"Error plot saved to: {output_path}")
        
        # Return error stats
        return {
            'total_genes_genome1': len(genome1_df),
            'total_genes_genome2': len(genome2_df),
            'common_genes': 0,
            'correlation': 0.0,
            'chromosomes_genome1': len(genome1_df['sequence'].unique()),
            'chromosomes_genome2': len(genome2_df['sequence'].unique())
        }
    
    # Rest of the original function for when common genes exist...
    # Filter to common genes only
    g1_common = genome1_df[genome1_df['busco_id'].isin(common_genes)].copy()
    g2_common = genome2_df[genome2_df['busco_id'].isin(common_genes)].copy()
    
    # Create mapping between genomes
    g1_positions = dict(zip(g1_common['busco_id'], g1_common['linear_position']))
    g2_positions = dict(zip(g2_common['busco_id'], g2_common['linear_position']))
    
    # Prepare data for plotting
    plot_data = []
    for gene in common_genes:
        if gene in g1_positions and gene in g2_positions:
            plot_data.append({
                'gene': gene,
                'genome1_pos': g1_positions[gene],
                'genome2_pos': g2_positions[gene]
            })
    
    plot_df = pd.DataFrame(plot_data)
    
    # Create the plot
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    
    # Main dotplot
    ax_main = axes[0, 0]
    scatter = ax_main.scatter(plot_df['genome1_pos'], plot_df['genome2_pos'], 
                             alpha=0.7, s=60, c='blue', edgecolors='darkblue', linewidth=0.5)
    
    ax_main.set_xlabel(f'{genome1_name} Linear Position')
    ax_main.set_ylabel(f'{genome2_name} Linear Position')
    ax_main.set_title(f'Linear Gene Order Comparison\n{genome1_name} vs {genome2_name}')
    ax_main.grid(True, alpha=0.3)
    
    # Add diagonal for perfect correlation
    max_pos = max(plot_df['genome1_pos'].max(), plot_df['genome2_pos'].max())
    ax_main.plot([0, max_pos], [0, max_pos], 'r--', alpha=0.5, linewidth=2, label='Perfect correlation')
    ax_main.legend()
    
    # Correlation statistics
    correlation = plot_df['genome1_pos'].corr(plot_df['genome2_pos'])
    ax_main.text(0.05, 0.95, f'Correlation: {correlation:.3f}', 
                transform=ax_main.transAxes, bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))
    
    # Chromosome distribution for genome1
    ax_chr1 = axes[0, 1]
    chr_counts1 = g1_common['sequence'].value_counts()
    ax_chr1.bar(range(len(chr_counts1)), chr_counts1.values, alpha=0.7, color='lightblue')
    ax_chr1.set_xticks(range(len(chr_counts1)))
    ax_chr1.set_xticklabels(chr_counts1.index, rotation=45, ha='right')
    ax_chr1.set_ylabel('Gene Count')
    ax_chr1.set_title(f'{genome1_name} Chromosome Distribution')
    
    # Chromosome distribution for genome2
    ax_chr2 = axes[1, 0]
    chr_counts2 = g2_common['sequence'].value_counts()
    ax_chr2.bar(range(len(chr_counts2)), chr_counts2.values, alpha=0.7, color='lightgreen')
    ax_chr2.set_xticks(range(len(chr_counts2)))
    ax_chr2.set_xticklabels(chr_counts2.index, rotation=45, ha='right')
    ax_chr2.set_ylabel('Gene Count')
    ax_chr2.set_title(f'{genome2_name} Chromosome Distribution')
    
    # Summary statistics
    ax_stats = axes[1, 1]
    ax_stats.axis('off')
    
    stats_text = f"""
COMPARISON SUMMARY

Total genes compared: {len(plot_data)}
{genome1_name}: {len(genome1_df)} genes, {len(genome1_df['sequence'].unique())} chromosomes
{genome2_name}: {len(genome2_df)} genes, {len(genome2_df['sequence'].unique())} chromosomes

Linear correlation: {correlation:.3f}

INTERPRETATION:
• Correlation > 0.7: High synteny conservation
• Correlation 0.3-0.7: Moderate synteny
• Correlation < 0.3: Low synteny / high rearrangement

COMMON GENES BY CHROMOSOME:
"""
    
    # Add chromosome overlap info
    for chr1 in sorted(g1_common['sequence'].unique()):
        chr1_genes = set(g1_common[g1_common['sequence'] == chr1]['busco_id'])
        
        for chr2 in sorted(g2_common['sequence'].unique()):
            chr2_genes = set(g2_common[g2_common['sequence'] == chr2]['busco_id'])
            overlap = len(chr1_genes & chr2_genes)
            
            if overlap > 0:
                stats_text += f"\n{chr1} ↔ {chr2}: {overlap} genes"
    
    ax_stats.text(0.05, 0.95, stats_text, transform=ax_stats.transAxes, fontsize=9, verticalalignment='top', fontfamily='monospace')

def main():
    parser = argparse.ArgumentParser(description='Create linear BUSCO dotplot comparison')
    parser.add_argument('genome1_tsv', help='First genome BUSCO TSV file')
    parser.add_argument('genome2_tsv', help='Second genome BUSCO TSV file')
    parser.add_argument('output_plot', help='Output plot file (PNG/PDF)')
    parser.add_argument('--name1', default='Genome1', help='Name for first genome')
    parser.add_argument('--name2', default='Genome2', help='Name for second genome')
    parser.add_argument('--use-linearised', action='store_true', help='Use enhanced linearised dotplot')
    parser.add_argument('--use-chromosome-sizes', action='store_true', help='Use chromosome sizes from Google Sheet')
    
    args = parser.parse_args()
    
    try:
        # Load genomes
        genome1_df = load_busco_tsv(args.genome1_tsv, args.name1)
        genome2_df = load_busco_tsv(args.genome2_tsv, args.name2)
        
        # Load chromosome sizes if requested
        chromosome_sizes_df = None
        if args.use_chromosome_sizes:
            try:
                chromosome_sizes_df = get_chromosome_sizes_from_sheet()
                print(f"Loaded chromosome sizes for {len(chromosome_sizes_df)} chromosomes")
            except Exception as e:
                print(f"Warning: Could not load chromosome sizes from Google Sheet: {e}")
                print("Continuing with estimated sizes from gene positions...")
        
        # Create dotplot
        if args.use_linearised:
            stats = create_linear_dotplot_with_linearised(genome1_df, genome2_df, 
                                        args.name1, args.name2, 
                                        args.output_plot, chromosome_sizes_df)
        else:
            stats = create_linear_dotplot(genome1_df, genome2_df, 
                                        args.name1, args.name2, 
                                        args.output_plot)
        
        print("\n" + "="*50)
        print("DOTPLOT COMPARISON COMPLETE")
        print("="*50)
        if args.use_linearised:
            print(f"Common genes: {stats.get('total_linearised_genes', 0)}")
            print(f"Linear correlation: {stats.get('linear_correlation', 0.0):.3f}")
        else:
            print(f"Common genes: {stats['common_genes']}")
            print(f"Linear correlation: {stats['correlation']:.3f}")
        print(f"Plot saved to: {args.output_plot}")
        
    except Exception as e:
        print(f"❌ Error: {e}")
        return 1
    
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main()), 
                 fontsize=12, verticalalignment='top', color='red')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Error plot saved to: {output_path}")
    
