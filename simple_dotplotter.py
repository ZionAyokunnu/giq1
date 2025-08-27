#!/usr/bin/env python3
"""
Independent simple BUSCO dotplot script for testing genome comparisons.
Creates clean diagonal dotplots using actual genomic coordinates.

Usage:
python3 simple_dotplot.py genome1.tsv genome2.tsv output.png [--name1 "Genome 1"] [--name2 "Genome 2"] [--chr-order1 OZ002744.1 chr2 chr1 chr4 chr3 chr5]

python3 simple_dotplotter.py \
    /Users/zionayokunnu/Documents/Giq/compare/root_agora_ancestral_genome.tsv \
    /Users/zionayokunnu/Documents/Bibionidae/busco-tables/Dioctria_rufipes.tsv \
    compare/experiment/simple/agora-chr-aware_vs_Dioctria_rufipes.png \
    --name1 "Agora -chr-aware " \
    --name2 "Dioctria rufipes" \
    --chr-order2 OZ002741.1 OZ002745.1 OZ002740.1 OZ002743.1 OZ002744.1 OZ002743.1 OZ002742.1
    
    
python3 simple_dotplotter.py \
    /Users/zionayokunnu/Documents/Bibionidae/busco-tables/Dioctria_rufipes.tsv \
    compare/Dioctria_linearis_linearized.tsv \
    compare/experiment/simple/Dioctria_rufipes_original_vs_Dioctria_linearis_linearised.png \
    --name1 "Dioctria rufipes_original" \
    --name2 "Dioctria linearis_linearised"
    
"""



import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse
import sys
from pathlib import Path


def load_busco_file(file_path, name):
    """Load BUSCO file handling both formats"""
    print(f"Loading {name}: {file_path}")
    
    try:
        # Check if file has comment headers
        with open(file_path, 'r') as f:
            first_line = f.readline().strip()
        
        if first_line.startswith('#'):
            # BUSCO format with comments
            df = pd.read_csv(file_path, sep='\t', comment='#', header=None, 
                             names=['busco_id', 'status', 'sequence', 'gene_start', 'gene_end', 'strand', 'score', 'length'])
            
            # Convert positions to numeric
            df['gene_start'] = pd.to_numeric(df['gene_start'], errors='coerce')
            df['gene_end'] = pd.to_numeric(df['gene_end'], errors='coerce')
            
            # Filter to complete genes only
            df = df[df['status'] == 'Complete'].copy()
            print(f"  Filtered to {len(df)} complete genes")
        else:
            # Clean TSV format
            df = pd.read_csv(file_path, sep='\t')
            df['gene_start'] = pd.to_numeric(df['gene_start'], errors='coerce')
            df['gene_end'] = pd.to_numeric(df['gene_end'], errors='coerce')
        
        print(f"  Loaded {len(df):,} genes from {len(df['sequence'].unique())} chromosomes")
        
        # Show chromosome distribution
        chr_counts = df['sequence'].value_counts().sort_index()
        for chr_name, count in chr_counts.items():
            print(f"    {chr_name}: {count:,} genes")
        
        return df
        
    except Exception as e:
        print(f"Error loading {name}: {e}")
        raise


def create_simple_busco_dotplot(genome1_df, genome2_df, genome1_name, genome2_name, output_path, chr_order1=None, chr_order2=None):
    """
    Create simple diagonal dotplot using linearized genomic coordinates
    Both genomes get linearized if they have multiple chromosomes
    
    Args:
        genome1_df, genome2_df: BUSCO DataFrames
        genome1_name, genome2_name: Display names
        output_path: Where to save the plot
        
    Returns:
        correlation: Pearson correlation coefficient
    """
    
    print(f"\nCreating simple BUSCO dotplot: {genome1_name} vs {genome2_name}")
    
    # Find common genes
    common_genes = set(genome1_df['busco_id']) & set(genome2_df['busco_id'])
    print(f"  Common genes: {len(common_genes):,}")
    
    if len(common_genes) == 0:
        print("  No common genes found!")
        return None
    
    # Filter to common genes only
    g1_common = genome1_df[genome1_df['busco_id'].isin(common_genes)].copy()
    g2_common = genome2_df[genome2_df['busco_id'].isin(common_genes)].copy()
    
    # Linearize genome1 if it has multiple chromosomes and store boundary info
    chr_boundaries_g1 = None
    if len(g1_common['sequence'].unique()) > 1:
        print(f"  Linearizing {genome1_name} ({len(g1_common['sequence'].unique())} chromosomes)")
        g1_order = get_chromosome_order(g1_common, chr_order1)
        g1_linearized = linearize_real_genome_coordinates(g1_common, g1_order)
        g1_positions = dict(zip(g1_linearized['busco_id'], g1_linearized['linearized_real_position']))
        
        # Calculate chromosome boundaries for axis labels
        chr_boundaries_g1 = {}
        cumulative_offset = 0
        
        g1_order = get_chromosome_order(g1_common, chr_order1)
        for chrom in g1_order:
            chrom_genes = g1_common[g1_common['sequence'] == chrom]
            chr_size = chrom_genes['gene_end'].max()
            chr_boundaries_g1[chrom] = {
                'start': cumulative_offset,
                'end': cumulative_offset + chr_size,
                'midpoint': cumulative_offset + (chr_size / 2)
            }
            cumulative_offset += chr_size
            
    else:
        print(f"  Using direct positions for {genome1_name} (single chromosome)")
        g1_positions = dict(zip(g1_common['busco_id'], g1_common['gene_start']))
    
    # Linearize genome2 if it has multiple chromosomes and store boundary info
    chr_boundaries_g2 = None
    if len(g2_common['sequence'].unique()) > 1:
        print(f"  Linearizing {genome2_name} ({len(g2_common['sequence'].unique())} chromosomes)")
        g2_order = get_chromosome_order(g2_common, chr_order2)
        g2_linearized = linearize_real_genome_coordinates(g2_common, g2_order)
        g2_positions = dict(zip(g2_linearized['busco_id'], g2_linearized['linearized_real_position']))
        
        # Calculate chromosome boundaries for axis labels
        chr_boundaries_g2 = {}
        cumulative_offset = 0
        
        g2_order = get_chromosome_order(g2_common, chr_order2)
        for chrom in g2_order:
            chrom_genes = g2_common[g2_common['sequence'] == chrom]
            chr_size = chrom_genes['gene_end'].max()
            chr_boundaries_g2[chrom] = {
                'start': cumulative_offset,
                'end': cumulative_offset + chr_size,
                'midpoint': cumulative_offset + (chr_size / 2)
            }
            cumulative_offset += chr_size
            
    else:
        print(f"  Using direct positions for {genome2_name} (single chromosome)")
        g2_positions = dict(zip(g2_common['busco_id'], g2_common['gene_start']))
    
    # Also create chromosome mappings for coloring
    g1_chromosomes = dict(zip(g1_common['busco_id'], g1_common['sequence']))
    g2_chromosomes = dict(zip(g2_common['busco_id'], g2_common['sequence']))
    
    # Prepare plotting data
    plot_data = []
    for gene in common_genes:
        if gene in g1_positions and gene in g2_positions:
            plot_data.append({
                'gene_id': gene,
                'x_pos': g1_positions[gene],
                'y_pos': g2_positions[gene], 
                'x_chr': g1_chromosomes[gene],
                'y_chr': g2_chromosomes[gene],
                'same_chr': g1_chromosomes[gene] == g2_chromosomes[gene]
            })
    
    plot_df = pd.DataFrame(plot_data)
    print(f"  Plotting {len(plot_df):,} gene pairs")
    
    # Calculate correlation
    correlation = 0.0
    if len(plot_df) > 1:
        correlation = plot_df['x_pos'].corr(plot_df['y_pos'])
        print(f"  Overall correlation: {correlation:.3f}")
    
    # Create figure
    plt.figure(figsize=(12, 10))
    
    # Separate same vs different chromosome genes
    same_chr = plot_df[plot_df['same_chr'] == True]
    diff_chr = plot_df[plot_df['same_chr'] == False]
    
    # Plot different chromosomes first (background)
    if len(diff_chr) > 0:
        plt.scatter(diff_chr['x_pos'], diff_chr['y_pos'], 
                   c='lightcoral', alpha=0.6, s=20, 
                   label=f'Different chromosomes ({len(diff_chr):,})',
                   edgecolors='none')
    
    # Plot same chromosomes on top (foreground)
    if len(same_chr) > 0:
        plt.scatter(same_chr['x_pos'], same_chr['y_pos'], 
                   c='steelblue', alpha=0.7, s=30,
                   label=f'Same chromosomes ({len(same_chr):,})',
                   edgecolors='none')
    
    # Add perfect correlation line
    min_pos = min(plot_df['x_pos'].min(), plot_df['y_pos'].min())
    max_pos = max(plot_df['x_pos'].max(), plot_df['y_pos'].max())
    plt.plot([min_pos, max_pos], [min_pos, max_pos], 
             'r--', alpha=0.5, linewidth=1.5, label='Perfect correlation')
    
    # Formatting
    plt.xlabel(f'{genome1_name} Linearized Position (bp)', fontsize=12)
    plt.ylabel(f'{genome2_name} Position (bp)', fontsize=12)
    plt.title(f'BUSCO Linearized Position Comparison\n{genome1_name} vs {genome2_name}\nR = {correlation:.3f}', 
              fontsize=14)
    
    # Format axes to show Mb
    plt.gca().xaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f'{x/1e6:.0f}M'))
    plt.gca().yaxis.set_major_formatter(plt.FuncFormatter(lambda y, p: f'{y/1e6:.0f}M'))
    
    # Add chromosome boundaries and labels for both genomes
    if chr_boundaries_g1:
        print(f"  Adding X-axis chromosome boundaries for {len(chr_boundaries_g1)} chromosomes")
        
        # Add vertical lines at chromosome boundaries
        for chrom, info in chr_boundaries_g1.items():
            if info['start'] > 0:  # Don't add line at position 0
                plt.axvline(x=info['start'], color='gray', linewidth=1, alpha=0.6, linestyle='--')
        
        # Add chromosome labels on top x-axis
        ax_top = plt.gca().secondary_xaxis('top')
        
        # Set chromosome labels at midpoints
        chr_positions = [info['midpoint'] for info in chr_boundaries_g1.values()]
        chr_labels = list(chr_boundaries_g1.keys())
        
        ax_top.set_xticks(chr_positions)
        ax_top.set_xticklabels(chr_labels, rotation=45, ha='left', fontsize=10)
        ax_top.tick_params(length=0)  # Remove tick marks
    
    if chr_boundaries_g2:
        print(f"  Adding Y-axis chromosome boundaries for {len(chr_boundaries_g2)} chromosomes")
        
        # Add horizontal lines at chromosome boundaries
        for chrom, info in chr_boundaries_g2.items():
            if info['start'] > 0:  # Don't add line at position 0
                plt.axhline(y=info['start'], color='gray', linewidth=1, alpha=0.6, linestyle='--')
        
        # Add chromosome labels on right y-axis
        ax_right = plt.gca().secondary_yaxis('right')
        
        # Set chromosome labels at midpoints
        chr_positions = [info['midpoint'] for info in chr_boundaries_g2.values()]
        chr_labels = list(chr_boundaries_g2.keys())
        
        ax_right.set_yticks(chr_positions)
        ax_right.set_yticklabels(chr_labels, rotation=0, ha='left', fontsize=10)
        ax_right.tick_params(length=0)  # Remove tick marks
    
    plt.grid(True, alpha=0.3)
    plt.legend(loc='upper left')
    
    # Add correlation text box
    plt.text(0.02, 0.98, f'Correlation: {correlation:.3f}', 
             transform=plt.gca().transAxes, fontsize=12,
             bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8),
             verticalalignment='top')
    
    # Save plot
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    
    print(f"  Simple dotplot saved to: {output_path}")
    
    # Calculate and display summary stats
    same_chr_pct = (len(same_chr) / len(plot_df)) * 100 if len(plot_df) > 0 else 0
    print(f"  Same chromosome: {same_chr_pct:.1f}% ({len(same_chr):,}/{len(plot_df):,})")
    
    return correlation



def linearize_real_genome_coordinates(normal_df, custom_order=None):
    """
    Linearize using REAL genomic coordinates only
    """
    linearized_df = normal_df.copy()
    
    # Group by chromosome and find real chromosome sizes
    chr_offsets = {}
    cumulative_offset = 0
    
    if custom_order:
        # Start with custom order, then add any missing chromosomes
        all_chroms = set(normal_df['sequence'].unique())
        chromosome_order = [chrom for chrom in custom_order if chrom in all_chroms]
        missing_chroms = sorted(all_chroms - set(chromosome_order))
        chromosome_order.extend(missing_chroms)
        if missing_chroms:
            print(f"  Warning: Adding missing chromosomes at end: {missing_chroms}")
    else:
        chromosome_order = sorted(normal_df['sequence'].unique())

    for chrom in chromosome_order:
        chr_offsets[chrom] = cumulative_offset
        chrom_genes = normal_df[normal_df['sequence'] == chrom]
        
        # Use REAL chromosome size (actual max position)
        chr_size = chrom_genes['gene_end'].max()
        cumulative_offset += chr_size
    
    # Add real linearized positions
    linearized_df['linearized_real_position'] = linearized_df.apply(
        lambda row: chr_offsets[row['sequence']] + row['gene_start'], axis=1
    )
    
    return linearized_df




def create_comparison_summary(genome1_df, genome2_df, genome1_name, genome2_name, correlation, output_dir):
    """Create text summary of the comparison"""
    
    summary_path = Path(output_dir) / f"{Path(output_dir).name}_comparison_summary.txt"
    
    # Find common genes
    common_genes = set(genome1_df['busco_id']) & set(genome2_df['busco_id'])
    g1_only = set(genome1_df['busco_id']) - common_genes
    g2_only = set(genome2_df['busco_id']) - common_genes
    
    summary_text = f"""
SIMPLE BUSCO DOTPLOT COMPARISON SUMMARY
{"="*60}

GENOMES COMPARED:
{genome1_name}: {len(genome1_df):,} genes, {len(genome1_df['sequence'].unique())} chromosomes
{genome2_name}: {len(genome2_df):,} genes, {len(genome2_df['sequence'].unique())} chromosomes

GENE OVERLAP:
Common genes: {len(common_genes):,}
{genome1_name} only: {len(g1_only):,}
{genome2_name} only: {len(g2_only):,}
Total unique genes: {len(common_genes | g1_only | g2_only):,}

CORRELATION ANALYSIS:
Overall correlation: {correlation:.3f}

INTERPRETATION:
• R > 0.9: Excellent synteny (very similar gene order)
• R 0.7-0.9: Good synteny (mostly similar gene order)
• R 0.3-0.7: Moderate synteny (some rearrangements)
• R < 0.3: Poor synteny (major rearrangements)

CHROMOSOME DISTRIBUTIONS:
{"-"*40}

{genome1_name}:
"""
    
    chr_counts1 = genome1_df['sequence'].value_counts().sort_index()
    for chr_name, count in chr_counts1.items():
        summary_text += f"  {chr_name}: {count:,} genes\n"
    
    summary_text += f"\n{genome2_name}:\n"
    
    chr_counts2 = genome2_df['sequence'].value_counts().sort_index()
    for chr_name, count in chr_counts2.items():
        summary_text += f"  {chr_name}: {count:,} genes\n"
    
    # Save summary
    with open(summary_path, 'w') as f:
        f.write(summary_text)
    
    print(f"\nComparison summary saved to: {summary_path}")
    return summary_path

def get_chromosome_order(df, custom_order=None):
    """Get chromosome ordering - custom if provided, otherwise default sort"""
    
    # Get all chromosomes in the dataset
    chroms = sorted(df['sequence'].unique())
    
    # If custom order provided via command line, use it
    if custom_order:
        return [chrom for chrom in custom_order if chrom in chroms]
    
    # Default: use original sorting
    return chroms

def main():
    parser = argparse.ArgumentParser(description='Simple BUSCO dotplot comparison')
    parser.add_argument('genome1', help='First BUSCO TSV file')
    parser.add_argument('genome2', help='Second BUSCO TSV file')
    parser.add_argument('output', help='Output PNG file')
    parser.add_argument('--name1', default='Genome 1', help='Name for first genome')
    parser.add_argument('--name2', default='Genome 2', help='Name for second genome')
    parser.add_argument('--chr-order1', nargs='*', help='Custom chromosome order for genome 1')
    parser.add_argument('--chr-order2', nargs='*', help='Custom chromosome order for genome 2')
    
    args = parser.parse_args()
    
    try:
        print("SIMPLE BUSCO DOTPLOT TOOL")
        print("="*50)
        
        # Load genomes
        genome1_df = load_busco_file(args.genome1, args.name1)
        genome2_df = load_busco_file(args.genome2, args.name2)
        
        # Create dotplot
        correlation = create_simple_busco_dotplot(
            genome1_df, genome2_df, 
            args.name1, args.name2, 
            args.output,
            args.chr_order1, args.chr_order2
        )
        
        if correlation is not None:
            # Create summary
            output_dir = Path(args.output).parent
            summary_path = create_comparison_summary(
                genome1_df, genome2_df,
                args.name1, args.name2,
                correlation, output_dir
            )
            
            print("\n" + "="*50)
            print("DOTPLOT COMPARISON COMPLETED")
            print("="*50)
            print(f"Correlation: {correlation:.3f}")
            print(f"Plot: {args.output}")
            print(f"Summary: {summary_path}")
            
            # Interpretation
            if correlation > 0.9:
                print("🟢 Excellent synteny - very similar gene order")
            elif correlation > 0.7:
                print("🟡 Good synteny - mostly similar gene order") 
            elif correlation > 0.3:
                print("🟠 Moderate synteny - some rearrangements")
            else:
                print("🔴 Poor synteny - major rearrangements")
        else:
            print("❌ No common genes found - comparison failed")
            return 1
        
        return 0
        
    except Exception as e:
        print(f"❌ Error: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())