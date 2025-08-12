#!/usr/bin/env python3
"""
Independent linear BUSCO dotplotter for comparing gene orders between two BUSCO TSV files.
Creates separate, clean plots instead of cluttered dashboards.
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


def create_simple_dotplot(genome1_df, genome2_df, genome1_name, genome2_name, output_path):
    """Create a clean, simple dotplot with just the essential information"""
    
    # Find common genes
    common_genes = set(genome1_df['busco_id']) & set(genome2_df['busco_id'])
    print(f"Common genes found: {len(common_genes)}")
    
    if len(common_genes) == 0:
        print("Warning: No common genes found between the two genomes!")
        return None
    
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
                'genome1_pos': g1_positions[gene],
                'genome2_pos': g2_positions[gene]
            })
    
    plot_df = pd.DataFrame(plot_data)
    
    # Calculate correlation
    correlation = 0.0
    if len(plot_df) > 1:
        correlation = plot_df['genome1_pos'].corr(plot_df['genome2_pos'])
    
    # Create clean plot
    plt.figure(figsize=(10, 8))
    
    # Scatter plot
    plt.scatter(plot_df['genome1_pos'], plot_df['genome2_pos'],
                alpha=0.6, s=50, c='steelblue', edgecolors='darkblue', linewidth=0.5)
    
    # Add diagonal for perfect correlation
    max_pos = max(plot_df['genome1_pos'].max(), plot_df['genome2_pos'].max())
    plt.plot([0, max_pos], [0, max_pos], 'r--', alpha=0.5, linewidth=2, label='Perfect correlation')
    
    # Formatting
    plt.xlabel(f'{genome1_name} Linear Position', fontsize=12)
    plt.ylabel(f'{genome2_name} Linear Position', fontsize=12)
    plt.title(f'Gene Order Comparison\n{genome1_name} vs {genome2_name}\nCorrelation: {correlation:.3f}', fontsize=14)
    plt.grid(True, alpha=0.3)
    plt.legend()
    
    # Save
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    
    print(f"Simple dotplot saved to: {output_path}")
    
    return {
        'common_genes': len(common_genes),
        'correlation': correlation,
        'total_genes_1': len(genome1_df),
        'total_genes_2': len(genome2_df),
        'chromosomes_1': len(genome1_df['sequence'].unique()),
        'chromosomes_2': len(genome2_df['sequence'].unique())
    }


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
           chr_size = chr_groups.loc[chr_name, 'max'] - chr_groups.loc[chr_name, 'min'] + 10000000  # Add 10Mb padding
           cumulative_offset += chr_size
  
   # Add linearised coordinates
   genome_linear['linear_position_mb'] = genome_linear.apply(
       lambda row: (chr_offsets.get(row['sequence'], 0) + row['start']) / 1e6, axis=1
   )
  
   return genome_linear, chr_offsets


def create_linearised_dotplot(genome1_df, genome2_df, genome1_name, genome2_name, output_path):
    """Create a clean linearised genome comparison plot"""
    
    # Find common genes
    common_genes = set(genome1_df['busco_id']) & set(genome2_df['busco_id'])
    print(f"Creating linearised plot with {len(common_genes)} common genes")
    
    if len(common_genes) == 0:
        return None
    
    # Filter to common genes only
    g1_common = genome1_df[genome1_df['busco_id'].isin(common_genes)].copy()
    g2_common = genome2_df[genome2_df['busco_id'].isin(common_genes)].copy()
    
    # Linearise coordinates for both genomes
    g1_linear, g1_offsets = linearise_genome_coordinates(g1_common)
    g2_linear, g2_offsets = linearise_genome_coordinates(g2_common)
    
    # Create mapping between genomes using linear coordinates
    g1_positions = dict(zip(g1_linear['busco_id'], g1_linear['linear_position_mb']))
    g2_positions = dict(zip(g2_linear['busco_id'], g2_linear['linear_position_mb']))
    
    # Prepare data for plotting
    plot_data = []
    for gene in common_genes:
        if gene in g1_positions and gene in g2_positions:
            plot_data.append({
                'genome1_linear': g1_positions[gene],
                'genome2_linear': g2_positions[gene]
            })
    
    plot_df = pd.DataFrame(plot_data)
    
    # Calculate correlation
    correlation = 0.0
    if len(plot_df) > 1:
        correlation = plot_df['genome1_linear'].corr(plot_df['genome2_linear'])
    
    # Create clean linearised plot
    plt.figure(figsize=(12, 10))
    
    # Scatter plot
    plt.scatter(plot_df['genome1_linear'], plot_df['genome2_linear'],
                alpha=0.6, s=30, c='red', edgecolors='darkred', linewidth=0.3)
    
    # Add chromosome boundaries
    for chr_name, offset in g1_offsets.items():
        if offset > 0:  # Don't add line at position 0
            plt.axvline(x=offset/1e6, color='grey', linewidth=0.8, alpha=0.5, linestyle='--')
    
    for chr_name, offset in g2_offsets.items():
        if offset > 0:  # Don't add line at position 0
            plt.axhline(y=offset/1e6, color='grey', linewidth=0.8, alpha=0.5, linestyle='--')
    
    # Formatting
    plt.xlabel(f'{genome1_name} Linearised Position (Mb)', fontsize=12)
    plt.ylabel(f'{genome2_name} Linearised Position (Mb)', fontsize=12)
    plt.title(f'Linearised Genome Comparison\n{genome1_name} vs {genome2_name}\nLinear R = {correlation:.3f}', fontsize=14)
    plt.grid(True, alpha=0.2)
    
    # Add chromosome labels
    xlim = plt.xlim()
    ylim = plt.ylim()
    
    # Add chromosome names as text annotations
    sorted_g1_chrs = sorted(g1_offsets.items(), key=lambda x: x[1])
    for chr_name, offset in sorted_g1_chrs:
        offset_mb = offset / 1e6
        if offset_mb < xlim[1]:
            plt.text(offset_mb + (xlim[1] * 0.01), ylim[1] * 0.98, chr_name,
                    rotation=90, ha='left', va='top', fontsize=8, alpha=0.7)
    
    sorted_g2_chrs = sorted(g2_offsets.items(), key=lambda x: x[1])
    for chr_name, offset in sorted_g2_chrs:
        offset_mb = offset / 1e6
        if offset_mb < ylim[1]:
            plt.text(xlim[1] * 0.98, offset_mb + (ylim[1] * 0.01), chr_name,
                    rotation=0, ha='right', va='bottom', fontsize=8, alpha=0.7)
    
    # Save
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    
    print(f"Linearised dotplot saved to: {output_path}")
    
    return {
        'linear_correlation': correlation,
        'total_linearised_genes': len(plot_df)
    }
    """Create a clean chromosome comparison plot"""
    
    # Find common genes
    common_genes = set(genome1_df['busco_id']) & set(genome2_df['busco_id'])
    
    if len(common_genes) == 0:
        return None
    
    # Filter to common genes only
    g1_common = genome1_df[genome1_df['busco_id'].isin(common_genes)].copy()
    g2_common = genome2_df[genome2_df['busco_id'].isin(common_genes)].copy()
    
    # Get chromosome counts
    chr_counts1 = g1_common['sequence'].value_counts().sort_index()
    chr_counts2 = g2_common['sequence'].value_counts().sort_index()
    
    # Create figure with subplots
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))
    
    # Plot 1: Genome1 chromosome distribution
    ax1.bar(range(len(chr_counts1)), chr_counts1.values, alpha=0.7, color='lightblue', edgecolor='darkblue')
    ax1.set_xticks(range(len(chr_counts1)))
    ax1.set_xticklabels(chr_counts1.index, rotation=45, ha='right')
    ax1.set_ylabel('Common Gene Count', fontsize=12)
    ax1.set_title(f'{genome1_name} - Common Genes per Chromosome', fontsize=14)
    ax1.grid(True, alpha=0.3, axis='y')
    
    # Plot 2: Genome2 chromosome distribution
    ax2.bar(range(len(chr_counts2)), chr_counts2.values, alpha=0.7, color='lightgreen', edgecolor='darkgreen')
    ax2.set_xticks(range(len(chr_counts2)))
    ax2.set_xticklabels(chr_counts2.index, rotation=45, ha='right')
    ax2.set_ylabel('Common Gene Count', fontsize=12)
    ax2.set_title(f'{genome2_name} - Common Genes per Chromosome', fontsize=14)
    ax2.grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    
    print(f"Chromosome comparison saved to: {output_path}")
    
    return True


def generate_detailed_summary(genome1_df, genome2_df, genome1_name, genome2_name, stats, output_path):
    """Generate a detailed text summary of the comparison"""
    
    if stats is None:
        return
    
    # Find common genes for detailed analysis
    common_genes = set(genome1_df['busco_id']) & set(genome2_df['busco_id'])
    g1_common = genome1_df[genome1_df['busco_id'].isin(common_genes)].copy()
    g2_common = genome2_df[genome2_df['busco_id'].isin(common_genes)].copy()
    
    # Generate comprehensive summary
    summary = f"""
BUSCO DOTPLOT COMPARISON SUMMARY
{"="*50}

DATASET INFORMATION:
{genome1_name}: {stats['total_genes_1']:,} total genes, {stats['chromosomes_1']} chromosomes
{genome2_name}: {stats['total_genes_2']:,} total genes, {stats['chromosomes_2']} chromosomes

COMPARISON RESULTS:
Common BUSCO genes: {stats['common_genes']:,}
Linear correlation: {stats['correlation']:.3f}

INTERPRETATION:
• Correlation > 0.7: High synteny conservation
• Correlation 0.3-0.7: Moderate synteny
• Correlation < 0.3: Low synteny / high rearrangement

CHROMOSOME DISTRIBUTIONS:
{"-"*30}

{genome1_name} chromosomes:
"""
    
    # Add chromosome details for genome1
    chr_counts1 = g1_common['sequence'].value_counts().sort_index()
    for chr_name, count in chr_counts1.items():
        summary += f"  {chr_name}: {count} genes\n"
    
    summary += f"\n{genome2_name} chromosomes:\n"
    
    # Add chromosome details for genome2
    chr_counts2 = g2_common['sequence'].value_counts().sort_index()
    for chr_name, count in chr_counts2.items():
        summary += f"  {chr_name}: {count} genes\n"
    
    summary += f"\nCHROMOSOME OVERLAP ANALYSIS:\n{'-'*30}\n"
    
    # Calculate chromosome-to-chromosome overlaps
    overlap_data = []
    for chr1 in sorted(g1_common['sequence'].unique()):
        chr1_genes = set(g1_common[g1_common['sequence'] == chr1]['busco_id'])
        
        for chr2 in sorted(g2_common['sequence'].unique()):
            chr2_genes = set(g2_common[g2_common['sequence'] == chr2]['busco_id'])
            overlap = len(chr1_genes & chr2_genes)
            
            if overlap > 0:
                overlap_data.append((chr1, chr2, overlap))
    
    # Sort by overlap count and show top matches
    overlap_data.sort(key=lambda x: x[2], reverse=True)
    
    for chr1, chr2, overlap in overlap_data[:20]:  # Show top 20
        summary += f"{chr1} ↔ {chr2}: {overlap} genes\n"
    
    if len(overlap_data) > 20:
        summary += f"... and {len(overlap_data) - 20} more chromosome pairs\n"
    
    # Save summary to file
    with open(output_path, 'w') as f:
        f.write(summary)
    
    print(f"Detailed summary saved to: {output_path}")


def main():
   parser = argparse.ArgumentParser(description='Create clean, separate BUSCO dotplot comparisons')
   parser.add_argument('genome1_tsv', help='First genome BUSCO TSV file')
   parser.add_argument('genome2_tsv', help='Second genome BUSCO TSV file')
   parser.add_argument('output_prefix', help='Output file prefix (will create multiple files)')
   parser.add_argument('--name1', default='Genome1', help='Name for first genome')
   parser.add_argument('--name2', default='Genome2', help='Name for second genome')
  
   args = parser.parse_args()
  
   try:
       # Load genomes
       genome1_df = load_busco_tsv(args.genome1_tsv, args.name1)
       genome2_df = load_busco_tsv(args.genome2_tsv, args.name2)
      
       # Create output paths
       base_path = Path(args.output_prefix)
       simple_dotplot = f"{base_path}_dotplot.png"
       linearised_plot = f"{base_path}_linearised.png"
       chromosome_plot = f"{base_path}_chromosomes.png"
       summary_file = f"{base_path}_summary.txt"
       
       print("\nCreating plots...")
       
       # Create simple dotplot
       stats = create_simple_dotplot(genome1_df, genome2_df, args.name1, args.name2, simple_dotplot)
       
       if stats is not None:
           # Create linearised dotplot
           linear_stats = create_linearised_dotplot(genome1_df, genome2_df, args.name1, args.name2, linearised_plot)
           
           # Create chromosome comparison
           create_chromosome_comparison(genome1_df, genome2_df, args.name1, args.name2, chromosome_plot)
           
           # Generate detailed summary
           generate_detailed_summary(genome1_df, genome2_df, args.name1, args.name2, stats, summary_file)
           
           print("\n" + "="*50)
           print("COMPARISON COMPLETE")
           print("="*50)
           print(f"Common genes: {stats['common_genes']}")
           print(f"Position correlation: {stats['correlation']:.3f}")
           if linear_stats:
               print(f"Linear correlation: {linear_stats['linear_correlation']:.3f}")
           print(f"\nFiles created:")
           print(f"  Dotplot: {simple_dotplot}")
           print(f"  Linearised: {linearised_plot}")
           print(f"  Chromosomes: {chromosome_plot}")
           print(f"  Summary: {summary_file}")
       else:
           print("❌ No common genes found - comparison cannot be completed")
      
   except Exception as e:
       print(f"❌ Error: {e}")
       return 1
  
   return 0


if __name__ == "__main__":
   import sys
   sys.exit(main())