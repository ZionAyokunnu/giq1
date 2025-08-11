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


def create_linear_dotplot(genome1_df, genome2_df, genome1_name, genome2_name, output_path):
    """Create linear dotplot comparing two genomes"""
    
    # Find common genes
    common_genes = set(genome1_df['busco_id']) & set(genome2_df['busco_id'])
    print(f"Common genes found: {len(common_genes)}")
    
    if len(common_genes) == 0:
        print("Warning: No common genes found between the two genomes!")
        return
    
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
    
    ax_stats.text(0.05, 0.95, stats_text, transform=ax_stats.transAxes, 
                 fontsize=9, verticalalignment='top', fontfamily='monospace')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Dotplot saved to: {output_path}")
    
    # Return summary statistics
    return {
        'total_genes_genome1': len(genome1_df),
        'total_genes_genome2': len(genome2_df),
        'common_genes': len(common_genes),
        'correlation': correlation,
        'chromosomes_genome1': len(genome1_df['sequence'].unique()),
        'chromosomes_genome2': len(genome2_df['sequence'].unique())
    }


def main():
    parser = argparse.ArgumentParser(description='Create linear BUSCO dotplot comparison')
    parser.add_argument('genome1_tsv', help='First genome BUSCO TSV file')
    parser.add_argument('genome2_tsv', help='Second genome BUSCO TSV file')
    parser.add_argument('output_plot', help='Output plot file (PNG/PDF)')
    parser.add_argument('--name1', default='Genome1', help='Name for first genome')
    parser.add_argument('--name2', default='Genome2', help='Name for second genome')
    
    args = parser.parse_args()
    
    try:
        # Load genomes
        genome1_df = load_busco_tsv(args.genome1_tsv, args.name1)
        genome2_df = load_busco_tsv(args.genome2_tsv, args.name2)
        
        # Create dotplot
        stats = create_linear_dotplot(genome1_df, genome2_df, 
                                    args.name1, args.name2, 
                                    args.output_plot)
        
        print("\n" + "="*50)
        print("DOTPLOT COMPARISON COMPLETE")
        print("="*50)
        print(f"Common genes: {stats['common_genes']}")
        print(f"Linear correlation: {stats['correlation']:.3f}")
        print(f"Plot saved to: {args.output_plot}")
        
    except Exception as e:
        print(f"❌ Error: {e}")
        return 1
    
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())


