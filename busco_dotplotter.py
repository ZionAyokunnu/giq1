#!/usr/bin/env python3
"""
Enhanced BUSCO dotplotter with chromosome content matching and three-color system.
Designed for AGORA (single linear chromosome) vs GIQ (multiple consolidated chromosomes) comparison.


script:
python3 enhanced_dotplotter.py agora.tsv giq_consolidated.tsv results \
    --window-size 75 --threshold 0.7 --agora-name "AGORA_Ancestral_genome" --giq-name "GIQ_ancestral_genome"

"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse
from pathlib import Path
import seaborn as sns
from collections import Counter, defaultdict


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
        
        # Add 'start' column for compatibility
        df['start'] = df['gene_start']
        
        # Sort by position within each chromosome
        df = df.sort_values(['sequence', 'gene_start'])
        
        # Add linear position for plotting
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


def segment_agora_by_giq_dominance(agora_df, giq_df, window_size=50, threshold=0.6, 
                                  output_segments_tsv=None, output_matched_tsv=None):
    """
    Segment AGORA's linear gene order based on GIQ chromosome dominance.
    
    Args:
        agora_df: AGORA BUSCO DataFrame (single chromosome)
        giq_df: GIQ BUSCO DataFrame (multiple chromosomes) 
        window_size: Number of consecutive AGORA genes to analyze
        threshold: Minimum percentage for confident assignment
        output_segments_tsv: Optional path to save segment analysis
        output_matched_tsv: Optional path to save final matched AGORA TSV
    
    Returns:
        matched_agora_df: AGORA DataFrame with updated chromosome assignments
        segments_df: DataFrame with segmentation details
    """
    
    print(f"\nSegmenting AGORA by GIQ chromosome dominance:")
    print(f"  Window size: {window_size} genes")
    print(f"  Confidence threshold: {threshold*100}%")
    
    # Create GIQ gene-to-chromosome mapping
    giq_gene_mapping = dict(zip(giq_df['busco_id'], giq_df['sequence']))
    giq_chromosomes = sorted(giq_df['sequence'].unique())
    
    print(f"  GIQ chromosomes: {giq_chromosomes}")
    print(f"  GIQ gene mapping: {len(giq_gene_mapping)} genes")
    
    # Get AGORA gene order (sorted by linear position)
    agora_ordered = agora_df.sort_values('linear_position').copy()
    agora_genes = agora_ordered['busco_id'].tolist()
    
    print(f"  AGORA gene order: {len(agora_genes)} genes")
    
    # Analyze windows for chromosome dominance
    segments = []
    
    for start_idx in range(0, len(agora_genes), window_size):
        end_idx = min(start_idx + window_size, len(agora_genes))
        window_genes = agora_genes[start_idx:end_idx]
        
        # Count genes per GIQ chromosome in this window
        chr_counts = Counter()
        mapped_genes = 0
        
        for gene in window_genes:
            if gene in giq_gene_mapping:
                chr_counts[giq_gene_mapping[gene]] += 1
                mapped_genes += 1
        
        # Find dominant chromosome
        if chr_counts:
            dominant_chr = chr_counts.most_common(1)[0][0]
            dominant_count = chr_counts[dominant_chr]
            confidence = dominant_count / len(window_genes)
            mapped_percentage = mapped_genes / len(window_genes)
        else:
            dominant_chr = 'chr_unassigned'
            dominant_count = 0
            confidence = 0.0
            mapped_percentage = 0.0
        
        # Assignment decision
        if confidence >= threshold:
            assigned_chr = f"{dominant_chr}_agora"
            assignment_type = 'confident'
        elif chr_counts:
            assigned_chr = f"{dominant_chr}_agora"
            assignment_type = 'best_guess'
        else:
            assigned_chr = 'chr_unassigned_agora'
            assignment_type = 'unassigned'
        
        segments.append({
            'segment_id': len(segments) + 1,
            'start_gene_idx': start_idx,
            'end_gene_idx': end_idx,
            'window_size': len(window_genes),
            'mapped_genes': mapped_genes,
            'mapped_percentage': mapped_percentage,
            'dominant_giq_chr': dominant_chr,
            'dominant_count': dominant_count,
            'confidence': confidence,
            'assigned_chr': assigned_chr,
            'assignment_type': assignment_type,
            'all_chr_counts': dict(chr_counts)
        })
        
        print(f"  Segment {len(segments):2d}: genes {start_idx:4d}-{end_idx:4d} → {assigned_chr} ({confidence:.2%})")
    
    segments_df = pd.DataFrame(segments)
    
    # Apply chromosome assignments to AGORA genes
    matched_agora_df = agora_ordered.copy()
    matched_agora_df['original_sequence'] = matched_agora_df['sequence']
    
    for _, segment in segments_df.iterrows():
        start_idx = segment['start_gene_idx'] 
        end_idx = segment['end_gene_idx']
        assigned_chr = segment['assigned_chr']
        
        # Update chromosome assignment for genes in this segment
        mask = (matched_agora_df.index >= agora_ordered.index[start_idx]) & \
               (matched_agora_df.index <= agora_ordered.index[end_idx-1])
        matched_agora_df.loc[mask, 'sequence'] = assigned_chr
    
    # Print summary
    print(f"\nSegmentation Summary:")
    assignment_counts = segments_df['assignment_type'].value_counts()
    for assignment_type, count in assignment_counts.items():
        print(f"  {assignment_type}: {count} segments")
    
    chr_assignments = matched_agora_df['sequence'].value_counts()
    print(f"\nFinal AGORA chromosome assignments:")
    for chr_name, count in chr_assignments.items():
        print(f"  {chr_name}: {count} genes")
    
    # Save intermediate files
    if output_segments_tsv:
        segments_df.to_csv(output_segments_tsv, sep='\t', index=False)
        print(f"\nSegmentation details saved to: {output_segments_tsv}")
    
    if output_matched_tsv:
        matched_agora_df.to_csv(output_matched_tsv, sep='\t', index=False)
        print(f"Matched AGORA TSV saved to: {output_matched_tsv}")
    
    return matched_agora_df, segments_df


def create_three_color_chromosome_dotplot(agora_df, giq_df, agora_name, giq_name, output_path):
    """
    Create three-color dotplot showing:
    - Blue: GIQ-only genes
    - Red: AGORA-only genes  
    - Green: Collision genes (both methods agree on same chromosome)
    """
    
    print(f"\nCreating three-color chromosome dotplot...")
    
    # Find all genes in both datasets
    agora_genes = set(agora_df['busco_id'])
    giq_genes = set(giq_df['busco_id'])
    
    # Create chromosome mappings
    agora_chr_map = dict(zip(agora_df['busco_id'], agora_df['sequence']))
    giq_chr_map = dict(zip(giq_df['busco_id'], giq_df['sequence']))
    
    agora_pos_map = dict(zip(agora_df['busco_id'], agora_df['linear_position']))
    giq_pos_map = dict(zip(giq_df['busco_id'], giq_df['linear_position']))
    
    # Categorize genes
    common_genes = agora_genes & giq_genes
    agora_only = agora_genes - giq_genes
    giq_only = giq_genes - agora_genes
    
    print(f"  Common genes: {len(common_genes)}")
    print(f"  AGORA-only genes: {len(agora_only)}")
    print(f"  GIQ-only genes: {len(giq_only)}")
    
    # For common genes, check if they're on matching chromosomes
    collision_genes = []
    mismatch_genes = []
    
    for gene in common_genes:
        agora_chr = agora_chr_map[gene]
        giq_chr = giq_chr_map[gene]
        
        # Check if chromosomes match (removing _agora/_giq suffixes)
        agora_chr_base = agora_chr.replace('_agora', '')
        giq_chr_base = giq_chr.replace('_giq', '') if '_giq' in giq_chr else giq_chr
        
        if agora_chr_base == giq_chr_base:
            collision_genes.append(gene)
        else:
            mismatch_genes.append(gene)
    
    print(f"  Collision genes (same chromosome): {len(collision_genes)}")
    print(f"  Mismatch genes (different chromosomes): {len(mismatch_genes)}")
    
    # Prepare plotting data
    plot_data = {
        'agora_pos': [],
        'giq_pos': [],
        'color': [],
        'gene_type': []
    }
    
    # Add collision genes (green)
    for gene in collision_genes:
        plot_data['agora_pos'].append(agora_pos_map[gene])
        plot_data['giq_pos'].append(giq_pos_map[gene])
        plot_data['color'].append('green')
        plot_data['gene_type'].append('collision')
    
    # Add mismatch genes (purple - show disagreement)
    for gene in mismatch_genes:
        plot_data['agora_pos'].append(agora_pos_map[gene])
        plot_data['giq_pos'].append(giq_pos_map[gene])
        plot_data['color'].append('purple')
        plot_data['gene_type'].append('mismatch')
    
    plot_df = pd.DataFrame(plot_data)
    
    # Create plot
    plt.figure(figsize=(12, 10))
    
    # Plot each category
    if len(collision_genes) > 0:
        collision_data = plot_df[plot_df['gene_type'] == 'collision']
        plt.scatter(collision_data['agora_pos'], collision_data['giq_pos'], 
                   c='green', alpha=0.7, s=30, label=f'Same chromosome ({len(collision_genes)})', 
                   edgecolors='darkgreen', linewidth=0.3)
    
    if len(mismatch_genes) > 0:
        mismatch_data = plot_df[plot_df['gene_type'] == 'mismatch']
        plt.scatter(mismatch_data['agora_pos'], mismatch_data['giq_pos'],
                   c='purple', alpha=0.7, s=30, label=f'Different chromosomes ({len(mismatch_genes)})', 
                   edgecolors='darkmagenta', linewidth=0.3)
    
    # Add chromosome boundaries and labels for AGORA
    agora_chr_boundaries = {}
    cumulative_pos = 0
    for chrom in sorted(agora_df['sequence'].unique()):
        chrom_genes = agora_df[agora_df['sequence'] == chrom]
        if len(chrom_genes) > 0:
            agora_chr_boundaries[chrom] = (cumulative_pos, cumulative_pos + len(chrom_genes))
            cumulative_pos += len(chrom_genes)
    
    # Add vertical lines for AGORA chromosome boundaries  
    for chrom, (start, end) in agora_chr_boundaries.items():
        if start > 0:
            plt.axvline(x=start, color='red', linewidth=1, alpha=0.5, linestyle='--')
        # Add chromosome labels
        mid_pos = (start + end) / 2
        plt.text(mid_pos, plt.ylim()[1] * 0.95, chrom.replace('_agora', ''), 
                rotation=45, ha='center', va='top', fontsize=8, color='red', alpha=0.8)
    
    # Add chromosome boundaries and labels for GIQ
    giq_chr_boundaries = {}
    cumulative_pos = 0
    for chrom in sorted(giq_df['sequence'].unique()):
        chrom_genes = giq_df[giq_df['sequence'] == chrom]
        if len(chrom_genes) > 0:
            giq_chr_boundaries[chrom] = (cumulative_pos, cumulative_pos + len(chrom_genes))
            cumulative_pos += len(chrom_genes)
    
    # Add horizontal lines for GIQ chromosome boundaries
    for chrom, (start, end) in giq_chr_boundaries.items():
        if start > 0:
            plt.axhline(y=start, color='blue', linewidth=1, alpha=0.5, linestyle='--')
        # Add chromosome labels
        mid_pos = (start + end) / 2
        plt.text(plt.xlim()[0] * 0.02, mid_pos, chrom, 
                rotation=0, ha='left', va='center', fontsize=8, color='blue', alpha=0.8)
    
    # Formatting
    plt.xlabel(f'{agora_name} Linear Gene Position', fontsize=12)
    plt.ylabel(f'{giq_name} Linear Gene Position', fontsize=12)
    plt.title(f'Three-Color Chromosome Synteny Analysis\n{agora_name} vs {giq_name}', fontsize=14)
    plt.grid(True, alpha=0.3)
    plt.legend(loc='upper right')
    
    # Calculate and display correlation
    if len(plot_df) > 1:
        correlation = plot_df['agora_pos'].corr(plot_df['giq_pos'])
        plt.text(0.02, 0.98, f'Correlation: {correlation:.3f}', 
                transform=plt.gca().transAxes, fontsize=10, 
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    
    print(f"Three-color dotplot saved to: {output_path}")
    
    # Return summary statistics
    return {
        'total_genes': len(agora_genes | giq_genes),
        'common_genes': len(common_genes),
        'agora_only': len(agora_only),
        'giq_only': len(giq_only),
        'collision_genes': len(collision_genes),
        'mismatch_genes': len(mismatch_genes),
        'synteny_agreement': len(collision_genes) / len(common_genes) if common_genes else 0
    }


def create_method_comparison_summary(agora_df, giq_df, segments_df, stats, agora_name, giq_name, output_path):
    """Create comprehensive summary of method comparison"""
    
    summary_text = f"""
AGORA vs GIQ METHOD COMPARISON SUMMARY
{"="*60}

DATASET INFORMATION:
{agora_name}: {len(agora_df):,} total genes, {len(agora_df['sequence'].unique())} chromosomes
{giq_name}: {len(giq_df):,} total genes, {len(giq_df['sequence'].unique())} chromosomes

SEGMENTATION RESULTS:
{"="*40}
Window size: 50 genes
Total segments: {len(segments_df)}

Segment confidence:
"""
    
    # Add segment confidence breakdown
    assignment_counts = segments_df['assignment_type'].value_counts()
    for assignment_type, count in assignment_counts.items():
        percentage = (count / len(segments_df)) * 100
        summary_text += f"  {assignment_type}: {count} segments ({percentage:.1f}%)\n"
    
    summary_text += f"\nCHROMOSOME ASSIGNMENTS:\n{'-'*40}\n"
    
    # Add AGORA chromosome assignments
    agora_chr_counts = agora_df['sequence'].value_counts()
    for chr_name, count in agora_chr_counts.items():
        summary_text += f"{chr_name}: {count} genes\n"
    
    summary_text += f"\nGIQ CHROMOSOMES:\n{'-'*40}\n"
    
    # Add GIQ chromosome counts
    giq_chr_counts = giq_df['sequence'].value_counts()
    for chr_name, count in giq_chr_counts.items():
        summary_text += f"{chr_name}: {count} genes\n"
    
    summary_text += f"""
SYNTENY ANALYSIS:
{"-"*40}
Total unique genes: {stats['total_genes']:,}
Common genes (both methods): {stats['common_genes']:,}
AGORA-only genes: {stats['agora_only']:,}
GIQ-only genes: {stats['giq_only']:,}

CHROMOSOME AGREEMENT:
Same chromosome assignment: {stats['collision_genes']:,} genes
Different chromosome assignment: {stats['mismatch_genes']:,} genes
Synteny agreement rate: {stats['synteny_agreement']:.1%}

INTERPRETATION:
• >80% agreement: High synteny conservation
• 60-80% agreement: Moderate synteny conservation  
• <60% agreement: High chromosomal rearrangement/uncertainty

CONFIDENCE ANALYSIS:
High confidence segments (≥60%): {assignment_counts.get('confident', 0)} 
Best guess segments (<60%): {assignment_counts.get('best_guess', 0)}
Unassigned segments: {assignment_counts.get('unassigned', 0)}

The segmentation approach successfully identified chromosome blocks in AGORA's 
linear gene order based on GIQ chromosome content dominance.
"""
    
    # Create text-only image
    fig, ax = plt.subplots(figsize=(12, 14))
    ax.axis('off')
    ax.text(0.05, 0.95, summary_text, transform=ax.transAxes, fontsize=10,
            verticalalignment='top', fontfamily='monospace')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    
    print(f"Method comparison summary saved to: {output_path}")


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


def create_main_linear_dotplot(genome1_df, genome2_df, genome1_name, genome2_name, output_path):
    """1. Main Linear Dotplot - Blue scatter plot with perfect correlation diagonal line"""
    
    # Find common genes
    common_genes = set(genome1_df['busco_id']) & set(genome2_df['busco_id'])
    
    if len(common_genes) == 0:
        return None
    
    # Filter to common genes only
    g1_common = genome1_df[genome1_df['busco_id'].isin(common_genes)].copy()
    g2_common = genome2_df[genome2_df['busco_id'].isin(common_genes)].copy()
    
    # Create mapping between genomes using linear_position
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
                alpha=0.6, s=50, c='blue', edgecolors='darkblue', linewidth=0.5)
    
    # Add diagonal for perfect correlation
    max_pos = max(plot_df['genome1_pos'].max(), plot_df['genome2_pos'].max())
    plt.plot([0, max_pos], [0, max_pos], 'r--', alpha=0.5, linewidth=2, label='Perfect correlation')
    
    # Formatting
    plt.xlabel(f'{genome1_name} Linear Position', fontsize=12)
    plt.ylabel(f'{genome2_name} Linear Position', fontsize=12)
    plt.title(f'Main Linear Gene Order Comparison\n{genome1_name} vs {genome2_name}\nCorrelation: {correlation:.3f}', fontsize=14)
    plt.grid(True, alpha=0.3)
    plt.legend()
    
    # Save
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    
    print(f"Main linear dotplot saved to: {output_path}")
    return correlation


def create_genome1_chromosome_distribution(genome1_df, genome2_df, genome1_name, genome2_name, output_path):
    """2. Genome1 Chromosome Distribution - Light blue bar chart"""
    
    # Find common genes
    common_genes = set(genome1_df['busco_id']) & set(genome2_df['busco_id'])
    
    if len(common_genes) == 0:
        return None
    
    # Filter to common genes only
    g1_common = genome1_df[genome1_df['busco_id'].isin(common_genes)].copy()
    
    # Get chromosome counts
    chr_counts1 = g1_common['sequence'].value_counts().sort_index()
    
    # Create figure
    plt.figure(figsize=(12, 8))
    
    # Bar chart
    plt.bar(range(len(chr_counts1)), chr_counts1.values, alpha=0.7, color='lightblue', edgecolor='darkblue')
    plt.xticks(range(len(chr_counts1)), chr_counts1.index, rotation=45, ha='right')
    plt.ylabel('Common Gene Count', fontsize=12)
    plt.title(f'{genome1_name} - Common Genes per Chromosome', fontsize=14)
    plt.grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    
    print(f"Genome1 chromosome distribution saved to: {output_path}")
    return True


def create_genome2_chromosome_distribution(genome1_df, genome2_df, genome1_name, genome2_name, output_path):
    """3. Genome2 Chromosome Distribution - Light green bar chart"""
    
    # Find common genes
    common_genes = set(genome1_df['busco_id']) & set(genome2_df['busco_id'])
    
    if len(common_genes) == 0:
        return None
    
    # Filter to common genes only
    g2_common = genome2_df[genome2_df['busco_id'].isin(common_genes)].copy()
    
    # Get chromosome counts
    chr_counts2 = g2_common['sequence'].value_counts().sort_index()
    
    # Create figure
    plt.figure(figsize=(12, 8))
    
    # Bar chart
    plt.bar(range(len(chr_counts2)), chr_counts2.values, alpha=0.7, color='lightgreen', edgecolor='darkgreen')
    plt.xticks(range(len(chr_counts2)), chr_counts2.index, rotation=45, ha='right')
    plt.ylabel('Common Gene Count', fontsize=12)
    plt.title(f'{genome2_name} - Common Genes per Chromosome', fontsize=14)
    plt.grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    
    print(f"Genome2 chromosome distribution saved to: {output_path}")
    return True


def create_all_summaries_panel(genome1_df, genome2_df, genome1_name, genome2_name, correlation, output_path):
    """4. All Summary Information - Combined statistics and chromosome overlap analysis"""
    
    # Find common genes for detailed analysis
    common_genes = set(genome1_df['busco_id']) & set(genome2_df['busco_id'])
    g1_common = genome1_df[genome1_df['busco_id'].isin(common_genes)].copy()
    g2_common = genome2_df[genome2_df['busco_id'].isin(common_genes)].copy()
    
    # Calculate chromosome-to-chromosome overlaps
    overlap_data = []
    for chr1 in sorted(g1_common['sequence'].unique()):
        chr1_genes = set(g1_common[g1_common['sequence'] == chr1]['busco_id'])
        
        for chr2 in sorted(g2_common['sequence'].unique()):
            chr2_genes = set(g2_common[g2_common['sequence'] == chr2]['busco_id'])
            overlap = len(chr1_genes & chr2_genes)
            
            if overlap > 0:
                overlap_data.append((chr1, chr2, overlap))
    
    # Sort by overlap count
    overlap_data.sort(key=lambda x: x[2], reverse=True)
    
    # Generate comprehensive summary text
    summary_text = f"""
BUSCO DOTPLOT COMPARISON SUMMARY
{"="*60}

DATASET INFORMATION:
{genome1_name}: {len(genome1_df):,} total genes, {len(genome1_df['sequence'].unique())} chromosomes
{genome2_name}: {len(genome2_df):,} total genes, {len(genome2_df['sequence'].unique())} chromosomes

COMPARISON RESULTS:
Common BUSCO genes: {len(common_genes):,}
Linear correlation: {correlation:.3f}

INTERPRETATION:
• Correlation > 0.7: High synteny conservation
• Correlation 0.3-0.7: Moderate synteny  
• Correlation < 0.3: Low synteny / high rearrangement

CHROMOSOME DISTRIBUTIONS:
{"-"*40}

{genome1_name} chromosomes:
"""
    
    # Add chromosome details for genome1
    chr_counts1 = g1_common['sequence'].value_counts().sort_index()
    for chr_name, count in chr_counts1.items():
        summary_text += f"  {chr_name}: {count} genes\n"
    
    summary_text += f"\n{genome2_name} chromosomes:\n"
    
    # Add chromosome details for genome2
    chr_counts2 = g2_common['sequence'].value_counts().sort_index()
    for chr_name, count in chr_counts2.items():
        summary_text += f"  {chr_name}: {count} genes\n"
    
    summary_text += f"\nCHROMOSOME OVERLAP ANALYSIS:\n{'-'*40}\n"
    
    # Add top chromosome overlaps
    for chr1, chr2, overlap in overlap_data[:30]:  # Show top 30
        summary_text += f"{chr1} ↔ {chr2}: {overlap} genes\n"
    
    if len(overlap_data) > 30:
        summary_text += f"... and {len(overlap_data) - 30} more chromosome pairs\n"
    
    # Create text-only image
    fig, ax = plt.subplots(figsize=(14, 16))
    ax.axis('off')
    ax.text(0.05, 0.95, summary_text, transform=ax.transAxes, fontsize=10, 
            verticalalignment='top', fontfamily='monospace', wrap=True)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    
    print(f"All summaries panel saved to: {output_path}")
    return True


def create_standard_gene_position_comparison(genome1_df, genome2_df, genome1_name, genome2_name, output_path):
    """5. Standard Gene Position Comparison - Blue dotplot using gene_start positions"""
    
    # Find common genes
    common_genes = set(genome1_df['busco_id']) & set(genome2_df['busco_id'])
    
    if len(common_genes) == 0:
        return None
    
    # Filter to common genes only
    g1_common = genome1_df[genome1_df['busco_id'].isin(common_genes)].copy()
    g2_common = genome2_df[genome2_df['busco_id'].isin(common_genes)].copy()
    
    # Create mapping between genomes using gene_start positions
    g1_positions = dict(zip(g1_common['busco_id'], g1_common['gene_start']))
    g2_positions = dict(zip(g2_common['busco_id'], g2_common['gene_start']))
    
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
                alpha=0.7, s=40, c='blue', edgecolors='darkblue', linewidth=0.5)
    
    # Formatting
    plt.xlabel(f'{genome1_name} Gene Start Position', fontsize=12)
    plt.ylabel(f'{genome2_name} Gene Start Position', fontsize=12)
    plt.title(f'Standard Gene Position Comparison\n{genome1_name} vs {genome2_name}\nR = {correlation:.3f}', fontsize=14)
    plt.grid(True, alpha=0.3)
    
    # Save
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    
    print(f"Standard gene position comparison saved to: {output_path}")
    return correlation


def create_linearised_genome_comparison(genome1_df, genome2_df, genome1_name, genome2_name, output_path):
    """6. Linearised Genome Comparison - Red scatter plot with chromosome boundaries and labels in Mb"""
    
    # Find common genes
    common_genes = set(genome1_df['busco_id']) & set(genome2_df['busco_id'])
    
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
    
    print(f"Linearised genome comparison saved to: {output_path}")
    return correlation


def create_side_by_side_chromosome_comparison(genome1_df, genome2_df, genome1_name, genome2_name, output_path):
    """7. Side-by-side Chromosome Comparison - Both genomes on same chart"""
    
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
    
    # Create side-by-side bar plot
    fig, ax = plt.subplots(figsize=(14, 8))
    
    x = np.arange(len(chr_counts1))
    width = 0.35
    
    ax.bar(x - width/2, chr_counts1.values, width, label=genome1_name, alpha=0.7, color='lightblue', edgecolor='darkblue')
    
    # Align genome2 chromosomes if they match
    chr_counts2_aligned = []
    for chr_name in chr_counts1.index:
        if chr_name in chr_counts2.index:
            chr_counts2_aligned.append(chr_counts2[chr_name])
        else:
            chr_counts2_aligned.append(0)
    
    ax.bar(x + width/2, chr_counts2_aligned, width, label=genome2_name, alpha=0.7, color='lightgreen', edgecolor='darkgreen')
    
    ax.set_xlabel('Chromosome', fontsize=12)
    ax.set_ylabel('Common Gene Count', fontsize=12)
    ax.set_title(f'Side-by-Side Chromosome Comparison\n{genome1_name} vs {genome2_name}', fontsize=14)
    ax.set_xticks(x)
    ax.set_xticklabels(chr_counts1.index, rotation=45, ha='right')
    ax.legend()
    ax.grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    
    print(f"Side-by-side chromosome comparison saved to: {output_path}")
    return True


def main():
    parser = argparse.ArgumentParser(description='Enhanced BUSCO dotplotter with chromosome matching and all visualization types')
    parser.add_argument('agora_tsv', help='AGORA BUSCO TSV file (single linear chromosome)')
    parser.add_argument('giq_tsv', help='GIQ BUSCO TSV file (multiple consolidated chromosomes)')
    parser.add_argument('output_prefix', help='Output file prefix (will create multiple files)')
    parser.add_argument('--agora-name', default='AGORA', help='Name for AGORA method')
    parser.add_argument('--giq-name', default='GIQ', help='Name for GIQ method')
    parser.add_argument('--window-size', type=int, default=50, help='Window size for chromosome segmentation (default: 50)')
    parser.add_argument('--threshold', type=float, default=0.6, help='Confidence threshold for chromosome assignment (default: 0.6)')
    
    args = parser.parse_args()
    
    try:
        print("="*60)
        print("ENHANCED BUSCO DOTPLOTTER WITH CHROMOSOME MATCHING")
        print("="*60)
        
        # Load datasets
        agora_df = load_busco_tsv(args.agora_tsv, args.agora_name)
        giq_df = load_busco_tsv(args.giq_tsv, args.giq_name)
        
        # Create output paths for all visualizations
        base_path = Path(args.output_prefix)
        
        # TSV outputs
        segments_tsv = f"{base_path}_segmentation_details.tsv"
        matched_agora_tsv = f"{base_path}_agora_matched.tsv"
        
        # New enhanced plots
        dotplot_png = f"{base_path}_three_color_dotplot.png"
        summary_png = f"{base_path}_method_comparison_summary.png"
        
        # Original 7 plots
        output1 = f"{base_path}_main_linear_dotplot.png"
        output2 = f"{base_path}_agora_chromosomes.png"
        output3 = f"{base_path}_giq_chromosomes.png"
        output4 = f"{base_path}_all_summaries.png"
        output5 = f"{base_path}_standard_comparison.png"
        output6 = f"{base_path}_linearised_comparison.png"
        output7 = f"{base_path}_side_by_side_chromosomes.png"
        
        # Step 1: Segment AGORA by GIQ chromosome dominance
        matched_agora_df, segments_df = segment_agora_by_giq_dominance(
            agora_df, giq_df, 
            window_size=args.window_size, 
            threshold=args.threshold,
            output_segments_tsv=segments_tsv,
            output_matched_tsv=matched_agora_tsv
        )
        
        print(f"\nCreating all 9 visualization plots...")
        
        # Step 2: Create NEW enhanced plots
        stats = create_three_color_chromosome_dotplot(
            matched_agora_df, giq_df, 
            args.agora_name, args.giq_name, 
            dotplot_png
        )
        
        create_method_comparison_summary(
            matched_agora_df, giq_df, segments_df, stats,
            args.agora_name, args.giq_name,
            summary_png
        )
        
        # Step 3: Create ORIGINAL 7 plots using matched AGORA data
        correlation1 = create_main_linear_dotplot(matched_agora_df, giq_df, args.agora_name, args.giq_name, output1)
        
        if correlation1 is not None:
            create_genome1_chromosome_distribution(matched_agora_df, giq_df, args.agora_name, args.giq_name, output2)
            create_genome2_chromosome_distribution(matched_agora_df, giq_df, args.agora_name, args.giq_name, output3)
            create_all_summaries_panel(matched_agora_df, giq_df, args.agora_name, args.giq_name, correlation1, output4)
            correlation5 = create_standard_gene_position_comparison(matched_agora_df, giq_df, args.agora_name, args.giq_name, output5)
            correlation6 = create_linearised_genome_comparison(matched_agora_df, giq_df, args.agora_name, args.giq_name, output6)
            create_side_by_side_chromosome_comparison(matched_agora_df, giq_df, args.agora_name, args.giq_name, output7)
        
        print("\n" + "="*60)
        print("COMPLETE ENHANCED DOTPLOT ANALYSIS FINISHED")
        print("="*60)
        print(f"Window size used: {args.window_size} genes")
        print(f"Confidence threshold: {args.threshold*100:.0f}%")
        print(f"Synteny agreement rate: {stats['synteny_agreement']:.1%}")
        
        print(f"\nFiles created:")
        print(f"  TSV OUTPUTS:")
        print(f"    1. Segmentation details: {segments_tsv}")
        print(f"    2. Matched AGORA TSV: {matched_agora_tsv}")
        print(f"  NEW ENHANCED PLOTS:")
        print(f"    3. Three-color dotplot: {dotplot_png}")
        print(f"    4. Method comparison: {summary_png}")
        print(f"  ORIGINAL 7 PLOTS:")
        print(f"    5. Main linear dotplot: {output1}")
        print(f"    6. AGORA chromosomes: {output2}")
        print(f"    7. GIQ chromosomes: {output3}")
        print(f"    8. All summaries: {output4}")
        print(f"    9. Standard comparison: {output5}")
        print(f"    10. Linearised comparison: {output6}")
        print(f"    11. Side-by-side chromosomes: {output7}")
        
        if correlation1: 
            print(f"\nCorrelation Results:")
            print(f"  Main linear correlation: {correlation1:.3f}")
            if correlation5: print(f"  Standard position correlation: {correlation5:.3f}")
            if correlation6: print(f"  Linearised correlation: {correlation6:.3f}")        # Step 1: Segment AGORA by GIQ chromosome dominance
        matched_agora_df, segments_df = segment_agora_by_giq_dominance(
            agora_df, giq_df, 
            window_size=args.window_size, 
            threshold=args.threshold,
            output_segments_tsv=segments_tsv,
            output_matched_tsv=matched_agora_tsv
        )
        
        print(f"\nCreating all 9 visualization plots...")
        
        # Step 2: Create NEW enhanced plots
        stats = create_three_color_chromosome_dotplot(
            matched_agora_df, giq_df, 
            args.agora_name, args.giq_name, 
            dotplot_png
        )
        
        create_method_comparison_summary(
            matched_agora_df, giq_df, segments_df, stats,
            args.agora_name, args.giq_name,
            summary_png
        )
        
        # Step 3: Create ORIGINAL 7 plots using matched AGORA data
        correlation1 = create_main_linear_dotplot(matched_agora_df, giq_df, args.agora_name, args.giq_name, output1)
        
        if correlation1 is not None:
            create_genome1_chromosome_distribution(matched_agora_df, giq_df, args.agora_name, args.giq_name, output2)
            create_genome2_chromosome_distribution(matched_agora_df, giq_df, args.agora_name, args.giq_name, output3)
            create_all_summaries_panel(matched_agora_df, giq_df, args.agora_name, args.giq_name, correlation1, output4)
            correlation5 = create_standard_gene_position_comparison(matched_agora_df, giq_df, args.agora_name, args.giq_name, output5)
            correlation6 = create_linearised_genome_comparison(matched_agora_df, giq_df, args.agora_name, args.giq_name, output6)
            create_side_by_side_chromosome_comparison(matched_agora_df, giq_df, args.agora_name, args.giq_name, output7)
        
        print("\n" + "="*60)
        print("COMPLETE ENHANCED DOTPLOT ANALYSIS FINISHED")
        print("="*60)
        print(f"Window size used: {args.window_size} genes")
        print(f"Confidence threshold: {args.threshold*100:.0f}%")
        print(f"Synteny agreement rate: {stats['synteny_agreement']:.1%}")
        
        print(f"\nFiles created:")
        print(f"  TSV OUTPUTS:")
        print(f"    1. Segmentation details: {segments_tsv}")
        print(f"    2. Matched AGORA TSV: {matched_agora_tsv}")
        print(f"  NEW ENHANCED PLOTS:")
        print(f"    3. Three-color dotplot: {dotplot_png}")
        print(f"    4. Method comparison: {summary_png}")
        print(f"  ORIGINAL 7 PLOTS:")
        print(f"    5. Main linear dotplot: {output1}")
        print(f"    6. AGORA chromosomes: {output2}")
        print(f"    7. GIQ chromosomes: {output3}")
        print(f"    8. All summaries: {output4}")
        print(f"    9. Standard comparison: {output5}")
        print(f"    10. Linearised comparison: {output6}")
        print(f"    11. Side-by-side chromosomes: {output7}")
        
        if correlation1: 
            print(f"\nCorrelation Results:")
            print(f"  Main linear correlation: {correlation1:.3f}")
            if correlation5: print(f"  Standard position correlation: {correlation5:.3f}")
            if correlation6: print(f"  Linearised correlation: {correlation6:.3f}")
        
        else:
            print("❌ No common genes found - comparison cannot be completed")
        
    except Exception as e:
        print(f"❌ Error: {e}")
        return 1
    
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())