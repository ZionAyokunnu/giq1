#!/usr/bin/env python3
"""
Enhanced BUSCO dotplotter with chromosome content matching and proper coordinate system.
Designed for AGORA (single linear chromosome) vs GIQ (multiple consolidated chromosomes) comparison.

Usage:
python3 busco_dotplotter.py \
    /Users/zionayokunnu/Documents/Giq/compare/Dioctria_linearis_linearized.tsv \
    /Users/zionayokunnu/Documents/Bibionidae/busco-tables/Dioctria_rufipes.tsv \
    compare/experiment/segmented/Lin-lin_vs_Ruf_w20t10 \
    --agora-name "Dioctria linearis Linearised" \
    --giq-name "Dioctria rufipes Normal" \
    --agora-coord-mode genomic \
    --giq-coord-mode genomic \
    --window-size 20 --threshold 0.1
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse
from pathlib import Path
import seaborn as sns
from collections import Counter, defaultdict


def load_busco_tsv(file_path: str, name: str):
    """Load and process BUSCO TSV file - handles both clean TSV and BUSCO format with # comments"""
    print(f"Loading {name} from: {file_path}")
    
    try:
        # First, check if file has BUSCO comment headers
        with open(file_path, 'r') as f:
            first_line = f.readline().strip()
        
        if first_line.startswith('#'):
            # BUSCO format with comments - skip comment lines
            df = pd.read_csv(file_path, sep='\t', comment='#', header=None, 
                 names=['busco_id', 'status', 'sequence', 'gene_start', 'gene_end', 'strand', 'score', 'length'])
            
            df['gene_start'] = pd.to_numeric(df['gene_start'], errors='coerce')
            df['gene_end'] = pd.to_numeric(df['gene_end'], errors='coerce')
            
            df = df[df['status'] == 'Complete'].copy()
            print(f"  After filtering to Complete: {len(df)} genes")
            
        else:
            # Regular TSV format
            df = pd.read_csv(file_path, sep='\t')
            print(f"  Loaded clean TSV format")
        
            df['gene_start'] = pd.to_numeric(df['gene_start'], errors='coerce')
            df['gene_end'] = pd.to_numeric(df['gene_end'], errors='coerce')
            
        print(f"  Loaded {len(df)} genes")
        
        # Handle different column naming conventions
        column_mapping = {
            'Busco id': 'busco_id',
            'Status': 'status', 
            'Sequence': 'sequence',
            'Gene Start': 'gene_start',
            'Gene End': 'gene_end',
            'Strand': 'strand'
        }
        
        # Rename columns if they exist
        df = df.rename(columns=column_mapping)
        
        # Ensure required columns exist
        required_cols = ['busco_id', 'sequence', 'gene_start', 'gene_end']
        missing_cols = [col for col in required_cols if col not in df.columns]
        
        if missing_cols:
            print(f"  Available columns: {list(df.columns)}")
            raise ValueError(f"Missing required columns: {missing_cols}")
        
        # Sort by position within each chromosome
        df = df.sort_values(['sequence', 'gene_start'])
        
        # Show chromosome distribution
        for chrom in sorted(df['sequence'].unique()):
            chrom_genes = df[df['sequence'] == chrom]
            print(f"  {chrom}: {len(chrom_genes)} genes")
        
        return df
        
    except Exception as e:
        print(f"Error loading {name}: {e}")
        raise


def linearize_genomic_coordinates(df):
    """
    Linearize genomic coordinates using cumulative chromosome lengths (like simple dotplotter)
    """
    linearized_df = df.copy()
    
    # Group by chromosome and find real chromosome sizes
    chr_offsets = {}
    cumulative_offset = 0
    
    for chrom in sorted(df['sequence'].unique()):
        chr_offsets[chrom] = cumulative_offset
        chrom_genes = df[df['sequence'] == chrom]
        
        # Use REAL chromosome size (actual max position)
        chr_size = chrom_genes['gene_end'].max()
        cumulative_offset += chr_size
    
    # Add linearized positions
    linearized_df['linearized_position'] = linearized_df.apply(
        lambda row: chr_offsets[row['sequence']] + row['gene_start'], axis=1
    )
    
    return linearized_df, chr_offsets


def get_genome_positions(df, coord_mode, name):
    """
    Get genome positions based on coordinate mode
    
    Args:
        df: BUSCO DataFrame
        coord_mode: 'artificial' or 'genomic'
        name: genome name for logging
    
    Returns:
        positions_dict: dict mapping busco_id to position
        chr_boundaries: dict with chromosome boundary info (None for artificial)
    """
    print(f"  Processing {name} coordinates in '{coord_mode}' mode")
    
    if coord_mode == 'artificial':
        # Use gene_start directly - don't linearize even if multiple chromosomes
        positions_dict = dict(zip(df['busco_id'], df['gene_start']))
        chr_boundaries = None
        print(f"    Using artificial coordinates directly ({len(positions_dict)} genes)")
        
    elif coord_mode == 'genomic':
        if len(df['sequence'].unique()) == 1:
            # Single chromosome - use direct positions
            positions_dict = dict(zip(df['busco_id'], df['gene_start']))
            chr_boundaries = None
            print(f"    Single chromosome - using direct positions ({len(positions_dict)} genes)")
        else:
            # Multiple chromosomes - linearize using cumulative lengths
            linearized_df, chr_offsets = linearize_genomic_coordinates(df)
            positions_dict = dict(zip(linearized_df['busco_id'], linearized_df['linearized_position']))
            
            # Calculate chromosome boundaries for plotting
            chr_boundaries = {}
            for chrom in sorted(df['sequence'].unique()):
                chrom_genes = df[df['sequence'] == chrom]
                chr_size = chrom_genes['gene_end'].max()
                chr_boundaries[chrom] = {
                    'start': chr_offsets[chrom],
                    'end': chr_offsets[chrom] + chr_size,
                    'midpoint': chr_offsets[chrom] + (chr_size / 2)
                }
            print(f"    Linearized {len(df['sequence'].unique())} chromosomes ({len(positions_dict)} genes)")
    
    else:
        raise ValueError(f"Invalid coord_mode: {coord_mode}. Must be 'artificial' or 'genomic'")
    
    return positions_dict, chr_boundaries


def regenerate_agora_coordinates(matched_agora_df):
    """
    Regenerate artificial coordinates for AGORA after segmentation.
    Each chromosome gets new continuous coordinates starting from 0.
    Preserves original gene lengths (gene_end - gene_start).
    """
    print(f"\nRegenerating AGORA coordinates after segmentation:")
    
    # Store original gene lengths
    matched_agora_df['original_gene_length'] = matched_agora_df['gene_end'] - matched_agora_df['gene_start']
    
    # Group by new chromosome assignments and regenerate coordinates
    regenerated_df = matched_agora_df.copy()
    
    for chrom in sorted(matched_agora_df['sequence'].unique()):
        chrom_mask = matched_agora_df['sequence'] == chrom
        chrom_genes = matched_agora_df[chrom_mask].copy()
        
        print(f"  {chrom}: {len(chrom_genes)} genes")
        
        # Generate new continuous coordinates for this chromosome
        new_positions = []
        current_pos = 0
        
        for _, gene in chrom_genes.iterrows():
            gene_length = gene['original_gene_length']
            
            # Assign new positions
            new_positions.append({
                'gene_start': current_pos,
                'gene_end': current_pos + gene_length
            })
            
            # Move to next position
            current_pos += gene_length
        
        # Update the dataframe with new coordinates
        for i, (idx, _) in enumerate(chrom_genes.iterrows()):
            regenerated_df.loc[idx, 'gene_start'] = new_positions[i]['gene_start']
            regenerated_df.loc[idx, 'gene_end'] = new_positions[i]['gene_end']
        
        print(f"    New coordinates: 0 to {current_pos - gene_length}")
    
    print(f"  Coordinate regeneration completed")
    return regenerated_df


def segment_agora_by_giq_dominance(agora_df, giq_df, window_size=50, threshold=0.6, 
                                  output_segments_tsv=None, output_matched_tsv=None):
    """
    Segment AGORA's linear gene order based on GIQ chromosome dominance.
    After segmentation, regenerates proper artificial coordinates per chromosome.
    
    Args:
        agora_df: AGORA BUSCO DataFrame (single chromosome)
        giq_df: GIQ BUSCO DataFrame (multiple chromosomes) 
        window_size: Number of consecutive AGORA genes to analyze
        threshold: Minimum percentage for confident assignment
        output_segments_tsv: Optional path to save segment analysis
        output_matched_tsv: Optional path to save final matched AGORA TSV
    
    Returns:
        matched_agora_df: AGORA DataFrame with updated chromosome assignments and coordinates
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
    
    # Get AGORA gene order (sorted by gene_start position)
    agora_ordered = agora_df.sort_values('gene_start').copy()
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
    
    # CRITICAL: Regenerate coordinates after segmentation
    matched_agora_df = regenerate_agora_coordinates(matched_agora_df)
    
    if output_matched_tsv is not None:
        debug_path = output_matched_tsv.replace('_agora_matched.tsv', '_debug_segmented_coords.tsv')
        matched_agora_df.to_csv(debug_path, sep='\t', index=False)
        print(f"DEBUG: Segmented coordinates saved to: {debug_path}")

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


def create_three_color_chromosome_dotplot(agora_df, giq_df, agora_name, giq_name, 
                                        agora_coord_mode, giq_coord_mode, output_path):
    """
    Create three-color dotplot showing:
    - Blue: GIQ-only genes
    - Red: AGORA-only genes  
    - Green: Collision genes (both methods agree on same chromosome)
    
    Now uses proper coordinate systems instead of linear_position
    """
    
    print(f"\nCreating three-color chromosome dotplot...")
    
    # Get positions using proper coordinate modes
    agora_positions, agora_boundaries = get_genome_positions(agora_df, agora_coord_mode, agora_name)
    giq_positions, giq_boundaries = get_genome_positions(giq_df, giq_coord_mode, giq_name)
    
    # Find all genes in both datasets
    agora_genes = set(agora_df['busco_id'])
    giq_genes = set(giq_df['busco_id'])
    
    # Create chromosome mappings
    agora_chr_map = dict(zip(agora_df['busco_id'], agora_df['sequence']))
    giq_chr_map = dict(zip(giq_df['busco_id'], giq_df['sequence']))
    
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
        if gene in agora_positions and gene in giq_positions:
            plot_data['agora_pos'].append(agora_positions[gene])
            plot_data['giq_pos'].append(giq_positions[gene])
            plot_data['color'].append('green')
            plot_data['gene_type'].append('collision')
    
    # Add mismatch genes (purple - show disagreement)
    for gene in mismatch_genes:
        if gene in agora_positions and gene in giq_positions:
            plot_data['agora_pos'].append(agora_positions[gene])
            plot_data['giq_pos'].append(giq_positions[gene])
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
    if agora_boundaries:
        for chrom, info in agora_boundaries.items():
            if info['start'] > 0:
                plt.axvline(x=info['start'], color='red', linewidth=1, alpha=0.5, linestyle='--')
            # Add chromosome labels
            plt.text(info['midpoint'], plt.ylim()[1] * 0.95, chrom.replace('_agora', ''), 
                    rotation=45, ha='center', va='top', fontsize=8, color='red', alpha=0.8)
    
    # Add chromosome boundaries and labels for GIQ
    if giq_boundaries:
        for chrom, info in giq_boundaries.items():
            if info['start'] > 0:
                plt.axhline(y=info['start'], color='blue', linewidth=1, alpha=0.5, linestyle='--')
            # Add chromosome labels
            plt.text(plt.xlim()[0] * 0.02, info['midpoint'], chrom, 
                    rotation=0, ha='left', va='center', fontsize=8, color='blue', alpha=0.8)
    
    # Formatting
    plt.xlabel(f'{agora_name} Position ({agora_coord_mode})', fontsize=12)
    plt.ylabel(f'{giq_name} Position ({giq_coord_mode})', fontsize=12)
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


def create_side_by_side_chromosome_comparison(agora_df, giq_df, agora_name, giq_name, output_path):
    """Side-by-side Chromosome Comparison - Both genomes on same chart with flexible matching"""
    
    # Find common genes
    common_genes = set(agora_df['busco_id']) & set(giq_df['busco_id'])
    
    if len(common_genes) == 0:
        return None
    
    # Filter to common genes only
    g1_common = agora_df[agora_df['busco_id'].isin(common_genes)].copy()
    g2_common = giq_df[giq_df['busco_id'].isin(common_genes)].copy()
    
    # FLEXIBLE MATCHING: Extract base chromosome names
    def get_base_chromosome_name(chr_name):
        """Remove suffixes like _agora, _giq to get base name"""
        return chr_name.replace('_agora', '').replace('_giq', '').replace('_matched', '')
    
    # Get all unique base chromosome names
    base_chromosomes = set()
    
    # Extract from both genomes
    for chr_name in g1_common['sequence'].unique():
        base_chromosomes.add(get_base_chromosome_name(chr_name))
    
    for chr_name in g2_common['sequence'].unique():
        base_chromosomes.add(get_base_chromosome_name(chr_name))
    
    # Count genes per base chromosome
    comparison_data = []
    for base_chr in sorted(base_chromosomes):
        # Count in genome1 (look for exact match or with suffix)
        g1_count = 0
        for chr_variant in [base_chr, f"{base_chr}_agora", f"{base_chr}_matched"]:
            if chr_variant in g1_common['sequence'].values:
                g1_count += len(g1_common[g1_common['sequence'] == chr_variant])
        
        # Count in genome2 (look for exact match or with suffix)  
        g2_count = 0
        for chr_variant in [base_chr, f"{base_chr}_giq", f"{base_chr}_ordinal"]:
            if chr_variant in g2_common['sequence'].values:
                g2_count += len(g2_common[g2_common['sequence'] == chr_variant])
        
        if g1_count > 0 or g2_count > 0:  # Only include if data exists
            comparison_data.append({
                'chromosome': base_chr,
                'genome1_count': g1_count,
                'genome2_count': g2_count
            })
    
    if not comparison_data:
        return None
    
    # Create side-by-side bar plot
    fig, ax = plt.subplots(figsize=(14, 8))
    
    chromosomes = [item['chromosome'] for item in comparison_data]
    g1_counts = [item['genome1_count'] for item in comparison_data]
    g2_counts = [item['genome2_count'] for item in comparison_data]
    
    x = np.arange(len(chromosomes))
    width = 0.35
    
    ax.bar(x - width/2, g1_counts, width, label=agora_name, alpha=0.7, color='lightblue', edgecolor='darkblue')
    ax.bar(x + width/2, g2_counts, width, label=giq_name, alpha=0.7, color='lightgreen', edgecolor='darkgreen')
    
    ax.set_xlabel('Chromosome', fontsize=12)
    ax.set_ylabel('Common Gene Count', fontsize=12)
    ax.set_title(f'Side-by-Side Chromosome Comparison\n{agora_name} vs {giq_name}', fontsize=14)
    ax.set_xticks(x)
    ax.set_xticklabels(chromosomes, rotation=45, ha='right')
    ax.legend()
    ax.grid(True, alpha=0.3, axis='y')
    
    # Add count labels on bars
    for i, (g1, g2) in enumerate(zip(g1_counts, g2_counts)):
        if g1 > 0:
            ax.text(i - width/2, g1 + max(g1_counts + g2_counts) * 0.01, str(g1), 
                   ha='center', va='bottom', fontsize=10)
        if g2 > 0:
            ax.text(i + width/2, g2 + max(g1_counts + g2_counts) * 0.01, str(g2), 
                   ha='center', va='bottom', fontsize=10)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    
    print(f"Side-by-side chromosome comparison saved to: {output_path}")
    print(f"  Matched chromosomes: {', '.join(chromosomes)}")
    print(f"  {agora_name} genes: {sum(g1_counts)}")
    print(f"  {giq_name} genes: {sum(g2_counts)}")
    
    return True


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
Window size: {segments_df.iloc[0]['window_size'] if len(segments_df) > 0 else 'N/A'} genes
Total segments: {len(segments_df)}

Segment confidence:
"""
    
    # Add segment confidence breakdown
    if len(segments_df) > 0:
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


def main():
    parser = argparse.ArgumentParser(description='Enhanced BUSCO dotplotter with chromosome matching and proper coordinates')
    parser.add_argument('agora_tsv', help='AGORA BUSCO TSV file (single linear chromosome)')
    parser.add_argument('giq_tsv', help='GIQ BUSCO TSV file (multiple consolidated chromosomes)')
    parser.add_argument('output_prefix', help='Output file prefix (will create multiple files)')
    parser.add_argument('--agora-name', default='AGORA', help='Name for AGORA method')
    parser.add_argument('--giq-name', default='GIQ', help='Name for GIQ method')
    parser.add_argument('--agora-coord-mode', choices=['artificial', 'genomic'], default='artificial', 
                        help='Coordinate mode for AGORA: artificial (use gene_start directly) or genomic (linearize if multi-chr)')
    parser.add_argument('--giq-coord-mode', choices=['artificial', 'genomic'], default='genomic',
                        help='Coordinate mode for GIQ: artificial (use gene_start directly) or genomic (linearize if multi-chr)')
    parser.add_argument('--window-size', type=int, default=50, help='Window size for chromosome segmentation (default: 50)')
    parser.add_argument('--threshold', type=float, default=0.6, help='Confidence threshold for chromosome assignment (default: 0.6)')
    
    args = parser.parse_args()
    
    try:
        print("="*60)
        print("ENHANCED BUSCO DOTPLOTTER WITH PROPER COORDINATES")
        print("="*60)
        
        # Load datasets
        agora_df = load_busco_tsv(args.agora_tsv, args.agora_name)
        giq_df = load_busco_tsv(args.giq_tsv, args.giq_name)
        
        # Create output paths for the 3 essential visualizations
        base_path = Path(args.output_prefix)
        
        # TSV outputs
        segments_tsv = f"{base_path}_segmentation_details.tsv"
        matched_agora_tsv = f"{base_path}_agora_matched.tsv"
        
        # The 3 essential plots
        dotplot_png = f"{base_path}_three_color_dotplot.png"
        side_by_side_png = f"{base_path}_side_by_side_chromosomes.png"
        summary_png = f"{base_path}_method_comparison_summary.png"
        
        # Step 1: Segment AGORA by GIQ chromosome dominance
        matched_agora_df, segments_df = segment_agora_by_giq_dominance(
            agora_df, giq_df, 
            window_size=args.window_size, 
            threshold=args.threshold,
            output_segments_tsv=segments_tsv,
            output_matched_tsv=matched_agora_tsv
        )
        
        print(f"\nCreating 3 essential visualization plots...")
        
        # Step 2: Create the 3 essential plots with proper coordinates
        stats = create_three_color_chromosome_dotplot(
            matched_agora_df, giq_df, 
            args.agora_name, args.giq_name,
            args.agora_coord_mode, args.giq_coord_mode,
            dotplot_png
        )
        
        create_side_by_side_chromosome_comparison(
            matched_agora_df, giq_df, 
            args.agora_name, args.giq_name,
            side_by_side_png
        )
        
        create_method_comparison_summary(
            matched_agora_df, giq_df, segments_df, stats,
            args.agora_name, args.giq_name,
            summary_png
        )
        
        print("\n" + "="*60)
        print("ENHANCED DOTPLOT ANALYSIS COMPLETED")
        print("="*60)
        print(f"Window size used: {args.window_size} genes")
        print(f"Confidence threshold: {args.threshold*100:.0f}%")
        print(f"AGORA coordinate mode: {args.agora_coord_mode}")
        print(f"GIQ coordinate mode: {args.giq_coord_mode}")
        print(f"Synteny agreement rate: {stats['synteny_agreement']:.1%}")
        
        print(f"\nFiles created:")
        print(f"  TSV OUTPUTS:")
        print(f"    1. Segmentation details: {segments_tsv}")
        print(f"    2. Matched AGORA TSV: {matched_agora_tsv}")
        print(f"  ESSENTIAL PLOTS:")
        print(f"    3. Three-color dotplot: {dotplot_png}")
        print(f"    4. Side-by-side chromosomes: {side_by_side_png}")
        print(f"    5. Method comparison summary: {summary_png}")
        
        print(f"\nResults Summary:")
        print(f"  Common genes: {stats['common_genes']:,}")
        print(f"  Collision genes (same chr): {stats['collision_genes']:,}")
        print(f"  Mismatch genes (diff chr): {stats['mismatch_genes']:,}")
        print(f"  Synteny agreement: {stats['synteny_agreement']:.1%}")
        
        if stats['synteny_agreement'] > 0.8:
            print("🟢 High synteny conservation")
        elif stats['synteny_agreement'] > 0.6:
            print("🟡 Moderate synteny conservation")
        else:
            print("🔴 High chromosomal rearrangement/uncertainty")
        
    except Exception as e:
        print(f"❌ Error: {e}")
        return 1

    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())