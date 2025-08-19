"""
BUSCO processing module for the Genome Inversion Quantifier II
Handles parsing, filtering, and more from BUSCO results
"""

import pandas as pd
import logging
import os
from pathlib import Path
from scipy.stats import pearsonr

from config import (
    CONFIG
)

logger = logging.getLogger(__name__)


def parse_busco_table(busco_path, config):
    
    # Create data directory for TSV outputs
    output_dir = config.get('output_dir', '.')
    data_dir = Path(output_dir) / 'data'
    data_dir.mkdir(parents=True, exist_ok=True)
    
    # Determine species name for file naming
    if busco_path == config.get('first_busco_path'):
        species_name = config.get('first_species_name', 'species1')
    else:
        species_name = config.get('second_species_name', 'species2')
    
    with open(busco_path, 'r') as f:
        lines = [line for line in f if not line.startswith('#')]
    
    busco_data = []
    parsing_errors = []
    duplicate_entries = []
    
    seen_entries = set()
    
    for line_num, line in enumerate(lines, 1):
        if line.strip():
            parts = line.strip().split('\t')
            if len(parts) >= 6:
                try:
                    entry_key = (parts[0], parts[2], parts[3], parts[4])
                    
                    if config.get('enable_duplicate_handling', False):
                        if entry_key in seen_entries:
                            duplicate_entries.append(line_num)
                            continue
                        seen_entries.add(entry_key)
                    
                    start_str = parts[3]
                    end_str = parts[4]
                    
                    if start_str != 'N/A' and end_str != 'N/A':
                        # Take coordinates exactly as written - no inference
                        gene_start = int(start_str)  # Parts[3] as-is
                        gene_end = int(end_str)      # Parts[4] as-is
                        
                        strand = parts[5]
                        
                        entry = {
                            'busco_id': parts[0],
                            'status': parts[1],
                            'sequence': parts[2],
                            'gene_start': gene_start,  # Exact from file
                            'gene_end': gene_end,      # Exact from file  
                            'strand': strand,          # Original strand
                            'score': float(parts[6]) if len(parts) > 6 and parts[6] != 'N/A' else None,
                            'length': int(parts[7]) if len(parts) > 7 and parts[7] != 'N/A' else abs(gene_end - gene_start),
                            'line_number': line_num,
                            'original_start': gene_start,  # Same as gene_start now
                            'original_end': gene_end       # Same as gene_end now
                        }
                        
                        busco_data.append(entry)
                    else:
                        parsing_errors.append(f"Line {line_num}: N/A coordinates")
                        
                except (ValueError, IndexError) as e:
                    parsing_errors.append(f"Line {line_num}: {e}")
            else:
                parsing_errors.append(f"Line {line_num}: Too few columns ({len(parts)})")
    
    busco_df = pd.DataFrame(busco_data)
    
    # STAGE 1: Save raw parsed data (before corrections)
    raw_parsed_file = data_dir / f'stage_01_raw_parsed_busco_{species_name}.tsv'
    busco_df.to_csv(raw_parsed_file, sep='\t', index=False)
    logger.info(f"Saved raw parsed data: {raw_parsed_file}")
    
    # Now check and fix any rows where start > end
    needs_correction = busco_df['gene_start'] > busco_df['gene_end']
    
    # Step 1: Find which rows need correction
    needs_correction = busco_df['gene_start'] > busco_df['gene_end']  # Boolean mask per row

    # Step 2: Only correct those specific rows
    busco_df.loc[needs_correction, ['gene_start', 'gene_end']] = \
        busco_df.loc[needs_correction, ['gene_end', 'gene_start']].values  # Swap only bad rows

    # Step 3: Flip strand only for corrected rows
    busco_df.loc[needs_correction & (busco_df['strand'] == '+'), 'strand'] = '-'  # + → -
    busco_df.loc[needs_correction & (busco_df['strand'] == '-'), 'strand'] = '+'  # - → +
    
    # STAGE 2: Save corrected data (after coordinate and strand corrections)
    corrected_file = data_dir / f'stage_02_strand_corrected_busco_{species_name}.tsv'
    busco_df.to_csv(corrected_file, sep='\t', index=False)
    logger.info(f"Saved strand corrected data: {corrected_file}")
    
    logger.info(f"    Total data lines: {len(lines)}")
    logger.info(f"    Parsing errors: {len(parsing_errors)}")
    
    return busco_df


def filter_busco_genes(busco_df, config, quality_info=None, species_name=None):
    
    # Create data directory for TSV outputs
    output_dir = config.get('output_dir', '.')
    data_dir = Path(output_dir) / 'data'
    data_dir.mkdir(parents=True, exist_ok=True)
    
    # Determine species name if not provided
    if species_name is None:
        species_name = "unknown_species"
    
    initial_count = len(busco_df)
    filtering_stats = {'initial': initial_count}
        
    status_filter = config.get('busco_status_filter', ['Complete'])
    filtered_df = busco_df[busco_df['status'].isin(status_filter)].copy()
    filtering_stats['status_filter'] = len(filtered_df)
    
    # STAGE 3: Save filtered data (after status filtering)
    filtered_file = data_dir / f'stage_03_post_filter_busco_{species_name}.tsv'
    filtered_df.to_csv(filtered_file, sep='\t', index=False)
    logger.info(f"Saved post-filter data for {species_name}: {filtered_file}")
    
    final_count = len(filtered_df)
    exclusion_rate = (initial_count - final_count) / initial_count if initial_count > 0 else 0
    
    # Warning about excessive exclusions
    if config.get('enable_exclusion_warnings', False):
        if exclusion_rate > 0.4:
            logger.warning(f"Very high exclusion rate: {exclusion_rate:.1%} of BUSCOs excluded")
        elif exclusion_rate > 0.3:
            logger.warning(f"High exclusion rate: {exclusion_rate:.1%} of BUSCOs excluded")
    
    logger.info(f"  Filtering summary: {initial_count} -> {final_count} ({exclusion_rate:.1%} excluded)")

    
    return filtered_df


def detect_flips(df1, df2, config):
    
    # Create data directory for TSV outputs
    output_dir = config.get('output_dir', '.')
    data_dir = Path(output_dir) / 'data'
    data_dir.mkdir(parents=True, exist_ok=True)
    
    df1 = df1.copy()
    df2 = df2.copy()
    
    species1_name = config.get('first_species_name', 'Species1')
    species2_name = config.get('second_species_name', 'Species2')
    

    
    print(f"{len(df1)} genes in {species1_name}, and {len(df2)} in {species2_name}")
    print('-' * 80)
    
    common_buscos = set(df1['busco_id']) & set(df2['busco_id'])
    joined_genes = []
    
    for busco_id in common_buscos:
        gene1 = df1[df1['busco_id'] == busco_id].iloc[0]
        gene2 = df2[df2['busco_id'] == busco_id].iloc[0]
        
        joined_genes.append({
            'busco_id': busco_id,
            'chr1': gene1['sequence'],
            'chr2': gene2['sequence'],
            'start1': gene1['gene_start'],
            'end1': gene1['gene_end'],
            'start2': gene2['gene_start'],
            'end2': gene2['gene_end'],
            'strand1': gene1['strand'],
            'strand2': gene2['strand'],
            'is_flipped': gene1['strand'] != gene2['strand']
        })
    
    joined_genes.sort(key=lambda x: x['busco_id'])
    joined_df = pd.DataFrame(joined_genes)
    
    # STAGE 4: Save joined genes dataframe
    joined_genes_file = data_dir / 'stage_04_joined_genes.tsv'
    joined_df.to_csv(joined_genes_file, sep='\t', index=False)
    logger.info(f"Saved joined genes data: {joined_genes_file}")
    
    total_flipped = joined_df['is_flipped'].sum()
    print(f"{len(joined_genes)} common genes with {total_flipped} flipped whole genes")
    
    chromosome_dict = {}
    for gene in joined_genes:
        chr_pair = (gene['chr1'], gene['chr2'])
        if chr_pair not in chromosome_dict:
            chromosome_dict[chr_pair] = []
        chromosome_dict[chr_pair].append(gene)
    
    print(f"{len(chromosome_dict)} chromosome pairs across {species1_name} and {species2_name}")
    
    inversion_results = []
    
    for chr_pair, genes in chromosome_dict.items():
        chr1, chr2 = chr_pair
        
        genes.sort(key=lambda x: x['start1'])
        
        total_genes = len(genes)
        flipped_genes = sum(1 for g in genes if g['is_flipped'])
        flip_rate = flipped_genes / total_genes if total_genes > 0 else 0
        
        # print(f"   {chr1} vs {chr2}: {total_genes} genes, {flipped_genes} flipped (rate: {flip_rate:.2f})")
        
        correlation = None
        p_value = None
        strand_consistency = (total_genes - flipped_genes) / total_genes if total_genes > 0 else 1.0
        
        # print(f"   Strand consistency is {strand_consistency:.2f} across the {species1_name} and {species2_name} genome")
        
        if total_genes >= 3:
            pos1_list = [g['start1'] for g in genes]
            pos2_list = [g['start2'] for g in genes]
            
            try:
                correlation, p_value = pearsonr(pos1_list, pos2_list)
                # print(f"   Correlation: {correlation:.3f} (p={p_value:.3f})")
            except:
                correlation, p_value = 0.0, 1.0
        
        inversion_type = "None"
        if flipped_genes > 0:
            if total_genes >= 3 and correlation is not None and correlation < 0:
                inversion_type = "Large Multigene Inversion within"
            elif flip_rate > 0.5:
                inversion_type = "Multiple Single Gene Flips"
            else:
                inversion_type = "Single Gene Inversions"
        
        inversion_results.append({
            'chr1': chr1,
            'chr2': chr2,
            'total_genes': total_genes,
            'flipped_genes': flipped_genes,
            'flip_rate': flip_rate,
            'correlation': correlation,
            'p_value': p_value,
            'strand_consistency': strand_consistency,
            'inversion_type': inversion_type,
            'gene_list': [g['busco_id'] for g in genes]
        })
    
    results_df = pd.DataFrame(inversion_results)
    
    # STAGE 5: Save inversion results dataframe
    inversion_results_file = data_dir / 'stage_05_inversion_results.tsv'
    results_df.to_csv(inversion_results_file, sep='\t', index=False)
    logger.info(f"Saved inversion results: {inversion_results_file}")
    
    print('*' * 80)
    print('*' * 80)
    print(f"   - Total chromosome pairs: {len(results_df)}")
    print(f"   - Pairs with inversions: {len(results_df[results_df['flipped_genes'] > 0])}")
    print(f"   - Total flipped genes: {results_df['flipped_genes'].sum()}")
    
    for i, gene in enumerate(joined_genes[:5]):
        print(f"  {gene['busco_id']}: strand1={gene['strand1']}, strand2={gene['strand2']}, flipped={gene['is_flipped']}")
    
    
    # Add this debug in detect_flips
    strand1_counts = {}
    strand2_counts = {}
    for gene in joined_genes:
        strand1_counts[gene['strand1']] = strand1_counts.get(gene['strand1'], 0) + 1
        strand2_counts[gene['strand2']] = strand2_counts.get(gene['strand2'], 0) + 1

    print(f"DEBUG - Strand distribution:")
    print(f"  Genome1 strands: {strand1_counts}")
    print(f"  Genome2 strands: {strand2_counts}")
    
    
    return results_df, joined_df