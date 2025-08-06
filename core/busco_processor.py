"""
BUSCO processing module for the Genome Inversion analyser
Handles parsing, filtering, and sequence extraction from BUSCO results
"""

import pandas as pd
import logging
from scipy.stats import pearsonr

from config import (
    CONFIG
)

logger = logging.getLogger(__name__)


def parse_busco_table(busco_path, config):
    """BUSCO table parsing. It handles negative strand genes"""
    logger.info(f"Parsing BUSCO table from {busco_path}")
    
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
                        coord1 = int(start_str)
                        coord2 = int(end_str)
                    
                        gene_start = min(coord1, coord2)
                        gene_end = max(coord1, coord2)
                        
                        strand = parts[5] if len(parts) > 5 else '+'
                        if coord1 > coord2 and strand == '+':
                            strand = '-' 
                        elif coord1 < coord2 and strand == '-':
                            strand = '+'
                        
                        entry = {
                            'busco_id': parts[0],
                            'status': parts[1],
                            'sequence': parts[2],
                            'gene_start': gene_start,
                            'gene_end': gene_end,
                            'strand': strand,
                            'score': float(parts[6]) if len(parts) > 6 and parts[6] != 'N/A' else None,
                            'length': int(parts[7]) if len(parts) > 7 and parts[7] != 'N/A' else (gene_end - gene_start),
                            'line_number': line_num,
                            'original_start': coord1,
                            'original_end': coord2
                        }
                        
                        busco_data.append(entry)
                    else:
                        parsing_errors.append(f"Line {line_num}: N/A coordinates")
                        
                except (ValueError, IndexError) as e:
                    parsing_errors.append(f"Line {line_num}: {e}")
            else:
                parsing_errors.append(f"Line {line_num}: Too few columns ({len(parts)})")
    
    busco_df = pd.DataFrame(busco_data)
    
    total_lines = len(lines)
    success_rate = len(busco_data) / total_lines * 100 if total_lines > 0 else 0
    
    logger.info(f"    Total data lines: {total_lines}")
    logger.info(f"    Parsing errors: {len(parsing_errors)}")
    
    return busco_df


def filter_busco_genes(busco_df, config, quality_info=None):
    """BUSCO filtering"""
    
    initial_count = len(busco_df)
    filtering_stats = {'initial': initial_count}
        
    status_filter = config.get('busco_status_filter', ['Complete'])
    filtered_df = busco_df[busco_df['status'].isin(status_filter)].copy()
    filtering_stats['status_filter'] = len(filtered_df)
    
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

def detect_inversions(df1, df2, config):
    """Detect inversions"""
    
    df1 = df1.copy()
    df2 = df2.copy()
    
    species1_name = config.get('first_species_name', 'Species1')
    species2_name = config.get('second_species_name', 'Species2')
    
    print(f"{len(df1)} genes in {species1_name}, and {len(df2)} in {species2_name}")
    
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
        
        print(f"   {chr1} vs {chr2}: {total_genes} genes, {flipped_genes} flipped (rate: {flip_rate:.2f})")
        
        correlation = None
        p_value = None
        strand_consistency = (total_genes - flipped_genes) / total_genes if total_genes > 0 else 1.0
        
        print(f"   Strand consistency is {strand_consistency:.2f} across the {species1_name} and {species2_name} genome")
        
        if total_genes >= 3:
            pos1_list = [g['start1'] for g in genes]
            pos2_list = [g['start2'] for g in genes]
            
            try:
                correlation, p_value = pearsonr(pos1_list, pos2_list)
                print(f"   Correlation: {correlation:.3f} (p={p_value:.3f})")
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
    
    print('*' * 80)
    print("ANALYSIS COMPLETE")
    print('*' * 80)
    print(f"   - Total chromosome pairs: {len(results_df)}")
    print(f"   - Pairs with inversions: {len(results_df[results_df['flipped_genes'] > 0])}")
    print(f"   - Total flipped genes: {results_df['flipped_genes'].sum()}")
    
    return results_df, joined_df


if __name__ == "__main__":
    file1 = CONFIG['first_busco_path']
    file2 = CONFIG['second_busco_path']
    species1_name = CONFIG['first_species_name']
    species2_name = CONFIG['second_species_name']
    
    df1 = parse_busco_table(file1, CONFIG)
    df2 = parse_busco_table(file2, CONFIG)
    

    df1 = df1[df1['status'].isin(CONFIG.get('busco_status_filter', ['Complete']))]
    df2 = df2[df2['status'].isin(CONFIG.get('busco_status_filter', ['Complete']))]
    
    inversion_results, joined_data = detect_inversions(df1, df2, CONFIG)
    
    print("\n" + "="*80)
    print("INVERSION ANALYSIS RESULTS")
    print("="*80)
    print(inversion_results[['chr1', 'chr2', 'total_genes', 'flipped_genes', 
                           'correlation', 'strand_consistency', 'inversion_type']])
    
    inversion_results.to_csv('inversion_analysis.csv', index=False)
    joined_data.to_csv('joined_gene_data_for_dotplot.csv', index=False)

else:
        print("No inversions found - end")
    