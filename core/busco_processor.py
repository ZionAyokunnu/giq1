

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
    output_dir = config.get('output_dir', '.')
    data_dir = Path(output_dir) / 'data'
    data_dir.mkdir(parents=True, exist_ok=True)
    
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
                        gene_start = int(start_str)  
                        gene_end = int(end_str)     
                        
                        strand = parts[5]
                        
                        entry = {
                            'busco_id': parts[0],
                            'status': parts[1],
                            'sequence': parts[2],
                            'gene_start': gene_start,  
                            'gene_end' : gene_end,
                            'strand': strand,          
                            'score': float(parts[6]),
                            'length': int(parts[7]) if len(parts) > 7 and parts[7] != 'N/A' else abs(gene_end - gene_start),
                            'line_number': line_num,
      
                        }
                        
                        busco_data.append(entry)
                    else:
                        parsing_errors.append(f"Line {line_num}: N/A coordinates")
                        
                except (ValueError, IndexError) as e:
                    parsing_errors.append(f"Line {line_num}: {e}")
            else:
                parsing_errors.append(f"Line {line_num}: Too few columns ({len(parts)})")
    
    busco_df = pd.DataFrame(busco_data)
    



    
    return busco_df


def correct_chromosome_orientation(df1, df2,joined_df, config):

    joined_genes = joined_df.to_dict('records')
    chromosome_dict = {}
    for gene in joined_genes:
        chr_pair = (gene['chr1'], gene['chr2'])
        if chr_pair not in chromosome_dict:
            chromosome_dict[chr_pair] = []
        chromosome_dict[chr_pair].append(gene)
    
    corrected_dfs = [df1.copy(), df2.copy()]
    
    for chr_pair, genes in chromosome_dict.items():
        chr1, chr2 = chr_pair
        
        if len(genes) >= 3: 
            genes.sort(key=lambda x: x['start1'])
            
            pos1_list = [g['start1'] for g in genes]
            pos2_list = [g['start2'] for g in genes]
            
            correlation, p_value = pearsonr(pos1_list, pos2_list)
            
            if correlation < -0.5 and p_value < 0.05:
                chr2_genes = corrected_dfs[1][corrected_dfs[1]['sequence'] == chr2]
                
   
                max_pos = chr2_genes['gene_start'].max()
                corrected_dfs[1].loc[corrected_dfs[1]['sequence'] == chr2, 'gene_start'] = \
                    max_pos - chr2_genes['gene_start']
                corrected_dfs[1].loc[corrected_dfs[1]['sequence'] == chr2, 'strand'] = \
                    chr2_genes['strand'].map({'+': '-', '-': '+'})
    
    return corrected_dfs[0], corrected_dfs[1]

def filter_busco_genes(busco_df, config, quality_info=None, species_name=None):
    

    output_dir = config.get('output_dir', '.')
    data_dir = Path(output_dir) / 'data'
    data_dir.mkdir(parents=True, exist_ok=True)
    

    if species_name is None:
        species_name = "unknown_species"
    
    initial_count = len(busco_df)
    filtering_stats = {'initial': initial_count}
        
    status_filter = config.get('busco_status_filter', ['Complete'])
    filtered_df = busco_df[busco_df['status'].isin(status_filter)].copy()
    filtering_stats['status_filter'] = len(filtered_df)
    

    
    final_count = len(filtered_df)
    exclusion_rate = (initial_count - final_count) / initial_count if initial_count > 0 else 0
    
    
    return filtered_df


def detect_flips(df1, df2, config):
   df1 = df1.copy()
   df2 = df2.copy()
   
   species1_name = config.get('first_species_name', 'Species1')
   species2_name = config.get('second_species_name', 'Species2')
   
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
   
   chromosome_dict = {}
   for gene in joined_genes:
       chr_pair = (gene['chr1'], gene['chr2'])
       if chr_pair not in chromosome_dict:
           chromosome_dict[chr_pair] = []
       chromosome_dict[chr_pair].append(gene)
   
   inversion_results = []
   all_inversion_events = []
   
   for chr_pair, genes in chromosome_dict.items():
       chr1, chr2 = chr_pair
       
       genes.sort(key=lambda x: x['start1'])
       
       total_genes = len(genes)
       flipped_genes = sum(1 for g in genes if g['is_flipped'])
       flip_rate = flipped_genes / total_genes if total_genes > 0 else 0
       
       correlation = None
       p_value = None
       strand_consistency = (total_genes - flipped_genes) / total_genes if total_genes > 0 else 1.0
       
       if total_genes >= 3:
           pos1_list = [g['start1'] for g in genes]
           pos2_list = [g['start2'] for g in genes]
           
           try:
               correlation, p_value = pearsonr(pos1_list, pos2_list)
           except:
               correlation, p_value = 0.0, 1.0
       
       inversion_events = []
       current_block = []
       
       for i, gene in enumerate(genes):
           if gene['is_flipped']:
               current_block.append({
                   'busco_id': gene['busco_id'],
                   'position': i,
                   'start1': gene['start1'],
                   'start2': gene['start2']
               })
           else:
               if current_block:
                   inversion_events.append({
                       'chr1': chr1,
                       'chr2': chr2,
                       'type': 'multi-gene' if len(current_block) > 1 else 'single-gene',
                       'gene_count': len(current_block),
                       'start_gene': current_block[0]['busco_id'],
                       'end_gene': current_block[-1]['busco_id'],
                       'start_pos1': current_block[0]['start1'],
                       'end_pos1': current_block[-1]['start1'],
                       'span_bp': abs(current_block[-1]['start1'] - current_block[0]['start1']),
                       'genes': [g['busco_id'] for g in current_block]
                   })
                   current_block = []
       
       if current_block:
           inversion_events.append({
               'chr1': chr1,
               'chr2': chr2,
               'type': 'multi-gene' if len(current_block) > 1 else 'single-gene',
               'gene_count': len(current_block),
               'start_gene': current_block[0]['busco_id'],
               'end_gene': current_block[-1]['busco_id'],
               'start_pos1': current_block[0]['start1'],
               'end_pos1': current_block[-1]['start1'],
               'span_bp': abs(current_block[-1]['start1'] - current_block[0]['start1']),
               'genes': [g['busco_id'] for g in current_block]
           })
       
       all_inversion_events.extend(inversion_events)
       
       single_gene_count = sum(1 for e in inversion_events if e['type'] == 'single-gene')
       multi_gene_count = sum(1 for e in inversion_events if e['type'] == 'multi-gene')
       largest_inversion = max([e['gene_count'] for e in inversion_events]) if inversion_events else 0
       
       inversion_type = "Syntenic (same orientation)"
       if flipped_genes > 0:
           if total_genes >= 3 and correlation is not None and correlation < 0:
               inversion_type = "Large-scale inversion (negative correlation)"
           elif flip_rate > 0.5:
               inversion_type = "High strand discordance (Many small inversions[>50%])"
           else:
               inversion_type = "Minor strand differences (<50%)"
       
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
           'single_gene_inversions': single_gene_count,
           'multi_gene_inversions': multi_gene_count,
           'largest_inversion': largest_inversion,
           'total_inversion_events': len(inversion_events),
           'gene_list': [g['busco_id'] for g in genes]
       })
   
   results_df = pd.DataFrame(inversion_results)
   
   config['all_inversion_events'] = all_inversion_events
   
   return results_df, joined_df


