import sys
import os
import subprocess
import logging
from pathlib import Path
import pandas as pd
import numpy as np
from typing import Dict, Optional, List, Tuple 
from collections import defaultdict

logger = logging.getLogger(__name__)


def run_ragtag_analysis(fasta1_path: str, fasta2_path: str, output_dir: str) -> Dict:

    ragtag_dir = os.path.join(output_dir, 'ragtag_alignment')
    os.makedirs(ragtag_dir, exist_ok=True)
    

    cmd = ['ragtag.py', 'scaffold', fasta1_path, fasta2_path, '-o', ragtag_dir]
    

    try:
        subprocess.run(cmd, capture_output=True, text=True, check=True)
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"RagTag failed: {e.stderr}")
    
    agp_file = os.path.join(ragtag_dir, 'ragtag.scaffold.agp')
    if not os.path.exists(agp_file):
        raise FileNotFoundError(f"Expected AGP file not found: {agp_file}")
    

    alignments = defaultdict(list)
    with open(agp_file, 'r') as f:
        for line in f:
            if line.startswith('#') or len(line.strip().split('\t')) < 9:
                continue
            
            parts = line.strip().split('\t')
            if parts[4] == 'W':
                ref_chr = parts[0].replace('_RagTag', '')
                query_chr = parts[5]
                alignment_length = int(parts[2]) - int(parts[1]) + 1
                alignments[query_chr].append((ref_chr, alignment_length))
    

    chromosome_mapping = {}
    for query_chr, alignment_list in alignments.items():
        if alignment_list:
            best_ref_chr = max(alignment_list, key=lambda x: x[1])[0]
            chromosome_mapping[query_chr] = best_ref_chr
    
    return {
        'chromosome_mapping': chromosome_mapping,
        'agp_file': agp_file
    }


def get_optimal_chromosome_pairs(df1: pd.DataFrame, df2: pd.DataFrame, 
                                chromosome_mapping: Dict = None) -> List[Tuple[str, str]]:

    chr_names_1 = set(df1['sequence'])
    chr_names_2 = set(df2['sequence'])
    

    common_chrs = chr_names_1 & chr_names_2
    compatibility_threshold = 0.1 * min(len(chr_names_1), len(chr_names_2))
    
    if len(common_chrs) >= compatibility_threshold:
        logger.info("Using direct chromosome name matching")
        return [(chr, chr) for chr in common_chrs]

    if chromosome_mapping:
        candidate_pairs = [
            (ref_chr, query_chr) for query_chr, ref_chr in chromosome_mapping.items()
            if query_chr in chr_names_2 and ref_chr in chr_names_1
        ]
    else:
        logger.warning("No chromosome mapping provided - trying all combinations")
        candidate_pairs = [(chr1, chr2) for chr1 in chr_names_1 for chr2 in chr_names_2]
    

    scored_pairs = []
    for chr1, chr2 in candidate_pairs:
        chr1_genes = set(df1[df1['sequence'] == chr1]['busco_id'])
        chr2_genes = set(df2[df2['sequence'] == chr2]['busco_id'])
        common_genes = len(chr1_genes & chr2_genes)
        
        if common_genes > 10: 
            scored_pairs.append((chr1, chr2, common_genes))
    

    scored_pairs.sort(key=lambda x: x[2], reverse=True)
    
    used_chr1 = set()
    used_chr2 = set()
    final_pairs = []
    
    high_confidence_threshold = 300
    
    for chr1, chr2, count in scored_pairs:
        if (chr1 not in used_chr1 and chr2 not in used_chr2) or count > high_confidence_threshold:
            final_pairs.append((chr1, chr2))
            if count <= high_confidence_threshold:
                used_chr1.add(chr1)
                used_chr2.add(chr2)
    
    return final_pairs


def extract_syntenic_data(df1: pd.DataFrame, df2: pd.DataFrame, 
                         chromosome_pairs: List[Tuple[str, str]]) -> Tuple[pd.DataFrame, pd.DataFrame]:

    syntenic_parts_1 = []
    syntenic_parts_2 = []
    
    for chr1, chr2 in chromosome_pairs:
        chr1_genes = df1[df1['sequence'] == chr1]
        chr2_genes = df2[df2['sequence'] == chr2]
        
        common_genes = set(chr1_genes['busco_id']) & set(chr2_genes['busco_id'])
        
        if len(common_genes) > 0:
            syntenic_chr1 = chr1_genes[chr1_genes['busco_id'].isin(common_genes)]
            syntenic_chr2 = chr2_genes[chr2_genes['busco_id'].isin(common_genes)]
            
            syntenic_parts_1.append(syntenic_chr1)
            syntenic_parts_2.append(syntenic_chr2)
    
    if not syntenic_parts_1:
        raise ValueError("No syntenic genes found between genomes")
    
    syntenic_df1 = pd.concat(syntenic_parts_1, ignore_index=True)
    syntenic_df2 = pd.concat(syntenic_parts_2, ignore_index=True)
    
    return syntenic_df1, syntenic_df2


def complete_chromosome_analysis(df1: pd.DataFrame, df2: pd.DataFrame, 
                                config: Dict) -> Tuple[pd.DataFrame, pd.DataFrame, List[Tuple[str, str]]]:

    chromosome_mapping = None
    if 'first_fasta_path' in config and 'second_fasta_path' in config:
        try:
            ragtag_results = run_ragtag_analysis(
                config['first_fasta_path'], 
                config['second_fasta_path'], 
                config.get('output_dir', '.')
            )
            chromosome_mapping = ragtag_results['chromosome_mapping']
            logger.info("RagTag analysis completed successfully")
        except Exception as e:
            logger.warning(f"RagTag analysis failed: {e}")
    

    chromosome_pairs = get_optimal_chromosome_pairs(df1, df2, chromosome_mapping)
    
    if not chromosome_pairs:
        raise ValueError("No suitable chromosome pairs found")
    

    syntenic_df1, syntenic_df2 = extract_syntenic_data(df1, df2, chromosome_pairs)
    
    
    return syntenic_df1, syntenic_df2, chromosome_pairs