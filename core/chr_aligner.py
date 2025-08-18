import sys
import os
import subprocess
import logging
from pathlib import Path
import pandas as pd
import numpy as np
from typing import Dict, Optional
from collections import defaultdict

logger = logging.getLogger(__name__)

def parse_fasta_headers(fasta_path):
    chromosomes = set()
    
    if not os.path.exists(fasta_path):
        raise FileNotFoundError(f"FASTA file not found: {fasta_path}")
    
    with open(fasta_path, 'r') as f:
        for line in f:
            if line.startswith('>'):
                chr_name = line.strip()[1:].split()[0]
                chromosomes.add(chr_name)
    
    return chromosomes


def parse_agp_alignments(agp_file_path):
    alignments = defaultdict(list)
    
    with open(agp_file_path, 'r') as f:
        for line in f:
            if line.startswith('#') or line.startswith('##'):
                continue
            
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            
            if parts[4] == 'W':
                ref_chr = parts[0].replace('_RagTag', '')
                query_chr = parts[5]
                start_pos = int(parts[1])
                end_pos = int(parts[2])
                alignment_length = end_pos - start_pos + 1
                
                alignments[query_chr].append((ref_chr, alignment_length))
    
    return alignments


def create_chromosome_mappings(agp_file_path, reference_fasta_path, query_fasta_path):
    alignments = parse_agp_alignments(agp_file_path)
    
    query_to_ref = {}
    for query_chr, alignment_list in alignments.items():
        if alignment_list:
            largest_alignment = max(alignment_list, key=lambda x: x[1])
            query_to_ref[query_chr] = largest_alignment[0]
    
    return {
        'query_to_reference': query_to_ref,
        'alignment_details': alignments
    }


def run_ragtag_scaffold(reference_fasta, query_fasta, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    
    cmd = [
        'ragtag.py', 'scaffold',
        reference_fasta, query_fasta,
        '-o', output_dir
    ]
    
    logger.info(f"Running RagTag: {os.path.basename(reference_fasta)} vs {os.path.basename(query_fasta)}")
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"RagTag failed: {e.stderr}")
    
    agp_file = os.path.join(output_dir, 'ragtag.scaffold.agp')
    if not os.path.exists(agp_file):
        raise FileNotFoundError(f"Expected AGP file not found: {agp_file}")
    
    return agp_file


def run_pairwise_ragtag(config):
    
    fasta1_path = config['first_fasta_path']
    fasta2_path = config['second_fasta_path']
    
    output_dir = os.path.join(config.get('output_dir', '.'), 'ragtag_alignment')
    
    # Run RagTag scaffold (fasta1 as reference, fasta2 as query)
    agp_file = run_ragtag_scaffold(fasta1_path, fasta2_path, output_dir)
    
    # Create chromosome mappings
    mappings = create_chromosome_mappings(agp_file, fasta1_path, fasta2_path)
    
    logger.info(f"RagTag mapped {len(mappings['query_to_reference'])} chromosomes")


    return mappings['query_to_reference']



def select_optimal_chromosome_pairs(chromosome_pairs_with_counts):
    """Select the best 1:1 chromosome pairing for diagonal synteny"""
    
    # Sort by syntenic gene count (highest first)
    sorted_pairs = sorted(chromosome_pairs_with_counts, key=lambda x: x[2], reverse=True)
    
    # For each genome1 chromosome, find best genome2 match
    chr1_to_best_chr2 = {}
    used_chr2 = set()
    
    for chr1, chr2, count in sorted_pairs:
        if chr1 not in chr1_to_best_chr2 and chr2 not in used_chr2:
            chr1_to_best_chr2[chr1] = chr2
            used_chr2.add(chr2)
    
    optimal_pairs = list(chr1_to_best_chr2.items())
    return optimal_pairs

def prepare_diagonal_synteny_data(joined_df, optimal_pairs):
    """Prepare data for diagonal synteny plot"""
    
    # Filter to only include optimal pairs
    diagonal_data = []
    
    for chr1, chr2 in optimal_pairs:
        pair_data = joined_df[
            (joined_df['chr1'] == chr1) & (joined_df['chr2'] == chr2)
        ].copy()
        diagonal_data.append(pair_data)
    
    if diagonal_data:
        diagonal_df = pd.concat(diagonal_data, ignore_index=True)
    else:
        diagonal_df = joined_df  # Fallback
    
    # Create ordered chromosome lists for plotting
    chr1_order = [pair[0] for pair in optimal_pairs]
    chr2_order = [pair[1] for pair in optimal_pairs] 
    
    return diagonal_df, chr1_order, chr2_order