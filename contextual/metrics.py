"""
BUSCO Contextual Metrics
"""

import pandas as pd
import numpy as np
import logging
from pathlib import Path
from typing import Dict, List, Any, Optional, Union
from collections import defaultdict

from config import (
    CONFIG
)

logger = logging.getLogger(__name__)


def compute_inversion_rate_per_mb_busco(inversion_df: pd.DataFrame, #Using Species A as the reference mb span
                                       joined_df: pd.DataFrame) -> Dict:

        
    total_inversions = len(inversion_df[inversion_df['flipped_genes'] > 0])
    total_flipped_genes = inversion_df['flipped_genes'].sum()
    total_genes = len(joined_df)
    
    genome_coverage_estimate = 0
    for _, row in inversion_df.iterrows():
        if row['total_genes'] > 1:
            chr_genes = joined_df[
                (joined_df['chr1'] == row['chr1']) & 
                (joined_df['chr2'] == row['chr2'])
            ]
            if len(chr_genes) > 0:
                chr1_span = chr_genes['end1'].max() - chr_genes['start1'].min()
                genome_coverage_estimate += chr1_span
    
    genome_size_mb = genome_coverage_estimate / 1_000_000 if genome_coverage_estimate > 0 else 0            #Print genome span for species 1 (reference Genome)
    
    if genome_size_mb == 0:
        logger.warning("Cannot estimate genome size from BUSCO data")
        return {}
    

    inversions_per_mb = total_inversions / genome_size_mb                           #Print
    flipped_genes_per_mb = total_flipped_genes / genome_size_mb                    #Print
    

    rates_by_type = {}
    if 'inversion_type' in inversion_df.columns:
        for inv_type in inversion_df['inversion_type'].unique():
            type_data = inversion_df[inversion_df['inversion_type'] == inv_type]
            type_inversions = len(type_data[type_data['flipped_genes'] > 0])
            rates_by_type[inv_type] = type_inversions / genome_size_mb              #Print
    

    rates_by_chr = {}
    for chrom in inversion_df['chr1'].unique():
        chr_data = inversion_df[inversion_df['chr1'] == chrom]
        chr_inversions = len(chr_data[chr_data['flipped_genes'] > 0])
        rates_by_chr[chrom] = chr_inversions / genome_size_mb if genome_size_mb > 0 else 0      #Print
    
    rate_metrics = {
        'total_chromosome_pairs_with_inversions': total_inversions,
        'total_flipped_genes': total_flipped_genes,
        'total_analysed_genes': total_genes,
        'genome_coverage_mb': genome_size_mb,
        'inversion_events_per_mb': inversions_per_mb,
        'flipped_genes_per_mb': flipped_genes_per_mb,
        'rates_by_type': rates_by_type,
        'rates_by_chromosome': rates_by_chr,
        'global_flip_rate': total_flipped_genes / total_genes if total_genes > 0 else 0
    }
    
    logger.info(f"Inversion rate: {inversions_per_mb:.3f} events per Mb, "
               f"{flipped_genes_per_mb:.3f} flipped genes per Mb")
    return rate_metrics




def _analyse_gene_density_correlation(inversion_df: pd.DataFrame, 
                                    joined_df: pd.DataFrame, config) -> Dict:
<<<<<<< HEAD
    """
    analyse correlation between gene density (per mb) and inversion density
    """
=======
>>>>>>> origin
    
    species1_name = config.get('first_species_name', 'Species A')
    species2_name = config.get('second_species_name', 'Species B')

    species_a_gene_density = {}  
    species_b_gene_density = {} 

    for chrom in joined_df['chr1'].unique():
        chr_genes = joined_df[joined_df['chr1'] == chrom]
        
        if len(chr_genes) > 1:
            chr_span = chr_genes['end1'].max() - chr_genes['start1'].min()
            density = len(chr_genes) / (chr_span / 1000000) if chr_span > 0 else 0
            species_a_gene_density[chrom] = density             #Print
    
    for chrom in joined_df['chr2'].unique():
        chr_genes = joined_df[joined_df['chr2'] == chrom]
        
        if len(chr_genes) > 1:
            chr_span = chr_genes['end2'].max() - chr_genes['start2'].min()
            density = len(chr_genes) / (chr_span / 1000000) if chr_span > 0 else 0
            species_b_gene_density[chrom] = density
    
    species_a_inversion_density = {}  # Inversions of A relative to B
    species_b_inversion_density = {}  # Inversions of B relative to A
    
    for chrom in species_a_gene_density.keys():
        chr_inv_pairs = inversion_df[inversion_df['chr1'] == chrom]
        total_flipped = chr_inv_pairs['flipped_genes'].sum()
        
        chr_genes = joined_df[joined_df['chr1'] == chrom]
        chr_span = chr_genes['end1'].max() - chr_genes['start1'].min()
        span_mb = chr_span / 1000000 if chr_span > 0 else 1
        
        species_a_inversion_density[chrom] = total_flipped / span_mb
    

    for chrom in species_b_gene_density.keys():
        chr_inv_pairs = inversion_df[inversion_df['chr2'] == chrom]
        total_flipped = chr_inv_pairs['flipped_genes'].sum()
        
        chr_genes = joined_df[joined_df['chr2'] == chrom]
        chr_span = chr_genes['end2'].max() - chr_genes['start2'].min()
        span_mb = chr_span / 1000000 if chr_span > 0 else 1
        
        species_b_inversion_density[chrom] = total_flipped / span_mb
    
    correlations = {}
    
    if len(species_a_gene_density) > 2:
        common_chroms_a = set(species_a_gene_density.keys()) & set(species_a_inversion_density.keys())
        if len(common_chroms_a) > 2:
            gene_vals_a = [species_a_gene_density[chrom] for chrom in common_chroms_a]
            inv_vals_a = [species_a_inversion_density[chrom] for chrom in common_chroms_a]
            
            if np.std(gene_vals_a) > 0 and np.std(inv_vals_a) > 0:
                corr_a = np.corrcoef(gene_vals_a, inv_vals_a)[0, 1]
                correlations['species_a_within_genome'] = {
                    'correlation': corr_a,
                    'chromosomes': list(common_chroms_a),
                    'description': '{species1_name} gene density vs {species1_name} inversion density'
                }
    
    if len(species_b_gene_density) > 2:
        common_chroms_b = set(species_b_gene_density.keys()) & set(species_b_inversion_density.keys())
        if len(common_chroms_b) > 2:
            gene_vals_b = [species_b_gene_density[chrom] for chrom in common_chroms_b]
            inv_vals_b = [species_b_inversion_density[chrom] for chrom in common_chroms_b]
            
            if np.std(gene_vals_b) > 0 and np.std(inv_vals_b) > 0:
                corr_b = np.corrcoef(gene_vals_b, inv_vals_b)[0, 1]
                correlations['species_b_within_genome'] = {
                    'correlation': corr_b,
                    'chromosomes': list(common_chroms_b),
                    'description': '{species2_name} gene density vs {species2_name} inversion density'
                }

    
    result = {
        'gene_density_species_a': species_a_gene_density,
        'gene_density_species_b': species_b_gene_density,
        'inversion_density_species_a': species_a_inversion_density,
        'inversion_density_species_b': species_b_inversion_density,
    }
    
    return result