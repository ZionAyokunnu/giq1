"""
Profiling with bin for Markov's model.

"""

import pandas as pd
import logging
from scipy.stats import pearsonr
from collections import defaultdict

from config import (
    CONFIG
)

logger = logging.getLogger(__name__)

def group_genomes_by_chromosome(corrected_genomes):
    """
    Group corrected genome data by chromosome.
    
    Args:
        corrected_genomes: dict {genome_id: corrected_df}
        
    Returns:
        dict: {genome_id: {chromosome: corrected_df}}
    """
    grouped_genomes = {}
    
    for genome_id, corrected_df in corrected_genomes.items():
        print(f"\nGrouping {genome_id} by chromosomes:")
        
        # Group by 'sequence' column (chromosome)
        chromosome_groups = corrected_df.groupby('sequence')
        
        grouped_genomes[genome_id] = {}
        for chromosome, chromosome_df in chromosome_groups:
            grouped_genomes[genome_id][chromosome] = chromosome_df.copy()
            print(f"  {chromosome}: {len(chromosome_df)} genes")
    
    return grouped_genomes


def calculate_gene_bin_overlaps(gene_start, gene_end, chromosome, bin_size_bp):
    """
    Calculate which bins a gene overlaps and the percentage overlap in each bin.
    
    Returns:
        list: [(bin_id, overlap_percentage), ...]
    """
    gene_length = gene_end - gene_start
    if gene_length <= 0:
        return []
    
    start_bin = gene_start // bin_size_bp
    end_bin = gene_end // bin_size_bp
    
    bin_overlaps = []
    
    for bin_num in range(start_bin, end_bin + 1):
        bin_start = bin_num * bin_size_bp
        bin_end = (bin_num + 1) * bin_size_bp
        
        overlap_start = max(gene_start, bin_start)
        overlap_end = min(gene_end, bin_end)
        overlap_length = max(0, overlap_end - overlap_start)
        
        if overlap_length > 0:
            overlap_percentage = (overlap_length / gene_length) * 100
            
            bin_id = f"{chromosome}_bin_{bin_num}"
            
            bin_overlaps.append((bin_id, overlap_percentage))
    
    return bin_overlaps


def assign_genes_to_bins(corrected_df, bin_size_kb=None):
    """
    Assign genes to genomic bins based on overlap percentages.
        
    Returns:
        dict: {busco_id: [(bin_id, overlap_percentage), ...]}
    """
    if bin_size_kb is None:
        bin_size_kb = CONFIG['position_bin_size_kb']
    
    bin_size_bp = bin_size_kb * 1000
    
    gene_bin_assignments = {}
    
    for _, gene in corrected_df.iterrows():
        busco_id = gene['busco_id']
        chromosome = gene['sequence']
        gene_start = gene['gene_start']
        gene_end = gene['gene_end']
        
        bin_overlaps = calculate_gene_bin_overlaps(
            gene_start, gene_end, chromosome, bin_size_bp
        )
        
        gene_bin_assignments[busco_id] = bin_overlaps
    
    return gene_bin_assignments


def process_genomes_binning(grouped_genomes, bin_size_kb=None):
    """
    Process multiple genomes for bin assignments, grouped by chromosome.
        
    Returns:
        dict: {genome_id: {chromosome: {busco_id: [(bin_id, overlap_percentage), ...]}}}
    """
    if bin_size_kb is None:
        bin_size_kb = CONFIG['position_bin_size_kb']
    
    genome_bin_assignments = {}
    
    for genome_id, chromosomes in grouped_genomes.items():
        genome_bin_assignments[genome_id] = {}
        
        for chromosome, corrected_df in chromosomes.items():
            print(f"  Processing {genome_id} - {chromosome}")
            gene_bins = assign_genes_to_bins(corrected_df, bin_size_kb)
            genome_bin_assignments[genome_id][chromosome] = gene_bins
            print(f"    Assigned {len(gene_bins)} genes to bins")
    
    return genome_bin_assignments


def build_markov_profile(genome_bin_assignments, calculation_method=None):
    """
    Build the Markov percentage profile from multiple genome bin assignments.
    Now handles chromosome-grouped data.
        
    Returns:
        dict: {bin_id: {busco_id: profile_data}}
    """
    if calculation_method is None:
        calculation_method = CONFIG['profile_calculation_method']
    
    bin_gene_data = defaultdict(lambda: defaultdict(list))
    
    # Handle nested structure: genome_id -> chromosome -> gene_bins
    for genome_id, chromosomes in genome_bin_assignments.items():
        for chromosome, gene_bins in chromosomes.items():
            for busco_id, bin_overlaps in gene_bins.items():
                for bin_id, overlap_percentage in bin_overlaps:
                    bin_gene_data[bin_id][busco_id].append(overlap_percentage)
    
    markov_profile = {}
    
    for bin_id, genes_data in bin_gene_data.items():
        markov_profile[bin_id] = {}
        
        for busco_id, percentages in genes_data.items():
            genome_count = len(percentages)
            total_genomes = len(genome_bin_assignments)
            
            avg_percentage = sum(percentages) / len(percentages)
            min_percentage = min(percentages)
            max_percentage = max(percentages)
            
            markov_profile[bin_id][busco_id] = {
                'average_percentage': round(avg_percentage, 2),
                'percentage_range': (round(min_percentage, 2), round(max_percentage, 2)),
                'genome_frequency': f"{genome_count}/{total_genomes}",
                'genome_count': genome_count,
                'total_genomes': total_genomes,
                'calculation_method': calculation_method
            }
    
    return markov_profile

#Needs more checking...
def get_profile_summary(markov_profile):
    """
    Get summary statistics of the Markov profile.
    
    Returns:
        dict: Summary statistics
    """
    total_bins = len(markov_profile)
    total_gene_bin_combinations = sum(len(genes) for genes in markov_profile.values())
    
    all_genes = set()
    for genes_data in markov_profile.values():
        all_genes.update(genes_data.keys())
    
    summary = {
        'total_bins': total_bins,
        'total_genes': len(all_genes),
        'total_gene_bin_combinations': total_gene_bin_combinations,
        'average_genes_per_bin': round(total_gene_bin_combinations / total_bins, 2) if total_bins > 0 else 0
    }
    
    return summary


def display_profile_sample(markov_profile, sample_bins=3):
    """
    Display a sample of the Markov profile for inspection.
    
    Args:
        markov_profile: Output from build_markov_profile()
        sample_bins: Number of bins to display
    """
    print("Markov Profile Sample:")
    print("=" * 50)
    
    bin_ids = list(markov_profile.keys())[:sample_bins]
    
    for bin_id in bin_ids:
        print(f"\n{bin_id}:")
        genes_data = markov_profile[bin_id]
        
        sample_genes = list(genes_data.keys())[:3]
        for gene_id in sample_genes:
            gene_data = genes_data[gene_id]
            print(f"  {gene_id}: {gene_data}")
        
        if len(genes_data) > 3:
            print(f"  ... and {len(genes_data) - 3} more genes")


if __name__ == "__main__":
    
    print("Profiling is working")
    
else:
    print("Profiling did not work")