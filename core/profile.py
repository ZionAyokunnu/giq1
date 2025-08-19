

import pandas as pd
import logging
from scipy.stats import pearsonr
from collections import defaultdict

from config import (
    CONFIG
)

logger = logging.getLogger(__name__)


def assign_ordinal_ranks_within_chromosome(corrected_df):
    """
    Assign ordinal ranks to genes within a chromosome based on start position.
    
    Args:
        corrected_df: DataFrame with gene information for a single chromosome
        
    Returns:
        DataFrame: corrected_df with added 'ordinal_rank' column
    """
    # Sort genes by start position within chromosome
    sorted_df = corrected_df.sort_values('gene_start').reset_index(drop=True)
    
    # Assign consecutive ordinal ranks (1, 2, 3, ...)
    sorted_df['ordinal_rank'] = range(1, len(sorted_df) + 1)
    
    return sorted_df


def calculate_gene_bin_overlaps_hybrid(gene_start, gene_end, chromosome, bin_size_bp, 
                                     ordinal_rank, chromosome_mappings=None, genome_id=None):

    gene_length = gene_end - gene_start
    if gene_length <= 0:
        return {
            'positional_bins': [],
            'ordinal_window': None,
            'gene_position': gene_start,
            'relative_position': 0.0
        }
    
    # Determine chromosome name to use
    if chromosome_mappings and genome_id:
        from core.chr_aligner import standardize_chromosome_name_unified
        standardized_chr = standardize_chromosome_name_unified(chromosome_mappings, genome_id, chromosome)
        bin_chromosome = standardized_chr
    else:
        bin_chromosome = chromosome
    
    # Calculate positional bins (existing logic)
    start_bin = gene_start // bin_size_bp
    end_bin = gene_end // bin_size_bp
    
    positional_bins = []
    for bin_num in range(start_bin, end_bin + 1):
        bin_start = bin_num * bin_size_bp
        bin_end = (bin_num + 1) * bin_size_bp
        
        overlap_start = max(gene_start, bin_start)
        overlap_end = min(gene_end, bin_end)
        overlap_length = max(0, overlap_end - overlap_start)
        
        if overlap_length > 0:
            overlap_percentage = (overlap_length / gene_length) * 100
            bin_id = f"{bin_chromosome}_bin_{bin_num}"
            positional_bins.append((bin_id, overlap_percentage))
    
    # Create ordinal window
    ordinal_window = f"{bin_chromosome}_rank_{ordinal_rank}"
    
    # Calculate relative position (we'll need chromosome length for this - placeholder for now)
    relative_position = 0.0  # TODO: Calculate actual relative position when chromosome length is available
    
    return {
        'positional_bins': positional_bins,
        'ordinal_window': ordinal_window,
        'gene_position': gene_start,
        'relative_position': relative_position
    }


def assign_genes_to_bins_hybrid(corrected_df, bin_size_kb=None, chromosome_mappings=None, genome_id=None):
    """
    Assign genes to both positional bins and ordinal windows.
    
    Args:
        corrected_df: DataFrame with gene information
        bin_size_kb: Bin size in kilobases
        chromosome_mappings: Optional chromosome mappings from align-chr
        genome_id: Genome identifier
        
    Returns:
        dict: {busco_id: {
            'positional_bins': [(bin_id, overlap_percentage), ...],
            'ordinal_window': ordinal_window_id,
            'gene_position': gene_start,
            'relative_position': relative_position_percentage
        }}
    """
    if bin_size_kb is None:
        bin_size_kb = CONFIG['position_bin_size_kb']
    
    bin_size_bp = bin_size_kb * 1000
    
    # First, assign ordinal ranks within this chromosome
    ranked_df = assign_ordinal_ranks_within_chromosome(corrected_df)
    
    gene_assignments = {}
    
    for _, gene in ranked_df.iterrows():
        busco_id = gene['busco_id']
        chromosome = gene['sequence']
        gene_start = gene['gene_start']
        gene_end = gene['gene_end']
        ordinal_rank = gene['ordinal_rank']
        
        # Calculate both positional and ordinal assignments
        hybrid_data = calculate_gene_bin_overlaps_hybrid(
            gene_start, gene_end, chromosome, bin_size_bp, ordinal_rank,
            chromosome_mappings, genome_id
        )
        
        gene_assignments[busco_id] = hybrid_data
    
    return gene_assignments

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
            
            grouped_genomes[genome_id][chromosome] = chromosome_df.reset_index(drop=True)
            print(f"  {chromosome}: {len(chromosome_df)} genes")
    
    return grouped_genomes


def process_genomes_binning_hybrid(grouped_genomes, bin_size_kb=None, chromosome_mappings=None):
    """
    Process multiple genomes for hybrid bin assignments, grouped by chromosome.
    
    Args:
        grouped_genomes: {genome_id: {chromosome: DataFrame}}
        bin_size_kb: Bin size in kilobases
        chromosome_mappings: Optional chromosome mappings from align-chr
        
    Returns:
        dict: {genome_id: {chromosome: {busco_id: hybrid_data}}}
    """
    if bin_size_kb is None:
        bin_size_kb = CONFIG['position_bin_size_kb']
    
    genome_assignments = {}
    
    for genome_id, chromosomes in grouped_genomes.items():
        genome_assignments[genome_id] = {}
        
        for chromosome, corrected_df in chromosomes.items():
            print(f"  Processing {genome_id} - {chromosome} (hybrid)")
            
            # Get hybrid assignments for this chromosome
            hybrid_assignments = assign_genes_to_bins_hybrid(
                corrected_df, bin_size_kb, chromosome_mappings, genome_id
            )
            
            genome_assignments[genome_id][chromosome] = hybrid_assignments
            print(f"    Assigned {len(hybrid_assignments)} genes to hybrid bins/ranks")
            
            # Debug: Show sample assignments
            if hybrid_assignments:
                sample_gene = next(iter(hybrid_assignments.keys()))
                sample_data = hybrid_assignments[sample_gene]
                print(f"    Sample positional: {sample_data['positional_bins'][0] if sample_data['positional_bins'] else 'None'}")
                print(f"    Sample ordinal: {sample_data['ordinal_window']}")
    
    return genome_assignments


def build_markov_profile_hybrid(genome_assignments, calculation_method=None):
    """
    Build hybrid Markov profile with both positional and ordinal data.
    
    Args:
        genome_assignments: Output from process_genomes_binning_hybrid()
        calculation_method: 'average' or 'range'
        
    Returns:
        dict: {
            'positional_profile': {bin_id: {busco_id: profile_data}},
            'ordinal_profile': {rank_id: {busco_id: profile_data}},
            'hybrid_summary': summary_stats
        }
    """
    if calculation_method is None:
        calculation_method = CONFIG['profile_calculation_method']
    
    # Separate data collection for positional and ordinal
    positional_data = defaultdict(lambda: defaultdict(list))
    ordinal_data = defaultdict(lambda: defaultdict(list))
    
    # Process each genome's assignments
    for genome_id, chromosomes in genome_assignments.items():
        for chromosome, gene_assignments in chromosomes.items():
            for busco_id, hybrid_data in gene_assignments.items():
                
                # Collect positional data
                for bin_id, overlap_percentage in hybrid_data['positional_bins']:
                    positional_data[bin_id][busco_id].append(overlap_percentage)
                
                # Collect ordinal data (100% presence in ordinal window)
                ordinal_window = hybrid_data['ordinal_window']
                if ordinal_window:
                    ordinal_data[ordinal_window][busco_id].append(100.0)  # 100% presence in rank
    
    # Build positional profile
    positional_profile = {}
    for bin_id, genes_data in positional_data.items():
        positional_profile[bin_id] = {}
        
        for busco_id, percentages in genes_data.items():
            genome_count = len(percentages)
            total_genomes = len(genome_assignments)
            
            avg_percentage = sum(percentages) / len(percentages)
            min_percentage = min(percentages)
            max_percentage = max(percentages)
            
            positional_profile[bin_id][busco_id] = {
                'average_percentage': round(avg_percentage, 2),
                'percentage_range': (round(min_percentage, 2), round(max_percentage, 2)),
                'genome_frequency': f"{genome_count}/{total_genomes}",
                'genome_count': genome_count,
                'total_genomes': total_genomes,
                'calculation_method': calculation_method,
                'profile_type': 'positional'
            }
    
    # Build ordinal profile
    ordinal_profile = {}
    for rank_id, genes_data in ordinal_data.items():
        ordinal_profile[rank_id] = {}
        
        for busco_id, presences in genes_data.items():
            genome_count = len(presences)
            total_genomes = len(genome_assignments)
            
            # For ordinal, we care about frequency of appearance at this rank
            frequency_percentage = (genome_count / total_genomes) * 100
            
            ordinal_profile[rank_id][busco_id] = {
                'frequency_percentage': round(frequency_percentage, 2),
                'genome_frequency': f"{genome_count}/{total_genomes}",
                'genome_count': genome_count,
                'total_genomes': total_genomes,
                'calculation_method': calculation_method,
                'profile_type': 'ordinal'
            }
    
    # Create hybrid summary
    all_genes_positional = set()
    for genes_data in positional_profile.values():
        all_genes_positional.update(genes_data.keys())
    
    all_genes_ordinal = set()
    for genes_data in ordinal_profile.values():
        all_genes_ordinal.update(genes_data.keys())
    
    hybrid_summary = {
        'total_genomes': len(genome_assignments),
        'positional_stats': {
            'total_bins': len(positional_profile),
            'total_genes': len(all_genes_positional),
            'gene_bin_combinations': sum(len(genes) for genes in positional_profile.values())
        },
        'ordinal_stats': {
            'total_ranks': len(ordinal_profile),
            'total_genes': len(all_genes_ordinal),
            'gene_rank_combinations': sum(len(genes) for genes in ordinal_profile.values())
        },
        'overlap_genes': len(all_genes_positional & all_genes_ordinal),
        'profile_type': 'hybrid'
    }
    
    print(f"Built hybrid profile:")
    print(f"  Positional: {hybrid_summary['positional_stats']['total_bins']} bins, {hybrid_summary['positional_stats']['total_genes']} genes")
    print(f"  Ordinal: {hybrid_summary['ordinal_stats']['total_ranks']} ranks, {hybrid_summary['ordinal_stats']['total_genes']} genes")
    print(f"  Overlapping genes: {hybrid_summary['overlap_genes']}")
    
    return {
        'positional_profile': positional_profile,
        'ordinal_profile': ordinal_profile,
        'hybrid_summary': hybrid_summary
    }


# Backward compatibility function
def build_markov_profile(genome_bin_assignments, calculation_method=None):
    """
    Enhanced build_markov_profile that detects input type and builds appropriate profile.
    
    Args:
        genome_bin_assignments: Either old format or new hybrid format
        calculation_method: 'average' or 'range'
        
    Returns:
        dict: Profile data (hybrid if input is hybrid, otherwise legacy format)
    """
    if calculation_method is None:
        calculation_method = CONFIG['profile_calculation_method']
    
    # Detect if input is hybrid format by checking structure
    sample_genome = next(iter(genome_bin_assignments.values()))
    sample_chromosome = next(iter(sample_genome.values()))
    sample_gene_data = next(iter(sample_chromosome.values()))
    
    # Check if it's hybrid format (contains 'positional_bins' and 'ordinal_window')
    if isinstance(sample_gene_data, dict) and 'positional_bins' in sample_gene_data:
        print("Detected hybrid input format - building hybrid profile")
        return build_markov_profile_hybrid(genome_bin_assignments, calculation_method)
    else:
        print("Detected legacy input format - building legacy profile")
        # Use existing legacy logic
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
                    'calculation_method': calculation_method,
                    'profile_type': 'legacy_positional'
                }
        
        return markov_profile


def get_profile_summary_hybrid(hybrid_profile):
    """
    Get summary statistics for hybrid profile.
    
    Args:
        hybrid_profile: Output from build_markov_profile_hybrid()
        
    Returns:
        dict: Summary statistics
    """
    if 'hybrid_summary' in hybrid_profile:
        return hybrid_profile['hybrid_summary']
    else:
        
        
        
        ########
        # Legacy profile
        return get_profile_summary(hybrid_profile)


def display_profile_sample_hybrid(hybrid_profile, sample_items=3):
    """
    Display sample of hybrid profile for inspection.
    
    Args:
        hybrid_profile: Output from build_markov_profile_hybrid()
        sample_items: Number of items to display for each type
    """
    print("Hybrid Profile Sample:")
    print("=" * 60)
    
    if 'positional_profile' in hybrid_profile:
        print("\nPOSITIONAL PROFILE:")
        print("-" * 30)
        positional_bins = list(hybrid_profile['positional_profile'].keys())[:sample_items]
        
        for bin_id in positional_bins:
            print(f"\n{bin_id}:")
            genes_data = hybrid_profile['positional_profile'][bin_id]
            sample_genes = list(genes_data.keys())[:3]
            for gene_id in sample_genes:
                gene_data = genes_data[gene_id]
                print(f"  {gene_id}: avg={gene_data['average_percentage']}%, freq={gene_data['genome_frequency']}")
        
        print("\nORDINAL PROFILE:")
        print("-" * 30)
        ordinal_ranks = list(hybrid_profile['ordinal_profile'].keys())[:sample_items]
        
        for rank_id in ordinal_ranks:
            print(f"\n{rank_id}:")
            genes_data = hybrid_profile['ordinal_profile'][rank_id]
            sample_genes = list(genes_data.keys())[:3]
            for gene_id in sample_genes:
                gene_data = genes_data[gene_id]
                print(f"  {gene_id}: freq={gene_data['frequency_percentage']}%, appears={gene_data['genome_frequency']}")
        
        print(f"\nHybrid Summary:")
        summary = hybrid_profile['hybrid_summary']
        print(f"  Total genomes: {summary['total_genomes']}")
        print(f"  Positional bins: {summary['positional_stats']['total_bins']}")
        print(f"  Ordinal ranks: {summary['ordinal_stats']['total_ranks']}")
        print(f"  Overlapping genes: {summary['overlap_genes']}")
    else:
     
     
     
     
     
     ########
        display_profile_sample(hybrid_profile, sample_items)


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


