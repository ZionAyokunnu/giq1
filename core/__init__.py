"""Core analysis modules for genome inversion analyser"""

# Import from quality_assessment
from .quality_assessment import (
    assess_assembly_quality,
    calculate_comprehensive_metrics,
    print_detailed_metrics,
    generate_iteration_report,
    save_pattern_analysis
    
)

# Import from busco_processor
from .busco_processor import (
    parse_busco_table,
    filter_busco_genes,
    detect_flips,
    correct_chromosome_orientation,
    track_strand_changes_per_iteration,
    generate_strand_debug_tsv
)

# Import from profile
from .profile import (
    assign_ordinal_ranks_within_chromosome,
    calculate_gene_bin_overlaps_hybrid,
    assign_genes_to_bins_hybrid,
    process_genomes_binning_hybrid,
    build_markov_profile_hybrid,
    build_markov_profile,
    get_profile_summary_hybrid,
    display_profile_sample_hybrid,
    display_profile_sample,
    get_profile_summary,
    group_genomes_by_chromosome,
    filter_chromosomes_for_genome
)

# Import from query_reference
from .query_reference import (
    extract_gene_distribution,
    calculate_distribution_stats,
    calculate_bit_scores,
    compare_query_genome_to_profile,
)

# Import from query_movement
from .query_movement import (
    extract_gene_ranges,
    extract_current_ranges,
    filter_same_chromosome,
    calculate_gene_movement,
    analyse_query_movements,
    get_movement_summary,
    convert_hybrid_to_legacy_format,
    extract_gene_distribution_hybrid_aware,
)

# Import from reverse
from .reverse import (
    extract_movement_sequence,
    detect_adjacency_inversions,
    detect_extended,
    detect_flip_in_pattern,
    apply_adjacency_inversion,
    apply_flip_inversion,
    iterative_detection,
    get_permutable_positions,
    calculate_position_probability,
    evaluate_inversion_step_probability,
    generate_gene_specific_pathways,
    generate_inversion_steps_for_gene,
    iterative_detection_gene_specific,
    evaluate_pathway_steps,
    get_gene_final_position,
    track_all_gene_final_positions,
    probability_weighted_inversion_analysis,
    check_events_iteration,
    create_pairwise_movement_sequence_per_chromosome,
    integrate_best_alternatives,
    extract_chromosome_movement_sequence,
    calculate_total_movement,
    find_non_overlapping_adjacencies,
    find_non_overlapping_flips,
    is_perfect_incremental,
)

# Import from chr_aligner
from .chr_aligner import (
    normalize_genome_name,
    get_genome_name_from_fasta,
    parse_fasta_headers,
    parse_agp_alignments,
    create_chromosome_mappings,
    run_ragtag_scaffold,
    run_all_pairwise_alignments,
    save_unified_mappings,
    find_connected_components,
    find_all_homologous_groups,
    get_total_gene_content,
    assign_hierarchical_standard_names,
    create_unified_mappings_multi_reference,
    load_chromosome_mappings,
    standardize_chromosome_name_unified,
    build_gene_content_homology_primary,
    enhance_homology_with_gene_content,
    filter_chromosomes_by_gene_count,
    filter_unmapped_by_gene_count,
)

# Import from transposition
from .transposition import (
    detect_and_apply_translocations,
    calculate_total_movement as calculate_total_movement_transposition,
    detect_translocation_patterns,
    group_into_translocation_blocks,
    find_optimal_insertion_point,
    apply_translocation,
)

# Import from profile_builder
from .profile_builder import (
    build_profile_command_hybrid,
    multi_reference_alignment_pipeline,
    save_stage_data,
    save_hybrid_profile_data,
    save_hybrid_assignments,
    standardize_chromosome_names_in_binning,
)

# Import from query_analyzer
from .query_analyzer import (
    analyze_query_command_hybrid,
)

# Import from convergence_analysis
from .convergence_analysis import (
    create_single_convergence_tsv,
    create_convergence_analysis_tsv,
    analyze_transposition_patterns,
)

from .formats import (
    run_algorithm_test,
)


# Define __all__ to control what gets imported with "from core import *"
__all__ = [
    # Quality assessment
    'assess_assembly_quality',

    # BUSCO processing
    'parse_busco_table',
    'filter_busco_genes',
    'detect_flips',
    'correct_chromosome_orientation',
    
    # Profile building
    'assign_ordinal_ranks_within_chromosome',
    'calculate_gene_bin_overlaps_hybrid',
    'assign_genes_to_bins_hybrid',
    'process_genomes_binning_hybrid',
    'build_markov_profile_hybrid',
    'build_markov_profile',
    'get_profile_summary_hybrid',
    'display_profile_sample_hybrid',
    'display_profile_sample',
    'get_profile_summary',
    'group_genomes_by_chromosome',
    'filter_chromosomes_for_genome',
    
    # Query reference
    'extract_gene_distribution',
    'calculate_distribution_stats',
    'calculate_bit_scores',
    'compare_query_genome_to_profile',
    
    # Query movement
    'extract_gene_ranges',
    'extract_current_ranges',
    'filter_same_chromosome',
    'calculate_gene_movement',
    'analyse_query_movements',
    'get_movement_summary',
    'convert_hybrid_to_legacy_format',
    'extract_gene_distribution_hybrid_aware',
    
    # Reverse/inversion analysis
    'extract_movement_sequence',
    'detect_adjacency_inversions',
    'detect_extended',
    'detect_flip_in_pattern',
    'apply_adjacency_inversion',
    'apply_flip_inversion',
    'iterative_detection',
    'get_permutable_positions',
    'calculate_position_probability',
    'evaluate_inversion_step_probability',
    'generate_gene_specific_pathways',
    'generate_inversion_steps_for_gene',
    'iterative_detection_gene_specific',
    'evaluate_pathway_steps',
    'get_gene_final_position',
    'track_all_gene_final_positions',
    'probability_weighted_inversion_analysis',
    'check_events_iteration',
    'create_pairwise_movement_sequence_per_chromosome',
    'integrate_best_alternatives',
    'extract_chromosome_movement_sequence',
    'calculate_total_movement',
    'find_non_overlapping_adjacencies',
    'find_non_overlapping_flips',
    'is_perfect_incremental',
    
    # Chromosome alignment
    'normalize_genome_name',
    'get_genome_name_from_fasta',
    'parse_fasta_headers',
    'parse_agp_alignments',
    'create_chromosome_mappings',
    'run_ragtag_scaffold',
    'run_all_pairwise_alignments',
    'save_unified_mappings',
    'find_connected_components',
    'find_all_homologous_groups',
    'get_total_gene_content',
    'assign_hierarchical_standard_names',
    'create_unified_mappings_multi_reference',
    'load_chromosome_mappings',
    'standardize_chromosome_name_unified',
    'build_gene_content_homology_primary',
    'enhance_homology_with_gene_content',
    'filter_chromosomes_by_gene_count',
    'filter_unmapped_by_gene_count',
    
    # Transposition
    'detect_and_apply_translocations',
    'calculate_total_movement_transposition',
    'detect_translocation_patterns',
    'group_into_translocation_blocks',
    'find_optimal_insertion_point',
    'apply_translocation',
    
    # Profile builder
    'build_profile_command_hybrid',
    'multi_reference_alignment_pipeline',
    'save_stage_data',
    'save_hybrid_profile_data',
    'save_hybrid_assignments',
    'standardize_chromosome_names_in_binning',
    
    # Query analyzer
    'analyze_query_command_hybrid',
    
    # Convergence analysis
    'create_single_convergence_tsv',
    'create_convergence_analysis_tsv',
    'analyze_transposition_patterns',
]