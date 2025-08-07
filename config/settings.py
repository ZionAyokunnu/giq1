"""
Configuration settings
"""

CONFIG = {

    'first_fasta_path': 'template.fna',
    'second_fasta_path': 'template.fna',
    'first_busco_path': 'template/full_table.tsv',
    'second_busco_path': 'template/full_table.tsv',
    'first_species_name': 'First_Species',
    'second_species_name': 'Second_Species',
    'enable_duplicate_handling': False ,
    'busco_status_filter': ['Complete', 'Duplicated'],

    'base_output_dir': 'results',
    'inversion_summary_csv': 'inversion_summary.csv',
    'chromosome_rearrangements_csv': 'chromosome_rearrangements.csv',
    
    
    'plot_width': 15,
    'plot_height': 10,
    'dpi': 300,
    'color_palette': 'viridis',
    'font_size': 12,
    'figure_format': 'png',
    
    'generate_dotplots': True,          
    'dotplot_show_labels': False,        
    'dotplot_by': 'busco',             
    'dotplot_size': (12, 10),            
    'synteny_color': '#1f77b4',          # Blue for syntenic regions
    'inversion_color': '#d62728',        # Red for inverted regions
    'confidence_alpha': True,            
    'show_synteny_blocks': True,         
    'show_breakpoints': True,          

    'use_synteny_plotter': True,
    'create_curved_ribbons': True,
    'create_chord_diagrams': True,
    
    'use_existing_diptera_tree': True,
    'diptera_tree_path': '/path/to/diptera.nwk',
    'prune_to_species': True,
    'annotate_with_inversions': True,
    
    'inversion_metrics': ['raw_count', 'normalized_rate', 'per_mb'],
    'edge_annotation_style': 'heatmap',  # 'heatmap', 'width', 'color'
    
    'create_publication_plots': True,
    'plot_formats': ['png', 'pdf', 'svg'],
    'high_dpi': 300,
    
    'external_tools': {
        'synteny_plotter': '/Users/za7/Documents/giq/external_tools/synteny_plotter'
    },
    
    'synteny_visualisation': {
        'enabled': True,
        'create_curved_ribbons': True,
        'create_straight_links': False,
        'options': {
            'filter_threshold': 3,  # Minimum BUSCOs per chromosome (R script -f parameter)
            'gap': 6,              # Gap between chromosomal sets (R script -g parameter)
            'alpha': 0,           # Transparency percentage (R script -alpha parameter)
        }
    },
    
    'tree_annotation': {
        'enabled': True,
        'source_tree_path': 'diptera_clean_20species.newick',  #to be changed with real tree
        'prune_to_target_species': True,
        'annotation_metrics': {
            'inversion_count': True,
            'inversion_rate_per_mb': True,
            'normalized_inversion_score': True
        },
        'visualisation': {
            'edge_annotation_style': 'heatmap',  # 'heatmap', 'width', 'color', 'labels'
            'node_support_values': True,
            'output_formats': ['png', 'pdf', 'newick'],
            'tree_layout': 'rectangular',  # 'rectangular', 'circular'
            'dpi': 300
        }
    },
    
    'publication_suite': {
        'enabled': True,
        'create_all_plots': True,
        'high_quality_output': True,
        'create_supplementary_data': True
    },

    'busco_phylogeny': {
        'enabled': True, 
        'approach': 'existing_tree',
        'source_tree_path': 'diptera_clean_20species.newick'
    },
    
    'syri_colors': {
    'syn': '#1f77b4',    # Blue for syntenic regions
    'inv': '#d62728',    # Red for inverted regions
    'trans': '#ff7f0e',  # Orange for translocations
    'dup': '#2ca02c',    # Green for duplications
    'del': '#9467bd'     # Purple for deletions
    },
    
    'position_bin_size_kb': 1.5,
    'profile_calculation_method': 'average',  # 'average', or 'range',
    'probability_threshold_for_target': 0.3,
    
    'permutable_positions_threshold': 0.5,  # ...Minimum probability for permutable positions
    'max_permutable_positions': 10  # ...Maximum number of permutable positions to consider
}
