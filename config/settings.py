"""
Configuration settings
"""

CONFIG = {
    'first_busco_path': 'template/full_table.tsv',
    'second_busco_path': 'template/full_table.tsv',
    'first_species_name': 'First_Species',
    'second_species_name': 'Second_Species',
    'busco_status_filter': ['Complete'],

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
    
    'inversion_metrics': ['raw_count', 'normalized_rate', 'per_mb'],
    
    'position_bin_size_kb': 1.5,
    'profile_calculation_method': 'average',  # 'average', or 'range',
    'probability_threshold_for_target': 0.3,
    
    'permutable_positions_threshold': 0.5,  # ...Minimum probability for permutable positions
    'max_permutable_positions': 10,  # ...Maximum number of permutable positions to consider
    
    'dotplot_config': {
        'figure_size': (12, 10),
        'dpi': 600,
        'point_size': 8,
        'point_alpha': 0.8,
        'chr_line_color': 'grey',
        'chr_line_width': 0.8,
        'font_size_labels': 14,
        'font_size_title': 18,
        'font_size_chr': 12,
        'show_legend': False,
        'margin_factor': 0.02,
        'color_palette': {
            "M1": "#1573af", "M2": "#e59d38", "M3": "#f0e354", 
            "M4": "#169e73", "M5": "#60b5e1", "M6": "#000000", 
            "unassigned": "#808080"
        }
    }
}
