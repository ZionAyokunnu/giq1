#!/usr/bin/env python3
"""
GIQ1 - Genome Inversion Quantifier
Main entry point.

Author: Zion Ayokunnu
Supervisors: Kamil Jaron, Sasha, Arif
Version: 1.0.0
"""

import sys
import logging
from pathlib import Path
import pandas as pd
import numpy as np
from typing import Dict, Optional, List, Tuple
from itertools import combinations
import argparse
from config import CONFIG
from utils import create_output_directory
from core import parse_busco_table, filter_busco_genes, detect_flips, complete_chromosome_analysis, correct_chromosome_orientation
from contextual import (
    compute_inversion_rate_per_mb_busco,
    _analyze_gene_density_correlation
)

from visualisation import (
    create_busco_dotplot,
    create_linearised_dotplot
)

logging.basicConfig(
    level=logging.INFO, 
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def parse_arguments() -> Tuple[List[str], Optional[List[str]]]:
    parser = argparse.ArgumentParser(description='GIQ1 - Genome Inversion Quantifier')
    parser.add_argument('busco_files', nargs='+', help='BUSCO table files (.tsv)')
    parser.add_argument('--fasta', nargs='+', help='FASTA files corresponding to BUSCO files (optional)')
    args = parser.parse_args()
    return args.busco_files, getattr(args, 'fasta', None)


def run_busco_inversion_analysis(config: Dict = None) -> Dict:
    if config is None:
        config = CONFIG
    
    output_dir = create_output_directory(config)
    logger.info(f"Output directory: {output_dir}")
    
    data_dir = output_dir / 'data'
    data_dir.mkdir(parents=True, exist_ok=True)
    
    species1_name = config.get('first_species_name', 'Species1')
    species2_name = config.get('second_species_name', 'Species2')
    
    try:
        busco_df1 = parse_busco_table(config['first_busco_path'], config)
        busco_df2 = parse_busco_table(config['second_busco_path'], config)
        
        filtered_df1 = filter_busco_genes(busco_df1, config, species_name=species1_name)
        filtered_df2 = filter_busco_genes(busco_df2, config, species_name=species2_name)
        
        syntenic_df1, syntenic_df2, chromosome_pairs = complete_chromosome_analysis(
            filtered_df1, filtered_df2, config
        )
        
        inversion_df, joined_df = detect_flips(syntenic_df1, syntenic_df2, config)
        
        corrected_df1, corrected_df2 = correct_chromosome_orientation(
            syntenic_df1, syntenic_df2, joined_df, config
        )
        
        final_inversion_df, final_joined_df = detect_flips(corrected_df1, corrected_df2, config)
        
        contextual_metrics = compute_contextual_metrics(
            final_inversion_df, final_joined_df, config
        )
        
        visualisation_results = create_analysis_visualisations(
            final_inversion_df, final_joined_df, output_dir, config
        )
        
        save_results(
            final_inversion_df, final_joined_df, contextual_metrics, output_dir, config
        )
        
        report_path = generate_analysis_report(
            final_inversion_df, final_joined_df, contextual_metrics, 
            output_dir, config
        )
        
        results = {
            'inversion_df': final_inversion_df,
            'joined_df': final_joined_df,
            'chromosome_pairs': chromosome_pairs,
            'contextual_metrics': contextual_metrics,
            'visualisation_results': visualisation_results,
            'output_dir': output_dir,
            'report_path': report_path,
            'config': config,
            'species_names': [species1_name, species2_name]
        }
        
        return results
        
    except Exception as e:
        logger.error(f"Analysis failed: {str(e)}")
        if config.get('enable_debug_output', False):
            import traceback
            traceback.print_exc()
        raise



def run_multispecies_analysis(busco_paths: List[str], fasta_paths: Optional[List[str]] = None) -> Dict:
    original_config = CONFIG.copy()
    
    base_output = Path(CONFIG.get('base_output_dir', './results'))
    base_output.mkdir(parents=True, exist_ok=True)
    
    pairs = list(combinations(busco_paths, 2))
    pairwise_results = {}
    
    for i, (path1, path2) in enumerate(pairs):
        species1_name = Path(path1).stem
        species2_name = Path(path2).stem
        
        try:
            pair_config = {}
            for key, value in CONFIG.items():
                if isinstance(value, dict):
                    pair_config[key] = value.copy()
                else:
                    pair_config[key] = value
            
            pair_config['first_busco_path'] = path1
            pair_config['second_busco_path'] = path2
            pair_config['first_species_name'] = species1_name
            pair_config['second_species_name'] = species2_name
            
            pair_dir = base_output / f'pairwise' / f'pair_{i}_{species1_name}_vs_{species2_name}'
            pair_dir.mkdir(parents=True, exist_ok=True)
            pair_config['base_output_dir'] = str(pair_dir)
            
            CONFIG.update({
                'first_species_name': species1_name,
                'second_species_name': species2_name,
                'base_output_dir': str(pair_dir)
            })
            
            if fasta_paths:
                fasta_dict = {Path(f).stem: f for f in fasta_paths}
                if species1_name in fasta_dict:
                    pair_config['first_fasta_path'] = fasta_dict[species1_name]
                if species2_name in fasta_dict:
                    pair_config['second_fasta_path'] = fasta_dict[species2_name]
            

            pairwise_results[(path1, path2)] = run_busco_inversion_analysis(pair_config)
            
        except Exception as e:
            logger.error(f"Pair {i} failed: {e}")
            continue
        
        finally:
            CONFIG.clear()
            CONFIG.update(original_config)
    

    # ancestral_results = reconstruct_ancestral_states(pairwise_results, busco_paths)
    # report_path = generate_multispecies_report(pairwise_results, ancestral_results, base_output)
    

    try:
        print("Starting ancestral reconstruction...")
        ancestral_results = reconstruct_ancestral_states(pairwise_results, busco_paths)
        print("✓ Ancestral reconstruction complete")
    except Exception as e:
        print(f"✗ Ancestral reconstruction failed: {e}")


    try:
        print("Starting report generation...")
        report_path = generate_multispecies_report(pairwise_results, ancestral_results, base_output)
        print("✓ Report complete")
    except Exception as e:
        print(f"✗ Report failed: {e}")
    
    return {
        'pairwise_results': pairwise_results,
        'ancestral_results': ancestral_results,
        'report_path': report_path,
        'output_dir': base_output
    }




def reconstruct_ancestral_states(pairwise_results: Dict, busco_paths: List[str]) -> Dict:
    all_species = [Path(p).stem for p in busco_paths]
    
    all_busco_ids = set()
    for result in pairwise_results.values():
        all_busco_ids.update(result['joined_df']['busco_id'])
    
    ancestral_orientations = {}
    confidence_scores = {}
    orientation_conflicts = {}  
    
    for busco_id in all_busco_ids:
        species_orientations = {} 
        
        for (path1, path2), result in pairwise_results.items():
            species1 = Path(path1).stem
            species2 = Path(path2).stem
            
            gene_data = result['joined_df'][result['joined_df']['busco_id'] == busco_id]
            if not gene_data.empty:
                row = gene_data.iloc[0]
                
                if species1 not in species_orientations:
                    species_orientations[species1] = []
                if species2 not in species_orientations:
                    species_orientations[species2] = []
                    
                species_orientations[species1].append(row['strand1'])
                species_orientations[species2].append(row['strand2'])
        
        final_orientations = {}
        conflicts = []
        
        for species, orientation_list in species_orientations.items():
            unique_orientations = set(orientation_list)
            
            if len(unique_orientations) > 1:

                conflicts.append(f"{species}: {orientation_list}")
                # Use majority rule
                from collections import Counter
                counts = Counter(orientation_list)
                most_common = counts.most_common(1)[0][0]
                final_orientations[species] = most_common
            else:

                final_orientations[species] = orientation_list[0]
        

        if conflicts:
            orientation_conflicts[busco_id] = conflicts
        
        if len(final_orientations) >= 2:

            strand_counts = pd.Series(list(final_orientations.values())).value_counts()
            most_common_strand = strand_counts.index[0]
            confidence = strand_counts.iloc[0] / len(final_orientations)
            
            ancestral_orientations[busco_id] = most_common_strand
            confidence_scores[busco_id] = confidence
    

    inversion_events = []
    for busco_id in ancestral_orientations:
        if busco_id in ancestral_orientations:
            ancestral_strand = ancestral_orientations[busco_id]
            
            for species in all_species:

                species_orientation = None
                for (path1, path2), result in pairwise_results.items():
                    species1 = Path(path1).stem
                    species2 = Path(path2).stem
                    
                    if species in [species1, species2]:
                        gene_data = result['joined_df'][result['joined_df']['busco_id'] == busco_id]
                        if not gene_data.empty:
                            row = gene_data.iloc[0]
                            if species == species1:
                                species_orientation = row['strand1']
                            else:
                                species_orientation = row['strand2']
                            break  
                
                if species_orientation and species_orientation != ancestral_strand:

                    chromosome = None
                    for (path1, path2), result in pairwise_results.items():
                        species1 = Path(path1).stem
                        species2 = Path(path2).stem
                        
                        if species in [species1, species2]:
                            gene_data = result['joined_df'][result['joined_df']['busco_id'] == busco_id]
                            if not gene_data.empty:
                                row = gene_data.iloc[0]
                                if species == species1:
                                    chromosome = row['chr1']
                                else:
                                    chromosome = row['chr2']
                                break
                    
                    inversion_events.append({
                        'busco_id': busco_id,
                        'ancestral_strand': ancestral_strand,
                        'inverted_species': species,
                        'chromosome': chromosome,
                        'confidence': confidence_scores[busco_id]
                    })
    
    return {
        'ancestral_orientations': ancestral_orientations,
        'confidence_scores': confidence_scores,
        'inversion_events': inversion_events,
        'species_list': all_species,
        'orientation_conflicts': orientation_conflicts 
    }


def generate_multispecies_report(pairwise_results: Dict, ancestral_results: Dict, output_dir: Path) -> Path:
    reports_dir = output_dir / 'reports'
    reports_dir.mkdir(exist_ok=True)
    
    ancestral_df = None
    inversions_df = None
    
    if ancestral_results['ancestral_orientations']:
        ancestral_df = pd.DataFrame([
            {
                'busco_id': busco_id,
                'ancestral_strand': strand,
                'confidence': ancestral_results['confidence_scores'].get(busco_id, 0),
                'confidence_category': 'high' if ancestral_results['confidence_scores'].get(busco_id, 0) > 0.8 
                                     else 'medium' if ancestral_results['confidence_scores'].get(busco_id, 0) > 0.6 
                                     else 'low'
            }
            for busco_id, strand in ancestral_results['ancestral_orientations'].items()
        ])
        ancestral_df = ancestral_df.sort_values('confidence', ascending=False)
        
        ancestral_tsv = reports_dir / 'ancestral_gene_orientations.tsv'
        ancestral_df.to_csv(ancestral_tsv, sep='\t', index=False)
    
    if ancestral_results['inversion_events']:
        inversions_df = pd.DataFrame(ancestral_results['inversion_events'])
        inversions_df = inversions_df.sort_values('confidence', ascending=False)
        
        inversions_tsv = reports_dir / 'lineage_specific_inversions.tsv'
        inversions_df.to_csv(inversions_tsv, sep='\t', index=False)
    
    
    report_path = reports_dir / 'multispecies_analysis_report.md'
    
    with open(report_path, 'w') as f:
        f.write("# GIQ1 - 3+ Species Report\n\n")
        f.write(f"**Generated:** {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        
        species_list = ancestral_results['species_list']
        f.write(f"**Species Analyzed:** {', '.join(species_list)}\n")
        f.write(f"**Number of Species:** {len(species_list)}\n")
        f.write(f"**Pairwise Comparisons:** {len(pairwise_results)}\n\n")
        
        f.write("# Ancestral State Reconstruction\n\n")
        f.write(f"- **Genes with ancestral state prediction:** {len(ancestral_results['ancestral_orientations'])}\n")
        f.write(f"- **High confidence predictions (>0.8):** {sum(1 for c in ancestral_results['confidence_scores'].values() if c > 0.8)}\n")
        f.write(f"- **Medium confidence predictions (0.6-0.8):** {sum(1 for c in ancestral_results['confidence_scores'].values() if 0.6 <= c <= 0.8)}\n")
        f.write(f"- **Low confidence predictions (<0.6):** {sum(1 for c in ancestral_results['confidence_scores'].values() if c < 0.6)}\n")
        f.write(f"- **Lineage-specific inversion events:** {len(ancestral_results['inversion_events'])}\n\n")
        

        if ancestral_df is not None and inversions_df is not None:

            combined_df = ancestral_df.copy()
            
            inversion_info = {}
            for _, inv_row in inversions_df.iterrows():
                busco_id = inv_row['busco_id']
                if busco_id not in inversion_info:
                    inversion_info[busco_id] = []
                inversion_info[busco_id].append({
                    'species': inv_row['inverted_species'],
                    'chromosome': inv_row['chromosome']
                })
            
            combined_df['inverted_species'] = combined_df['busco_id'].apply(
                lambda x: ', '.join([info['species'] for info in inversion_info.get(x, [])])
            )
            combined_df['inverted_chromosomes'] = combined_df['busco_id'].apply(
                lambda x: ', '.join([info['chromosome'] for info in inversion_info.get(x, [])])
            )
            combined_df['has_inversions'] = combined_df['inverted_species'] != ''
        
            
            f.write("# Ancestral Reconstruction Results (Top 20)\n\n")
            f.write("| BUSCO ID ---| AncStr | Conf | Categ | Inv. Sp. | Inv. Chr. |\n")
            f.write("|-------------|--------|------|-------|----------|-----------|\n")
            
            for _, row in combined_df.head(20).iterrows():
                inverted_species = row['inverted_species'] if row['inverted_species'] else 'None'
                inverted_chrs = row['inverted_chromosomes'] if row['inverted_chromosomes'] else 'None'
                
                f.write(f"| {row['busco_id']} | {row['ancestral_strand']} | {row['confidence']:.3f} | "
                    f"{row['confidence_category']} | {inverted_species} | {inverted_chrs} |\n")
            
            f.write(f"\n{len(combined_df)} total genes, "
                f"{len(inversions_df)} with lineage-specific inversions*\n\n")

        f.write("# Pairwise Summary\n\n")
        f.write("| Species 1 | Species 2 | Total Genes | Flipped Genes | Flip Rate |\n")
        f.write("|-----------|-----------|-------------|---------------|----------|\n \n \n")

        processed_pairs = set()
        for (path1, path2), result in pairwise_results.items():
            species1 = Path(path1).stem
            species2 = Path(path2).stem
            
            pair_key = tuple(sorted([species1, species2]))
            if pair_key in processed_pairs:
                continue
            processed_pairs.add(pair_key)
            
            total_genes = len(result['joined_df'])
            flipped_genes = result['joined_df']['is_flipped'].sum()
            flip_rate = flipped_genes / total_genes if total_genes > 0 else 0
            
            f.write(f"| {species1} | {species2} | {total_genes} | {flipped_genes} | {flip_rate:.3f} |\n")
            

            if ancestral_results['confidence_scores']:
                confidences = list(ancestral_results['confidence_scores'].values())
                f.write("# Confidence Score Statistics\n\n")
                f.write(f"- **Mean confidence:** {np.mean(confidences):.3f}\n")
                f.write(f"- **Median confidence:** {np.median(confidences):.3f}\n")
                f.write(f"- **Min confidence:** {min(confidences):.3f}\n")
                f.write(f"- **Max confidence:** {max(confidences):.3f}\n\n")
        
    
    return report_path


def compute_contextual_metrics(inversion_df: pd.DataFrame, 
                             joined_df: pd.DataFrame, 
                             config: Dict) -> Dict:
    
    contextual_metrics = {}
    
    rate_metrics = compute_inversion_rate_per_mb_busco(inversion_df, joined_df)
    contextual_metrics['inversion_rates'] = rate_metrics
    
    gene_density_corr = _analyze_gene_density_correlation(inversion_df, joined_df, config)
    contextual_metrics['gene_density_correlation'] = gene_density_corr

    
    return contextual_metrics


def create_analysis_visualisations(inversion_df: pd.DataFrame,
                                 joined_df: pd.DataFrame,
                                 output_dir: Path,
                                 config: Dict) -> Dict:
    
    plots_dir = output_dir / 'plots'
    plots_dir.mkdir(exist_ok=True)
    

    data_dir = output_dir / 'data'
    data_dir.mkdir(parents=True, exist_ok=True)
    
    visualisation_results = {}
    
    species1_name = config.get('first_species_name', 'Species1')
    species2_name = config.get('second_species_name', 'Species2')
    
    import matplotlib.pyplot as plt
    plt.style.use('default')
    plt.rcParams['figure.dpi'] = config.get('dpi', 300)
    plt.rcParams['font.size'] = config.get('font_size', 12)
    
    
    try:
        try:
            create_busco_dotplot(joined_df, plots_dir, config)
            visualisation_results['dotplot'] = plots_dir / 'busco_dotplot.png'
        except Exception as e:
            logger.warning(f"create_busco_dotplot failed: {e}")
            
        try:
           
            fig, ax = create_linearised_dotplot(
                joined_df=joined_df, 
                plots_dir=plots_dir,
                config=config['dotplot_config']  # Pass the nested dotplot config
            )
            visualisation_results['linearised_dotplot'] = plots_dir / 'linearised_busco_dotplot.png'
            
        except Exception as e:
            logger.warning(f"linearised dotplot failed: {e}")
        
    except Exception as e:
        logger.error(f"visualisation creation failed: {e}")
    
    return visualisation_results

def save_results(inversion_df: pd.DataFrame,
                joined_df: pd.DataFrame,
                contextual_metrics: Dict,
                output_dir: Path,
                config: Dict):
    
    data_dir = output_dir / 'data'
    data_dir.mkdir(exist_ok=True)
    

    
    import json
    with open(data_dir / 'contextual_metrics.json', 'w') as f:
        json.dump(contextual_metrics, f, indent=2, default=str)
    



def generate_analysis_report(inversion_df: pd.DataFrame,
                         joined_df: pd.DataFrame,
                         contextual_metrics: Dict,
                         output_dir: Path,
                         config: Dict) -> Path:
  
  reports_dir = output_dir / 'reports'
  reports_dir.mkdir(exist_ok=True)
  
  report_path = reports_dir / 'analysis_report.md'
  
  species1_name = config.get('first_species_name', 'Species1')
  species2_name = config.get('second_species_name', 'Species2')
  all_inversion_events = config.get('all_inversion_events', [])
  
  with open(report_path, 'w') as f:
      f.write("# GIQ1 - Genome Inversion Quantification Report\n\n")
      f.write(f"**Generated:** {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
      f.write(f"**Species Comparison:** {species1_name} vs {species2_name}\n\n")
      
      f.write("# Analysis Summary\n\n")
      f.write(f"- **Total chromosome pairs analyzed:** {len(inversion_df)}\n")
      f.write(f"- **Pairs with inversions:** {len(inversion_df[inversion_df['flipped_genes'] > 0])}\n")
      f.write(f"- **Total genes compared:** {len(joined_df)}\n")
      f.write(f"- **Total flipped genes:** {joined_df['is_flipped'].sum()}\n")
      f.write(f"- **Overall flip rate:** {joined_df['is_flipped'].mean():.1%}\n\n")
      
      f.write("# Inversion Event Analysis\n\n")
      single_gene_events = [e for e in all_inversion_events if e['type'] == 'single-gene']
      multi_gene_events = [e for e in all_inversion_events if e['type'] == 'multi-gene']
      largest_event = max([e['gene_count'] for e in all_inversion_events]) if all_inversion_events else 0
      
      f.write(f"- **Total inversion events detected:** {len(all_inversion_events)}\n")
      f.write(f"- **Single-gene inversions:** {len(single_gene_events)}\n")
      f.write(f"- **Multi-gene inversions:** {len(multi_gene_events)}\n")
      f.write(f"- **Largest inversion:** {largest_event} genes\n\n")
      
      f.write("# Inversion Types Distribution\n\n")
      inversion_types = inversion_df['inversion_type'].value_counts()
      for inv_type, count in inversion_types.items():
          f.write(f"- **{inv_type}:** {count} pairs\n")
      f.write("\n")
      
      f.write("# Contextual Metrics\n\n")
      if 'inversion_rates' in contextual_metrics:
          rates = contextual_metrics['inversion_rates']
          f.write(f"- **Genome coverage (Mb):** {rates.get('genome_coverage_mb', 0):.2f}\n")
          f.write(f"- **Inversions per Mb:** {rates.get('inversion_events_per_mb', 0):.3f}\n")
          f.write(f"- **Flipped genes per Mb:** {rates.get('flipped_genes_per_mb', 0):.3f}\n")
      f.write("\n")
      
      f.write("# Chromosome Pair Details\n\n")
      f.write("| Chr1 | Chr2 | Total Genes | Flipped | Rate | Single | Multi | Largest | Type |\n")
      f.write("|------|------|-------------|---------|------|--------|-------|---------|------|\n")
      
      for _, row in inversion_df.iterrows():
          f.write(f"| {row['chr1']} | {row['chr2']} | {row['total_genes']} | "
                 f"{row['flipped_genes']} | {row['flip_rate']:.2f} | "
                 f"{row['single_gene_inversions']} | {row['multi_gene_inversions']} | "
                 f"{row['largest_inversion']} | {row['inversion_type']} |\n")
      
      if len(all_inversion_events) > 0:
          f.write("\n# Detailed Inversion Events\n\n")
          f.write(f"| -------------Chr Pair----------- | -------Type-------- | -Genes- | --------Start Gene------ | --------End Gene-------- | ------{species1_name} Coordinates--- | -----{species2_name} Coordinates------ | -----Span (bp)---------- |\n")
          f.write("|-----------------------------------|---------------------|---------|--------------------------|--------------------------|--------------------------------------|----------------------------------------|--------------------------|\n")
          
          for event in sorted(all_inversion_events, key=lambda x: x['gene_count'], reverse=True)[:50]:
              coord1 = f"{event['start_pos1']:,} - {event['end_pos1']:,}"
              coord2 = f"{event['start_pos2']:,} - {event['end_pos2']:,}"
              f.write(f"| {event['chr1']}-{event['chr2']} | {event['type']} | {event['gene_count']} | "
                     f" {event['start_gene']} | {event['end_gene']} | {coord1} | {coord2} | {event['span_bp']:,} |\n")
          
          if len(all_inversion_events) > 50:
              f.write(f"\n*Showing top 50 events. Total: {len(all_inversion_events)} events.*\n")
      
      f.write(f"\n---\n*Report generated by GIQ1 v1.0.0*\n")
  
  logger.info(f"Analysis report saved to {report_path}")
  
  return report_path


def main():
    try:
        busco_paths, fasta_paths = parse_arguments()
    except:
        # Fallback to hardcoded paths for testing
        CONFIG.update({
        'first_busco_path': '/Users/zionayokunnu/Documents/Bibionidae/busco-tables/Dioctria_linearis.tsv',
        'second_busco_path': '/Users/zionayokunnu/Documents/Bibionidae/busco-tables/Dioctria_rufipes.tsv',
        'first_species_name': 'Dioctria_linearis',
        'second_species_name': 'Dioctria_rufipes',
        'first_fasta_path' : '/Users/zionayokunnu/Documents/Bibionidae/fasta/Dioctria_linearis.fna',
        'second_fasta_path' : '/Users/zionayokunnu/Documents/Bibionidae/fasta/Dioctria_rufipes.fna'
        })
        
        import random
        random.seed(42)
        np.random.seed(42)
        
        try:
            results = run_busco_inversion_analysis()
            
            print(f"Species analyzed: {' vs '.join(results['species_names'])}")
            print(f"Total inversions detected: {results['inversion_df']['flipped_genes'].sum()}")
            print(f"Output directory: {results['output_dir']}")
            print(f"Report: {results['report_path']}")
            print("=" * 80)
            
            return 0
            
        except Exception as e:
            logger.error(f"Analysis failed: {str(e)}")
            print(f"\nERROR: {str(e)}")
            
            if CONFIG.get('enable_debug_output', False):
                import traceback
                traceback.print_exc()
            
            return 1
    
    import random
    random.seed(42)
    np.random.seed(42)
    
    try:
        if len(busco_paths) == 2:
            CONFIG.update({
                'first_busco_path': busco_paths[0],
                'second_busco_path': busco_paths[1],
                'first_species_name': Path(busco_paths[0]).stem,
                'second_species_name': Path(busco_paths[1]).stem,
            })
            
            if fasta_paths and len(fasta_paths) >= 2:
                CONFIG.update({
                    'first_fasta_path': fasta_paths[0],
                    'second_fasta_path': fasta_paths[1]
                })
            
            results = run_busco_inversion_analysis()
            print(f"Species analyzed: {' vs '.join(results['species_names'])}")
            print(f"Total inversions detected: {results['inversion_df']['flipped_genes'].sum()}")
            print(f"Output directory: {results['output_dir']}")
            print(f"Report: {results['report_path']}")
            
        elif len(busco_paths) >= 3:
            results = run_multispecies_analysis(busco_paths, fasta_paths)
            species_names = [Path(p).stem for p in busco_paths]
            print(f"Multi-species analysis: {', '.join(species_names)}")
            print(f"Pairwise comparisons: {len(results['pairwise_results'])}")
            print(f"Output directory: {results['output_dir']}")
            print(f"Report: {results['report_path']}")
            
        else:
            raise ValueError("Need at least 2 BUSCO files")
            
        print("=" * 80)
        return 0
        
    except Exception as e:
        logger.error(f"Analysis failed: {str(e)}")
        print(f"\nERROR: {str(e)}")
        return 1


if __name__ == "__main__":
    sys.exit(main())