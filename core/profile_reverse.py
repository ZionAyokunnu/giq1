import math
from config.settings import CONFIG
from .query_reference import extract_gene_distribution
from .reverse import detect_extended, apply_adjacency_inversion, apply_flip_inversion, detect_adjacency_inversions, extract_movement_sequence, iterative_detection

def get_permutable_positions(gene_distribution, min_probability=None, max_positions=None):
    """
    Get permutable positions for a gene based on probability thresholds.
    
    Args:
        gene_distribution: {bin_id: normalized_probability}
        min_probability: Minimum probability threshold (uses CONFIG if None)
        max_positions: Maximum positions to return (uses CONFIG if None)
        
    Returns:
        list: [(bin_number, probability)] sorted by probability (descending)
    """
    if min_probability is None:
        min_probability = CONFIG['permutable_positions_threshold']
    if max_positions is None:
        max_positions = CONFIG['max_permutable_positions']
    
    permutable_positions = []
    
    for bin_id, probability in gene_distribution.items():
        if probability >= min_probability:
            try:
                bin_number = int(bin_id.split('_bin_')[1])
                permutable_positions.append((bin_number, probability))
            except (IndexError, ValueError):
                continue
    
    permutable_positions.sort(key=lambda x: (-x[1], x[0]))
    
    return permutable_positions[:max_positions]


def calculate_position_probability(gene_id, target_position, markov_profile):
    """
    Calculate probability of a gene being at a specific position.
    
    Args:
        gene_id: Gene identifier
        target_position: Target bin number
        markov_profile: Markov profile
        
    Returns:
        float: Probability of gene being at target position
    """
    gene_distribution = extract_gene_distribution(markov_profile, gene_id)
    
    for bin_id, probability in gene_distribution.items():
        try:
            bin_number = int(bin_id.split('_bin_')[1])
            if bin_number == target_position:
                return probability
        except (IndexError, ValueError):
            continue
    
    return 0.0 


def evaluate_inversion_step_probability(inversion_step, markov_profile, threshold=None):
    """
    Evaluate if an inversion step meets probability thresholds.
    
    Args:
        inversion_step: Inversion event record with genes and positions
        markov_profile: Markov profile
        threshold: Probability threshold (uses CONFIG if None)
        
    Returns:
        dict: Evaluation results
    """
    if threshold is None:
        threshold = CONFIG['permutable_positions_threshold']
    
    gene_probabilities = {}
    passed_genes = []
    failed_genes = []
    
    genes = inversion_step['genes']
    positions = inversion_step['positions']
    
    for i, gene_id in enumerate(genes):
        final_position = positions[i]
        probability = calculate_position_probability(gene_id, final_position, markov_profile)
        
        gene_probabilities[gene_id] = {
            'position': final_position,
            'probability': probability,
            'passed': probability >= threshold
        }
        
        if probability >= threshold:
            passed_genes.append(gene_id)
        else:
            failed_genes.append(gene_id)
    
    overall_probability = 1.0
    for prob_data in gene_probabilities.values():
        overall_probability *= prob_data['probability']
    
    bit_score = -math.log2(overall_probability) if overall_probability > 0 else float('inf')
    
    return {
        'gene_probabilities': gene_probabilities,
        'passed_genes': passed_genes,
        'failed_genes': failed_genes,
        'overall_probability': overall_probability,
        'bit_score': bit_score,
        'threshold_passed': len(failed_genes) == 0
    }


def generate_gene_specific_pathways(target_gene, current_position, markov_profile, other_genes_positions, initial_sequence):
    """
    Generate alternative pathways for moving a specific gene to its permutable positions.
    
    Args:
        target_gene: Gene that needs alternative pathway
        current_position: Current position of the target gene
        markov_profile: Markov profile
        other_genes_positions: {gene_id: current_position} for other genes
        initial_sequence: Initial movement sequence [(gene_id, position, movement)]
        
    Returns:
        list: Alternative pathways ranked by combined bit score
    """
    # Get permutable positions for target gene
    gene_distribution = extract_gene_distribution(markov_profile, target_gene)
    permutable_positions = get_permutable_positions(gene_distribution)
    
    alternative_pathways = []
    
    for target_position, target_probability in permutable_positions:
        if target_position == current_position:
            continue  # Already at this position
        
        # Generate inversion steps to move target_gene from current to target position
        pathway_steps = generate_inversion_steps_for_gene(
            target_gene, current_position, target_position, other_genes_positions
        )
        
        # Evaluate pathway including target gene and all affected genes
        pathway_evaluation = evaluate_pathway_steps(pathway_steps, markov_profile, initial_sequence)
        
        alternative_pathways.append({
            'target_gene': target_gene,
            'target_position': target_position,
            'target_probability': target_probability,
            'pathway_steps': pathway_steps,
            'evaluation': pathway_evaluation
        })
    
    # Sort by combined bit score (lower is better)
    alternative_pathways.sort(key=lambda x: x['evaluation']['combined_bit_score'])
    
    return alternative_pathways


def generate_inversion_steps_for_gene(gene_id, current_pos, target_pos, other_genes_positions):
    """
    Generate inversion steps to move a specific gene from current to target position.
    Uses the full iterative detection algorithm with gene-specific constraint.
    
    Args:
        gene_id: Gene to move
        current_pos: Current position
        target_pos: Target position
        other_genes_positions: {gene_id: position} for other genes
        
    Returns:
        list: Inversion steps that involve the target gene
    """
    movement_sequence = []
    
    target_movement = target_pos - current_pos
    movement_sequence.append((gene_id, current_pos, target_movement))
    
    for other_gene, pos in other_genes_positions.items():
        movement_sequence.append((other_gene, pos, 0.0)) 
    
    movement_sequence.sort(key=lambda x: x[1])
    
    gene_specific_events = iterative_detection_gene_specific(movement_sequence, gene_id, target_pos)
    
    return gene_specific_events


def iterative_detection_gene_specific(movement_sequence, target_gene, target_position, max_iterations=100):
    """
    Run iterative detection but only apply inversions that involve the target gene.
    
    Args:
        movement_sequence: [(gene_id, position, movement)]
        target_gene: Gene that must be involved in all inversions
        target_position: Target position for the gene
        max_iterations: Maximum iterations
        
    Returns:
        list: Inversion events that involve target gene
    """
    current_sequence = movement_sequence.copy()
    gene_specific_events = []
    iteration = 0
    
    while iteration < max_iterations:
        iteration += 1
        applied_inversion = False
        
        flip_patterns = detect_extended(current_sequence)
        if flip_patterns:
  
            for start_idx, end_idx, flip_size in flip_patterns:
                segment_genes = [current_sequence[i][0] for i in range(start_idx, end_idx + 1)]
                
                if target_gene in segment_genes:
                    current_sequence, inversion_record = apply_flip_inversion(
                        current_sequence, start_idx, end_idx, flip_size
                    )
                    inversion_record['iteration'] = iteration
                    inversion_record['gene_specific_target'] = target_gene
                    gene_specific_events.append(inversion_record)
                    applied_inversion = True
                    break
        

        if not applied_inversion:
            adjacency_inversions = detect_adjacency_inversions(current_sequence)
            if adjacency_inversions:

                for index1, index2 in adjacency_inversions:
                    gene1 = current_sequence[index1][0]
                    gene2 = current_sequence[index2][0]
                    
                    if target_gene == gene1 or target_gene == gene2:
                        current_sequence, inversion_record = apply_adjacency_inversion(
                            current_sequence, index1, index2
                        )
                        inversion_record['iteration'] = iteration
                        inversion_record['gene_specific_target'] = target_gene
                        gene_specific_events.append(inversion_record)
                        applied_inversion = True
                        break
        

        target_gene_current_pos = None
        for gene_name, pos, movement in current_sequence:
            if gene_name == target_gene:
                target_gene_current_pos = pos
                break
        
        if target_gene_current_pos == target_position:
            break 
        

        if not applied_inversion:
            break
    
    return gene_specific_events


def evaluate_pathway_steps(pathway_steps, markov_profile, initial_sequence):
    """
    Evaluate all steps in a pathway, including target gene and all affected genes.
    
    Args:
        pathway_steps: List of inversion steps
        markov_profile: Markov profile
        initial_sequence: Initial movement sequence [(gene_id, position, movement)]
        
    Returns:
        dict: Pathway evaluation
    """
    step_evaluations = []
    total_bit_score = 0.0
    all_affected_genes = set()
    
    for step in pathway_steps:
        step_eval = evaluate_inversion_step_probability(step, markov_profile)
        step_evaluations.append(step_eval)
        total_bit_score += step_eval['bit_score']
        
        all_affected_genes.update(step['genes'])

    final_positions = track_all_gene_final_positions(pathway_steps, initial_sequence)
    

    combined_probability = 1.0
    affected_gene_probabilities = {}
    
    for gene_id in all_affected_genes:
        final_position = final_positions.get(gene_id)
        if final_position is not None:
            gene_probability = calculate_position_probability(gene_id, final_position, markov_profile)
            
            affected_gene_probabilities[gene_id] = {
                'final_position': final_position,
                'probability': gene_probability
            }
            
            combined_probability *= gene_probability
    
    combined_bit_score = -math.log2(combined_probability) if combined_probability > 0 else float('inf')
    
    return {
        'step_evaluations': step_evaluations,
        'total_bit_score': total_bit_score,
        'combined_bit_score': combined_bit_score,
        'affected_genes': list(all_affected_genes),
        'affected_gene_probabilities': affected_gene_probabilities,
        'final_positions': final_positions,
        'pathway_valid': combined_probability > 0
    }


def get_gene_final_position(gene_id, pathway_steps, initial_sequence):
    """
    Get the final position of a gene after all pathway steps by applying each step sequentially.
    
    Args:
        gene_id: Gene identifier
        pathway_steps: List of inversion steps
        initial_sequence: Initial movement sequence [(gene_id, position, movement)]
        
    Returns:
        int: Final position of the gene
    """
    current_sequence = initial_sequence.copy()
    
    for step in pathway_steps:
        if step['type'] == 'flip':
            # Find the positions in current sequence
            start_pos = min(step['positions'])
            end_pos = max(step['positions'])
            

            start_idx = None
            end_idx = None
            for i, (gene, pos, movement) in enumerate(current_sequence):
                if pos == start_pos and start_idx is None:
                    start_idx = i
                if pos == end_pos:
                    end_idx = i
            
            if start_idx is not None and end_idx is not None:
                current_sequence, _ = apply_flip_inversion(
                    current_sequence, start_idx, end_idx, step.get('flip_size', 1)
                )
        
        elif step['type'] == 'adjacency':

            genes_in_step = step['genes']
            if len(genes_in_step) >= 2:
                gene1, gene2 = genes_in_step[0], genes_in_step[1]
                
                index1 = None
                index2 = None
                for i, (gene, pos, movement) in enumerate(current_sequence):
                    if gene == gene1:
                        index1 = i
                    elif gene == gene2:
                        index2 = i
                
                if index1 is not None and index2 is not None:
                    current_sequence, _ = apply_adjacency_inversion(
                        current_sequence, index1, index2
                    )
    
    for gene, pos, movement in current_sequence:
        if gene == gene_id:
            return pos
    
    return None 


def track_all_gene_final_positions(pathway_steps, initial_sequence):
    """
    Track final positions of all genes after pathway steps.
    
    Args:
        pathway_steps: List of inversion steps
        initial_sequence: Initial movement sequence [(gene_id, position, movement)]
        
    Returns:
        dict: {gene_id: final_position}
    """
    current_sequence = initial_sequence.copy()
    
    for step in pathway_steps:
        if step['type'] == 'flip':
            start_pos = min(step['positions'])
            end_pos = max(step['positions'])
            
            start_idx = None
            end_idx = None
            for i, (gene, pos, movement) in enumerate(current_sequence):
                if pos == start_pos and start_idx is None:
                    start_idx = i
                if pos == end_pos:
                    end_idx = i
            
            if start_idx is not None and end_idx is not None:
                current_sequence, _ = apply_flip_inversion(
                    current_sequence, start_idx, end_idx, step.get('flip_size', 1)
                )
        
        elif step['type'] == 'adjacency':
            genes_in_step = step['genes']
            if len(genes_in_step) >= 2:
                gene1, gene2 = genes_in_step[0], genes_in_step[1]
                
                index1 = None
                index2 = None
                for i, (gene, pos, movement) in enumerate(current_sequence):
                    if gene == gene1:
                        index1 = i
                    elif gene == gene2:
                        index2 = i
                
                if index1 is not None and index2 is not None:
                    current_sequence, _ = apply_adjacency_inversion(
                        current_sequence, index1, index2
                    )
    
    final_positions = {}
    for gene, pos, movement in current_sequence:
        final_positions[gene] = pos
    
    return final_positions

def probability_weighted_inversion_analysis(movement_results, markov_profile):
    """
    Perform probability-weighted inversion analysis with alternative pathways.
    
    Args:
        movement_results: Output from movement analysis
        markov_profile: Markov profile
        
    Returns:
        dict: Complete probability-weighted analysis
    """
    standard_analysis = check_events_iteration(movement_results)
    
    # Evaluate each inversion step for probability
    evaluated_steps = []
    problematic_genes = []
    
    for event in standard_analysis['inversion_events']:
        step_eval = evaluate_inversion_step_probability(event, markov_profile)
        
        evaluated_steps.append({
            'original_event': event,
            'probability_evaluation': step_eval
        })
        
        if not step_eval['threshold_passed']:
            problematic_genes.extend(step_eval['failed_genes'])
    
    alternative_pathways = {}
    
    movement_sequence = extract_movement_sequence(movement_results)
    current_positions = {gene_id: pos for gene_id, pos, _ in movement_sequence}
    
    for gene_id in set(problematic_genes):
        if gene_id in current_positions:
            other_positions = {g: p for g, p in current_positions.items() if g != gene_id}
            
            alternatives = generate_gene_specific_pathways(
                gene_id, current_positions[gene_id], markov_profile, other_positions, movement_sequence
            )
            
            alternative_pathways[gene_id] = alternatives
            
    total_standard_bit_score = sum(step['probability_evaluation']['bit_score'] 
                                  for step in evaluated_steps)
    
    integrated_analysis = integrate_best_alternatives(
        standard_analysis, 
        alternative_pathways, 
        evaluated_steps
    )
        
    return {
        'standard_analysis': standard_analysis,
        'integrated_analysis': integrated_analysis,
        'evaluated_steps': evaluated_steps,
        'total_standard_bit_score': total_standard_bit_score,
        'problematic_genes': list(set(problematic_genes)),
        'alternative_pathways': alternative_pathways,
        'has_alternatives': len(alternative_pathways) > 0
    }


def integrate_best_alternatives(standard_analysis, alternative_pathways, evaluated_steps):
    """Replace failed inversions with best alternatives"""
    
    final_events = standard_analysis['inversion_events'].copy()
    
    # For each failed event, replace with best alternative
    for i, step in enumerate(evaluated_steps):
        if not step['probability_evaluation']['threshold_passed']:

            failed_genes = step['probability_evaluation']['failed_genes']
            
            # Replace with alternatives for these genes
            for gene_id in failed_genes:
                if gene_id in alternative_pathways and alternative_pathways[gene_id]:
                    best_alternative = alternative_pathways[gene_id][0] 
  
                    final_events[i] = best_alternative['pathway_steps'][0]  # Use best pathway
    
    return {
        'final_sequence': final_events, 
        'inversion_events': final_events,
        'total_events': len(final_events),
        'converged': True 
    }


def check_events_iteration(movement_results):
    """
    Run inversion detection separately for each chromosome
    """
    all_chromosome_results = {}
    combined_events = []
    
    for chromosome, gene_results in movement_results.items():
        chromosome_sequence = extract_chromosome_movement_sequence(gene_results)
        
        chromosome_analysis = iterative_detection(chromosome_sequence)
        
        all_chromosome_results[chromosome] = chromosome_analysis
        combined_events.extend(chromosome_analysis['inversion_events'])
    
    return {
        'inversion_events': combined_events,
        'chromosome_results': all_chromosome_results,
        'total_events': len(combined_events),
        'total_gene_inversions': sum(event['gene_inversions'] for event in combined_events),
        'adjacency_events': sum(1 for event in combined_events if event['type'] == 'adjacency'),
        'flip_events': sum(1 for event in combined_events if event['type'] == 'flip'),
        'converged': all(r['converged'] for r in all_chromosome_results.values()),
        'iterations': max(r.get('iterations', 0) for r in all_chromosome_results.values())
    }

def extract_chromosome_movement_sequence(gene_results):
    """Extract movement sequence for a single chromosome"""
    movement_sequence = []
    
    for gene_id, result in gene_results.items():
        if result['current_ranges']:
            current_position = result['current_ranges'][0][0]
            mean_movement = result['movement_analysis']['mean_movement']
            movement_sequence.append((gene_id, current_position, mean_movement))
    
    movement_sequence.sort(key=lambda x: x[1])
    return movement_sequence
