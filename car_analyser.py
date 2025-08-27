
import re
from collections import defaultdict, Counter
import csv
import sys

def parse_car_data(car_file):

    cars = {}
    with open(car_file, 'r') as f:
        for line in f:
            line = line.strip()
            if ':' in line:
                parts = line.split(': ')
                car_id = parts[0]
                gene_ids = [int(x) for x in parts[1].split()]
                cars[car_id] = gene_ids
    return cars

def parse_dioctria_busco(busco_file):

    busco_to_chr = {}
    chr_stats = defaultdict(int)
    
    with open(busco_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('#') or not line:
                continue
            
            parts = line.split('\t')
            if len(parts) >= 3:
                busco_id = parts[0]
                status = parts[1]
                chromosome = parts[2]
                
                if status == 'Complete':
                    busco_to_chr[busco_id] = chromosome
                    chr_stats[chromosome] += 1
    
    return busco_to_chr, dict(chr_stats)

def load_ancgenes_mapping(ancgenes_file):
    """
    Load the ancGenes mapping file that converts numerical indices to BUSCO IDs
    Handles AGORA format: 504100.00.1 species|busco_id species|busco_id
    """
    import bz2
    
    index_to_busco = {}
    try:

        if ancgenes_file.endswith('.bz2'):
            file_opener = lambda f: bz2.open(f, 'rt', encoding='utf-8')
        else:
            file_opener = lambda f: open(f, 'r', encoding='utf-8')
        
        with file_opener(ancgenes_file) as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if line and not line.startswith('#'):
                    parts = line.split()
                    
                    if len(parts) >= 2:
                        try:
                            # Extract index from first column (e.g., "504100.00.1" -> 1)
                            ancestor_index = parts[0]
                            if '.' in ancestor_index:
                                index = int(ancestor_index.split('.')[-1])
                            else:
                                index = int(parts[0])
                            
                            # Extract BUSCO ID from second column (e.g., "species|100019at7147" -> "100019at7147")
                            species_busco = parts[1]
                            if '|' in species_busco:
                                busco_id = species_busco.split('|')[1]
                            else:
                                busco_id = species_busco
                            
                            index_to_busco[index] = busco_id
                            
                        except (ValueError, IndexError) as e:
                            print(f"Warning: Could not parse line {line_num}: {line}")
                            continue
        
        # print(f"Loaded {len(index_to_busco)} index-to-BUSCO mappings")
        if len(index_to_busco) > 0:
            
            sample_items = list(index_to_busco.items())[:3]
            # print(f"Sample mappings: {sample_items}")
        
    except FileNotFoundError:
        print(f"Warning: ancGenes mapping file '{ancgenes_file}' not found.")
        
        return None
    except Exception as e:
        print(f"Error reading ancGenes file: {e}")
        return None
    
    return index_to_busco

def analyse_car_chromosome_alignment(cars, busco_to_chr, index_to_busco=None):

    results = []
    
    for car_id, gene_indices in cars.items():
        car_analysis = {
            'car_id': car_id,
            'total_genes': len(gene_indices),
            'mapped_genes': 0,
            'unmapped_genes': 0,
            'chromosome_counts': Counter(),
            'chromosome_percentages': {},
            'dominant_chromosome': None,
            'dominant_percentage': 0.0,
            'purity_score': 0.0,
            'mapped_busco_ids': [],
            'unmapped_indices': []
        }
        
        # Map indices to BUSCO IDs and then to chromosomes
        for gene_index in gene_indices:
            busco_id = None
            chromosome = None
            
            # Convert index to BUSCO ID
            if index_to_busco and gene_index in index_to_busco:
                busco_id = index_to_busco[gene_index]
                car_analysis['mapped_busco_ids'].append(busco_id)
                
                # Convert BUSCO ID to chromosome
                if busco_id in busco_to_chr:
                    chromosome = busco_to_chr[busco_id]
                    car_analysis['chromosome_counts'][chromosome] += 1
                    car_analysis['mapped_genes'] += 1
                else:
                    car_analysis['unmapped_genes'] += 1
            else:
                car_analysis['unmapped_indices'].append(gene_index)
                car_analysis['unmapped_genes'] += 1
        

        if car_analysis['mapped_genes'] > 0:
            total_mapped = car_analysis['mapped_genes']
            
            for chromosome, count in car_analysis['chromosome_counts'].items():
                percentage = (count / total_mapped) * 100
                car_analysis['chromosome_percentages'][chromosome] = percentage
            

            if car_analysis['chromosome_counts']:
                dominant_chr = car_analysis['chromosome_counts'].most_common(1)[0]
                car_analysis['dominant_chromosome'] = dominant_chr[0]
                car_analysis['dominant_percentage'] = (dominant_chr[1] / total_mapped) * 100
                car_analysis['purity_score'] = car_analysis['dominant_percentage'] / 100
        
        results.append(car_analysis)
    
    return results

def calculate_global_metrics(results):

    total_cars = len(results)
    high_purity_cars = sum(1 for r in results if r['purity_score'] >= 0.8)
    medium_purity_cars = sum(1 for r in results if 0.5 <= r['purity_score'] < 0.8)
    low_purity_cars = sum(1 for r in results if r['purity_score'] < 0.5)
    
    all_chromosomes = set()
    for result in results:
        all_chromosomes.update(result['chromosome_counts'].keys())
    
    avg_purity = sum(r['purity_score'] for r in results) / total_cars if total_cars > 0 else 0
    avg_car_size = sum(r['total_genes'] for r in results) / total_cars if total_cars > 0 else 0
    
    return {
        'total_cars': total_cars,
        'high_purity_cars': high_purity_cars,
        'medium_purity_cars': medium_purity_cars,
        'low_purity_cars': low_purity_cars,
        'average_purity': avg_purity,
        'average_car_size': avg_car_size,
        'chromosomes_found': sorted(all_chromosomes)
    }

def write_results(results, global_metrics, output_file):

    with open(output_file, 'w', newline='') as csvfile:

        csvfile.write("# CAR to Chromosome Alignment Analysis\n")
        csvfile.write(f"# Total CARs: {global_metrics['total_cars']}\n")
        csvfile.write(f"# High purity (what percentage of genes within that CAR come from the same chromosome?) CARs (≥80%): {global_metrics['high_purity_cars']}\n")
        csvfile.write(f"# Medium purity CARs (50-79%): {global_metrics['medium_purity_cars']}\n")
        csvfile.write(f"# Low purity CARs (<50%): {global_metrics['low_purity_cars']}\n")
        csvfile.write(f"# Average purity: {global_metrics['average_purity']:.2%}\n")
        csvfile.write(f"# Average CAR size: {global_metrics['average_car_size']:.1f} genes\n")
        csvfile.write(f"# Chromosomes found: {', '.join(global_metrics['chromosomes_found'])}\n")
        csvfile.write("#\n")
        
        # Write detailed results
        writer = csv.writer(csvfile)
        writer.writerow([
            'CAR_ID', 'Total_Genes', 'Mapped_Genes', 'Unmapped_Genes',
            'Dominant_Chromosome', 'Dominant_Percentage', 'Purity_Score',
            'All_Chromosome_Percentages'
        ])
        
        for result in results:
            chr_percentages = '; '.join([
                f"{chr}:{pct:.1f}%" 
                for chr, pct in sorted(result['chromosome_percentages'].items())
            ])
            
            writer.writerow([
                result['car_id'],
                result['total_genes'],
                result['mapped_genes'],
                result['unmapped_genes'],
                result['dominant_chromosome'] or 'None',
                f"{result['dominant_percentage']:.1f}%",
                f"{result['purity_score']:.3f}",
                chr_percentages
            ])

def print_summary(results, global_metrics):

    print("\n" + "="*60)
    print("CAR TO CHROMOSOME ALIGNMENT ANALYSIS")
    print("="*60)
    
    print(f"Total CARs analysed: {global_metrics['total_cars']}")
    print(f"Average CAR size: {global_metrics['average_car_size']:.1f} genes")
    print(f"Average purity: {global_metrics['average_purity']:.2%}")
    
    print(f"\nPurity Distribution:")
    print(f"  High purity (≥80%): {global_metrics['high_purity_cars']} CARs")
    print(f"  Medium purity (50-79%): {global_metrics['medium_purity_cars']} CARs")
    print(f"  Low purity (<50%): {global_metrics['low_purity_cars']} CARs")
    
    print(f"\nChromosomes found: {', '.join(global_metrics['chromosomes_found'])}")
    
    print(f"\nTop 10 CARs by purity:")
    sorted_results = sorted(results, key=lambda x: x['purity_score'], reverse=True)
    for i, result in enumerate(sorted_results[:10]):
        print(f"  {result['car_id']}: {result['purity_score']:.2%} "
              f"({result['dominant_chromosome']}, {result['total_genes']} genes)")
    
    print(f"\nLargest CARs:")
    sorted_by_size = sorted(results, key=lambda x: x['total_genes'], reverse=True)
    for i, result in enumerate(sorted_by_size[:5]):
        print(f"  {result['car_id']}: {result['total_genes']} genes, "
              f"{result['purity_score']:.2%} purity ({result['dominant_chromosome']})")




    """
    Extract gene IDs from each CAR with the below:
    
    bzcat agora_results/scaffolded_blocks.504100.00.list.bz2 | awk '{
    printf "CAR_%d: ", NR
    for(i=3; i<=(2+$2); i++) printf "%s ", $i
    print ""
}' > agora_results/car_gene_assignments_scaffolded.txt
    
    
    bzcat agora_results/scaffolded_blocks.504100.00.list.bz2 | awk '{
    printf "CAR_%d: ", NR
    for(i=3; i<=(2+$2); i++) printf "%s ", $i
    print ""
}' > agora_results/car_gene_assignments_scaffolded-chr-aware.txt

python3 car_analyser.py
    
    """
def main():

    car_file = "/Users/zionayokunnu/Documents/Giq/agora_results/car_gene_assignments_scaffolded.txt"
    dioctria_busco_file = "/Users/zionayokunnu/Documents/Bibionidae/busco-tables/Dioctria_rufipes.tsv"
    ancgenes_file = "/Users/zionayokunnu/Documents/Giq/reformatted/ancGenes.504100.00.list.bz2"
    output_file = "/Users/zionayokunnu/Documents/Giq/agora_results/car_chromosome_analysis-anchor6.csv"

    
    print("Loading CAR data...")
    cars = parse_car_data(car_file)
    print(f"Loaded {len(cars)} CARs")
    
    print("Loading Dioctria BUSCO data...")
    busco_to_chr, chr_stats = parse_dioctria_busco(dioctria_busco_file)
    print(f"Loaded {len(busco_to_chr)} BUSCO-to-chromosome mappings")
    print(f"Chromosomes in Dioctria: {list(chr_stats.keys())}")
    
    print("Loading ancGenes mapping...")
    index_to_busco = load_ancgenes_mapping(ancgenes_file)
    if index_to_busco:
        print(f"Loaded {len(index_to_busco)} index-to-BUSCO mappings")
    else:
        print("No ancGenes mapping available - will report unmapped genes")
    
    print("Analyzing CAR-chromosome alignments...")
    results = analyse_car_chromosome_alignment(cars, busco_to_chr, index_to_busco)
    
    print("Calculating global metrics...")
    global_metrics = calculate_global_metrics(results)
    
    print(f"Writing results to {output_file}...")
    write_results(results, global_metrics, output_file)
    
    print_summary(results, global_metrics)
    print(f"\nDetailed results saved to: {output_file}")

if __name__ == "__main__":
    main()