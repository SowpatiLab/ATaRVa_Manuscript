# main.py
#!/usr/bin/env python3
import random
from collections import defaultdict
from tqdm import tqdm
import sys
import argparse

from interrupted_TR_simulation import (
    mutate_sequence_with_motif_spacing,
    update_tracking_with_substitutions
)
from compound_TR_simulation import generate_complex_repeat_with_tracking

def parse_args():
    parser = argparse.ArgumentParser(description="Complex Tandem Repeat Simulator")
    
    parser.add_argument('-l', '--num-locations', type=int, default=100, 
                       help='Number of repeat locations')
    parser.add_argument('-o', '--out-prefix', type=str, default='', 
                       help='Output file prefix')
    
    # Sequence type flags
    parser.add_argument('--pure', action='store_true', 
                       help='Generate only pure repeats (no fragments, no mutations)')
    parser.add_argument('--point-mutation', action='store_true', 
                       help='Generate repeats with point mutations only')
    parser.add_argument('--fragmented', action='store_true', 
                       help='Generate fragmented repeats (default)')
    parser.add_argument('--fragmented-with-mutations', action='store_true',
                       help='Generate fragmented repeats with point mutations')
    
    # Mutation parameters
    parser.add_argument('--purity', type=float, default=0.90,
                       help='Purity level for point mutations (0.0-1.0)')
    
    # Existing parameters
    parser.add_argument('--min-motif-size', type=int, default=3, 
                       help='Minimum motif length')
    parser.add_argument('--max-motif-size', type=int, default=10, 
                       help='Maximum motif length')
    parser.add_argument('--min-repeat-length', type=int, default=50, 
                       help='Minimum repeat sequence length')
    parser.add_argument('--max-repeat-length', type=int, default=200, 
                       help='Maximum repeat sequence length')
    parser.add_argument('--no-fragment-fraction', type=float, default=0.2, 
                       help='Fraction of sequences with no fragments (0.0-1.0)')
    
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()
    
    return parser.parse_args()

def main():
    args = parse_args()
    
    if args.out_prefix == '':
        args.out_prefix = "%06x" % random.randint(0, 0xFFFFFFFFFF)
    
    # Determine sequence type
    if args.pure:
        sequence_type = "pure"
    elif args.point_mutation:
        sequence_type = "point_mutation" 
    elif args.fragmented_with_mutations:
        sequence_type = "fragmented_with_mutations"
    else:
        sequence_type = "fragmented"
    
    print(f"File prefix: {args.out_prefix}")
    print(f"Sequence type: {sequence_type}")
    if sequence_type in ["point_mutation", "fragmented_with_mutations"]:
        print(f"Mutation purity: {args.purity}")
    print(f"Motif size range: {args.min_motif_size}-{args.max_motif_size}")
    print(f"Repeat length range: {args.min_repeat_length}-{args.max_repeat_length}")
    if sequence_type in ["fragmented", "fragmented_with_mutations"]:
        print(f"No-fragment fraction: {args.no_fragment_fraction}")
    
    # Load motifs
    rmotifs = defaultdict(list)
    try:
        with open('./filtered-2-HG38_2-10_motifs_d2d.tsv') as fh:
            for line in fh:
                motif, kmer = line.strip().split('\t')
                kmer = int(kmer)
                if args.min_motif_size <= kmer <= args.max_motif_size:
                    rmotifs[kmer].append(motif)
    except FileNotFoundError:
        print("Error: Motif file not found.")
        sys.exit(1)
    
    available_sizes = [k for k in rmotifs.keys() if rmotifs[k]]
    if not available_sizes:
        print(f"Error: No motifs found in size range {args.min_motif_size}-{args.max_motif_size}")
        sys.exit(1)
    
    print(f"Available motif sizes: {sorted(available_sizes)}")
    
    # Open output files
    bed_file = open(f'sim_{args.out_prefix}.bed', 'w')
    fasta_file = open(f'sim_{args.out_prefix}.fa', 'w')
    
    idx_digits = len(str(args.num_locations))
    
    print("Generating complex repeats...")
    for idx in tqdm(range(args.num_locations)):
        repeat_id = f'R{idx:0{idx_digits}d}'
        
        # Choose random motif size and motif
        motif_size = random.choice(available_sizes)
        primary_motif = random.choice(rmotifs[motif_size])
        
        # Determine repeat length
        repeat_length = random.randint(args.min_repeat_length, args.max_repeat_length)
        
        # Handle different sequence types
        mutations_info = []
        
        if sequence_type == "fragmented_with_mutations":
            # Generate fragmented sequence - FIXED UNPACKING
            result = generate_complex_repeat_with_tracking(
                primary_motif, repeat_length, sequence_type, fragment_mutation_rate=0.95
            )
            # Take only first 3 values
            repeat_seq, structure_summary, num_fragments = result[0], result[1], result[2]
            
            # Apply substitutions
            repeat_seq, mutations_info = mutate_sequence_with_motif_spacing(
                repeat_seq, primary_motif, args.purity
            )
            
            # Update tracking to reflect substitutions
            structure_summary = update_tracking_with_substitutions(
                structure_summary, primary_motif, mutations_info, repeat_seq
            )
            
        elif sequence_type == "point_mutation":
            # Generate pure sequence - FIXED UNPACKING
            result = generate_complex_repeat_with_tracking(
                primary_motif, repeat_length, "pure"
            )
            # Take only first 3 values
            repeat_seq, structure_summary, num_fragments = result[0], result[1], result[2]
            
            # Apply substitutions
            repeat_seq, mutations_info = mutate_sequence_with_motif_spacing(
                repeat_seq, primary_motif, args.purity
            )
            
            # Update tracking
            structure_summary = update_tracking_with_substitutions(
                structure_summary, primary_motif, mutations_info, repeat_seq
            )
            
        else:
            # Original behavior - FIXED UNPACKING
            result = generate_complex_repeat_with_tracking(
                primary_motif, repeat_length, sequence_type
            )
            # Take only first 3 values
            repeat_seq, structure_summary, num_fragments = result[0], result[1], result[2]
        
        # Write to FASTA
        print(f'>{repeat_id};{primary_motif};{len(primary_motif)}', file=fasta_file)
        print(repeat_seq, file=fasta_file)
        
        # Write to BED
        bed_columns = [
            repeat_id,
            '0',
            str(len(repeat_seq)),
            str(motif_size),
            primary_motif,
            structure_summary,
            str(num_fragments),
            str(len(repeat_seq))
        ]
        
        if mutations_info:
            mutations_str = ';'.join(['|'.join(mut) for mut in mutations_info])
            bed_columns.append(mutations_str)
        
        print('\t'.join(bed_columns), file=bed_file)
    
    bed_file.close()
    fasta_file.close()
    
    print(f"Simulation completed!")
    print(f"Generated {args.num_locations} {sequence_type} repeat sequences")
    print(f"Output files: sim_{args.out_prefix}.bed, sim_{args.out_prefix}.fa")

if __name__ == "__main__":
    main()