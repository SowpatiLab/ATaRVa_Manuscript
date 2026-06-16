import argparse
import sys
from Bio import pairwise2

def expand_format1(decomp):
    """Expand format 1 decomposition to sequence with X breaks"""
    parts = decomp.split('|')
    expanded = []
    for part in parts:
        if ':' in part:
            seq, count = part.split(':')
            expanded.append(seq * int(count))
        else:
            # Handle single elements without colon
            expanded.append(part)
    return 'X'.join(expanded)  # Keep X as gap markers

def expand_format2(decomp):
    """Expand format 2 decomposition to sequence with X breaks"""
    parts = decomp.replace('(', '').replace(')', '').split('-')
    expanded = []
    for part in parts:
        if part[-1].isdigit():
            for i, char in enumerate(reversed(part)):
                if not char.isdigit():
                    break
            seq_end = len(part) - i
            seq = part[:seq_end]
            count = int(part[seq_end:])
            expanded.append(seq * int(count))
        else:
            expanded.append(part)
    return 'X'.join(expanded)  # Keep X as gap markers

def calculate_alignment_metrics(seq1, seq2, match_score=2, mismatch_penalty=-1, gap_open_penalty=-2, gap_extend_penalty=-0.5):
    """Calculate alignment metrics with gap penalties"""
    
    # Replace X with gaps for alignment (X becomes - in alignment)
    seq1_clean = seq1.replace('X', '-')
    seq2_clean = seq2.replace('X', '-')
    
    # Perform global alignment with affine gap penalties
    alignments = pairwise2.align.globalms(
        seq1_clean, 
        seq2_clean,
        match_score,          # Score for matches
        mismatch_penalty,     # Penalty for mismatches
        gap_open_penalty,     # Penalty for opening a gap
        gap_extend_penalty,   # Penalty for extending a gap
        one_alignment_only=True
    )
    
    if alignments:
        best_alignment = alignments[0]
        raw_score = best_alignment[2]  # Raw alignment score
        
        # Calculate normalized score
        # Use the shorter sequence length for normalization
        seq1_bases = seq1.replace('X', '')
        seq2_bases = seq2.replace('X', '')
        min_length = min(len(seq1_bases), len(seq2_bases))
        
        if min_length > 0:
            # Maximum possible = perfect alignment of the shorter sequence
            max_possible = min_length * match_score
            percent_match = (raw_score / max_possible) * 100
        else:
            percent_match = 0
        
        return raw_score, percent_match, best_alignment
    return 0, 0, None

def read_bed_file(filename):
    """Read BED file and return lines"""
    with open(filename, 'r') as f:
        return [line.strip() for line in f if line.strip()]

def process_bed_files_with_gap_penalty(simulated_bed, decomposed_bed, output_bed, 
                                     match_score=2, mismatch_penalty=-1, 
                                     gap_open_penalty=-2, gap_extend_penalty=-0.5,
                                     debug=False, save_mismatches=False):
    """Process BED files with configurable gap penalties"""
    
    # Read files
    sim_lines_raw = read_bed_file(simulated_bed)
    decomp_lines_raw = read_bed_file(decomposed_bed)
    
    # Parse lines
    sim_lines = [line.split('\t') for line in sim_lines_raw]
    decomp_lines = [line.split('\t') for line in decomp_lines_raw]
    
    # Ensure we have the same number of rows
    min_rows = min(len(sim_lines), len(decomp_lines))
    
    # Statistics
    stats = {
        'total': 0,
        'length_mismatches': 0,
        'base_length_mismatches': 0,
        'errors': 0
    }
    
    # Store mismatched lines
    mismatched_lines = []
    
    # Open output files
    out = open(output_bed, 'w')
    if save_mismatches:
        mismatch_out = open(output_bed.replace('.bed', '_mismatches.bed'), 'w')
        mismatch_out.write("line_number\tchrom\tstart\tend\tsim_decomp\tdecomp_decomp\tsim_length\tdecomp_length\tsim_bases\tdecomp_bases\tdiff_bases\tsim_raw_line\tdecomp_raw_line\n")
    
    # Write header
    out.write("# Alignment parameters: match={}, mismatch={}, gap_open={}, gap_extend={}\n".format(
        match_score, mismatch_penalty, gap_open_penalty, gap_extend_penalty))
    out.write("chrom\tstart\tend\tsim_motif\tdecomp_motif\traw_score\tpercent_match\tsim_decomp\tdecomp_decomp\tsim_seq_len\tdecomp_seq_len\tsim_sequence\tdecomp_sequence\n")
    
    processed = 0
    total_raw_score = 0
    total_percent_match = 0
    
    for i in range(min_rows):
        sim_parts = sim_lines[i]
        decomp_parts = decomp_lines[i]
        
        if len(sim_parts) >= 6 and len(decomp_parts) >= 6:
            chrom, start, end = sim_parts[0], sim_parts[1], sim_parts[2]
            sim_motif = sim_parts[3]
            sim_decomp = sim_parts[5]
            decomp_motif = decomp_parts[3]
            decomp_decomp = decomp_parts[5] if len(decomp_parts) > 5 else decomp_parts[4]
            
            try:
                # Expand sequences (keeping X as gap markers)
                seq1 = expand_format1(sim_decomp)
                seq2 = expand_format2(decomp_decomp)
                
                stats['total'] += 1
                
                # Check length match
                len1 = len(seq1)
                len2 = len(seq2)
                
                # Check base length (without X)
                seq1_bases = seq1.replace('X', '')
                seq2_bases = seq2.replace('X', '')
                base_len1 = len(seq1_bases)
                base_len2 = len(seq2_bases)
                
                mismatch_detected = False
                mismatch_info = None
                
                if base_len1 != base_len2:
                    stats['base_length_mismatches'] += 1
                    mismatch_detected = True
                    
                    # Store mismatch info
                    mismatch_info = {
                        'line_number': i + 1,
                        'chrom': chrom,
                        'start': start,
                        'end': end,
                        'sim_decomp': sim_decomp,
                        'decomp_decomp': decomp_decomp,
                        'sim_length': len1,
                        'decomp_length': len2,
                        'sim_bases': base_len1,
                        'decomp_bases': base_len2,
                        'diff_bases': abs(base_len1 - base_len2),
                        'sim_raw_line': sim_lines_raw[i],
                        'decomp_raw_line': decomp_lines_raw[i]
                    }
                    
                    mismatched_lines.append(mismatch_info)
                    
                    if debug:
                        print(f"\n{'='*80}")
                        print(f"[DEBUG] BASE LENGTH MISMATCH at line {i+1}:")
                        print(f"{'='*80}")
                        print(f"Chromosome: {chrom}")
                        print(f"Region: {start}-{end}")
                        print(f"\nSimulated decomposition: {sim_decomp}")
                        print(f"Decomposed decomposition: {decomp_decomp}")
                        print(f"\nSimulated expanded length: {len1} chars ({base_len1} bases without X)")
                        print(f"Decomposed expanded length: {len2} chars ({base_len2} bases without X)")
                        print(f"Base length difference: {abs(base_len1 - base_len2)} bases")
                        print(f"\nSimulated raw BED line: {sim_lines_raw[i]}")
                        print(f"Decomposed raw BED line: {decomp_lines_raw[i]}")
                        
                        # Show first 100 chars of each expanded sequence
                        print(f"\nFirst 100 chars of expanded sequences:")
                        print(f"Simulated:  {seq1[:100]}")
                        print(f"Decomposed: {seq2[:100]}")
                        
                        # Show pattern differences
                        print(f"\nPattern analysis:")
                        print(f"Simulated pattern counts: {sim_decomp.count(':') + sim_decomp.count('|') + 1} segments")
                        print(f"Decomposed pattern counts: {decomp_decomp.count(')')} explicit segments")
                        print(f"{'='*80}")
                
                if len1 != len2 and base_len1 == base_len2:
                    stats['length_mismatches'] += 1
                    if debug:
                        print(f"\n[DEBUG] Line {i+1}: X-count mismatch only")
                        print(f"  Same base length: {base_len1}")
                        print(f"  Different total length due to X separators")
                        print(f"  Simulated X count: {seq1.count('X')}")
                        print(f"  Decomposed X count: {seq2.count('X')}")
                
                # Calculate alignment with gap penalties
                raw_score, percent_match, alignment = calculate_alignment_metrics(
                    seq1, seq2, match_score, mismatch_penalty, gap_open_penalty, gap_extend_penalty
                )
                
                # Write output
                out.write(f"{chrom}\t{start}\t{end}\t{sim_motif}\t{decomp_motif}\t")
                out.write(f"{raw_score:.2f}\t{percent_match:.2f}\t")
                out.write(f"{sim_decomp}\t{decomp_decomp}\t{len1}\t{len2}\t{seq1}\t{seq2}\n")
                
                processed += 1
                total_raw_score += raw_score
                total_percent_match += percent_match
                
            except Exception as e:
                stats['errors'] += 1
                if debug:
                    print(f"\n[ERROR] Processing line {i+1}: {e}")
                    print(f"  Simulated: {sim_parts}")
                    print(f"  Decomposed: {decomp_parts}")
                continue
    
    # Close output files
    out.close()
    
    # Save mismatches to file
    if save_mismatches and mismatched_lines:
        for mismatch in mismatched_lines:
            mismatch_out.write(f"{mismatch['line_number']}\t{mismatch['chrom']}\t{mismatch['start']}\t{mismatch['end']}\t")
            mismatch_out.write(f"{mismatch['sim_decomp']}\t{mismatch['decomp_decomp']}\t")
            mismatch_out.write(f"{mismatch['sim_length']}\t{mismatch['decomp_length']}\t")
            mismatch_out.write(f"{mismatch['sim_bases']}\t{mismatch['decomp_bases']}\t{mismatch['diff_bases']}\t")
            mismatch_out.write(f"{mismatch['sim_raw_line']}\t{mismatch['decomp_raw_line']}\n")
        mismatch_out.close()
        print(f"\nSaved {len(mismatched_lines)} mismatched lines to: {output_bed.replace('.bed', '_mismatches.bed')}")
    
    # Print summary
    print(f"\n{'='*80}")
    print(f"PROCESSING SUMMARY")
    print(f"{'='*80}")
    print(f"Total sequences: {stats['total']}")
    print(f"Successfully processed: {processed}")
    print(f"Base length mismatches (without X): {stats['base_length_mismatches']} ({stats['base_length_mismatches']/stats['total']*100:.1f}%)")
    print(f"X-count mismatches only: {stats['length_mismatches']} ({stats['length_mismatches']/stats['total']*100:.1f}%)")
    print(f"Errors: {stats['errors']}")
    
    if processed > 0:
        print(f"\nALIGNMENT RESULTS")
        print(f"{'='*80}")
        print(f"Average raw alignment score: {total_raw_score/processed:.2f}")
        print(f"Average percent match: {total_percent_match/processed:.2f}%")
        print(f"Scoring parameters: match={match_score}, mismatch={mismatch_penalty}, gap_open={gap_open_penalty}, gap_extend={gap_extend_penalty}")
    
    # Print top mismatches if any
    if mismatched_lines:
        print(f"\n{'='*80}")
        print(f"TOP 10 BASE LENGTH MISMATCHES (by difference):")
        print(f"{'='*80}")
        
        # Sort by difference (largest first)
        sorted_mismatches = sorted(mismatched_lines, key=lambda x: x['diff_bases'], reverse=True)
        
        for j, mismatch in enumerate(sorted_mismatches[:10]):
            print(f"\n#{j+1} Line {mismatch['line_number']}: {mismatch['chrom']}:{mismatch['start']}-{mismatch['end']}")
            print(f"  Simulated: {mismatch['sim_bases']} bases")
            print(f"  Decomposed: {mismatch['decomp_bases']} bases")
            print(f"  Difference: {mismatch['diff_bases']} bases")
            print(f"  Simulated pattern: {mismatch['sim_decomp']}")
            print(f"  Decomposed pattern: {mismatch['decomp_decomp']}")
    
    if not processed:
        print("No sequences processed successfully")

def main():
    parser = argparse.ArgumentParser(description='Compare BED files with sequence decomposition and calculate alignment scores with gap penalties')
    parser.add_argument('simulated_bed', help='Input BED file from simulated data (format 1)')
    parser.add_argument('decomposed_bed', help='Input BED file from decomposition script (format 2)')
    parser.add_argument('output_bed', help='Output BED file with alignment scores')
    parser.add_argument('--match', type=float, default=1.0, help='Match score (default: 1.0)')
    parser.add_argument('--mismatch', type=float, default=-1.0, help='Mismatch penalty (default: -1.0)')
    parser.add_argument('--gap_open', type=float, default=0.0, help='Gap opening penalty (default: 0.0)')
    parser.add_argument('--gap_extend', type=float, default=0.0, help='Gap extension penalty (default: 0.0)')
    parser.add_argument('--debug', '-d', action='store_true', help='Enable debug output for length mismatches')
    parser.add_argument('--save_mismatches', '-s', action='store_true', help='Save mismatched lines to separate file')
    parser.add_argument('--verbose', '-v', action='store_true', help='Verbose output')
    
    args = parser.parse_args()
    
    if args.verbose:
        print(f"Processing files:")
        print(f"  Simulated data: {args.simulated_bed}")
        print(f"  Decomposed data: {args.decomposed_bed}")
        print(f"  Output: {args.output_bed}")
        print(f"  Scoring: match={args.match}, mismatch={args.mismatch}, gap_open={args.gap_open}, gap_extend={args.gap_extend}")
        if args.debug:
            print(f"  Debug mode: ON")
        if args.save_mismatches:
            print(f"  Save mismatches: ON")
    
    try:
        process_bed_files_with_gap_penalty(
            args.simulated_bed, 
            args.decomposed_bed, 
            args.output_bed,
            args.match,
            args.mismatch,
            args.gap_open,
            args.gap_extend,
            args.debug,
            args.save_mismatches
        )
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
