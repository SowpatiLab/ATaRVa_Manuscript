from collections import Counter
from bitarray import bitarray
import regex as re
import argparse

divisor_dict = {2:[1], 3:[1], 4:[1,2], 5:[1], 6:[1,2,3], 7:[1], 8:[1,2,4], 9:[1,3], 10:[1,2,5]}

# Dictionary to store first appearance of cyclic motifs for the current sequence
cyclic_motif_registry = {}

def get_cyclic_variants(motif):
    """Get all cyclic rotations of a motif"""
    n = len(motif)
    return [motif[i:] + motif[:i] for i in range(n)]

def get_canonical_motif(motif):
    """Get the canonical (lexicographically smallest) cyclic variant"""
    return min(get_cyclic_variants(motif))

def register_and_get_motif(motif):
    """Register a motif and return the first variant for its cyclic class"""
    if not motif or len(motif) == 0:
        return motif
    
    canonical = get_canonical_motif(motif)
    
    if canonical not in cyclic_motif_registry:
        # First time seeing this cyclic class, register the actual motif
        cyclic_motif_registry[canonical] = motif
        return motif
    else:
        # Return the first registered variant for this cyclic class
        return cyclic_motif_registry[canonical]

def convert_to_bitset(seq):
    lbit = {'A': '0', 'C': '0', 'G': '1', 'T': '1', 'N': '1'}
    rbit = {'A': '0', 'C': '1', 'G': '0', 'T': '1', 'N': '1'}
    
    lbitseq = bitarray()
    rbitseq = bitarray()
    
    for s in seq:
        lbitseq.extend(lbit.get(s, '1'))
        rbitseq.extend(rbit.get(s, '1'))
    
    return lbitseq, rbitseq

def shift_and_match(seq):
    shift_list = []
    
    for shift in range(1, 11):
        lbitseq, rbitseq = convert_to_bitset(seq)
        lmatch = ~(lbitseq ^ (lbitseq >> shift))
        rmatch = ~(rbitseq ^ (rbitseq >> shift))
        match = lmatch & rmatch
        shift_list.append(match)

    return shift_list

def kmp_search_non_overlapping(text, pattern):
    def compute_lps(pattern):
        lps = [0] * len(pattern)
        length = 0
        i = 1
        while i < len(pattern):
            if pattern[i] == pattern[length]:
                length += 1
                lps[i] = length
                i += 1
            else:
                if length != 0:
                    length = lps[length - 1]
                else:
                    lps[i] = 0
                    i += 1
        return lps

    lps = compute_lps(pattern)
    result = []
    i = 0
    j = 0  

    while i < len(text):
        if pattern[j] == text[i]:
            i += 1
            j += 1

        if j == len(pattern):  
            result.append(i - j)
            j = 0  
        elif i < len(text) and pattern[j] != text[i]:
            if j != 0:
                j = lps[j - 1]
            else:
                i += 1

    return result

def get_most_frequent_motif(sequence, motif_size, motif):
    """Get the most frequent motif, using the first-appearing cyclic variant"""
    if len(sequence) < motif_size:  
        return sequence, None
    
    # Get all motifs of the given size
    repeating_units = [sequence[i:i + motif_size] for i in range(len(sequence) - motif_size + 1)]
    motif_counts = Counter(repeating_units)
    
    if not motif_counts:  
        return sequence, None
    
    # If a specific motif is provided in the input, use it as the reference
    if motif:
        # Register the provided motif first
        registered_motif = register_and_get_motif(motif)
        
        # Check if this motif or its cyclic variants exist in the sequence
        for m, count in motif_counts.items():
            if get_canonical_motif(m) == get_canonical_motif(registered_motif):
                # Use the registered variant (first-appearing cyclic form)
                return registered_motif, None
    
    # Find the most frequent motif, considering cyclic variations
    best_motif = None
    best_count = 0
    
    # Group motifs by their canonical form and sum their counts
    motif_groups = {}
    for m, count in motif_counts.items():
        canonical = get_canonical_motif(m)
        if canonical not in motif_groups:
            motif_groups[canonical] = {'total': 0, 'variants': []}
        motif_groups[canonical]['total'] += count
        motif_groups[canonical]['variants'].append((m, count))
    
    # For each canonical group, find the best candidate
    for canonical, group_info in motif_groups.items():
        total_count = group_info['total']
        
        if total_count > best_count:
            best_count = total_count
            
            # Choose the variant to use
            if canonical in cyclic_motif_registry:
                # Use the already registered variant
                best_motif = cyclic_motif_registry[canonical]
            else:
                # This is the first time seeing this cyclic class
                # Use the most frequent variant, or the first one if tied
                variants = group_info['variants']
                # Sort by frequency, then by first appearance
                max_freq = max(count for _, count in variants)
                top_variants = [m for m, count in variants if count == max_freq]
                
                # Register and use the first top variant
                best_motif = register_and_get_motif(top_variants[0])
    
    return best_motif, None

def max_match(shift_list, gap_regions):
    gap_wise_shift = []
    
    for each_gap in gap_regions:
        max_matches = 0
        best_shift = -1
        for idx, each_shift in enumerate(shift_list):
            if each_shift == 0: 
                continue
            each_shift = each_shift[each_gap[0]: each_gap[1]]
            count = each_shift.count()
            
            if count > max_matches and count > (idx + 1):
                max_matches = count
                best_shift = idx + 1
        gap_wise_shift.append(best_shift)
    return gap_wise_shift

slide_threshold = {1:1, 2:2, 3:3, 4:4, 5:5, 6:6, 7:7, 8:8}

def window_scan(shift_list, motif_size, sequence, sequential_decomp, sequential_part, overall_boundary, current_gap):
    shift_seq = shift_list[motif_size - 1]
    slide_size = motif_size
    
    i = current_gap[0]
    start = i
    end = current_gap[1]
    start_track = start
    initial = True
    
    while i < (current_gap[1] - slide_size + 1):
        words = shift_seq[i:i + slide_size]

        if (words.count() / len(words)) >= 0.9:
            if initial:
                calc_start = i - motif_size if i - motif_size >= start_track else start_track
                if calc_start > -1:
                    for b in overall_boundary:
                        if not (b[0] <= calc_start < b[1]):
                            pass
                        else:
                            start = b[1]
                            break
                    else:
                        start = calc_start
                else:
                    start = i
                initial = False
            end = i + motif_size
            i += motif_size
            continue
        else:
            i += 1
            if not initial:
                calc_end = end + motif_size
                for b in overall_boundary:
                    if (start < calc_end <= b[0]) or (calc_end > start >= b[1]):
                        pass
                    else:
                        end = b[0]
                        break
                else:
                    end = calc_end    
                if start != end:
                    start_track = decomposer([start, end], sequence, motif_size, sequential_decomp, sequential_part, shift_list, overall_boundary)
                initial = True
    
    if not initial:
        calc_end = end + motif_size
        for b in overall_boundary:
            if (start < calc_end <= b[0]) or (calc_end > start >= b[1]):
                pass
            else:
                end = b[0]
                break
        else:
            end = calc_end 
        if start != end:
            start_track = decomposer([start, end], sequence, motif_size, sequential_decomp, sequential_part, shift_list, overall_boundary)

def shift_decomp(seq, motif_size, motif, boundary, state):
    decomposed_parts = []
    count = 1  
    
    # Get the primary motif (will use registered cyclic variant if available)
    primary_motif, _ = get_most_frequent_motif(seq, motif_size, motif)
    
    # Search for the registered motif pattern
    positions = kmp_search_non_overlapping(seq, primary_motif)
    
    if not positions:
        return [seq], boundary

    seq = seq[positions[0] : positions[-1] + motif_size]
    
    if state:
        b1 = boundary[0]
        boundary = [b1 + positions[0], b1 + positions[-1] + motif_size]
        
    for i in range(1, len(positions)):
        if positions[i] == positions[i - 1] + len(primary_motif):
            count += 1
        else:
            decomposed_parts.append(f"({primary_motif}){count}")
            interspersed = seq[positions[i - 1] + len(primary_motif):positions[i]]
            if interspersed:
                secondary_decomp, boundary = shift_decomp(interspersed, motif_size, '', boundary, False)
                if secondary_decomp:
                    decomposed_parts.extend(secondary_decomp)
            count = 1
    
    decomposed_parts.append(f"({primary_motif}){count}")
    last_motif_end = positions[-1] + len(primary_motif)
    leftover_sequence = seq[last_motif_end:]
    
    if leftover_sequence:
        decomposed_parts.append(leftover_sequence)
    
    return decomposed_parts, boundary

def decomposer(tuple_bound, sequence, best_shift, sequential_decomp, sequential_part, shift_list, overall_boundary):
    processed_seq, tuple_bound = shift_decomp(sequence[tuple_bound[0]:tuple_bound[1]], best_shift, '', tuple_bound, True)
    b1 = tuple_bound[0]
    b2 = tuple_bound[1]
    sequential_decomp[b1] = processed_seq
    sequential_part[b1:b2] = 0
    
    for id in range(len(shift_list)):
        if shift_list[id] == 0: 
            continue
        shift_list[id][b1:b2] = 0
    
    overall_boundary.append(tuple_bound)
    return b2

def gap_boundaries(overall_boundary, seq_len, chunk_state):
    overall_boundary = sorted(overall_boundary)
    
    if chunk_state:
        chunk_boundary = overall_boundary[0][0]
        chunk_boundary_end = overall_boundary[-1][1]
    else:
        chunk_boundary = 0
        chunk_boundary_end = seq_len
        
    last_start = overall_boundary[-1][0]
    gap_regions = []
    
    for id, occupied_reg in enumerate(overall_boundary):
        if id == 0:
            if occupied_reg[0] != chunk_boundary:
                gap_regions.append([chunk_boundary, occupied_reg[0]])
            elif len(overall_boundary) == 1:
                if occupied_reg[1] <= (seq_len - 1):
                    gap_regions.append([occupied_reg[1], seq_len])
            current_end = occupied_reg[1]
        else:
            if occupied_reg[0] != current_end:
                gap_regions.append([current_end, occupied_reg[0]])
            current_end = occupied_reg[1]
            
    if current_end != chunk_boundary_end:
        gap_regions.append([current_end, chunk_boundary_end])
        
    return gap_regions
    

def motif_decomposition(sequence, motif_size):
    # Reset cyclic motif registry for each sequence
    global cyclic_motif_registry
    cyclic_motif_registry = {}
    
    shift_list = shift_and_match(sequence)
    overall_boundary = []
    sequential_decomp = {}
    sequential_part = bitarray([1] * len(sequence))
    seq_len = len(sequence)
    gap_regions = [[0, seq_len]]
    rounds = 0
    
    while any(sequential_part):
        if rounds == 1:
            shift_list[0][:] = 0
            
        if (rounds == 0) and (motif_size == 1):
            gap_wise_shift = [1]
        else:
            gap_wise_shift = max_match(shift_list, gap_regions)

        for index, best_shift in enumerate(gap_wise_shift):
            current_gap = gap_regions[index]
            
            if best_shift != -1:
                previous_ovbound = overall_boundary.copy()
                window_scan(shift_list, best_shift, sequence, sequential_decomp, sequential_part, overall_boundary, current_gap)
                
                if previous_ovbound == overall_boundary:
                    c1 = current_gap[0]
                    c2 = current_gap[1]
                    covered_chunks = [i for i in overall_boundary if (i[0] >= c1 and i[1] <= c2) or ((i[1] == c1) or (i[0] == c2))]
                    
                    if covered_chunks == []:
                        sequential_part[c1 : c2] = 0
                        sequential_decomp[c1] = [sequence[c1 : c2]]
                        overall_boundary.append([c1, c2])
                        continue

                    gap_in_chunks = gap_boundaries(covered_chunks, c2, True)
                    if not gap_in_chunks:
                        if (c1 == 0) or (c2 == seq_len):
                            gap_in_chunks = [current_gap]
                    
                    for each_gap in gap_in_chunks:
                        sequential_part[each_gap[0] : each_gap[1]] = 0
                        for shift_id in range(len(shift_list)):
                            shift_list[shift_id][each_gap[0] : each_gap[1]] = 0
                        sequential_decomp[each_gap[0]] = [sequence[each_gap[0] : each_gap[1]]]
                        overall_boundary.append([each_gap[0], each_gap[1]])
            else:
                if sequence[current_gap[0] : current_gap[1]]:
                    sequential_decomp[current_gap[0]] = [sequence[current_gap[0] : current_gap[1]]]
                sequential_part[current_gap[0] : current_gap[1]] = 0
                overall_boundary.append([current_gap[0], current_gap[1]])

        gap_regions = gap_boundaries(overall_boundary, seq_len, False)
        rounds += 1
        
        if rounds >= 20:
            break
    
    #sequential_decomp = dict([(k,v) for k,v in sequential_decomp.items() if v!=''])
    fseq = [i for k, v in sorted(sequential_decomp.items(), key=lambda x: x[0]) for i in v]
    
    if len(fseq) == 1:
        if '(' in fseq[0]:
            non_rep_percent = 0
        else:
            non_rep_percent = 1
    else:
        fseq, non_rep_percent = refine_decomposition(fseq, motif_size, len(sequence))
    
    # Return in the original format: "-".join(fseq)
    return ["-".join(fseq), non_rep_percent]

def refine_decomposition(split_seq, bed_motif, seq_len, depth=0):
    new_seq_list = []
    non_repeat = 0
    
    for i in split_seq:
        if '(' in i:
            end_point = i.index(')')
            count = int(i[end_point + 1:])
            tmp_motif = i[1:end_point]
            motif_len = len(tmp_motif)
            
            if len(set(tmp_motif)) == 1:
                new_seq_list.append(f'({tmp_motif[0]}){motif_len * count}')
            elif (motif_len % 2 == 0):
                mid_point = int(motif_len / 2)
                if tmp_motif[0: mid_point] == tmp_motif[mid_point:]:
                    check_motif = tmp_motif[0: mid_point]
                    new_seq_list.append(f'({check_motif}){count * 2}')
                else:
                    new_seq_list.append(i)
            elif motif_len <= bed_motif:
                new_seq_list.append(i)
            elif (motif_len != 1) & (motif_len > bed_motif) & (motif_len % bed_motif == 0):
                for div in divisor_dict.get(motif_len, [1]):
                    check_motif = tmp_motif[: div]
                    pattern = "(" + check_motif + "){e<=" + str(0) + "}"
                    matches = re.finditer(pattern, tmp_motif, overlapped=False)
                    tot_rep = sum(1 for _ in matches)
                    if tot_rep * len(check_motif) == motif_len:
                        motif_count = tot_rep * count
                        new_seq_list.append(f'({check_motif}){motif_count}')
                        break
                else:
                    new_seq_list.append(i)
            else:
                new_seq_list.append(i)
        else:
            new_seq_list.append(i)
            non_repeat += len(i)
    
    non_rep_percent = round(non_repeat / seq_len, 2)
    
    if split_seq == new_seq_list:
        return [new_seq_list, non_rep_percent]
    else:
        return refine_decomposition(new_seq_list, bed_motif, seq_len, depth + 1)

def read_fasta(filename):
    sequences = {}
    with open(filename, 'r') as f:
        current_id = None
        current_seq = []
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_id is not None:
                    sequences[current_id] = {
                        'sequence': ''.join(current_seq),
                        'primary_motif': primary_motif,
                        'motif_size': motif_size
                    }
                header_parts = line[1:].split(';')
                current_id = header_parts[0]
                primary_motif = header_parts[1] if len(header_parts) > 1 else ""
                motif_size = int(header_parts[2]) if len(header_parts) > 2 else 0
                current_seq = []
            else:
                current_seq.append(line)
        
        if current_id is not None:
            sequences[current_id] = {
                'sequence': ''.join(current_seq),
                'primary_motif': primary_motif,
                'motif_size': motif_size
            }
    
    return sequences

def process_fasta_to_bed(fasta_file, output_bed):
    sequences = read_fasta(fasta_file)
    
    with open(output_bed, 'w') as bed_out:
        # Original format: ID, Start, End, primary_motif, motif_size, Decomp_seq
        for seq_id, seq_info in sequences.items():
            sequence = seq_info['sequence']
            primary_motif = seq_info['primary_motif']
            motif_size = seq_info['motif_size']
            
            decomp_result = motif_decomposition(sequence, motif_size)
            decomp_seq = decomp_result[0]
            non_rep_percent = decomp_result[1]
            
            # Keep the original BED format
            bed_line = f"{seq_id}\t0\t{len(sequence)}\t{primary_motif}\t{motif_size}\t{decomp_seq}\n"
            bed_out.write(bed_line)

def main():
    parser = argparse.ArgumentParser(description='Motif Decomposition Tool - Process FASTA sequences and output BED format')
    parser.add_argument('-i', '--input', required=True, help='Input FASTA file')
    parser.add_argument('-o', '--output', required=True, help='Output BED file')
    
    args = parser.parse_args()
    
    print(f"Processing {args.input}...")
    process_fasta_to_bed(args.input, args.output)
    print(f"Results saved to {args.output}")

if __name__ == "__main__":
    main()