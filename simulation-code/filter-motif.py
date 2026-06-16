import re

def has_consecutive_repeats(motif):
    """Check if motif has consecutive repeated nucleotides"""
    for i in range(len(motif) - 1):
        if motif[i] == motif[i + 1]:
            return True
    return False

def has_repeated_pattern(motif, min_pattern_length=2):
    """Check if motif contains repeated substring patterns"""
    n = len(motif)
    
    # Check for repeated substrings of different lengths
    for length in range(min_pattern_length, n//2 + 1):
        for i in range(n - length*2 + 1):
            substring1 = motif[i:i+length]
            substring2 = motif[i+length:i+2*length]
            if substring1 == substring2:
                return True
    return False

# Process file
input_file = "filtered-HG38_2-10_motifs_d2d.tsv"
output_file = "filtered-2-HG38_2-10_motifs_d2d.tsv"

with open(input_file, 'r') as f_in, open(output_file, 'w') as f_out:
    for line in f_in:
        line = line.strip()
        if not line:
            continue
            
        parts = line.split()
        if len(parts) >= 1:
            motif = parts[0]
            
            # Keep only if NO consecutive repeats
            if not has_repeated_pattern(motif):
                f_out.write(line + '\n')
