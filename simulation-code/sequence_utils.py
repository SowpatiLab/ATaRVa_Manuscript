# sequence_utils.py
def recalculate_tracking_from_sequence(sequence, primary_motif):
    """
    Recalculate tracking by actually scanning the sequence for motifs
    Returns accurate tracking that matches the actual sequence
    """
    motif_len = len(primary_motif)
    i = 0
    tracking_parts = []
    
    # Generate all possible fragments from the primary motif
    possible_motifs = []
    
    # Primary motif
    possible_motifs.append((primary_motif, motif_len))
    
    # All fragments (size 2 to motif_len-1)
    for frag_size in range(2, motif_len):
        for start in range(motif_len - frag_size + 1):
            fragment = primary_motif[start:start+frag_size]
            possible_motifs.append((fragment, frag_size))
    
    # Sort by length (longest first) for greedy matching
    possible_motifs.sort(key=lambda x: x[1], reverse=True)
    
    # Also consider monomers
    monomers = ['A', 'C', 'G', 'T']
    
    while i < len(sequence):
        matched = False
        
        # Try all possible motifs (longest first)
        for motif, motif_len_val in possible_motifs:
            if i + motif_len_val <= len(sequence):
                # Check for consecutive repeats of this motif
                count = 0
                j = i
                while (j + motif_len_val <= len(sequence) and 
                       sequence[j:j+motif_len_val] == motif):
                    count += 1
                    j += motif_len_val
                
                if count > 0:
                    tracking_parts.append(f"{motif}:{count}")
                    i = j
                    matched = True
                    break
        
        if not matched:
            # Try monomers
            current_base = sequence[i]
            if current_base in monomers:
                count = 0
                j = i
                while j < len(sequence) and sequence[j] == current_base:
                    count += 1
                    j += 1
                
                if count > 0:
                    tracking_parts.append(f"{current_base}:{count}")
                    i = j
                    matched = True
        
        if not matched:
            # Skip one base and try again
            i += 1
    
    # Combine consecutive identical motifs
    if tracking_parts:
        final_parts = []
        current_motif = None
        current_count = 0
        
        for part in tracking_parts:
            motif, count_str = part.split(":")
            count = int(count_str)
            
            if motif == current_motif:
                current_count += count
            else:
                if current_motif is not None:
                    final_parts.append(f"{current_motif}:{current_count}")
                current_motif = motif
                current_count = count
        
        if current_motif is not None:
            final_parts.append(f"{current_motif}:{current_count}")
        
        return "|".join(final_parts)
    
    return ""

def verify_sequence_vs_tracking(sequence, tracking_str, repeat_id=""):
    """Verify tracking matches sequence"""
    reconstructed = ""
    parts = tracking_str.split("|")
    
    for part in parts:
        if ":" in part:
            motif, count_str = part.split(":")
            count = int(count_str)
            reconstructed += motif * count
        else:
            reconstructed += part
    
    if reconstructed == sequence:
        return True
    else:
        print(f"ERROR {repeat_id}: Tracking mismatch!")
        print(f"  Sequence length: {len(sequence)}")
        print(f"  Reconstructed length: {len(reconstructed)}")
        return False
