import random
from interrupted-TR-simulation import *


def generate_complex_repeat_with_tracking(primary_motif, total_length, sequence_type="fragmented", fragment_mutation_rate=0.8):
    """
    Generate complex repeat with tracking that EXACTLY matches the final sequence
    sequence_type: "pure", "point_mutation", or "fragmented"
    fragment_mutation_rate: probability of mutating a fragment when added (0.0-1.0)
    """
    motif_len = len(primary_motif)
    sequence = ""
    
    # Track everything we try to insert
    insertion_sequence = []
    fragment_count = 0
    
    # Start with 3 repeats of the primary motif
    initial_repeats = 3
    sequence += primary_motif * initial_repeats
    insertion_sequence.append(f"{primary_motif}:{initial_repeats}")
    
    if sequence_type == "pure":
        # Just fill with primary motif
        remaining_repeats = (total_length - len(sequence)) // motif_len
        if remaining_repeats > 0:
            sequence += primary_motif * remaining_repeats
            insertion_sequence.append(f"{primary_motif}:{remaining_repeats}")
    
    elif sequence_type == "fragmented":
        # Available fragment sizes: 2 to motif_len - 1
        if motif_len > 2:
            valid_fragment_sizes = list(range(2, motif_len))
        else:
            valid_fragment_sizes = []
        
        # Add loop safety counter
        max_iterations = total_length * 3
        iteration_count = 0
        
        while len(sequence) < total_length and iteration_count < max_iterations:
            iteration_count += 1
            remaining_length = total_length - len(sequence)
            
            # More conservative break condition
            if remaining_length < motif_len:
                break
                
            choice = random.random()
            made_progress = False
            
            if choice < 0.7:  # Primary motif (70%)
                max_repeats = remaining_length // motif_len
                if max_repeats > 0:
                    repeats = random.randint(1, min(4, max_repeats))
                    new_segment = primary_motif * repeats
                    # Only add if it fits completely
                    if len(sequence) + len(new_segment) <= total_length:
                        sequence += new_segment
                        insertion_sequence.append(f"{primary_motif}:{repeats}")
                        made_progress = True
                
            elif choice < 0.9:  # Fragment (20%)
                if valid_fragment_sizes and remaining_length >= 6:  # Need space for at least 3 copies
                    fragment_size = random.choice(valid_fragment_sizes)
                    max_start_pos = motif_len - fragment_size
                    
                    if max_start_pos >= 0:
                        start_pos = random.randint(0, max_start_pos)
                        fragment = primary_motif[start_pos:start_pos + fragment_size]
                        max_repeats = remaining_length // len(fragment)
                        
                        # ENSURE AT LEAST 3 COPIES
                        if max_repeats >= 3:
                            repeats = random.randint(3, min(6, max_repeats))  # Min 3 copies
                            new_segment = fragment * repeats
                            # Only add if it fits completely
                            if len(sequence) + len(new_segment) <= total_length:
                                
                                # CHECK IF WE SHOULD MUTATE THIS FRAGMENT (80% chance)
                                should_mutate =  fragment_mutation_rate
                                if should_mutate:
                                    # Mutate the fragment using cyclic variation
                                    mutated_fragment = apply_cyclic_mutation(fragment)
                                    new_segment = mutated_fragment * repeats
                                    # Update tracking to show mutation
                                    insertion_sequence.append(f"{mutated_fragment}:{repeats}")
                                else:
                                    insertion_sequence.append(f"{fragment}:{repeats}")
                                
                                sequence += new_segment
                                fragment_count += 1
                                made_progress = True
                
            else:  # Mixed fragments (10%)
                if valid_fragment_sizes and remaining_length >= 9:  # Need space for at least 3 copies
                    fragment_size = random.choice(valid_fragment_sizes)
                    max_start_pos = motif_len - fragment_size
                    if max_start_pos >= 0:
                        start_pos = random.randint(0, max_start_pos)
                        fragment = primary_motif[start_pos:start_pos + fragment_size]
                        segment_length = random.randint(motif_len + 1, motif_len * 3)
                        segment_length = min(segment_length, remaining_length)
                        repeats = segment_length // len(fragment)
                        
                        # ENSURE AT LEAST 3 COPIES
                        if repeats >= 3:
                            new_segment = fragment * repeats
                            new_segment = new_segment[:segment_length]
                            # Only add if it fits completely
                            if len(sequence) + len(new_segment) <= total_length:
                                
                                # CHECK IF WE SHOULD MUTATE THIS FRAGMENT (80% chance)
                                should_mutate = random.random() < fragment_mutation_rate
                                if should_mutate:
                                    # Mutate the fragment using cyclic variation
                                    mutated_fragment = apply_cyclic_mutation(fragment)
                                    new_segment = mutated_fragment * repeats
                                    new_segment = new_segment[:segment_length]
                                    insertion_sequence.append(f"{mutated_fragment}:{repeats}")
                                else:
                                    insertion_sequence.append(f"{fragment}:{repeats}")
                                
                                sequence += new_segment
                                fragment_count += 1
                                made_progress = True
            
            # If we didn't make progress for a while, break to avoid infinite loop
            if not made_progress and remaining_length < (motif_len * 2):
                break
    
    # Final padding with complete motifs only if we're close
    if len(sequence) < total_length:
        remaining = total_length - len(sequence)
        if remaining >= motif_len:
            fits_repeats = remaining // motif_len
            if fits_repeats >= 1:
                final_segment = primary_motif * fits_repeats
                sequence += final_segment
                
                # Update tracking
                if insertion_sequence and insertion_sequence[-1].startswith(primary_motif + ":") and not "mutated" in insertion_sequence[-1]:
                    last_entry = insertion_sequence[-1]
                    current_repeats = int(last_entry.split(":")[1])
                    insertion_sequence[-1] = f"{primary_motif}:{current_repeats + fits_repeats}"
                else:
                    insertion_sequence.append(f"{primary_motif}:{fits_repeats}")
    
    original_length = len(sequence)
    sequence = sequence[:total_length]
    trimmed_length = original_length - len(sequence)
    
    # ADJUST TRACKING TO MATCH ACTUAL SEQUENCE LENGTH
    if trimmed_length > 0:
        adjusted_sequence = adjust_tracking_to_length(insertion_sequence, trimmed_length, primary_motif)
    else:
        adjusted_sequence = insertion_sequence
    
    # Combine consecutive primary motif entries and simplify monomers
    combined_sequence = []
    current_primary_count = 0
    current_monomer = None
    current_monomer_count = 0
    
    for entry in adjusted_sequence:
        if entry.startswith(primary_motif + ":") and not "mutated" in entry:
            repeats = int(entry.split(":")[1])
            current_primary_count += repeats
        else:
            # Check if this is a monomer entry
            if ":" in entry and not "mutated" in entry:
                parts = entry.split(":")
                motif = parts[0]
                count = int(parts[1])
                
                # Check if it's a single nucleotide (monomer)
                if len(motif) == 1:
                    if motif == current_monomer:
                        current_monomer_count += count
                    else:
                        if current_monomer_count > 0:
                            combined_sequence.append(f"{current_monomer}:{current_monomer_count}")
                        current_monomer = motif
                        current_monomer_count = count
                    continue
            
            # Handle non-monomer entries
            if current_primary_count > 0:
                combined_sequence.append(f"{primary_motif}:{current_primary_count}")
                current_primary_count = 0
            
            if current_monomer_count > 0:
                combined_sequence.append(f"{current_monomer}:{current_monomer_count}")
                current_monomer = None
                current_monomer_count = 0
            
            combined_sequence.append(entry)
    
    # Don't forget the last accumulated entries
    if current_primary_count > 0:
        combined_sequence.append(f"{primary_motif}:{current_primary_count}")
    if current_monomer_count > 0:
        combined_sequence.append(f"{current_monomer}:{current_monomer_count}")
    
    return sequence, "|".join(combined_sequence), fragment_count


def apply_cyclic_mutation(fragment):
    """
    Apply cyclic mutation to a fragment (A→C→G→T→A)
    """
    CYCLIC_MAP = {'A': 'C', 'C': 'G', 'G': 'T', 'T': 'A'}
    mutated_fragment = ""
    for nucleotide in fragment:
        mutated_fragment += CYCLIC_MAP.get(nucleotide, nucleotide)
    return mutated_fragment


def adjust_tracking_to_length(insertion_sequence, trimmed_length, primary_motif):
    """
    Adjust tracking by trimming from the end to match actual sequence length
    """
    if trimmed_length <= 0 or not insertion_sequence:
        return insertion_sequence
    
    adjusted_sequence = insertion_sequence.copy()
    remaining_trim = trimmed_length
    
    while remaining_trim > 0 and adjusted_sequence:
        last_entry = adjusted_sequence[-1]
        
        if ":" in last_entry and "mutated" not in last_entry:
            parts = last_entry.split(":")
            motif = parts[0]
            count = int(parts[1])
            
            # Calculate total length represented by this entry
            if len(motif) == 1:  # Monomer
                total_length = count
            else:
                total_length = len(motif) * count
            
            if total_length <= remaining_trim:
                # Remove entire entry
                adjusted_sequence.pop()
                remaining_trim -= total_length
            else:
                # Trim this entry
                keep_length = total_length - remaining_trim
                
                if len(motif) == 1:  # Monomer
                    adjusted_sequence[-1] = f"{motif}:{keep_length}"
                else:
                    complete_repeats = keep_length // len(motif)
                    partial_bases = keep_length % len(motif)
                    
                    if complete_repeats > 0:
                        if partial_bases > 0:
                            adjusted_sequence[-1] = f"{motif}:{complete_repeats}+{partial_bases}"
                        else:
                            adjusted_sequence[-1] = f"{motif}:{complete_repeats}"
                    else:
                        if partial_bases > 0:
                            adjusted_sequence[-1] = f"{motif[:partial_bases]}"
                        else:
                            adjusted_sequence.pop()
                
                remaining_trim = 0
        else:
            # Remove non-counted entries
            adjusted_sequence.pop()
            remaining_trim -= 1
    
    return adjusted_sequence