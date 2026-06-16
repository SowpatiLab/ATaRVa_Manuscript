# point_mutation.py
import random
from config import NUCLEOTIDES

def mutate_sequence_substitutions_only(sequence, purity=0.90):
    """
    Add only substitution mutations to a sequence
    """
    seq_length = len(sequence)
    num_mutations = int((1 - purity) * seq_length)
    
    if num_mutations <= 0:
        return sequence, []
    
    mutations_info = []
    seq_list = list(sequence)
    
    # Select unique positions for mutations
    positions = random.sample(range(seq_length), min(num_mutations, seq_length))
    
    for pos in sorted(positions):
        original_base = seq_list[pos]
        possible_bases = list(NUCLEOTIDES - {original_base})
        if possible_bases:
            new_base = random.choice(possible_bases)
            seq_list[pos] = new_base
            mutations_info.append(['S', str(pos), f'{original_base}/{new_base}'])
    
    return ''.join(seq_list), mutations_info

def mutate_sequence_with_motif_spacing(sequence, primary_motif, purity=0.90):
    """
    Add substitutions with spacing of approximately 1 mutation per 4 motifs
    PROTECTS first 3 motifs from mutations
    """
    motif_len = len(primary_motif)
    seq_length = len(sequence)
    
    if motif_len == 0 or seq_length == 0:
        return sequence, []
    
    # Calculate number of motifs in sequence
    num_motifs = seq_length // motif_len
    
    # Target: 1 mutation per 4 motifs (excluding first 3)
    motifs_available_for_mutation = max(0, num_motifs - 3)  # Exclude first 3 motifs
    target_mutations = max(1, motifs_available_for_mutation // 4)
    
    # Adjust based on purity if it results in fewer mutations
    purity_mutations = int((1 - purity) * seq_length)
    num_mutations = min(target_mutations, purity_mutations)
    
    if num_mutations <= 0:
        return sequence, []
    
    mutations_info = []
    seq_list = list(sequence)
    
    # Distribute mutations evenly across motifs (excluding first 3)
    mutation_positions = []
    
    if motifs_available_for_mutation > 0:
        # Motif indices start from 0, so exclude motifs 0, 1, 2 (first 3 motifs)
        available_motif_indices = list(range(3, num_motifs))
        
        # Calculate how many motifs we can actually mutate
        num_motifs_to_mutate = min(num_mutations, len(available_motif_indices))
        
        if num_motifs_to_mutate > 0:
            # Randomly select which motifs get mutations
            motifs_to_mutate = random.sample(available_motif_indices, num_motifs_to_mutate)
            
            for motif_index in motifs_to_mutate:
                # Pick a random position within this motif
                start_pos = motif_index * motif_len
                end_pos = min(start_pos + motif_len, seq_length)
                
                if start_pos < end_pos:
                    # Pick a random position within this motif
                    pos_in_motif = random.randint(0, min(motif_len - 1, end_pos - start_pos - 1))
                    mutation_pos = start_pos + pos_in_motif
                    mutation_positions.append(mutation_pos)
    
    # If we need more mutations or no motifs found, add random ones (still excluding first 3 motifs)
    if len(mutation_positions) < num_mutations:
        additional_needed = num_mutations - len(mutation_positions)
        
        # Calculate positions that are NOT in the first 3 motifs
        protected_region_end = 3 * motif_len
        
        # All available positions (excluding protected region)
        available_positions = []
        for pos in range(seq_length):
            # Skip positions in first 3 motifs
            if pos >= protected_region_end:
                available_positions.append(pos)
        
        # Remove positions already selected for mutation
        available_positions = [p for p in available_positions if p not in mutation_positions]
        
        if available_positions:
            additional_positions = random.sample(available_positions, 
                                               min(additional_needed, len(available_positions)))
            mutation_positions.extend(additional_positions)
    
    # Apply substitutions at selected positions
    for pos in sorted(mutation_positions):
        original_base = seq_list[pos]
        possible_bases = list(NUCLEOTIDES - {original_base})
        if possible_bases:
            new_base = random.choice(possible_bases)
            seq_list[pos] = new_base
            mutations_info.append(['S', str(pos), f'{original_base}/{new_base}'])
    
    return ''.join(seq_list), mutations_info

def update_tracking_with_substitutions(tracking_str, primary_motif, mutations_info, sequence):
    """
    Update tracking to reflect substitution mutations in the sequence
    Example: TACG:10|TAC:3|TACG:1|CG:4|TACG:2 + S|66|A/G → TACG:10|TAC:3|TACG:1|CG:4|TGCG:1|TACG:1
    """
    if not mutations_info:
        return tracking_str

    
    
    # Convert mutations to dict by position
    mutation_dict = {}
    for mut_type, pos_str, info in mutations_info:
        if mut_type == 'S':  # Only handle substitutions
            pos = int(pos_str)
            mutation_dict[pos] = info
    
    if not mutation_dict:
        return tracking_str
    
    # Parse tracking parts
    tracking_parts = tracking_str.split("|")
    
    # Build position map and reconstruct motifs
    motif_blocks = []
    current_pos = 0
    
    for part in tracking_parts:
        if ":" in part:
            motif, count_str = part.split(":")
            count = int(count_str)
            motif_len = len(motif)
            
            for i in range(count):
                start = current_pos
                end = current_pos + motif_len
                
                # Check if this motif instance has any mutations
                mutated_motif = list(motif)
                has_mutation = False
                
                for pos in range(start, end):
                    if pos in mutation_dict:
                        # Apply the substitution
                        original, new = mutation_dict[pos].split('/')
                        pos_in_motif = pos - start
                        if pos_in_motif < len(mutated_motif) and mutated_motif[pos_in_motif] == original:
                            mutated_motif[pos_in_motif] = new
                            has_mutation = True
                
                # Store the motif (original or mutated)
                final_motif = ''.join(mutated_motif)
                motif_blocks.append({
                    'motif': final_motif,
                    'start': start,
                    'end': end
                })
                current_pos = end
        else:
            # Single character (monomer)
            motif_blocks.append({
                'motif': part,
                'start': current_pos,
                'end': current_pos + 1
            })
            current_pos += 1
    
    # Combine consecutive identical motifs
    new_tracking_parts = []
    current_motif = None
    current_count = 0
    
    for block in motif_blocks:
        if block['motif'] == current_motif:
            current_count += 1
        else:
            if current_motif is not None:
                if current_count > 1:
                    new_tracking_parts.append(f"{current_motif}:{current_count}")
                else:
                    new_tracking_parts.append(current_motif)
            current_motif = block['motif']
            current_count = 1
    
    # Add the last motif
    if current_motif is not None:
        if current_count > 1:
            new_tracking_parts.append(f"{current_motif}:{current_count}")
        else:
            new_tracking_parts.append(current_motif)
    
    return "|".join(new_tracking_parts)

def is_monomer_repeat(motif):
    """Check if motif is repeated single nucleotides (AA, CCC, etc.)"""
    if len(motif) <= 1:
        return None
    
    # More Pythonic implementation
    return motif[0] if all(char == motif[0] for char in motif) else None