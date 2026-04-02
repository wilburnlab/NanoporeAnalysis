def find_seq_matches(query_seq, target_seq, max_edit_distance, query_ID = 'None', min_length = 0, skip_reverse = False, block_overlap = True) :
    """
    Finds matches between target_seq and query_seq using edlib.align in both forward and reverse direction, then translates the result into a list of location matches.
    
    Args :
        query_seq (str) : sequence to be used as query. This is usually the smaller sequence and within the target.
        target_seq (str) : sequence to be used as target. This is usually the larger sequence, and the query is meant to be within it.
        max_edit_distance (int) : the maximum allowable editDistance from an edlib alignment. Anything higher will be discarded. The edit distance is a count of how many characters have to be changed to turn the query into the target.
        query_ID (str) : the name of the query sequence. Defaults to 'None'.
        min_length (int) : the minmum allowed length of an alignment. Anything shorter will be discarded.
        skip_reverse (bool) : whether or not to skip the reverse alignment. Only useful for when aligning two things that have known orientations, ie the distal SSP after assigning the barcode. Defaults to False.
        
    Returns :
        matches (list) : a list of dictionaries for each alignment that meets the requirments.
    """
    matches = []
    target_seq_len = len(target_seq)
    query_seq_len = len(query_seq)
    forward_alignment = edlib.align(query_seq, target_seq, mode='HW', task='locations')
    if forward_alignment['editDistance'] <= max_edit_distance :
        for location in forward_alignment['locations'] :
            location_overlapped = False
            if block_overlap :
                for match in matches :
                    if abs( np.mean(location) - np.mean([ match['start'], match['end'] ]) ) < 10 :
                        location_overlapped = True
            if not location_overlapped :
                if location[0] == None :
                    alignment_end = 0
                else :
                    alignment_end = location[0]
                alignment_length = abs( alignment_end - location[1] ) + 1
                if alignment_length >= min_length :
                    score = ( alignment_length - forward_alignment['editDistance'] ) / query_seq_len
                    matches.append({
                        'query_ID' : query_ID,
                        'edit_distance' : forward_alignment['editDistance'],
                        'edit_score' : score,
                        'start' : alignment_end,
                        'end' : location[1],
                        'direction' : 'forward'
                    })
    if not skip_reverse :
        rev_alignment = edlib.align(utils.reverse_complement(query_seq), target_seq, mode='HW', task='locations')
        if rev_alignment['editDistance'] <= max_edit_distance :
            for location in rev_alignment['locations'] :
                location_overlapped = False
                if block_overlap :
                    for match in matches :
                        if abs( np.mean(location) - np.mean([ match['start'], match['end'] ]) ) < 10 :
                            location_overlapped = True
                if not location_overlapped :
                    if location[0] == None :
                        alignment_length = abs( 0 - location[1] ) + 1
                    else :
                        alignment_length = abs( location[0] - location[1] ) + 1
                    if alignment_length >= min_length :
                        score = ( alignment_length - rev_alignment['editDistance'] ) / query_seq_len
                        matches.append({
                            'query_ID' : query_ID,
                            'edit_distance' : rev_alignment['editDistance'],
                            'edit_score' : score,
                            'start' : location[0],
                            'end' : location[1],
                            'direction' : 'reverse'
                        })
    return matches

def merge_overlapped_indices(index_pairs: list, tolerated_mismatches: int) -> list:
    """
    Provided a list of (start,end) index pairs, combine overlapping pairs into single pairs that span the full alignment range.
    
    Args :
        index_pairs (list) : list of [start (int),end (int)] tuples.
        tolerated_mismatches (int) : number of allowed gap characters between alignments that will still be grouped together.
        
    Returns :
        reduced_pairs (list) : list of condensed [start (int),end (int)] tuples.
    """
    reduced_pairs = []
    for i, pair in enumerate(index_pairs):
        if i == 0:
            current_start, current_end = pair
        start, end = pair
        if end <= current_end+1+tolerated_mismatches:
            current_end = end
        else:
            align_len = current_end - current_start + 1
            reduced_pairs.append({'start':current_start,'end':current_end,'length':align_len})
            current_start = start
            current_end = end
    align_len = current_end - current_start + 1
    reduced_pairs.append({'start':current_start,'end':current_end,'length':align_len})
    return reduced_pairs
    
def find_polyX(sequence: str, X: str, N: int, tolerated_mismatches: int, min_len: int) -> list:
    """
    Find runs of single nucleotides (X), using edlib align with search query X*N with tolerated_mismatches in the search
    
    Args :
        sequence (str) : sequence to be searched.
        X (str) : single character to be used. Must be just 1 character long, not case sensitive.
        N (int) : number of characters to be considered a polyX sequence.
        tolerated_mismatches (int) : number of allowed gap characters between alignments that will still be grouped together.
        min_len (int) : minimum length of the final polyX region(s).
        
    Returns :
        sorted_indices (list) : list of [start (int),end (int)] tuples sorted by length.
    """
    assert len(X)==1, f'More than one nt provided as query ({X}) for find_polyX'
    query = X*N
    alignment = edlib.align(query, sequence, mode='HW', task='locations')
    align_indices = alignment['locations']
    if len(align_indices) == 0:
        return []
    else:
        merged_indices = merge_overlapped_indices(align_indices, tolerated_mismatches)
        filtered_indices = [ i for i in merged_indices if i['length'] >= min_len]
        sorted_indices = sorted(filtered_indices, key=lambda x:x['length'], reverse=True)
        return sorted_indices
    
def find_polyA_bidirectional(sequence: str, N: int, tolerated_mismatches: int, min_len: int) -> list:
    """
    Finds any polyA sequences in the sequence in either forward or reverse direction that meets the minimum length.
    
    Args :
        sequence (str) : sequence to be searched.
        N (int) : number of characters to be considered a polyX sequence.
        tolerated_mismatches (int) : number of allowed gap characters between alignments that will still be grouped together.
        min_len (int) : minimum length of the final polyX region(s).
        
    Returns :
        indices (list) : list of [start (int),end (int)] tuples.
    """
    for_indices = find_polyX(sequence, 'A', N, tolerated_mismatches, min_len)
    rev_indices = find_polyX(sequence, 'T', N, tolerated_mismatches, min_len)
    indices = []
    for index_set in for_indices :
        index_set['direction'] = 'forward'
        index_set['query_ID'] = 'polyA'
        indices.append(index_set)
    for index_set in rev_indices :
        index_set['direction'] = 'reverse'
        index_set['query_ID'] = 'polyA'
        indices.append(index_set)
    return indices