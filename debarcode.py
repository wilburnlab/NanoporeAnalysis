from pathlib import Path
import concurrent
import concurrent.futures
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq
from pyarrow import dataset
from pyarrow import csv
import numpy as np
import edlib
from NanoporeAnalysis import utils

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
    for_alignment = edlib.align(query_seq, target_seq, mode='HW', task='locations')
    if for_alignment['editDistance'] <= max_edit_distance and for_alignment['editDistance'] >= 0 :
        for location in for_alignment['locations'] :
            location_overlapped = False
            if block_overlap :
                for match in matches :
                    if location[0] < match['end'] and location[1] > match['start'] :
                        location_overlapped = True
            if not location_overlapped :
                if location[0] == None :
                    alignment_start = 0
                else :
                    alignment_start = location[0]
                alignment_length = abs( alignment_start - location[1] ) + 1
                if alignment_length >= min_length :
                    score = ( alignment_length - for_alignment['editDistance'] ) / query_seq_len
                    matches.append({
                        'query_ID' : query_ID,
                        'edit_distance' : for_alignment['editDistance'],
                        'edit_score' : score,
                        'start' : alignment_start,
                        'end' : location[1],
                        'direction' : 'forward'
                    })
    if not skip_reverse :
        rev_alignment = edlib.align(utils.reverse_complement(query_seq), target_seq, mode='HW', task='locations')
        if rev_alignment['editDistance'] <= max_edit_distance and rev_alignment['editDistance'] >= 0 :
            for location in rev_alignment['locations'] :
                location_overlapped = False
                if block_overlap :
                    for match in matches :
                        if location[0] < match['end'] and location[1] > match['start'] :
                            location_overlapped = True
                if not location_overlapped :
                    if location[0] == None :
                        alignment_start = 0
                    else :
                        alignment_start = location[0]
                    alignment_length = abs( alignment_start - location[1] ) + 1
                    if alignment_length >= min_length :
                        score = ( alignment_length - rev_alignment['editDistance'] ) / query_seq_len
                        matches.append({
                            'query_ID' : query_ID,
                            'edit_distance' : rev_alignment['editDistance'],
                            'edit_score' : score,
                            'start' : alignment_start,
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
        if end <= current_end + 1 + tolerated_mismatches:
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

def evaluate_polyA_UMI_SSP(polyA, UMI_match, SSP_match, polyA_UMI_SSP_set, max_gap = 3) :
    """
    Judge a set of polyA, UMI, and SSP alignments based on whether they're in the correct order (polyA, then UMI, then SSP) and are close enough to each other via max_gap. Uses polyA_UMI_SSP_set as a blank template and fills it in. First calculates gaps, then checks if directions of UMI and polyA are matching and their gap is allowable, then fills out the template polyA_UMI_SSP_set from them. 
    
    Args :
        polyA (list): the indices of a polyA sequence from find_polyA_bidirectional(). This is a tuple of the start and end indices [start, end].
        UMI_match (dict) : the alignment for a UMI, output by find_seq_matches(). This is a dictionary of the alignment attributes.
        SSP_match (dict) : the alignment for the SSP sequence, output by find_seq_matches(). This is a dictionary of the alignment attributes.
        polyA_UMI_SSP_set (dict) : a template dict with default values for the keys. The default values are intended to indicate an invalid match, as failing to meet the requirements of this function will result in the default values being returned.
        max_gap (int) : the maximum gap between the UMI and the polyA and SSP. Defaults to 5.
    
    Returns :
        polyA_UMIT_SSP_set (dict) : the template dict with changed values, corresponding to whether the three alignments form a valid set representing a true barcode.
    """
    polyA_UMI_gap = abs(polyA['end'] - UMI_match['location'][0])
    UMI_SSP_gap = abs(SSP_match['location'][0] - UMI_match['location'][1])
    if UMI_match['direction'] == polyA['direction'] and polyA_UMI_gap <= max_gap :
        polyA_UMI_SSP_set['barcode_ID'] = UMI_match['query_ID']
        polyA_UMI_SSP_set['barcode_polyA_start'] = polyA['start']
        polyA_UMI_SSP_set['barcode_polyA_end'] = polyA['end']
        polyA_UMI_SSP_set['barcode_polyA_len'] = polyA['length']
        polyA_UMI_SSP_set['barcode_UMI_start'] = UMI_match['location'][0]
        polyA_UMI_SSP_set['barcode_UMI_end'] = UMI_match['location'][1]
        polyA_UMI_SSP_set['barcode_UMI_edit_distance'] = UMI_match['edit_distance']
        polyA_UMI_SSP_set['barcode_direction'] = UMI_match['direction']
        polyA_UMI_SSP_set['barcode_flag'][0] = 1
        polyA_UMI_SSP_set['barcode_flag'][1] = 1
    if UMI_match['direction'] == SSP_match['direction'] and UMI_SSP_gap <= max_gap :
        polyA_UMI_SSP_set['barcode_ID'] = UMI_match['query_ID']
        polyA_UMI_SSP_set['barcode_UMI_start'] = UMI_match['location'][0]
        polyA_UMI_SSP_set['barcode_UMI_end'] = UMI_match['location'][1]
        polyA_UMI_SSP_set['barcode_UMI_edit_distance'] = UMI_match['edit_distance']
        polyA_UMI_SSP_set['barcode_SSP_start'] = SSP_match['location'][0]
        polyA_UMI_SSP_set['barcode_SSP_end'] = SSP_match['location'][1]
        polyA_UMI_SSP_set['barcode_SSP_edit_distance'] = SSP_match['edit_distance']
        polyA_UMI_SSP_set['barcode_direction'] = UMI_match['direction']
        polyA_UMI_SSP_set['barcode_flag'][1] = 1
        polyA_UMI_SSP_set['barcode_flag'][2] = 1
        polyA_UMI_SSP_set['barcode_score'] = polyA_UMI_SSP_set['barcode_UMI_edit_distance'] + polyA_UMI_SSP_set['barcode_SSP_edit_distance']
    return polyA_UMI_SSP_set
    
def reverse_matches(matches, seq_len) :
    """
    Takes match outputs from find_seq_matches and reverses the orientation by subtracting the location indices from the length of the sequence. All matches must at least correspond to the same sequence length, usually the same sequence as well.
    
    Args :
        matches (list) : list of the matches to be reversed.
        seq_len (int) : length of the sequence that the matches correspond to.
    
    Returns :
        reversed_matches (list) : list of matches after being reversed. They are in the same order as given in matches.
    """
    reversed_matches = []
    for match in matches :
        reversed_matches.append({
        'query_ID' : match['query_ID'],
        'edit_distance' : match['edit_distance'],
        'edit_score' : match['edit_score'],
        'location' : [ seq_len - 1 - match['location'][1], seq_len - 1 - match['location'][0] ],
        'direction' : match['direction']
    })
    return reversed_matches

def pick_best_UMI(seq, UMIs, UMI_max_score, UMI_min_len, direction = 'forward') :
    UMI_matches = []
    for UMI in UMIs :
        if direction == 'reverse' :
            UMI_match = find_seq_matches(utils.reverse_complement(UMI[1]), seq, UMI_max_score, UMI[0], min_length = UMI_min_len, skip_reverse = True)
        else :
            UMI_match = find_seq_matches(UMI[1], seq, UMI_max_score, UMI[0], min_length = UMI_min_len, skip_reverse = True)
        if len(UMI_match) > 0 :
            UMI_matches.append(UMI_match[0])
    if len(UMI_matches) > 1 :
        UMI_matches.sort( key = lambda x : x['edit_score'], reverse = True )
        if UMI_matches[0]['edit_score'] - UMI_matches[1]['edit_score'] < 0.05 :
            return False
    elif len(UMI_matches) == 0 :
        return False
    match_dict =  {
        'barcode_ID' : UMI_matches[0]['query_ID'],
        'barcode_UMI_start' : UMI_matches[0]['start'],
        'barcode_UMI_end' : UMI_matches[0]['end'],
        'barcode_UMI_edit_distance' : UMI_matches[0]['edit_distance'],
        'barcode_direction' : direction,
        'barcode_flag' : [0,1,0],
    }
    return match_dict

def parse_polyA_UMI_SSP(seq, UMIs, SSP, UMI_max_score, SSP_max_score, max_gap = 5, UMI_min_len = 0, SSP_min_len = 0, polyA_N = 4, polyA_tolerated_mismatches = 1, polyA_min_len = 10) :
    """
    Finds potential unique molecular identifiers (UMIs), strand switching primer (SSP), and polyA sequences that meet certain minimum quality criteria, then determines if a UMI can be paired with an SSP and/or polyA alignment to represent a true barcode sequence. 
    
    Args :
        seq (str) : the sequence to be checked.
        UMIs (list) : list of [name, sequence] pairs for UMIs.
        SSP (str) : sequence for the SSP.
        UMI_max_score (int) : maximum edit_distance for the UMI alignments. Any alignments with a higher edit_distance are discarded.
        SSP_max_score (int) : maximum edit_distance for the SSP alignments. Any alignments with a higher edit_distance are discarded.
        max_gap (int) : the maximum gap between the UMI and the polyA or SSP. Defaults to 5
        UMI_min_len (int) : minimum length for the UMI alignments.
        SSP_min_len (int) : minimum length for the SSP alignments.
        polyA_N (int) : number of characters to be considered a polyA sequence. Defaults to 4.
        polyA_tolerated_mismatches (int) : number of allowed gap characters between polyA alignments that will be grouped together into a larger polyA. Defaults to 1.
        polyA_min_len (int) : minimum length of the final polyX region(s) that will be considered a true polyA. Defaults to 10.
    
    Returns :
        barcode_dict (dict) : dictionary with the final decided attributes for the seq and whether it contains an identifiable UMI/barcode.
    """
    max_UMI_len = max( [len(x[1]) for x in UMIs] )
    seq_len = len(seq)
    SSP_matches = find_seq_matches(SSP, seq, SSP_max_score, 'SSP', min_length = SSP_min_len)
    polyAs = find_polyA_bidirectional(seq, N = polyA_N, tolerated_mismatches = polyA_tolerated_mismatches, min_len = polyA_min_len)
    all_features = SSP_matches + polyAs
    all_features.sort( key = lambda x : x['start'] )
    possible_barcodes = []
    possible_distal_SSPs = []
    last_feature = {'query_ID' : 'None', 'direction' : None}
    for feature in all_features :
        if feature['query_ID'] + last_feature['query_ID'] in [ 'SSPpolyA', 'polyASSP' ] and feature['direction'] == last_feature['direction'] and feature['start'] - last_feature['end'] <= max_UMI_len + max_gap + max_gap :
            if feature['query_ID'] == 'SSP' :
                possible_barcodes[-1]['barcode_SSP_start'] = feature['start']
                possible_barcodes[-1]['barcode_SSP_end'] = feature['end']
                possible_barcodes[-1]['barcode_SSP_edit_distance'] = feature['edit_distance']
                possible_barcodes[-1]['barcode_flag'][2] = 1
            elif feature['query_ID'] == 'polyA' :
                possible_barcodes[-1]['barcode_polyA_start'] = feature['start']
                possible_barcodes[-1]['barcode_polyA_end'] = feature['end']
                possible_barcodes[-1]['barcode_polyA_len'] = feature['length']
                possible_barcodes[-1]['barcode_flag'][0] = 1
            last_feature = feature
        elif feature['query_ID'] == 'SSP' :
            if feature['direction'] == 'forward' :
                best_UMI = pick_best_UMI(seq[ feature['start'] - max_gap - max_UMI_len : feature['start'] + max_gap ], UMIs, UMI_max_score, UMI_min_len, direction = 'forward')
                if best_UMI :
                    best_UMI['barcode_SSP_start'] = feature['start']
                    best_UMI['barcode_SSP_end'] = feature['end']
                    best_UMI['barcode_SSP_edit_distance'] = feature['edit_distance']
                    best_UMI['barcode_flag'][2] = 1
                    best_UMI['barcode_UMI_start'] += feature['start'] - max_gap - max_UMI_len
                    best_UMI['barcode_UMI_end'] += feature['start'] - max_gap - max_UMI_len
                    possible_barcodes.append(best_UMI)
                    last_feature = feature
                else :
                    possible_distal_SSPs.append(feature)
            if feature['direction'] == 'reverse' :
                best_UMI = pick_best_UMI(seq[ feature['end'] - max_gap : feature['end'] + max_gap + max_UMI_len ], UMIs, UMI_max_score, UMI_min_len, direction = 'reverse')
                if best_UMI :
                    best_UMI['barcode_SSP_start'] = feature['start']
                    best_UMI['barcode_SSP_end'] = feature['end']
                    best_UMI['barcode_SSP_edit_distance'] = feature['edit_distance']
                    best_UMI['barcode_flag'][2] = 1
                    best_UMI['barcode_UMI_start'] += feature['end'] - max_gap
                    best_UMI['barcode_UMI_end'] += feature['end'] - max_gap
                    possible_barcodes.append(best_UMI)
                    last_feature = feature
                else :
                    possible_distal_SSPs.append(feature)
        elif feature['query_ID'] == 'polyA' :
            if feature['direction'] == 'forward' :
                best_UMI = pick_best_UMI(seq[ feature['end'] - max_gap : feature['end'] + max_gap + max_UMI_len ], UMIs, UMI_max_score, UMI_min_len, direction = 'forward')
                if best_UMI :
                    best_UMI['barcode_polyA_start'] = feature['start']
                    best_UMI['barcode_polyA_end'] = feature['end']
                    best_UMI['barcode_polyA_len'] = feature['length']
                    best_UMI['barcode_flag'][0] = 1
                    best_UMI['barcode_UMI_start'] += feature['end'] - max_gap
                    best_UMI['barcode_UMI_end'] += feature['end'] - max_gap
                    possible_barcodes.append(best_UMI)
                    last_feature = feature
            if feature['direction'] == 'reverse' :
                best_UMI = pick_best_UMI(seq[ feature['start'] - max_gap - max_UMI_len : feature['start'] + max_gap ], UMIs, UMI_max_score, UMI_min_len, direction = 'reverse')
                if best_UMI :
                    best_UMI['barcode_polyA_start'] = feature['start']
                    best_UMI['barcode_polyA_end'] = feature['end']
                    best_UMI['barcode_polyA_len'] = feature['length']
                    best_UMI['barcode_flag'][0] = 1
                    best_UMI['barcode_UMI_start'] += feature['start'] - max_gap - max_UMI_len
                    best_UMI['barcode_UMI_end'] += feature['start'] - max_gap - max_UMI_len
                    possible_barcodes.append(best_UMI)
                    last_feature = feature
    if len(possible_barcodes) == 0 :
        return {'barcode_ID': 'none matched', 'barcode_flag' : [0,0,0]}
    elif len(possible_barcodes) > 1 :
        return {'barcode_ID':'multiple', 'barcode_flag' : [0,0,0]}
    else :
        barcode_dict = possible_barcodes[0]
        left_barcode_indices = [ barcode_dict['barcode_UMI_start'] ]
        if 'barcode_polyA_start' in barcode_dict :
            left_barcode_indices.append(barcode_dict['barcode_polyA_start'])
        if 'barcode_SSP_start' in barcode_dict :
            left_barcode_indices.append(barcode_dict['barcode_SSP_start'])
        left_barcode_index = min(left_barcode_indices)
        right_barcode_indices = [ barcode_dict['barcode_UMI_end'] ]
        if 'barcode_polyA_end' in barcode_dict :
            left_barcode_indices.append(barcode_dict['barcode_polyA_end'])
        if 'barcode_SSP_end' in barcode_dict :
            left_barcode_indices.append(barcode_dict['barcode_SSP_end'])
        right_barcode_index = max(right_barcode_indices)
        if barcode_dict['barcode_direction'] == 'forward' :
            barcode_dict['barcode_biological_seq_indices'] = [ 1, left_barcode_index - 1 ]
        else :
            barcode_dict['barcode_biological_seq_indices'] = [ right_barcode_index + 1, seq_len - 1 ]
        for possible_distal_SSP in possible_distal_SSPs :
            if barcode_dict['barcode_direction'] == 'forward' and possible_distal_SSP['direction'] == 'reverse' :
                if possible_distal_SSP['end'] < left_barcode_index :
                    barcode_dict['barcode_biological_seq_indices'][0] = possible_distal_SSP['end'] + 1
                    barcode_dict['barcode_distal_SSP_start'] = possible_distal_SSP['start']
                    barcode_dict['barcode_distal_SSP_end'] = possible_distal_SSP['end']
                    barcode_dict['barcode_distal_SSP_edit_distance'] = possible_distal_SSP['edit_distance']
            elif barcode_dict['barcode_direction'] == 'reverse' and possible_distal_SSP['direction'] == 'forward' :
                if possible_distal_SSP['start'] > right_barcode_index :
                    barcode_dict['barcode_biological_seq_indices'][1] = possible_distal_SSP['start'] - 1
                    barcode_dict['barcode_distal_SSP_start'] = possible_distal_SSP['start']
                    barcode_dict['barcode_distal_SSP_end'] = possible_distal_SSP['end']
                    barcode_dict['barcode_distal_SSP_edit_distance'] = possible_distal_SSP['edit_distance']
        return barcode_dict

def debarcode_table(table, UMIs, SSP, UMI_max_score, SSP_max_score, max_gap = 5, UMI_min_len = 0, SSP_min_len = 0, polyA_N = 4, polyA_tolerated_mismatches = 1, polyA_min_len = 10, sample_dict = None, include_duplex_parents = False) :
    """
    Takes a pyarrow table with one sequences in each row under a column titled 'seq' and runs parse_polyA_UMI_SSP() on all the reads within.
    
    Args :
        table (pyarrow table) : table with a column called 'seq' with the sequences to be analyzed. Seqs must be strings.
        UMIs (list) : list of [name, sequence] pairs for UMIs.
        SSP (str) : sequence for the SSP.
        UMI_max_score (int) : maximum edit_distance for the UMI alignments. Any alignments with a higher edit_distance are discarded.
        SSP_max_score (int) : maximum edit_distance for the SSP alignments. Any alignments with a higher edit_distance are discarded.
        max_gap (int) : the maximum gap between the UMI and the polyA or SSP. Defaults to 5
        UMI_min_len (int) : minimum length for the UMI alignments.
        SSP_min_len (int) : minimum length for the SSP alignments.
        polyA_N (int) : number of characters to be considered a polyA sequence. Defaults to 4.
        polyA_tolerated_mismatches (int) : number of allowed gap characters between polyA alignments that will be grouped together into a larger polyA. Defaults to 1.
        polyA_min_len (int) : minimum length of the final polyX region(s) that will be considered a true polyA. Defaults to 10.
    
    Returns :
        table (pyarrow table) : the input table with the output from parse_polyA_UMI_SSP() added as new columns.
    """
    seqs = table.column('seq').to_pylist()
    if 'dx:i' in table.column_names :
        duplex_flags = table.column('dx:i').to_pylist()
    else :
        duplex_flags = pa.nulls(len(seqs)).fill_nulls(0).to_pylist()
    parsed_seqs_dict = {
        'barcode_ID' : [],
        'barcode_polyA_start' : [],
        'barcode_polyA_end' : [],
        'barcode_polyA_len' : [],
        'barcode_UMI_start' : [],
        'barcode_UMI_end' : [],
        'barcode_UMI_edit_distance' : [],
        'barcode_SSP_start' : [],
        'barcode_SSP_end' : [],
        'barcode_SSP_edit_distance' : [],
        'barcode_direction' : [],
        'barcode_flag' : [],
        'barcode_biological_seq_indices' : [],
        'barcode_distal_SSP_start' : [],
        'barcode_distal_SSP_end' : [],
        'barcode_distal_SSP_edit_distance' : []
    }
    barcode_sample_IDs = []
    for seq, duplex_flag in zip( seqs, duplex_flags ) :
        if duplex_flag != -1 or include_duplex_parents :
            barcode_dict = parse_polyA_UMI_SSP(seq, UMIs, SSP, UMI_max_score, SSP_max_score, max_gap, UMI_min_len, SSP_min_len, polyA_N, polyA_tolerated_mismatches, polyA_min_len)
        else :
            barcode_dict = {}
        for key in parsed_seqs_dict :
            if key in barcode_dict :
                parsed_seqs_dict[key].append(barcode_dict[key])
            else :
                parsed_seqs_dict[key].append(None)
        if barcode_dict['barcode_ID'] not in ['multiple', 'none matched'] and sample_dict != None :
            barcode_sample_IDs.append(sample_dict[barcode_dict['barcode_ID']])
        else :
            barcode_sample_IDs.append('none')
    for key in parsed_seqs_dict :
        table = table.append_column(key, [parsed_seqs_dict[key]])
    table = table.append_column('barcode_sample_ID', [barcode_sample_IDs])
    return table

def debarcode_table_from_file(file, UMIs, SSP, UMI_max_score, SSP_max_score, max_gap = 5, UMI_min_len = 0, SSP_min_len = 0, polyA_N = 4, polyA_tolerated_mismatches = 1, polyA_min_len = 10, sample_dict = None, resume=False, overwrite=False) :
    """
    Load a pyarrow table from a parquet file, then runs passes it to debarcode_table() which runs parse_polyA_UMI_SSP() which attempts to find a barcode in all the sequences within the file. This then saves the table back to file, overwriting the original.
    
    Args :
        UMIs (list) : list of [name, sequence] pairs for UMIs.
        SSP (str) : sequence for the SSP.
        UMI_max_score (int) : maximum edit_distance for the UMI alignments. Any alignments with a higher edit_distance are discarded.
        SSP_max_score (int) : maximum edit_distance for the SSP alignments. Any alignments with a higher edit_distance are discarded.
        max_gap (int) : the maximum gap between the UMI and the polyA or SSP. Defaults to 5
        UMI_min_len (int) : minimum length for the UMI alignments.
        SSP_min_len (int) : minimum length for the SSP alignments.
        polyA_N (int) : number of characters to be considered a polyA sequence. Defaults to 4.
        polyA_tolerated_mismatches (int) : number of allowed gap characters between polyA alignments that will be grouped together into a larger polyA. Defaults to 1.
        polyA_min_len (int) : minimum length of the final polyX region(s) that will be considered a true polyA. Defaults to 10.
        resume (bool) : whether to resume debarcoding from a paused or broken run. This will only debarcode files that don't already have barcode information in them, so it can't be used to continue a re-debarcoding session, as the files will all still have the original barcode data. Defaults to False
        overwrite (bool) : whether to allow for overwriting debarcoded files. If True, it will remove any existing barode data and continue with normal debarcoding. If False, it will skip over any files with barcode data.
    """
    table = pq.read_table(file)
    if 'barcode_ID' in table.column_names :
        if resume == True :
            print("skip resume", file)
            print(table.column_names)
            del table
            return
        elif overwrite == True :
            print("removing old barcode data...")
            columns_to_drop = [ column for column in table.column_names if 'SSP_' in column or 'barcode_' in column ]
            if len(columns_to_drop) > 0 :
                table = table.drop_columns(columns_to_drop)
        else :
            del table
            raise ValueError("Error: some of this data may have already been debarcoded. Please set resume = True if resuming or set overwrite = True if you'd like to overwrite any existing barcode information.")
    print("debarcoding file ", Path(file).stem)
    table = debarcode_table(table, UMIs, SSP, UMI_max_score, SSP_max_score, max_gap, UMI_min_len, SSP_min_len, polyA_N, polyA_tolerated_mismatches, polyA_min_len, sample_dict)
    pq.write_table(table, file)
    del table
    return

def debarcode(path_dataset, UMIs, SSP, UMI_max_score = 3, SSP_max_score = 5, max_gap = 5, UMI_min_len = 12, SSP_min_len = 15, polyA_N = 4, polyA_tolerated_mismatches = 1, polyA_min_len = 10, workers = 1, sample_dict = None, resume=False, overwrite=False) :
    """
    Finds all parquet files within the given dataset directory under path_dataset, then pushes them through debarcode_table_from_file() to attempt to find barcodes.
    
    Args :
        path_dataset (str) : the path to the folder containing the parquet files to be debarcoded. Should be path_out/pa_dataset/ where path_out is the same as what was used in build_parquet_dataset_from_sam.
        UMIs (str) : list of [name, sequence] pairs for UMIs (unique molecular identifiers). These are the unique sequences of the barcodes. These should be given as the reverse complement of the actual sequence of the barcode primers,
            such that the full barcode would be (polyA sequence)-(UMI)-(SSP sequence). This is crucial for proper identification of barcodes.
        SSP (str) : sequence for the SSP (strand-switching primer), which should be the same sequence as the constant region on the barcodes.
        UMI_max_score (int) : maximum edit_distance for the UMI alignments. Any alignments with a higher edit_distance are discarded.
        SSP_max_score (int) : maximum edit_distance for the SSP alignments. Any alignments with a higher edit_distance are discarded.
        max_gap (int) : the maximum gap between the UMI and the polyA or SSP. Defaults to 5
        UMI_min_len (int) : minimum length for the UMI alignments.
        SSP_min_len (int) : minimum length for the SSP alignments.
        polyA_N (int) : number of characters to be considered a polyA sequence. Defaults to 4.
        polyA_tolerated_mismatches (int) : number of allowed gap characters between polyA alignments that will be grouped together into a larger polyA. Defaults to 1.
        polyA_min_len (int) : minimum length of the final polyX region(s) that will be considered a true polyA. Defaults to 10.
        workers (int) : number of parallel debarcoding processes to use. Defaults to 4.
        resume (bool) : whether to resume debarcoding from a paused or broken run. This will only debarcode files that don't already have barcode information in them, so it can't be used to continue a re-debarcoding session, as the files will all still have the original barcode data. Defaults to False
        overwrite (bool) : whether to allow for overwriting debarcoded files. If True, it will remove any existing barode data and continue with normal debarcoding. If False, it will skip over any files with barcode data.
    """
    files = [x for x in Path(path_dataset).iterdir() if x.is_file()]
    with concurrent.futures.ProcessPoolExecutor( max_workers=workers ) as executor :
        futures = [ executor.submit( debarcode_table_from_file, file, UMIs, SSP, UMI_max_score, SSP_max_score, max_gap, UMI_min_len, SSP_min_len, polyA_N, polyA_tolerated_mismatches, polyA_min_len, sample_dict, resume, overwrite ) for file in files ]
        concurrent.futures.wait( futures )
    print("Finished debarcoding!")
    return