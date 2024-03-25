'''
Tools for performing and manipulating sequence alignments
'''

import scipy.stats
from skbio.alignment import StripedSmithWaterman
from .utils import reverse_complement
from cigar import Cigar

import numpy as np
import edlib
from NanoporeAnalysis.utils import reverse_complement as rc

def align(query: str,
          subject: str,
          max_score: int,
          score_only: bool = False):
    '''
    Align sequences using edlib and compute the edge distance, i.e. steps to make
    the step sequences match. If score_only, only return the edge distance
    '''
    #print( query, subject, max_score, score_only )
    task = 'distance' if score_only else 'locations'
    alignment = edlib.align(query, subject, mode='HW', task=task, k=max_score)
    return alignment


def merge_overlapped_indices(index_pairs: list,
                             tolerated_mismatches: int) -> list:
    '''
    Provided a list of (start,end) index pairs, combine overlapping pairs into
    single pairs that span the full alignment range
    '''
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
 
    
def find_polyX(sequence: str,
               X: str,
               N: int,
               tolerated_mismatches: int,
               min_len: int) -> list:
    '''
    Find runs of single nucleotides (X), using edlib align with search query X*N
    with tolerated_mismatches in the search
    '''
    assert len(X)==1, f'More than one nt provided as quuery ({X}) for find_polyX'
    query = X*N
    alignment = align(query, sequence, max_score=0 )
    align_indices = alignment['locations']
    if len(align_indices) == 0:
        return []
    else:
        merged_indices = merge_overlapped_indices(align_indices, tolerated_mismatches)
        filtered_indices = [ i for i in merged_indices if i['length'] >= min_len]
        sorted_indices = sorted(filtered_indices, key=lambda x:x['length'], reverse=True)
        return sorted_indices
    
    
def find_best_polyA(sequence: str,
                    search_N: int=4,
                    allowed_mismatches: int=1,
                    min_len: int=10):
    '''
    Find best polyA tail
    CURRENT VERSION: use longest tail, consider adding some 3' filter in
    '''
    polyA_tails = find_polyX(sequence, 'A', search_N, allowed_mismatches, min_len )
    if len(polyA_tails) > 0:
        return polyA_tails[0]
    else:
        return None
    


def split_sequence_by_polyA(sequence: str,
                            ):
    '''
    Find the best polyA tail, and split sequence accordingly
    Returns the pre-polyA sequence, post-polyA sequence, and the polyA tail itself
    '''
    polyA = find_best_polyA(sequence)
    if polyA:
        pre_polyA = sequence[:polyA['start']]
        post_polyA = sequence[polyA['end']+1:]
    else:
        pre_polyA = sequence
        post_polyA = None
    return pre_polyA, post_polyA, polyA
        
    
    
def parse_adapter5_seq(sequence:str,
                       ssp:str,
                       max_score:int):
    '''
    Find 5' SSP and split the read accordingly. If the SSP is not found, then
    the entire sequence defaults as the transcript sequence.
    Returns the transcript sequence, adapter sequence, and SSP alignment score
    '''
    ssp_alignment = align(ssp,sequence,max_score)
    if len(ssp_alignment['locations']) > 0:
        score = (len(ssp)-ssp_alignment['editDistance'])/len(ssp)
        ssp_locations = ssp_alignment['locations']
        max_end = np.max([l[1] for l in ssp_locations])
        adapter_seq = sequence[:max_end+1]
        transcript_seq = sequence[max_end+1:]
    else:
        score = 0.5/len(ssp)
        adapter_seq = None
        transcript_seq = sequence
    return transcript_seq, adapter_seq, score
    
        
        
        
def score_umi(putative_umi: str,
              umi_dict: dict) -> dict:
    '''
    Provided a dict with {name:umi}, return the scores per umi sorted by rank
    '''
    scores = {}
    for name, umi in umi_dict.items():
        umi_align = align(rc(umi), putative_umi, max_score=10, score_only=False)
        #print(umi_align)
        if umi_align['editDistance'] == -1:
            umi_align['editDistance'] = len(umi)
        if len(umi_align['locations']) != 0 :
            umi_alignment_length = np.abs(umi_align['locations'][0][0] - umi_align['locations'][0][1])
        else :
            umi_alignment_length = 0
        umi_score = np.max([ umi_alignment_length - umi_align['editDistance'], 0.5 ])/len(umi)
        scores[name] = umi_score
    scores = dict(sorted(scores.items(), key=lambda x:x[1], reverse=True))
    return scores


def parse_adapter3_seq(sequence: str,
                       primer: str,
                       umi_dict: dict,
                       max_align_score: int = 10):
    '''
    Identify sequence characteristics in the 3' adapter sequence, including any potential UMIs
    '''
    primer_alignment = align(primer, sequence, max_score=max_align_score)
    primer_locations = primer_alignment['locations']
    #print( sequence, primer, primer_alignment )
    if len(primer_locations) > 0:
        primer_score = (len(primer)-primer_alignment['editDistance'])/len(primer)
        post_start = np.min([l[0] for l in primer_locations])
    else:
        primer_score = 0.5/len(primer)
    post_start = len(sequence) # Retain full sequence for UMI alignment
    putative_umi = sequence[:post_start]
    umi_scores = score_umi(putative_umi, umi_dict)
    return primer_score, umi_scores
    

def annotate_read(read: dict, 
                  ssp: str, 
                  primer3: str, 
                  umis: dict,
                  max_polyA_score_len: float = 20.0):
    '''
    Search a read for a PolyA tail and then score the 5'/3' ends accordingly
    
    *** POSSIBLE PATCH FOR LATER ***
    Currently max_polyA_score_len is arbitrary, consider changing later
    ***
    
    '''
    sequence = read['sequence']
    pre_polyA, post_polyA, polyA = split_sequence_by_polyA(sequence)
    if post_polyA:
        # Score potential UMIs in 3' adapter sequence
        primer_3_score, umi_scores = parse_adapter3_seq(post_polyA, primer3, umis)
        polyA_score = np.min([polyA['length']/max_polyA_score_len, 1.0])
    else:
        primer_3_score = 0.5/len(primer3)
        umi_scores = dict([(umi,0.5/len(umi)) for umi in umis])
        polyA_score = 0.5/max_polyA_score_len
    transcript_seq, adapter5_seq, primer_5_score = parse_adapter5_seq(pre_polyA, ssp, 10) # HARD CODED MAX SCORE ATM, CHANGE LATER
    read_annotations = { '3\' primer':                 primer3,
                         '3\' primer alignment score': primer_3_score,
                         '3\' adapter sequence':       post_polyA,
                         'Strand switch primer':       ssp,
                         'SSP alignment score':        primer_5_score,
                         '5\' adapter sequence':       adapter5_seq,
                         'UMI Scores':                 umi_scores,
                         'PolyA coordinates':          polyA,
                         'PolyA score':                polyA_score,
                         'Transcript sequence':        transcript_seq,
                       }
    return read | read_annotations
    
    
def process_read(read: dict,
                 ssp: str,
                 primer3: str,
                 umis: dict,
                 min_polyA_score: float = 0.5,
                 min_ssp_score: float = 0.5,
                 min_umi_score: float = 0.5):
    '''
    Annotate reads and, if possible, re-orient along forward strand
    
    *** POSSIBLE PATCH FOR LATER ***
    Currently min_polyA_score and min_ssp_score are *very* arbitrarily chosen
    based on examination of real data; usually real scores are substantially
    better than 0.5, but seemed like a good enough starting point. Consider
    changing based on further testing.
    ***
    '''
    forward_analysis = annotate_read(read, ssp, primer3, umis)
    rc_read = dict(read)
    rc_read['sequence'] = rc(rc_read['sequence'])
    reverse_analysis = annotate_read(rc_read, ssp, primer3, umis)
    
    if np.max([list(forward_analysis['UMI Scores'].values())[0],list(reverse_analysis['UMI Scores'].values())[0]]) > min_umi_score:
        # Found a PolyA tail
        delta_polyA_score = forward_analysis['PolyA score'] - reverse_analysis['PolyA score']
        if max([forward_analysis['PolyA score'], reverse_analysis['PolyA score']]) > min_polyA_score:
            analysis = forward_analysis if delta_polyA_score > 0 else reverse_analysis
            annotation = 'Full-length cDNA' if analysis['SSP alignment score'] > min_ssp_score else "3\' only"
        else:
            analysis = forward_analysis # Need to actually solve for the correct one
            annotation = "Duplexed 3\' only"
        best_umi, second_umi = list(analysis['UMI Scores'])[:2]
        best_umi_score = analysis['UMI Scores'][best_umi]
        delta_umi_score = best_umi_score - analysis['UMI Scores'][second_umi]
        
    else:
        # Did not find a PolyA tail, check for 5' features
        best_umi, best_umi_score, delta_umi_score = [None]*3
        max_ssp_score = np.max([forward_analysis['SSP alignment score'],reverse_analysis['SSP alignment score']])
        if max_ssp_score > min_ssp_score:
            delta_ssp_score = forward_analysis['SSP alignment score']-reverse_analysis['SSP alignment score']
            if np.abs(delta_ssp_score) > min_ssp_score:
                analysis = forward_analysis if delta_ssp_score > 0 else reverse_analysis
                annotation = "5\' only"
            else:
                analysis = forward_analysis # Need to solve for which one is better
                annotation = "Duplexed 5\' only"
        else:
            analysis = forward_analysis # Unknown, no better guess
            annotation = "Uncharacterized sequence"
    analysis['annotation'] = annotation 
    analysis['Best UMI'] = best_umi
    analysis['Best UMI Score'] = best_umi_score
    analysis['Delta UMI Score'] = delta_umi_score
    return analysis, forward_analysis, reverse_analysis

def align_to_dict(align_result):
    """Convert a skbio.alignment object into a simple dict.
    
    Args:
        align_result (class 'skbio.alignment._ssw_wrapper.AlignmentStructure'):
            A computed alignment returned by one of several alignment
            functions available in the skbio.alignment package.
        
    Returns:
        align_dict (dict): A dictionary of all public attributes of align_result.
    """
    
    align_dict ={}
    
    for x in dir(align_result):
        if x[0] != '_' and type(align_result[x]) in [str,int]:
            align_dict[x] = align_result[x]
    
    return align_dict

def reverse_alignment(alignment):
    """
    Takes a scikitbio alignment after being passed through align_to_dict and gives the reverse complement of all sequences, as well as changes the beginning and ending indexes to reflect the reversal. Note that the target_end_suboptimal cannot be changed to reflect the suboptimal alignment in the reversed direction, but its indices are still reversed here.
    
    Returns:
        Reverse-complemented dict
    """
    
    alignment_rev = alignment.copy()
    alignment_rev['aligned_query_sequence'] = reverse_complement(alignment['aligned_query_sequence'])
    alignment_rev['aligned_target_sequence'] = reverse_complement(alignment['aligned_target_sequence'])
    alignment_rev['query_sequence'] = reverse_complement(alignment['query_sequence'])
    alignment_rev['target_sequence'] = reverse_complement(alignment['target_sequence'])
    alignment_rev['query_begin'] = len(alignment['query_sequence']) - 1 - alignment['query_end']
    alignment_rev['query_end'] = len(alignment['query_sequence']) - 1 - alignment['query_begin']
    alignment_rev['target_begin'] = len(alignment['target_sequence']) - 1 - alignment['target_end_optimal']
    alignment_rev['target_end_optimal'] = len(alignment['target_sequence']) - 1 - alignment['target_begin']
    alignment_rev['target_end_suboptimal'] = len(alignment['target_sequence']) - 1 - alignment['target_end_suboptimal']
    alignment_rev['cigar'] = Cigar(alignment['cigar'])._reverse_cigar()
    return alignment_rev

def get_best_barcode(target_sequence, primers, strand_switch_primer, filter_strand_switch_score, seq_id='ID_not_specified', tags = None, score_filter=None, distance_filter=None):
    """
    Get the best alignment for a barcode sequence and a matching strand-switching primer against a target sequence
    First align the strand switching primer, then align the barcodes to the resulting local regions of the target sequence. Returns a scikitbio stripedsmithwaterman alignment with the alignment parameters modified to reflect an alignment done against the full target, ie adjust the alignment start/end points to accurately reflect the position of the alignment in the whole target instead of the local region derived from the initial primer alignment.
    
    Args:
        target_sequence (str): the sequence to be aligned to.
        primers (dict): a dict of barcodes with the keys as barcode names and the values as the DNA sequences with only characters ATCGatcg.
        strand_switch_primer (str): the primer used for library prep which is reverse complement to a concensus part of all barcodes.
        filter_strand_switch_score (int): filter for the strand-switch primer alignment. If no alignment is found above this filter, the barcode is aligned against the entire target sequence.
        seq_id='ID_not_specified' (str): ID for the target sequence. Usually a NGS read ID
        score_filter=None (int): score filter for barcode alignment. Speeds up alignment up to ~10-20%, but filters >90% of the maximum filter for the given barcode length.
        distance_filter=None (int): filter for the minimum aligned length for the barcodes. Similar to but less effective than score_filter.
        
    Returns:
        A dict with seven items: 
            'alignment' = dict of barcode alignment after being re-indexed,
            'direction' = direction of barcode alignment
            'primer_name' = name of best aligned barcode
            'score' = barcode alignment score
            'strand_switch_primer_align' = alignment dict for strand_switch_primer
            'seq_id' = name of target sequence
            'biological_sequence' = the region of target sequence between the strand_switch_primer and the barcode or just everything up to the barcode if the filter_strand_switch_score wasn't reached.
    """
    
    max_score = 0
    result = {}
    
    query = StripedSmithWaterman(reverse_complement(strand_switch_primer))
    fwd_ss_align_dict = align_to_dict(query(target_sequence))
    fwd_start = max(0, fwd_ss_align_dict['target_begin'] - 50)
    fwd_end = max(0, fwd_ss_align_dict['target_end_optimal'] + 50)
    rev_ss_align_dict = align_to_dict(query(reverse_complement(target_sequence)))
    rev_start = max(0, rev_ss_align_dict['target_begin'] - 50)
    rev_end = max(0, rev_ss_align_dict['target_end_optimal'] + 50)
    
    if fwd_ss_align_dict['optimal_alignment_score'] >= filter_strand_switch_score:
        fwd_seq_region = target_sequence[fwd_start:fwd_end]
    else:
        fwd_seq_region = target_sequence
        
    if rev_ss_align_dict['optimal_alignment_score'] >= filter_strand_switch_score:
        rev_seq_region = reverse_complement(target_sequence)[rev_start:rev_end]
    else:
        rev_seq_region = reverse_complement(target_sequence)
    
    for name, seq in primers.items():
        #try aligning in both directions
        query = StripedSmithWaterman(seq, score_filter=score_filter, distance_filter=distance_filter)
        
        fwd = query(fwd_seq_region)
        rev = query(rev_seq_region)

        if fwd.optimal_alignment_score > rev.optimal_alignment_score:
            alignment = fwd
            direction = 'forward'
        else:
            alignment = rev
            direction = 'reverse'
    
        score = alignment.optimal_alignment_score
        
        if score > max_score or max_score == 0:
            result['alignment'] = align_to_dict(alignment)
            result['direction'] = direction
            result['primer_name'] = name
            result['score'] = score
            max_score = score
    if result['direction'] == 'forward':
        result['strand_switch_primer_align'] = reverse_alignment(rev_ss_align_dict)
        initial_ss_align_score = fwd_ss_align_dict['optimal_alignment_score']
    else :
        result['strand_switch_primer_align'] = reverse_alignment(fwd_ss_align_dict)
        initial_ss_align_score = rev_ss_align_dict['optimal_alignment_score']
        
    if initial_ss_align_score >= filter_strand_switch_score :
        if result['direction'] == 'forward':
            result['alignment']['target_begin'] = result['alignment']['target_begin'] + fwd_start
            result['alignment']['target_end_optimal'] = result['alignment']['target_end_optimal'] + fwd_start
            result['alignment']['target_end_suboptimal'] = result['alignment']['target_end_suboptimal'] + fwd_start
            result['alignment']['target_sequence'] = target_sequence
        else :
            result['alignment']['target_begin'] = result['alignment']['target_begin'] + rev_start
            result['alignment']['target_end_optimal'] = result['alignment']['target_end_optimal'] + rev_start
            result['alignment']['target_end_suboptimal'] = result['alignment']['target_end_suboptimal'] + rev_start
            result['alignment']['target_sequence'] = reverse_complement(target_sequence)
    
    if result['strand_switch_primer_align']['optimal_alignment_score'] >= filter_strand_switch_score :
        result['biological_sequence'] = result['alignment']['target_sequence'][result['strand_switch_primer_align']['target_end_optimal'] + 1 : result['alignment']['target_begin']]
    else :
        result['biological_sequence'] = result['alignment']['target_sequence'][ 0 : result['alignment']['target_begin']]
    result['barcode_score'] = result['alignment']['optimal_alignment_score']
    result['strand_switch_score'] = result['strand_switch_primer_align']['optimal_alignment_score']
    result['read_length'] = len(target_sequence)
    result['seq_id'] = seq_id.replace('@','')
    result['tags'] = tags
    
    return(result)

def get_best_alignment(target_sequence, primers, seq_id='ID_not_specified'):
        
    """Attempt to align every primer to target_sequence and return the best result.
    
    Uses the striped Smith-Waterman algorithm as implemented in the
    skbio.alignment package.
    
    Args:
        target_sequence (str): sequence to which the primers are aligned.
        primers (dict): primers to try aligning, with primer names as keys
            and corresponding sequences as values.
        seq_id (str): optional parameter to uniquely identify the target sequence.
            Defaults to 'ID_not_specified'.
        
    Returns:
        result (dict): the highest-scoring alignment, with the following keys:
            'alignment' :  a dict of the best-scoring alignment, as returned by
                the align_to_dict function.
            'direction' : direction of the best-scoring alignment; can be 
                'forward' or 'reverse'.
            'primer_name' : name of the best-scoring primer.
            'score' : the score of the best-scoring alignment.
            'seq_id' : unique identifier of the target sequence.
    """

    max_score = 0
    result = {}
    
    for name, seq in primers.items():
        #try aligning in both directions
        query = StripedSmithWaterman(seq)
        fwd = query(target_sequence)
        rev = query(reverse_complement(target_sequence))

        if fwd.optimal_alignment_score > rev.optimal_alignment_score:
            alignment = fwd
            direction = 'forward'
        else:
            alignment = rev
            direction = 'reverse'
    
        score = alignment.optimal_alignment_score
        
        if score > max_score or max_score == 0:
            result['alignment'] = align_to_dict(alignment)
            result['direction'] = direction
            result['primer_name'] = name
            result['score'] = score
            max_score = score
    
    result['seq_id'] = seq_id

    return(result)


def show_alignment_simple(alignments, line_length=100):
    """
    Print a basic text representation of one or more alignments. This will account for indels by creating a global reference sequence from all alignments. Aligned portions of each query are shown in uppercase letters and the unaligned regions are shown in lowercase, cut off at the ends of the target sequence to make visualization easy.
    
    Args:
        alignments (dict or list of dicts): One or more alignments, as
            returned by the align_to_dict function. If more than one
            alignment is given, they must have the same target sequence (reverse complement is also accepted).
            This can be useful, e.g., for looking at the forward and
            reverse primers aligned to the same sequence.
        line_length (int): Number of characters to print in each line.
            Defaults to 100.
        
    Returns:
        None
        
    """
    
    # Assert the input alignments as a list
    alignments = [alignments] if not isinstance(alignments, list) else alignments
    ref_target = alignments[0]["target_sequence"]
    
    # Set up a variable to hold the aligned queries and targets
    aligned_strings = []
    
    # Go through each alignment and set up the query and target as aligned strings.
    for alignment in alignments:
        
        # Prepare alignment view
        # Decide if the query sequence should have a gap, where the target extends beyond the query, or extend the query to the start of the target.
        if alignment['target_begin'] >= alignment['query_begin'] :
            align_view_query = ' ' * ( alignment['target_begin'] - alignment['query_begin'] ) + alignment["query_sequence"][:alignment["query_begin"]].lower() + alignment["aligned_query_sequence"].upper()
        elif alignment['target_begin'] < alignment['query_begin'] :
            align_view_query = alignment["query_sequence"][alignment['query_begin'] - alignment['target_begin']:alignment["query_begin"]].lower() + alignment["aligned_query_sequence"].upper() + alignment["query_sequence"][ alignment["query_end"] : alignment["query_end"] + len(alignment['target_sequence']) - alignment['target_end_optimal'] ].lower()
            
        # Add the unaligned beginning of the target sequence, if needed.
        align_view_query_target = alignment['target_sequence'][:alignment["target_begin"]] + alignment["aligned_target_sequence"]
        
        # Add a trailing space if needed to make all sequences a similar length.
        trailing_char = ' ' * ( len(ref_target) - alignment['target_end_optimal'] )
        align_view_query = align_view_query + trailing_char
        align_view_query_target = align_view_query_target + trailing_char
        
        # Take reverse complement of the alignment view if needed, and check if target is the same for every alignment.    
        if alignment["target_sequence"] == reverse_complement(ref_target):
            align_view_query = reverse_complement(align_view_query)
            align_view_query_target = reverse_complement(align_view_query_target)
        elif alignment["target_sequence"].lower() != ref_target.lower():
            raise Exception("Alignments have different target sequences!")
            
        # Append both the query and target to the strings holder.
        aligned_strings.append([align_view_query, align_view_query_target])
    
    # Begin setting up a global alignment. Start a residue counter as an index. The aligner will continue through the aligned strings until it reaches the length of the 'blank' ref_target, which will get indels inserted anywhere that an alignment shows one.
    residue = 0
    while residue < len(ref_target) :
        for i in range(len(aligned_strings)) :
            # If there's an insert at the current index position for any of the aligned target sequences, insert a '-' into the ref_target and any aligned query that doesn't have one there.
            if aligned_strings[i][1][residue] == '-' :
                if ref_target[residue] != '-' :
                    ref_target = ref_target[:residue] + '-' + ref_target[residue:]
                for j in range(len(aligned_strings)) :
                    if aligned_strings[j][1][residue] != '-' :
                        aligned_strings[j][1] = aligned_strings[j][1][:residue] + '-' + aligned_strings[j][1][residue:]
                        aligned_strings[j][0] = aligned_strings[j][0][:residue] + '-' + aligned_strings[j][0][residue:]
        residue += 1
                        
    # Print the reference sequence and aligned targets:
    for i in range( (len(ref_target) // line_length ) + 1) :
        #Print the global target sequence.
        print(ref_target[i*line_length:(i+1)*line_length])
        
        # Set up a blank boolean to mark whether or not the current line has any actual aligned sequences or if no sequences (included unaligned regions) show up yet. Also create a print holder variable to hold the seq until after everything has been tested to be blank or not.
        blank = True
        aligned_prints = []
        for aligned_string in aligned_strings :
            line = aligned_string[0][i*line_length:(i+1)*line_length]
            aligned_prints.append(line)
            # If any alignment has something in this line, mark blank as False.
            if line.replace(' ', '') != str('') :
                blank = False
        if blank == False : 
            for line in aligned_prints :
                print(line)
        # Print a gap
        print(' '*line_length)
        
    return