'''

Tools for performing and manipulating sequence alignments

'''

import scipy.stats
from skbio.alignment import StripedSmithWaterman
from .utils import reverse_complement
from cigar import Cigar

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
    """Print a basic text representation of one or more alignments.
    
    Args:
        alignments (dict or list of dicts): One or more alignments, as
            returned by the align_to_dict function. If more than one
            alignment is given, they must have the same target sequence.
            This can be useful, e.g., for looking at the forward and
            reverse primers aligned to the same sequence.
        line_length (int): Number of characters to print in each line.
            Defaults to 100.
        
    Returns:
        None
        
    """
    
    alignments = [alignments] if not isinstance(alignments, list) else alignments 
    
    ref_target = alignments[0]["target_sequence"] #reference target sequence for all alignments
    
    aligned_strings = []
    
    for alignment in alignments:
        
        # Prepare alignment view
        align_view = '-' * alignment["target_begin"] + alignment["aligned_query_sequence"]
        num_trailing_char = len(ref_target) - len(align_view)
        align_view = align_view + '-' * num_trailing_char
        
        # Take reverse complement of the alignment view if needed        
        if alignment["target_sequence"] == reverse_complement(ref_target):
            align_view = reverse_complement(align_view)
        elif alignment["target_sequence"] != ref_target:
            raise Exception("Alignments have different target sequences!")
        aligned_strings.append(align_view)

    # Print the reference sequence and aligned targets:
    for i in range (len(ref_target) // line_length + 1):
        print(ref_target[i*line_length:(i+1)*line_length])
        
        for aligned_string in aligned_strings:
            print(aligned_string[i*line_length:(i+1)*line_length])
        print(' '*line_length)
    
    return