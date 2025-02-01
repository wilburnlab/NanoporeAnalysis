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
               min_len: int,
               max_len) -> list:
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
        filtered_indices = [ i for i in merged_indices if i['length'] >= min_len and i['length'] <= max_len]
        sorted_indices = sorted(filtered_indices, key=lambda x:x['start'], reverse=True) # Put the last one first
        return sorted_indices
    
    
def find_best_polyA(sequence: str,
                    search_N: int=4,
                    allowed_mismatches: int=1,
                    min_len: int=15,
                    max_len: int=40): # Patch to min length = 15, max = 40
    '''
    Find best polyA tail
    CURRENT VERSION: use longest tail, consider adding some 3' filter in
    '''
    polyA_tails = find_polyX(sequence, 'A', search_N, allowed_mismatches, min_len, max_len)
    if len(polyA_tails) > 0:
        return polyA_tails[0]
    else:
        return None
    


def split_sequence_by_polyA(sequence: str):
    '''
    Find the best polyA tail, and split sequence accordingly
    Returns the pre-polyA sequence, post-polyA sequence, and the polyA tail itself
    '''
    polyA = find_best_polyA(sequence)
    if polyA is not None:
        pre_polyA = sequence[:polyA['start']]
        post_polyA = sequence[polyA['end']+1:]
        polyA_stats = dict([(f"PolyA {k}", v) for k,v in polyA.items()])
    else:
        pre_polyA = sequence
        post_polyA = None
        polyA_stats = dict()
    return pre_polyA, post_polyA, polyA_stats
        
    
    
def parse_adapter5_seq(sequence:str,
                       ssp:str,
                       max_score:int = 2):
    '''
    Find 5' SSP and split the read accordingly. If the SSP is not found, then
    the entire sequence defaults as the transcript sequence.
    Returns the transcript sequence, adapter sequence, and SSP alignment score
    '''
    if sequence is None:
        score, trim_idx, adapter_seq = [None]*3
        transcript_seq = sequence
    else:
        ssp_alignment = align(ssp,sequence,max_score)
        if len(ssp_alignment['locations']) > 0:
            score = ssp_alignment['editDistance'] #(len(ssp)-ssp_alignment['editDistance'])/len(ssp)
            ssp_locations = ssp_alignment['locations']
            max_end = np.max([l[1] for l in ssp_locations])
            trim_idx = max_end+1
            adapter_seq = sequence[:trim_idx]
            transcript_seq = sequence[trim_idx:]
        else:
            score, trim_idx, adapter_seq = [None]*3
            transcript_seq = sequence
            #score = None #0.5/len(ssp)
            #trim_idx = None
        #adapter_seq = None if trim_idx is None else sequence[:trim_idx]
        #transcript_seq = sequence[trim_idx:]
    return transcript_seq, adapter_seq, score, trim_idx
    
        
        
        
def score_umi(putative_umi: str,
              umi_dict: dict,
              max_score: int = 2) -> dict:
    '''
    Provided a dict with {name:umi}, return the scores per umi sorted by rank
    '''
    scores = {}
    if putative_umi is None:
        scores = dict([(name,None) for name in umi_dict])
        for k in ['Best UMI', 'Best UMI Score', 'Delta UMI Score']:
            scores[k] = None
    else:
        for name, umi in umi_dict.items():
            umi_align = align(rc(umi), putative_umi, max_score=max_score, score_only=False) # Change max_score to 2
            scores[name] = umi_align['editDistance'] #umi_score
        
        select_scores = [i for i in scores.items() if i[1] is not None]
        if len(select_scores) == 0:
            best_umi, best_umi_score, delta_umi_score = [None]*3
        sorted_scores = sorted(select_scores, key=lambda x:x[1], reverse=True)
        best_umi, best_umi_score = sorted_scores[0]
        if len(sorted_scores) > 1:
            delta_umi_score = best_umi_score - sorted_scores[1][1]
        else:
            delta_umi_score = None
        scores['Best UMI'] = best_umi
        scores['Best UMI Score'] = best_umi_score
        scores['Delta UMI Score'] = delta_umi_score

    return scores


def parse_adapter3_seq(sequence: str,
                       primer: str,
                       umi_dict: dict,
                       max_align_score: int = 2): # Restrict max edits to 2 for adapter
    '''
    Identify sequence characteristics in the 3' adapter sequence, including any potential UMIs
    '''
    if sequence is None:
        primer_score = None
        putative_umi = None
    else:
        primer_alignment = align(primer, sequence, max_score=max_align_score)
        primer_locations = primer_alignment['locations']
        #primer_score = primer_alignment['editDistance']
        #print( sequence, primer, primer_alignment )
        if len(primer_locations) > 0:
            primer_score = primer_alignment['editDistance']  #(len(primer)-primer_alignment['editDistance'])/len(primer)
            post_start = np.min([l[0] for l in primer_locations])
            putative_umi = sequence[:post_start]
        else:
            primer_score = None #0.0 #0.5/len(primer)
            #post_start = len(sequence) # Retain full sequence for UMI alignment
            putative_umi = None
    umi_scores = score_umi(putative_umi, umi_dict)
    return primer_score, umi_scores
    

def compute_read_statistics(read: dict):
    read_length = len(read['sequence'])
    mean_quality = np.mean([ord(x) for x in read['QUAL']])-33.0
    read_stats = {'read_length':read_length,
                  'mean_quality':mean_quality}
    return read_stats


def annotate_read(read: dict, 
                  ssp: str, 
                  primer3: str, 
                  umis: dict,
                  min_transcript_len: int = 200):
    '''
    Search a read for a PolyA tail and then score the 5'/3' ends accordingly
    
    *** POSSIBLE PATCH FOR LATER ***
    Currently max_polyA_score_len is arbitrary, consider changing later
    ***
    
    '''

    ## Read sequence, identify polyA tail if present, process 3' and 5' ends
    sequence = read['sequence']
    read_stats = compute_read_statistics(read)
    pre_polyA, post_polyA, polyA_stats = split_sequence_by_polyA(sequence)
    primer_3_score, umi_scores = parse_adapter3_seq(post_polyA, primer3, umis)
    transcript_seq, adapter5_seq, primer_5_score, primer5_trim_idx = parse_adapter5_seq(pre_polyA, ssp)
    
    ## Analyze the identified transcript sequence
    if transcript_seq is not None:
        transcript_len = len(transcript_seq)
        if primer5_trim_idx is not None:
            transcript_qual = read['QUAL'][primer5_trim_idx:primer5_trim_idx+transcript_len]
        else:
            transcript_qual = None 
    else:
        transcript_len = None
        transcript_qual = None



    t_scores = np.array([primer_5_score, primer_3_score, umi_scores['Best UMI Score']])
    if np.all(t_scores == None): # Garbage read
        cDNA_status = 'Unknown'
    else: # Something was detected
        if np.all(t_scores != None): # Everyone is good :D
            cDNA_status = 'Complete'
        else: # More complicated
            if np.all(t_scores[1:] != None): # 3' end intact with barcode - good for quant
                cDNA_status = '3\' Only'
            elif np.all(t_scores[:2] != None): # Supposedly complete with no barcode
                cDNA_status = 'Missing Barcode'
            elif np.all(t_scores[1:] == None): # Only 5' adapter detected
                cDNA_status = '5\' Only'
            else: # Not sure what this would be...
                cDNA_status = 'Unknown' 
    if transcript_len is not None and transcript_len < min_transcript_len:
        cDNA_status += ' Fragment'
    '''
        if np.all(end_scores != None):
            cDNA_status = 'Complete'
        else:
            cDNA_status = '5\' Only' if end_scores[1] == None else '3\' Only'
        if transcript_len < min_transcript_len:
            cDNA_status += ' Fragment'

    if primer_3_score is not None and primer_5_score is not None:
        cDNA_status = 'Complete'
    else:
        if primer_3_score is not None:
            cDNA_status = '3\' Only'
        elif primer_5_score is not None:
            cDNA_status = '5\' Only'
        else:
            cDNA_status = 'Unknown'
    '''

    ## Collate the annotations
    read_annotations = { '3\' primer':                 primer3,
                         '3\' primer alignment score': primer_3_score,
                         '3\' adapter sequence':       post_polyA,
                         'Strand switch primer':       ssp,
                         'SSP alignment score':        primer_5_score,
                         '5\' adapter trim':           primer5_trim_idx,
                         '5\' adapter sequence':       adapter5_seq,
                         'Transcript sequence':        transcript_seq,
                         'Transcript quality':         transcript_qual,
                         'Transcript length':          transcript_len,
                         'cDNA status':                cDNA_status
                       }
    
    return read | read_stats | polyA_stats | umi_scores | read_annotations
    

    '''
    if polyA_stats is not None: #post_polyA:
        # Potential polyA tail found, verify its reasonable by finding 3' primer
        primer_3_score, umi_scores = parse_adapter3_seq(post_polyA, primer3, umis)
        if primer_3_score is not None:
            ## Something?
        else:
            ## Something else?

        # Score potential UMIs in 3' adapter sequence
        primer_3_score, umi_scores = parse_adapter3_seq(post_polyA, primer3, umis)
        polyA_score = polyA_stats['PolyA length'] #np.min([polyA['PolyA length']/max_polyA_score_len, 1.0])
        polyA_start_idx, polyA_end_idx = polyA['PolyA start'], polyA['PolyA end']
        polyA_length = polyA['PolyA length']
    else:
        primer_3_score = 0.0 #0.5/len(primer3)
        umi_scores = 0.0 #dict([(umi,0.5/len(umi)) for umi in umis])
        polyA_score = 0.0 #0.5/max_polyA_score_len
        polyA_start_idx, polyA_end_idx, polyA_length = [None]*3
    transcript_seq, adapter5_seq, primer_5_score, primer5_trim_idx = parse_adapter5_seq(pre_polyA, ssp, 2) # HARD CODED MAX SCORE ATM, CHANGE LATER - Adjust max to 2 edits from 10 for greater stringency
    transcript_len = len(transcript_seq)
    #transcript_qual = read['QUAL'][primer5_trim_idx:polyA_start_idx]
    transcript_qual = read['QUAL'][primer5_trim_idx:primer5_trim_idx+transcript_len]


    read_annotations = { '3\' primer':                 primer3,
                         '3\' primer alignment score': primer_3_score,
                         '3\' adapter sequence':       post_polyA,
                         'Strand switch primer':       ssp,
                         'SSP alignment score':        primer_5_score,
                         '5\' adapter trim':           primer5_trim_idx,
                         '5\' adapter sequence':       adapter5_seq,
                         #'UMI Scores':                 umi_scores,
                         #'PolyA coordinates':          polyA,
                         'PolyA score':                polyA_score,
                         'PolyA start':                polyA_start_idx,
                         'PolyA end':                  polyA_end_idx,
                         'PolyA length':               polyA_length,
                         'Transcript sequence':        transcript_seq,
                         'Transcript quality':         transcript_qual,
                         'Transcript length':          transcript_len
                       }
    return read | read_stats | read_annotations | polyA | umi_scores
    '''
    
def process_read(read: dict,
                 ssp: str,
                 primer3: str,
                 umis: dict,
                 min_polyA_score: float = 0.0,
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
    
    # Analyze read in both forward and reverse directions
    f_analysis = annotate_read(read, ssp, primer3, umis)
    rc_read = dict(read)
    rc_read['sequence'] = rc(rc_read['sequence'])
    r_analysis = annotate_read(rc_read, ssp, primer3, umis)
    
    # Select forward or reverse direction based on analysis
    cDNA = np.array([a['cDNA status'] for a in [f_analysis,r_analysis]])
    if np.all((cDNA == 'Unknown') | (cDNA == 'Unknown Fragment')):
        analysis = f_analysis # Default to forward read, since we're guessing
    else:
        if np.any((cDNA == 'Unknown') | (cDNA == 'Unknown Fragment')): # Only one good option
            analysis = f_analysis if cDNA[1] == 'Unknown' else r_analysis
        else: # Multiple options
            if np.any(cDNA == 'Complete'):
                # Assume only one possible good read here
                analysis = f_analysis if cDNA[0] == 'Complete' else r_analysis
            else:
                # Return forward analysis with a modified label
                analysis = f_analysis
                analysis['cDNA status'] = 'Complex' #f"{f_analysis['cDNA status']}_{r_analysis['cDNA status']}"

    return analysis
        




    '''
    #if np.max([list(forward_analysis['UMI Scores'].values())[0],list(reverse_analysis['UMI Scores'].values())[0]]) > min_umi_score:
    #if np.max([d[umi] for d in [forward_analysis, reverse_analysis] for umi in umis]) > min_umi_score:
    if np.max([a['PolyA score'] for a in [for_analysis, rev_analysis]]) > min_polyA_score:
        # Found a PolyA tail
        delta_polyA_score = for_analysis['PolyA score'] - rev_analysis['PolyA score']
        if delta_polyA_score > 0: #max([for_analysis['PolyA score'], rev_analysis['PolyA score']]) > min_polyA_score:
            analysis = for_analysis if delta_polyA_score > 0 else rev_analysis
            annotation = 'Full-length cDNA' if analysis['SSP alignment score'] > min_ssp_score else "3\' only"
        else:
            analysis = for_analysis # Need to actually solve for the correct one
            annotation = "Duplexed 3\' only"
        umi_scores = [(umi, analysis[umi]) for umi in umis]
        umi_dict = dict(sorted(umi_scores, key=lambda x:x[1], reverse=True))
        best_umi, second_umi = list(umi_dict)[:2]
        best_umi_score = umi_dict[best_umi]
        delta_umi_score = best_umi_score - umi_dict[second_umi]
        
    else:
        # Did not find a PolyA tail, check for 5' features
        best_umi, best_umi_score, delta_umi_score = [None]*3
        max_ssp_score = np.max([for_analysis['SSP alignment score'],rev_analysis['SSP alignment score']])
        if max_ssp_score > min_ssp_score:
            delta_ssp_score = for_analysis['SSP alignment score']-rev_analysis['SSP alignment score']
            if np.abs(delta_ssp_score) > min_ssp_score:
                analysis = for_analysis if delta_ssp_score > 0 else rev_analysis
                annotation = "5\' only"
            else:
                analysis = for_analysis # Need to solve for which one is better
                annotation = "Duplexed 5\' only"
        else:
            analysis = for_analysis # Unknown, no better guess
            annotation = "Uncharacterized sequence"
    analysis['annotation'] = annotation 
    analysis['Best UMI'] = best_umi
    analysis['Best UMI Score'] = best_umi_score
    analysis['Delta UMI Score'] = delta_umi_score
    return analysis #, forward_analysis, reverse_analysis
    '''

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