'''

Tools for performing and manipulating sequence alignments

'''

from skbio.alignment import StripedSmithWaterman
from .utils import reverse_complement

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