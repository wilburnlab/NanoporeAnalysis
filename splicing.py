import re
import copy
import subprocess
import sys
import time
from pathlib import Path
import concurrent
import concurrent.futures
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq
from pyarrow import dataset
from pyarrow import csv
from scipy import stats
from scipy import optimize
from scipy.optimize import differential_evolution
import numpy as np
if '/users/PAS1669/smith12380/py' not in sys.path :
    sys.path.append('/users/PAS1669/smith12380/py')
from NanoporeAnalysis import utils
from NanoporeAnalysis import local_io
from NanoporeAnalysis import analysis
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar

def parse_cigar(cig, ignore_I = False, minimum_region_length = 0 ) :
    """
    Parses the CIGAR string of a sequence alignment into a list of [] items.
    
    Args :
        cig (str) : the CIGAR string to be parsed. Should only contain letters A-Z and numbers.
        ignore_I (bool) : whether or not to ignore the letter I, which represents inserts in the alignment. If True, inserts will simply be ignored.
        minimum_region_length (int) : the minimum length of items to return. If greater than 0, smaller items will be added to the preceeding item.
            This can be useful if only looking for larger features, as reducing the total number of items significantly speeds up future processing.
    
    Returns :
        cig_parsed (list of [length (int), value (char)] items) : the parsed list of items.
    """
    cig_split = re.findall("[0-9]+|[a-zA-Z]+", cig)
    cig_numbers = cig_split[0::2]
    cig_letters = cig_split[1::2]
    cig_parsed = []
    previous_letter = None
    for number_str, letter in zip(cig_numbers, cig_letters) :
        if letter == 'I' :
            if not ignore_I :
                if letter == previous_letter :
                    cig_parsed[-1][0] += int(number_str)
                elif int(number_str) >= minimum_region_length :
                    cig_parsed.append([ int(number_str), letter ])
                    previous_letter = letter
        elif ( int(number_str) < minimum_region_length or letter == previous_letter ) and previous_letter != None :
            cig_parsed[-1][0] += int(number_str)
        else :
            cig_parsed.append([ int(number_str), letter ])
            previous_letter = letter
    return cig_parsed

def clip_cigar(cig, clip) :
    """
    This is a utility function for the next function. Clips the start of a parsed CIGAR by the given number of nucleotides.
        Ignores any insert regions, as this is meant to be used while aligning multiple alignments to a common template.
        Can be used with non-cigar inputs that also have the same [length,...] format.
        
    Args :
        cig (list of [length (int), value (char)] items) : the parsed CIGAR. Must be parsed, cannot be the raw string.
        clip (int) : length of cig to remove from the start of the cig.
        
    Returns :
        cig_copy (list of [length (int), value (char)] items) : the clipped cig.
    """
    if type(cig) == str :
        raise TypeError("Error: this cig doesn't appear to be parsed. Must be a parsed cigar in the form (list of [length (int), value (char)] items).")
    clip_index = 0
    cig_copy = copy.deepcopy(cig)
    while clip_index < clip :
        clip_remaining = clip - clip_index
        if len(cig_copy) == 0 :
            return None
        elif cig_copy[0][1] == 'I' :
            del cig_copy[0]
        elif cig_copy[0][0] > clip_remaining :
            clip_index += clip_remaining
            cig_copy[0][0] -= clip_remaining
        elif cig_copy[0][0] <= clip_remaining :
            clip_index += cig_copy[0][0]
            del cig_copy[0]
    return cig_copy

def clip_cigar_ends(cig, left_clip, right_clip, return_clipped_seq_lengths = False) :
    """
    Clips the ends of a parsed CIGAR sequence by the prescribed lengths. Ignores insert regions.

    Args :
        cig (list of [length (int), value (char)] items) : the parsed CIGAR. Must be parsed, cannot be the raw string.
        left_clip (int) : length of cig to remove from the left end.
        right_clip (int) : length of cig to remove from the right end.
        return_clipped_seq_lengths (bool) : whether or not to return the length of the cig's underlying sequence that is clipped. This is used to re-align the
            sequence to the clipped cig.

    Returns :
        cig_parsed (list of [length, value] pairs) : the clipped CIGAR sequence.
    """
    if type(cig) == str :
        raise TypeError("Error: this cig doesn't appear to be parsed. Must be a parsed cigar in the form (list of [length (int), value (char)] items).")
    left_clipped_length = 0
    right_clipped_length = 0
    if left_clip != 0 :
        clipped_cig = clip_cigar( cig, left_clip )
        if return_clipped_seq_lengths :
            left_clipped_length = get_total_cigar_length(cig, chars_to_ignore=['N', 'D', ' ']) - get_total_cigar_length(clipped_cig, chars_to_ignore=['N', 'D', ' '])
        cig = clipped_cig
    if cig != None and right_clip != 0 :
        cig.reverse()
        clipped_cig = clip_cigar( cig, right_clip )
        if return_clipped_seq_lengths :
            right_clipped_length = get_total_cigar_length(cig, chars_to_ignore=['N', 'D', ' ']) - get_total_cigar_length(clipped_cig, chars_to_ignore=['N', 'D', ' '])
        cig = clipped_cig
        if cig != None :
            cig.reverse()
    if return_clipped_seq_lengths :
        return cig, left_clipped_length, right_clipped_length
    else :
        return cig

def pad_cigar(cig, r_st, r_en, pileup_indices, index_0, ignore_ends = False) :
    """
    Pads parsed CIGAR sequence by adding 'blank' items onto the ends. Does so based on given r_st (start) and r_en (end)
        values that denote the point along the common reference sequence that each CIGAR sequence starts and ends at, respectively, and
        pileup_indices, which denotes the minimum and maximum r_st and r_en of all cigs being aligned.
        This is meant to be used with process_cigar_strings, which takes care of aligning CIGARs.

    Args :
        cig (list of [length, value] pairs) : the parsed CIGAR. Must be parsed, cannot be the raw string.
        r_st (int) : the 0-indexed point along the reference sequence that the CIGAR starts at.
        r_en (int) : the 0-indexed point (exclusive) along the reference sequence that the CIGAR ends at.
        pileup_indices (list of [int, int]) : the minimum r_st and maximum r_en of all cigs being aligned.
        index_0 (bool) : whether or not to pad the cig out to 0 on the left end.
        ignore_ends (bool) : whether or not to remove the first and last items of the cig.

    Returns :
        cig (list of [length, value] pairs) : the padded CIGAR.
    """
    if index_0 :
        left_gap = r_st
    else :
        left_gap = max( r_st - pileup_indices[0], 0 )
    right_gap = max( pileup_indices[1] - r_en, 0 )
    if ignore_ends :
        cig[0][1] = ' '
        cig[-1][1] = ' '
    if left_gap > 0 :
        cig.insert(0, [left_gap, ' '])
    if right_gap > 0 :
        cig.append([right_gap, ' '])
    return cig

def process_cigar_strings(cigs, r_sts, r_ens, index_0 = False, left_clip = 0, right_clip = 0, minimum_region_length = 0, ignore_I = False, return_null_cigs = True, ignore_ends = False, return_passed_cigs_indices = False) :
    """
    Process CIGAR strings and align them together. These should all be CIGARs from alignments to a common reference sequence.

    Args :
        cigs (list of str) : the CIGAR strings. Refer to parse_cigar() for more information.
        r_sts (list of int) : the 0-indexed point along the reference sequence that each CIGAR starts at.
        r_en (list of int) : the 0-indexed point (exclusive) along the reference sequence that each CIGAR ends at.
        index_0 (bool) : whether or not to pad the cig out to 0 on the left end.
        left_clip (int) : length of cig to remove from the left end.
        right_clip (int) : length of cig to remove from the right end.
        minimum_region_length (int) : the minimum length of items to return. If greater than 0, smaller items will be added to the preceeding item.
            This can be useful if only looking for larger features, as reducing the total number of items significantly speeds up future processing.
        ignore_I (bool) : whether or not to ignore the letter I, which represents inserts in the alignment. If True, inserts will simply be ignored.
        return_null_cigs (bool) : whether or not to return cigs that get evaluated to a blank sequence at some point in the process,
            usually due to clipping or some error in the CIGAR string. This can be set to True when the returned cigs need to be matched with some other list.
        ignore_ends (bool) : whether or not to remove the first and last items of the cig.

    Returns :
        cigs_gapped (list of cig (list of [length (int), value (char)] items)) : the processed and aligned cigs.
    """
    cigs_parsed = [ parse_cigar( cig, ignore_I = ignore_I, minimum_region_length = minimum_region_length ) for cig in cigs ]
    r_sts_passed = [ x[0] for x in zip(r_sts, cigs_parsed) if x[1] != None]
    r_ens_passed = [ x[0] for x in zip(r_ens, cigs_parsed) if x[1] != None]
    if len(r_sts_passed) == 0 :
        return []
    pileup_indices = [ min(r_sts_passed), max(r_ens_passed) ]
    cigs_gapped = []
    indices = range(len(cigs_parsed))
    indices_passed = []
    for cig_parsed, r_st, r_en, i in zip(cigs_parsed, r_sts, r_ens, indices) :
        if cig_parsed != None :
            cig_parsed = pad_cigar(cig_parsed, r_st, r_en, pileup_indices, index_0, ignore_ends)
            cig_parsed = clip_cigar_ends(cig_parsed, left_clip, right_clip)
        if cig_parsed == None :
            if return_null_cigs :
                cigs_gapped.append(cig_parsed)
        else :
            cigs_gapped.append(cig_parsed)
            indices_passed.append(i)
    if return_passed_cigs_indices :
        return cigs_gapped, indices_passed
    else :
        return cigs_gapped

def get_total_cigar_length(cig, chars_to_ignore = []) :
    """
    Evaluates the length of the parsed cig while potentially ignoring inserts.

    Args:
        cig (list of [length (int), value (char)] items) : parsed cigar.
        chars_to_ignore (list of cigar values (usually characters)) : which values of the cigar to ignore. Useful to exclude things like inserts.
        
    Returns :
        total (int) : length of the cig.
    """
    total = 0
    for item in cig :
        if item[1] not in chars_to_ignore:
            total += item[0]
    return total

def cigar_to_str(cig, chars_to_ignore = []) :
    """
    Converts a parsed CIGAR sequence into a fully expanded CIGAR string.

    Args:
        cig (list of [length (int), value (char)] items) : parsed cigar.
        chars_to_ignore (list of char) : any characters/values to be ignored in the cigar.
        
    Returns :
        cig_str (str) : converted CIGAR string.
    """
    cig_parsed = parse_cigar(cig) if type(cig) != list else cig
    cig_str = ''
    for item in cig_parsed :
        if item[1] not in chars_to_ignore :
            cig_str += item[0] * item[1]
    return cig_str

def compress_cigar_str(cigar_str) :
    """
    Compresses an expanded alignment string into a CIGAR-type sequence. Not currently used, but this is a good utility function.

    Args :
        cigar_str (str) : a string of characters, each representing a matching, missing, inserted, or mismatched position along an alignment of strings.
        
    Returns :
        compressed_cig (list of [length (int), value (char)] items) : the resulting CIGAR.
    """
    compressed_cig = []
    current_value = 0
    for i, position in enumerate(cigar_str) :
        if position != current_value :
            compressed_cig.append([ 1, position ])
            current_value = position
        else :
            compressed_cig[ len(compressed_cig) - 1 ][0] += 1
    return compressed_cig

def return_to_cigar(cig, chars_to_ignore = []) :
    """
    Converts a processed cigar back into a CIGAR format string.

    Args :
        cig (list of [length (int), value (char)] items) : a list of (length, value) items representing a sequence alignment, ie the result of parse_cigar().
        
    Returns :
        cig_str (str) : CIGAR format of the cig.
    """
    cig_str = ''
    for item in cig :
        if item[1] not in chars_to_ignore :
            cig_str += str(item[0])
            cig_str += str(item[1])
    return cig_str

def pad_compact_seqs(compact_seqs, pad_value) :
    """
    Adds 'blank' items to the right end of compact seqs (lists of [length, values...]) to create equal-length seqs.

    Args :
        compact_seqs (list[ list[length, values...] ]) : the compact seqs (CIGAR-style packaging of sequence/alignment data).
        pad_value (list[values...]) : the value to be added to the right end of the compact seqs. This should be the same format as the compact_seq items
            but with the length value removed, and it should represent 'blank' space (ie likely not 'D' for a deletion, since there is no underlying sequence
            alignment beyond the cigar)
            
    Returns :
        compact_seqs (list[ list[length, values...] ]) : the padded compact seqs.
    """
    lengths = [ get_total_cigar_length(x) for x in compact_seqs ]
    max_length = max(lengths)
    for compact_seq, length in zip (compact_seqs, lengths) :
        if length < max_length :
            compact_seq.append([max_length - length] + pad_value)
    return compact_seqs

def align_compact_seqs(compact_seqs, return_numpy = True) :
    """
    Splits the items within compact seqs (lists of [length, values...]) such that each seq has the same number of items and matching length values.
        i.e. turns [[2,3,1], [5,2,1]] and [[7,2,1]] into [[2,3,1], [5,2,1]] and [[2,2,1], [5,2,1]]. 
        Note that the input seqs MUST have the same total length (sum of the first values of each compact seq item).

    Args :
        compact_seqs (list[ list[length, values...] ]) : the compact seqs (CIGAR-style packaging of sequence/alignment data). MUST all be the same total length.
        return_numpy (bool) : whether or not to return the seqs as a numpy array, which is mainly a problem if your seqs have non-numerical data types. If
            False, the resulting numpy array will all be of one type (ie if there are strings in the seqs, all values (including lengths) will be strings).
        
    Returns :
        aligned_seqs (list of [length, values...]) : the aligned seqs. Will be a numpy array if return_numpy = True.
    """
    split_indices = [ np.cumsum( [x[0] for x in compact_seq] ) for compact_seq in compact_seqs ]
    unique_splits = np.unique(np.concatenate(split_indices)).tolist() + [0]
    unique_splits.sort()
    diffs = np.diff(unique_splits)
    sort_indices = np.array([ compact_seq_indices.searchsorted(unique_splits, side='right')[:-1] for compact_seq_indices in split_indices ])
    aligned_seqs = np.array([ np.insert(np.take_along_axis(np.array(compact_seq)[:,1:], np.expand_dims(indices, 1), axis=0), 0, diffs, axis=1) for compact_seq, indices in zip(compact_seqs, sort_indices) ])
    if not return_numpy :
        lengths_arrays = aligned_seqs[:, :, :1].astype(np.int64()).tolist()
        content_arrays = aligned_seqs[:, :, 1:].tolist()
        aligned_seqs = [ [ length + content for length, content in zip(lengths_array, content_array) ] for lengths_array, content_array in zip(lengths_arrays, content_arrays) ]
    return aligned_seqs

def merge_pileups(pileups) :
    """
    Aligns and adds together two pileups.

    Args :
        pileups (list[ list[length (int), num aligned (int), num matching (int)] ]) : the pileups to merge. Look at cigs_to_pileup() for more info on pileups.
    
    Returns :
        pileup (list[length (int), num aligned (int), num matching (int)]) : the merged pileup.
    """
    if len(pileups) == 0 :
        return None
    elif len(pileups) == 1 :
        return pileups[0]
    else :
        pileups_padded = pad_compact_seqs(pileups, [0,0])
        aligned_pileups = align_compact_seqs(pileups_padded, return_numpy = True)
        merged = np.concatenate([aligned_pileups[0,:,:1], np.sum(aligned_pileups[:,:,1:], axis=0)], axis=1).tolist()
        return merged

def cigs_to_pileup(cigs, r_sts, r_ens, minimum_region_length = 5, ignore_ends = False) :
    """
    Converts CIGAR strings into an alignment pileup to find the total occupancy along the reference sequence. Uses the r_sts and r_ens values to
        align the CIGAR sequnces. See the functions process_cigar_strings() and align_compact_seqs() for more details. Treats M in the CIGAR as an aligned
        match, and treats N, D as gaps, and ignores I (inserts make pileups very complicated). Meant to be used with CIGARs that refer to a common reference.

    Args : 
        cigs (list of str) : the CIGAR strings to be piled up.
        r_sts (list of int) : the 0-indexed point along the reference sequence that each CIGAR starts at.
        r_en (int) : the 0-indexed point (exclusive) along the reference sequence that each CIGAR ends at.
        minimum_region_length (int) : the minimum length of items to return. If greater than 0, smaller items will be added to the preceeding item.
            This can be useful if only looking for larger features, as reducing the total number of items significantly speeds up future processing.
        ignore_ends (bool) : whether or not to remove the first and last items of the cig.

    Returns:
        pileup (list of [length, num aligned CIGARs, num matching CIGARs]) : the resulting pileup in a compact seq format, where each item represents a region
            of the common reference sequence. The first value of each item is the length of that region, the second value is the number of CIGARs that
            covered that region, and the third is the number that matched that region.
    """
    cigs_grouped = [ cigs[ i : i+100 ] for i in range(0, len(cigs), 100) ]
    r_sts_grouped = [ r_sts[ i : i+100 ] for i in range(0, len(r_sts), 100) ]
    r_ens_grouped = [ r_ens[ i : i+100 ] for i in range(0, len(r_ens), 100) ]
    pileups = []
    for cigs_group, r_st_group, r_en_group in zip( cigs_grouped, r_sts_grouped, r_ens_grouped ) :
        cigs_parsed = process_cigar_strings( cigs_group, r_st_group, r_en_group, index_0 = True, ignore_I = True, minimum_region_length = minimum_region_length, return_null_cigs = False, ignore_ends = ignore_ends )
        if len(cigs_parsed) != 0 :
            cigs_aligned = align_compact_seqs(cigs_parsed)
            Ms = np.count_nonzero(cigs_aligned[:,:,1] == 'M', axis=0)
            NDs = np.count_nonzero((cigs_aligned[:,:,1] == 'N') | (cigs_aligned[:,:,1] == 'D'), axis=0)
            pileup = np.insert( np.stack( [cigs_aligned[0,:,0].astype('<i4'), Ms + NDs], axis=1), 2, Ms, axis=1 ).tolist()
            pileups.append(pileup)
    del cigs_grouped
    del r_sts_grouped
    del r_ens_grouped
    if len(pileups) == 0 :
        return [None]
    if len(pileups) > 1 :
        merged = merge_pileups(copy.deepcopy(pileups))
        del pileups
        return merged
    else :
        return pileups[0]

def remove_outlier_alignments(table) :
    """
    Removes sequence alignments from a pyarrow table that are distant outliers. Uses a naive method of making histogram bins and removing the lower and upper
        bins that include less than 5% of values.

    Args :
        table (pyarrow table) : the table with alignments.
        
    Returns :
        table (pyarrow table) : the table with outliers removed (if there are any).
    """
    num_alignments = table.num_rows
    r_sts = [ int(x) for x in table.column('minimap2_r_st').to_pylist() ]
    r_st_counts, r_st_bins = np.histogram( r_sts, bins='doane')
    r_st_bins_and_counts = [ x for x in zip(r_st_counts, r_st_bins)]
    r_st_bins_and_counts.sort(key = lambda x : x[1])
    r_st_bins_inclusions = [ sum([ x[0] for x in r_st_bins_and_counts[i:] ]) / num_alignments for i in range(len(r_st_bins_and_counts)) ]
    r_st_bins_valid = [ x for x, y in zip(r_st_bins_and_counts, r_st_bins_inclusions) if y >= 0.95 ]
    r_st_bound = max( r_st_bins_valid, key = lambda x : x[0])[1]
    r_ens = [ int(x) for x in table.column('minimap2_r_en').to_pylist() ]
    r_en_counts, r_en_bins = np.histogram( r_ens, bins='doane')
    r_en_bins_and_counts = [ x for x in zip(r_en_counts, r_en_bins[1:])]
    r_en_bins_and_counts.sort(key = lambda x : x[1], reverse=True)
    r_en_bins_inclusions = [ sum([ x[0] for x in r_en_bins_and_counts[i:] ]) / num_alignments for i in range(len(r_en_bins_and_counts)) ]
    r_en_bins_valid = [ x for x, y in zip(r_en_bins_and_counts, r_en_bins_inclusions) if y >= 0.95 ]
    r_en_bound = max( r_en_bins_valid, key = lambda x : x[0])[1]
    return table.filter( (pc.field('minimap2_r_st').cast(pa.float64()) >= r_st_bound) & (pc.field('minimap2_r_en').cast(pa.float64()) <= r_en_bound) )

def table_to_pileup_by_gene(table, minimum_region_length = 5, max_match = 1, ignore_ends = False, quiet = True, minimum_num_reads = 5) :
    """
    Creates a pileup per gene for all reads in a table. See the function cigs_to_pileup() for more details on the pileups.

    Args :
        table (pyarrow table) : the table with alignments.
        minimum_region_length (int) : the minimum length of items to return. If greater than 0, smaller items will be added to the preceeding item.
            This can be useful if only looking for larger features, as reducing the total number of items significantly speeds up future processing.
        max_match (float between 0 and 1) : the maximum allowed fraction of aligned reference sequence that is matched by each read. This is useful to exclude reads that
            are completely or nearly matched to the reference, which can be bad alignments or weird transcripts, if you're working with mRNA that is supposed
            to have introns. Defaults to 1, which does not filter anything out.
        ignore_ends (bool) : whether or not to remove the first and last items of the cig.
        quiet (bool) : whether to print the progress of this function as it goes.
        minimum_num_reads (int) : minimum number of reads required per gene to create a pileup. This speeds things up a bit by skipping low-count genes
            that you are going to want to ignore in later steps anyways.

    Returns : 
        pileup_dict (dict) : a dictionary of gene:pileup pairings. Genes with too few reads are given a pileup of [None].
    """
    pileup_dict = {}
    genes = table.column('minimap2_gene_ID').unique().to_pylist()
    for i, gene in enumerate(genes) :
        gene_table = table.filter( pc.field('minimap2_gene_ID') == gene )
        if max_match < 1 :
            gene_table = gene_table.filter(pc.less_equal( pc.divide( pc.field('minimap2_mlen').cast(pa.float32()), pc.abs( pc.subtract( pc.field('minimap2_r_en'), pc.field('minimap2_r_st') ) ) ), pc.scalar(max_match) ))
        if gene_table.num_rows >= minimum_num_reads :
            gene_table = remove_outlier_alignments(gene_table)
            if gene_table.num_rows >= minimum_num_reads :
                ctgs = gene_table.column('minimap2_ctg')
                if len(ctgs.unique().to_pylist()) > 1 :
                    ctg = pa.table( gene_table.column('minimap2_ctg').value_counts() ).sort_by( [('counts', 'descending')] ).column('values').to_pylist()[0]
                    gene_table = gene_table.filter( pc.field('minimap2_ctg') == ctg )
                cigs = gene_table.column('minimap2_cigar').to_pylist()
                r_sts = [ int(x) for x in gene_table.column('minimap2_r_st').to_pylist() ]
                r_ens = [ int(x) for x in gene_table.column('minimap2_r_en').to_pylist() ]
                pileup_dict[gene] = [cigs_to_pileup( cigs, r_sts, r_ens, minimum_region_length, ignore_ends )]
            else :
                pileup_dict[gene] = [None]
        else :
            pileup_dict[gene] = [None]
        if not quiet :
            if i % 100 == 0 and i != 0 :
                print(i, len(genes))
    return pileup_dict

def create_pileups(path_dataset, path_out_pileup = None, minimum_region_length = 15, workers = 4, max_match = 1, ignore_ends = False, quiet = True, minimum_num_reads = 5) :
    """
    Creates sequence depth pileups against the genomic reference per gene for each barcode, creates a table with each gene as a column, and a row per barcode
        with the sample ID annotated, in the 'barcode_ID' and 'sample' columns, respectively. Pileups are a list of compressed, cigar-like regions in the form
        of [ number of bases long, number of reads spanning the region, number of reads matching the region. ]. The overall process goes as follows:
        gather all reads for one gene for one barcode -> clip CIGARS to the region -> align all their CIGAR sequences together using their alignments to the
        reference genome while keeping them collapsed as [ length of region, data for region ] items -> go through all CIGARs one region at a time and
        accumulate them into the depth pileups -> repeat for each gene, for each barcode, and accumulate into table. By virtue of how minimap creates the
        CIGARs during alignment, they always begin and end with bases matched to the reference, but will 'ignore' small mismatches within the alignment,
        instead marking them as matches. Indels are *never* ignored, but are sometimes labels as N in the CIGAR instead of an indel... IDK why. This function
        is fairly well validated, I have ironed out kinks and it consistently aligns and makes pileups that are accurate to the actual data. 
    
    Args:
        path_datase (str, or list of str) : the path(s) to the dataset. Can be given a list of paths if the data isn't all in one.
        path_out_pileup (str) : path to output the pileup table, if desired. Can be left blank to not write to file.
        minimum_region_length (int) : minimum region length to allow within the pileups. Anything smaller than this will be ignored and added to the preceeding region, if it's not an indel, or just cut out, if its an indel (doing so preserves inter-cigar alignment). This is only done to limit the computational load, which is primarily in the handling of the cigar regions. Fewer regions makes for much faster work and smaller memory needs.
        workers (int) : number of parallel tasks to run. I've been using as many as possible, as this is a massive process on the entire dataset.
            This also doesn't use a ton of memory, since it's working with CIGAR strings.
        max_match (float between 0 and 1) : the maximum allowed fraction of aligned reference sequence that is matched by each read. This is useful to exclude reads that
            are completely or nearly matched to the reference, which can be bad alignments or weird transcripts, if you're working with mRNA that is supposed
            to have introns. Defaults to 1, which does not filter anything out.
        ignore_ends (bool) : whether or not to remove the first and last items of the cig.
        quiet (bool) : whether to print the progress of this function as it goes.
        minimum_num_reads (int) : minimum number of reads required per gene to create a pileup. This speeds things up a bit by skipping low-count genes
            that you are going to want to ignore in later steps anyways.
        
    Returns None. Outputs to path_out_pileup.
    """
    print("loading dataset")
    if type(path_dataset) == list :
        ds = dataset.dataset([ dataset.dataset(x) for x in path_dataset ])
    else :
        ds = dataset.dataset(path_dataset)
    print("cleaning dataset")
    ds = ds.filter(~pc.field('barcode_sample_ID').isin(pa.array([None], type=pa.string())))
    sample_IDs = ds.to_table(columns=['barcode_sample_ID']).column('barcode_sample_ID').unique().to_pylist()
    sample_ID_list = []
    barcode_ID_list = []
    pileup_table = False
    for sample_ID in sample_IDs :
        print("filtering dataset")
        ds_filtered_sample = ds.filter( pc.field('barcode_sample_ID') == sample_ID ).filter(~pc.field('minimap2_gene_ID').isin(pa.array([None], type=pa.string())))
        sample_barcodes = ds_filtered_sample.to_table(columns=['barcode_ID']).column('barcode_ID').unique().to_pylist()
        for barcode in sample_barcodes :
            print(barcode)
            sample_ID_list.append(sample_ID)
            barcode_ID_list.append(barcode)
            print("loading barcode table")
            barcode_filtered_ds = ds_filtered_sample.filter(pc.field('barcode_ID') == barcode)
            pileup_dict = {}
            genes = barcode_filtered_ds.to_table(columns=['minimap2_gene_ID']).column('minimap2_gene_ID').unique().to_pylist()
            gene_chunks = np.array_split( np.array(genes), workers * 10 )
            futures = []
            print("starting pileup")
            with concurrent.futures.ProcessPoolExecutor( max_workers = workers ) as executor :
                for gene_chunk in gene_chunks :
                    gene_table = barcode_filtered_ds.filter( pc.field('minimap2_gene_ID').isin(gene_chunk) ).scanner(columns = [ 'barcode_ID', 'barcode_sample_ID', 'minimap2_ctg', 'minimap2_cigar','minimap2_r_st','minimap2_r_en', 'minimap2_gene_ID', 'minimap2_mlen' ])
                    futures.append(executor.submit( table_to_pileup_by_gene, gene_table,
                                                   minimum_region_length = minimum_region_length, max_match = max_match,
                                                   ignore_ends = ignore_ends, quiet = quiet, minimum_num_reads = minimum_num_reads ))
                concurrent.futures.wait( futures )
            for future in futures :
                pileup_dict.update(future.result())
            if not pileup_table :
                pileup_table = pa.table(pileup_dict)
            else :
                pileup_table = pa.concat_tables([ pileup_table, pa.table(pileup_dict) ], promote_options="default")
    pileup_table = pileup_table.add_column(0, 'barcode_ID', [barcode_ID_list])
    pileup_table = pileup_table.add_column(0, 'barcode_sample_ID', [sample_ID_list])
    if path_out_pileup != None :
        pq.write_table(pa.table(pileup_table), path_out_pileup)
    return

def get_regions(sample_1_pileups, sample_2_pileups, min_average_read_depth = 20, min_read_depth = 5, min_fraction_change = 0.1) :
    """
    Takes pileups for two samples (should have multiple pileups per sample from replicates), and separates then into distinct regions based on
        sudden changes in the sequence coverage along the pileups. This function aligns the pileups, then goes along them and aggregates each item/sub-region
        in the pileups into a larger region while tracking the global fraction of matched:aligned reads (fraction spliced in). Whenever a item/sub-region is
        reached with a fraction spliced in that differs from that of the current region, it starts a new region. This usually happens at the boundaries of
        introns/exons.

    Args :
        sample_1_pileups (list of lists of [length, num aligned, num matched]) : the pileups for sample 1. See create_pileups() for more details on pileups.
        sample_2_pileups (list of lists of [length, num aligned, num matched]) : the pileups for sample 2. See create_pileups() for more details on pileups.
        min_average_read_depth (int) : the minimum number of average reads aligned or covering each region. Regions that don't meet this will be ignored.
            Averaged across all pileups. Recommend setting to at LEAST 1.
        min_read_depth (int) : the minimum number of reads per pileup aligned or covering each region. Regions that don't meet this will be ignored.
            Not averaged across pileups, so all pileups must meet this minimum. Recommend setting to at LEAST 1.
        min_fraction_change (float between 0 and 1) : The minimum change in the fraction of matched:aligned reads (fraction spliced in) to establish a new region.

    Returns :
        regions (list of [ 
            [region start location, region end location], [
                [ subregion length, [sample 1 num matched, sample 1 num aligned, sample 2 num matched, sample 2 num aligned] ]...
            ] ]) : The resulting regions with all data preserved. These regions are given in a really weird data structure, but I couldn't be bothered to
            make it better. Each region is a list with the boundary indices as a 0-indexed [start, end], then a list of the sub-regions so no data is lost.
    """
    pileups = sample_1_pileups + sample_2_pileups
    pileups = pad_compact_seqs(pileups, [0, 0])
    aligned_pileups = align_compact_seqs(pileups)
    num_pileups_1 = len(sample_1_pileups)
    sample_1_pileups_grouped = [ [ pileup[j] for pileup in aligned_pileups[ : num_pileups_1 ] ] for j in range( len(aligned_pileups[0]) ) ]
    num_pileups_2 = len(sample_2_pileups)
    sample_2_pileups_grouped = [ [ pileup[j] for pileup in aligned_pileups[ num_pileups_1 : ] ] for j in range( len(aligned_pileups[0]) ) ]
    position_index = 1
    last_region_mean_fractions = 0
    last_region_recorded = False
    last_region_length = 0
    regions = []
    for sample_1_set, sample_2_set in zip( sample_1_pileups_grouped, sample_2_pileups_grouped ) :
        sample_1_mean_depth = np.mean([ x[1] for x in sample_1_set ])
        sample_2_mean_depth = np.mean([ x[1] for x in sample_2_set ])
        combined_min_depth = min([ x[1] for x in sample_1_set ] + [ x[1] for x in sample_2_set ])
        if sample_1_mean_depth >= min_average_read_depth and sample_2_mean_depth >= min_average_read_depth and combined_min_depth >= min_read_depth :
            sample_1_spliced_ins = [ max(item[2], 0.1) for item in sample_1_set ]
            sample_1_totals = [ item[1] for item in sample_1_set ]
            sample_2_spliced_ins = [ max(item[2], 0.1) for item in sample_2_set ]
            sample_2_totals = [ item[1] for item in sample_2_set ]
            region_mean_fractions = [ np.mean(np.divide( sample_1_spliced_ins, sample_1_totals)), np.mean(np.divide( sample_2_spliced_ins, sample_2_totals )) ]
            if last_region_recorded != False :
                means_diffs = np.absolute(np.subtract( region_mean_fractions, last_region_mean_fractions ))
                if np.max(means_diffs) >= min_fraction_change and last_region_recorded == True :
                    regions.append([ [ position_index, position_index + sample_1_set[0][0] ], [ [ sample_1_set[0][0], [ sample_1_spliced_ins, sample_1_totals, sample_2_spliced_ins, sample_2_totals ] ] ] ])
                    last_region_mean_fractions = region_mean_fractions
                    last_region_length = sample_1_set[0][0]
                else :
                    regions[-1] = [ [ regions[-1][0][0], position_index + sample_1_set[0][0] ], regions[-1][1] + [ [ sample_1_set[0][0], [ sample_1_spliced_ins, sample_1_totals, sample_2_spliced_ins, sample_2_totals ] ] ] ]
                    last_region_mean_fractions = np.divide( np.add( np.multiply( last_region_mean_fractions, last_region_length ), np.multiply( region_mean_fractions, sample_1_set[0][0] ) ), last_region_length + sample_1_set[0][0] )
                    last_region_length += sample_1_set[0][0]
            else :
                regions.append([ [ position_index, position_index + sample_1_set[0][0] ], [ [ sample_1_set[0][0], [ sample_1_spliced_ins, sample_1_totals, sample_2_spliced_ins, sample_2_totals ] ] ] ])
                last_region_mean_fractions = region_mean_fractions
                last_region_length = sample_1_set[0][0]
            last_region_recorded = True
        else :
            last_region_mean_fractions = 0
            last_region_recorded = False
            last_region_length = 0
        position_index += sample_1_set[0][0]
    return regions

def splicing_region_neg_binom(region, gene = None, min_mean_fraction_spliced_in = 0.05, min_length = 15) :
    """
    Fits a negative binomial model to the data of a region for two samples from get_regions(). Note that sample 1 (the first one given to get_regions()) is
        considered the numerator or dependent variable here, so fraction_change will represent fraction spliced in of sample_2 - sample_1.

    Args :
        region (super nested list) : the data for the region. See get_regions() for details on how this looks.
        gene (str) : the name of the gene.
        min_mean_fraction_spliced_in (float between 0 and 1) : the minimum fraction spliced in (matched/aligned reads) to include. This is useful for excluding
            introns, which makes future FDR correction better, since introns probably aren't what you're hypothesis testing and this analysis pipeline can be
            rather liberal with the splitting of regions. Setting to 0 will disable this filter.
        min_length (int) : minimum length of region to include.
    
    Returns
        region_dict (dict) : dictionary detailing the gene, indices, and statistical values for the region:
            pvalue is calculated based on the log-likelihood functions of the negative binomial fits for both samples and a null hypothesis where
                all pileups are treated as one sample.
            fraction_change is the fraction spliced in of sample_2 - sample_1, so this can range from -1 to 1.
    """
    region_dict = {}
    length = region[0][1] - region[0][0]
    if length >= min_length :
        sample_1_spliced_ins = np.array([ np.array(subregion[1][0]) * subregion[0] for subregion in region[1] ]).sum(axis=0) / length
        sample_1_totals = np.array([ np.array(subregion[1][1]) * subregion[0] for subregion in region[1] ]).sum(axis=0) / length
        sample_2_spliced_ins = np.array([ np.array(subregion[1][2]) * subregion[0] for subregion in region[1] ]).sum(axis=0) / length
        sample_2_totals = np.array([ np.array(subregion[1][3]) * subregion[0] for subregion in region[1] ]).sum(axis=0) / length
        mean_spliced_in = np.mean(np.concatenate([ sample_1_spliced_ins, sample_2_spliced_ins ]))
        mean_fraction_spliced_in = np.mean(np.divide( np.concatenate([ sample_1_spliced_ins, sample_2_spliced_ins ]), np.concatenate([ sample_1_totals, sample_2_totals ]) ))
        if mean_fraction_spliced_in >= min_mean_fraction_spliced_in :
            sample_1_spliced_ins = np.maximum( np.round(sample_1_spliced_ins), 1)
            sample_1_totals = np.round(sample_1_totals)
            sample_1_spliced_outs = np.maximum( np.subtract( sample_1_totals, sample_1_spliced_ins ), 1)
            sample_2_spliced_ins = np.maximum( np.round(sample_2_spliced_ins), 1)
            sample_2_totals = np.round(sample_2_totals)
            sample_2_spliced_outs = np.maximum( np.subtract( sample_2_totals, sample_2_spliced_ins ), 1)
            combined_spliced_ins = np.concatenate([ sample_1_spliced_ins, sample_2_spliced_ins ])
            combined_totals = np.concatenate([sample_1_totals, sample_2_totals])
            combined_spliced_outs = np.concatenate([sample_1_spliced_outs, sample_2_spliced_outs])
            def sample_1_nllf(params) :
                return -np.sum([ stats.nbinom.logpmf(k, n, params[0]) for k, n in zip( sample_1_spliced_outs, sample_1_spliced_ins ) ])
            def sample_2_nllf(params) :
                return -np.sum([ stats.nbinom.logpmf(k, n, params[0]) for k, n in zip( sample_2_spliced_outs, sample_2_spliced_ins ) ])
            def combined_nllf(params) :
                return -np.sum([ stats.nbinom.logpmf(k, n, params[0]) for k, n in zip( combined_spliced_outs, combined_spliced_ins ) ])
            sample_1_minimized = optimize.differential_evolution(sample_1_nllf, bounds=[(0,1)], strategy='best1bin', x0=[min(1, np.mean(np.divide( sample_1_spliced_ins, sample_1_totals )))], maxiter = 100, tol=0.01, updating='deferred')
            sample_2_minimized = optimize.differential_evolution(sample_2_nllf, bounds=[(0,1)], strategy='best1bin', x0=[min(1, np.mean(np.divide( sample_2_spliced_ins, sample_2_totals )))], maxiter = 100, tol=0.01, updating='deferred')
            combined_minimized = optimize.differential_evolution(combined_nllf, bounds=[(0,1)], strategy='best1bin', x0=[min(1, np.mean(np.divide( combined_spliced_ins, combined_totals )))], maxiter = 100, tol=0.01, updating='deferred')
            if sample_1_minimized.success and sample_2_minimized.success and combined_minimized.success :
                p_value = np.exp(-combined_nllf([combined_minimized.x[0]])) / np.exp(-np.sum([ sample_1_nllf([sample_1_minimized.x[0]]), sample_2_nllf([sample_2_minimized.x[0]]) ]))
                region_dict['gene'] = gene
                region_dict['indices'] = region[0]
                region_dict['pvalue'] = p_value
                region_dict['fraction_change'] = sample_2_minimized.x[0] - sample_1_minimized.x[0]
                region_dict['avg_sample_1_counts'] = np.mean(sample_1_totals)
                region_dict['avg_sample_2_counts'] = np.mean(sample_2_totals)
    return region_dict

def pileups_to_regions(sample_1_dict, sample_2_dict, min_average_read_depth = 20, min_fraction_change = 0.1, min_read_depth = 5, min_mean_fraction_spliced_in = 0.05, min_length = 15) :
    """
    This is a utility function for compare_pileups that analyzes pileups per gene to find interesting regions. See get_regions() and splicing_region_neg_binom()
        for more details on how this works.
    Args :
        sample_1_dict (dict) : gene:[pileup] pairings for sample 1. Should contain replicates, i.e. multiple pileups per gene.
        sample_2_dict (dict) : gene:[pileup] pairings for sample 2. Should contain replicates, i.e. multiple pileups per gene.
        min_average_read_depth (int) : the minimum number of average reads aligned or covering each region. Regions that don't meet this will be ignored.
            Averaged across all pileups. Recommend setting to at LEAST 1.
        min_fraction_change (float between 0 and 1) : The minimum change in the fraction of matched:aligned reads (fraction spliced in) to establish a new region.
        min_read_depth (int) : the minimum number of reads per pileup aligned or covering each region. Regions that don't meet this will be ignored.
            Not averaged across pileups, so all pileups must meet this minimum. Recommend setting to at LEAST 1.
        min_mean_fraction_spliced_in (float between 0 and 1) : the minimum fraction spliced in (matched/aligned reads) to include. This is useful for excluding
            introns, which makes future FDR correction better, since introns probably aren't what you're hypothesis testing and this analysis pipeline can be
            rather liberal with the splitting of regions. Setting to 0 will disable this filter.
        min_length (int) : minimum length of region to include.
        
    """
    processed_regions_dict = {'gene' : [], 'indices' : [], 'pvalue' : [], 'fraction_change' : [], 'avg_sample_1_counts' : [], 'avg_sample_2_counts' : []}
    for gene in list(sample_1_dict.keys()) :
        if (gene not in [None, 'none', 'None']) and (gene in sample_2_dict) :
            if None not in sample_1_dict[gene] + sample_2_dict[gene] :
                regions = get_regions( sample_1_dict[gene], sample_2_dict[gene], min_average_read_depth=min_average_read_depth,
                                      min_fraction_change=min_fraction_change, min_read_depth=min_read_depth )
                for region in regions :
                    region_dict = splicing_region_neg_binom(region, gene, min_mean_fraction_spliced_in=min_mean_fraction_spliced_in, min_length=min_length)
                    for key in region_dict :
                        processed_regions_dict[key].append(region_dict[key])
    return processed_regions_dict

def compare_pileups(samples, path_pileup_table, workers = 4, print_results = True, min_average_read_depth = 20, min_fraction_change = 0.1, min_read_depth = 5, min_mean_fraction_spliced_in = 0.05, min_length = 15) :
    """
    Looks in the pileup table generated by create_pileups and compares the pileups between the indicated samples, finding any regions in the pileups that
        differ between them. This is done by aligning the pileups together, then going through each region and comparing between them by looking at coverage
        to determine difference between samples and overall count number to normalize and also get significance values. Refer to pileups_to_regions(),
        get_regions(), and splicing_region_neg_binom() for details.
    
    Args:
        samples (list of str) : the samples to compare.
        path_pileup_table (str) : path to the pileup table output by create_pileups.
        workers (int) : number of parallel processes to run. This is fairly memory efficient and computationally complex, so use as many processes as possible.
        print_results (bool) : whether to print the resulting table of regions.
        min_average_read_depth (int) : the minimum number of average reads aligned or covering each region. Regions that don't meet this will be ignored.
            Averaged across all pileups. Recommend setting to at LEAST 1.
        min_fraction_change (float between 0 and 1) : The minimum change in the fraction of matched:aligned reads (fraction spliced in) to establish a new region.
        min_read_depth (int) : the minimum number of reads per pileup aligned or covering each region. Regions that don't meet this will be ignored.
            Not averaged across pileups, so all pileups must meet this minimum. Recommend setting to at LEAST 1.
        min_mean_fraction_spliced_in (float between 0 and 1) : the minimum fraction spliced in (matched/aligned reads) to include. This is useful for excluding
            introns, which makes future FDR correction better, since introns probably aren't what you're hypothesis testing and this analysis pipeline can be
            rather liberal with the splitting of regions. Setting to 0 will disable this filter.
        min_length (int) : minimum length of region to include.
    
    Returns :
        regions_table (dict) : table of the selected regions with columns for gene, indices, pvalue, padj, fraction_change, and average counts for the samples.
            each row is a region, and the FDR correction for padj is Benjamini-Hochberg.
    """
    pileup_ds = dataset.dataset(path_pileup_table)
    sample_1_table = pileup_ds.filter( pc.field('barcode_sample_ID') == samples[0] ).to_table().drop_columns(['barcode_sample_ID', 'barcode_ID'])
    sample_2_table = pileup_ds.filter( pc.field('barcode_sample_ID') == samples[1] ).to_table().drop_columns(['barcode_sample_ID', 'barcode_ID'])
    print("finished prep")
    with concurrent.futures.ProcessPoolExecutor( max_workers = workers ) as executor :
        futures = []
        column_chunks = [ [ int(y) for y in x ] for x in np.array_split( np.linspace(0, sample_1_table.num_columns - 1, num = sample_1_table.num_columns, dtype=int), workers ) ]
        for column_chunk in column_chunks :
            futures.append(executor.submit( pileups_to_regions,
                                           sample_1_table.select(column_chunk).to_pydict(), sample_2_table.select(column_chunk).to_pydict(),
                                           min_average_read_depth=min_average_read_depth, min_fraction_change=min_fraction_change, min_read_depth=min_read_depth,
                                           min_mean_fraction_spliced_in=min_mean_fraction_spliced_in, min_length=min_length
                                          ))
        print("finished sending functions")
        concurrent.futures.wait( futures )
    regions_dict = {'gene' : [], 'indices' : [], 'pvalue' : [], 'fraction_change' : [], 'avg_sample_1_counts' : [], 'avg_sample_2_counts' : []}
    print("finished futures")
    for future in futures :
        for key in future.result() :
            regions_dict[key] = regions_dict[key] + future.result()[key]
    regions_table = pa.table(regions_dict).sort_by('pvalue')
    padj = []
    N = regions_table.num_rows
    for rank, pvalue in zip(range(1,N+1), regions_table.column('pvalue').to_pylist()) :
        if pvalue == None :
            padj.append(1)
        else :
            padj.append(min(1, pvalue * N / rank))
    regions_table = regions_table.append_column('padj', [padj])
    if print_results :
        display_table = regions_table.sort_by('pvalue').slice(0,100)
        print(regions_table.num_rows)
        print( "{:>20}{:<25}{:15}{:15}{:15}{:20}{:20}".format('gene|', 'indices', 'pvalue', 'padj', 'fraction_change', 'avg_sample_1_counts', 'avg_sample_2_counts') )
        for row in display_table.to_pylist() :
            print( "{:>20}{:<25}{:<15.5}{:<15.5}{:<15.5}{:<20.0f}{:<20.0f}".format( row['gene'], str(row['indices']), row['pvalue'], row['padj'], row['fraction_change'], row['avg_sample_1_counts'], row['avg_sample_2_counts'] ) )
    return regions_table

def show_pileup_for_gene(path_out_pileup, gene, samples = None, colors_dict = None, target_indices = None, average_replicates = True, show_x_labels = True, scale_bar = False, save_fig = None, index_start = 0, reverse = False ) :
    """
    Displays the pileups for a given gene within a pileup file from create_pileups() as a coverage plot.

    Args :
        path_out_pileup (str) : the path for the pileup file from create_pileups().
        gene (str) : the gene to pull data for.
        samples (list of sample names (str)) : the names of samples to show. If None, this will show ALL samples.
        colors_dict (dict) : pairings of sample:color to use for the plot. The color can be given in any way that matplotlib accepts. If left as None, this
            will use the default TABLEAU_COLORS of matplotlib. Specifying your colors is the best way to know what line on the plot is what sample.
        target_indices (list of [start, end]) : the 0-indexed, right-exclusive indices to display. If None, this will show the full region spanned by the
            pileups of the gene.
        average_replicates (bool) : whether to average the pileups for each sample into one line.
        show_x_labels (bool) : whether to show the x_label indices.
        scale_bar (False or int) : the size of scale bar to show, in number of bases. Leave False to disable.
        save_fig (str) : path to save the figure to. Leave False to disable.
        index_start (int) : the position along the reference to re-index the data to. This can be nice to show the x axis as starting from the end of a gene.
        reverse (bool) : whether to display the plot as reversed.
    Returns : None
    """
    if samples != None :
        gene_table = pq.read_table( path_out_pileup, columns=['barcode_sample_ID', gene], filters = pc.field('barcode_sample_ID').isin(samples) ).drop_null()
    else :
        gene_table = pq.read_table( path_out_pileup, columns=['barcode_sample_ID', gene] ).drop_null()
        samples = gene_table.column('barcode_sample_ID').unique().to_pylist()
    if colors_dict == None :
        colors_dict = {}
        default_color_list = list(mpl.colors.TABLEAU_COLORS.values())
        for i, sample in enumerate(samples) :
            colors_dict[sample] = default_color_list[i]
    if target_indices != None and type(target_indices[0]) != list :
        target_indices = [target_indices]
    gene_table = gene_table.set_column(1, gene, [pad_compact_seqs(gene_table.column(gene).to_pylist(), [0, 0])])
    first_covered_index = get_total_cigar_length(gene_table.column(gene).to_pylist()[0])
    for compact_seq in gene_table.column(gene).to_pylist() :
        seq_index = 0
        for x in compact_seq :
            if seq_index < first_covered_index :
                if x[1] > 0 :
                    first_covered_index = seq_index
                else :
                    seq_index += x[0]
            else :
                break
    fig, ax = plt.subplots(1, 1)
    for sample in samples :
        sample_table = gene_table.filter( pc.field('barcode_sample_ID') == sample )
        pileups = sample_table.column(gene).to_pylist()
        pileup_length = get_total_cigar_length(pileups[0])
        if target_indices != None :
            clipped_pileups = []
            for indices in target_indices :
                right_clip = pileup_length - indices[1]
                clipped_pileups.append([ clip_cigar_ends(pileup, indices[0] - 1, right_clip) for pileup in copy.deepcopy(pileups) ])
            num_targets = len(target_indices)
            for i in range(len(clipped_pileups[0])) :
                concat_pileup = []
                for j in range(num_targets) :
                    concat_pileup += clipped_pileups[j][i]
                pileups[i] = concat_pileup
        else :
            pileups = [ clip_cigar_ends(pileup, first_covered_index, 0) for pileup in pileups ]
        pileups_aligned = align_compact_seqs(pileups)
        if average_replicates :
            summed_pileups = np.sum(pileups_aligned[:,:, 1:], axis=0)
            pileups_aligned = np.expand_dims(np.insert(summed_pileups, 0, pileups_aligned[0,:,0], axis=1), 0)
        for pileup in pileups_aligned :
            fill_vals = np.divide( pileup[:,2], pileup[:,1], out = np.full(len(pileup), None), where = pileup[:,1] > 0 )
            coverages = np.repeat( fill_vals, pileup[:,0] )
            percents = np.multiply( coverages, 100, out = np.full(len(coverages), None), where = np.not_equal(coverages, np.full(len(coverages), None)) )
            if not reverse :
                ax.plot( range(len(percents)), percents, color = colors_dict[sample] )
            else :
                ax.plot( range(len(percents), 0, -1), percents, color = colors_dict[sample] )
    if target_indices != None :
        display_indices = []
        axis_break_indices = [0]
        for indices in target_indices :
            index_list = list(range(indices[0], indices[1]+1))
            display_indices += index_list
            axis_break_indices.append(axis_break_indices[-1] + len(index_list))
        display_length = len(display_indices)
        labels = np.absolute(np.take(display_indices, np.linspace(0, display_length - 1, num=10, dtype=np.int64())) - index_start)
        if len(axis_break_indices) > 2 :
            kwargs = dict(marker=[(-0.5, -1), (0.5, 1)], markersize=6, linestyle="none", color='k', mec='k', mew=1, clip_on=False)
            for break_index in axis_break_indices[1:-1] :
                ax.plot( [break_index / display_length], [0], transform=ax.transAxes, **kwargs )
    else :
        display_indices = [first_covered_index, pileup_length]
        display_length = pileup_length - first_covered_index
        labels = np.absolute(np.linspace(display_indices[0], display_indices[1], num=10, dtype=np.int64()) - index_start)
    ax.set_xlim(0, display_length - 1)
    ticks = np.linspace(0, display_length - 1, num=10, dtype=np.int64())
    if reverse :
        ticks = np.flip(ticks)
    if show_x_labels :
        ax.set_xticks(ticks, labels = labels, rotation=45 )
    else :
        ax.axes.get_xaxis().set_visible(False)
    if scale_bar :
        sb_artist = AnchoredSizeBar( ax.transAxes, scale_bar / display_length, str(scale_bar) + ' bp', loc='lower left', bbox_to_anchor=(1, 0), bbox_transform=ax.transAxes , pad=0.1, borderpad=0.2, sep=5, frameon=False )
        ax.add_artist(sb_artist)
    ax.grid(axis='y', zorder=0)
    ax.set_axisbelow(True)
    ax.set_ylim(-1,101)
    if save_fig != None :
        fig.savefig(save_fig)
    print('done!')
    return

def align_inserts(compact_seqs, return_num_inserts = False, min_insert_len = 3) :
    """
    Aligns the insert regions of CIGAR sequences and similarly formatted data. This is meant to be used with show_alignments(). This finds the inserts
        within the seqs and adds gap regions into all other seqs, so that all seqs are aligned and the same length.

    Args :
        compact_seqs (list of [length, values...]) : the compact seqs (CIGAR-style packaging of sequence/alignment data). These need to be equal lengths,
            ignoring inserts, so you should process CIGARs with process_cigar_strings().
        return_num_inserts (bool) : whether to return the total number of inserts along with the compact_seqs.
        min_insert_len (int) : minimum length of insert to include. Shorter inserts are ignored, which can be nice to remove noise artifacts.
    Returns :
        compact_seqs (list of [length, values...]) : the compact_seqs aligned.
    """
    inserts = []
    for compact_seq in compact_seqs :
        non_i_index = 0
        i = 0
        while i < len(compact_seq) :
            if compact_seq[i][1] == 'I' :
                if compact_seq[i][0] >= min_insert_len :
                    inserts.append([non_i_index, compact_seq[i][0]])
                    i += 1
                else :
                    del compact_seq[i]
            else :
                non_i_index += compact_seq[i][0]
                i += 1
    inserts.sort(key = lambda x : x[1], reverse=True)
    inserts.sort(key = lambda x : x[0])
    inserts_non_redundant = []
    non_i_index = 0
    max_index = get_total_cigar_length(compact_seqs[0], chars_to_ignore=['I'])
    for insert in inserts :
        if insert[0] > non_i_index :
            inserts_non_redundant.append(insert)
            non_i_index = insert[0]
    total_inserts = sum([ x[1] for x in inserts_non_redundant ])
    if max_index < non_i_index :
        insert_seq.append([max_index - non_i_index, ' '])
    for compact_seq in compact_seqs :
        non_i_index = 0
        i = 0
        for insert in inserts_non_redundant :
            while non_i_index < insert[0] :
                if compact_seq[i][0] <= insert[0] - non_i_index :
                    non_i_index += compact_seq[i][0]
                else :
                    compact_seq.insert(i+1, [ non_i_index + compact_seq[i][0] - insert[0], compact_seq[i][1] ])
                    compact_seq[i][0] = insert[0] - non_i_index
                    non_i_index = insert[0]
                i += 1
            if compact_seq[i][1] == 'I' :
                if compact_seq[i][0] < insert[1] :
                    compact_seq.insert(i+1, [ insert[1] - compact_seq[i][0], ' ' ])
                    i += 1
            else :
                compact_seq.insert(i, [ insert[1], ' ' ])
            i += 1
    if return_num_inserts :
        return compact_seqs, total_inserts
    else :
        return compact_seqs

def load_alignment_table(path_dataset, genes = None, ref_indices = None, ref_ctg = None, samples_to_compare = None) :
    """
    Convenience function to load the relevant data from a pyarrow dataset in order to look at the alignments.

    Args :
        path_dataset (str or list of str) : the path(s) to the dataset(s).
        genes (str or list of str) : the gene(s) to include.
        ref_indices (list of int [start, end]) : the 0-indexed indices along the reference to filter the data by.
            Only selects reads that overlap with this region. MUST be used with ref_ctg, otherwise the data will be garbage.
        ref_ctg (str) : the contig to filter by.
        samples_to_compare (list of str) : the samples to select.
    
    Returns :
        gene_table (pyarrow Table) : the filtered data.
    """
    print("loading data")
    if type(path_dataset) == list :
        ds = dataset.dataset([ dataset.dataset(x) for x in path_dataset ])
    else :
        ds = dataset.dataset(path_dataset)
    if samples_to_compare != None :
        ds = ds.filter(pc.field('barcode_sample_ID').isin(samples_to_compare))
    if genes != None :
        if type(genes) != list :
            genes = [genes]
        ds = ds.filter(pc.field('minimap2_gene_ID').isin(genes))
    if ref_indices != None :
        ds = ds.filter((pc.field('minimap2_ctg') == ref_ctg) & (pc.field('minimap2_r_st') <= ref_indices[1]) & (pc.field('minimap2_r_en') >= ref_indices[0]))
    gene_table = ds.to_table( columns = [ 'barcode_ID', 'barcode_sample_ID', 'minimap2_cigar', 'minimap2_r_st', 'minimap2_r_en', 'minimap2_gene_ID', 'minimap2_ctg', 'minimap2_mlen', 'minimap2_strand', 'barcode_direction' ])
    print("done loading")
    return gene_table
    
def align_cig_and_seq(cig, seq, del_inserts = False, add_gaps = True) :
    """
    Aligns a sequence with a corresponding CIGAR alignment representation. This can remove inserts (parts of the sequence that don't fit into the alignment 
        reference) and add gaps into the sequence where there are deletions (parts of the alignment reference that the sequence doesn't have).
        
    Args :
        cig (list of [length (int), value (char)]) : the parsed CIGAR. Must be parsed, cannot be the raw string.
        seq (str) : the sequence that the CIGAR represents.
        del_inserts (bool) : whether or not to delete inserts from the sequence (indicated by 'I' in the CIGAR).
        add_gaps (bool) : whether or not to add spaces into the sequence to represent deletions in the CIGAR (indicated by 'D', 'N', or ' ').
        
    Returns :
        algined_seq (str) : the edited sequence.
    """
    if type(cig) == str :
        raise TypeError("Error: this cig doesn't appear to be parsed. Must be a parsed cigar in the form list of [length (int), value (char)].")
    i = 0
    aligned_seq = ''
    for item in cig :
        if item[1] == 'M' :
            aligned_seq += seq[ i : i + item[0] ]
            i += item[0]
        elif item[1] in ['D', 'N', ' '] :
            if add_gaps :
                aligned_seq += ' ' * item[0]
        elif item[1] == 'I' :
            if not del_inserts :
                aligned_seq += seq[ i : i + item[0] ]
            i += item[0]
    return aligned_seq

def process_table_cigs(table, ignore_I = False, index_0 = False, minimum_region_length = 0, ignore_ends = False, max_match = 1, remove_outliers = False, left_clip = 0, right_clip = 0) :
    """
    Takes a table of mapped reads and processes the CIGARs together to give processed cigs and the table of data for reads with valid cigs.

    Args :
        table (pyarrow.Table) : table with the reads. All included reads will be processed together, so this should only include desired reads.
        ignore_I (bool) : whether or not to remove inserts in the cigs.
        index_0 (bool) : whether or not to include a 'gap' at the start of cigs such that all cigs start at index 0 on the alignment reference sequence.
        minimum_region_length (int) : the minimum length of items to return. If greater than 0, smaller items will be added to the preceeding item.
            This can be useful if only looking for larger features, as reducing the total number of items significantly speeds up future processing.
        ignore_ends (bool) : whether or not to ignore the first and last region of alignments. This can just simplify the visual.
        max_match (float between 0 and 1) : the maximum allowed fraction of aligned reference sequence that is matched by each read. This is useful to exclude reads that
            are completely or nearly matched to the reference, which can be bad alignments or weird transcripts, if you're working with mRNA that is supposed
            to have introns. Set to 1 to not filter anything out.
        remove_outliers (bool) : whether or not to remove outlier alignments with remove_outlier_alignments().
        left_clip (int) : length of cig to remove from the left end.
        right_clip (int) : length of cig to remove from the right end.
        
    Returns :
        parsed_cigs ( list of [list of [length (int), value (char)]] ) : the processed CIGAR sequences clipped and padded as specified.
        table (pyarrow.Table) : the same data as the input table, filtered to only include the reads that had valid CIGAR sequences after processing.
    """
    ctgs = table.column('minimap2_ctg')
    if len(ctgs.unique().to_pylist()) > 1 :
        print("Warning: multiple contig IDs found for this gene. Check alignments to verify that they are all correctly matched. Moving forward with the most abundant contig ID. Contigs: ", ctgs.value_counts().to_pylist())
        ctg = pa.table( table.column('minimap2_ctg').value_counts() ).sort_by( [('counts', 'descending')] ).column('values').to_pylist()[0]
        table = table.filter( pc.field('minimap2_ctg') == ctg )
    else :
        ctg = ctgs.unique().to_pylist()[0]
    if max_match != 1 :
        table = table.filter(pc.less_equal( pc.divide( pc.field('minimap2_mlen').cast(pa.float32()), pc.abs( pc.subtract( pc.field('minimap2_r_en'), pc.field('minimap2_r_st') ) ) ), pc.scalar(max_match) ))
    if remove_outliers :
        table = remove_outlier_alignments(table)
    table_dict = table.to_pydict()
    parsed_cigs, passed_indices = process_cigar_strings(
        table_dict['minimap2_cigar'], table_dict['minimap2_r_st'], table_dict['minimap2_r_en'], index_0 = index_0, ignore_I = ignore_I,
        minimum_region_length = minimum_region_length, return_null_cigs = False, ignore_ends = ignore_ends, return_passed_cigs_indices = True,
        left_clip = left_clip, right_clip = right_clip
    )
    return parsed_cigs, table.take(passed_indices)

def show_alignments_visual(table, target_indices = None, remove_outliers = False, minimum_region_length = 0, max_match = 1, ignore_ends = False, ignore_I = True, reduce_introns = False, save_fig = None) :
    """
    Shows the alignments for given data, aligned to each other.
    
    Args :
        table (pyarrow Table) : the data including alignments that you'd like to show. Can be from load_alignment_table, or from custom filtering of a
            pyarrow dataset that has been processed through the map_count functions.
        target_indices (list of int [start, end]) : the 0-indexed region to limit the display to. If None, this will show the full width of all alignments.
        remove_outliers (bool) : whether or not to remove outlier alignments with remove_outlier_alignments().
        max_match (float between 0 and 1) : the maximum allowed fraction of aligned reference sequence that is matched by each read. This is useful to exclude reads that
            are completely or nearly matched to the reference, which can be bad alignments or weird transcripts, if you're working with mRNA that is supposed
            to have introns. Set to 1 to not filter anything out.
        ignore_ends (bool) : whether or not to ignore the first and last region of alignments. This can just simplify the visual.
        ignore_I (bool) : whether or not to ignore inserts in the alignments.
        reduce_introns (int or False) : the length to shorten long introns to. Introns are considered regions that have a global average fraction spliced-in
            of less than 0.05. Set to False to disable.
        save_fig (str) : path to save the figure to.
        
    Returns :
        None
    """
    print("starting gap and parse")
    r_sts = table.column('minimap2_r_st').to_pylist()
    r_ens = table.column('minimap2_r_en').to_pylist()
    if target_indices != None :
        cigs_parsed, table_passed = process_table_cigs(
            table, ignore_I = ignore_I, minimum_region_length = minimum_region_length, ignore_ends = ignore_ends, max_match = max_match, index_0 = True,
            remove_outliers = remove_outliers, left_clip = target_indices[0], right_clip = max(r_ens) - target_indices[1]
        )
    else :
        cigs_parsed, table_passed = process_table_cigs(
            table, ignore_I = ignore_I, minimum_region_length = minimum_region_length, ignore_ends = ignore_ends, max_match = max_match, index_0 = False,
            remove_outliers = remove_outliers
        )
    sample_IDs = table_passed.column('barcode_sample_ID').to_pylist()
    barcode_IDs = table_passed.column('barcode_ID').to_pylist()
    r_sts = table_passed.column('minimap2_r_st').to_pylist()
    r_ens = table_passed.column('minimap2_r_en').to_pylist()
    if not ignore_I :
        cigs_parsed, total_inserts = align_inserts(cigs_parsed, return_num_inserts = True)
    time_parse = time.time()
    print("finished gap and parse")
    cig_value_dict = {'M':0, 'D':1, 'N':1, ' ':2, 'I':2}
    identified_cigs = list(zip( cigs_parsed, sample_IDs, r_sts, r_ens ))
    identified_cigs.sort( key = lambda x : sum([ cig_value_dict[y[1]] * y[0] for y in x[0] ]) )
    r_sts_passed = []
    r_ens_passed = []
    rows_by_sample = {}
    for cig, sample_ID, r_st, r_en in identified_cigs :
        if ''.join([ x[1] for x in cig ]).count(' ') < len(cig) :
            r_sts_passed.append(r_st)
            r_ens_passed.append(r_en)
            if sample_ID not in rows_by_sample :
                rows_by_sample[sample_ID] = []
            color_row = np.array([])
            for item in cig :
                if item[1] == 'M' :
                    color_row = np.concatenate((color_row, np.full(item[0], np.int8(1))))
                elif item[1] == ' ' :
                    color_row = np.concatenate((color_row, np.full(item[0], np.int8(0))))
                elif item[1] == 'I' :
                    color_row = np.concatenate((color_row, np.full(item[0], np.int8(2))))
                else :
                    color_row = np.concatenate((color_row, np.full(item[0], np.int8(-1))))
            rows_by_sample[sample_ID].append(color_row)
    print("finished making colors")
    if target_indices != None :
        indices = target_indices
    else :
        indices = [min(r_sts_passed), max(r_ens_passed)]
    if not ignore_I :
        indices[1] += total_inserts
    print(min(r_sts_passed), max(r_ens_passed))
    fig, axes = plt.subplots(len(rows_by_sample), sharex=True)
    cmap = mpl.colors.ListedColormap(['lightsalmon','white','darkblue','green'])
    if reduce_introns :
        colors_array = np.concatenate([ np.array(rows) for rows in rows_by_sample.values() ])
        valid_indices = []
        num_conceq_intron_cols = 0
        for i, col in enumerate(colors_array.mean( axis=0, where= colors_array != 0 )) :
            if col < -0.95 :
                num_conceq_intron_cols += 1
                if num_conceq_intron_cols < reduce_introns :
                    valid_indices.append(i)
            else :
                valid_indices.append(i)
                num_conceq_intron_cols = 0
        for sample in rows_by_sample :
            rows_by_sample[sample] = np.take(rows_by_sample[sample], valid_indices, axis=1)
        display_indices = np.take(range(indices[0], indices[1]+1), valid_indices)
        labels = np.take(display_indices, np.linspace(0, len(display_indices) - 1, num=50, dtype=np.int64()))
    else :
        labels = np.linspace(indices[0], indices[1], num=50, dtype=np.int64())
    ticks = np.linspace(0, len(list(rows_by_sample.values())[0][0]), num=50, dtype=np.int64())
    if len(rows_by_sample) == 1 :
        for sample in rows_by_sample :
            axes.imshow(rows_by_sample[sample], interpolation = 'none', aspect = 'auto', cmap = cmap, vmin=-1, vmax=2)
            axes.set_ylabel(sample)
        axes.set_xticks(ticks, labels, rotation='vertical' )
    else :
        for sample, ax in zip(rows_by_sample, axes) :
            ax.imshow(rows_by_sample[sample], interpolation = 'none', aspect = 'auto', cmap = cmap, vmin=-1, vmax=2)
            ax.set_ylabel(sample)
        ax.set_xticks(ticks, labels, rotation='vertical' )
    if save_fig != None :
        fig.savefig(save_fig)
    return

def show_read_alignments(table, reference_seq = None, target_indices = None, ignore_I = False, minimum_region_length = 0, ignore_ends = False, max_match = 1, remove_outliers = False, print_match_indicators = True, print_cigars = False) :
    """
    Prints the given reads aligned together, including an optional reference sequence.
    
    Args :
        table (pyarrow Table) : the data including alignments that you'd like to show. Can be from load_alignment_table, or from custom filtering of a
            pyarrow dataset that has been processed through the map_count functions.
        reference_seq (str) : the sequence that the reads were aligned to.
        target_indices (list of [start (int), end (int)]) : the 0-indexed region to limit the display to. If None, this will show the full width of all alignments.
        ignore_I (bool) : whether or not to ignore inserts in the alignments.
        minimum_region_length (int) : the minimum length of items to return. If greater than 0, smaller items will be added to the preceeding item.
            This can be useful if only looking for larger features, as reducing the total number of items significantly speeds up future processing.
        ignore_ends (bool) : whether or not to ignore the first and last region of alignments. This can just simplify the visual.
        max_match (float between 0 and 1) : the maximum allowed fraction of aligned reference sequence that is matched by each read. This is useful to exclude reads that
            are completely or nearly matched to the reference, which can be bad alignments or weird transcripts, if you're working with mRNA that is supposed
            to have introns. Set to 1 to not filter anything out.
        remove_outliers (bool) : whether or not to remove outlier alignments with remove_outlier_alignments().
        print_match_indicators (bool) : whether or not to print the match indicators where (relative to the reference seq), '|' means a match, 'X' means a mismatch,
            '^' means an insert, '-' means a deletion, and ' ' means a gap where the read has no data (for inserts from other reads or before/after the read).
            Only prints if reference_seq is provided to compare against. Note that the 'M' in the cigar sequences can also indicate a small mismatch due to a
            quirk with minimap2 alignment.
        print_cigars (bool) : whether or not to print the cigar sequences.
        
    Returns :
        None
    """
    cigs_passed, table_passed = process_table_cigs(
        table, ignore_I = ignore_I, minimum_region_length = minimum_region_length, ignore_ends = ignore_ends, max_match = 1,
        remove_outliers = False, index_0 = True
    )
    if len(cigs_passed) != 0 :
        seqs = table_passed.column('seq').to_pylist()
        q_sts = table_passed.column('minimap2_q_st').to_pylist()
        q_ens = table_passed.column('minimap2_q_en').to_pylist()
        bio_seq_indices = table_passed.column('barcode_biological_seq_indices').to_pylist()
        directions = table_passed.column('barcode_direction').to_pylist()
        strands = table_passed.column('minimap2_strand').to_pylist()
        r_sts = table_passed.column('minimap2_r_st').to_pylist()
        r_ens = table_passed.column('minimap2_r_en').to_pylist()
        seqs_passed = []
        for seq, q_st, q_en, bio_seq_index, direction, strand in zip(seqs, q_sts, q_ens, bio_seq_indices, directions, strands) :
            if direction == 'reverse' :
                seq_selected = utils.reverse_complement(seq[bio_seq_index[0] : bio_seq_index[1]])[q_st : q_en]
            else :
                seq_selected = seq[bio_seq_index[0] : bio_seq_index[1]][q_st : q_en]
            if strand == -1 :
                seq_selected = utils.reverse_complement(seq_selected)
            seqs_passed.append(seq_selected)
        if target_indices == None :
            target_indices = [min(r_sts), max(r_ens)]
        cig_clip_results = [ clip_cigar_ends( cig, target_indices[0], max( 0, max(r_ens) - target_indices[1] ), return_clipped_seq_lengths= True ) for cig in cigs_passed ]
        cigs_clipped = [ x[0] for x in cig_clip_results ]
        clips = [ (x[1], x[2]) for x in cig_clip_results ]
        seqs_passed = [ x[clip[0]:] for x, r_st, clip in zip(seqs_passed, r_sts, clips) ]
        if reference_seq != None :
            cigs_clipped = [[[target_indices[1] - target_indices[0], 'M']]] + cigs_clipped
            seqs_passed = [ reference_seq[target_indices[0] : target_indices[1]+1].upper() ] + seqs_passed
        if ignore_I :
            cigs_aligned = align_compact_seqs(cigs_clipped, return_numpy = False)
            seqs_passed = [ align_cig_and_seq(cig, seq, del_inserts = True, add_gaps = True) for cig, seq in zip(cigs_aligned, seqs_passed) ]
        else :
            cigs_aligned = align_inserts(cigs_clipped, min_insert_len = 1)
            seqs_passed = [ align_cig_and_seq(cig, seq, del_inserts = False, add_gaps = True) for cig, seq in zip(cigs_aligned, seqs_passed) ]
        cigs_expanded = [ cigar_to_str(cig) for cig in cigs_aligned ]
        step_size = 250
        num_inserts = 0
        for i in range(0, len(seqs_passed[0]), step_size) :
            left_index_statement = 'V ' + str(target_indices[0] + i - num_inserts)
            line_length = len(cigs_expanded[0][i:i+step_size])
            num_inserts_in_line = len(re.findall("[I ]", cigs_expanded[0][i:i+step_size]))
            num_inserts += num_inserts_in_line
            right_index_statement = str(target_indices[0] + i + line_length - num_inserts) + ' V'
            index_statement_gap_len = line_length - len(left_index_statement) - len(right_index_statement)
            if index_statement_gap_len >= 5 :
                index_statement_gap = index_statement_gap_len * ' '
                print(left_index_statement + index_statement_gap + right_index_statement)
            elif index_statement_gap_len + len(right_index_statement) >= 0 :
                index_statement_gap = (line_length - len(right_index_statement)) * ' '
                print(index_statement_gap + right_index_statement)
            else :
                index_statement_gap = (line_length - 1) * ' '
                right_index_statement = 'V ' + str(target_indices[0] + i + line_length - num_inserts)
                print(index_statement_gap + right_index_statement)
                
            for j, cig, seq in zip(range(len(cigs_expanded)), cigs_expanded, seqs_passed) :
                if reference_seq != None and j == 0 :
                    ref_seq_selected = seq[i:i+step_size]
                    print(ref_seq_selected)
                else :
                    seq_selected = seq[i:i+step_size]
                    cig_selected = cig[i:i+step_size]
                    if print_match_indicators and reference_seq :
                        gaps = np.equal(list(cig_selected), ' ')
                        matches = np.logical_and( np.equal(list(ref_seq_selected), list(seq_selected)), np.logical_not(gaps) )
                        inserts = np.equal(list(cig_selected), 'I')
                        dels = np.logical_or( np.equal(list(cig_selected), 'N'), np.equal(list(cig_selected), 'D') )
                        indicators = np.where( matches, np.full(len(matches), '|'), np.full(len(matches), 'X') )
                        indicators = np.where( gaps, np.full(len(matches), ' '), indicators )
                        indicators = np.where( inserts, np.full(len(matches), '^'), indicators )
                        indicators = np.where( dels, np.full(len(matches), '-'), indicators )
                        print(''.join(indicators))
                    if print_cigars :
                        print(cig_selected)
                    print(seq_selected)
            print()
    return

def get_PSI_in_gene(path_pileup, gene, target_indices, reverse=False, quiet=False) :
    """
    Gets the percent spliced in (PSI) in a gene within given indices, based on the pileups generated with create_pileups(). Returns and prints the data.

    Args :
        path_pileup (str) : path to the pileup file generated with create_pileups().
        gene (str) : the gene to select.
        target_indices ([start (int), end (int)]) : the 0-indexed indices to get the PSI for.
        reverse (bool) : whether to provide values as percent spliced rather than percent spliced in.

    Returns :
        results_table (pyarrow Table) : the resulting PSI values per sample with error values (+,- 2 standard errors) and extra data.
    """
    gene_table = pq.read_table( path_pileup, columns=['barcode_sample_ID', gene] ).drop_null()
    gene_table = gene_table.set_column(1, gene, [pad_compact_seqs(gene_table.column(gene).to_pylist(), [0, 0])])
    pileups = gene_table.column(gene).to_pylist()
    pileup_length = get_total_cigar_length(pileups[0])
    right_clip = pileup_length - target_indices[1]
    gene_table = gene_table.set_column(1, gene, [[ clip_cigar_ends(pileup, target_indices[0], right_clip) for pileup in pileups ]])
    samples = gene_table.column('barcode_sample_ID').unique().to_pylist()
    results_dict = {
        'sample' : [],
        'average_PSI' : [],
        'std_dev_errors' : [],
        'counts' : [],
        'PSI_values' : []
    }
    for sample in samples :
        sample_table = gene_table.filter( pc.field('barcode_sample_ID') == sample )
        pileups = sample_table.column(gene).to_pylist()
        pileups_aligned = align_compact_seqs(pileups)
        PSIs = []
        counts = []
        for pileup in pileups_aligned :
            fill_vals = np.divide( pileup[:,2], pileup[:,1], out = np.full(len(pileup), None), where = pileup[:,1] > 0 )
            coverages = np.repeat( fill_vals, pileup[:,0] )
            if np.all(coverages) :
                PSIs.append(max(0.01, min(0.99, np.average(coverages))))
                counts.append(np.max( pileup[:, 1] ))
        if len(PSIs) > 0 :
            if reverse :
                PSIs = np.subtract(1, PSIs)
            PSIs_logit = np.log(np.divide(PSIs, np.subtract(1, PSIs)))
            PSIs_mean_logit = np.mean(PSIs_logit)
            PSIs_mean = 1 / ( 1 + np.exp(-PSIs_mean_logit) )
            PSIs_stddev_logit = np.abs(np.std(PSIs_logit, ddof=1))
            PSIs_pm_2_stders_logit = [ PSIs_mean_logit + PSIs_stddev_logit, PSIs_mean_logit - PSIs_stddev_logit ]
            PSIs_pm_2_stders = np.divide(1, np.add(1, np.exp(np.subtract(0, PSIs_pm_2_stders_logit))))
            PSIs_pm_2_stders = np.abs(np.subtract(PSIs_mean, PSIs_pm_2_stders))
            results_dict['sample'].append(sample)
            results_dict['average_PSI'].append(100 * PSIs_mean)
            results_dict['std_dev_errors'].append(100*PSIs_pm_2_stders)
            results_dict['PSI_values'].append(np.multiply(PSIs, 100))
            results_dict['counts'].append(counts)
    results_table = pa.table(results_dict).sort_by('sample')
    if not quiet :
        print( "%15s%20s%20s%25s%25s" % ('sample', 'average_PSI', 'std_dev_errors', 'PSI_values', 'counts') )
        for row in results_table.to_pylist() :
            print('%(sample)15s%(average_PSI)15.5f%(std_dev_errors)20s%(PSI_values)25s%(counts)25s' % row)
    return results_table