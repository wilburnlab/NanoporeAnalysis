"""
Tools for performing and manipulating sequence alignments
"""

import edlib
import numpy as np

from NanoporeAnalysis.utils import return_best_orf
from NanoporeAnalysis.utils import reverse_complement as rc


def align(query: str, subject: str, max_score: int, score_only: bool = False):
    """
    Align sequences using edlib and compute the edge distance, i.e. steps to make
    the step sequences match. If score_only, only return the edge distance
    """
    # print( query, subject, max_score, score_only )
    task = "distance" if score_only else "locations"
    alignment = edlib.align(query, subject, mode="HW", task=task, k=max_score)
    return alignment


def merge_overlapped_indices(index_pairs: list, tolerated_mismatches: int) -> list:
    """
    Provided a list of (start,end) index pairs, combine overlapping pairs into
    single pairs that span the full alignment range
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
            reduced_pairs.append({"start": current_start, "end": current_end, "length": align_len})
            current_start = start
            current_end = end
    align_len = current_end - current_start + 1
    reduced_pairs.append({"start": current_start, "end": current_end, "length": align_len})
    return reduced_pairs


def find_polyX(sequence: str, X: str, N: int, tolerated_mismatches: int, min_len: int, max_len) -> list:
    """
    Find runs of single nucleotides (X), using edlib align with search query X*N
    with tolerated_mismatches in the search
    """
    assert len(X) == 1, f"More than one nt provided as quuery ({X}) for find_polyX"
    query = X * N
    alignment = align(query, sequence, max_score=0)
    align_indices = alignment["locations"]
    if len(align_indices) == 0:
        return []
    else:
        merged_indices = merge_overlapped_indices(align_indices, tolerated_mismatches)
        filtered_indices = [i for i in merged_indices if i["length"] >= min_len and i["length"] <= max_len]
        sorted_indices = sorted(filtered_indices, key=lambda x: x["start"], reverse=True)  # Put the last one first
        return sorted_indices


def find_best_polyA(
    sequence: str, search_N: int = 4, allowed_mismatches: int = 1, min_len: int = 15, max_len: int = 40
):  # Patch to min length = 15, max = 40
    """
    Find best polyA tail
    CURRENT VERSION: use longest tail, consider adding some 3' filter in
    """
    polyA_tails = find_polyX(sequence, "A", search_N, allowed_mismatches, min_len, max_len)
    if len(polyA_tails) > 0:
        return polyA_tails[0]
    else:
        return None


def split_sequence_by_polyA(sequence: str):
    """
    Find the best polyA tail, and split sequence accordingly
    Returns the pre-polyA sequence, post-polyA sequence, and the polyA tail itself
    """
    polyA = find_best_polyA(sequence)
    if polyA is not None:
        pre_polyA = sequence[: polyA["start"]]
        post_polyA = sequence[polyA["end"] + 1 :]
        polyA_stats = dict([(f"PolyA {k}", v) for k, v in polyA.items()])
    else:
        pre_polyA = sequence
        post_polyA = None
        polyA_stats = dict()
    return pre_polyA, post_polyA, polyA_stats


def parse_adapter5_seq(sequence: str, ssp: str, max_score: int = 2):
    """
    Find 5' SSP and split the read accordingly. If the SSP is not found, then
    the entire sequence defaults as the transcript sequence.
    Returns the transcript sequence, adapter sequence, and SSP alignment score
    """
    if sequence is None:
        score, trim_idx, adapter_seq = [None] * 3
        transcript_seq = sequence
    else:
        ssp_alignment = align(ssp, sequence, max_score)
        if len(ssp_alignment["locations"]) > 0:
            score = ssp_alignment["editDistance"]  # (len(ssp)-ssp_alignment['editDistance'])/len(ssp)
            ssp_locations = ssp_alignment["locations"]
            max_end = np.max([loc[1] for loc in ssp_locations])
            trim_idx = max_end + 1
            adapter_seq = sequence[:trim_idx]
            transcript_seq = sequence[trim_idx:]
        else:
            score, trim_idx, adapter_seq = [None] * 3
            transcript_seq = sequence
            # score = None #0.5/len(ssp)
            # trim_idx = None
        # adapter_seq = None if trim_idx is None else sequence[:trim_idx]
        # transcript_seq = sequence[trim_idx:]
    return transcript_seq, adapter_seq, score, trim_idx


def score_umi(putative_umi: str, umi_dict: dict, max_score: int = 2) -> dict:
    """
    Provided a dict with {name:umi}, return the scores per umi sorted by rank
    """
    scores = {}
    if putative_umi is None:
        scores = dict([(name, None) for name in umi_dict])
        for k in ["Best UMI", "Best UMI Score", "Delta UMI Score"]:
            scores[k] = None
    else:
        for name, umi in umi_dict.items():
            umi_align = align(rc(umi), putative_umi, max_score=max_score, score_only=False)  # Change max_score to 2
            scores[name] = umi_align["editDistance"]  # umi_score

        select_scores = [i for i in scores.items() if i[1] is not None]
        if len(select_scores) == 0:
            best_umi, best_umi_score, delta_umi_score = [None] * 3
        sorted_scores = sorted(select_scores, key=lambda x: x[1], reverse=True)
        best_umi, best_umi_score = sorted_scores[0]
        delta_umi_score = best_umi_score - sorted_scores[1][1] if len(sorted_scores) > 1 else None
        scores["Best UMI"] = best_umi
        scores["Best UMI Score"] = best_umi_score
        scores["Delta UMI Score"] = delta_umi_score

    return scores


def parse_adapter3_seq(
    sequence: str, primer: str, umi_dict: dict, max_align_score: int = 2
):  # Restrict max edits to 2 for adapter
    """
    Identify sequence characteristics in the 3' adapter sequence, including any potential UMIs
    """
    if sequence is None:
        primer_score = None
        putative_umi = None
    else:
        primer_alignment = align(primer, sequence, max_score=max_align_score)
        primer_locations = primer_alignment["locations"]
        # primer_score = primer_alignment['editDistance']
        # print( sequence, primer, primer_alignment )
        if len(primer_locations) > 0:
            primer_score = primer_alignment[
                "editDistance"
            ]  # (len(primer)-primer_alignment['editDistance'])/len(primer)
            post_start = np.min([loc[0] for loc in primer_locations])
            putative_umi = sequence[:post_start]
        else:
            primer_score = None  # 0.0 #0.5/len(primer)
            # post_start = len(sequence) # Retain full sequence for UMI alignment
            putative_umi = None
    umi_scores = score_umi(putative_umi, umi_dict)
    return primer_score, umi_scores


def compute_read_statistics(read: dict):
    read_length = len(read["sequence"])
    mean_quality = np.mean([ord(x) for x in read["QUAL"]]) - 33.0
    read_stats = {"read_length": read_length, "mean_quality": mean_quality}
    return read_stats


def annotate_read(read: dict, ssp: str, primer3: str, umis: dict, min_transcript_len: int = 200):
    """
    Search a read for a PolyA tail and then score the 5'/3' ends accordingly

    *** POSSIBLE PATCH FOR LATER ***
    Currently max_polyA_score_len is arbitrary, consider changing later
    ***

    """

    ## Read sequence, identify polyA tail if present, process 3' and 5' ends
    sequence = read["sequence"]
    read_stats = compute_read_statistics(read)
    pre_polyA, post_polyA, polyA_stats = split_sequence_by_polyA(sequence)
    primer_3_score, umi_scores = parse_adapter3_seq(post_polyA, primer3, umis)
    transcript_seq, adapter5_seq, primer_5_score, primer5_trim_idx = parse_adapter5_seq(pre_polyA, ssp)

    # Construct sample label
    sample_label = f"{read['Source']}_{umi_scores['Best UMI']}"

    ## Analyze the identified transcript sequence
    if transcript_seq is not None:
        transcript_len = len(transcript_seq)
        if primer5_trim_idx is not None:
            transcript_qual = read["QUAL"][primer5_trim_idx : primer5_trim_idx + transcript_len]
        else:
            transcript_qual = None
    else:
        transcript_len = None
        transcript_qual = None

    t_scores = np.array([primer_5_score, primer_3_score, umi_scores["Best UMI Score"]])
    if np.all(t_scores == None):  # noqa: E711 - numpy element-wise comparison
        cDNA_status = "Unknown"
    else:  # Something was detected
        if np.all(t_scores != None):  # noqa: E711 - numpy element-wise comparison
            cDNA_status = "Complete"
        else:  # More complicated
            if np.all(t_scores[1:] != None):  # noqa: E711 - numpy element-wise comparison
                cDNA_status = "3' Only"
            elif np.all(t_scores[:2] != None):  # noqa: E711 - numpy element-wise comparison
                cDNA_status = "Missing Barcode"
            elif np.all(t_scores[1:] == None):  # noqa: E711 - numpy element-wise comparison
                cDNA_status = "5' Only"
            else:  # Not sure what this would be...
                cDNA_status = "Unknown"
    if transcript_len is not None and transcript_len < min_transcript_len:
        cDNA_status += " Fragment"

    ## Collate the annotations
    read_annotations = {
        "Sample_ID": sample_label,
        "3' primer": primer3,
        "3' primer alignment score": primer_3_score,
        "3' adapter sequence": post_polyA,
        "Strand switch primer": ssp,
        "SSP alignment score": primer_5_score,
        "5' adapter trim": primer5_trim_idx,
        "5' adapter sequence": adapter5_seq,
        "Transcript sequence": transcript_seq,
        "Transcript quality": transcript_qual,
        "Transcript length": transcript_len,
        "cDNA status": cDNA_status,
    }

    ## Add ORF search results
    transcript_for_orf = transcript_seq if cDNA_status in ["Complete", "3' Only"] else None
    best_orf = return_best_orf(transcript_for_orf)

    return read | read_stats | polyA_stats | umi_scores | read_annotations | best_orf


def process_read(
    read: dict,
    ssp: str,
    primer3: str,
    umis: dict,
    min_polyA_score: float = 0.0,
    min_ssp_score: float = 0.5,
    min_umi_score: float = 0.5,
):
    """
    Annotate reads and, if possible, re-orient along forward strand

    *** POSSIBLE PATCH FOR LATER ***
    Currently min_polyA_score and min_ssp_score are *very* arbitrarily chosen
    based on examination of real data; usually real scores are substantially
    better than 0.5, but seemed like a good enough starting point. Consider
    changing based on further testing.
    ***
    """

    # Analyze read in both forward and reverse directions
    f_analysis = annotate_read(read, ssp, primer3, umis)
    rc_read = dict(read)
    rc_read["sequence"] = rc(rc_read["sequence"])
    r_analysis = annotate_read(rc_read, ssp, primer3, umis)

    # Select forward or reverse direction based on analysis
    cDNA = np.array([a["cDNA status"] for a in [f_analysis, r_analysis]])
    if np.all((cDNA == "Unknown") | (cDNA == "Unknown Fragment")):
        analysis = f_analysis  # Default to forward read, since we're guessing
    else:
        if np.any((cDNA == "Unknown") | (cDNA == "Unknown Fragment")):  # Only one good option
            analysis = f_analysis if cDNA[1] == "Unknown" else r_analysis
        else:  # Multiple options
            if np.any(cDNA == "Complete"):
                # Assume only one possible good read here
                analysis = f_analysis if cDNA[0] == "Complete" else r_analysis
            else:
                # Return forward analysis with a modified label
                analysis = f_analysis
                analysis["cDNA status"] = "Complex"  # f"{f_analysis['cDNA status']}_{r_analysis['cDNA status']}"

    # Encode movemap if present in the read data
    if "mv:B" in analysis:
        analysis["ts"] = np.int32(analysis["ts:i"])
        movemap = analysis["mv:B"].split(",")
        analysis["stride"] = np.int16(movemap[1])
        analysis["movemap"] = np.asarray(movemap[2:], "uint8")

    return analysis
