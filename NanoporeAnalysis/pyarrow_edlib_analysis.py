import time
import os
from pathlib import Path
import math
import concurrent
import concurrent.futures
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq
from pyarrow import dataset
import matplotlib.pyplot as plt
import numpy as np
from NanoporeAnalysis import utils
from NanoporeAnalysis import local_io
import edlib

def pod5_split() :
    
    #generate the index of all reads and their channels
    if not os.path.isfile(Path(path_out + "/view.txt")) or reset_pod5_view == True :
        pod5_view.view_pod5([Path(path_data)], Path(path_out), include = "read_id, channel", force_overwrite=True, threads=threads)
    else :
        view = pd.read_table(Path(path_out + "/view.txt"), sep='\t')
    
    channels = view['channel'].unique()
    channels_per_worker = np.array_split(np.array(channels), workers)
    
    for i in range(workers) :
        
        view.query("channel in " + str(channels_per_worker[i]))['read_id'].to_csv(Path(path_out + "/view_current.txt"), index = False, sep=' ', header=False)
        pod5_filter.filter_pod5([Path(path_data)], Path(path_out + '/split_pod5s_' + str(i)), threads = threads, force_overwrite = True, ids = Path(path_out + "/view_current.txt"))
    return

def sam_to_parquet(file, path_dataset, basename_template = None, sam_or_bam = 'sam') :
    table_dict = {
        'ID' : [],
        'seq' : [],
        'seq_len' : [],
        'qual' : []
    }
    table_dict_keys = table_dict.keys()
    tag_names = []
    with open(file, 'r') as handle :
        i = 0
        for line in handle.readlines() :
            if line[0] != '@' :
                line_split = line.split('\t')
                line_dict = dict.fromkeys(table_dict_keys)
                
                line_dict['ID'] = line_split[0]
                line_dict['seq'] = line_split[9]
                line_dict['seq_len'] = len(line_dict['seq'])
                line_dict['qual'] = line_split[10]
                
                for tag in line_split[11:] :
                    tag_name = tag[:4]
                    tag_content = tag[5:]
                    line_dict[tag_name] = tag_content
                
                for key in line_dict.keys() :
                    if key not in table_dict_keys :
                        table_dict[key] = pa.nulls(i).to_pylist()
                    table_dict[key].append(line_dict[key])
                i+=1
    table = pa.table(table_dict)
    dataset.write_dataset(table, Path(path_out + '/pa_dataset'), format='parquet', basename_template = basename_template, max_rows_per_file = 200000, max_rows_per_group = 200000, existing_data_behavior='overwrite_or_ignore')
    return 'done'

def build_parquet_dataset_from_sam(sam_folder, path_dataset) :
    files = [x for x in Path(sam_folder).iterdir() if x.is_file()]
    for file in files :
        sam_to_parquet(file, path_dataset, basename_template = str(file.stem + '_part-{i}.parquet'))
    return 'done'

def find_seq_matches(target_seq, query_seq, max_edit_distance, query_ID, min_length = 0, skip_reverse = False) :
    
    matches = []
    target_seq_len = len(target_seq)
    query_seq_len = len(query_seq)
    for_alignment = edlib.align(query_seq, target_seq, mode='HW', task='locations')
    for location in for_alignment['locations'] :
        alignment_length = abs(location[0] - location[1])
        if for_alignment['editDistance'] <= max_edit_distance and alignment_length >= min_length :
            score = ( for_alignment['editDistance'] + query_seq_len - alignment_length ) / query_seq_len
            matches.append({
                'query_ID' : query_ID,
                'edit_distance' : for_alignment['editDistance'],
                'edit_score' : score,
                'location' : location,
                'direction' : 'forward'
            })
    if not skip_reverse :
        rev_alignment = edlib.align(query_seq, utils.reverse_complement(target_seq), mode='HW', task='locations')
        for location in rev_alignment['locations'] :
            alignment_length = abs(location[0] - location[1])
            if rev_alignment['editDistance'] <= max_edit_distance and alignment_length >= min_length :
                score = ( rev_alignment['editDistance'] + query_seq_len - alignment_length ) / query_seq_len
                matches.append({
                    'query_ID' : query_ID,
                    'edit_distance' : rev_alignment['editDistance'],
                    'edit_score' : score,
                    'location' : location,
                    'direction' : 'reverse'
                })
    return matches

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
    alignment = edlib.align(query, sequence, mode='HW', task='locations')
    align_indices = alignment['locations']
    if len(align_indices) == 0:
        return []
    else:
        merged_indices = merge_overlapped_indices(align_indices, tolerated_mismatches)
        filtered_indices = [ i for i in merged_indices if i['length'] >= min_len]
        sorted_indices = sorted(filtered_indices, key=lambda x:x['length'], reverse=True)
        return sorted_indices
    
def find_polyA_bidirectional(sequence: str,
                             N: int,
                             tolerated_mismatches: int,
                             min_len: int) -> list:
    for_indices = find_polyX(sequence, 'A', N, tolerated_mismatches, min_len)
    rev_indices = find_polyX(utils.reverse_complement(sequence), 'A', N, tolerated_mismatches, min_len)
    indices = []
    for index_set in for_indices :
        index_set['direction'] = 'forward'
        indices.append(index_set)
    for index_set in rev_indices :
        index_set['direction'] = 'reverse'
        indices.append(index_set)
    return indices

def pick_best_match(matches, min_match_location = [0, 0]) :
    best_match = {
        'query_ID' : 'None matched',
        'edit_distance' : -1,
        'edit_score' : 1000,
        'location' : [-1,-1],
        'direction' : 'None'
    }
    for match in matches :
        if min_match_location[0] != 0 :
            match_start = min_match_location[1] - 1  - match['location'][1] if match['direction'] == 'reverse' else match['location'][0]
        if match_start >= min_match_location[0] :
            if match['edit_score'] < best_match['edit_score'] :
                best_match = match
            elif match['edit_score'] == best_match['edit_score'] :
                best_match = {
                    'query_ID' : 'Multiple',
                    'edit_distance' : 1000,
                    'edit_score' : match['edit_score'],
                    'location' : [-1, -1],
                    'direction' : 'None'
                }
    return best_match

def evaluate_polyA_UMI_SSP(polyA, UMI_match, SSP_match, polyA_UMI_SSP_set, max_gap = 5) :
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
        polyA_UMI_SSP_set['direction'] = UMI_match['direction']
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
        polyA_UMI_SSP_set['direction'] = UMI_match['direction']
        polyA_UMI_SSP_set['barcode_flag'][1] = 1
        polyA_UMI_SSP_set['barcode_flag'][2] = 1
        polyA_UMI_SSP_set['barcode_score'] = polyA_UMI_SSP_set['barcode_UMI_edit_distance'] + polyA_UMI_SSP_set['barcode_SSP_edit_distance']
    return polyA_UMI_SSP_set
    
def reverse_matches(matches, seq_len) :
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

def parse_polyA_UMI_SSP(seq, UMIs, SSP, UMI_max_score, SSP_max_score, max_gap = 5, UMI_min_len = 0, SSP_min_len = 0, polyA_N = 4, polyA_tolerated_mismatches = 1, polyA_min_len = 10) :
#     seq_len = len(seq)
#     polyAs = find_polyA_bidirectional(seq, N = polyA_N, tolerated_mismatches = polyA_tolerated_mismatches, min_len = polyA_min_len)
#     UMI_matches = []
#     SSP_matches = find_seq_matches(seq, SSP, SSP_max_score, 'SSP', min_length = SSP_min_len)
#     UMI_seq_regions = []
#     for polyA in polyAs :
#         region = [ polyA['end'] , min(polyA['end'] + 50, seq_len - 1) , polyA['direction'] ]
#         if region[1] - region[0] > 0 :
#             UMI_seq_regions.append(region)
#     for SSP_match in SSP_matches :
#         region = [ max(SSP_match['location'][0] - 50, 0) , SSP_match['location'][0] , SSP_match['direction'] ]
#         if region[1] - region[0] > 0 :
#             UMI_seq_regions.append(region)
#     if len(SSP_matches) == 0 :
#         SSP_matches = [{ 'location' : [-1,-1], 'direction' : 'None', 'edit_distance' : 1000 }]
#     for UMI in UMIs :
#         for region in UMI_seq_regions :
#             if region[2] == 'forward' :
#                 tmp_UMI_matches = find_seq_matches(seq[region[0] : region[1]], UMI[1], UMI_max_score, UMI[0], min_length = UMI_min_len, skip_reverse = True)
#             else :
#                 tmp_UMI_matches = find_seq_matches(utils.reverse_complement(seq)[region[0] : region[1]], UMI[1], UMI_max_score, UMI[0], min_length = UMI_min_len, skip_reverse = True)
#             if len(tmp_UMI_matches) > 0 :
#                 for UMI_match in tmp_UMI_matches :
#                     UMI_match['location'] = [ UMI_match['location'][0] + region[0], UMI_match['location'][1] + region[0] ]
#                     UMI_matches.append(UMI_match)
    seq_len = len(seq)
    polyAs = find_polyA_bidirectional(seq, N = polyA_N, tolerated_mismatches = polyA_tolerated_mismatches, min_len = polyA_min_len)
    SSP_matches = find_seq_matches(seq, SSP, SSP_max_score, 'SSP', min_length = SSP_min_len)
    if len(SSP_matches) == 0 :
        SSP_matches = [{ 'location' : [-1,-1], 'direction' : 'None', 'edit_distance' : 1000 }]
    UMI_matches = []
    for UMI in UMIs :
        tmp_UMI_matches = find_seq_matches(seq, UMI[1], UMI_max_score, UMI[0], min_length = UMI_min_len)
        if len(tmp_UMI_matches) > 0 :
            for UMI_match in tmp_UMI_matches :
                UMI_matches.append(UMI_match)
    polyA_UMI_SSP_default = {
        'barcode_ID' : 'None matched',
        'barcode_polyA_start' : -1,
        'barcode_polyA_end' : -1,
        'barcode_polyA_len' : -1,
        'barcode_UMI_start' : -1,
        'barcode_UMI_end' : -1,
        'barcode_UMI_edit_distance' : -1,
        'barcode_SSP_start' : -1,
        'barcode_SSP_end' : -1,
        'barcode_SSP_edit_distance' : -1,
        'barcode_score' : 1000,
        'direction' : 'None',
        'barcode_flag' : [0,0,0]
    }
    polyA_UMI_SSP_set = polyA_UMI_SSP_default.copy()
    for UMI_match in UMI_matches :
        for SSP_match in SSP_matches :
            for polyA in polyAs :
                barcode_evaluation = evaluate_polyA_UMI_SSP(polyA, UMI_match, SSP_match, polyA_UMI_SSP_default, max_gap)
                if barcode_evaluation['barcode_flag'][1] == 1 :
                    if barcode_evaluation['barcode_score'] < polyA_UMI_SSP_set['barcode_score'] :
                        polyA_UMI_SSP_set = barcode_evaluation
                    elif barcode_evaluation['barcode_score'] == polyA_UMI_SSP_set['barcode_score'] :
                        if barcode_evaluation['barcode_polyA_len'] > polyA_UMI_SSP_set['barcode_polyA_len'] :
                            polyA_UMI_SSP_set = barcode_evaluation
                        elif barcode_evaluation['barcode_ID'] != polyA_UMI_SSP_set['barcode_ID'] :
                            polyA_UMI_SSP_set['barcode_ID'] = 'Multiple'
                            polyA_UMI_SSP_set['barcode_score'] = barcode_evaluation['barcode_score']
#     if polyA_UMI_SSP_set['barcode_ID'] == 'None matched' and len(UMI_matches) > 0 :
#         best_UMI = pick_best_match(UMI_matches, [ 0.8 * seq_len, seq_len ])
#         polyA_UMI_SSP_set['barcode_ID'] = best_UMI['query_ID']
#         polyA_UMI_SSP_set['barcode_UMI_start'] = best_UMI['location'][0]
#         polyA_UMI_SSP_set['barcode_UMI_end'] = best_UMI['location'][1]
#         polyA_UMI_SSP_set['barcode_UMI_edit_distance'] = best_UMI['edit_distance']
#         polyA_UMI_SSP_set['direction'] = best_UMI['direction']
        
    if SSP_matches[0]['direction'] != 'None' :
        SSP_matches = reverse_matches(SSP_matches, seq_len)
    if polyA_UMI_SSP_set['direction'] == 'reverse' :
        seq = utils.reverse_complement(seq)
        tmp_UMI_start = polyA_UMI_SSP_set['barcode_UMI_start']
        polyA_UMI_SSP_set['barcode_UMI_start'] = seq_len - 1 - polyA_UMI_SSP_set['barcode_UMI_end']
        polyA_UMI_SSP_set['barcode_UMI_end'] = seq_len - 1 - tmp_UMI_start
        if polyA_UMI_SSP_set['barcode_SSP_start'] != -1 :
            polyA_UMI_SSP_set['barcode_SSP_start'] = seq_len - 1 - polyA_UMI_SSP_set['barcode_SSP_end']
            polyA_UMI_SSP_set['barcode_SSP_end'] = seq_len - 1 - polyA_UMI_SSP_set['barcode_SSP_start']
    distal_SSP = {
        'SSP_start' : -1,
        'SSP_end' : -1,
        'SSP_edit_distance' : 1000,
    }
    max_SSP_start = polyA_UMI_SSP_set['barcode_UMI_start'] if polyA_UMI_SSP_set['barcode_UMI_start'] != -1 else seq_len - 1
    for SSP_match in SSP_matches :
        if SSP_match['direction'] != polyA_UMI_SSP_set['direction'] and SSP_match['location'][1] <= max_SSP_start :
            if SSP_match['edit_distance'] < distal_SSP['SSP_edit_distance'] :
                distal_SSP['SSP_start'] = SSP_match['location'][0]
                distal_SSP['SSP_end'] = SSP_match['location'][1]
                distal_SSP['SSP_edit_distance'] = SSP_match['edit_distance']
            elif SSP_match['edit_distance'] == distal_SSP['SSP_edit_distance'] and SSP_match['location'][1] < distal_SSP['SSP_end'] :
                distal_SSP['SSP_start'] = SSP_match['location'][0]
                distal_SSP['SSP_end'] = SSP_match['location'][1]
                distal_SSP['SSP_edit_distance'] = SSP_match['edit_distance']
    
    parsed = polyA_UMI_SSP_set
    parsed.update(distal_SSP)
    if parsed['barcode_UMI_start'] - parsed['SSP_end'] > seq_len * 0.5 and parsed['SSP_end'] != -1 and parsed['barcode_flag'][1] == 1 :
        parsed['biological_seq_indices'] = [ parsed['SSP_end'], parsed['barcode_UMI_start'] ]
    else :
        parsed['biological_seq_indices'] = [ 0, seq_len - 1 ]
    return parsed

def debarcode_table(table, UMIs, SSP, UMI_max_score, SSP_max_score, max_gap = 5, UMI_min_len = 0, SSP_min_len = 0, polyA_N = 4, polyA_tolerated_mismatches = 1, polyA_min_len = 10) :
    seqs = table.column('seq').to_pylist()
    parsed_seqs = {
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
        'barcode_score' : [],
        'direction' : [],
        'barcode_flag' : [],
        'SSP_start' : [],
        'SSP_end' : [],
        'SSP_edit_distance' : [],
        'biological_seq_indices' : []
    }
    for i in range(len(seqs)) :
        parsed_seq = parse_polyA_UMI_SSP(seqs[i], UMIs, SSP, UMI_max_score, SSP_max_score, max_gap, UMI_min_len, SSP_min_len, polyA_N, polyA_tolerated_mismatches, polyA_min_len)
        for key in parsed_seqs :
            parsed_seqs[key].append(parsed_seq[key])
    for key in parsed_seqs :
        table = table.append_column(key, [parsed_seqs[key]])
    return table

def debarcode_table_from_file(file, UMIs, SSP, UMI_max_score, SSP_max_score, max_gap = 5, UMI_min_len = 0, SSP_min_len = 0, polyA_N = 4, polyA_tolerated_mismatches = 1, polyA_min_len = 10, resume=False, overwrite=False) :
    table = pq.read_table(file)
    if 'barcode_ID' in table.column_names :
        if resume == True :
            print('skip resume', file)
            print(table.column_names)
            del table
            return
        elif overwrite == True :
            for column_name in [
                'barcode_ID',
                'barcode_polyA_start',
                'barcode_polyA_end',
                'barcode_polyA_len',
                'barcode_UMI_start',
                'barcode_UMI_end',
                'barcode_UMI_edit_distance',
                'barcode_SSP_start',
                'barcode_SSP_end',
                'barcode_SSP_edit_distance',
                'barcode_score',
                'direction',
                'barcode_flag',
                'SSP_start',
                'SSP_end',
                'SSP_edit_distance',
                'biological_seq_indices'
            ] :
                if column_name in table.column_names :
                    table = table.drop_columns([column_name])
        else :
            print('skip no overwrite', file)
            del table
            return
    table = debarcode_table(table, UMIs, SSP, UMI_max_score, SSP_max_score, max_gap, UMI_min_len, SSP_min_len, polyA_N, polyA_tolerated_mismatches, polyA_min_len)
    pq.write_table(table, file)
    del table
    print(file)
    return

def load_barcodes(path) :
    barcodes = []
    with open(path, 'r') as handle :
        for line in handle.readlines() :
            line_split = line.split(',')
            barcodes.append([line_split[0], utils.reverse_complement(line_split[2])])
    return barcodes

def debarcode(dataset_dir, UMIs_path, SSP, UMI_max_score, SSP_max_score, max_gap = 5, UMI_min_len = 0, SSP_min_len = 0, polyA_N = 4, polyA_tolerated_mismatches = 1, polyA_min_len = 10, workers = 4, resume=False, overwrite=False ) :
    UMIs = load_barcodes(UMIs_path)
    files = [x for x in Path(dataset_dir).iterdir() if x.is_file()]
    split = math.ceil( len(files) / workers )
    files_chunks = np.array_split(np.array(files), split)
    for chunk in files_chunks : 
        with concurrent.futures.ProcessPoolExecutor( max_workers=workers ) as executor :
            futures = [ executor.submit( debarcode_table_from_file, file, UMIs, SSP, UMI_max_score, SSP_max_score, max_gap, UMI_min_len, SSP_min_len, polyA_N, polyA_tolerated_mismatches, polyA_min_len, resume, overwrite ) for file in chunk]
            concurrent.futures.wait( futures )
    return

def minimap2_table(table, path_ref, preset='splice') :
    seqs = table.column('seq').to_pylist()
    biological_seq_indices = table.column('biological_seq_indices').to_pylist()
    combined_scores = table.column('combined_score').to_pylist()
    SSP_edit_distances = table.column('SSP_edit_distance').to_pylist()
    barcode_IDs = table.column('barcode_ID').to_pylist()
    aligner = mappy.Aligner(path_ref, preset=preset)
    alignments_core = {
        'minimap2_q_st' : [],
        'minimap2_q_en' : [],
        'minimap2_strand' : [],
        'minimap2_ctg' : [],
        'minimap2_ctg_len' : [],
        'minimap2_r_st' : [],
        'minimap2_r_en' : [],
        'minimap2_mlen' : [],
        'minimap2_blen' : [],
        'minimap2_mapq' : []
    }
    alignments_tags = {}
    alignments_core_keys = alignments_core.keys()
    alignments_tags_keys = alignments_tags.keys()
    i=0
    for i in range(len(seqs)) :
        hit_core_dict = dict.fromkeys(alignments_core_keys)
        hit_tags_dict = dict.fromkeys(alignments_tags_keys)
        if combined_scores[i] <= 20 and barcode_IDs[i] != 'Multiple' and SSP_edit_distances[i] <= 5 :
            indices = biological_seq_indices[i]
            for hit in aligner.map(seqs[i][indices[0] : indices[1]]) :
                if hit.is_primary :
                    hit_split = str(hit).split('\t')
                    j = 0
                    for key in alignments_core_keys :
                        hit_core_dict[key] = hit_split[j] 
                        j += 1
                    for tag in hit_split[j:] :
                        tag_name = 'minimap2_' + tag[:4]
                        tag_content = tag[5:]
                        hit_tags_dict[tag_name] = tag_content
        
        for key in alignments_core_keys :
            alignments_core[key].append(hit_core_dict[key])
        for key in hit_tags_dict.keys() :
            if key not in alignments_tags_keys :
                alignments_tags[key] = pa.nulls(len(alignments_core['minimap2_q_st']) - 1).to_pylist()
            alignments_tags[key].append(hit_tags_dict[key])
    for key in alignments_core :
        table = table.append_column(key, [alignments_core[key]])
    for key in alignments_tags :
        table = table.append_column(key, [alignments_tags[key]])
    del aligner
    return table

def minimap2_table_from_file(file, path_ref, preset='splice', resume=False, overwrite=False) :
    table = pq.read_table(file)
    if 'minimap2_q_st' in table.column_names :
        if resume == True :
            print('skip resume', file)
            del table
            return
        elif overwrite == True :
            columns_to_drop = [ name for name in table.column_names if 'minimap2' in name ]
            table = table.drop_columns(columns_to_drop)
        else :
            print('skip no overwrite', file)
            del table
            return
    table = minimap2_table(table, path_ref, preset=preset)
    pq.write_table(table, file)
    del table
    print(file)
    return

def minimap2(dataset_dir, path_ref, preset='splice', resume=False, overwrite=False, workers = 4) :
    files = [x for x in Path(dataset_dir).iterdir() if x.is_file()]
    split = math.ceil( len(files) / workers )
    files_chunks = np.array_split(np.array(files), split)
    for chunk in files_chunks :
        with concurrent.futures.ProcessPoolExecutor( max_workers=workers ) as executor :
            futures = [ executor.submit( minimap2_table_from_file, file, path_ref, preset=preset, resume=resume, overwrite=overwrite ) for file in chunk]
            concurrent.futures.wait( futures )
    return

def count_mapped_reads(dataset_dir, path_ref, path_out_csv, sample_dict) :
    ref = local_io.read_fastx(path_ref)
    contigs_to_gene_name = {}
    gene_names = []
    
    for key in ref.keys() :
        comma_split = key.split(', ')
        contig_name_split = comma_split[0].split(' ')
        contig_id = contig_name_split[0]
        gene_name = ' '.join(contig_name_split[1:])
        contigs_to_gene_name[contig_id] = gene_name
        gene_names.append(gene_name)
    
    table = pq.read_table(dataset_dir, columns=['barcode_ID', 'minimap2_ctg', 'minimap2_mapq'])
    barcodes = table.column('barcode_ID').unique().to_pylist()
    counts = False
    for barcode in barcodes :
        if barcode not in ['None matched', 'Multiple'] :
            counts_in_barcode_array = table.filter( pc.field('barcode_ID') == barcode ).column('minimap2_ctg').drop_null().value_counts()
            counts_in_barcode_dict = dict.fromkeys(gene_names, 0)
            total_counts = 0
            for count in counts_in_barcode_array :
                gene_name = contigs_to_gene_name[str(count[0])]
                counts_in_barcode_dict[gene_name] += count[1].as_py()
                total_counts += count[1].as_py()
            for key in counts_in_barcode_dict :
                counts_in_barcode_dict[key] = [ counts_in_barcode_dict[key] ]#/ total_counts ]
            
            counts_in_barcode_table = pa.table(counts_in_barcode_dict).add_column(0, 'barcode_ID', [[barcode]]).add_column(0, 'sample', [[sample_dict[barcode]]])
            if not counts :
                counts = counts_in_barcode_table
            else :
                counts = pa.concat_tables([counts, counts_in_barcode_table])
    
    csv.write_csv(counts, path_out_csv)
    print('done')
    return

def qc_metrics(path_dataset) :
    table = pq.read_table(path_dataset, columns = ['seq_len', 'combined_score', 'barcode_ID', 'SSP_edit_distance'])
    qc = {}
    hist_max = 5000
    fig1, qc['alignment_scores_hist'] = plt.subplots(2,2, sharex=True, sharey=True)

    group_1 = table.filter(pc.field("combined_score") <= 20).filter(pc.field('barcode_ID') != 'Multiple').filter(pc.field('SSP_edit_distance') <= 5)
    qc['alignment_scores_hist'][0,0].hist(group_1.column('seq_len').to_pylist(), bins=50, range=(0, hist_max))
    qc['alignment_scores_hist'][0,0].set_title('Passed Barcode and Strand Switch', {'fontsize': 10})
    group_2 = table.filter(pc.field("combined_score") <= 20).filter(pc.field('barcode_ID') != 'Multiple').filter(pc.field('SSP_edit_distance') > 5)
    qc['alignment_scores_hist'][0,1].hist(group_2.column('seq_len').to_pylist(), bins=50, range=(0, hist_max))
    qc['alignment_scores_hist'][0,1].set_title('Passed Barcode, failed Strand Switch', {'fontsize': 10})
    group_3 = table.filter(pc.field("combined_score") > 20).filter(pc.field('SSP_edit_distance') <= 5)
    qc['alignment_scores_hist'][1,0].hist(group_3.column('seq_len').to_pylist(), bins=50, range=(0, hist_max))
    qc['alignment_scores_hist'][1,0].set_title('Failed Barcode, passed Strand Switch', {'fontsize': 10})
    group_4 = table.filter(pc.field("combined_score") > 20).filter(pc.field('SSP_edit_distance') > 5)
    qc['alignment_scores_hist'][1,1].hist(group_4.column('seq_len').to_pylist(), bins=50, range=(0, hist_max))
    qc['alignment_scores_hist'][1,1].set_title('Failed Barcode and Strand Switch', {'fontsize': 10})

    barcodes = table.column('barcode_ID').unique().to_pylist()
    for barcode in barcodes :
        print(barcode + ' passed barcode and SSP: ' + str(group_1.filter(pc.field('barcode_ID') == barcode).num_rows))
        print(barcode + ' passed barcode, failed SSP: ' + str(group_2.filter(pc.field('barcode_ID') == barcode).num_rows))
        print(barcode + ' failed barcode, passed SSP: ' + str(group_3.filter(pc.field('barcode_ID') == barcode).num_rows))
        print(barcode + ' failed barcode and SSP: ' + str(group_4.filter(pc.field('barcode_ID') == barcode).num_rows))
    return qc