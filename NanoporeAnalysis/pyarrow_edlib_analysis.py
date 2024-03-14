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

def find_seq_matches(target_seq, query_seq, max_edit_distance, query_ID, min_length=0) :
    
    matches = []
    target_seq_len = len(target_seq)
    query_seq_len = len(query_seq)
    for_alignment = edlib.align(query_seq, target_seq, mode='HW', task='path')
    rev_alignment = edlib.align(utils.reverse_complement(query_seq), target_seq, mode='HW', task='path')
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

def evaluate_barcode_SSP_pair(barcode_match, SSP_match, max_gap = 5) :
    barcode_start = barcode_match['location'][0]
    barcode_end = barcode_match['location'][1]
    SSP_start = SSP_match['location'][0]
    SSP_end = SSP_match['location'][1]
    barcode_SSP_pair = None
    pair_gap = abs(SSP_start - barcode_end)
    if pair_gap <= max_gap and barcode_match['direction'] == SSP_match['direction'] :
        barcode_SSP_pair = {
            'barcode_ID' : barcode_match['query_ID'],
            'barcode_UMI_start' : barcode_start,
            'barcode_UMI_end' : barcode_end,
            'barcode_UMI_edit_distance' : barcode_match['edit_distance'],
            'barcode_SSP_start' : SSP_start,
            'barcode_SSP_end' : SSP_end,
            'barcode_SSP_edit_distance' : SSP_match['edit_distance'],
            'combined_score' : barcode_match['edit_score'] + SSP_match['edit_score'],
            'direction' : barcode_match['direction']
        }
    return barcode_SSP_pair

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

def parse_SSP_and_barcode(seq, barcodes, SSP, max_barcode_score, max_SSP_score, max_gap = 5, min_length_barcode = 0, min_length_SSP = 0) :
    seq_len = len(seq)
    SSP_matches = find_seq_matches(seq, SSP, max_SSP_score, 'SSP', min_length = min_length_SSP)
    barcode_matches = []
    for barcode in barcodes :
        tmp_barcode_matches = find_seq_matches(seq, barcode[1], max_barcode_score, barcode[0], min_length = min_length_barcode)
        if len(tmp_barcode_matches) > 0 :
            for barcode_match in tmp_barcode_matches :
                barcode_matches.append(barcode_match)
                
    barcode_SSP_pair_default = {
        'barcode_ID' : 'None matched',
        'barcode_UMI_start' : -1,
        'barcode_UMI_end' : -1,
        'barcode_UMI_edit_distance' : -1,
        'barcode_SSP_start' : -1,
        'barcode_SSP_end' : -1,
        'barcode_SSP_edit_distance' : -1,
        'combined_score' : 1000,
        'direction' : 'None'
    }
    barcode_SSP_pair = barcode_SSP_pair_default
    for barcode_match in barcode_matches :
        for SSP_match in SSP_matches :
            barcode_SSP_evaluation = evaluate_barcode_SSP_pair(barcode_match, SSP_match, max_gap)
            if barcode_SSP_evaluation != None :
                if barcode_SSP_evaluation['combined_score'] < barcode_SSP_pair['combined_score'] :
                    barcode_SSP_pair = barcode_SSP_evaluation
                elif barcode_SSP_evaluation['combined_score'] == barcode_SSP_pair['combined_score'] :
                    barcode_SSP_pair = barcode_SSP_pair_default
                    barcode_SSP_pair['barcode_ID'] = 'Multiple'
                    barcode_SSP_pair['combined_score'] = barcode_SSP_evaluation['combined_score']
    
    if barcode_SSP_pair['barcode_ID'] == 'None matched' and len(barcode_matches) > 0 :
        best_barcode = pick_best_match(barcode_matches, [ 0.8 * seq_len, seq_len ])
        barcode_SSP_pair['barcode_ID'] = best_barcode['query_ID']
        barcode_SSP_pair['barcode_UMI_start'] = best_barcode['location'][0]
        barcode_SSP_pair['barcode_UMI_end'] = best_barcode['location'][1]
        barcode_SSP_pair['barcode_UMI_edit_distance'] = best_barcode['edit_distance']
        barcode_SSP_pair['direction'] = best_barcode['direction']
        
    if barcode_SSP_pair['direction'] == 'reverse' :
        seq = utils.reverse_complement(seq)
        SSP_matches = reverse_matches(SSP_matches, seq_len)
        tmp_barcode_UMI_start = barcode_SSP_pair['barcode_UMI_start']
        barcode_SSP_pair['barcode_UMI_start'] = seq_len - 1 - barcode_SSP_pair['barcode_UMI_end']
        barcode_SSP_pair['barcode_UMI_end'] = seq_len - 1 - tmp_barcode_UMI_start
        if barcode_SSP_pair['barcode_SSP_start'] != -1 :
            barcode_SSP_pair['barcode_SSP_start'] = seq_len - 1 - barcode_SSP_pair['barcode_SSP_end']
            barcode_SSP_pair['barcode_SSP_end'] = seq_len - 1 - barcode_SSP_pair['barcode_SSP_start']
    
    distal_SSP = {
        'SSP_start' : -1,
        'SSP_end' : -1,
        'SSP_edit_distance' : 1000,
    }
    max_SSP_start = barcode_SSP_pair['barcode_UMI_start'] if barcode_SSP_pair['barcode_UMI_start'] != -1 else seq_len - 1
    for SSP_match in SSP_matches :
        if SSP_match['direction'] != barcode_SSP_pair['direction'] and SSP_match['location'][1] <= max_SSP_start :
            if SSP_match['edit_distance'] < distal_SSP['SSP_edit_distance'] :
                distal_SSP['SSP_start'] = SSP_match['location'][0]
                distal_SSP['SSP_end'] = SSP_match['location'][1]
                distal_SSP['SSP_edit_distance'] = SSP_match['edit_distance']
            elif SSP_match['edit_distance'] == distal_SSP['SSP_edit_distance'] and SSP_match['location'][1] < distal_SSP['SSP_end'] :
                distal_SSP['SSP_start'] = SSP_match['location'][0]
                distal_SSP['SSP_end'] = SSP_match['location'][1]
                distal_SSP['SSP_edit_distance'] = SSP_match['edit_distance']
    
    parsed = barcode_SSP_pair
    parsed.update(distal_SSP)
    if parsed['barcode_UMI_start'] - parsed['SSP_end'] > seq_len * 0.5 and parsed['barcode_SSP_start'] != -1 :
        parsed['biological_seq_indices'] = [ parsed['SSP_end'], parsed['barcode_UMI_start'] ]
    else :
        parsed['biological_seq_indices'] = [ 0, seq_len - 1 ]
    return parsed

def debarcode_table(table, barcodes, SSP, max_barcode_score, max_SSP_score, max_gap=5, min_length_barcode=0, min_length_SSP=0) :
    seqs = table.column('seq')
    parsed_seqs = {
        'barcode_ID' : [],
        'barcode_UMI_start' : [],
        'barcode_UMI_end' : [],
        'barcode_UMI_edit_distance' : [],
        'barcode_SSP_start' : [],
        'barcode_SSP_end' : [],
        'barcode_SSP_edit_distance' : [],
        'combined_score' : [],
        'direction' : [],
        'SSP_start' : [],
        'SSP_end' : [],
        'SSP_edit_distance' : [],
        'biological_seq_indices' : []
    }
    for i in range(len(seqs)) :
        parsed_seq = parse_SSP_and_barcode(str(seqs[i]), barcodes, SSP, max_barcode_score, max_SSP_score, max_gap, min_length_barcode, min_length_SSP)
        for key in parsed_seqs :
            parsed_seqs[key].append(parsed_seq[key])
    for key in parsed_seqs :
        table = table.append_column(key, [parsed_seqs[key]])
    return table

def debarcode_table_from_file(file, barcodes, SSP, max_barcode_score, max_SSP_score, max_gap=5, min_length_barcode=0, min_length_SSP=0, resume=False, overwrite=False) :
    table = pq.read_table(file)
    if 'barcode_ID' in table.column_names :
        if resume == True :
            print('skip resume', file)
            print(table.column_names)
            del table
            return
        elif overwrite == True :
            table = table.drop_columns([
                'barcode_ID',
                'barcode_UMI_start',
                'barcode_UMI_end',
                'barcode_UMI_edit_distance',
                'barcode_SSP_start',
                'barcode_SSP_end',
                'barcode_SSP_edit_distance',
                'combined_score',
                'direction',
                'SSP_start',
                'SSP_end',
                'SSP_edit_distance',
                'biological_seq_indices'
            ])
        else :
            print('skip no overwrite', file)
            del table
            return
    table = debarcode_table(table, barcodes, SSP, max_barcode_score, max_SSP_score, max_gap, min_length_barcode, min_length_SSP)
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

def debarcode(dataset_dir, barcode_path, SSP, max_barcode_score, max_SSP_score, max_gap=5, min_length_barcode=0, min_length_SSP=0, workers=4, resume=False, overwrite=False ) :
    barcodes = load_barcodes(barcode_path)
    files = [x for x in Path(dataset_dir).iterdir() if x.is_file()]
    split = math.ceil( len(files) / workers )
    files_chunks = np.array_split(np.array(files), split)
    for chunk in files_chunks : 
#         debarcode_table_from_file( file, barcodes, SSP, max_barcode_score, max_SSP_score, max_gap, min_length_barcode, min_length_SSP, resume, overwrite )
        with concurrent.futures.ProcessPoolExecutor( max_workers=workers ) as executor :
            futures = [ executor.submit( debarcode_table_from_file, file, barcodes, SSP, max_barcode_score, max_SSP_score, max_gap, min_length_barcode, min_length_SSP, resume, overwrite ) for file in chunk]
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