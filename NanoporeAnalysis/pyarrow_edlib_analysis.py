def pod5_split() :
    
    #generate the index of all reads and their channels
    if not os.path.isfile(Path(path_out + "/view.txt")) or reset_pod5_view == True :
        pod5_view.view_pod5([Path(path_data)], Path(path_out), include = "read_id, channel", force_overwrite=True, threads=threads)
    else :
        view = pd.read_table(Path(path_out + "/view.txt"), sep='\t')
    
    channels = view['channel'].unique()
    channels_per_worker = np.array_split(np.array(channels), workers)
    
    for i in range(workers) :
        
        view.query("channel in " str(channels_per_worker[i]))['read_id'].to_csv(Path(path_out + "/view_current.txt"), index = False, sep=' ', header=False)
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
            score = for_alignment['editDistance'] + query_seq_len - alignment_length / query_seq_len
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
            score = rev_alignment['editDistance'] + query_seq_len - alignment_length / query_seq_len
            matches.append({
                'query_ID' : query_ID,
                'edit_distance' : rev_alignment['editDistance'],
                'edit_score' : score,
                'location' : [ target_seq_len - location[0], target_seq_len - location[1] ],
                'direction' : 'reverse'
            })
    return matches

def pick_best_match(matches) :
    best_match = {
        'edit_score' : 1000,
    }
    for match in matches :
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

def parse_SSP_and_barcode(seq, barcodes, SSP, max_barcode_score, max_SSP_score, max_gap = 5, min_length_barcode = 0, min_length_SSP = 0) :
    seq_len = len(seq)
    SSP_matches = find_seq_matches(seq, SSP, max_SSP_score, 'SSP', min_length = min_length_SSP)
    barcode_matches = []
    for barcode in barcodes :
        tmp_barcode_matches = find_seq_matches(seq, barcode[1], max_barcode_score, barcode[0], min_length = min_length_barcode)
        if len(tmp_barcode_matches) > 0 :
            barcode_matches.append(pick_best_match(tmp_barcode_matches))
    barcode_SSP_pair = {
        'barcode_ID' : 'None matched',
        'barcode_UMI_start' : 0,
        'barcode_UMI_end' : 0,
        'barcode_UMI_edit_distance' : 0,
        'barcode_SSP_start' : 0,
        'barcode_SSP_end' : 0,
        'barcode_SSP_edit_distance' : 0,
        'combined_score' : 1000,
        'direction' : 'None'
    }
    for barcode_match in barcode_matches :
        for SSP_match in SSP_matches :
            barcode_SSP_evaluation = evaluate_barcode_SSP_pair(barcode_match, SSP_match, max_gap)
            if barcode_SSP_evaluation != None :
                if barcode_SSP_evaluation['combined_score'] < barcode_SSP_pair['combined_score'] :
                    barcode_SSP_pair = barcode_SSP_evaluation
                elif barcode_SSP_evaluation['combined_score'] == barcode_SSP_pair['combined_score'] :
                    barcode_SSP_pair['barcode_ID'] = 'Multiple'
                    barcode_SSP_pair['barcode_ID'] = barcode_SSP_evaluation['combined_score']
    
    if barcode_SSP_pair['direction'] == 'reverse' :
        seq = utils.reverse_complement(seq)
    elif barcode_SSP_pair['barcode_ID'] == 'None matched' and len(barcode_matches) > 0 :
        best_barcode = pick_best_match(barcode_matches)
        barcode_SSP_pair['barcode_ID'] = best_barcode['query_ID']
        barcode_SSP_pair['barcode_UMI_start'] = best_barcode['location'][0]
        barcode_SSP_pair['barcode_UMI_end'] = best_barcode['location'][1]
        barcode_SSP_pair['barcode_UMI_edit_distance'] = best_barcode['edit_distance']
        barcode_SSP_pair['direction'] = best_barcode['direction']
    
    distal_SSP = {
        'SSP_start' : 0,
        'SSP_end' : 0,
        'SSP_edit_distance' : 1000,
    }
    for SSP_match in SSP_matches :
        if SSP_match['direction'] != barcode_SSP_pair['direction'] :
            if SSP_match['edit_score'] < distal_SSP['SSP_edit_distance'] :
                distal_SSP['SSP_start'] = SSP_match['location'][0]
                distal_SSP['SSP_end'] = SSP_match['location'][1]
                distal_SSP['SSP_edit_distance'] = SSP_match['edit_distance']
            elif SSP_match['edit_score'] == distal_SSP['SSP_edit_distance'] and SSP_match['location'][1] < distal_SSP['SSP_end'] :
                distal_SSP['SSP_start'] = SSP_match['location'][0]
                distal_SSP['SSP_end'] = SSP_match['location'][1]
                distal_SSP['SSP_edit_distance'] = SSP_match['edit_distance']
    
    
    parsed = barcode_SSP_pair
    parsed.update(distal_SSP)
    if parsed['barcode_UMI_start'] - parsed['SSP_end'] > seq_len * 0.5 :
        parsed['biological_seq_indices'] = [ parsed['SSP_end'], parsed['barcode_UMI_start'] ]
    else :
        parsed['biological_seq_indices'] = [ 0, seq_len - 1 ]
    return parsed

def debarcode_table(table, barcodes, SSP, max_barcode_score, max_SSP_score, max_gap=5, min_length_barcode=0, min_length_SSP=0) :
    seqs = table.column('seq')
    duplex_tags = table.column('dx:i')
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

def load_barcodes(path) :
    barcodes = []
    with open(path, 'r') as handle :
        for line in handle.readlines() :
            line_split = line.split(',')
            barcodes.append([line_split[0], utils.reverse_complement(line_split[2])])
    return barcodes

def debarcode(dataset_dir, barcode_path, SSP, max_barcode_score, max_SSP_score, max_gap=5, min_length_barcode=0, min_length_SSP=0) :
    barcodes = load_barcodes(barcode_path)
    files = [x for x in Path(dataset_dir).iterdir() if x.is_file()]
    for file in files :
        table = pq.read_table(file).slice(0,1000)
        table = debarcode_table(table, barcodes, SSP, max_barcode_score, max_SSP_score, max_gap, min_length_barcode, min_length_SSP)
    return table