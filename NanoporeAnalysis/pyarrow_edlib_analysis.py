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
        i = 0
        for line in handle.readlines() :
            if line[0] != '@' :
                line_split = line.split('\t')
                line_dict = dict.fromkeys(table_dict_keys)
                
                line_dict['ID'] = line_split[0]
                line_dict['seq'] = line_split[9]
                line_dict['seq_len'] = len(seq)
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

def find_seq_matches(target_seq, query_seq, max_edit_distance, min_length=0) :
    
    matches = []
    for_alignment = edlib.align(query_seq, target_seq, mode='HW', task='path')
    rev_alignment = edlib.align(utils.reverse_complement(query_seq), target_seq, mode='HW', task='path')
    for location in for_alignment['locations'] :
        length = abs(location[0] - location[1])
        if for_alignment['editDistance'] <= max_edit_distance and length >= min_length :
            matches.append({
                'edit_distance' : for_alignment['editDistance'],
                'edit_score' : for_alignment['editDistance'] / length,
                'location' : location, 
                'direction' : 'forward'
            })
    for location in rev_alignment['locations'] :
        length = abs(location[0] - location[1])
        if rev_alignment['editDistance'] <= max_edit_distance and length >= min_length :
            matches.append({
                'edit_distance' : rev_alignment['editDistance'],
                'edit_score' : for_alignment['editDistance'] / length,
                'location' : location, 
                'direction' : 'reverse'
            })
    return matches
def find_barcode_SSP_pairs(barcode_matches, barcode_ID, SSP_match, max_gap = 5) :
    barcode_SSP_pairs = []
    for barcode_match in barcode_matches :
        
        barcode_start = barcode_match['location'][0]
        barcode_end = barcode_match['location'][1]
        SSP_start = SSP_match['location'][0]
        SSP_end = SSP_match['location'][1]
        
        forward_pair_gap = abs(SSP_start - barcode_end  )
        reverse_pair_gap = abs( barcode_start - SSP_end)
        if forward_pair_gap <= max_gap and barcode_match['direction'] == 'forward' and SSP_match['direction'] == 'forward' :
            barcode_SSP_pairs.append({
                'barcode_ID' : barcode_ID,
                'barcode_start' : barcode_start,
                'barcode_end' : barcode_end,
                'barcode_score' : barcode_match['edit_distance'],
                'SSP_start' : SSP_start,
                'SSP_end' : SSP_end,
                'SSP_score' : SSP_match['edit_distance'],
                'combined_score' : barcode_match['edit_score'] + SSP_match['edit_score'],
                'direction' : 'forward'
            })
        elif reverse_pair_gap <= max_gap and barcode_match['direction'] == 'reverse' and SSP_match['direction'] == 'reverse' :
            barcode_SSP_pairs.append({
                'barcode_ID' : barcode_ID,
                'barcode_start' : barcode_start,
                'barcode_end' : barcode_end,
                'barcode_score' : barcode_match['edit_distance'],
                'SSP_start' : SSP_start,
                'SSP_end' : SSP_end,
                'SSP_score' : SSP_match['edit_distance'],
                'combined_score' : barcode_match['edit_score'] + SSP_match['edit_score'],
                'direction' : 'reverse'
            })
    
    return barcode_SSP_pairs

def parse_SSP_and_barcode(seq, barcodes, SSP, max_barcode_score, max_SSP_score, max_gap = 5, min_length_barcode=0, min_length_SSP=0) :
    
    seq_len = len(seq)
    barcode_SSP_pairs = []
    distal_SSP = {
        'edit_distance' : 'None',
        'location' : ['None', 'None'],
        'direction' : 'None'
    }
    SSP_matches = find_seq_matches(seq, SSP, max_SSP_score, min_length = min_length_SSP)
    if SSP_matches == [] :
        SSP_matches.append({
            'edit_distance' : 100,
            'location' : [0, seq_len], 
            'direction' : 'reverse'
        })
    for barcode in barcodes :
        barcode_matches = find_seq_matches(seq, barcode[1], max_barcode_score, min_length = min_length_barcode)
        for SSP_match in SSP_matches :
            for pair in find_barcode_SSP_pairs(barcode_matches, barcode[0], SSP_match, max_gap) :
                barcode_SSP_pairs.append(pair)
                
    barcode_SSP_pair = {
        'combined_score' : 100,
        'barcode_ID' : 'None'
    }
    for pair in barcode_SSP_pairs :
        if pair['combined_score'] < barcode_SSP_pair['combined_score'] :
            barcode_SSP_pair = pair
        elif pair['combined_score'] == barcode_SSP_pair['combined_score'] :
            barcode_SSP_pair['barcode_ID'] = 'Multiple'
            
    if barcode_SSP_pair['barcode_ID'] not in ['Multiple', 'None']  :
        if barcode_SSP_pair['direction'] == 'reverse' :
            seq = utils.reverse_complement(seq)
            barcode_start = seq_len - barcode_SSP_pair['barcode_end']
            SSP_end = seq_len - barcode_SSP_pair['SSP_start']
            
            for possible_distal_SSP_match in SSP_matches :
                if possible_distal_SSP_match['direction'] == 'forward' and possible_distal_SSP_match['location'][0] > barcode_SSP_pair['barcode_end'] :
                    location = [ seq_len - possible_distal_SSP_match['location'][1], seq_len - possible_distal_SSP_match['location'][0] ]
                    possible_distal_SSP_match = {
                        'edit_distance' : possible_distal_SSP_match['edit_distance'],
                        'location' : location,
                        'direction' : 'forward'
                    }
                    if distal_SSP['edit_distance'] == 'None' :
                        distal_SSP = possible_distal_SSP_match
                    elif distal_SSP['edit_distance'] > possible_distal_SSP_match['edit_distance'] :
                        distal_SSP = possible_distal_SSP_match
                    elif distal_SSP['location'][0] > possible_distal_SSP_match['location'][0] :
                        distal_SSP = possible_distal_SSP_match
        else :
            barcode_start = barcode_SSP_pair['barcode_start']
            SSP_end = barcode_SSP_pair['SSP_end']
            
            for possible_distal_SSP_match in SSP_matches :
                if possible_distal_SSP_match['direction'] == 'reverse' and possible_distal_SSP_match['location'][1] < barcode_SSP_pair['barcode_start'] :
                    if distal_SSP['edit_distance'] == 'None' :
                        distal_SSP = possible_distal_SSP_match
                    elif distal_SSP['edit_distance'] > possible_distal_SSP_match['edit_distance'] :
                        distal_SSP = possible_distal_SSP_match
                    elif distal_SSP['location'][0] > possible_distal_SSP_match['location'][0] :
                        distal_SSP = possible_distal_SSP_match
                        
        parsed = {
            'barcode_ID' : str(barcode_SSP_pair['barcode_ID']),
            'barcode_SSP_pair_start' : str(barcode_start),
            'barcode_SSP_pair_end' :str( SSP_end),
            'barcode_SSP_pair_score' : str(barcode_SSP_pair['combined_score']),
            'SSP_start' : str(distal_SSP['location'][0]),
            'SSP_end' : str(distal_SSP['location'][1]),
            'SSP_score' : str(distal_SSP['edit_distance']),
            'direction' : str(barcode_SSP_pair['direction'])
        }
    elif barcode_SSP_pair['barcode_ID'] == 'None' :
        parsed = {
            'barcode_ID' : 'None matched',
            'barcode_SSP_pair_start' : 'None',
            'barcode_SSP_pair_end' : 'None',
            'barcode_SSP_pair_score' : 'None',
            'SSP_start' : 'None',
            'SSP_end' : 'None',
            'SSP_score' : 'None',
            'direction' : 'None'
        }
    elif barcode_SSP_pair['barcode_ID'] == 'Multiple' :
        parsed = {
            'barcode_ID' : 'Multiple',
            'barcode_SSP_pair_start' : 'None',
            'barcode_SSP_pair_end' : 'None',
            'barcode_SSP_pair_score' : 'None',
            'SSP_start' : 'None',
            'SSP_end' : 'None',
            'SSP_score' : 'None',
            'direction' : 'None'
        }
        
    return parsed

def debarcode_table(table, barcodes, SSP, max_barcode_score, max_SSP_score, max_gap=5, min_length_barcode=0, min_length_SSP=0) :
    seqs = table.column('seq')
    duplex_tags = table.column('dx:i')
    parsed_seqs = {
        'barcode_ID' : [],
        'barcode_SSP_pair_start' : [],
        'barcode_SSP_pair_end' : [],
        'barcode_SSP_pair_score' : [],
        'SSP_start' : [],
        'SSP_end' : [],
        'SSP_score' : [],
        'direction' : []
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
        table = pq.read_table(file)
        table = debarcode_table(table, barcodes, SSP, max_barcode_score, max_SSP_score, max_gap, min_length_barcode, min_length_SSP)
    return table