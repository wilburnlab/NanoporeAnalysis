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
from pyarrow import csv
import matplotlib.pyplot as plt
import mappy
import numpy as np
from NanoporeAnalysis import utils
from NanoporeAnalysis import local_io
import edlib

def pod5_split(path_out, path_pod5, python_env, account, mail, mail_type = 'None', workers=5) :
    """
    Runs slurm jobs to split .pod5 files into separate folders based on sequencing channel. This is needed for running dorado duplexing on the OSC with multiple workers. .pod5 files are split and saved under path_out/split_pod5/{channel}, and the job scripts and .out logs are saved under path_out/logs_split_pod5. Currently requires a specified python environment with the pod5 package from ONT installed.
    
    Args:
        path_out (str) : directory to be used as the output folder. This will be created if it doesn't exist.
        path_pod5 (str) : directory containing .pod5 files. The function will recursively find the files in this folder.
        python_env (str) : name of your conda environment with pod5 installed.
        account (str) : OSC account to be billed for the compute time.
        mail (str) : email for any status messages from the jobs.
        mail_type (str) : type of messages to recieve for the jobs ie 'Start', 'All', 'Error'. Defaults to 'None'
        workers (int) : number of jobs to run for the pod5 splitting. Defaults to 5. Rough recommendation is to do at most 5M reads per worker ie use >= 10 workers for 50M reads. Doing more just runs the risk of memory and time issues, and the OSC is hardly ever out of CPU space.
    """
    Path(path_out).mkdir(parents=True, exist_ok=True)
    if not os.path.isfile(Path(path_out + "/view.txt")) :
        pod5_view.view_pod5([Path(path_pod5)], Path(path_out + '/view.txt'), include = "read_id, channel", force_overwrite=True, threads=workers)
    view = csv.read_csv(path_out + "/view.txt", parse_options=csv.ParseOptions(delimiter='\t'))
    channels = view.column('channel').unique().to_pylist()
    if Path(path_out + '/split_pod5/').is_dir() :
        shutil.rmtree(Path(path_out + '/split_pod5/'))
    if Path(path_out + '/logs_split_pod5/').is_dir() :
        shutil.rmtree(Path(path_out + '/logs_split_pod5/'))
    Path(path_out + '/split_pod5/').mkdir(parents=True, exist_ok=True)
    Path(path_out + '/logs_split_pod5/').mkdir(parents=True, exist_ok=True)
    for channel in channels :
        Path(path_out + '/split_pod5/' + str(channel)).mkdir(parents=True, exist_ok=True)
    files = [str(x) for x in Path(path_pod5).iterdir()]
    file_chunks = np.array_split(files, workers )
    print("done with prep work")
    for i in range(workers) :
        script = [
            "#!/bin/bash\n",
            "#SBATCH --account=" + account + "\n",
            "#SBATCH --job-name=pod5_split_" + str(i) + "\n",
            "#SBATCH --nodes=1 --ntasks-per-node=1\n",
            "#SBATCH --cpus-per-task=4\n",
            "#SBATCH --output=" + path_out + "/logs_split_pod5/" + str(i) + ".out" + "\n",
            "#SBATCH --mail-type=" + mail_type + "\n",
            "#SBATCH --time=04:00:00\n",
            "#SBATCH --mail-user=" + mail + "\n",
            "module load miniconda3/23.3.1-py310\n",
            "conda activate " + python_env + "\n",
            "pod5 subset " + ' '.join(file_chunks[i]) + " --table " + path_out + "/view.txt -o " + path_out + "/split_pod5/ --threads 4 --template '{channel}/{channel}_" + str(i) + ".pod5' -M --columns channel"
        ]
        with open(str(path_out + "/logs_split_pod5/" + str(i) + ".sh"), 'w') as handle:
            handle.writelines(script)
        args = ["sbatch", str(path_out + "/logs_split_pod5/" + str(i) + ".sh")]
        subprocess.run(args)
    return

def copy_split_pod5s(path_out, workers_dorado = 1, workers_copying = 10) :
    """
    Copies the split pod5 files from pod5_split into folders for each dorado worker in run_dorado_job. Finds pod5s in path_out/split_pod5/ and creates a folder for each dorado worker under path_out/grouped_pod5/. For the sake of speed, this can run in parallel. 
    
    Args:
        path_out (str) : directory to be used as the output folder, same as what was used for pod5_split.
        workers_dorado (int) : number of folder to split the pod5s into. Should match the number of dorado jobs that you want to run.
        workers_copying (int) : number of parallel workers to use for the copying process. Defaults to 10.
    """
    if Path(path_out + '/grouped_pod5/').is_dir() :
        shutil.rmtree(Path(path_out + '/grouped_pod5/'))
    Path(path_out + "/grouped_pod5/").mkdir(parents=True, exist_ok=True)
    channel_folders = [x for x in Path(path_out + "/split_pod5/").iterdir()]
    channel_folder_chunks = np.array_split( channel_folders, workers_dorado )
    for i in range(len(channel_folder_chunks)) :
        print("started moving group ", i)
        Path( path_out + "/grouped_pod5/" + str(i) ).mkdir(parents=True, exist_ok=True)
        with concurrent.futures.ProcessPoolExecutor( max_workers = workers_copying ) as executor :
            futures = [ executor.submit( shutil.copytree, folder, Path( path_out + "/grouped_pod5/" + str(i) + "/" + folder.stem ) ) for folder in channel_folder_chunks[i] ]
            concurrent.futures.wait( futures )
        print("finished moving group ", i)
    return

def check_pod5_jobs(user) :
    """
    Checks squeue to see if any pod5_split jobs are running.
    
    Args :
        user (str) : the OSC username used for the pod5_split jobs.
    Returns :
        bool, True if any jobs are found.
    """
    args = ['squeue', '-u', user, '-v', '-r', '--states=BF,CF,CG,DL,F,NF,OOM,PD,PR,R,RD,RS,RV,ST,S,TO', '--format="%.18i %.9P %.20j %.8u %.2t %.10M %.6D %R"']
    squeue = subprocess.run(args, capture_output=True)
    stdout = str(squeue.stdout)
    squeue_nl_split = stdout.split('\\n')
    for item in squeue_nl_split :
        if 'pod5_split' in item :
            return True
    return False

def run_dorado_job(path_out, path_dorado, account, mail, mail_type = 'None', i = 0) :
    """
    Runs a slurm job for dorado duplexing, sending the output to path_out/basecalled/{i}.sam. Creates a script file in path_out/logs_dorado/, which is also where the .out file from the job will be. Currently is hardcoded to request 1 GPU, 12 CPU cores, 80G memory, and 12 hours.
    
    Args :
        path_out (str) : directory to be used as the output folder. This will fail if it doesn't exist already.
        path_dorado (str) : path to the dorado executable.
        account (str) : OSC account to be billed for the compute time.
        mail (str) : email for any status messages from the jobs.
        mail_type (str) : type of messages to recieve for the jobs ie 'Start', 'All', 'Error'. Defaults to 'None'
        i (int) : the worker_ID to be appended to the job script and .sam output file.
    """
    script = [
        "#!/bin/bash\n",
        "#SBATCH --account=" + account + "\n",
        "#SBATCH --job-name=dorado_" + str(i) + "\n",
        "#SBATCH --nodes=1\n",
        "#SBATCH --ntasks=1\n",
        "#SBATCH --cpus-per-task=12\n",
        "#SBATCH --mem=80G\n",
        "#SBATCH --gres=gpu:1\n",
        "#SBATCH --time=12:00:00\n",
        "#SBATCH --output=" + path_out + "/logs_dorado/dorado_" + str(i) + ".out" + "\n",
        "#SBATCH --mail-type=" + mail_type + "\n",
        "#SBATCH --mail-user=" + mail + "\n",
        str(path_dorado) + " duplex sup --emit-sam -r -v " + str(path_out) + "/grouped_pod5/" + str(i) + "/ > " + str(path_out) + "/basecalled/" + str(i) + ".sam"
    ]
    with open(path_out + '/logs_dorado/dorado_' + str(i) + '.sh', 'w') as handle:
        handle.writelines(script)
    args = ['sbatch', str(path_out + '/logs_dorado/dorado_' + str(i) + '.sh')]
    subprocess.run(args)
    return

def dorado_slurm(path_out, path_dorado, account, mail, mail_type = 'None', python_env = None, workers_dorado = 1, workers_pod5_copy = 10, skip_split = False, skip_copy = False) :
    """
    Run dorado duplex basecalling through slurm jobs. First will split the .pod5 files by sequencing channel using pod5_split(), unless skip_split is set to True. Puts split files under path_out/split_pod5/. Uses check_pod5_jobs() to make sure no pod5_split jobs are running under the current OSC account and will not run if there are, otherwise things will break. Uses copy_split_pod5(0) to copy the split pod5s into a folder for each dorado worker under path_out/grouped_pod5/. Then runs parallel dorado jobs requesting GPUs, several CPUs, and sufficient RAM. Note that this will only run on NVIDIA GPUs with Tensor cores, which are only in the Voltair line and newer. The Owens cluster does not have these, and dorado will not run on it.
    
    Args :
        path_out (str) : directory to be used as the output folder, same as what was used in pod5_split. This will be created if it doesn't exist.
        path_dorado (str) : path to the dorado executable.
        account (str) : OSC account to be billed for the compute time.
        mail (str) : email for any status messages from the jobs.
        mail_type (str) : type of messages to recieve for the jobs ie 'Start', 'All', 'Error'. Defaults to 'None'
        python_env (str) : name of your conda environment with pod5 installed. If your files are already split, this is optional.
        workers_dorado (int) : number of jobs to run for dorado. If you already ran pod5_split, this should be the SAME number of workers as used for that. You'll either leave some reads un-basecalled or run unnecessary jobs if not.
        workers_pod5_copy (int) : number of workers to use when copying the pod5 files from pod5_split. This just helps speed things up a bit if you're running on multiple cores, it's not super critical.
        skip_split (bool) : whether to skip over running pod5_split(). Set this to true if you've already run pod5_split or if you're just re-running this function. If False, this will overwrite any existing files.
        skip_copy (bool) : whether to skip over running copy_split_pod5s(). Set this to true if you're just re-running this function. If False, this will overwrite any existing files.
    """
    if not skip_split :
        pod5_split(path_out, path_pod5, python_env, account, mail, mail_type = 'None', workers = workers_dorado)
    while check_pod5_jobs() :
        print("pod5 spit still running")
        time.sleep(180)
    if not skip_copy :
        copy_split_pod5s(path_out, workers_dorado, workers_pod5_copy)
    if Path(path_out + '/basecalled/').is_dir() :
        shutil.rmtree(Path(path_out + '/basecalled/'))
    Path(path_out + "/basecalled/").mkdir(parents=True, exist_ok=True)
    if Path(path_out + '/logs_dorado/').is_dir() :
        shutil.rmtree(Path(path_out + '/logs_dorado/'))
    Path(path_out + "/logs_dorado/").mkdir(parents=True, exist_ok=True)
    for i in range(workers_dorado) :
        run_dorado_job(path_out, path_dorado, account, mail, i)
    return

def check_dorado_jobs(path_out, rerun_broken_jobs = True, path_dorado = None, account = None, mail = None, mail_type = None) :
    """
    Checks the .out logs in path_out/logs_dorado/ to make sure that all dorado jobs ran successfully. Checks for any lines in the .out with 'Basecalled' in them, which appears to only be present if dorado has completed basecalling. Unless rerun_broken_jobs is set to False, this will also rerun the broken jobs via run_dorado_job(). path_out, path_dorado, account, mail, and mail_type must all be set to appropriate things in order to the rerun to work. They are only given defaults so that this can be run easily with rerun_broken_jobs = False.
    
    Args :
        path_out (str) : directory used as the path_out folder for dorado.
        rerun_broken_jobs (bool) : whether or not to rerun any broken jobs that are found. Defaults to True.
        path_dorado (str) : path to the dorado executable. 
        account (str) : OSC account to be billed for the compute time.
        mail (str) : email for any status messages from the jobs.
        mail_type (str) : type of messages to recieve for the jobs ie 'Start', 'All', 'Error'. Defaults to 'None'
    """
    dorado_logs = [ x for x in Path(path_out + "/logs_dorado/").iterdir() if x.suffix == '.out' ]
    broken_jobs = []
    for dorado_log in dorado_logs :
        broken = True
        with open(dorado_log, 'r') as handle :
            for line in handle.readlines() :
                if 'Basecalled' in line :
                    broken = False
                    break
        if broken == True :
            name = dorado_log.stem
            i = name.split('_')[1]
            broken_jobs.append(i)
    if rerun_broken_jobs :
        for broken_job in broken_jobs :
            run_dorado_job(path_out, path_dorado, account, mail, broken_job)
    else :
        print("Broken jobs were found : ", broken_jobs)
    return

def sam_to_parquet(file, path_out, basename_template = None) :
    """
    Transfers data from .sam files to an Apache parquet database by going line-by-line through the .sam file and building a pyarrow table. Generates ID, seq, seq_len, and qual tags at minimum, then automatically creates tags for any other content in the .sam file assuming the usual format where of abcd:xyz... where xxxx is the name and anything x and beyond is the value. Also, this skips over all tags in the first ten columns of the .sam, other than the three explicitly pulled tags: ID, seq, and qual. The others do not seem to be used by dorado. Outputs to path_out/pa_dataset/ following the basename_template. Current hardcoded to limit the number of entries in the files to 200000.
    
    Args :
        file (str or Path) : the .sam file to be converted.
        path_out (str) : the directory being used for output. The /pa_dataset/ folder will be under this directory. Must already exist.
        basename_template (str) : the template to be used for writing .parquet files. Requires {i} to be present, which will be replaced by a number for each different file created ie 'dataset_A_part_{i}.parquet' where {i} will be 0 for the first file, 1, 2, etc. for subsequent files if the dataset is too large for one file.
    """
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

def build_parquet_dataset_from_sam(sam_folder, path_out) :
    """
    Takes all .sam files in sam_folder and builds an Apache parquet database under path_out/pa_dataset. This function should be parallelized with multiple instances of sam_to_parquet.
    
    Args :
        sam_folder (str) : directory containing .sam files. Does not currently search recursively.
        path_out (str) : the directory being used for output. The /pa_dataset/ folder will be under this directory. Must already exist.
    """
    files = [x for x in Path(sam_folder).iterdir() if x.is_file()]
    for file in files :
        print("moving file ", file.stem)
        sam_to_parquet(file, path_out, basename_template = str(file.stem + '_part-{i}.parquet'))
    return 'done'

def find_seq_matches(query_seq, target_seq, max_edit_distance, query_ID = 'None', min_length = 0, skip_reverse = False) :
    """
    Finds matches between target_seq and query_seq using edlib.align in both forward and reverse direction, then translates the result into a list of location matches.
    
    Args :
        query_seq (str) : sequence to be used as query. This is usually the smaller sequence and within the target.
        target_seq (str) : sequence to be used as target. This is usually the larger sequence, and the query is meant to be within it.
        max_edit_distance (int) : the maximum allowable editDistance from an edlib alignment. Anything higher will be discarded. The edit distance is a count of how many characters have to be changed to turn the query into the target.
        query_ID (str) : the name of the query sequence. Defaults to 'None'.
        min_length (int) : the minmum allowed length of an alignment. Anything shorter will be discarded.
        skip_reverse (bool) : whether or not to skip the reverse alignment. Only useful for when aligning two things that have known orientations, ie the distal SSP after assigning the barcode. Defaults to False.
        
    Returns :
        matches (list) : a list of dictionaries for each alignment that meets the requirments.
    """
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

def merge_overlapped_indices(index_pairs: list, tolerated_mismatches: int) -> list:
    """
    Provided a list of (start,end) index pairs, combine overlapping pairs into single pairs that span the full alignment range.
    
    Args :
        index_pairs (list) : list of [start (int),end (int)] tuples.
        tolerated_mismatches (int) : number of allowed gap characters between alignments that will still be grouped together.
        
    Returns :
        reduced_pairs (list) : list of condensed [start (int),end (int)] tuples.
    """
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
    
def find_polyX(sequence: str, X: str, N: int, tolerated_mismatches: int, min_len: int) -> list:
    """
    Find runs of single nucleotides (X), using edlib align with search query X*N with tolerated_mismatches in the search
    
    Args :
        sequence (str) : sequence to be searched.
        X (str) : single character to be used. Must be just 1 character long, not case sensitive.
        N (int) : number of characters to be considered a polyX sequence.
        tolerated_mismatches (int) : number of allowed gap characters between alignments that will still be grouped together.
        min_len (int) : minimum length of the final polyX region(s).
        
    Returns :
        sorted_indices (list) : list of [start (int),end (int)] tuples sorted by length.
    """
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
    
def find_polyA_bidirectional(sequence: str, N: int, tolerated_mismatches: int, min_len: int) -> list:
    """
    Finds any polyA sequences in the sequence in either forward or reverse direction that meets the minimum length.
    
    Args :
        sequence (str) : sequence to be searched.
        N (int) : number of characters to be considered a polyX sequence.
        tolerated_mismatches (int) : number of allowed gap characters between alignments that will still be grouped together.
        min_len (int) : minimum length of the final polyX region(s).
        
    Returns :
        indices (list) : list of [start (int),end (int)] tuples.
    """
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
    """
    Selects the best alignment match from a list. Not currently used... also can't remember how the min_match_location was used.
    """
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
    """
    Judge a set of polyA, UMI, and SSP alignments based on whether they're in the correct order (polyA, then UMI, then SSP) and are close enough to each other via max_gap. Uses polyA_UMI_SSP_set as a blank template and fills it in. First calculates gaps, then checks if directions of UMI and polyA are matching and their gap is allowable, then fills out the template polyA_UMI_SSP_set from them. 
    
    Args :
        polyA (list): the indices of a polyA sequence from find_polyA_bidirectional(). This is a tuple of the start and end indices [start, end].
        UMI_match (dict) : the alignment for a UMI, output by find_seq_matches(). This is a dictionary of the alignment attributes.
        SSP_match (dict) : the alignment for the SSP sequence, output by find_seq_matches(). This is a dictionary of the alignment attributes.
        polyA_UMI_SSP_set (dict) : a template dict with default values for the keys. The default values are intended to indicate an invalid match, as failing to meet the requirements of this function will result in the default values being returned.
        max_gap (int) : the maximum gap between the UMI and the polyA and SSP. Defaults to 5.
    
    Returns :
        polyA_UMIT_SSP_set (dict) : the template dict with changed values, corresponding to whether the three alignments form a valid set representing a true barcode.
    """
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
    """
    Takes match outputs from find_seq_matches and reverses the orientation by subtracting the location indices from the length of the sequence. All matches must at least correspond to the same sequence length, usually the same sequence as well.
    
    Args :
        matches (list) : list of the matches to be reversed.
        seq_len (int) : length of the sequence that the matches correspond to.
    
    Returns :
        reversed_matches (list) : list of matches after being reversed. They are in the same order as given in matches.
    """
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
    """
    Finds potential unique molecular identifiers (UMIs), strand switching primer (SSP), and polyA sequences that meet certain minimum quality criteria, then determines if a UMI can be paired with an SSP and/or polyA alignment to represent a true barcode sequence. 
    
    Args :
        seq (str) : the sequence to be checked.
        UMIs (list) : list of [name, sequence] pairs for UMIs.
        SSP (str) : sequence for the SSP.
        UMI_max_score (int) : maximum edit_distance for the UMI alignments. Any alignments with a higher edit_distance are discarded.
        SSP_max_score (int) : maximum edit_distance for the SSP alignments. Any alignments with a higher edit_distance are discarded.
        max_gap (int) : the maximum gap between the UMI and the polyA or SSP. Defaults to 5
        UMI_min_len (int) : minimum length for the UMI alignments.
        SSP_min_len (int) : minimum length for the SSP alignments.
        polyA_N (int) : number of characters to be considered a polyA sequence. Defaults to 4.
        polyA_tolerated_mismatches (int) : number of allowed gap characters between polyA alignments that will be grouped together into a larger polyA. Defaults to 1.
        polyA_min_len (int) : minimum length of the final polyX region(s) that will be considered a true polyA. Defaults to 10.
    
    Returns :
        parsed (dict) : dictionary with the final decided attributes for the seq and whether it contains an identifiable UMI/barcode.
    """
    seq_len = len(seq)
    polyAs = find_polyA_bidirectional(seq, N = polyA_N, tolerated_mismatches = polyA_tolerated_mismatches, min_len = polyA_min_len)
    SSP_matches = find_seq_matches(SSP, seq, SSP_max_score, 'SSP', min_length = SSP_min_len)
    if len(SSP_matches) == 0 :
        SSP_matches = [{ 'location' : [-1,-1], 'direction' : 'None', 'edit_distance' : 1000 }]
    UMI_matches = []
    for UMI in UMIs :
        tmp_UMI_matches = find_seq_matches(UMI[1], seq, UMI_max_score, UMI[0], min_length = UMI_min_len)
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
    """
    Takes a pyarrow table with one sequences in each row under a column titled 'seq' and runs parse_polyA_UMI_SSP() on all the reads within.
    
    Args :
        table (pyarrow table) : table with a column called 'seq' with the sequences to be analyzed. Seqs must be strings.
        UMIs (list) : list of [name, sequence] pairs for UMIs.
        SSP (str) : sequence for the SSP.
        UMI_max_score (int) : maximum edit_distance for the UMI alignments. Any alignments with a higher edit_distance are discarded.
        SSP_max_score (int) : maximum edit_distance for the SSP alignments. Any alignments with a higher edit_distance are discarded.
        max_gap (int) : the maximum gap between the UMI and the polyA or SSP. Defaults to 5
        UMI_min_len (int) : minimum length for the UMI alignments.
        SSP_min_len (int) : minimum length for the SSP alignments.
        polyA_N (int) : number of characters to be considered a polyA sequence. Defaults to 4.
        polyA_tolerated_mismatches (int) : number of allowed gap characters between polyA alignments that will be grouped together into a larger polyA. Defaults to 1.
        polyA_min_len (int) : minimum length of the final polyX region(s) that will be considered a true polyA. Defaults to 10.
    
    Returns :
        table (pyarrow table) : the input table with the output from parse_polyA_UMI_SSP() added as new columns.
    """
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
    """
    Load a pyarrow table from a parquet file, then runs passes it to debarcode_table() which runs parse_polyA_UMI_SSP() which attempts to find a barcode in all the sequences within the file. This then saves the table back to file, overwriting the original.
    
    Args :
        UMIs (list) : list of [name, sequence] pairs for UMIs.
        SSP (str) : sequence for the SSP.
        UMI_max_score (int) : maximum edit_distance for the UMI alignments. Any alignments with a higher edit_distance are discarded.
        SSP_max_score (int) : maximum edit_distance for the SSP alignments. Any alignments with a higher edit_distance are discarded.
        max_gap (int) : the maximum gap between the UMI and the polyA or SSP. Defaults to 5
        UMI_min_len (int) : minimum length for the UMI alignments.
        SSP_min_len (int) : minimum length for the SSP alignments.
        polyA_N (int) : number of characters to be considered a polyA sequence. Defaults to 4.
        polyA_tolerated_mismatches (int) : number of allowed gap characters between polyA alignments that will be grouped together into a larger polyA. Defaults to 1.
        polyA_min_len (int) : minimum length of the final polyX region(s) that will be considered a true polyA. Defaults to 10.
        resume (bool) : whether to resume debarcoding from a paused or broken run. This will only debarcode files that don't already have barcode information in them, so it can't be used to continue a re-debarcoding session, as the files will all still have the original barcode data. Defaults to False
        overwrite (bool) : whether to allow for overwriting debarcoded files. If True, it will remove any existing barode data and continue with normal debarcoding. If False, it will skip over any files with barcode data.
    """
    table = pq.read_table(file)
    if 'barcode_ID' in table.column_names :
        if resume == True :
            print("skip resume", file)
            print(table.column_names)
            del table
            return
        elif overwrite == True :
            print("removing old barcoding data...")
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
            print("skip no overwrite", file)
            del table
            return
    print("debarcoding file ", file.stem)
    table = debarcode_table(table, UMIs, SSP, UMI_max_score, SSP_max_score, max_gap, UMI_min_len, SSP_min_len, polyA_N, polyA_tolerated_mismatches, polyA_min_len)
    pq.write_table(table, file)
    del table
    return

def load_barcodes(path) :
    """
    Loads barcodes from a csv where each row is a barcode in the shape of Name,SSP/constant sequence,unique molecular identifier (UMI),polyT-VN sequence. This is not a super critical function, it will be better to just load the barcodes into debarcode() as a list, as the only important things are the UMIs and names.
    """
    barcodes = []
    with open(path, 'r') as handle :
        for line in handle.readlines() :
            line_split = line.split(',')
            barcodes.append([line_split[0], utils.reverse_complement(line_split[2])])
    return barcodes

def debarcode(dataset_dir, UMIs_path, SSP, UMI_max_score, SSP_max_score, max_gap = 5, UMI_min_len = 0, SSP_min_len = 0, polyA_N = 4, polyA_tolerated_mismatches = 1, polyA_min_len = 10, workers = 4, resume=False, overwrite=False ) :
    """
    Finds all parquet files within the given dataset directory under dataset_dir, then pushes them through debarcode_table_from_file() to attempt to find barcodes.
    
    Args :
        dataset_dir (str) : the path to the folder containing the parquet files to be debarcoded. Should be path_out/pa_dataset/ where path_out is the same as what was used in build_parquet_dataset_from_sam.
        UMIs_path (str) : path pointing to a csv to feed to load_barcodes(). This should be changed to just accept a list, instead.
        SSP (str) : sequence for the SSP.
        UMI_max_score (int) : maximum edit_distance for the UMI alignments. Any alignments with a higher edit_distance are discarded.
        SSP_max_score (int) : maximum edit_distance for the SSP alignments. Any alignments with a higher edit_distance are discarded.
        max_gap (int) : the maximum gap between the UMI and the polyA or SSP. Defaults to 5
        UMI_min_len (int) : minimum length for the UMI alignments.
        SSP_min_len (int) : minimum length for the SSP alignments.
        polyA_N (int) : number of characters to be considered a polyA sequence. Defaults to 4.
        polyA_tolerated_mismatches (int) : number of allowed gap characters between polyA alignments that will be grouped together into a larger polyA. Defaults to 1.
        polyA_min_len (int) : minimum length of the final polyX region(s) that will be considered a true polyA. Defaults to 10.
        workers (int) : number of parallel debarcoding processes to use. Defaults to 4.
        resume (bool) : whether to resume debarcoding from a paused or broken run. This will only debarcode files that don't already have barcode information in them, so it can't be used to continue a re-debarcoding session, as the files will all still have the original barcode data. Defaults to False
        overwrite (bool) : whether to allow for overwriting debarcoded files. If True, it will remove any existing barode data and continue with normal debarcoding. If False, it will skip over any files with barcode data.
    """
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
    """
    Runs minimap2 on the biological sequences in a pyarrow table. 
    
    Args :
        table (pyarrow table) : contains sequencing reads as rows. This should only be run after running debarcode(), as it relies on some barcode information.
        path_ref (str) : path to the sequence reference to be used.
        preset (str) : the preset to use for minimap2. Defaults to 'splice', which is designed to be splice-aware.
    
    Returns :
        table (pyarrow table) : the input table with minimap alignments appended as new columns.
    """
    seqs = table.column('seq').to_pylist()
    biological_seq_indices = table.column('biological_seq_indices').to_pylist()
    barcode_scores = table.column('barcode_score').to_pylist()
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
        if barcode_scores[i] <= 20 and barcode_IDs[i] != 'Multiple' and SSP_edit_distances[i] <= 5 :
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
    """
    Opens a parquet file and feeds the table within to minimap2_table().
    
    Args :
        file (str or Path) : path to the parquet file to be mapped.
        path_ref (str) : path to the sequence reference to be used.
        preset (str) : the preset to use for minimap2. Defaults to 'splice', which is designed to be splice-aware.
        resume (bool) : whether to resume mapping from a paused or broken run. This will only map files that don't already have minimap information in them, so it can't be used to continue a re-mapping session, as the files will all still have the original minimap data. Defaults to False.
        overwrite (bool) : whether to allow for overwriting mapped files. If True, it will remove any existing minimap data and continue with normal mapping. If False, it will skip over any files with minimap data.
    """
    print(file, '                  ')
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
    """
    Goes through the parquet dataset in dataset_dir and runs everything through minimap2 to assign sequencing reads to genes. Opens all files individually and runs minimap2_table_from_file(), which also saves the results to disk.
    
    Args :
        dataset_dir (str) : the path to the folder containing the parquet files to be mapped. Should be path_out/pa_dataset/ where path_out is the same as what was used in build_parquet_dataset_from_sam.
        path_ref (str) : path to the sequence reference to be used. Currently intended to work with a transcripomic reference, as count_mapped_reads will only work with that.
        preset (str) : the preset to use for minimap2. Defaults to 'splice', which is designed to be splice-aware.
        resume (bool) : whether to resume mapping from a paused or broken run. This will only map files that don't already have minimap information in them, so it can't be used to continue a re-mapping session, as the files will all still have the original minimap data. Defaults to False.
        overwrite (bool) : whether to allow for overwriting mapped files. If True, it will remove any existing minimap data and continue with normal mapping. If False, it will skip over any files with minimap data.
        workers (int) : number of parallel mapping processes to run. Note that minimap2 has high memory requirements, appearing to need about 8-12G per process to be stable. Using less will possisbly result in broken runs which are not currently set to re-run. Defaults to 4.
    """
    files = [x for x in Path(dataset_dir).iterdir() if x.is_file()]
    split = math.ceil( len(files) / workers )
    files_chunks = np.array_split(np.array(files), split)
    for chunk in files_chunks :
        with concurrent.futures.ProcessPoolExecutor( max_workers=workers ) as executor :
            futures = [ executor.submit( minimap2_table_from_file, file, path_ref, preset=preset, resume=resume, overwrite=overwrite ) for file in chunk]
            concurrent.futures.wait( futures )
    return

def count_mapped_reads(dataset_dir, path_ref, path_out_csv, sample_dict, run_label = "None") :
    """
    Counts the mapped reads in the dataset_dir and tallies up gene counts, collapsing gene and transcript variants into one count per gene ID. Gene IDs are defined by the name in parentheses in the minimap reference file. Saves the resulting counts as csv.
    
    Args :
        dataset_dir (str) : the path to the folder containing the parquet files to be counted. Should be path_out/pa_dataset/ where path_out is the same as what was used in build_parquet_dataset_from_sam.
        path_ref (str) : path to the sequence reference to be used. Currently intended to work with a transcripomic reference, essentially in the format of a fastq file where each row is a transcript with a contig ID, then sequence, then some tags. The critical information is the gene ID in parentheses, which must be the last parentheses in the first field of each row. ie XC001.4 Gene Name Etc. (version 1) (Gene ID), sequence, other information.
        path_out_csv (str) : the full path to where the results should be saved as csv.
        sample_dict (dict) : a dictionary to define which barcodes belong to what sample names. Must be in the form of 'barcode ID' : 'sample name'.
        run_label (str) : a label that is assigned to all the reads in this dataset to potentially differentiate them from reads from other runs that might share the same barcode and sample names. Useful for combining reads from two runs with the same samples or if a run fails and is restarted.
    """
    ref = local_io.read_fastx(path_ref)
    contigs_to_gene_IDs = {}
    gene_IDs = []
    for key in ref.keys() :
        comma_split = key.split(', ')
        space_split = comma_split[0].split(' ')
        left_para_split = key.split('(')
        right_para_split = left_para_split[ len(left_para_split) - 1 ].split(')')
        gene_ID = right_para_split[0]
        contig_ID = space_split[0]
        contigs_to_gene_IDs[contig_ID] = gene_ID
        gene_IDs.append(gene_ID)
    
    table = pq.read_table(dataset_dir, columns=['barcode_ID', 'minimap2_ctg', 'minimap2_mapq'])
    barcodes = table.column('barcode_ID').unique().to_pylist()
    counts = False
    for barcode in barcodes :
        if barcode not in ['None matched', 'Multiple'] :
            counts_in_barcode_array = table.filter( pc.field('barcode_ID') == barcode ).column('minimap2_ctg').drop_null().value_counts()
            counts_in_barcode_dict = dict.fromkeys(gene_IDs, 0)
#             total_counts = 0
            for count in counts_in_barcode_array :
                gene_ID = contigs_to_gene_IDs[str(count[0])]
                counts_in_barcode_dict[gene_ID] += count[1].as_py()
#                 total_counts += count[1].as_py()
            for key in counts_in_barcode_dict :
                counts_in_barcode_dict[key] = [ counts_in_barcode_dict[key] ]
            
            counts_in_barcode_table = pa.table(counts_in_barcode_dict).add_column(0, 'barcode_ID', [[barcode]]).add_column(0, 'sample', [[sample_dict[barcode]]]).add_column(0, 'run_label', [[run_label]])
            if not counts :
                counts = counts_in_barcode_table
            else :
                counts = pa.concat_tables([counts, counts_in_barcode_table])
    
    csv.write_csv(counts, path_out_csv)
    print('done')
    return

def qc_metrics(path_dataset) :
    table = pq.read_table(path_dataset, columns = ['seq_len', 'barcode_score', 'barcode_ID', 'SSP_edit_distance'])
    qc = {}
    hist_max = 5000
    fig1, qc['alignment_scores_hist'] = plt.subplots(2,2, sharex=True, sharey=True)

    group_1 = table.filter(pc.field('barcode_score') <= 20).filter(pc.field('barcode_ID') != 'Multiple').filter(pc.field('SSP_edit_distance') <= 5)
    qc['alignment_scores_hist'][0,0].hist(group_1.column('seq_len').to_pylist(), bins=50, range=(0, hist_max))
    qc['alignment_scores_hist'][0,0].set_title('Passed Barcode and Strand Switch', {'fontsize': 10})
    group_2 = table.filter(pc.field('barcode_score') <= 20).filter(pc.field('barcode_ID') != 'Multiple').filter(pc.field('SSP_edit_distance') > 5)
    qc['alignment_scores_hist'][0,1].hist(group_2.column('seq_len').to_pylist(), bins=50, range=(0, hist_max))
    qc['alignment_scores_hist'][0,1].set_title('Passed Barcode, failed Strand Switch', {'fontsize': 10})
    group_3 = table.filter(pc.field('barcode_score') > 20).filter(pc.field('SSP_edit_distance') <= 5)
    qc['alignment_scores_hist'][1,0].hist(group_3.column('seq_len').to_pylist(), bins=50, range=(0, hist_max))
    qc['alignment_scores_hist'][1,0].set_title('Failed Barcode, passed Strand Switch', {'fontsize': 10})
    group_4 = table.filter(pc.field('barcode_score') > 20).filter(pc.field('SSP_edit_distance') > 5)
    qc['alignment_scores_hist'][1,1].hist(group_4.column('seq_len').to_pylist(), bins=50, range=(0, hist_max))
    qc['alignment_scores_hist'][1,1].set_title('Failed Barcode and Strand Switch', {'fontsize': 10})

    barcodes = table.column('barcode_ID').unique().to_pylist()
    for barcode in barcodes :
        print(barcode + ' passed barcode and SSP: ' + str(group_1.filter(pc.field('barcode_ID') == barcode).num_rows))
        print(barcode + ' passed barcode, failed SSP: ' + str(group_2.filter(pc.field('barcode_ID') == barcode).num_rows))
        print(barcode + ' failed barcode, passed SSP: ' + str(group_3.filter(pc.field('barcode_ID') == barcode).num_rows))
        print(barcode + ' failed barcode and SSP: ' + str(group_4.filter(pc.field('barcode_ID') == barcode).num_rows))
    return qc