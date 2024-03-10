import sys
import os.path
import shutil
import subprocess
import time
import math
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import concurrent
import concurrent.futures
import pickle
import pysam
from bioinfokit import analys, visuz
import pydeseq2
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from NanoporeAnalysis import local_io
from NanoporeAnalysis import utils
from NanoporeAnalysis import align
from pod5.tools import pod5_view
from pod5.tools import pod5_subset
from pod5.tools import pod5_filter
from importlib import reload
reload(align)

"""
Todo:
Give Dorado options. ie Simplex calling
Filter out tRNA and rRNA? Also find the % of total reads that are t and r RNA
Find better p value calculation for DEG analysis
Combine barcode and SS alignment? align a combined sequence while setting the gap penalty to 0
drop gap penalty for alignment
FSA package
parasailer for alignment?
edlib for a crude alignment?
isONclust for denovo clustering
hidden markov model for tRNA and rRNA detection
Clean up the dorado function(s)
barcode alignment gives full scoring to the polyA, which means that a polyA region without a properly identifiable barcode can pass the threshold
change mappy to run once. The index building may be the limiting factor.
generalize functions a bit more

QC metrics:
    plot barcode total pass/fail numers per barcode or sample in a stacked bar plot.
    Compare between barcodes of a sample
    Add more options to the alignment viewer
    Read depth per gene
    Do a pileup map for a gene to see where reads lie
    count reads starting/stopping in the middle of a gene
    
Analysis functions:
    Diff exp between samples/barcodes for all features referenced or for specified features.
    Diff exp of mRNA variants? Idea: use CIGAR ouput from mappy alignment if possible. Check for large deletions. Spladder uses a similar method 
    PCA
"""

def run_dorado_slurm(path_dorado, path_model, path_out, path_data, account, mail, workers=1, threads=1, overwrite=False, skip_split=False) :
    """
    For running dorado through SLURM cluster job scheduling.
    """
    if path_data != None:
        if not ( Path(path_data).is_dir() or os.listdir(path_data) ) :
            raise FileNotFoundError("Error: Data path is invalid or doesn't exist. Please point to an existing directory containing data.")
    
    Path(path_out + '/pod5_sam').mkdir(parents=True, exist_ok=True)
    Path(path_out + '/fastq').mkdir(parents=True, exist_ok=True)
    Path(path_out + '/logs').mkdir(parents=True, exist_ok=True)

    if (os.listdir(path_out + '/pod5s_by_channel') or os.listdir(path_out + '/fastq')) and overwrite == False :
        raise FileExistsError("Error: output path (fastq or pod5_sam) may contain data. Please point to a fresh directory or specify overwrite=True")
    elif overwrite == True :
        shutil.rmtree(path_out + '/fastq')
        shutil.rmtree(path_out + '/pod5_sam')
        Path(path_out + '/fastq').mkdir(parents=True, exist_ok=True)
        Path(path_out + '/pod5_sam').mkdir(parents=True, exist_ok=True)

    if skip_split == False :
        if os.listdir(path_out + '/pod5s_by_channel') and overwrite == False :
            raise FileExistsError("Error: output path (split_pod5s) may contain data. Please point to a fresh directory or specify overwrite=True")
        elif overwrite == True :
            shutil.rmtree(path_out + '/pod5s_by_channel')
            Path(path_out + '/pod5s_by_channel').mkdir(parents=True, exist_ok=True)
        print("started making the table")
        pod5_view.view_pod5([Path(path_data)], Path(path_out), include = "read_id, channel", force_overwrite=True, threads=threads)
        print("finished making the table")
        
        pod5_view = pd.read_table(Path(path_out + "/view.txt"), sep='\t')
#         for channel in pod5_view['channel'].unique() :
#             pod5_view[pod5_view['channel'] == channel]['read_id'].to_csv(Path(path_out + "/view_current.txt"), index = False, sep=' ', header=False)
#             pod5_filter.filter_pod5([Path(path_data)], Path(path_out + '/split_pod5s'), threads = threads, force_overwrite = True, ids = Path(path_out + "/view.txt"))
        
#         pod5_subset.subset_pod5([Path(path_data)], Path(path_out + '/split_pod5s'), threads = threads, force_overwrite = True, columns = ["channel"], table = Path(path_out + "/view.txt"))
        print("finished subsetting")
        
        channel_folders = [x for x in Path(path_out + '/pod5s_by_channel').iterdir() if x.is_dir()]
        arrays = np.array_split(np.array(channel_folders), workers)
        for i in range(workers) :
            Path(path_out + '/split_pod5s_byworker/' + str(i)).mkdir(parents=True, exist_ok=True)
            for path in arrays[i] :
                for file in [x for x in path.glob('*') if x.is_file()]:
                    shutil.copy(str(file), str(path_out + '/split_pod5s_byworker/' + str(i)))
#     channel_folders = [x for x in Path(path_out + '/pod5s_by_channel').iterdir() if x.is_dir()]
    for i in range(workers) :
        script = [
            '#!/bin/bash\n',
            str('#SBATCH --account=' + account + '\n'),
            str('#SBATCH --job-name=dorado_' + str(i) + '\n'),
            '#SBATCH --nodes=1\n',
            '#SBATCH --ntasks=1\n',
            '#SBATCH --cpus-per-task=2\n',
            '#SBATCH --mem=20G\n',
            '#SBATCH --gres=gpu:2\n',
            '#SBATCH --time=02:00:00\n',
            str('#SBATCH --output=' + path_out + '/logs' + '/dorado_' + str(i) + '.out' + '\n'),
            '#SBATCH --mail-type=None\n',
            str('#SBATCH --mail-user=' + mail + '\n'),
            str(str(path_dorado) + " duplex " + path_model + ' --emit-sam ' + str(path_out + '/split_pod5s_byworker/' + str(i)) + '.pod5 > ' + str(path_out) + '/pod5_sam' + '/' + str(i) + '.sam')
        ]
        with open(str(path_out + '/logs' + '/dorado_' + str(i) + '.sh'), 'w') as handle:
            handle.writelines(script)
        args = ['sbatch', str(path_out + '/logs' + '/dorado_' + str(i) + '.sh')]
        subprocess.run(args)

    return

def run_dorado(path_dorado, path_model, path_out, path_data, workers=1, skip_split=False, overwrite=False) :
    """
    Run Dorado duplex basecalling on provided data.
    Args :
        path_dorado (str) : path to dorado executable
        path_model (str) : path to dorado model
        workers (int) : number of dorado processes to spawn. In this function, this just serves to split the data up for manageability and stability reasons.
        skip_split (boolean) : whether or not to skip splitting the pod5 files. Set to true if the pod5s have previously been split for duplex calling.
    Output: None
    """
    if path_data != None:
        if not ( Path(path_data).is_dir() or os.listdir(path_data) ) :
            raise FileNotFoundError("Error: Data path is invalid or doesn't exist. Please point to an existing directory containing data.")

    Path(path_out + '/pod5_sam').mkdir(parents=True, exist_ok=True)
    Path(path_out + '/fastq').mkdir(parents=True, exist_ok=True)
    Path(path_out + '/split_pod5s').mkdir(parents=True, exist_ok=True)

    if (os.listdir(path_out + '/split_pod5s') or os.listdir(path_out + '/fastq')) and overwrite == False :
        raise FileExistsError("Error: output path (fastq or pod5_sam) may contain data. Please point to a fresh directory or specify overwrite=True")
    elif overwrite == True :
        shutil.rmtree(path_out + '/fastq')
        shutil.rmtree(path_out + '/pod5_sam')
        Path(path_out + '/fastq').mkdir(parents=True, exist_ok=True)
        Path(path_out + '/pod5_sam').mkdir(parents=True, exist_ok=True)

    if skip_split == False :
        if os.listdir(path_out + '/split_pod5s') and overwrite == False :
            raise FileExistsError("Error: output path (split_pod5s) may contain data. Please point to a fresh directory or specify overwrite=True")
        elif overwrite == True :
            shutil.rmtree(path_out + '/split_pod5s')
            Path(path_out + '/split_pod5s').mkdir(parents=True, exist_ok=True)
        pod5_view.view_pod5([Path(path_data)], Path(path_out), include = "read_id, channel", force_overwrite=True)
        pod5_subset.subset_pod5([Path(path_data)], Path(path_out + '/split_pod5s'), columns = ["channel"], table = Path(path_out + "/view.txt"))
        files = [x for x in Path(path_out + '/split_pod5s').iterdir() if x.is_file()]
        arrays = np.array_split(np.array(files), workers)
        for i in range(workers) :
            Path(path_out + '/split_pod5s_byworker/' + str(i)).mkdir(parents=True, exist_ok=True)
            for path in arrays[i] :
                shutil.copy(str(path), str(path_out + '/split_pod5s_byworker/' + str(i)))

    for i in range(workers) :
        args = [str(path_dorado + " duplex " + path_model + ' ' + str(path_out + '/split_pod5s_byworker/' + str(i)) + ' > ' + path_out + '/pod5_sam' + '/' + str(i) + '.bam')]
        subprocess.run(args, shell=True)

    files = [x for x in Path(path_out + '/pod5_sam').iterdir() if x.is_file()]
    Path(path_out + '/pod5_sam' + '/sorted_by_name/').mkdir(parents=True, exist_ok=True)
    
    for file in files : 
        print("converting ", file, " to fastq")
        pysam.sort('-o', str(path_out + '/pod5_sam' + '/sorted_by_name/' + file.stem + '_sorted.bam'), '-n', str(file) )
        Path(str(path_out + '/fastq' + '/'+ file.stem + '.fq')).touch()
        pysam.fastq('-0', str(path_out + '/fastq' + '/'+ file.stem + '.fq'), '-T', 'dx', str(path_out + '/pod5_sam' + '/sorted_by_name/' + file.stem + '_sorted.bam'))

    return

def debarcode(path_barcodes, strand_switch_primer, path_out, path_fastq=None, overwrite=False, workers=None, filter_barcode_score = 0, filter_barcode_distance = 0, filter_strand_switch_score = None, threshold_barcode = None, threshold_strand_switch = None) :
    """
    aligns barcodes and strand_switch primer, then saves to file by barcode and whether the fliter was met.

    args:
        path_barcodes (str) : file path for a csv containing a column defining barcode primer names and a second column defining their sequences
        strand_switch_primer (str) : sequence for the strand_switch_primer
        path_fastq (str) : file path to a folder with fastx files. For use with data basecalled on machine/elsewhere
        filter_barcode_score (int) : filter for barcode alignment. slightly speeds up alignment by 10-20% if set to ~90% of the max alignment score
        filter_barcode_distance (int) : similar to filter_barcode_score, but for min aligned distance. Less effective
        filter_strand_switch_score : minimum score for the strand_switch_primer alignment. Below this filter, a full alignment is done for the barcode. By default, this is set to half the maximum score for the given sequence
        threshold_barcode : minimum score for barcode alignment when sorting into files. Reads not meeting this filter get put into 'unclassified.fq'. By default, this is set to half the maximum score for the given sequence.
        threshold_strand_switch : minimum score for ss primer alignment when sorting into files. Reads not meeting this filter get put into a '*_incomplete.fq' file for their given barcode. By default, this is set to half the maximum score for the given sequence.
    Output:
        None: Writes to file.
    """
    if path_fastq == None : 
        path_fastq = path_out + '/fastq'
    if os.listdir(path_fastq) == False :
        raise FileNotFoundError("No fastq files found. Please run basecalling to generate fastq files or set path_fastq to a folder containing fastq files.")

    Path(path_out + "/temp_aligned_reads").mkdir(parents=True, exist_ok=True)
    Path(path_out + '/debarcoded').mkdir(parents=True, exist_ok=True)
    Path(path_out + '/fasta').mkdir(parents=True, exist_ok=True)
    if ( os.listdir(path_out + "/temp_aligned_reads") or os.listdir(path_out + '/debarcoded') or os.listdir(path_out + '/fasta') ) and overwrite == False :
        raise FileExistsError("Error: output path may contain data. Please point to a fresh directory or specify overwrite=True")
    elif overwrite == True :
        shutil.rmtree(path_out + "/temp_aligned_reads")
        shutil.rmtree(path_out + '/debarcoded')
        shutil.rmtree(path_out + '/fasta')
        Path(path_out + "/temp_aligned_reads").mkdir(parents=True, exist_ok=True)
        Path(path_out + '/debarcoded').mkdir(parents=True, exist_ok=True)
        Path(path_out + '/fasta').mkdir(parents=True, exist_ok=True)

    if filter_strand_switch_score == None : filter_strand_switch_score = len(strand_switch_primer) * 0.5 * 2
    if threshold_strand_switch == None : threshold_strand_switch = len(strand_switch_primer) * 0.5 * 2

    #load barcodes
    barcodes = pd.read_csv(path_barcodes, names=['name', 'seq'])
    barcodes_clean = barcodes.copy()
    barcodes_clean['seq'] = barcodes['seq'].str.replace('[^ATCGatcg]', '', regex=True).apply(utils.reverse_complement)
    barcodes_clean = barcodes_clean.set_index('name').T.to_dict('records')[0]

    #find data
    files = [x for x in Path(path_fastq).iterdir() if x.is_file()]

    start_time = time.time()
    aligned_reads = []

    dump_index = 1
    read_count = 0

    for file in files:

        file_name = str(Path(file).stem)

        print('Debarcoding file ' + file_name + ';   File ', files.index(file) + 1, 'out of ', len(files))
        data = local_io.read_fastx(file)

        executor = concurrent.futures.ProcessPoolExecutor( max_workers=workers )
        if 'Tags' in data[list(data.keys())[0]].keys() :
            futures = [ executor.submit( align.get_best_barcode, data[read_id]['Sequence'], barcodes_clean, strand_switch_primer, filter_strand_switch_score, read_id, data[read_id]['Tags'], filter_barcode_score, filter_barcode_distance) for read_id in data if not 'dx:i:-1' in data[read_id]['Tags']]
        else :
            futures = [ executor.submit( align.get_best_barcode, data[read_id]['Sequence'], barcodes_clean, strand_switch_primer, filter_strand_switch_score, read_id, None, filter_barcode_score, filter_barcode_distance) for read_id in data]
        concurrent.futures.wait( futures )
        for f in futures: aligned_reads.append(f.result())

        if sys.getsizeof(aligned_reads) > 100000 or files.index(file) + 1 >= len(files) :

            print('Saved data to temp file ' + str(dump_index))
            with open(path_out + "/temp_aligned_reads" +'/saved_read_alignments_' + str(dump_index) + '.txt', 'wb') as handle:
                pickle.dump(aligned_reads, handle)
            dump_index += 1
            read_count += len(aligned_reads)
            aligned_reads = []

    print("Debarcoded ", read_count, " reads in: ", time.time()-start_time, " sec")
    
    # Classify each sequence based on how well the barcode and opposite primer are aligned
    files_aligned = [str(x) for x in Path(path_out + "/temp_aligned_reads").iterdir() if x.is_file()]
    
    print("Sorting and saving debarcoded reads.")
    
    for i in range(len(files_aligned)) :
        
        seq_by_barcode = {'unclassified' : {}}
        data_by_barcode = {'unclassified' : {}}

        for primer_name in barcodes_clean.keys():
            seq_by_barcode[primer_name] = {}
            seq_by_barcode[primer_name + '_incomplete'] = {}
            data_by_barcode[primer_name] = {}
            data_by_barcode[primer_name + '_incomplete'] = {}

        with open(files_aligned[i], 'rb') as handle:
            aligned_reads = pickle.load(handle)
        print("sorting ", i)

        for a in aligned_reads:

            score_barcode = a['alignment']['optimal_alignment_score']
            score_strand_switch = a['strand_switch_primer_align']['optimal_alignment_score']

            primer_name = a['primer_name']

            key, bio_seq = a['seq_id'], a['biological_sequence']
            del a['seq_id']

            if threshold_barcode == None : 
                tmp_threshold_barcode = len(a['alignment']['query_sequence']) * 0.5 * 2 
            else :
                tmp_threshold_barcode = threshold_barcode

            a['threshold_barcode'] = tmp_threshold_barcode
            a['threshold_strand_switch'] = threshold_strand_switch

            if score_barcode > tmp_threshold_barcode:
                if score_strand_switch > threshold_strand_switch:
                    seq_by_barcode[primer_name][key] = bio_seq
                    data_by_barcode[primer_name][key] = a
                else:
                    seq_by_barcode[primer_name + '_incomplete'][key] = bio_seq
                    data_by_barcode[primer_name + '_incomplete'][key] = a

            else:
                seq_by_barcode['unclassified'][key] = bio_seq
                data_by_barcode['unclassified'][key] = a

        # Save data into multiple FASTA files

        for category in seq_by_barcode.keys():
            fasta_file = Path(path_out + '/fasta') / (category + '_' + str(i) + '.fa')
            local_io.write_fasta(fasta_file, seq_by_barcode[category], append=False)
            df = pd.DataFrame.from_dict(data_by_barcode[category], orient = 'index')
            df.to_csv(str(path_out + '/debarcoded' + '/' + category + '_' + str(i) + '.csv'), header=True, index_label='read_id')
        print("Saved ", i, " to disk")
        
    return

def minimap2(path_minimap2, path_ref, path_out, workers = 3, overwrite=False, finish=False) :
    """
    Calls Minimap2 to map all reads in given path based on a given reference. Saves mapped reads into new folder.

    Args:
        path_minimap2 (str) = path to minimap2 executable
        path_ref (str) = path to reference file

    Output: None
    """
    Path(path_out + '/combined_sam').mkdir(parents=True, exist_ok=True)
    Path(path_out + '/combined_bam').mkdir(parents=True, exist_ok=True)
    if (os.listdir(path_out + '/combined_sam') or os.listdir(path_out + '/combined_bam')) and overwrite == False and finish == False :
        raise FileExistsError("Error: output path may contain data. Please point to a fresh directory or specify overwrite=True")
    elif overwrite == True :
        shutil.rmtree(path_out + '/combined_sam')
        shutil.rmtree(path_out + '/combined_bam')
        Path(path_out + '/combined_sam').mkdir(parents=True, exist_ok=True)
        Path(path_out + '/combined_bam').mkdir(parents=True, exist_ok=True)

    #Minimap2 alignment
    if finish == True :
        files_ran = [x.stem for x in Path(path_out + '/combined_bam').iterdir() if x.is_file() and x.suffix == '.bam' ]
        files_debarcoded = [str(x) for x in Path(path_out + '/combined_fasta').iterdir() if x.is_file() and x.suffix == '.fasta' and 'unclassified' not in x.stem and x.stem not in files_ran ]
    else :
        files_debarcoded = [str(x) for x in Path(path_out + '/combined_fasta').iterdir() if x.is_file() and x.suffix == '.fasta' and 'unclassified' not in x.stem ]
    
    executor = concurrent.futures.ProcessPoolExecutor( max_workers = workers )
    futures = [ executor.submit( minimap2_util, path_minimap2, path_ref, file, path_out ) for file in files_debarcoded]
    concurrent.futures.wait( futures )
    
#     for file in files_debarcoded :
#         print('mapping file ', Path(file).stem)
#         args = [str(path_minimap2 + ' -ax splice ' + path_ref + ' ' + file + ' > ' + path_out + '/combined_sam' + '/' + Path(file).stem + '.sam')]
#         subprocess.run(args, shell=True)
#         pysam.sort("-o", str(path_out + '/combined_bam' + '/' + Path(file).stem + '.bam'), path_out + '/combined_sam' + '/' + Path(file).stem + '.sam')
#         pysam.index(str(path_out + '/combined_bam' + '/' + Path(file).stem + '.bam'))

    print("Done mapping")
    
    return

def minimap2_util(path_minimap2, path_ref, file, path_out) :
    print('mapping file ', Path(file).stem)
    args = [str(path_minimap2 + ' -ax splice ' + path_ref + ' ' + file + ' > ' + path_out + '/combined_sam/' + Path(file).stem + '.sam')]
    subprocess.run(args, shell=True)
    pysam.sort("-o", str(path_out + '/combined_bam/' + Path(file).stem + '.bam'), path_out + '/combined_sam/' + Path(file).stem + '.sam')
    pysam.index(str(path_out + '/combined_bam/' + Path(file).stem + '.bam'))
    return

def count_reads(path_ref_file, samples, path_out, workers = 4, count = True, overwrite=False) :

    """
    Makes read counts for mapped reads in path_out against a provided .bed file. Outputs to a new folder.

    Args:
        path_bed_file (str) = path to bed file
        samples (dict) = dictionary defining the barcodes for each sample. Follow format {sample name : [barcodes]...}
        count (boolean) = whether or not to perform the read counting. Set to False if counting has already been done previously, ie to re-define sample barcodes.
    Outputs: None
    """

    if count == True :
        bam_files = [ x for x in Path(path_out + '/combined_bam').iterdir() if x.is_file() and x.suffix == '.bam' and x.stem != 'unclassified' and 'incomplete' not in x.stem ]
        start_time = time.time()
        counts = pd.DataFrame()
        seqs = pd.DataFrame()
        
        executor = concurrent.futures.ProcessPoolExecutor( max_workers = workers )
        futures = [ executor.submit( count_util, file, path_out, path_ref_file, overwrite ) for file in bam_files ]
        concurrent.futures.wait( futures )
        
        for future in futures :
            result = future.result()
            if counts.empty :
                counts = pd.DataFrame(data = result[0], index = [result[1]] ).rename_axis('barcode')
            else :
                counts = pd.concat([counts, pd.DataFrame(data = result[0], index = [result[1]] ).rename_axis('barcode')])
        
        counts.to_csv(str(path_out + '/counts.csv'), header=True, index_label='barcode')
        
        print('finished counting in ' + str( time.time()-start_time ) + 'sec')
    else:
        counts = pd.read_csv(str(path_out + '/counts.csv'), index_col = 'barcode')
    
    counts_by_sample = counts
    counts_by_sample.insert(0,'sample', '~')
    for sample in samples :
        for barcode in samples[sample] :
            for index in counts_by_sample.index :
                if barcode == index :
                    counts_by_sample.at[index, 'sample'] = sample
#     print(counts_by_sample)
    counts_by_sample.reset_index().set_index(['sample', 'barcode']).to_csv(str(path_out + '/counts_by_sample.csv'), header=True, index_label=['sample', 'barcode'])

    print('done')

    return

def count_util(file, path_out, path_ref_file, overwrite) :
    print("counting " + file.stem)
    bam_file = pysam.AlignmentFile(file, 'rb')
    sample_name = file.stem
    with open(path_out + '/combined_metadata' + '/' + sample_name + '.csv') as handle :
        total_reads = sum(1 for row in handle) - 1
    """
    meta_data = pd.read_csv(path_out + '/combined_metadata' + '/' + sample_name + '.csv', dtype={'minimap_alignment' : str}).set_index('read_id') #note: take out dtype argument after rerunning
    if 'minimap_alignment' in meta_data.columns and overwrite==False :
        raise ValueError("Counting has already been done for this file. Please set overwrite=True if you'd like to continue.")
    elif 'minimap_alignment' in meta_data.columns :
        meta_data.drop(columns='minimap_alignment', inplace=True)
        meta_data.insert(len(meta_data.columns), 'minimap_alignment', 'Not mapped')
    else :
        meta_data.insert(len(meta_data.columns), 'minimap_alignment', 'Not mapped')"""
        
#     total_reads = len(meta_data.index)
    
    ref = local_io.read_fastx(path_ref_file)
    contigs = []
    counts = {}

    for key in ref.keys() :
        comma_split = key.split(', ')
        contig_name_split = comma_split[0].split(' ')
        contig_id = contig_name_split[0]
        gene_name = ' '.join(contig_name_split[1:])
        contigs.append( [ contig_id, gene_name ] )
        counts[gene_name] = 0

    for contig in contigs :
        bam_iter = bam_file.fetch(contig[0])
        read_count = 0
        for x in bam_iter:
            if not ( x.is_secondary or x.is_supplementary ) :
                try :
#                     meta_data.at[x.query_name, 'minimap_alignment'] = [contig[0] + ' ' + contig[1], x.reference_start, x.cigarstring, x.reference_end]
                    read_count += 1 #some attempts at adding to meta_data raise an error. If these events should still result in an increase to read_count, move this line outside of the try block
                except:
                    pass
        
        
        RPM_multiplier =  total_reads# / 1000000
        read_count = read_count / RPM_multiplier
        counts[contig[1]] += read_count
#         counts_by_feature.append([contig[0], read_count])
        
#     with open(path_bed_file, 'r') as bed_in:
#         bed_lines = bed_in.readlines()
#     for line in bed_lines :
#         cols = line.split('\t')
#         bam_iter = bam_file.fetch(cols[0], int(cols[1]), int(cols[2]))
        
#         read_count = 0
#         read_seq = []
#         for x in bam_iter:
#             if not ( x.is_secondary or x.is_supplementary ) :
#                 try :
#                     meta_data.at[x.query_name, 'minimap_alignment'] = [cols[3], x.reference_start, x.cigarstring, x.reference_end]
# #                     seq = x.query_sequence
# #                     if isinstance(seq, str): # why are some empty?..
# #                         read_seq.append(' '.join([x.query_name, seq]))
#                     read_count += 1 #some attempts at adding to meta_data raise an error. If these events should still result in an increase to read_count, move this line outside of the try block
#                 except:
#                     pass
        
#         read_seq_str = '\n'.join(read_seq)
#         seqs_by_feature.append(read_seq_str)
    
    bam_file.close()
#     meta_data.to_csv(str(path_out + '/combined_metadata' + '/' + sample_name + '.csv'), header=True, index_label='read_id')
    return [counts, sample_name]

def qc_metrics(path_out, on='library', samples=None, hist_max=3000) :
    """
    Calculate quality control metrics.
    Args:
        alignment_scores_hist_max (int) : set the max x-axis value for the read length histograms for legibility. Default is 3000 bp.
    Output:
        None, displays histograms for read lengths according to passing barcode and strand_switch thresholds. Uses the thresholds set during debarcoding. 
    """
    Path(path_out + '/qc').mkdir(parents=True, exist_ok=True)

    qc = {}
    files = [str(x) for x in Path(path_out + '/combined_metadata').iterdir() if x.is_file()]
    df = pd.DataFrame()
    for file in files :
        df = pd.concat([df, pd.read_csv(file).set_index('read_id')], axis = 0)
    qc['read_count'] = len(df.index)
    #print(df)

    if on == 'library' : 
        qc = metrics_util(df, hist_max, title='Library')
    if on == 'sample' :
        if samples == None :
            raise ValueError("Please provide a 'samples' list")
        for sample in sample :
            df_tmp = df[df['primer_name'] in sample]
            qc[sample] = metrics_util(df_tmp, hist_max, title=sample)
    if on == 'barcodes' :
        for barcode in df['primer_name'].drop_duplicates() :
            df_tmp = df[df['primer_name'] == barcode]
            qc[barcode] = metrics_util(df_tmp, hist_max, title=barcode)

    return

def metrics_util(path_out, df, hist_max, title) :
    qc = {}
    fig1, qc['alignment_scores_hist'] = plt.subplots(2,2, sharex=True)
    qc['alignment_scores_hist'][0,0].hist(df[(df['barcode_score'] >= df['threshold_barcode']) & (df['strand_switch_score'] >= df['threshold_strand_switch'])]['read_length'], bins=50, range=(0, hist_max))
    qc['alignment_scores_hist'][0,0].set_title('Passed Barcode and Strand Switch', {'fontsize': 10})
    qc['alignment_scores_hist'][0,1].hist(df[(df['barcode_score'] >= df['threshold_barcode']) & (df['strand_switch_score'] < df['threshold_strand_switch'])]["read_length"], bins=50, range=(0, hist_max))
    qc['alignment_scores_hist'][0,1].set_title('Passed Barcode, failed Strand Switch', {'fontsize': 10})
    qc['alignment_scores_hist'][1,0].hist(df[(df['barcode_score'] < df['threshold_barcode']) & (df['strand_switch_score'] >= df['threshold_strand_switch'])]["read_length"], bins=50, range=(0, hist_max))
    qc['alignment_scores_hist'][1,0].set_title('Failed Barcode, passed Strand Switch', {'fontsize': 10})
    qc['alignment_scores_hist'][1,1].hist(df[(df['barcode_score'] < df['threshold_barcode']) & (df['strand_switch_score'] < df['threshold_strand_switch'])]["read_length"], bins=50, range=(0, hist_max))
    qc['alignment_scores_hist'][1,1].set_title('Failed Barcode and Strand Switch', {'fontsize': 10})
    fig1.suptitle(title)

    qc['passed_barcode_passed_strand_switch'] = len(df[(df['barcode_score'] >= df['threshold_barcode']) & (df['strand_switch_score'] >= df['threshold_strand_switch'])])
    qc['passed_barcode_failed_strand_switch'] = len(df[(df['barcode_score'] >= df['threshold_barcode']) & (df['strand_switch_score'] < df['threshold_strand_switch'])])
    qc['failed_barcode_passed_strand_switch'] = len(df[(df['barcode_score'] < df['threshold_barcode']) & (df['strand_switch_score'] >= df['threshold_strand_switch'])])
    qc['failed_barcode_failed_strand_switch'] = len(df[(df['barcode_score'] < df['threshold_barcode']) & (df['strand_switch_score'] < df['threshold_strand_switch'])])

    return qc

def diff_exp(path_out, samples, meta_data, design_factor) :
    """
    Do a differential expression analysis, display a volcano plot, and save the diff exp gene results.
    Args:
        samples (list of str) : List of sample names to include in the analysis
        meta_data (dataframe) : meta data for the samples with indexes of sample names and named columns referring to the experimental conditions that differentiate them, ie cell type, growth condition, treatment.
        design_factor (str) : the experimental condition to do diff exp on. Must match one of the column names in meta_data.
    Output:
        Shows resulting volcano plot and saves analysis to file under path_out with a name referring to the samples and design_factor provided.
    """
    counts_by_sample = pd.read_csv(path_out + '/counts_by_sample.csv', index_col = ['sample', 'barcode']).reset_index('sample')
    volcano_data = counts_by_sample[counts_by_sample['sample'].isin(samples)].drop('sample', axis=1)
    volcano_data = volcano_data.T[volcano_data.min(axis=0) >= 20].T
    volcano_meta = meta_data.loc[samples, :]
#     print(volcano_meta)
    dds = DeseqDataSet(
        counts=volcano_data,
        metadata=volcano_meta,
        design_factors=design_factor,
        quiet=False
    )
    dds.deseq2()
    stat_res = DeseqStats(dds, quiet=False)
    stat_res.summary()
    df = stat_res.results_df.reset_index().dropna().copy()
    pd.DataFrame.from_dict(df).to_csv(str(path_out + '/diff_seq_' + str(samples) + '_' + design_factor + '.csv'), header=True)
    df = df.drop(df['pvalue'].idxmin(axis=0))

    visuz.GeneExpression.volcano(df=df, lfc='log2FoldChange', pv='pvalue', show=True, lfc_thr=(0.5,0.5), pv_thr=(0.05,0.05))
    return

def view_alignments(path_out, barcodes=None, thresholds='any', num=10, funcs=None, reset=False, metadata=None) :

    assert thresholds in ['any', 'neither', 'both', 'barcode', 'strand_switch'], "Please specify thresholds = 'any', 'neither', 'both', 'barcode', or 'strand_switch'"

    if meta_data != None or reset == True:
        meta_data = pd.DataFrame()
        files = [str(x) for x in Path(path_out + '/debarcoded').iterdir() if x.is_file()]
        print("Importing metadata...")
        for file in files :
            meta_data = pd.concat([meta_data, pd.read_csv(file, dtype={'minimap_alignment' : str}).set_index('read_id')], axis = 0) #note: take out dtype argument after rerunning
    else :
        print("Reusing held metadata. Specify reset = True to reload metadata.")

    df = meta_data

    if barcodes != None :
        df = df[df['primer_name'].map(lambda x : x in barcodes)]
        if df.empty :
            raise ValueError("None of specified barcodes found in the metadata")

    if thresholds == 'neither' :
        df = df[(df['barcode_score'] < df['threshold_barcode']) & (df['strand_switch_score'] < df['threshold_strand_switch'])]
    elif thresholds == 'barcode' :
        df = df[(df['barcode_score'] >= df['threshold_barcode']) & (df['strand_switch_score'] < df['threshold_strand_switch'])]
    elif thresholds == 'strand_switch' :
        df = df[(df['barcode_score'] < df['threshold_barcode']) & (df['strand_switch_score'] >= df['threshold_strand_switch'])]
    elif thresholds == 'both' :
        df = df[(df['barcode_score'] >= df['threshold_barcode']) & (df['strand_switch_score'] >= df['threshold_strand_switch'])]

    if funcs != None :
        for func in funcs :
            df = df[func(df)]

    if df.empty :
        raise ValueError("No data fits the specified requirements.")


    if num > len(df.index) :
        num = len(df.index)
        print('num is greater than the number of available reads. Capping at ', num)
    df = df.sample(num, axis=0)
    for read in df.index :
        alignments = [eval(df.at[read, 'alignment']), eval(df.at[read, 'strand_switch_primer_align'])]
        print("Barcode and strand_switch_primer alignment for ", read)
        align.show_alignment_simple(alignments)
    return meta_data

def barcode_variation(path_out, samples) :
    
    counts = pd.read_csv(str(path_out + '/counts.csv'), index_col = ['barcode'])
    counts_by_sample = pd.DataFrame(index=counts.columns)
    variation = pd.DataFrame(index=counts.columns)
    
    for sample in samples :
        variation.insert(0, sample, counts.loc[samples[sample]].std(axis=0) / counts.loc[samples[sample]].mean(axis=0))
        counts_by_sample.insert(0, sample, counts.loc[samples[sample]].sum(axis=0))
        plt.scatter(counts_by_sample[sample], variation[sample])
    return

    """
    data structures:
    
    debarcoded _.csv :
        csv w/ header row and a row for each read.
        cols: read_id (str), alignment (str representation of dict, see below), direction (str), primer_name (str), score (int), strand_switch_primer_align (str representation of dict, see below), biological_sequence (str), barcode_score (int), strand_switch_score (int), read_length (int), tags (str), threshold_barcode (int), threshold_strand_switch (int), minimap_alignment (if counts() has been run: None or str representation of list, see below. Else, doesn't exist)
        
        alignment dicts:
            {'aligned_query_sequence': str, 
            'aligned_target_sequence': str, 
            'cigar': str, 
            'optimal_alignment_score': int, 
            'query_begin': int, 
            'query_end': int, 
            'query_sequence': str, 
            'suboptimal_alignment_score': int, 
            'target_begin': int, 
            'target_end_optimal': int, 
            'target_end_suboptimal': int, 
            'target_sequence': str}
        minimap_alignment list:
            None or [reference_name, reference_start, cigarstring, reference_end]
    
    """
    """
    Changes:
    added sample aligment viewer view_alignments()
    changed debarcode() to allow for fastqs without tags
    Added a file check for fastq files to debarcode()
    Added ability for minimap2() to continue from a previous run by setting finish=True
    Dropped object architecture: set all functions to take path_out instead. Enabled qc functions to return and take meta_data in order to not have to re-import the meta data every time.
    Reworked count_reads() to use a better data structure and to use RPMK
    Reworked other functions to use the updated count structure
    """