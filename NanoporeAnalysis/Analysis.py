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
from numbers import Number
from pod5.tools import pod5_view
from pod5.tools import pod5_subset
from importlib import reload
reload(align)

"""
condense/cleanup code
progress reporter for alignment
Give Dorado options. ie Simplex calling
Move more path declarations into __init__
Skip over fastq conversion by using bams and pysam

QC metrics:
    plot Read lengths by barcode or sample.
    plot Biological sequence lengths by barcode or sample.
    Compare between barcodes of a sample
    View sample alignments for any relevant steps.
    
Analysis functions:
    Diff exp between samples/barcodes for all features referenced or for specified features.
    Diff exp of mRNA variants
    PCA    
    
Check alignments from the unclassified file. See if any can be properly barcoded.
Add method for chopping up input data into smaller files
make  get_best_alignment return something special for when the ss primer filter isnt met

"""
class ONTanalysis :
    def __init__(self, path_data, path_out) :
        self.path_data = str(Path(path_data).absolute())
        if Path(self.path_data).is_dir() == False :
            print("Error: Data path is invalid or doesn't exist. Please point to an existing directory containing Pod5 files.")
        self.path_out = str(Path(path_out).absolute())
        try :
            Path(self.path_out).mkdir(parents=True, exist_ok=True)
        except FileExistsError :
            if Path(self.path_out + "/temp_aligned_reads/").is_dir() :
                print("Error: Ouput path may already contain data. Please point to an empty or non-existing directory")
        self.path_split_pod5s = self.path_out + '/split_pod5s'
        Path(self.path_split_pod5s).mkdir(parents=True, exist_ok=True)
        self.path_logs = self.path_out + '/logs'
        Path(self.path_logs).mkdir(parents=True, exist_ok=True)
        self.path_pod5_bams = self.path_out + '/pod5_bams'
        Path(self.path_pod5_bams).mkdir(parents=True, exist_ok=True)
        self.path_tmp_aligned = self.path_out + "/temp_aligned_reads/"
        Path(self.path_tmp_aligned).mkdir(parents=True, exist_ok=True)
        self.path_debarcoded = self.path_out + '/debarcoded'
        Path(self.path_debarcoded).mkdir(parents=True, exist_ok=True)
        self.path_sams = self.path_out + '/sams'
        Path(self.path_sams).mkdir(parents=True, exist_ok=True)
        self.path_bams = self.path_out + '/bams'
        Path(self.path_bams).mkdir(parents=True, exist_ok=True)
        self.path_fastqs = self.path_out + '/fastq'
        Path(self.path_fastqs).mkdir(parents=True, exist_ok=True)
        self.path_fastas = self.path_out + '/fasta'
        Path(self.path_fastas).mkdir(parents=True, exist_ok=True)
        self.path_basecalls = self.path_out + '/basecalls'
        Path(self.path_fastqs).mkdir(parents=True, exist_ok=True)
        self.path_qc = self.path_out + '/qc.'
    
    def run_dorado_slurm(self, path_dorado, path_model, account, mail, workers=1) :
        """
        
        
        """
        pod5_view.view_pod5([Path(self.path_data)], Path(self.path_out), include = "read_id, channel", force_overwrite=True)
        pod5_subset.subset_pod5([Path(self.path_data)], self.path_split_pod5s, columns = ["channel"], table = Path(self.path_out + "/view.txt"))
        files = [x for x in Path(self.path_split_pod5s).iterdir() if x.is_file()]
        arrays = np.array_split(np.array(files), workers)
        for i in range(workers) :
            shutil.move(arrays[i], str(self.path_out + '/split_pod5s_byworker' + str(i)))
            script = [
                '#!/bin/bash\n',
                str('#SBATCH --account=' + account + '\n'),
                str('#SBATCH --job-name=dorado_' + str(i) + '\n'),
                '#SBATCH --nodes=1 --ntasks-per-node=1 --gpus-per-node=2\n',
                str('#SBATCH --output=' + self.path_logs + 'dorado_' + str(i) + '.out' + '\n'),
                '#SBATCH --mail-type=FAIL\n',
                str('#SBATCH --mail-user=' + mail + '\n'),
                str('cd ' + self.path_basecalls + '\n'),
                str(path_dorado + " duplex " + path_model + self.path_split_pod5s + '/' + i + ' > ' + self.path_pod5_bams + i + '.bam'),
            ]
            with open(str(self.path_logs + 'dorado_' + str(i) + '.sh'), 'w') as handle:
                handle.writelines(script)
            args = ['sbatch', str(self.path_logs + 'dorado_' + str(i) + '.sh')]
            subprocess.run(args)
        return
        
    def run_dorado(self, path_dorado, path_model, workers=1, skip_split=False) :
        if skip_split == False :
            pod5_view.view_pod5([Path(self.path_data)], Path(self.path_out), include = "read_id, channel", force_overwrite=True)
            pod5_subset.subset_pod5([Path(self.path_data)], Path(self.path_split_pod5s), columns = ["channel"], table = Path(self.path_out + "/view.txt"))
        files = [x for x in Path(self.path_split_pod5s).iterdir() if x.is_file()]
        arrays = np.array_split(np.array(files), workers)
        for i in range(workers) :
            Path(self.path_out + '/split_pod5s_byworker/' + str(i)).mkdir(parents=True, exist_ok=True)
            for path in arrays[i] :
                shutil.copy(str(path), str(self.path_out + '/split_pod5s_byworker/' + str(i)))
            args = [str(path_dorado + " duplex " + path_model + ' ' + str(self.path_out + '/split_pod5s_byworker/' + str(i)) + ' > ' + self.path_pod5_bams + '/' + str(i) + '.bam')]
            subprocess.run(args, shell=True)
        return
    
    def bam_to_fastq(self) :
        files = [x for x in Path(self.path_pod5_bams).iterdir() if x.is_file()]
        Path(self.path_pod5_bams + '/sorted_by_name/').mkdir(parents=True, exist_ok=True)
        
        for file in files : 
            print(file)
            pysam.sort('-o', str(self.path_pod5_bams + '/sorted_by_name/' + file.stem + '_sorted.bam'), '-n', str(file) )
            Path(str(self.path_fastqs + '/'+ file.stem + '.fq')).touch()
            pysam.fastq('-0', str(self.path_fastqs + '/'+ file.stem + '.fq'), '-T', 'dx', str(self.path_pod5_bams + '/sorted_by_name/' + file.stem + '_sorted.bam'))
            
        return
                
    def debarcode(self, path_barcodes, strand_switch_primer, filter_barcode_score = 0, filter_barcode_distance = 0, filter_strand_switch_score = None, threshold_barcode = None, threshold_strand_switch = None) :
        """
        aligns barcodes and strand_switch primer, then saves to file by barcode and whether the fliter was met.

        args:
            path_barcodes (str) : file path for a csv containing the barcodes
            strand_switch_primer (str) : sequence for the strand_switch_primer
            filter_barcode_score (int) : filter for barcode alignment. slightly speeds up alignment by 10-20% if set to ~90% of the max alignment score
            filter_barcode_distance (int) : similar to filter_barcode_score, but for min aligned distance. Less effective
            filter_strand_switch_score : minimum score for the strand_switch_primer alignment. Below this filter, a full alignment is done for the barcode
            threshold_barcode : minimum score for barcode alignment when sorting into files. Reads not meeting this filter get put into 'unclassified.fastq'
            threshold_strand_switch : minimum score for ss primer alignment when sorting into files. Reads not meeting this filter get put into a '*_incomplete.fastq' file for their given barcode
        Output:
            None: Writes to file.
        """
        if filter_strand_switch_score == None : filter_strand_switch_score = len(strand_switch_primer) * 0.5 * 2
        if threshold_strand_switch == None : threshold_strand_switch = len(strand_switch_primer) * 0.5 * 2
        
        reload(align)

        #load barcodes
        barcodes = pd.read_csv(path_barcodes, names=['name', 'seq'])
        barcodes_clean = barcodes.copy()
        barcodes_clean['seq'] = barcodes['seq'].str.replace('[^ATCGatcg]', '', regex=True).apply(utils.reverse_complement)
        barcodes_clean = barcodes_clean.set_index('name').T.to_dict('records')[0]
        
        #find data
        files = [x for x in Path(self.path_fastqs).iterdir() if x.is_file()][:1]

        start_time = time.time()
        aligned_reads = []
        dump_index = 1
        read_count = 0

        for file in files:

            file_name = str(Path(file).stem)

            print('Debarcoding file ' + file_name + ';   File ', files.index(file) + 1, 'out of ', len(files))
            data = local_io.read_fastx(file)
            
            executor = concurrent.futures.ProcessPoolExecutor()
            futures = [ executor.submit( align.get_best_barcode, data[read_id]['Sequence'], barcodes_clean, strand_switch_primer, filter_strand_switch_score, read_id, data[read_id]['Tags'], filter_barcode_score, filter_barcode_distance) for read_id in data if not 'dx:i:-1' in data[read_id]['Tags']]
            concurrent.futures.wait( futures )
            for f in futures: aligned_reads.append(f.result())
            
            if sys.getsizeof(aligned_reads) > 100000 or files.index(file) + 1 >= len(files) :
                
                print('Saved data to temp file ' + str(dump_index))
                with open(self.path_tmp_aligned +'saved_read_alignments_' + str(dump_index) + '.txt', 'wb') as handle:
                    pickle.dump(aligned_reads, handle)
                dump_index += 1
                read_count += len(aligned_reads)
                aligned_reads = []
        
        print("Debarcoded ", read_count, " reads in: ", time.time()-start_time, " sec")
        
        print("Saving debarcoded reads.")
        
        start_time = time.time()
        
        # Initialize data structure for storing results
        seq_by_barcode = {'unclassified' : {}}
        data_by_barcode = {'unclassified' : {}}
        
        for primer_name in barcodes_clean.keys():
            seq_by_barcode[primer_name] = {}
            seq_by_barcode[primer_name + '_incomplete'] = {}
            data_by_barcode[primer_name] = {}
            data_by_barcode[primer_name + '_incomplete'] = {}
        
        # Classify each sequence based on how well the barcode and opposite primer are aligned
        files_aligned = [str(x) for x in Path(self.path_tmp_aligned).iterdir() if x.is_file()]
        
        for file in files_aligned :
            with open(file, 'rb') as handle:
                aligned_reads = pickle.load(handle)
                
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
                
        start_time = time.time()
        
        # Save data into multiple FASTA files

        for category in seq_by_barcode.keys():
            fasta_file = Path(self.path_fastas) / (category + '.fa')
            local_io.write_fasta(fasta_file, seq_by_barcode[category], mode='w')
            df = pd.DataFrame.from_dict(data_by_barcode[category], orient = 'index')
            df.to_csv(str(self.path_debarcoded + '/' + category + '.csv'), header=True, index_label='read_id')
        
        for item in Path(self.path_tmp_aligned).glob('*') :
            if item.is_file():
                item.unlink()
            elif item.is_dir():
                shutil.rmtree(item)
        
        print("Saved to disk in: ", time.time()-start_time, " sec")
        return
    
    def SLURM_Minimap2(self, path_minimap2, path_ref, account, mail) :
        """
        Calls Minimap2 to map all reads in given path based on a given reference. Saves mapped reads into new folder.

        Args:
            path_ref (str) = path to reference file
            path_minimap2 (str) = path to minimap2 executable
            account (str) = OSC account
            mail (str) = email account to send updates to

        Output: None
        """

        #Minimap2 alignment via SLURM
        
        files_debarcoded = [str(x) for x in Path(self.path_debarcoded).iterdir() if x.is_file() ] #and x.name != 'unclassified.fa']
        for file in files_debarcoded :
            script = [
                '#!/bin/bash\n',
                str('#SBATCH --account=' + account + '\n'),
                str('#SBATCH --job-name=minimap' + Path(file).stem + '\n'),
                '#SBATCH --nodes=1\n',
                str('#SBATCH --output=' + str(path_tmp_mapped + Path(file).stem + '.out') + '\n'),
                '#SBATCH --mail-type=FAIL\n',
                '#SBATCH --mem=32gb\n',
                str('#SBATCH --mail-user=' + mail + '\n'),
                str(path_minimap2 + ' -ax splice ' + path_ref + ' ' + file + ' > ' + self.path_sams + '/' + Path(file).stem + '.sam'),
            ]
            with open(str(path_tmp_mapped + Path(file).stem + '.sh'), 'w') as handle:
                handle.writelines(script)
                
            args = ['sbatch', str(path_tmp_mapped + Path(file).stem + '.sh')]
            
            subprocess.run(args)
        return
    
    def minimap2(self, path_minimap2, path_ref) :
        """
        Calls Minimap2 to map all reads in given path based on a given reference. Saves mapped reads into new folder.

        Args:
            path_ref (str) = path to reference file
            path_minimap2 (str) = path to minimap2 executable

        Output: None
        """

        #Minimap2 alignment
        files_debarcoded = [str(x) for x in Path(self.path_fastas).iterdir() if x.is_file() and x.suffix == '.fa' ]
        for file in files_debarcoded :
            print('mapping file ', Path(file).stem)
            args = [str(path_minimap2 + ' -ax splice ' + path_ref + ' ' + file + ' > ' + self.path_sams + '/' + Path(file).stem + '.sam')]
            subprocess.run(args, shell=True)
            pysam.sort("-o", str(self.path_bams + '/' + Path(file).stem + '.bam'), self.path_sams + '/' + Path(file).stem + '.sam')
            pysam.index(str(self.path_bams + '/' + Path(file).stem + '.bam'))
            
            
        print("Done mapping")
        return
    
    def count_reads(self, path_bed_file, samples, count = True) :
        
        """
        Makes read counts for mapped reads in self.path_out against a provided .bed file. Outputs to a new folder.
        
        Args:
            path_bed_file (str) = path to bed file
        Outputs: None
        """
        
        if count == True :
            bam_files = [str(x) for x in Path(self.path_bams).iterdir() if x.is_file() and x.suffix == '.bam']
            start_time = time.time()
            counts = pd.DataFrame()
            for file in bam_files:
                print("counting " + str(Path(file).stem))
                bam_file = pysam.AlignmentFile(file, 'rb')
                sample_name = Path(file).stem
                counts_by_feature = []
                with open(path_bed_file, 'r') as bed_in:
                    bed_lines = bed_in.readlines()
                for line in bed_lines :
                    cols = line.split('\t')
                    bam_iter = bam_file.fetch(cols[0], int(cols[1]), int(cols[2]))

                    read_count = 0
                    read_seq = []
                    for x in bam_iter:
                        name = x.query_name
                        seq = x.query_sequence
                        #if not isinstance(seq, str): seq = "EMPTY"
                        #if isinstance(seq, str):
                            #read_seq.append(' '.join([name, seq]))
                        read_count += 1
                    #read_seq_str = '\n'.join(read_seq)
                    feature_name = cols[3]
                    counts_by_feature.append([feature_name, read_count])
                if counts.empty :
                    counts = pd.DataFrame(data = counts_by_feature, columns = ['feature_name', sample_name]).set_index('feature_name')
                else :
                    counts = pd.concat([counts, pd.DataFrame(data = counts_by_feature, columns = ['feature_name', sample_name]).set_index('feature_name')], axis=1)
            counts.to_csv(str(self.path_out + '/counts.csv'), header=True, index_label='feature_name')
                
            print('finished counting in ' + str( time.time()-start_time ) + 'sec')
        else:
            counts = pd.read_csv(str(self.path_out + '/counts.txt')).set_index('feature_name')
        
        counts_by_sample = pd.DataFrame()
        
        for sample in samples :
            
            for barcode in samples[sample] :
                if counts_by_sample.empty :
                    counts_by_sample = pd.DataFrame(data=counts[barcode].values, index=counts[barcode].index, columns=[sample])
                else :
                    counts_by_sample = counts_by_sample.combine(pd.DataFrame(data=counts[barcode].values, index=counts[barcode].index, columns=[sample]), lambda x,y : x+y, fill_value=0, overwrite=True)
        
        counts_by_sample.to_csv(str(self.path_out + '/counts_by_sample.csv'), header=True, index_label='feature_name')
        
        print('done')
        
        return
    
    def count_reads_util(self, bed_line, bam_file_path, sample_name):
        """Defunct. Return the # of reads from bam file that overlap an interval given by bed file line."""
        
        cols = bed_line.split('\t')
        
        bamfile = pysam.AlignmentFile(bam_file_path, "rb")
        bam_iter = bamfile.fetch(cols[0], int(cols[1]), int(cols[2]))
        
        read_count = 0
        read_seq = []
        for x in bam_iter:
            name = x.query_name
            seq = x.query_sequence
            #if not isinstance(seq, str): seq = "EMPTY"
            #if isinstance(seq, str):
                #read_seq.append(' '.join([name, seq]))
            read_count += 1
        #read_seq_str = '\n'.join(read_seq)
        feature_name = cols[3]
        
        
        return [sample_name, feature_name, read_count]
    
    def qc_metrics(self, alignment_scores_hist_max=None) :
        qc = {}
        files = [str(x) for x in Path(self.path_debarcoded).iterdir() if x.is_file()]
        df = pd.DataFrame()
        for file in files :
            df = pd.concat([df, pd.read_csv(file).set_index('read_id')], axis = 0)
        qc['read_count'] = len(df.index)
        #print(df)
        if alignment_scores_hist_max == None : alignment_scores_hist_max = 3000
        fig1, qc['alignment_scores_hist'] = plt.subplots(2,2, sharex=True)
        qc['alignment_scores_hist'][0,0].hist(df[(df['barcode_score'] >= 60) & (df['strand_switch_score'] >= 20)]['read_length'], bins=50, range=(0, alignment_scores_hist_max))
        qc['alignment_scores_hist'][0,0].set_title('Passed Barcode and Strand Switch', {'fontsize': 10})
        qc['alignment_scores_hist'][0,1].hist(df[(df['barcode_score'] >= 60) & (df['strand_switch_score'] < 20)]["read_length"], bins=50, range=(0, alignment_scores_hist_max))
        qc['alignment_scores_hist'][0,1].set_title('Passed Barcode, failed Strand Switch', {'fontsize': 10})
        qc['alignment_scores_hist'][1,0].hist(df[(df['barcode_score'] < 60) & (df['strand_switch_score'] >= 20)]["read_length"], bins=50, range=(0, alignment_scores_hist_max))
        qc['alignment_scores_hist'][1,0].set_title('Failed Barcode, passed Strand Switch', {'fontsize': 10})
        qc['alignment_scores_hist'][1,1].hist(df[(df['barcode_score'] < 60) & (df['strand_switch_score'] < 20)]["read_length"], bins=50, range=(0, alignment_scores_hist_max))
        qc['alignment_scores_hist'][1,1].set_title('Failed Barcode and Strand Switch', {'fontsize': 10})
        
        
        qc['passed_barcode_passed_strand_switch'] = len(df[(df['barcode_score'] >= 60) & (df['strand_switch_score'] >= 20)])
        qc['passed_barcode_failed_strand_switch'] = len(df[(df['barcode_score'] >= 60) & (df['strand_switch_score'] < 20)])
        qc['failed_barcode_passed_strand_switch'] = len(df[(df['barcode_score'] < 60) & (df['strand_switch_score'] >= 20)])
        qc['failed_barcode_failed_strand_switch'] = len(df[(df['barcode_score'] < 60) & (df['strand_switch_score'] < 20)])
        
        for barcode in df['primer_name'].drop_duplicates() :
            qc[str(barcode + '_count')] = len(df[df['primer_name'] == barcode])
        return
    
    def diff_exp(self, samples, meta_data, design_factors) :
        counts_by_sample = pd.read_csv(self.path_out + '/counts_by_sample.csv', index_col='feature_name')
        
        volcano_data = counts_by_sample[counts_by_sample.min(axis=1) >= 5].loc[:,samples].T
        volcano_meta = meta_data.loc[samples, :]
        
        dds = DeseqDataSet(
            counts=volcano_data,
            metadata=volcano_meta,
            design_factors=design_factors,
            quiet=True
        )
        dds.deseq2()
        stat_res = DeseqStats(dds, quiet=True)
        stat_res.summary()
        df = stat_res.results_df.reset_index().dropna().copy()
        pd.DataFrame.from_dict(df).to_csv(str(self.path_out + '/diff_seq_' + str(samples) + '_' + design_factors + '.csv'), header=True)
        df = df.drop(df['pvalue'].idxmin(axis=0))
        
        visuz.GeneExpression.volcano(df=df, lfc='log2FoldChange', pv='pvalue', show=True)
        return