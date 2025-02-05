from pathlib import Path
import math
import subprocess
import concurrent
import concurrent.futures
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq
from pyarrow import dataset
from pyarrow import csv
import matplotlib.pyplot as plt
import shutil
from scipy import stats
import numpy as np
from NanoporeAnalysis import local_io

def generate_priors(path_csv, samples, funcs_to_check = None, fields_to_check = None, any_or_all = 'all', field_to_sort_by = None, N = None, min_avg_count = 5) :
    """
    Uses compare_counts() to generate a list of genes to use as priors for another compare_counts(), meant to make the FDR more useful. This is able to be customized to restrict based on any field in the table output by compare_counts().
    
    Args :
        path_csv (str) : path to the counts csv output by count_reads(). Must have a column per gene plus columns for barcode_ID, gene, and run_label. Each row denotes a different barcode/sample/run.
        samples (list) : a list of 2 strings ['sample_1', 'sample_2'] denoting the samples to be compared. The fold changes will be defined as the first sample divided by the second.
        funcs_to_check (list) : a list of functions that accept the corresponding value in the table under fields_to_check and outputs a boolean. ie lambda x : if x >= 2. Defaults to None, in which case the function will include all N genes (see below for details on N).
        fields_to_check (list) : the fields in the compare_counts() table to be used to select genes. Should be made in parallel with funcs_to_check and should be equal length.
        any_or_all (str) : whether all funcs in funcs_to_check need to be True in order to add the gene as a prior or if a single one will suffice. Must be either 'any' or 'all'.
        field_to_sort_by (str) : the field in the compare_counts() table to be used to sort the reads prior to selection. Only useful if using the N argument to select the first N genes. Defaults to None, in which case the order in the counts csv is preserved, which is sorted by ascending pvalue when originally created by count_mapped_reads().
        N (int) : the number of genes to select. These will be selected from the gene list after sorting by field_to_sort_by. ie if N = 100 and field_to_sort_by = 'padj', this will output the 100 genes with the lowest padj. Defaults to None, in which case all genes will be evaluated, if func_to_check is set.
       
    Returns :
        priors (list) : list of gene IDs to be used as priors for compare_counts().
    comparison_table = compare_counts(path_csv, samples, make_plots = False)
    """
    comparison_table = compare_counts(path_csv, samples, make_plots = False, print_results = False, min_avg_count = min_avg_count)
    if N == None or N > comparison_table.num_rows :
        N = comparison_table.num_rows
    if field_to_sort_by != None :
        comparison_table = comparison_table.sort_by(field_to_sort_by)
    rows = comparison_table.to_pylist()
    priors = {}
    for row in rows :
        if funcs_to_check != None :
            gene_bool = True
            for func, field in zip( funcs_to_check, fields_to_check )  :
                if func( row[field] ) :
                    if any_or_all == 'any' :
                        gene_bool = True
                        break
                else :
                    gene_bool = False
                    if any_or_all == 'all' :
                        break
            if gene_bool == True :
                priors[ row['gene'] ] = row['log_2_fold_change']
        else :
            priors[ row['gene'] ] = row['log_2_fold_change']
        if len(priors) == N :
            break
    print("Found ", len(priors), " genes to use as priors.")
    return priors

def compare_counts(path_csv, samples, priors = None, neg_priors = None, gene_names = None, make_plots = True, print_results = True, sort_table_by = None, display_table_length = 100, min_avg_count = 10, use_padj = True, save_fig = False) :
    """
    Compares the gene counts in the counts csv between two samples. Converts counts to a fraction of the total number of reads, converts this to logit, uses walsch's t test with unequal variance to determine pvalue, finds the log_2 of the fold change between them, performs p adjustment, and outputs a volcano plot (if set).
    
    Args :
        path_csv (str) : path to the counts csv output by count_reads(). Must have a column per gene plus columns for barcode_ID, gene, and run_label. Each row denotes a different barcode/sample/run.
        samples (list) : a list of strings ['sample_1', 'sample_2'] denoting the samples to be compared. The fold changes will be defined as the first sample(s) divided by the second.
        priors (list) : list of genes to restrict the analysis to. This can be the output from generate_priors(). Defaults to None, in which case all genes are considered.
        make_plots (bool) : whether or not to make the volcano plot. Uses log_2_fold_change for the x axis and log_10_pvalue for the y axis and uses the colors column in the comparison_table to define the colors of the data points, which denote points with padj < 0.05. Defaults to True.
        
    Returns :
        comparison_table (pyarrow table) : a table with each gene as a row and columns for gene, log_2_fold_change, pvalue, log_10_pvalue, padj, and color.
    """
    counts = csv.read_csv(path_csv, read_options = csv.ReadOptions(block_size = 10000000)).to_pydict()
    gene_IDs_proper = []
    gene_IDs_degenerate = []
    with open("/fs/ess/PAS2506/Users/Vlad/Ref_seqs/GRCh38_latest_genomic_v1.bed", 'r') as bed_handle :
        for line in bed_handle.readlines() :
            bed_split = line.split('\t')
            if bed_split[0][0:2] == 'NC' :
                gene_IDs_proper.append(bed_split[3])
            else :
                gene_IDs_degenerate.append(bed_split[3])
    for gene in gene_IDs_degenerate :
        if gene in counts :
            if '-' in gene[5:] :
                last_dash = gene.rfind('-')
                try :
                    end_num = int(gene[last_dash+1:])
                    gene_base = gene[:last_dash]
                    if gene_base in gene_IDs_proper and gene_base in counts :
                        degen_values = counts[gene]
                        proper_values = counts[gene_base]
                        new_values = [ x+y for x, y in zip(degen_values, proper_values) ]
                        del counts[gene]
                        counts[gene_base] = new_values
                except :
                    pass
    counts = pa.table(counts)
    comparison_dict = {'gene' : [], 'log_2_fold_change' : [], 'pvalue' : [], 'neg_log_10_pvalue' : [], 'sample_1_values' : [], 'sample_2_values' : [], 'sample_1_values_norm' : [], 'sample_2_values_norm' : [], 'log_bounds' : []}
    counts_sample_1 = counts.filter( pc.field('sample') == samples[0] ).drop_columns(['barcode_ID', 'sample', 'run_label'])
    counts_sample_2 = counts.filter( pc.field('sample') == samples[1] ).drop_columns(['barcode_ID', 'sample', 'run_label'])
    count_totals_sample_1 = []
    count_totals_sample_2 = []
    for row in counts_sample_1.to_pylist() :
        count_totals_sample_1.append(sum(list(row.values())))
    for row in counts_sample_2.to_pylist() :
        count_totals_sample_2.append(sum(list(row.values())))
    barcodes_sample_1 = counts.filter( pc.field('sample') == samples[0] ).column('barcode_ID').to_pylist()
    barcodes_sample_2 = counts.filter( pc.field('sample') == samples[1] ).column('barcode_ID').to_pylist()
    print(count_totals_sample_1, count_totals_sample_2)
    print(barcodes_sample_1, barcodes_sample_2)
    if gene_names != None :
        columns = [ x for x in gene_names if x in counts.column_names ]
        counts_sample_1 = counts_sample_1.select(columns)
        counts_sample_2 = counts_sample_2.select(columns)
    counts_sample_1 = counts_sample_1.to_pydict()
    counts_sample_2 = counts_sample_2.to_pydict()
    genes = list(counts_sample_1.keys())
    for gene in genes :
        sample_1_values = counts_sample_1[gene]
        sample_2_values = counts_sample_2[gene]
        if np.mean(sample_1_values) >= min_avg_count and np.mean(sample_2_values) >= min_avg_count :
            sample_1_normalized = [ max(sample_1_value, 0.1) / count_totals_sample_1 for sample_1_value, count_totals_sample_1 in zip(sample_1_values, count_totals_sample_1) ]
            sample_2_normalized = [ max(sample_2_value, 0.1) / count_totals_sample_2 for sample_2_value, count_totals_sample_2 in zip(sample_2_values, count_totals_sample_2) ]
            sample_1_logit = [ np.log( value ) - np.log( 1 - value ) for value in sample_1_normalized ]
            sample_2_logit = [ np.log( value ) - np.log( 1 - value ) for value in sample_2_normalized ]
            t_test = stats.ttest_ind(sample_1_logit, sample_2_logit, equal_var = False)
            mean_sample_1 = np.mean(sample_1_normalized)
            mean_sample_2 = np.mean(sample_2_normalized)
            diff_ratio = mean_sample_2 / mean_sample_1
            log_2_fold_change = np.log2( diff_ratio )
            sample_1_stdev = np.std( sample_1_normalized, ddof = len(sample_1_normalized) - 1 )
            sample_2_stdev = np.std( sample_2_normalized, ddof = len(sample_2_normalized) - 1 )
            diff_error = diff_ratio * np.sqrt( ((sample_1_stdev/mean_sample_1)**2) + ((sample_2_stdev/mean_sample_2)**2) )
            if diff_error < diff_ratio :
                lower_log_bound = np.log2(diff_ratio - diff_error)
            else :
                lower_log_bound = None
            upper_log_bound = np.log2(diff_ratio + diff_error)
            comparison_dict['gene'].append(gene)
            comparison_dict['log_2_fold_change'].append(log_2_fold_change)
            comparison_dict['pvalue'].append(t_test.pvalue)
            comparison_dict['neg_log_10_pvalue'].append(-np.log10(t_test.pvalue))
            comparison_dict['sample_1_values'].append( sample_1_values )
            comparison_dict['sample_2_values'].append( sample_2_values )
            comparison_dict['sample_1_values_norm'].append( sample_1_normalized )
            comparison_dict['sample_2_values_norm'].append( sample_2_normalized )
            comparison_dict['log_bounds'].append( [lower_log_bound, upper_log_bound] )
    comparison_table = pa.table(comparison_dict).sort_by('pvalue')
    padjs = []
    colors = []
    sizes = []
    p_values = comparison_table.column('pvalue').to_pylist()
    log_2_fold_change = comparison_table.column('log_2_fold_change').to_pylist()
    significant = True
    N = comparison_table.num_rows if priors == None else len(priors)
    rank = 1
    for row in comparison_table.to_pylist() :
        skip = False
        if priors == None :
            if neg_priors != None :
                if row['gene'] in neg_priors :
                    if np.sign(row['log_2_fold_change']) == np.sign(neg_priors[row['gene']]) :
                        padjs.append(1)
                        colors.append('black')
                        sizes.append(1)
                        skip = True
            if not skip :
                padj = (row['pvalue'] * N) / rank
                padjs.append( padj )
                if ( padj <= 0.05 and significant == True ) or use_padj == False :
                    if row['log_2_fold_change'] > 0 :
                        colors.append('red')
                    elif row['log_2_fold_change'] < 0 :
                        colors.append('blue')
                    sizes.append(20)
                else :
                    colors.append('grey')
                    sizes.append(1)
                    significant = False
                rank += 1
        elif priors != None :
            if neg_priors != None :
                if row['gene'] in neg_priors :
                    if np.sign(row['log_2_fold_change']) == np.sign(neg_priors[row['gene']]) :
                        padjs.append(1)
                        colors.append('black')
                        sizes.append(1)
                        skip = True
            if not skip :
                if row['gene'] in priors :
                    padjs.append( (row['pvalue'] * N) / rank )
                    if priors[row['gene']] > 0 :
                        colors.append( 'red' )
                    elif priors[row['gene']] < 0 :
                        colors.append( 'blue' )
                    sizes.append(20)
                    rank += 1
                else :
                    padjs.append(1)
                    colors.append('grey')
                    sizes.append(1)
    comparison_table = comparison_table.append_column('padj', [padjs])
    comparison_table = comparison_table.append_column('color', [colors])
    comparison_table = comparison_table.append_column('size', [sizes])
    if make_plots :
        fig, ax = plt.subplots(1)
        ax.scatter(comparison_table.column('log_2_fold_change').to_pylist(), comparison_table.column('neg_log_10_pvalue').to_pylist(), c = colors, s = sizes)
        ax.set_xlim(-3, 3)
        ax.set_ylim(0,4)
        ax.set_xlabel(r"log$_2$( fold change )", fontsize=15)
        ax.set_ylabel(r"-log$_{10}$( p-value )", fontsize=15)
        plt.show()
        if save_fig :
            fig.savefig(save_fig, bbox_inches='tight')
    if sort_table_by != None :
        comparison_table = comparison_table.sort_by(sort_table_by)
    if priors != None and use_padj == True :
        comparison_table = comparison_table.filter( ~ (pc.field('color').isin(['grey', 'black'])) )
    if print_results :
        print( "%25s%15s%15s%15s%7s%20s%20s%15s" % ('gene', 'log_2_fold_change', 'pvalue', 'padj', 'color', 'sample_1_values', 'sample_2_values', 'log_bounds') )
        for row in comparison_table.to_pylist()[:display_table_length] :
            print('%(gene)25s%(log_2_fold_change)15.10f%(pvalue)15.10f%(padj)15.10f%(color)7s %(sample_1_values)20s %(sample_2_values)20s %(log_bounds)15s' % row)
    return comparison_table

def qc_metrics(path_dataset) :
    ds = pa.dataset.dataset([pa.dataset.dataset(x) for x in Path(path_dataset).iterdir()])
    qc = {}
    hist_max = 5000
    fig1, qc['alignment_scores_hist'] = plt.subplots(3,1, sharex=True, sharey=True)
    fig2, qc['barcode_bar'] = plt.subplots()
    ds_pass_barcode_pass_SSP = ds.filter( (~(pc.field('barcode_ID').isin(pa.array(['none matched', 'multiple', None], type=pa.string())))) & (~(pc.field('barcode_distal_SSP_edit_distance').isin(pa.array([None, -1], type=pa.float64())))) )
    ds_pass_barcode_fail_SSP = ds.filter( (~(pc.field('barcode_ID').isin(pa.array(['none matched', 'multiple', None], type=pa.string())))) & ((pc.field('barcode_distal_SSP_edit_distance').isin(pa.array([None, -1], type=pa.float64())))) )
    ds_fail_barcode = ds.filter( (pc.field('barcode_ID').isin(pa.array(['none matched', 'multiple', None], type=pa.string()))) )
    table_pass_barcode_pass_SSP = ds_pass_barcode_pass_SSP.to_table(columns = ['seq_len'])
    table_pass_barcode_fail_SSP = ds_pass_barcode_fail_SSP.to_table(columns = ['seq_len'])
    table_fail_barcode = ds_fail_barcode.to_table(columns = ['seq_len'])
    qc['alignment_scores_hist'][0].hist(table_pass_barcode_pass_SSP.column('seq_len').to_pylist(), bins=50, range=(0, hist_max))
    qc['alignment_scores_hist'][0].set_title('+ Barcode + SSP: ' + str(table_pass_barcode_pass_SSP.num_rows), {'fontsize': 8})
    qc['alignment_scores_hist'][1].hist(table_pass_barcode_fail_SSP.column('seq_len').to_pylist(), bins=50, range=(0, hist_max))
    qc['alignment_scores_hist'][1].set_title('+ Barcode - SSP: ' + str(table_pass_barcode_fail_SSP.num_rows), {'fontsize': 8})
    qc['alignment_scores_hist'][2].hist(table_fail_barcode.column('seq_len').to_pylist(), bins=50, range=(0, hist_max))
    qc['alignment_scores_hist'][2].set_title('- Barcode: ' + str(table_fail_barcode.num_rows), {'fontsize': 8})
    barcode_counts = []
    barcode_counts_strings = []
    barcode_counts_passed_SSP = ds_pass_barcode_pass_SSP.to_table(columns=['barcode_ID']).column('barcode_ID').value_counts().to_pylist()
    barcode_counts_failed_SSP = ds_pass_barcode_fail_SSP.to_table(columns=['barcode_ID']).column('barcode_ID').value_counts().to_pylist()
    # barcode_counts_passed_SSP.sort( key = lambda x : x['values'] )
    # barcode_counts_failed_SSP.sort( key = lambda x : x['values'] )
    for passed, failed in zip(barcode_counts_passed_SSP, barcode_counts_failed_SSP) :
        barcode_counts.append(passed['counts'])
        barcode_counts.append(failed['counts'])
        barcode_counts_strings.append(passed['values'] + ' ++: ' + str('%15i' % passed['counts']))
        barcode_counts_strings.append(failed['values'] + ' +-: ' + str('%15i' % failed['counts']))
    qc['barcode_bar'].barh(barcode_counts_strings, barcode_counts)
    return