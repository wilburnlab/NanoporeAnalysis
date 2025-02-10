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
        min_avg_count (any number not null or inf) : only genes with a total average counts value above min_avg_counts will be included. 
       
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

def compare_counts(path_csv, samples, priors = None, neg_priors = None, make_plots = True, print_results = True, sort_table_by = None, display_table_length = 100, min_avg_count = 10, use_padj = True, save_fig = False) :
    """
    Compares the gene counts in the counts csv between two samples. Converts counts to a fraction of the total number of reads per barcode, converts this to logit space,
        uses walsch's t test with unequal variance to determine pvalue, finds the log_2 of the fold change between them, performs benjamini-hochberg pvalue adjustment,
        and outputs a volcano plot (if set) and prints the resulting comparison table (if set). This comparison can be restricted to certain genes through the priors,
        neg_priors, and min_avg_count options, but this won't change the calculation for total reads per barcode. The propogated error is also calculated for the fold_changes.
        This error can be applied to the log_2_fold_change by raising 2 to the log_2_fold_change value, then adding or subtracting the propogated_error. This produces an undefined
        value when the error is greater than the fold change, which usually indicates that the fold change isn't very accurate.
    
    Args :
        path_csv (str) : path to the counts csv output by count_reads(). Must have a column per gene plus columns for barcode_ID, gene, and run_label.Each row denotes a different barcode/sample/run.
        samples (list) : a list of strings ['sample_1', 'sample_2'] denoting the samples to be compared. The fold changes will be defined as the first sample divided by the second.
        priors (list or dict) : list of genes to restrict the analysis to, or a dict of gene:log_2_fold_change. This can be the output from generate_priors().
            If this is set, then only genes included in priors will be shown in color in plots and included in the pvalue adjustment.
            If priors is a dict (such as the output of generate_priors()), then the colors of the resulting plot will represent the values in the priors dict.
        neg_priors (list or dict) : list of genes to block from analysis. This can be the output from generate_priors().
            If this is a dict, then genes will only be blocked if the sign of the log_2_fold_change of the gene matches the sign of the value in the neg_priors dict.
            Genes in neg_priors will be shown in black on the plot if neg_priors is a list or if it's a dict and the signs match.
        make_plots (bool) : whether or not to make the volcano plot. Uses log_2_fold_change for the x axis and -log_10(pvalue) for the y axis and uses the colors
            column in the comparison_table to define the colors of the data points. Defaults to True.
        print_results (bool) : whether or not to print the table of comparison values.
        sort_table_by (list of tuples of (field, direction) pairs) : the fields to sort the results table by. Must be a list of tuples, such as [('pvalue', 'ascending')].
        display_table_length (int) : number of rows of the results table to print. Defaults to 100 so your screen doesn't get totally filled.
        min_avg_count (any number not null or inf) : only genes with a total average counts value above min_avg_counts will be included in this analysis. 
        use_padj (bool) : whether or not to color the plot based on the padj instead of pvalue.
        save_fig (str) : path to save the plot to. If None, the plot won't be saved. Recommended to save as a .svg file.
        
    Returns :
        comparison_table (pyarrow table) : a table with each gene as a row and columns for gene, log_2_fold_change, pvalue, padj, color, sample_1_values, sample_2_values, and propogated_error.
    """
    counts = csv.read_csv(path_csv, read_options = csv.ReadOptions(block_size = 10000000)).filter(pc.field('sample').isin(samples))
    comparison_dict = {'gene' : [], 'log_2_fold_change' : [], 'pvalue' : [], 'sample_1_values' : [], 'sample_2_values' : [], 'propogated_error' : []}
    sample_IDs = counts.column('sample').to_pylist()
    counts_values = counts.drop_columns(['barcode_ID', 'sample', 'run_label'])
    per_row_totals = pa.array([ pc.sum(list(x.values())) for x in counts_values.to_pylist() ])
    genes_above_min_avg = [ x for x in counts_values.column_names if pc.mean(counts_values.column(x)).as_py() >= min_avg_count ]
    counts_passed = counts_values.select(genes_above_min_avg)
    counts_passed_labeled = counts_passed.append_column('sample', [sample_IDs])
    comparison_dict['sample_1_values'] = counts_passed_labeled.filter(pc.field('sample') == samples[0]).drop_columns(['sample']).to_pydict().values()
    comparison_dict['sample_2_values'] = counts_passed_labeled.filter(pc.field('sample') == samples[1]).drop_columns(['sample']).to_pydict().values()
    comparison_dict['gene'] = counts_passed.column_names
    fill_array = pa.array([ 0.1 for x in range(counts_passed.num_rows) ])
    counts_filled = pa.Table.from_arrays( [ pc.max_element_wise(x, fill_array) for x in counts_passed.itercolumns() ], counts_passed.column_names )
    counts_normalized = pa.Table.from_arrays( [ pc.divide( x, per_row_totals ) for x in counts_filled.itercolumns() ], counts_filled.column_names )
    counts_normalized_labeled = counts_normalized.append_column('sample', [sample_IDs])
    counts_logits = pa.Table.from_arrays( [ pc.subtract( pc.ln(x), pc.ln(pc.subtract( 1, x ))) for x in counts_normalized.itercolumns() ], counts_normalized.column_names )
    counts_logits_labeled = counts_logits.append_column('sample', [sample_IDs])
    sample_1_means = [ pc.mean(x) for x in counts_normalized_labeled.filter(pc.field('sample') == samples[0]).drop_columns(['sample']).itercolumns() ]
    sample_2_means = [ pc.mean(x) for x in counts_normalized_labeled.filter(pc.field('sample') == samples[1]).drop_columns(['sample']).itercolumns() ]
    sample_1_logits = counts_logits_labeled.filter(pc.field('sample') == samples[0]).drop_columns(['sample']).to_pydict().values()
    sample_2_logits = counts_logits_labeled.filter(pc.field('sample') == samples[1]).drop_columns(['sample']).to_pydict().values()
    comparison_dict['pvalue'] = [ stats.ttest_ind(x, y, equal_var = False).pvalue for x, y in zip(sample_1_logits, sample_2_logits) ]
    fold_changes = pc.divide(sample_2_means, sample_1_means)
    comparison_dict['log_2_fold_change'] = pc.log2(fold_changes)
    sample_1_stddevs = [ pc.stddev( x, ddof = len(x) - 1 ) for x in counts_normalized_labeled.filter(pc.field('sample') == samples[0]).drop_columns(['sample']).itercolumns() ]
    sample_2_stddevs = [ pc.stddev( x, ddof = len(x) - 1 ) for x in counts_normalized_labeled.filter(pc.field('sample') == samples[1]).drop_columns(['sample']).itercolumns() ]
    comparison_dict['propogated_error'] = pc.multiply( fold_changes, pc.sqrt(pc.add( pc.power( pc.divide(sample_1_stddevs, sample_1_means), 2 ), pc.power(pc.divide(sample_2_stddevs, sample_2_means), 2 ) )))
    comparison_table = pa.table(comparison_dict).sort_by('pvalue')
    black_genes = []
    grey_genes = []
    if type(neg_priors) == list :
        black_genes += neg_priors
    elif type(neg_priors) == dict :
        for x in neg_priors :
            if x in comparison_table.column('gene').to_pylist() :
                if np.sign(neg_priors[x]) == np.sign(comparison_table.filter(pc.field('gene') == x).to_pylist()[0]['log_2_fold_change']) :
                    black_genes.append(x)
    if type(priors) == list :
        grey_genes += [ x for x in comparison_table.column('gene').to_pylist() if x not in priors + black_genes ]
    elif type(priors) == dict :
        grey_genes += [ x for x in comparison_table.column('gene').to_pylist() if x not in list(priors.keys()) + black_genes ]
    color_genes = [ x for x in comparison_table.column('gene').to_pylist() if x not in grey_genes + black_genes ]
    N = len(color_genes)
    significant = True
    rank = 1
    padjs = []
    colors = []
    sizes = []
    for row in comparison_table.to_pylist() :
        if row['gene'] in color_genes :
            padj = (row['pvalue'] * N) / rank
            padjs.append( padj )
            if type(priors) == dict :
                if priors[row['gene']] > 0 :
                    sizes.append(20)
                    colors.append('red')
                else :
                    sizes.append(20)
                    colors.append('blue')
            elif ( padj <= 0.05 and significant == True ) or use_padj == False :
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
        elif row['gene'] in grey_genes :
            padjs.append(1)
            colors.append('grey')
            sizes.append(1)
        elif row['gene'] in black_genes :
            padjs.append(1)
            colors.append('black')
            sizes.append(1)
    comparison_table = comparison_table.append_column('padj', [padjs])
    comparison_table = comparison_table.append_column('color', [colors])
    comparison_table = comparison_table.append_column('size', [sizes])
    comparison_table = comparison_table.sort_by([('size', 'descending')])
    if make_plots :
        fig, ax = plt.subplots(1)
        ax.scatter(comparison_table.filter(pc.field('color')=='grey').column('log_2_fold_change').to_pylist(),
                   pc.negate(pc.log10(comparison_table.filter(pc.field('color')=='grey').column('pvalue'))).to_pylist(),
                   c = comparison_table.filter(pc.field('color')=='grey').column('color').to_pylist(),
                   s = comparison_table.filter(pc.field('color')=='grey').column('size').to_pylist())
        ax.scatter(comparison_table.filter(pc.field('color')!='grey').column('log_2_fold_change').to_pylist(),
                   pc.negate(pc.log10(comparison_table.filter(pc.field('color')!='grey').column('pvalue'))).to_pylist(),
                   c = comparison_table.filter(pc.field('color')!='grey').column('color').to_pylist(),
                   s = comparison_table.filter(pc.field('color')!='grey').column('size').to_pylist())
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
        print( "%25s%15s%15s%15s%7s%20s%20s%15s" % ('gene', 'log_2_fold_change', 'pvalue', 'padj', 'color', 'sample_1_values', 'sample_2_values', 'propogated_error') )
        for row in comparison_table.to_pylist()[:display_table_length] :
            print('%(gene)25s%(log_2_fold_change)15.10f%(pvalue)15.10f%(padj)15.10f%(color)7s %(sample_1_values)20s %(sample_2_values)20s %(propogated_error)15s' % row)
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