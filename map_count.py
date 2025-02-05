from pathlib import Path
import concurrent
import concurrent.futures
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq
from pyarrow import dataset as ds
from pyarrow import csv
import numpy as np
import mappy
from NanoporeAnalysis import utils

def assign_gene_ID_genomic(table, path_bed) :
    bed_dict = {
        'bed_ctg' : [],
        'bed_start' : [],
        'bed_end' : [],
        'bed_gene' : [],
        'bed_strand' : [],
    }
    with open(path_bed, 'r') as bed_handle :
        for line in bed_handle.readlines() :
            bed_split = line.split('\t')
            bed_dict['bed_ctg'].append(bed_split[0])
            bed_dict['bed_start'].append(int(bed_split[1]))
            bed_dict['bed_end'].append(int(bed_split[2]))
            bed_dict['bed_gene'].append(bed_split[3])
            bed_dict['bed_strand'].append(bed_split[5].rstrip('\n'))
    bed_table = pa.table(bed_dict)
    minimap2_gene_IDs = []
    for row in table.select(['minimap2_r_st', 'minimap2_r_en', 'minimap2_strand', 'minimap2_ctg']).to_pylist() :
        if row['minimap2_ctg'] == None :
            minimap2_gene_IDs.append('None')
        else :
            potential_matches = bed_table.filter( (pc.field('bed_ctg') == row['minimap2_ctg']) & (pc.field('bed_start') <= row['minimap2_r_en']) & (pc.field('bed_end') >= row['minimap2_r_st']) ).to_pylist()
            if len(potential_matches) == 1 :
                minimap2_gene_IDs.append(potential_matches[0]['bed_gene'])
            elif len(potential_matches) > 1 :
                matches_scored = []
                for match in potential_matches :
                    match_score = (min(int(match['bed_end']), row['minimap2_r_en']) - max(int(match['bed_start']), row['minimap2_r_st'])) / max( row['minimap2_r_en'] - row['minimap2_r_st'], int(match['bed_end']) - int(match['bed_start']) )
                    if match['bed_strand'] == row['minimap2_strand'] :
                        match_score += 1
                    matches_scored.append([match, match_score])
                minimap2_gene_IDs.append(max( matches_scored, key = lambda x : x[1] )[0]['bed_gene'])
            else :
                minimap2_gene_IDs.append('none')
    return table.append_column('minimap2_gene_ID', [minimap2_gene_IDs])

def assign_gene_ID_rna(table, path_ref) :
    ref = local_io.read_fastx(path_ref)
    contigs_to_gene_IDs = {}
    gene_IDs = []
    for key in ref.keys() :
        comma_split = key.split(', ')
        for item in comma_split :
            if 'transcript variant' not in item and '(' in item and ')' in item :
                left_para_split = item.split('(')
                right_para_split = left_para_split[ len(left_para_split) - 1 ].split(')')
                gene_ID = right_para_split[0]
        space_split = comma_split[0].split(' ')
        contig_ID = space_split[0]
        contigs_to_gene_IDs[contig_ID] = gene_ID
        gene_IDs.append(gene_ID)
    minimap2_gene_IDs = []
    for ctg in table.column(['minimap2_ctg']).to_pylist() :
        if ctg == None :
            minimap2_gene_IDs.append('None')
        else :
            minimap2_gene_IDs.append(contigs_to_gene_IDs[ctg])
    return table.append_column('minimap2_gene_ID', [minimap2_gene_IDs])

def minimap2_table(table, path_ref, preset='splice', path_bed=None, skip_duplex_parents=True) :
    """
    Runs minimap2 on the biological sequences in a pyarrow table. 
    
    Args :
        table (pyarrow table) : contains sequencing reads as rows. This should only be run after running debarcode(), as it relies on some barcode information.
        path_ref (str) : path to the sequence reference to be used.
        preset (str) : the preset to use for minimap2. Defaults to 'splice', which is designed to be splice-aware.
    
    Returns :
        table (pyarrow table) : the input table with minimap alignments appended as new columns.
    """
    if skip_duplex_parents and 'dx:i' in table.column_names :
        skip_table = table.filter( (pc.field('dx:i') == '-1') & ( pc.field('barcode_ID').isin(['none matched', 'multiple']) ) )
        valid_table = table.filter( (pc.field('dx:i') != '-1') & ( ~ pc.field('barcode_ID').isin(['none matched', 'multiple']) ) )
    else :
        skip_table = table.filter( pc.field('barcode_ID').isin(['none matched', 'multiple']) )
        valid_table = table.filter( ~ pc.field('barcode_ID').isin(['none matched', 'multiple']) )
    seqs = valid_table.column('seq').to_pylist()
    biological_seq_indices = valid_table.column('barcode_biological_seq_indices').to_pylist()
    barcode_flags = valid_table.column('barcode_flag').to_pylist()
    directions = valid_table.column('barcode_direction').to_pylist()
    aligner = mappy.Aligner(path_ref, preset=preset, best_n=1)
    alignments_dict = {
        'minimap2_q_st' : [],
        'minimap2_q_en' : [],
        'minimap2_strand' : [],
        'minimap2_ctg' : [],
        'minimap2_ctg_len' : [],
        'minimap2_r_st' : [],
        'minimap2_r_en' : [],
        'minimap2_mlen' : [],
        'minimap2_blen' : [],
        'minimap2_mapq' : [],
        'minimap2_mapq' : [],
        'minimap2_cigar' : [],
        'minimap2_trans_strand' : []
    }
    for seq, flag, indices, direction in zip(seqs, barcode_flags, biological_seq_indices, directions) :
        minimap2_gene_ID = None
        alignment = dict.fromkeys(alignments_dict.keys())
        if flag[1] == 1 :
            if direction == 'forward' :
                hits = aligner.map(seq[indices[0] : indices[1]])
            else :
                hits = aligner.map(utils.reverse_complement(seq[indices[0] : indices[1]]))
            for hit in hits :
                if hit.is_primary :
                    alignment['minimap2_q_st'] = hit.q_st
                    alignment['minimap2_q_en'] = hit.q_en
                    alignment['minimap2_strand'] = hit.strand
                    alignment['minimap2_ctg'] = hit.ctg
                    alignment['minimap2_ctg_len'] = hit.ctg_len
                    alignment['minimap2_mlen'] = hit.mlen
                    alignment['minimap2_blen'] = hit.blen
                    alignment['minimap2_mapq'] = hit.mapq
                    alignment['minimap2_cigar'] = hit.cigar_str
                    alignment['minimap2_trans_strand'] = hit.trans_strand
                    if hit.r_st < hit.r_en :
                        alignment['minimap2_r_st'] = hit.r_st
                        alignment['minimap2_r_en'] = hit.r_en
                    else :
                        alignment['minimap2_r_st'] = hit.r_en
                        alignment['minimap2_r_en'] = hit.r_st
                    break
        for key in alignments_dict :
            alignments_dict[key].append(alignment[key])
    del aligner
    for key in alignments_dict :
        valid_table = valid_table.append_column(key, [alignments_dict[key]])
    mapped_table = pa.concat_tables([valid_table, skip_table], promote_options='default')
    if path_bed != None :
        mapped_table = assign_gene_ID_genomic(mapped_table, path_bed)
    else :
        mapped_table = assign_gene_ID_rna(mapped_table, path_ref)
    return mapped_table

def minimap2_table_from_file(file, path_ref, preset='splice', resume=False, overwrite=False, path_bed = None, skip_duplex_parents=True) :
    """
    Opens a parquet file and feeds the table within to minimap2_table().
    
    Args :
        file (str or Path) : path to the parquet file to be mapped.
        path_ref (str) : path to the sequence reference to be used.
        preset (str) : the preset to use for minimap2. Defaults to 'splice', which is designed to be splice-aware.
        resume (bool) : whether to resume mapping from a paused or broken run. This will only map files that don't already have minimap information in them, so it can't be used to continue a re-mapping session, as the files will all still have the original minimap data. Defaults to False.
        overwrite (bool) : whether to allow for overwriting mapped files. If True, it will remove any existing minimap data and continue with normal mapping. If False, it will skip over any files with minimap data.
    """
    table = pq.read_table(file)
    if 'minimap2_q_st' in table.column_names :
        if resume == True :
            print('skipping ', file, ', already mapped ')
            del table
            return
        elif overwrite == True :
            columns_to_drop = [ name for name in table.column_names if 'minimap2' in name ]
            table = table.drop_columns(columns_to_drop)
        else :
            raise ValueError("Error, this data has already been mapped with minimap2. Please set resume=True or overwrite=True if you'd like resume mapping a dataset or overwrite existing data.")
    print("started mapping ", file)
    table = minimap2_table( table, path_ref, preset=preset, path_bed=path_bed, skip_duplex_parents=skip_duplex_parents)
    pq.write_table(table, file)
    del table
    print("finished mapping ", file)
    return

def minimap2(dataset_dir, path_ref, preset='splice', resume=False, overwrite=False, workers = 4, path_bed = None) :
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
    files = [ x for x in Path(dataset_dir).iterdir() if x.is_file() ]
    with concurrent.futures.ProcessPoolExecutor( max_workers=workers ) as executor :
        futures = [ executor.submit( minimap2_table_from_file, file, path_ref, preset=preset, resume=resume, overwrite=overwrite, path_bed=path_bed ) for file in files ]
        concurrent.futures.wait( futures )
        print(futures)
    print("finished minimapping!")
    return

def strip_minimap_data(dataset_dir) :
    files = [x for x in Path(dataset_dir).iterdir() if x.is_file()]
    for file in files :
        print(file)
        table = pq.read_table(file)
        if 'minimap2_q_st' in table.column_names :
            print("stripping")
            columns_to_drop = [ name for name in table.column_names if 'minimap2' in name ]
            table = table.drop_columns(columns_to_drop)
            pq.write_table(table, file)
    print('Done')
    return

def count_mapped_reads(dataset_dir, path_out_csv, sample_dict, run_label = "None", max_match = 1) :
    """
    Counts the mapped reads in the dataset_dir and tallies up gene counts, collapsing gene and transcript variants into one count per gene ID. Gene IDs are defined by the name in parentheses in the minimap reference file. Saves the resulting counts as csv.
    
    Args :
        dataset_dir (str) : the path to the folder containing the parquet files to be counted. Should be path_out/pa_dataset/ where path_out is the same as what was used in build_parquet_dataset_from_sam.
        path_ref (str) : path to the sequence reference to be used. Currently intended to work with a transcripomic reference, essentially in the format of a fastq file where each row is a transcript with a contig ID, then sequence, then some tags. The critical information is the gene ID in parentheses, which must be the last parentheses in the first field of each row. ie XC001.4 Gene Name Etc. (version 1) (Gene ID), sequence, other information.
        path_out_csv (str) : the full path to where the results should be saved as csv. Will have a column per gene plus columns for barcode_ID, gene, and run_label. Each row denotes a different barcode/sample/run.
        sample_dict (dict) : a dictionary to define which barcodes belong to what sample names. Must be in the form of 'barcode ID' : 'sample name'.
        run_label (str) : a label that is assigned to all the reads in this dataset to potentially differentiate them from reads from other runs that might share the same barcode and sample names. Useful for combining reads from two runs with the same samples or if a run fails and is restarted.
    """
    dataset = ds.dataset(dataset_dir).filter( (~pc.field('minimap2_gene_ID').isin(['none', 'None', None])) & (~pc.field('barcode_ID').isin(['none matched', 'multiple'])) )
    if max_match != 1 :
        dataset = dataset.filter(pc.less_equal( pc.divide( pc.field('minimap2_mlen').cast(pa.float32()), pc.abs( pc.subtract( pc.field('minimap2_r_en'), pc.field('minimap2_r_st') ) ) ), pc.scalar(max_match) ) )
    if 'dx:i' in dataset.schema.names :
        dataset = dataset.filter( pc.field('dx:i') != '-1' )
        table = dataset.to_table( columns = ['barcode_ID', 'minimap2_gene_ID', 'dx:i', 'minimap2_ctg'] ).drop_null()
        table = table.drop_columns(['dx:i']).append_column('dx:i', table.column('dx:i').cast(pa.int8()))
        table_grouped = table.group_by(['barcode_ID', 'minimap2_gene_ID']).aggregate([ ([], 'count_all'), ('dx:i', 'sum') ])
        table_grouped = table_grouped.drop_columns(['count_all', 'dx:i_sum']).append_column('count_all', pc.add(table_grouped.column('count_all'), table_grouped.column('dx:i_sum')))
    else :
        table = dataset.to_table( columns = ['barcode_ID', 'minimap2_gene_ID', 'minimap2_ctg'] )
        table_grouped = table.group_by(['barcode_ID', 'minimap2_gene_ID']).aggregate([ ([], 'count_all') ])
    tables_by_barcode = []
    barcodes = table_grouped.column('barcode_ID').unique().to_pylist()
    for barcode in barcodes :
        row_values = table_grouped.filter(pc.field('barcode_ID') == barcode).column('count_all').to_pylist() + [run_label, sample_dict[barcode], barcode]
        gene_names = table_grouped.filter(pc.field('barcode_ID') == barcode).column('minimap2_gene_ID').to_pylist() + ['run_label', 'sample', 'barcode_ID']
        row_dict = {}
        for val, name in zip(row_values, gene_names) :
            row_dict[name] = [val]
        tables_by_barcode.append(pa.table(row_dict))
    counts = pa.concat_tables(tables_by_barcode, promote_options='default')
    counts = utils.fill_table_nulls(counts, 0)
    csv.write_csv(counts, path_out_csv)
    print('done')
    return

# def combine_degenerate_counts(counts) :
#     gene_IDs_proper = []
#     gene_IDs_degenerate = []
#     with open("/fs/ess/PAS2506/Users/Vlad/Ref_seqs/GRCh38_latest_genomic_v1.bed", 'r') as bed_handle :
#         for line in bed_handle.readlines() :
#             bed_split = line.split('\t')
#             if bed_split[0][0:2] == 'NC' :
#                 gene_IDs_proper.append(bed_split[3])
#             else :
#                 gene_IDs_degenerate.append(bed_split[3])
#     for gene in gene_IDs_degenerate :
#         if gene in counts :
#             if '-' in gene[5:] :
#                 last_dash = gene.rfind('-')
#                 try :
#                     end_num = int(gene[last_dash+1:])
#                     gene_base = gene[:last_dash]
#                     if gene_base in gene_IDs_proper and gene_base in counts :
#                         degen_values = counts[gene]
#                         proper_values = counts[gene_base]
#                         new_values = [ x+y for x, y in zip(degen_values, proper_values) ]
#                         del counts[gene]
#                         counts[gene_base] = new_values
#                 except :
#                     pass
#     return counts

# def count_mapped_reads_genomic(dataset_dir, path_out_csv, sample_dict, run_label = "None", max_match = 1) :
#     dataset = ds.dataset(dataset_dir)
#     if 'dx:i' in dataset.schema.names :
#         dataset = dataset.filter( pc.field('dx:i') != '-1' )
#     if max_match != 1 :
#         dataset = dataset.filter(pc.less_equal( pc.divide( pc.field('minimap2_mlen').cast(pa.float32()), pc.abs( pc.subtract( pc.field('minimap2_r_en'), pc.field('minimap2_r_st') ) ) ), pc.scalar(max_match) ) )
#     if 'dx:i' in dataset.schema.names :
#         table = dataset.to_table( columns = ['barcode_ID', 'minimap2_genomic_gene_ID', 'dx:i', 'minimap2_ctg'] )
#         table = table.drop_columns(['dx:i']).append_column('dx:i', table.column('dx:i').cast(pa.int8()))
#     else :
#         table = dataset.to_table( columns = ['barcode_ID', 'minimap2_genomic_gene_ID', 'minimap2_ctg'] )
#         table = table.append_column('dx:i', pa.nulls(table.num_rows).fill_nulls(0))
#     table_grouped = table.group_by(['barcode_ID', 'minimap2_genomic_gene_ID']).aggregate([ ([], 'count_all'), ('dx:i', 'sum') ])
#     barcodes = table.column('barcode_ID').unique().to_pylist()
#     genes = table.column('minimap2_genomic_gene_ID').unique().drop_null().to_pylist()
#     counts_dict = {}
#     for barcode in barcodes:
#         counts_dict[barcode] = {}
#         for gene in genes :
#             counts_dict[barcode][gene] = [0]
#     for row in table_grouped.to_pylist() :
#         if row['minimap2_genomic_gene_ID'] != None :
#             counts_dict[row['barcode_ID']][row['minimap2_genomic_gene_ID']] = [row['count_all'] + row['dx:i_sum']]
#     counts = False
#     for barcode in barcodes :
#         counts_in_barcode_table = pa.table(counts_dict[barcode]).add_column(0, 'barcode_ID', [[barcode]]).add_column(0, 'sample', [[sample_dict[barcode]]]).add_column(0, 'run_label', [[run_label]])
#         if not counts :
#             counts = counts_in_barcode_table
#         else :
#             counts = pa.concat_tables([counts, counts_in_barcode_table], promote_options='default')
#     # counts = combine_degenerate_counts(counts)
#     csv.write_csv(counts, path_out_csv)
#     print('done!')
#     return

def combine_counts(path_csvs, path_out_csv) :
    """
    Combine two counts csv files into a master table in a new location without modifying any counts. Intended to be the point where different runs are brought together for analysis. The resulting csv can be used in compare_counts().
    
    Args :
        path_csvs (list) : list of two paths to the csv files to be combined. Should be given as strings.
        path_out_csv (str) : path to where you want the resulting table to go.
    """
    table_1 = csv.read_csv(path_csvs[0], read_options = csv.ReadOptions(block_size = 10000000))
    table_2 = csv.read_csv(path_csvs[1], read_options = csv.ReadOptions(block_size = 10000000))
    table_combined = pa.concat_tables([table_1, table_2], promote_options='default')
    table_combined = pa.table([x.fill_null(0) for x in table_combined.itercolumns()], names=table_combined.column_names)
    csv.write_csv(table_combined, path_out_csv)
    print("Done combining counts")
    return

def add_counts(path_csvs, path_out_csv, new_run_label) :
    """
    Add together the counts within multiple counts csvs and generate a new csv table. This functions by grouping together values based on barcode and sample then adding them together per gene. Currently ignores any run_label tags and applies a new tag (new_run_label) to the new aggregated table. This is meant to combine counts from separate runs that are from the same samples ie restarting a broken sequencing run or generating additional data. 
    
    Args :
        path_csvs (list) : list of two paths to the csv files to be combined. Should be given as strings.
        path_out_csv (str) : path to where you want the resulting table to go.
        new_run_label (str) : the new run_label tag to apply to the resulting counts csv.
    """
    tables = [ csv.read_csv(path_csv, read_options = csv.ReadOptions(block_size = 10000000)) for path_csv in path_csvs ]
    new_table = pa.concat_tables(tables, promote_options='default')
    aggregations = [ (column_name, 'sum') for column_name in new_table.column_names if column_name not in ['sample', 'barcode_ID', 'run_label', 'none'] ]
    aggregated_table = pa.TableGroupBy(new_table, ['sample', 'barcode_ID'], use_threads = False).aggregate(aggregations)
    aggregated_table = aggregated_table.add_column( 0, 'run_label', [[ new_run_label for x in range(aggregated_table.num_rows) ]] )
    edited_names = aggregated_table.column_names[:3] + [ x[:-4] for x in aggregated_table.column_names[3:] ]
    aggregated_table = aggregated_table.rename_columns(edited_names)
    aggregated_table = utils.fill_table_nulls(aggregated_table, 0)
    csv.write_csv(aggregated_table, path_out_csv)
    print("Done adding counts")
    return