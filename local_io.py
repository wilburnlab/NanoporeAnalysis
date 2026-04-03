'''
Local I/O Functions
'''
import gzip
from pathlib import Path
import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
import pyarrow.compute as pc
from pyarrow import dataset
import concurrent
import concurrent.futures

class FileContentsError(RuntimeError):
    pass

class fastx_reader :
    def __init__(self, file, first_word = False, return_tags = False) :
        self.file = Path(file)
        self.first_word = first_word
        self.return_tags = return_tags
        file_opener = gzip.open if self.file.suffix == '.gz' else open
        self.reader = file_opener(self.file, 'rt')
        first_line = self.reader.readline()
        self.mode = 'a' if first_line[0] == '>' else 'q'
        self.reader.seek(0)
    def __iter__(self) :
        return self
    def __enter__(self) :
        return self
    def __exit__(self, exc_type, exc_value, traceback) :
        self.reader.close()
        return
    def __next__(self) :
        if self.mode == 'a': # FASTA
            seq = False
            while True :
                line = self.reader.readline()
                if line == '' :
                    if seq :
                        return name, seq
                    else :
                        raise StopIteration #end the iterator if blank line is found, which should only happen after reading the whole file
                if line[0] == '>': # New seq
                    if seq :
                        self.reader.seek(self.reader.tell() - len(line))
                        if self.return_tags :
                            return name, seq, tags
                        else :
                            return name, seq
                    line = line.rstrip()
                    if self.first_word :
                        split = line[1:].split(' ')
                        name = split[0]
                        if self.return_tags :
                            if len(split) > 1 :
                                tags = split[1:]
                            else :
                                tags = None
                    else :
                        name = line[1:]
                    seq = ''
                else:
                    if seq != False :
                        seq += line
        else: # FASTQ
            while True:
                ID_line = fastx.readline().rstrip()
                if len(ID_line) == 0:
                    raise StopIteration #end the iterator if blank line is found, which should only happen after reading the whole file
                if ID_line[0] == '@' :
                    ID_line = ID_line[1:]
                else :
                    raise Exception("Error: no '@' found at start of read ID line. The file may not match normal fastq format.")
                sequence = fastx.readline().rstrip()
                fastx.readline()
                score = decode_phred( fastx.readline().rstrip() )
                if self.first_word :
                    split = ID_line[1:].split(' ')
                    name = split[0]
                    if self.return_tags :
                        if len(split) > 1 :
                            tags = split[1:]
                        else :
                            tags = None
                else :
                    name = ID_line
                if self.return_tags :
                    return name, seq, score, tags
                else :
                    return name, seq, score

def decode_phred(score_str: str) -> np.array:
    '''
    Decode a quality score string into a numpy array
    '''
    #return np.array(list(bytes(score_str,'ascii')))
    return np.array( [ord(x) for x in score_str] )

def read_fastx(file_name, first_word = True, return_tags = False):
    '''
    Return dict with sequence data where keys are sequence names
    If FASTA: Values are sequences
    If FASTQ: Values are dicts with keys "Sequence" and "Score"
    '''
    seq_dict = {}
    with fastx_reader(file_name, first_word = first_word, return_tags = return_tags) as reader :
        if reader.mode == 'a' :
            if return_tags :
                seq_dict = { read[0] : (read[1], read[2]) for read in reader }
            else :
                seq_dict = { read[0] : read[1] for read in reader }
        elif reader.mode == 'q' :
            if return_tags :
                seq_dict = { read[0] : {
                    'Sequence' : read[1], 'Score' : read[2], 'Tags' : read[3]
                } for read in reader }
            else :
                seq_dict = { read[0] : {
                    'Sequence' : read[1], 'Score' : read[2]
                } for read in reader }
    # Ensure there are sequences in the file
    if not seq_dict:
        raise FileContentsError(f'No sequences found in FAST{mode.upper()} file: {file_name}')
    return seq_dict

def get_read_from_fastx(file, read_ID, first_word = True) :
    """
    Finds a specified read in a fasta or fastq file.

    Args:
        file (str or Path) : file path to fasta or fastq file.
        read_ID (str) : the read ID to find.
        first_word (bool) : whether to trim the fastx file read IDs to just the first 'word', ie trim everything past any space character.

    Returns :
        ctg_seq (str) : if found, the read sequence.
    """
    ctg_seq = None
    with fastx_reader(file, first_word = True) as reader :
        for contig in reader :
            if contig[0] == read_ID :
                ctg_seq = contig[1]
                return ctg_seq
    raise Exception("Error: matched contig not in reference file")

def write_fastx(file_name: str | Path, 
                seq_dict: dict,
                chars_per_line: int = None,
                append: bool = False,
                reduced_name: bool = False):
    '''
    Write a dict of sequences to a multipart FASTA file.
    
    Args:
        file_name (str or Path): destination file.
        seq_dict (dict): Keys are sequence names and values are the 
            corresponding sequences.
        chars_per_line (int): max number of characters on every line of the FASTA file.
            Defaults to None (no split per line).
        append (bool): sets parameter to pass to the open() function. If False, mode='w' and overwrites
            the file; if True, then mode='a' and appends to an existing file. Defaults to False.  
    
    Returns:
        None
    '''
    mode = 'a' if append else 'w'
    file_name = Path(file_name)
    with open(file_name, mode) as fout:
        for name, content in seq_dict.items():
            if reduced_name:
                name = name.split(' ')[0]
            if type(content) == str :
                if chars_per_line:
                    seq = '\n'.join([content[i:i+chars_per_line] for i in range(0,len(content),chars_per_line)])
                fout.write(f'>{name}\n{seq}\n')
            else :
                seq = content['Sequence']
                if chars_per_line:
                    seq = '\n'.join([content[i:i+chars_per_line] for i in range(0,len(content),chars_per_line)])
                score = content['Score']
                if 'Tags' in content :
                    tags = content['Tags'].join(' ')
                    fout.write(f'@{name} {tags}\n{seq}\n+\n{score}\n')
                else :
                    fout.write(f'@{name}\n{seq}\n+\n{score}\n')
    return None

def sam_to_parquet(file, path_out, basename_template = None) :
    """
    Transfers data from .sam files to an Apache parquet database by going line-by-line through the .sam file and building a pyarrow table. 
        Generates ID, seq, seq_len, and qual tags at minimum, then automatically creates tags for any other content in the .sam file assuming
        the format abcd:xyz... where abcd is the name and anything x and beyond is the value. Also, this skips over all tags in the first ten
        columns of the .sam, other than the three explicitly pulled tags: ID, seq, and qual. The others do not seem to be used by Dorado.
        Outputs to path_out following the basename_template. Currently hardcoded to limit the number of entries in the files to 200000.
    
    Args :
        file (str or Path) : the .sam file to be converted.
        path_out (str) : the directory being used for output.
        basename_template (str) : the optional template to be used for naming parquet files. Requires {i} to be present, which will be replaced by a number
            for each different file created ie 'dataset_A_part_{i}.parquet' where {i} will be 0 for the first file, 1, 2, etc. for subsequent files if the dataset is too large for one file.
    """
    Path(path_out).mkdir(parents=True, exist_ok=True)
    table_dict = {
        'ID' : [],
        'seq' : [],
        'seq_len' : [],
        'qual' : []
    }
    file = Path(file)
    if basename_template == None :
        basename_template = str(file.stem + '_part-{i}.parquet')
    print("moving file ", file.stem)
    table_dict_keys = table_dict.keys()
    with open(file, 'r') as handle_read :
        i = 0
        j = 0
        handle_write = False
        for line in handle_read.readlines() :
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
                        table_dict[key] = pa.nulls(len(table_dict['ID'])).to_pylist()
                        if handle_write :
                            handle_write.close()
                        handle_reader = pq.ParquetFile(path_out + basename_template.replace('{i}', str(j)))
                        table = pa.table(table_dict)
                        handle_write = pq.ParquetWriter(path_out + basename_template.replace( '{i}', str(j) + 'temp' ), table.schema)
                        for batch in handle_reader.iter_batches(batch_size=1000) :
                            handle_write.write(batch.append_column(key, pa.nulls(batch.num_rows).to_pylist()))
                        handle_reader.close()
                        handle_write.close()
                        os.remove(path_out + basename_template.replace('{i}', str(j)))
                        os.rename(path_out + basename_template.replace( '{i}', str(j) + 'temp' ), path_out + basename_template.replace('{i}', str(j)))
                    table_dict[key].append(line_dict[key])
                i += 1
            if i != 0 and i % 1000 == 0 :
                table = pa.table(table_dict)
                if not handle_write :
                    handle_write = pq.ParquetWriter(path_out + basename_template.replace('{i}', str(j)), table.schema)
                handle_write.write( table )
                del table
                table_dict = {
                    'ID' : [],
                    'seq' : [],
                    'seq_len' : [],
                    'qual' : []
                }
            if i == max_file_len :
                table = pa.table(table_dict)
                handle_write.write( table )
                handle_write.close()
                del table
                table_dict = {
                    'ID' : [],
                    'seq' : [],
                    'seq_len' : [],
                    'qual' : []
                }
                i = 0
                j += 1
                handle_write = pq.ParquetWriter(path_out + basename_template.replace('{i}', str(j)))
        table = pa.table(table_dict)
        if table.num_rows > 0 :
            handle_write.write( table )
            handle_write.close()
            del table
        elif j == 0 :
            raise FileContentsError(f'Sam file {file} appears to be empty. Dorado may not have finished running or ran improperly.')
    return

def build_parquet_dataset_from_sam(path_sam, path_out, delete_sam_folder = False) :
    """
    Takes all .sam files in path_sam and builds an Apache parquet database under path_out.
    
    Args :
        path_out (str) : the directory being used for output.
        path_sam (str) : directory containing .sam files. Does not search recursively.
        delete_sam_folder (bool) : whether or not to delete the folder at path_sam after building the parquet dataset. Defaults to False for safety.
    """
    files = [x for x in Path(path_sam).iterdir() if x.is_file()]
    Path(path_out).mkdir(parents=True, exist_ok=True)
    for file in files :
        if file.suffix == '.sam' :
            sam_to_parquet( file, path_out, basename_template = str(file.stem + '_part-{i}.parquet') )
    if delete_sam_folder :
        shutil.rmtree(path_sam)
    return

def fastq_to_parquet(path_fastq, path_out, additional_tags = False, max_file_len = 200000, omit_qual = False) :
    Path(path_out).mkdir(parents=True, exist_ok=True)
    seq_dict = {
        'ID' : [],
        'seq' : [],
        'seq_len' : [],
        'qual' : [],
    }
    with open(path_fastq, 'rt') as handle_open:
        j = 0
        k = 0
        handle_write = False
        while True:
            name = handle_open.readline().rstrip()
            if len(name) == 0: break
            if name[0] == '@' :
                name = name[1:]
            sequence = handle_open.readline().rstrip()
            handle_open.readline()
            score = handle_open.readline().rstrip()
            seq_dict['ID'].append(name)
            seq_dict['seq'].append(sequence)
            seq_dict['seq_len'].append(len(sequence))
            seq_dict['qual'].append(score)
            if j != 0 and j % 100000 == 0 :
                if additional_tags :
                    for tag in additional_tags :
                        seq_dict[tag] = np.full(len(seq_dict['ID']), additional_tags[tag])
                if omit_qual :
                    del seq_dict['qual']
                table = pa.table(seq_dict)
                if not handle_write :
                    handle_write = pq.ParquetWriter( path_out + '/' + path_fastq.stem + '_' + str(k), table.schema )
                handle_write.write_table(table)
                del table
                seq_dict = {
                    'ID' : [],
                    'seq' : [],
                    'seq_len' : [],
                    'qual' : []
                }
            j += 1
            if j == max_file_len :
                if additional_tags :
                    for tag in additional_tags :
                        seq_dict[tag] = np.full(len(seq_dict['ID']), additional_tags[tag])
                if omit_qual :
                    del seq_dict['qual']
                table = pa.table(seq_dict)
                if not handle_write :
                    handle_write = pq.ParquetWriter( path_out + '/' + path_fastq.stem + '_' + str(k), table.schema )
                handle_write.write_table(table)
                handle_write.close()
                seq_dict = {
                    'ID' : [],
                    'seq' : [],
                    'seq_len' : [],
                    'qual' : []
                }
                k += 1
                handle_write = pq.ParquetWriter( path_out + '/' + path_fastq.stem + '_' + str(k), table.schema )
                del table
                j = 0
        if additional_tags :
            for tag in additional_tags :
                seq_dict[tag] = np.full(len(seq_dict['ID']), additional_tags[tag])
        if omit_qual :
            del seq_dict['qual']
        table = pa.table(seq_dict)
        if not handle_write :
            handle_write = pq.ParquetWriter( path_out + '/' + path_fastq.stem + '_' + str(k), table.schema )
        handle_write.write_table( table )
        handle_write.close()
        del table
    return

def build_parquet_dataset_from_fastq(path_fastq, path_out, workers = 1, additional_tags = False, max_file_len = 200000, omit_qual = False, fill_pre_map_tags = False) :
    """
    Reads fastq file(s) to build a parquet dataset in the specified path_out folder.
    Args :
        path_fastx (str or list(str,)) : path(s) pointing to the fastq file(s). Must be fastq format, not fasta. Must point to files, not folders.
        path_out (str) : path pointing to the folder where you'd like the parquet dataset to be built. This folder will be made if it doesn't exist.
        
    """
    Path(path_out).mkdir(parents=True, exist_ok=True)
    if type(path_fastq) == list :
        path_fastq = [ Path(x) for x in path_fastq ]
    else :
        path_fastq = [ Path(path_fastq) ]
    with concurrent.futures.ProcessPoolExecutor( max_workers=workers ) as executor :
        futures = []
        for i, path in enumerate(path_fastq) :
            if additional_tags :
                this_additional_tags = { key : vals[i] for key, vals in additional_tags.items() }
                futures.append(executor.submit( fastq_to_parquet, path, path_out, additional_tags = this_additional_tags, max_file_len = max_file_len, omit_qual = omit_qual ))
            else :
                futures.append(executor.submit( fastq_to_parquet, path, path_out, additional_tags = False, max_file_len = max_file_len, omit_qual = omit_qual ))
        concurrent.futures.wait( futures )
        print(futures[0].result())
    return

def export_table_to_fastq(table, path_output, append = False) :
    """
    Exports reads in the given table to a fastq file. Table must have columns: 'ID' for the read ID, 'seq' for the sequence, and 'qual' for the read quality.
    Args :
        table (pyarrow.Table) : table with reads to be output. Must have columns: 'ID' for the read ID, 'seq' for the sequence, and 'qual' for the read quality.
        path_output (str or Path) : file path to write to.
        append (bool) : whether or not to append the reads to the output file, if it already exists. Should not be set to True unless the file already exists.
    Returns None
    """
    read_dict = {}
    for read in table.to_pylist() :
        read_dict[read['ID']] = {
            'Sequence' : read['seq'],
            'Score' : read['qual']
        }
    write_fastx(path_output, read_dict, append = append)
    return

def export_fastq(data, path_output, aggregations = None, write_descriptive_file_names = True) :
    """
    Exports sequences into fastq format files. This function can take data (in the same table format used throughout this package) in the form of a path (to 
        a parquet file or dataset directory, such as that made by build_parquet_dataset_from_sam), a pyarrow table, or a pyarrow dataset, or a list of any of 
        these (lists cannot contain pyarrow InMemoryDataset). The data can be split into different files based on
        values of table columns given through aggregations: all unique combinations of the specified columns will be given separate fastq files. To be valid,
        data must be a tabular format with at least these columns: 'ID' for the read ID, 'seq' for the read sequence, and 'qual' for the read quality.
    Args :
        data (str, Path, pyarrow.Table, pyarrow.dataset.dataset, or a list of any of these) : path(s) pointing to valid parquet file(s) or pyarrow dataset(s),
            or pre-loaded table(s) or dataset(s). Can be given as a list containing any one or more of these formats (except InMemoryDataset).
        path_output (str or Path) : path pointing to the desired fastq output file, OR to a directory if using aggregations to split the data. Will throw an
            error if a file path is given while using aggregations. Giving a directory path without using aggregations will just create a file with no suffix.
        aggregations (list( str, )) : list of the column names within the data that you'd like the resulting fastq to be split by. Each unique combination of
            values for these columns will be given a unique fastq file. i.e. if this table is given:
            'A' | 'B' | 'ID' | 'seq'  | 'qual'
             1  | 'a' | '1'  | 'cat'  | '!!!'
             2  | 'b' | '2'  | 'dog'  | '@@@'
             3  | 'b' | '3'  | 'fish' | '@#$!'
             3  | 'b' | '4'  | 'ant'  | '$$$'
            and aggregations = ['A', 'B'], this function will create the files '1_a.fastq' with read 1, '2_a.fastq' with read 2, and '3_b' with reads 3 and 4.
            The order of items listed in aggregations will determine the order of values in the file name, so it'd be 'a_1.fastq' etc. if it were flipped in
            the given example. These long names are only given if write_descriptive_file_names = True.
        write_descriptive_file_names (bool) : whether or not to write the fastq file names with the unique values of the columns given to aggregations. If
            False, the files will be names numerically. This has no effect if aggregations = None.
    Returns None
    """
    if type(data) == list :
        ds = dataset.dataset([ dataset.dataset(path) for path in data ])
    elif type(data) == str or type(data) == pa.Table :
        ds = dataset.dataset(data)
    elif type(data) == dataset.dataset or type(data) == dataset.InMemoryDataset :
        ds = data
    else :
        raise TypeError("Error: input 'data' variable is not a list, pyarrow Table, or pyarrow Dataset")
    if aggregations != None :
        Path(path_output).mkdir(parents=True, exist_ok=True)
        table_aggregates = ds.to_table(columns = aggregations)
        grouped = table_aggregates.group_by(aggregations).aggregate([])
        for i, row in enumerate(grouped.to_pylist()) :
            append_fastq_file = False
            if write_descriptive_file_names :
                path_output_row = path_output + '/' + '_'.join(list(row.values())) + '.fastq'
            else :
                path_output_row = path_output + '/' + str(i) + '.fastq'
            ds_filtered = ds
            for tag in list(row.keys()) :
                ds_filtered = ds_filtered.filter( pc.field(tag) == row[tag] )
            for batch in ds_filtered.to_batches(columns=['ID', 'seq', 'qual']) :
                export_table_to_fastq(batch, path_output_row, append = append_fastq_file)
                append_fastq_file = True
    else :
        append_fastq_file = False
        for batch in ds.to_batches(columns=['ID', 'seq', 'qual']) :
            export_table_to_fastq(batch, path_output, append = append_fastq_file)
            append_fastq_file = True
    return


