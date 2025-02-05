'''
Local I/O Functions
'''
import gzip
from pathlib import Path
import numpy as np
import pyarrow as pa
from pyarrow import dataset

class FileContentsError(RuntimeError):
    pass
        
def decode_phred(score_str: str) -> np.array:
    '''
    Decode a quality score string into a numpy array
    '''
    #return np.array(list(bytes(score_str,'ascii')))
    return np.array( [ord(x) for x in score_str] )

def read_fastx(file_name: str | Path,
               first_word: bool = False) -> dict:
    '''
    Return dict with sequence data where keys are sequence names
    If FASTA: Values are sequences
    If FASTQ: Values are dicts with keys "Sequence" and "Score"
    '''
    seq_dict = { }
    file_name = Path(file_name)
    file_opener = gzip.open if file_name.suffix == '.gz' else open
    with file_opener(file_name, 'rt') as fastx:
        first_line = fastx.readline()
        assert first_line[0] in ['>','@'], 'Invalid file '+file_name
        mode = 'a' if first_line[0] == '>' else 'q'
        fastx.seek(0)
        if mode == 'a': # FASTA
            seq = False
            for line in fastx:
                line = line.rstrip()
                if line != '' :
                    if line[0] == '>': # New seq
                        if seq :
                            seq_dict[name] = seq
                        name = line[1:].split(' ')[0] if first_word else line[1:]
                        seq = ''
                    else:
                        seq += line
            seq_dict[name] = seq
        else: # FASTQ
            while True:
                name = fastx.readline().rstrip()
                sequence = fastx.readline().rstrip()
                fastx.readline()
                score = decode_phred( fastx.readline().rstrip() )
                if len(name) == 0: break
                if '\t' in name :
                    split = name.split('\t')
                    seq_dict[split[0]] = { 'Sequence' : sequence, 'Score' : score , 'Tags' : split[1:]}
                else : 
                    seq_dict[name] = { 'Sequence' : sequence, 'Score' : score }
                    
    # Ensure there are sequences in the file
    if not seq_dict:
        raise FileContentsError(f'No sequences found in FAST{mode.upper()} file: {file_name}')
    
    return seq_dict

def write_fasta(file_name: str | Path, 
                seq_dict: dict,
                chars_per_line: int = None,
                append: bool = False,
                reduced_name: bool = False):
    # DBW comment 11/18/23
    # The "mode" variable was replaced with "append" to be more independently descriptive and not rely on the syntax 
    # of the open function. I had to change line 283 in Analysis.py that used this. Also, nchars has been replaced 
    # with chars_per_line and is now also togglable, so that we have the option to produce fasta files with no line 
    # breaks in the sequence. The 80 bases per line norm is a relic of when the FASTA format was originally defined 
    # back in the ~1980s and terminal widths were fixed size, and I don't know of any modern software packages that
    # require fixed width data. It's a good feature to keep for flexibility and visualization. I also added a
    # "reduced_name" (often useful when dealing with Uniprot-style data)
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
        for name, seq in seq_dict.items():
            if reduced_name:
                name = name.split(' ')[0]
            if chars_per_line:
                seq = '\n'.join([seq[i:i+chars_per_line] for i in range(0,len(seq),chars_per_line)])
            fout.write(f'>{name}\n{seq}\n')
    return None

def sam_to_parquet(file, path_out, basename_template = None) :
    """
    Transfers data from .sam files to an Apache parquet database by going line-by-line through the .sam file and building a pyarrow table. 
        Generates ID, seq, seq_len, and qual tags at minimum, then automatically creates tags for any other content in the .sam file assuming
        the format abcd:xyz... where abcd is the name and anything x and beyond is the value. Also, this skips over all tags in the first ten
        columns of the .sam, other than the three explicitly pulled tags: ID, seq, and qual. The others do not seem to be used by Dorado.
        Outputs to path_out following the basename_template. Current hardcoded to limit the number of entries in the files to 200000.
    
    Args :
        file (str or Path) : the .sam file to be converted.
        path_out (str) : the directory being used for output.
        basename_template (str) : the optional template to be used for naming parquet files. Requires {i} to be present, which will be replaced by a number
            for each different file created ie 'dataset_A_part_{i}.parquet' where {i} will be 0 for the first file, 1, 2, etc. for subsequent files if the dataset is too large for one file.
    """
    table_dict = {
        'ID' : [],
        'seq' : [],
        'seq_len' : [],
        'qual' : []
    }
    print("moving file ", file.stem)
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
    dataset.write_dataset(table, Path(path_out), format='parquet', basename_template = basename_template, max_rows_per_file = 200000, max_rows_per_group = 200000, existing_data_behavior='overwrite_or_ignore')
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
    for file in files :
        sam_to_parquet( file, path_out, basename_template = str(file.stem + '_part-{i}.parquet') )
    if delete_sam_folder :
        shutil.rmtree(path_sam)
    return