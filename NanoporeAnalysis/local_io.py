'''
Local I/O Functions
'''
import gzip
from pathlib import Path

class FileContentsError(RuntimeError):
    pass
        

def read_fastx(file_name: str | Path,
               first_word: bool) -> dict:
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
        for line in fastx:
            if mode == 'a': # FASTA
                if line[0] == '>':  # New seq
                    name = line[1:].rstrip()
                    if first_word:
                        name = name.split(' ')[0]  # Uniprot style
                    seq_dict[name] = ''
                else:
                    seq_dict[name] += line.rstrip()
            else: # FASTQ
                if line[0] == '@': # New seq
                    name = line[1:].rstrip()
                    if first_word:
                        name = name.spplit(' ')[0]
                    mode = 'Sequence'
                    seq_dict[name] = {'Sequence':'', 'Score':''}
                elif line[0] == '+':
                    mode = 'Score'
                else:
                    seq_dict[name][mode] = line.rstrip()

    # Ensure there are sequences in the file
    if not seq_dict:
        raise FileContentsError(f'No sequences found in FASTA file: "{file_name}"')

    return seq_dict
