'''
Local I/O Functions
'''
import gzip
from pathlib import Path
import numpy as np

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
            for line in fastx:
                line = line.rstrip()
                if line[0] == '>': # New seq
                    name = line[1:].split(' ')[0] if first_word else line[1:]
                    seq_dict[name] = ''
                else:
                    seq_dict[name] += line
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
    fout = open( file_name, mode )
    with open(file_name, mode) as fout:
        for name, seq in seq_dict.items():
            if reduced_name:
                name = name.split(' ')[0]
            if chars_per_line:
                seq = '\n'.join([seq[i:i+chars_per_line] for i in range(0,len(seq),chars_per_line)])
            fout.write(f'>{name}\n{seq}\n')
    return None
            
'''           
            
            if mode == 'a': # FASTA
                if line[0] == '>':  # New seq
                    name = line[1:].split(' ')[0] if first_word else line[1:]
                    seq_dict[name] = ''
                else:
                    seq_dict[name] += line
            else: # FASTQ
                if line[0] == '@': # New seq
                    name = line[1:].split(' ')[0] if first_word else line[1:]
                    if name == '0**(&\'&$#%#$%#$%#%.0102-&%$+++,:45<7,A<==?CDAA96?A1-.,/6577:9:8:32**227::2/6.3569=@CFHJKEEFIJM5436464660123==7400/0,++<@JIJC?=?8.-0399;7>85CJ{{=A?;6306<:72114/--=98FCB>=+*,::+-9<B==./<C?<=0/409:8889>:998,,/(((9688D60%*)3302:30\',1;A?<8;0,\'%&2;J{{FDF=,+/>9?FQH?54375D=BDGDKJEA779+()+2*())+4859;:::BDF<;=975549867>>:/)&&\'\'(%(%(#$$$$%&$%%&\'$#"$%,*.,03,,,:>:75270\'(**:A1*)*/,,&)$%()9=<?==<<D436F{:KL628?N90)((4)EG{999==---/1*(%+7>/..$##$)\'(\'\'%\'\'((138?;?<=8-)*();C@=<<()@DKG=>BIAEFFFPADIK;::B9C3AA?<8=;;I?>MGIIAB4448;GHM-*)*733=32-***((##"%(;A9:5=399/-182-\'&&\')04\'512()()(\'()(\'&548;125/6622.14A)(31:)9560**,23)&%&%#&9.*+55/3./:21DCBA==2/.+)(),-800*))*//*\',)(,,0.&&)+*2;>621.$$#%\'29?;-\'%"%&)-\'%&117803:431.6879>6=<<:7:FKI94;877))+*<3-**0.,%)*).0..897,+-126+*,.33767/\'(((+':
                        print( 'ACK' )
                        print( )
                    mode = 'Sequence'
                    seq_dict[name] = {}
                elif line[0] == '+':
                    mode = 'Score'
                else:
                    seq_dict[name][mode] = line if mode == 'Sequence' else decode_phred(line)
                    #if mode == 'Sequence':
                    #    
                    #else:
                    #    seq_dict[name][mode] = decode_phred(line)
                        #seq_dict[name][mode] += [ ord(c) for c in line.rstrip() ]
                        #for c in line.rstrip(): seq_dict[name][mode].append(ord(c))

    

    return seq_dict
'''