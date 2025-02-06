'''
Basic utility functions
'''


import re
import regex
import time
import datetime
import numpy as np

from NanoporeAnalysis.reference_data import NUCLEOTIDES, RESIDUES, CODONS, CODON_TO_RESIDUE


def timer( start, ):
    return str( datetime.timedelta( seconds=round( time.time() - start ) ) )

def return_count_dict(records_input: dict | list,
                      key: str,
                      descending = True):
    '''
    For the provided record_dict whose keys are names and values are records,
    compute a count dictionary of the provided key within the records, and
    return the sorted dictionary by count (default is highest values first)
    '''
    #print(list(record_dict.items())[:10])
    if type(records_input) == type({}): # Need to get records from values
        records = list(records_input.values())
    else: # Assume list or other type of array
        records = list(records_input)
    record_values = [r.get(key,None) for r in records]
    unique_values = set(record_values)
    n_per_unique = dict(zip(unique_values, [0]*len(unique_values)))
    for v in record_values:
        n_per_unique[v] += 1
    n_per_unique = dict(sorted(n_per_unique.items(), 
                               key=lambda x:x[1], 
                               reverse=descending))
    return n_per_unique

def reverse_complement(sequence):
    """Return the reverse complement of a DNA nucleotide string"""
    rc_trans = str.maketrans('ATCGatcg', 'TAGCtagc')
    seq_trans = sequence.translate(rc_trans)
    rc = seq_trans[::-1]
    return rc

def identify_alphabet(sequence: str) -> str:
    '''
    Return the likely alphabet of a sequence as DNA, Protein, or Unknown
    '''
    characters = set(sequence)
    dna_check = characters <= set(NUCLEOTIDES)
    protein_check = characters <= set(RESIDUES)

    if dna_check:
        return 'DNA'
    elif protein_check:
        return 'Protein'
    else:
        return 'Unknown'

def decode_phred(score_str: str) -> np.array:
    '''
    Decode a quality score string into a numpy array
    '''
    return np.array([ord(x)-33 for x in score_str])


def encode_phred(phred_array: np.ndarray) -> str:
    # Convert each integer in the array back to a character
    return ''.join(chr(score+33) for score in phred_array)


def translate(sequence: str) -> str:
    '''
    Translate DNA sequence to protein
    '''

    observed_alphabet = identify_alphabet(sequence)
    assert observed_alphabet == 'DNA', 'Attempted to translate non-DNA sequence: ' + sequence

    dna_sequence = sequence.strip().upper()
    if len(dna_sequence) < 3:
        return ''

    obs_codons = re.findall('...', dna_sequence)
    obs_residues = [CODON_TO_RESIDUE.get(c, 'X') for c in obs_codons]
    return ''.join(obs_residues)


def orf_check(prot: str) -> bool:
    '''
    Check if a protein sequence is a perfect start->stop with no ambiguous characters
    '''
    return prot[0] == 'M' and prot[-1] == '.' and prot.count('.') == 1 and prot.count('X') == 0


def len_check(sequence: str,
              min_len: int,
              max_len: int) -> bool:
    '''
    Check if a sequence is within a given min/max len
    '''
    seq_len = len(sequence)
    return seq_len >= min_len and seq_len <= max_len

def return_null_orf():
    return dict([(k,None) for k in 
                 ['ORF','Protein','ORF_start','ORF_end',
                  'ORF_strand','Protein_length','n']])


def orf_searcher(dna_sequence: str,
                 min_length: int = 30,
                 both_strands: bool = False,
                 longest_only: bool = True,
                 dir=None) -> list:
    '''
    Generate a list of records describing potential ORFs in a given DNA sequence
    '''

    if dna_sequence is None:
        return []

    ORF_search = regex.compile('ATG(?:...){2,}?(?:TAG|TAA|TGA|TRA|TAR)')

    records = []
    for strand in range(both_strands + 1):
        if strand == 0:  # For strand
            sequence = dna_sequence
            strand_label = '+'
        else:           # Rev strand
            sequence = reverse_complement(dna_sequence)
            strand_label = '-'

        if type(sequence) != type(''):
            print(sequence, type(sequence))
        orfs = ORF_search.findall(sequence, overlapped=True)

        for orf in set(orfs):
            translation = translate(orf)
            if translation.count('.') != 1:
                continue  # Limit to 1 stop codon
            protein = translation[:-1]
            if len(protein) <= min_length:
                continue  # Ensure # of AAs >= min length

            if strand_label == '+':
                start_pos = dna_sequence.find(orf) + 1
            else:
                start_pos = dna_sequence.find(reverse_complement(orf)) + 1
            orf_length = len(orf)
            end_pos = start_pos + orf_length - 1

            record = return_null_orf()
            record['ORF'] = orf
            record['Protein'] = protein
            record['ORF_start'] = start_pos
            record['ORF_end'] = end_pos
            record['ORF_strand'] = strand_label
            record['Protein_length'] = len(protein)
            record['n'] = 1
            records.append(record)

    # Prune records of internal ORFs
    if longest_only:
        longest_orfs = dict(
            [((r['ORF_strand'], r['ORF_end']), {'Protein_length': 0}) for r in records])
        for record in records:
            key = (record['ORF_strand'], record['ORF_end'])
            if record['Protein_length'] > longest_orfs[key]['Protein_length']:
                longest_orfs[key] = record
        records = sorted(longest_orfs.values(), 
                         key=lambda x:x['Protein_length'], 
                         reverse=True)

    if not dir:
        return records
    else:
        # Dump the records to the provided directory
        out_filename = sequence_name.replace('/', '_').replace('\\', '_') + '.pkl'
        out_path = os.path.join(dir, out_filename)
        pickle.dump(records, open(out_path, 'wb'))
        return None


def return_best_orf(dna_sequence: str,
                    min_length: int = 30,
                    both_strands: bool = False,
                    longest_only: bool = True):
    candidate_orfs = orf_searcher(dna_sequence,
                                  min_length,
                                  both_strands,
                                  longest_only)
    if len(candidate_orfs) > 0:
        sorted_orfs = sorted(candidate_orfs, 
                             key=lambda x:x['Protein_length'], 
                             reverse=True)
        best_orf = sorted_orfs[0]
    else:
        best_orf = return_null_orf()
    return best_orf