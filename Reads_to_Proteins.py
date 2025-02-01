'''
Script to process DNA reads into likely predicted proteins
'''

## Basic steps
#  1. Read FASTA/FASTX file
#  2. Predict most likely protein per sequence
#  3. More analysis...which let's make sure the first two work ok

# Libraries
import sys, time, argparse, pickle
from concurrent.futures import ProcessPoolExecutor

from NanoporeAnalysis.local_io import read_fastx
from NanoporeAnalysis.utils import orf_searcher, return_count_dict



def parse_args(args):
    parser = argparse.ArgumentParser(usage = 'Reads_to_proteins.py [optional arguments] input.fast[a/q] output.fasta',
                                    description = 'Predict likely proteome from de-multiplexed reads',
                                    formatter_class = argparse.RawTextHelpFormatter, )
    parser.add_argument('input_fastx',
                        type = str,
                        help = 'Path to input FASTX file')
    parser.add_argument('output_fasta',
                        type = str,
                        help = 'Path to output protein sequences (FASTA)',
                        default = 'output.fasta' )
    parser.add_argument('--n_cores',
                        type = int,
                        help = 'Number of cores to use in multiprocessing',
                        default = 1)
    return parser.parse_args(args)



def return_best_orf(sequence: str):
    candidate_orfs = orf_searcher(sequence,
                                  min_length = 30,
                                  both_strands = False,
                                  longest_only = True)
    if len(candidate_orfs) > 0:
        sorted_orfs = sorted(candidate_orfs, 
                             key=lambda x:x['Protein_length'], 
                             reverse=True)
        best_orf = sorted_orfs[0]
    else:
        best_orf = {}
    return best_orf



def main():
    start_time = time.time()
    args = parse_args(sys.argv[1:])

    reads = read_fastx(args.input_fastx)
    n_reads = len(reads)
    print(f"{n_reads} reads successfully read from {args.input_fastx} ({time.time()-start_time:.1f} s)")
    
    names, records = zip(*reads.items())
    sequences = [r['Sequence'] for r in records]
    batch_size = int(1e6)
    with ProcessPoolExecutor(args.n_cores) as executor:
        for i in range(0, len(sequences), batch_size):
            orfs = executor.map(return_best_orf, 
                                sequences[i:i+batch_size], 
                                chunksize=100)
            for name, orf in zip(names[i:i+batch_size], orfs):
                reads[name].update(orf)

    n_per_orf = return_count_dict(reads, 'ORF')
    n_orfs = n_reads - n_per_orf[None]
    f_orfs = n_orfs/n_reads
    print(f"{n_orfs} extracted ({f_orfs:.2%}; {time.time()-start_time:.1f}) s")


    n_per_protein = return_count_dict(reads, 'Protein')

    with open(f"{args.output_fasta}.n_per_orf.pkl", 'wb') as fout:
        pickle.dump(n_per_orf, fout)
    with open(f"{args.output_fasta}.n_per_protein.pkl", 'wb') as fout:
        pickle.dump(n_per_protein, fout)


    

if __name__ == "__main__":
    main()
