    '''
From parquet files generate ORF and Protein quant tables
'''

import sys
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import pyarrow.dataset as ds

from NanoporeAnalysis.utils import return_count_dict
from NanoporeAnalysis.local_io import write_fastx




def process_pqts(pqt_paths: list,
                 out_prefix: str):
    
    # Extract relevant records from the dataset
    columns = ['Sample_ID', 'ORF', 'Protein']
    dataset = ds.dataset(pqt_paths)
    records = dataset.to_table(filter=((ds.field('cDNA status') == 'Complete') &
                                        ds.field('Protein').is_valid()), 
                               columns=columns).to_pylist()
    
    # Construct counts of each data type and build Protein x Sample array
    #print(records)
    counts = dict([(c,return_count_dict(records,c)) for c in columns])

    # Filter singleton proteins
    counts['Protein'] = dict([i for i in counts['Protein'].items() if i[1] > 1])
    proteins = list(counts['Protein'])
    n_proteins = len(proteins)

    # Generate Protein labels
    protein_labels = [f"P{i}" for i in range(n_proteins)]
    protein_to_label = dict(zip(proteins, protein_labels))
    label_to_protein = dict([i[::-1] for i in protein_to_label.items()])

    # Write out fasta of protein sequences
    protein_fasta_file = f"{out_prefix}_proteins.fasta"
    fasta_dict = dict([(l,{'Sequence':p, 'Quality':None, 'Tags':None})
                       for l, p in label_to_protein.items()])
    write_fastx(protein_fasta_file, fasta_dict)

    # Generate counts per protein per sample
    samples = list(counts['Sample_ID'])
    n_samples = len(samples)
    count_array = np.zeros((n_proteins,n_samples))

    protein_to_idx = dict(zip(proteins, range(n_proteins)))
    sample_to_idx = dict(zip(samples, range(n_samples)))
    for r in records:
        i = protein_to_idx.get(r['Protein'], None)
        j = sample_to_idx.get(r['Sample_ID'], None)
        if i is not None:
            count_array[i,j] += 1

    # Export counts to CSV
    protein_df = pd.DataFrame(count_array, columns=samples)
    protein_df['Protein'] = protein_labels
    protein_df['Sequence'] = proteins
    protein_df = protein_df[['Protein']+samples+['Sequence']] # Re-order the columns
    protein_tsv = f"{out_prefix}_protein-counts.tab"
    protein_df.to_csv(protein_tsv, sep='\t')

    


if __name__ == "__main__":
    #if len(sys.argv) != 3:
    #    print("Usage: process_parquet.py <out_prefix>")
    #    sys.exit(1)
    
    out_prefix = sys.argv[1]
    pqt_files = sys.argv[2:]

    process_pqts(pqt_files, out_prefix)


