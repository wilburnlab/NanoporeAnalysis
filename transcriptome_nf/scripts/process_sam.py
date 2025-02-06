'''
Read in a sam file, parse the read, return a parquet file
'''

# Libraries
import sys
import pyarrow as pa
import pyarrow.parquet as pq

from NanoporeAnalysis.local_io import read_fastx
from NanoporeAnalysis.utils import reverse_complement, return_best_orf
from NanoporeAnalysis.align import process_read


def process_sam_file(sam_file: str, 
                     primer_file: str,
                     output_file: str):
    '''
    Process sam file, return output
    '''
    
    ## Process primer file
    primers = read_fastx(primer_file)
    ssp = primers['SSP']['Sequence']
    primer3 = reverse_complement(primers[list(primers)[0]]['Sequence'][:25])
    p3_len = len(primer3)
    umis = dict([(f"BC{p[7:9]}",primers[p]['Sequence'][p3_len:p3_len+16]) 
                 for p in primers if p[:2] == 'dT'])

    # Process SAM file
    sam_fields = ['QNAME', 'FLAG', 'RNAME', 'POS', 'MAPQ', 'CIGAR', 
                  'RNEXT', 'PNEXT', 'TLEN', 'sequence', 'QUAL']
    records = []
    for line in open(sam_file):
        if line[0] == '@': continue
        items = line.rstrip().split('\t')
        sam_dict = dict(zip(sam_fields, items[:11]))
        for x in items[11:]:
            sam_dict[x[:4]] = x[5:]
        sam_dict['Source'] = sam_file[:sam_file.rfind('_')] # Remove batch ID
        sam_dict = process_read(sam_dict, ssp, primer3, umis)
        records.append(sam_dict)
        
    table = pa.Table.from_pylist(records)
    pq.write_table(table, output_file)
    


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: process_sam.py <input_sam> <primers file> <output_parquet>")
        sys.exit(1)
    
    input_sam = sys.argv[1]
    primer_file = sys.argv[2]
    output_parquet = sys.argv[3]
    
    process_sam_file(input_sam, primer_file, output_parquet)


