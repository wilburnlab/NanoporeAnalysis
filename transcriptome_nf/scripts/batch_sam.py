'''
Split the sam file into multiple files for batching
'''

import os
import sys

def split_sam_file(input_sam: str,
                   batch_size: int,
                   output_dir: str):
    '''
    Split an individual sam file
    '''
    sam_fields = ['QNAME', 'FLAG', 'RNAME', 'POS', 'MAPQ', 'CIGAR', 
                    'RNEXT', 'PNEXT', 'TLEN', 'sequence', 'QUAL']
    output_path = output_dir #os.path.join(os.environ['workDir'], output_dir)
    with open(input_sam) as fin:
        header_lines = [fin.readline() for _ in range(3)]
        k = 0
        for i, line in enumerate(fin):
            if i%batch_size == 0:
                if k > 0:
                    fout.close()
                out_filename = f"{input_sam.split('/')[-1][:-4]}_{k}.sam"
                file_out = os.path.join(output_path, out_filename)
                #file_out = f"{output_dir}/{out_filename}"
                fout = open(file_out, 'w')
                for hline in header_lines:
                    fout.write(hline)
                k+=1
            fout.write(line)
        fout.close()


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: split_sam.py <input_sam> <batch_size> <output_dir>")
        sys.exit(1)
    
    input_sam = sys.argv[1]
    batch_size = int(sys.argv[2])
    output_dir = sys.argv[3]

    #assert os.path.isdir(output_dir), f"{output_dir} doesn't exist"
    
    if os.path.isdir(input_sam): # Folder of multiple sam files
        sam_files = [f for f in os.listdir(input_sam) if f[-4:] == '.sam']
        for sam_file in sam_files:
            split_sam_file(f"{input_sam}/{sam_file}", 
                           batch_size, 
                           output_dir)
    else:
        split_sam_file(input_sam, batch_size, output_dir)



