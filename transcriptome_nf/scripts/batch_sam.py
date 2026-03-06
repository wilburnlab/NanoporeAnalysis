"""
Split the sam file into multiple files for batching
"""

import sys
from pathlib import Path


def split_sam_file(input_sam: str, batch_size: int, output_dir: str):
    """
    Split an individual sam file
    """
    input_sam = Path(input_sam)
    output_path = Path(output_dir)
    fout = None
    with open(input_sam) as fin:
        header_lines = [fin.readline() for _ in range(3)]
        k = 0
        for i, line in enumerate(fin):
            if i % batch_size == 0:
                if fout is not None:
                    fout.close()
                out_filename = f"{input_sam.stem}_{k}.sam"
                file_out = output_path / out_filename
                fout = open(file_out, "w")  # noqa: SIM115
                for hline in header_lines:
                    fout.write(hline)
                k += 1
            fout.write(line)
        if fout is not None:
            fout.close()


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: split_sam.py <input_sam> <batch_size> <output_dir>")
        sys.exit(1)

    input_sam = sys.argv[1]
    batch_size = int(sys.argv[2])
    output_dir = sys.argv[3]

    input_path = Path(input_sam)
    if input_path.is_dir():  # Folder of multiple sam files
        sam_files = [f for f in input_path.iterdir() if f.suffix == ".sam"]
        for sam_file in sam_files:
            split_sam_file(str(sam_file), batch_size, output_dir)
    else:
        split_sam_file(input_sam, batch_size, output_dir)
