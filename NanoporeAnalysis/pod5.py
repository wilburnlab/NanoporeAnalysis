#!/usr/bin/env python3

"""
pod5_split.py
-------------
Split a single POD5 file into multiple smaller POD5 files
approx. by file size.

Usage:
  ./pod5_split.py --input big_file.pod5 --chunk-size 2G --out-prefix chunk_
"""

import h5py
import os
import argparse
import math

def parse_size(size_str):
    """
    Parse a human-readable size string like 2G, 500M, 100K into bytes.
    """
    size_str = size_str.strip().upper()
    unit_map = {'K': 1024, 'M': 1024**2, 'G': 1024**3}
    if size_str[-1] in unit_map:
        return int(float(size_str[:-1]) * unit_map[size_str[-1]])
    else:
        return int(size_str)

def copy_group(src_group, dst_group):
    """
    Recursively copy one group from source to destination.
    Because reads are stored in groups, you want to copy the entire group.
    """
    for key, item in src_group.items():
        if isinstance(item, h5py.Group):
            # Create a new group in the destination
            new_group = dst_group.create_group(key)
            copy_group(item, new_group)
        elif isinstance(item, h5py.Dataset):
            # Copy the dataset
            dst_group.create_dataset(key, data=item[()])

def split_pod5(input_file, chunk_size, out_prefix):
    """
    Split the given POD5 file by target chunk size (bytes).
    """
    with h5py.File(input_file, 'r') as f_in:
        # Typically, the reads are located under something like f_in['Reads'].
        # Confirm your file structure (for example, might be f_in['Raw/Reads']).
        if 'Reads' not in f_in:
            raise RuntimeError(f"No 'Reads' group found in {input_file}. Check file structure.")

        read_group = f_in['Reads']
        read_names = list(read_group.keys())  # a list of read IDs, can be huge
        read_count = len(read_names)
        print(f"Total reads found in {input_file}: {read_count}")

        chunk_index = 1
        out_file = None
        f_out = None

        # We will open a new chunk file, copy reads until we are near chunk_size,
        # then close that file and proceed with the next chunk, etc.
        for i, read_id in enumerate(read_names, start=1):
            # If we don't have an open chunk file (or if it's grown too large),
            # open a new chunk file
            if not f_out:
                out_file = f"{out_prefix}{chunk_index}.pod5"
                f_out = h5py.File(out_file, 'w')
                # Create matching top-level structure
                f_out.create_group('Reads')

            # Copy the read group from input to the chunk
            src_group = read_group[read_id]
            dst_group = f_out['Reads'].create_group(read_id)
            copy_group(src_group, dst_group)

            # Check current chunk file size
            f_out.flush()  # make sure data is written
            current_size = os.path.getsize(out_file)
            if current_size >= chunk_size:
                # Close this chunk and move on
                print(f"Chunk {chunk_index} reached ~{current_size} bytes with {i} reads copied total.")
                f_out.close()
                f_out = None
                chunk_index += 1

        # If there's an open file at the end, close it
        if f_out:
            f_out.close()

def main():
    parser = argparse.ArgumentParser(description="Split a POD5 file by approximate file size.")
    parser.add_argument("--input", required=True, help="Input POD5 file path.")
    parser.add_argument("--chunk-size", default="2G", help="Target chunk size (e.g., 2G, 500M).")
    parser.add_argument("--out-prefix", default="chunk_", help="Prefix for output POD5 files.")
    args = parser.parse_args()

    chunk_size_bytes = parse_size(args.chunk_size)

    print(f"Splitting {args.input} into chunks of ~{chunk_size_bytes} bytes.")
    split_pod5(args.input, chunk_size_bytes, args.out_prefix)

if __name__ == "__main__":
    main()
