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


from pathlib import Path
from pod5 import Reader
import polars as pl


import uuid
from operator import attrgetter
from contextlib import suppress
import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq

def return_pod5_schema(id_encoding: str = "uuid16",          
                       calibration: str = "scale_shift"):
    '''
    Return a pyarrow schema for processing POD5 files

    id_encoding is one of ("uuid16" | "binary" | "string")
    calibration is one of ("scale_shift" | "range_digitisation_offset")
    '''

    # Map encodings to Arrow dtypes (only the dtype varies)
    id_type_map = {"uuid16": pa.binary(16),
                   "binary": pa.binary(),
                   "string": pa.string()}
    try:
        id_type = id_type_map[id_encoding]
    except KeyError:
        raise ValueError(f"id_encoding must be one of {list(id_type_map)}, got {id_encoding!r}")

    # Map calibration mode to (name, dtype) pairs
    calib_cols_map = {"scale_shift": [("scale", pa.float32()),
                                      ("shift", pa.float32())],
                      "range_digitisation_offset": [("range", pa.float32()),
                                                    ("digitisation", pa.float32()),
                                                    ("offset", pa.float32())]}
    try:
        calib_cols = calib_cols_map[calibration]
    except KeyError:
        raise ValueError(f"calibration must be one of {list(calib_cols_map)}, got {calibration!r}")

    # Build the schema fields
    fields = [pa.field("read_id", id_type, nullable=False),
              #pa.field("read_id_str", pa.string(), nullable=False),
              pa.field("sample_rate", pa.int32(), nullable=False),
              pa.field("length",      pa.int32(),  nullable=False),
              *(pa.field(name, dtype, nullable=False) for name, dtype in calib_cols),
              pa.field("source_file", pa.string(), nullable=False),
              pa.field("read_index",  pa.int64(),  nullable=False),
              pa.field("signal",      pa.large_list(pa.int16()), nullable=False)]

    return pa.schema(fields)




def _to_uuid(obj) -> uuid.UUID:
    """Strictly coerce a POD5 read_id to uuid.UUID."""
    if isinstance(obj, uuid.UUID):
        return obj
    if isinstance(obj, (bytes, bytearray)) and len(obj) == 16:
        return uuid.UUID(bytes=bytes(obj))
    if isinstance(obj, str):
        return uuid.UUID(obj)  # canonical UUID text
    raise TypeError(f"read_id must be UUID, 16-byte UUID, or UUID string, got {type(obj)}")


def encode_read_id(read_or_id, id_encoding: str = "uuid16"):
    """
    Encode a POD5 read_id with minimal branching.
    id_encoding: 'uuid16' | 'binary' | 'string'
    """
    rid = getattr(read_or_id, "read_id", read_or_id)
    u = _to_uuid(rid)
    match id_encoding:
        case "uuid16" | "binary":
            return u.bytes
        case "string":
            return str(u)
        case _:
            raise ValueError("id_encoding must be 'uuid16' | 'binary' | 'string'")


#def read_calibration(read, mode: str = "scale_shift"):
#    if mode == "scale_shift":
#        return float(read.scale), float(read.shift)
#    elif mode == "range_digitisation_offset":
#        return float(read.range), float(read.digitisation), float(read.offset)
#    else:
#        raise ValueError("mode must be 'scale_shift' or 'range_digitisation_offset'")



def get_first_attr(obj, *paths):
    """
    Return the first non-None attribute among dotted paths, or None if none found.
    """
    for p in paths:
        with suppress(AttributeError):
            val = attrgetter(p)(obj)  # supports "calibration.range", etc.
            if val is not None:
                return val
    return None

def read_calibration(read, mode: str = "scale_shift"):
    if mode == "scale_shift":
        scale = get_first_attr(read,
                               "predicted_scaling.scale",
                               "tracked_scaling.scale",
                               "calibration.scale",
                               "scale")
        shift = get_first_attr(read,
                               "predicted_scaling.shift",
                               "tracked_scaling.shift",
                               "calibration.shift",
                               "offset",  # some builds expose offset equivalent
                               "shift")
        if scale is not None and shift is not None:
            return float(scale), float(shift)

        # fallback to R/D/O → scale/shift
        R = get_first_attr(read, "calibration_range", "calibration.range")
        D = get_first_attr(read, "calibration_digitisation", "calibration.digitisation")
        O = get_first_attr(read, "calibration.offset")
        if None not in (R, D, O):
            scale = float(R) / float(D)
            shift = float(O) * scale
            return scale, shift

        raise AttributeError("Missing scale/shift (and cannot derive from range/digitisation/offset).")

    elif mode == "range_digitisation_offset":
        R = get_first_attr(read, "calibration_range", "calibration.range")
        D = get_first_attr(read, "calibration_digitisation", "calibration.digitisation")
        O = get_first_attr(read, "calibration.offset", "calibration.shift")  # some APIs use 'shift'
        if None in (R, D, O):
            raise AttributeError("Missing range/digitisation/offset on this POD5 version.")
        return float(R), float(D), float(O)

    raise ValueError("mode must be 'scale_shift' or 'range_digitisation_offset'")


def read_to_record(read, 
                   i: int, 
                   source_name: str,
                   id_encoding: str = "uuid16",
                   calibration: str = "scale_shift"):
    '''
    Process read from pod5 into record
    '''
    #print(help(read))
    scalars = {'read_id': encode_read_id(read, id_encoding),
               #'read_id_str' : read.read_id,
               'sample_rate' : int(read.run_info.sample_rate),
               'signal': read.signal,
               'length': int(read.signal.size),
               'source_file': source_name,
               'read_index': int(i),
               }
    cal = dict(zip(calibration.split('_'), 
                   read_calibration(read, calibration)))
    return scalars|cal








def records_to_table(records, 
                     schema):
    '''
    Process records into pyarrow table, write to parquet if out_path is provided
    '''
    # Build scalar arrays from records (everything except 'signal')
    scalar_cols = [n for n in schema.names if n != "signal"]
    arrays = {k: pa.array([r[k] for r in records]) for k in scalar_cols}
    # Build the list column directly from list of NumPy arrays
    arrays["signal"] = pa.array([r["signal"] for r in records],
                        type=pa.large_list(pa.int16()))
    # Write the table and reset records
    table = pa.table(arrays, schema=schema)
    return table






def pod5_to_parquet(input_pod5,
                    out_dir,
                    id_encoding: str = "uuid16",
                    calibration: str = "scale_shift",
                    batch_size: int = 4096,
                    compression: str = "zstd",
                    compression_level: int = 4,
                    split_by_batch: bool = False,   # False: single file; True: part-00000.pqt, ...
                    file_prefix: str | None = None):   # override output filename prefix
    """
    Stream a POD5 into a Parquet dataset with metadata + signal in one table.

    Produces either:
      - one Parquet file kept open and appended with row groups (split_by_batch=False), or
      - many part files (split_by_batch=True), one per batch.

    Returns: Path to the written file or directory.
    """
    input_pod5 = Path(input_pod5)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    schema = return_pod5_schema(id_encoding, calibration)
    source_name = input_pod5.name
    prefix = file_prefix or source_name.rsplit(".", 1)[0]

    part_idx = 0
    writer = None
    if split_by_batch:
        out_path = out_dir/f"{prefix}.part-{part_idx:05d}.parquet"
    else:
        out_path = out_dir/f"{prefix}.parquet"
        writer = pq.ParquetWriter(out_path, 
                                  schema, 
                                  compression=compression, 
                                  compression_level=compression_level, 
                                  use_dictionary=True)

    def write_table():
        nonlocal writer, part_idx
        if not records:
            return
        table = records_to_table(records,schema)
        if split_by_batch:
            pq.write_table(table, 
                           out_path,
                           compression=compression,
                           compression_level=compression_level,
                           use_dictionary=True,
                           row_group_size=len(table))
        else:
            writer.write_table(table, row_group_size=len(table)) # Append to existing table
                    
    records = []
    with Reader(input_pod5) as r:
        for i, read in enumerate(r.reads()):
            records.append(read_to_record(read,i,source_name,id_encoding,calibration))
            if len(records) >= batch_size:
                write_table()
                records.clear()
                if split_by_batch:
                    part_idx += 1
                    out_path = out_dir/f"{prefix}.part-{part_idx:05d}.parquet"
        # Write out final batch
        if len(records) > 0:
            write_table()
            records.clear()

    if writer is not None:
        writer.close()




























def pod5_read_to_record(read):
    return {'rid' : str(read.read_id),
            'sr' : read.run_info.sample_rate, # Hz
            'adc' : read.signal, # raw ADC counts (int16)
            'pa' : read.signal_pa, # already applies per-read calibration
            't' : np.arange(read.signal_pa.size) / read.run_info.sample_rate}


def read_pod5(file_name: str|Path,
              out_file: str|Path = None):
    '''
    Read a pod5 file to a polars dataframe
    '''
    file_name = Path(file_name)
    records = {'rid':[], 'sr':[], 'size':[], 'adc':[], 'pa':[]}
    with Reader(file_name) as r:
        for read in r.reads():
            records['rid'].append(str(read.read_id)) # name
            records['sr'].append(read.run_info.sample_rate) # Hz
            records['size'].append(read.signal_pa.size) # N
            records['adc'].append(read.signal) # raw ADC counts (int16)
            records['pa'].append(read.signal_pa) # applies per-read calibration, picoAmp

    if out_file is not None:
        out_file = Path(out_file)
        df.write_parquet(out_file)
    return pl.DataFrame(records)



    #        pod5_dict = pod5_read_to_record(read)
    #        records.get(pod5_dict['rid'],{}).update(pod5_dict)









##
## Old code from Chat that I don't think we need anymore
## Consider deleting
##

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
