"""
POD5 signal processing and Parquet conversion
"""

import uuid
from contextlib import suppress
from operator import attrgetter
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
from pod5 import Reader


def return_pod5_schema(id_encoding: str = "uuid16", calibration: str = "scale_shift"):
    """
    Return a pyarrow schema for processing POD5 files

    id_encoding is one of ("uuid16" | "binary" | "string")
    calibration is one of ("scale_shift" | "range_digitisation_offset")
    """

    # Map encodings to Arrow dtypes (only the dtype varies)
    id_type_map = {"uuid16": pa.binary(16), "binary": pa.binary(), "string": pa.string()}
    try:
        id_type = id_type_map[id_encoding]
    except KeyError as err:
        raise ValueError(f"id_encoding must be one of {list(id_type_map)}, got {id_encoding!r}") from err

    # Map calibration mode to (name, dtype) pairs
    calib_cols_map = {
        "scale_shift": [("scale", pa.float32()), ("shift", pa.float32())],
        "range_digitisation_offset": [
            ("range", pa.float32()),
            ("digitisation", pa.float32()),
            ("offset", pa.float32()),
        ],
    }
    try:
        calib_cols = calib_cols_map[calibration]
    except KeyError as err:
        raise ValueError(f"calibration must be one of {list(calib_cols_map)}, got {calibration!r}") from err

    # Build the schema fields
    fields = [
        pa.field("read_id", id_type, nullable=False),
        # pa.field("read_id_str", pa.string(), nullable=False),
        pa.field("sample_rate", pa.int32(), nullable=False),
        pa.field("length", pa.int32(), nullable=False),
        *(pa.field(name, dtype, nullable=False) for name, dtype in calib_cols),
        pa.field("source_file", pa.string(), nullable=False),
        pa.field("read_index", pa.int64(), nullable=False),
        pa.field("signal", pa.large_list(pa.int16()), nullable=False),
    ]

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
    if id_encoding in ("uuid16", "binary"):
        return u.bytes
    elif id_encoding == "string":
        return str(u)
    else:
        raise ValueError("id_encoding must be 'uuid16' | 'binary' | 'string'")


# def read_calibration(read, mode: str = "scale_shift"):
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
        scale = get_first_attr(read, "predicted_scaling.scale", "tracked_scaling.scale", "calibration.scale", "scale")
        shift = get_first_attr(
            read,
            "predicted_scaling.shift",
            "tracked_scaling.shift",
            "calibration.shift",
            "offset",  # some builds expose offset equivalent
            "shift",
        )
        if scale is not None and shift is not None:
            return float(scale), float(shift)

        # fallback to range/digitisation/offset → scale/shift
        cal_range = get_first_attr(read, "calibration_range", "calibration.range")
        cal_digit = get_first_attr(read, "calibration_digitisation", "calibration.digitisation")
        cal_offset = get_first_attr(read, "calibration.offset")
        if None not in (cal_range, cal_digit, cal_offset):
            scale = float(cal_range) / float(cal_digit)
            shift = float(cal_offset) * scale
            return scale, shift

        raise AttributeError("Missing scale/shift (and cannot derive from range/digitisation/offset).")

    elif mode == "range_digitisation_offset":
        cal_range = get_first_attr(read, "calibration_range", "calibration.range")
        cal_digit = get_first_attr(read, "calibration_digitisation", "calibration.digitisation")
        cal_offset = get_first_attr(read, "calibration.offset", "calibration.shift")  # some APIs use 'shift'
        if None in (cal_range, cal_digit, cal_offset):
            raise AttributeError("Missing range/digitisation/offset on this POD5 version.")
        return float(cal_range), float(cal_digit), float(cal_offset)

    raise ValueError("mode must be 'scale_shift' or 'range_digitisation_offset'")


def read_to_record(read, i: int, source_name: str, id_encoding: str = "uuid16", calibration: str = "scale_shift"):
    """
    Process read from pod5 into record
    """
    # print(help(read))
    scalars = {
        "read_id": encode_read_id(read, id_encoding),
        #'read_id_str' : read.read_id,
        "sample_rate": int(read.run_info.sample_rate),
        "signal": read.signal,
        "length": int(read.signal.size),
        "source_file": source_name,
        "read_index": int(i),
    }
    cal = dict(zip(calibration.split("_"), read_calibration(read, calibration)))
    return scalars | cal


def records_to_table(records, schema):
    """
    Process records into pyarrow table, write to parquet if out_path is provided
    """
    # Build scalar arrays from records (everything except 'signal')
    scalar_cols = [n for n in schema.names if n != "signal"]
    arrays = {k: pa.array([r[k] for r in records]) for k in scalar_cols}
    # Build the list column directly from list of NumPy arrays
    arrays["signal"] = pa.array([r["signal"] for r in records], type=pa.large_list(pa.int16()))
    # Write the table and reset records
    table = pa.table(arrays, schema=schema)
    return table


def pod5_to_parquet(
    input_pod5,
    out_dir,
    id_encoding: str = "uuid16",
    calibration: str = "scale_shift",
    batch_size: int = 4096,
    compression: str = "zstd",
    compression_level: int = 4,
    split_by_batch: bool = False,  # False: single file; True: part-00000.pqt, ...
    file_prefix: str | None = None,
):  # override output filename prefix
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
        out_path = out_dir / f"{prefix}.part-{part_idx:05d}.parquet"
    else:
        out_path = out_dir / f"{prefix}.parquet"
        writer = pq.ParquetWriter(
            out_path, schema, compression=compression, compression_level=compression_level, use_dictionary=True
        )

    def write_table():
        nonlocal writer, part_idx
        if not records:
            return
        table = records_to_table(records, schema)
        if split_by_batch:
            pq.write_table(
                table,
                out_path,
                compression=compression,
                compression_level=compression_level,
                use_dictionary=True,
                row_group_size=len(table),
            )
        else:
            writer.write_table(table, row_group_size=len(table))  # Append to existing table

    records = []
    with Reader(input_pod5) as r:
        for i, read in enumerate(r.reads()):
            records.append(read_to_record(read, i, source_name, id_encoding, calibration))
            if len(records) >= batch_size:
                write_table()
                records.clear()
                if split_by_batch:
                    part_idx += 1
                    out_path = out_dir / f"{prefix}.part-{part_idx:05d}.parquet"
        # Write out final batch
        if len(records) > 0:
            write_table()
            records.clear()

    if writer is not None:
        writer.close()
