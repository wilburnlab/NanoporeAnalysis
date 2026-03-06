"""
Local I/O Functions
"""

import gzip
from pathlib import Path

import polars as pl
from pod5 import Reader


class FileContentsError(RuntimeError):
    pass


def read_fastx(file_name: str | Path, first_word: bool = False) -> dict:
    """
    Return dict with sequence data where keys are sequence names
    If FASTA: Values are sequences
    If FASTQ: Values are dicts with keys "Sequence" and "Score"
    """
    seq_dict = {}
    file_name = Path(file_name)
    file_opener = gzip.open if file_name.suffix == ".gz" else open
    with file_opener(file_name, "rt") as fastx:
        # Determine file type
        first_line = fastx.readline()
        assert first_line[0] in [">", "@"], f"Invalid file {file_name}"
        mode = "a" if first_line[0] == ">" else "q"
        fastx.seek(0)

        if mode == "a":  # FASTA
            for line in fastx:
                line = line.rstrip()
                if len(line) > 0 and line[0] == ">":  # New seq
                    name = line[1:].split(" ")[0] if first_word else line[1:]
                    seq_dict[name] = {"Sequence": "", "Quality": None, "Tags": None}
                else:
                    seq_dict[name]["Sequence"] += line
        else:  # FASTQ
            while True:
                header = fastx.readline().rstrip()
                if len(header) == 0:
                    break
                sequence = fastx.readline().rstrip()
                fastx.readline()
                qual = fastx.readline().rstrip()

                if header.find("\t") >= 0:
                    parts = header.split("\t")
                    name = parts[0][1:]
                    tags = parts[1:]
                else:
                    name = header[1:]
                    tags = []
                seq_dict[name] = {"Sequence": sequence, "Quality": qual, "Tags": tags}

    # Ensure there are sequences in the file
    if not seq_dict:
        raise FileContentsError(f"No sequences found in FAST{mode.upper()} file: {file_name}")

    return seq_dict


def write_fastx(
    file_name: str | Path, seq_dict: dict, chars_per_line: int = None, append: bool = False, reduced_name: bool = False
):
    file_name = Path(file_name)  # Ensure the path is Path object
    extension = file_name.suffix[1:]  # split('.')[-1]
    assert extension in ["fasta", "fa", "fastq", "fq"], f"Invalid extension for {file_name}"
    out_mode = "a" if extension in ["fasta", "fa"] else "q"

    write_mode = "a" if append else "w"
    file_name = Path(file_name)
    with open(file_name, write_mode) as fout:
        if out_mode == "q":  # FASTQ
            for name, r in seq_dict.items():
                fout.write(f"@{name}\n{r['Sequence']}\n+\n{r['Quality']}\n")
        else:  # FASTA
            for name, r in seq_dict.items():
                seq = r["Sequence"]
                if reduced_name:
                    name = name.split(" ")[0]
                if chars_per_line:
                    seq = "\n".join([seq[i : i + chars_per_line] for i in range(0, len(seq), chars_per_line)])
                fout.write(f">{name}\n{seq}\n")

    return None


def write_fasta(
    file_name: str | Path, seq_dict: dict, chars_per_line: int = None, append: bool = False, reduced_name: bool = False
):
    # DBW comment 11/18/23
    # The "mode" variable was replaced with "append" to be more independently descriptive and not rely on the syntax
    # of the open function. I had to change line 283 in Analysis.py that used this. Also, nchars has been replaced
    # with chars_per_line and is now also togglable, so that we have the option to produce fasta files with no line
    # breaks in the sequence. The 80 bases per line norm is a relic of when the FASTA format was originally defined
    # back in the ~1980s and terminal widths were fixed size, and I don't know of any modern software packages that
    # require fixed width data. It's a good feature to keep for flexibility and visualization. I also added a
    # "reduced_name" (often useful when dealing with Uniprot-style data)
    """
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
    """
    mode = "a" if append else "w"
    file_name = Path(file_name)
    with open(file_name, mode) as fout:
        for name, seq in seq_dict.items():
            if reduced_name:
                name = name.split(" ")[0]
            if chars_per_line:
                seq = "\n".join([seq[i : i + chars_per_line] for i in range(0, len(seq), chars_per_line)])
            fout.write(f">{name}\n{seq}\n")
    return None


def sam_to_fastx(sam_file: str | Path, output_file: str | Path):
    """
    Process sam file, return output
    """
    sam_file = Path(sam_file)
    output_file = Path(output_file)
    # Process SAM file
    sam_fields = ["QNAME", "FLAG", "RNAME", "POS", "MAPQ", "CIGAR", "RNEXT", "PNEXT", "TLEN", "SEQ", "QUAL"]
    with open(output_file, "w") as fout:
        mode = "q" if output_file.suffix in [".fq", ".fastq"] else "a"
        with open(sam_file) as fin:
            for line in fin:
                if line[0] == "@":
                    continue
                items = line.rstrip().split("\t")
                sam_dict = dict(zip(sam_fields, items[:11]))
                if mode == "a":
                    fout.write(f">{sam_dict['QNAME']}\n{sam_dict['SEQ']}\n")
                elif mode == "q":
                    fout.write(f"@{sam_dict['QNAME']}\n{sam_dict['SEQ']}\n+\n{sam_dict['QUAL']}\n")


def read_pod5(file_name: str | Path, out_file: str | Path = None):
    """
    Read a pod5 file to a polars dataframe
    """
    file_name = Path(file_name)
    records = {"rid": [], "sr": [], "size": [], "adc": [], "pa": []}
    with Reader(file_name) as r:
        for read in r.reads():
            records["rid"].append(str(read.read_id))  # name
            records["sr"].append(read.run_info.sample_rate)  # Hz
            records["size"].append(read.signal_pa.size)  # N
            records["adc"].append(read.signal)  # raw ADC counts (int16)
            records["pa"].append(read.signal_pa)  # applies per-read calibration, picoAmp

    df = pl.DataFrame(records)
    if out_file is not None:
        out_file = Path(out_file)
        df.write_parquet(out_file)
    return df
