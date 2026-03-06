# NanoporeAnalysis

## Purpose and Scientific Domain

NanoporeAnalysis is a Python library for processing Oxford Nanopore Technologies (ONT) long-read sequencing data, specifically designed for **cDNA sequencing workflows**. The pipeline covers:

1. Raw signal handling (POD5 files)
2. Read annotation: polyA-tail detection, adapter/SSP identification, UMI scoring, and read orientation
3. ORF prediction and protein extraction from transcript sequences
4. Parallelized batch processing via Nextflow for HPC/cluster environments

The library targets lab-specific cDNA library prep protocols that use custom barcodes, a strand-switch primer (SSP), and UMIs, where reads must be classified by identity and orientation, and translated into predicted protein sequences.

## Module Structure

```
NanoporeAnalysis/
  __init__.py          # Empty
  align.py             # Sequence alignment and read annotation (edlib-based; legacy skbio code at bottom)
  local_io.py          # FASTA/FASTQ/SAM I/O with gzip support, POD5 reading
  pod5.py              # POD5 signal processing and Parquet conversion
  reference_data.py    # Biological constants (codon table, nucleotides, residues)
  utils.py             # Sequence utilities: reverse complement, translation, ORF search, phred, counting
  Analysis.py          # LEGACY -- student code, slated for removal. Do not build on this.
  gff2bed.py           # LEGACY -- vendored GFF3-to-BED converter, slated for removal.

Reads_to_Proteins.py   # Standalone CLI for ORF prediction from FASTX files (multiprocessing)

transcriptome_nf/      # Nextflow (DSL2) pipeline for parallelized read processing
  main.nf              # Pipeline definition: batch -> process -> extract_proteins
  nextflow.config      # Execution profiles (local, SLURM) and parameters
  scripts/
    batch_sam.py       # Split SAM files into batches for parallel processing
    process_sam.py     # Per-batch read annotation via align.process_read(), outputs Parquet
    extract_proteins.py # Aggregate Parquet files, filter complete cDNAs, produce protein FASTA + count matrix
    merge_parquets.py  # Parquet merge utility (currently unused)
    primers.txt        # Primer/UMI configuration file
  Dockerfile           # Container definition
  environment.yml      # Conda environment spec
```

### align.py (the active alignment engine)
Contains the **edlib-based alignment layer** (the intended path forward) and legacy skbio code:

**Active edlib-based functions** (top of file):
- `align()` -- Edit-distance alignment via edlib (HW mode, semi-global)
- `merge_overlapped_indices()` -- Combines overlapping alignment location pairs
- `find_polyX()` / `find_best_polyA()` -- Detect homopolymer runs (polyA tails, min 15 nt, max 40 nt). Sorted by position (last occurrence first) rather than length.
- `split_sequence_by_polyA()` -- Split reads at the polyA boundary, returns pre-polyA, post-polyA, and polyA stats dict
- `parse_adapter5_seq()` -- Find 5' SSP (max 2 edits) and extract transcript sequence with trim index
- `parse_adapter3_seq()` -- Identify 3' adapter features including UMIs (max 2 edits)
- `score_umi()` -- Score putative UMI sequences against a known UMI dictionary. Returns edit distances plus Best UMI / Best UMI Score / Delta UMI Score.
- `compute_read_statistics()` -- Read length and mean quality score
- `annotate_read()` -- Full read annotation: polyA, SSP, UMI, adapter sequences, transcript extraction, ORF prediction, cDNA status classification (Complete / 3' Only / 5' Only / Missing Barcode / Unknown, with Fragment suffix for short transcripts). Constructs a Sample_ID label from source + best UMI.
- `process_read()` -- Analyzes both forward and reverse-complement orientations, selects the best based on cDNA status classification. Also encodes movemap data if present in SAM tags.

**Legacy skbio SSW functions** (bottom of file, to be removed):
- `align_to_dict()`, `reverse_alignment()`, `get_best_barcode()`, `get_best_alignment()`, `show_alignment_simple()` -- Used by the legacy Analysis.py. These depend on `skbio` and `cigar`. The skbio imports are already commented out.

### local_io.py
- `read_fastx()` -- Reads FASTA or FASTQ (auto-detected, gzip-aware), returns dict keyed by read ID. All entries use a consistent format: `{'Sequence': str, 'Quality': str|None, 'Tags': list|None}`.
- `write_fastx()` -- Writes FASTA or FASTQ based on file extension, with optional line wrapping.
- `write_fasta()` -- Original FASTA-only writer (still present alongside `write_fastx`).
- `sam_to_fastx()` -- Convert SAM file to FASTA or FASTQ.
- `read_pod5()` -- Read POD5 file into a Polars DataFrame with signal data.
- `FileContentsError` -- Custom exception for empty files.

### pod5.py
Handles raw ONT signal data (POD5 format) with conversion to Parquet for efficient storage and querying:
- `pod5_to_parquet()` -- Stream POD5 into Parquet with configurable batching, compression (zstd), and optional file splitting. Stores read metadata + raw signal in a single table.
- `return_pod5_schema()` -- Configurable PyArrow schema for POD5 data (UUID encoding options, calibration modes).
- `read_to_record()` / `records_to_table()` -- Record-level processing helpers.
- `read_calibration()` -- Robust calibration extraction supporting multiple POD5 API versions (scale/shift or range/digitisation/offset).
- `read_pod5()` -- Simple POD5 reader to Polars DataFrame (also in local_io.py).
- Legacy `split_pod5()` via h5py at bottom of file.

### reference_data.py
Biological constants used by the translation machinery:
- `NUCLEOTIDES`, `RESIDUES` -- Valid character sets for DNA and protein sequences
- `CODONS`, `CODON_TO_RESIDUE` -- Standard genetic code translation table (stop codons as `.`)

### utils.py
Expanded beyond basic helpers to include core bioinformatics utilities:
- `reverse_complement()` -- DNA reverse complement
- `timer()` -- Elapsed time formatting
- `decode_phred()` / `encode_phred()` -- Phred quality score conversion (ASCII offset 33)
- `identify_alphabet()` -- Classify sequence as DNA, Protein, or Unknown
- `translate()` -- DNA to protein translation using the standard genetic code
- `orf_searcher()` -- Find all ORFs in a DNA sequence using regex (ATG to stop codon). Supports both strands, minimum length filtering, and pruning of internal ORFs sharing the same stop codon.
- `return_best_orf()` -- Wrapper returning the longest ORF
- `orf_check()` / `len_check()` -- Validation helpers
- `return_count_dict()` -- Count occurrences of a field across a list/dict of records

### transcriptome_nf/ (Nextflow pipeline)
A Nextflow (DSL2) pipeline that parallelizes the read-processing workflow for transcriptomics data. This is the intended deployment mechanism for the core Python library on cluster/HPC environments. Three stages:

1. **`generate_batches`** -- Splits input SAM files into batches of N reads (default 100,000) via `scripts/batch_sam.py`. Enables parallel processing of large sequencing runs.
2. **`process_reads`** -- Calls `scripts/process_sam.py` on each batch, which parses SAM records, loads primers from a config file, and runs `NanoporeAnalysis.align.process_read()` on every read. Outputs per-batch Parquet files.
3. **`extract_proteins`** -- Collects all Parquet files, filters for complete cDNAs with valid proteins, generates a protein FASTA and a Protein x Sample count matrix.

Supports local execution (36 CPUs / 64 GB default) and SLURM profiles. Dorado basecalling process is stubbed out (commented). Intermediate data uses Parquet format (via PyArrow).

### Reads_to_Proteins.py (standalone CLI)
Command-line tool for ORF prediction from FASTX files with multiprocessing support (`ProcessPoolExecutor`). Reads sequences, finds best ORF per read, and outputs per-ORF and per-protein count dictionaries as pickle files.

### Analysis.py (LEGACY -- do not extend)
Contains a student's pipeline functions (Dorado basecalling, debarcoding via skbio SSW, minimap2 mapping, read counting, DESeq2 differential expression, QC). This file is **not the intended path forward** and should eventually be removed. It has known bugs and stale imports.

### gff2bed.py (LEGACY -- do not extend)
Vendored GFF3-to-BED converter (MIT license, Anthony Aylward). Slated for removal as part of dependency minimization.

## Core Architectural Decisions

- **Minimize dependencies**: Prefer lightweight, well-maintained libraries. `edlib` is the chosen alignment backend. `skbio` is considered poorly designed and should not be used in new code. The `cigar`, `bioinfokit`, `pydeseq2`, and `scipy` imports are legacy dependencies.
- **edlib for alignment**: All alignment work uses `edlib` (edit-distance based, HW mode for semi-global alignment). The `align()` wrapper in align.py is the standard entry point. Scores are reported as raw edit distances (not normalized).
- **Read-as-dict model**: Reads are represented as dicts that get progressively enriched with annotation keys via dict merging (`read | annotations`). SAM fields, alignment scores, polyA stats, UMI scores, ORF predictions, and cDNA status all accumulate on the same dict.
- **Parquet as intermediate format**: The active pipeline (Nextflow) uses PyArrow Parquet files for intermediate and output data, replacing the CSV/pickle approach of the legacy code.
- **Nextflow for parallelism**: Batch processing is handled by Nextflow, which splits SAM files and fans out `process_read()` calls across batches. The Python library itself stays single-threaded per invocation.
- **POD5 signal handling**: Raw nanopore signal is read via the `pod5` library and can be converted to Parquet with configurable schemas and compression for downstream analysis.

## What Is Complete vs. Incomplete

### Complete / functional
- edlib-based alignment engine (align, polyA detection, adapter parsing, UMI scoring)
- Read annotation and orientation pipeline (`annotate_read`, `process_read`) with cDNA status classification
- ORF prediction and protein extraction (`orf_searcher`, `return_best_orf`, integrated into `annotate_read`)
- FASTA/FASTQ I/O (read and write, gzip-aware, consistent dict format)
- SAM-to-FASTX conversion
- POD5 reading and POD5-to-Parquet conversion
- Nextflow pipeline for batch processing (SAM -> annotated Parquet -> protein FASTA + count matrix)
- Standalone CLI for reads-to-proteins (`Reads_to_Proteins.py`)
- Biological reference data (codon table, translation)

### Incomplete / in progress
- **Dorado basecalling integration** -- stubbed out in the Nextflow config (commented process), no active non-legacy basecalling wrapper
- **UMI deduplication** -- scoring exists but no dedup or UMI-aware counting step
- **Differential expression** -- no current non-legacy DE implementation
- **QC and visualization** -- no current non-legacy QC
- **Alignment tuning** -- `max_score` defaults (2 edits for SSP/UMI/primer) and polyA length bounds (15-40) are functional but may need adjustment per experiment

## Observed Conventions

- **Path handling**: Mix of string concatenation and `pathlib.Path`. No consistent standard yet.
- **Naming**: Functions use `snake_case`. Module files use `snake_case` (except legacy `Analysis.py`).
- **Docstrings**: Present on most functions, loose Google-style format.
- **Alignment scores**: edlib layer reports raw edit distances (lower = better). `None` indicates no alignment found (replaces the old `0.5/len` floor).
- **Data format**: Active code uses dicts for reads, Parquet for persistence. Legacy code used CSV and pickle.
- **License**: Apache 2.0.

## Testing Conventions and Tooling

- **Test framework**: pytest, configured in `pyproject.toml` under `[tool.pytest.ini_options]`.
- **Test location**: `tests/` directory. Test files named `test_<module>.py`.
- **Running tests**: `pytest` from the project root. Use `pytest -v` for verbose output.
- **Linter**: ruff, configured in `pyproject.toml` under `[tool.ruff]`. Rules: E, F, I, UP, B, SIM. Line length 120. Target Python 3.9.
- **Run linter**: `ruff check NanoporeAnalysis/ tests/` (add `--fix` to auto-fix).
- **Pre-commit**: `.pre-commit-config.yaml` runs ruff lint + format checks on commit via `ruff-pre-commit`. Install with `pre-commit install`.
- **CI**: GitHub Actions (`.github/workflows/ci.yml`) runs ruff lint and pytest on Python 3.9/3.11/3.12 for pushes and PRs to `main`.
- **Package install**: `pip install -e ".[dev]"` installs the package in editable mode with dev dependencies (pytest, ruff, pre-commit).
- **Dev environment**: `environment.yml` at project root defines a conda environment (`nanopore-dev`) with all runtime and dev dependencies.

### Design decisions
- **ruff over flake8/black/isort**: Single tool replaces linter + formatter + import sorter. Fast, minimal config. The `ruff-pre-commit` hook is used rather than running ruff via a `local` hook, so pre-commit manages its own ruff version independently.
- **pyproject.toml as single config source**: All project metadata, build config, pytest settings, and ruff settings live in `pyproject.toml`. No `setup.py`, `setup.cfg`, `tox.ini`, or separate `.flake8`/`.isort.cfg` files.
- **Ruff rule selection**: E (pycodestyle), F (pyflakes), I (isort), UP (pyupgrade), B (bugbear), SIM (simplify). E501 (line length) is ignored since the formatter handles wrapping. These rules catch real bugs and style issues without being noisy.
- **CI matrix**: Tests run on Python 3.9 (minimum supported), 3.11, and 3.12. Linting runs once on latest Python. CI installs via `pip install -e ".[dev]"` rather than conda to keep the workflow simple and fast.
- **regex as a runtime dependency**: Added to `pyproject.toml` since `utils.py` imports it (used by `orf_searcher` for possessive quantifiers not available in the stdlib `re` module).
