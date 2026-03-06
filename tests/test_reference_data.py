"""Smoke tests for NanoporeAnalysis.reference_data."""

from NanoporeAnalysis.reference_data import CODON_TO_RESIDUE, CODONS, NUCLEOTIDES


def test_codon_table_complete():
    # 64 standard codons + gap entries
    assert len(CODONS) >= 64


def test_all_codons_map_to_residue():
    for codon in CODONS:
        assert codon in CODON_TO_RESIDUE


def test_nucleotides():
    assert {"A", "T", "C", "G"}.issubset(set(NUCLEOTIDES))


def test_stop_codons():
    stops = [c for c, r in CODON_TO_RESIDUE.items() if r == "."]
    assert set(stops) == {"TAA", "TAG", "TGA"}
