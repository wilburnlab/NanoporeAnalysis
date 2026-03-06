"""Smoke tests for NanoporeAnalysis.utils."""

from NanoporeAnalysis.utils import (
    decode_phred,
    encode_phred,
    identify_alphabet,
    reverse_complement,
    translate,
)


def test_reverse_complement():
    assert reverse_complement("ATCG") == "CGAT"
    assert reverse_complement("") == ""


def test_phred_roundtrip():
    scores = [0, 10, 20, 30, 40]
    encoded = encode_phred(scores)
    decoded = decode_phred(encoded)
    assert list(decoded) == scores


def test_identify_alphabet_dna():
    assert identify_alphabet("ATCGATCG") == "DNA"


def test_identify_alphabet_protein():
    assert identify_alphabet("MVLSPADKTNVK") == "Protein"


def test_translate_basic():
    # ATG = M, TAA = stop (.)
    assert translate("ATGTAA") == "M."
