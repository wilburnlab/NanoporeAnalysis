#===============================================================================
# gff2bed.py

"""
Copyright (c) 2022 Anthony Aylward

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

#===============================================================================

"""Convert GFF3 data to BED format"""

# Imports ======================================================================

import gzip
from pathlib import PosixPath




# Functions ====================================================================

def parse_gff_attributes(attr: str, graceful: bool = False):
    """Parse an entry from the "attr" column of a GFF3 file and return it as
    a dict

    Parameters
    ----------
    attr : str
        feature attribute string
    graceful : bool
        if False, throw an error when a malformed tag-value pair is encountered
        if True, ignore malformed pairs gracefully

    Returns
    ------
    dict
        attr entries as a dict
    """

    return dict(pair.split('=') for pair in attr.split(';')
                if ('=' in pair) or not graceful)


def parse(gff_file, type: str = 'gene', parse_attr: bool = True,
          graceful: bool = False):
    """Parse a GFF3 file and yield its lines as tuples

    Parameters
    ----------
    gff_file
        String, PosixPath, or file-like object representing a GFF3 file
    type
        string indicating feature type to include, or None to include all
        features [gene]
    parse_attr : bool
        if False, do not parse attributes [True]
    graceful : bool
        if False, throw an error when a malformed tag-value pair is encountered
        if True, ignore malformed pairs gracefully [False]
    
    Yields
    ------
    seqid, start, end, strand, attr
        coordinates of a feature
    """

    with (
        gff_file if not isinstance(gff_file, (str, PosixPath))
        else gzip.open(gff_file, 'rt') if str(gff_file).endswith('.gz')
        else open(gff_file, 'r')
    ) as f:
        for line in f:
             if not line.startswith('#'):
                seqid, _, t, start, end, _, strand, _, attr = line.rstrip().split('\t')
                if ((t == type) or (type is None)):
                    if parse_attr:
                        yield (seqid, int(start), int(end), strand,
                            parse_gff_attributes(attr, graceful=graceful))
                    else:
                        yield seqid, int(start), int(end), strand, '.'


def convert(gff_data, tag: str = 'ID'):
    """Convert rows of GFF3 data to BED data. This involves reordering the
    columns to conform with BED format and shifting the coordinates to 0-based
    half-open values.

    Parameters
    ----------
    gff_data
        iterable of data from gff2bed.parse()
    tag : str
        GFF3 attribute tag to parse [ID]

    Yields
    ------
    tuple
        a row of BED data
    """

    for seqid, start, end, strand, attr in gff_data:
        yield seqid, start - 1, end, attr[tag], 0, strand
