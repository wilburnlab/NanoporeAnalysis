'''
Basic utility functions
'''

import time, datetime


def timer( start, ):
    return str( datetime.timedelta( seconds=round( time.time() - start ) ) )

def reverse_complement(sequence):
    """Return the reverse complement of a DNA nucleotide string"""
    rc_trans = str.maketrans('ATCGatcg', 'TAGCtagc')
    seq_trans = sequence.translate(rc_trans)
    rc = seq_trans[::-1]
    return rc