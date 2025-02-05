'''
Basic utility functions
'''

import time, datetime
import pyarrow as pa

def timer( start, ):
    return str( datetime.timedelta( seconds=round( time.time() - start ) ) )

def reverse_complement(sequence):
    """Return the reverse complement of a DNA nucleotide string"""
    rc_trans = str.maketrans('ATCGatcg', 'TAGCtagc')
    seq_trans = sequence.translate(rc_trans)
    rc = seq_trans[::-1]
    return rc

def fill_table_nulls(table, fill_value) :
    table = pa.table([x.fill_null( pa.array([fill_value]).cast(x.type).to_pylist()[0] ) for x in table.itercolumns()], names=table.column_names)
    return table