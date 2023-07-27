'''
Basic utility functions
'''

import time, datetime


def timer( start, ):
    return str( datetime.timedelta( seconds=round( time.time() - start ) ) )
