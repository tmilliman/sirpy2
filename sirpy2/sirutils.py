#!/usr/bin/env python3

# script to convert BYU .sir files to geotiffs

import re
import datetime as dt
import numpy as np


def parseFilename(filename, verbose=False):

    """
    Parse BYU filenames and return dictionary.
    Filenames are of the format:

    msfa-a-NAm18-267-271.sir
    quev-a-NAm99-201-204.sir

    Data Product

    ascat
    =====
    msfa = "all passes"
    mafa = "ascending passes only"
    mdfa = "descending passes only"
    mmfa = "morning passes only"
    mnfa = "mid-day passes only"
    mefa = "evening passes only"
    
    qscat
    =====
    quev = "all passes, egg reconstruction, v-polarization


    Image classes:

    a = sigma-0 value in dB normalized to 40 degree (ref.) incidence
        angle

    V = standard deviation of reconstruction error in dB

    Region:

    three-letter code

    Date

    Year-DOY start-DOY end

    """

    fname_patt = "(?P<product>)-(?P<itype>[aV])"
    fname_patt = fname_patt + "-(?P<region>\w\w\w)(?P<year>\d\d)"
    fname_patt = fname_patt + "-(?P<doy_start>\d\d\d)"
    fname_patt = fname_patt + "-(?P<doy_end>\d\d\d)\."

    m = re.match(fname_patt, filename)
    if m is not None:
        year2 = int(m.group('year'))
        if (year2 > 80):
            year4 = year2 + 1900
        else:
            year4 = 2000 + year2
        fparts = {'product': m.group('product'),
                  'itype': m.group('itype'),
                  'region': m.group('region'),
                  'year': year4,
                  'doy_start': int(m.group('doy_start')),
                  'doy_end': int(m.group('doy_end'))}
        return fparts
    else:
        return None


def fn2dt(filename):
    """
    Parse filename and return a datetime object for the middle
    of the doy range.
    """

    fndict = parseFilename(filename)
    if fndict is not None:
        center_doy = np.mean([float(fndict['doy_start']),
                             float(fndict['doy_end'])])
        fndt = dt.datetime(fndict['year'], 1, 1, 0, 0, 0) + \
            dt.timedelta(days=center_doy-1)
    else:
        fndt = None

    return fndt
    
if __name__ == '__main__':

    """
    If invoked on the command line run tests."
    """

    filename = 'msfa-a-NAm07-001-005.sir'
    fparts = parseFilename(filename)

    print(fparts)
