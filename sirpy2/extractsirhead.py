#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Created on Dec 10, 2018

@author: TEM
"""
from __future__ import division
from .loadsir import loadsir
from .ease2helper import ease2_map_info
import numpy as np
import sys

REGIONS = {
    1: 'Ama',
    204: 'CAm',
    205: 'NAm',
    206: 'SAm',
    207: 'NAf',
    208: 'SAf',
    209: 'Sib',
    210: 'Eur',
    211: 'SAs',
    212: 'ChJ',
    213: 'Ind',
    214: 'Aus'}


def sirstring(a,cnt):
    """sirstring(a,cnt)
    decodes and returns a string starting at a that is cnt bytes long
    """
    string=np.zeros(2*cnt,dtype='i1')
    for i in range(cnt):
        j=(i-1)*2+1;
        string[j-1]=(np.mod(a[i-1],256))
        string[j]=np.floor(a[i-1]/256)
    res=(str(bytearray(string).decode("utf-8")).rstrip(' \n\0\1\r'))

    return res
    
def extractsirhead(head):
    """
    printsirhead(head)

    extract information from the SIR file header info array from
    loadsir and return a dictionary of values

    INPUTS:
        head:       scaled header information block
    """

    region = np.int(head[17])
    region_type = np.int(head[18])
    nsx    = np.int(head[0])
    nsy    = np.int(head[1])
    iopt   = np.int(head[16])
    xdeg   = head[2]
    ydeg   = head[3]
    ascale = head[5]
    bscale = head[6]
    a0     = head[7]
    b0     = head[8]

    ioff=np.int(head[9])
    iscale=np.int(head[10])
    year = np.int(head[11])
    doy_start = np.int(head[12])
    doy_end = np.int(head[13])
    nhead = np.int(head[40])
    nhtype = np.int(head[4])
    if nhtype == 1:
        nhead = 1

            

    ndes = np.int(head[41])
    ldes = np.int(head[42])
    nia = np.int(head[43])
    idatatype=np.int(head[47])
    
    iscale_sc = np.int(head[30])
    ixdeg_off = np.int(head[126])
    iydeg_off = np.int(head[127])
    ideg_sc   = np.int(head[168])
    ia0_off   = np.int(head[189])
    ib0_off   = np.int(head[240])
    i0_sc     = np.int(head[255])

    proj4_string = None
    if iopt == 2:
        # create proj4 string for projection
        semimajor_radius = 6378135.0
        f = 298.260
        semiminor_radius = semimajor_radius * (1.0 - 1.0/f)
        latitude_of_projection_origin = ydeg
        longitude_of_projection_origin = xdeg
        latitude_of_true_scale = ydeg
        iascale = 1./ascale
        ibscale = 1./bscale
        xgrid_min = a0 * 1000.
        xgrid_max = (a0 + nsx/ascale) * 1000.
        ygrid_min = b0 * 1000
        ygrid_max = (b0 + nsy/bscale) * 1000.
        proj4_string_1 = "datum=wgs84 +proj=laea +lat_0={} "
        proj4_string_1 = proj4_string_1.format(latitude_of_projection_origin)
        proj4_string_2 = "+lon_0={} +k=1 +x_0=0 +y_0=0 "
        proj4_string_2 = proj4_string_2.format(longitude_of_projection_origin)
        proj4_string_3 = "+a={:.3f} +b={:.3f} +rf={:.3f} +units=m +no_defs "
        proj4_string_3 = proj4_string_3.format(semimajor_radius,
                                               semiminor_radius, 50000000.)
        proj4_string_4 = "+a_ullr {:.3f} {:.3f} {:.3f} {:.3f}"
        proj4_string_4 = proj4_string_4.format(xgrid_min, ygrid_min,
                                               xgrid_max, ygrid_max)
        proj4_string = proj4_string_1 + proj4_string_2 + \
            proj4_string_3 + proj4_string_4

    sensor = sirstring(head[19:39], 20)
    datatype = sirstring(head[57:79], 22)

    # make a dictionary of extracted values
    header_dict = {
        'region': REGIONS[region],
        'region_type': region_type,
        'year': year,
        'doy_start': doy_start,
        'doy_end': doy_end,
        'nhtype': nhtype,
        'nsx': nsx,
        'nsy': nsy,
        'iopt': iopt,
        'xdeg': xdeg,
        'ydeg': ydeg,
        'ascale': ascale,
        'bscale': bscale,
        'iascale': iascale,
        'ibscale': ibscale,
        'a0': a0,
        'b0': b0,
        'ioff': ioff,
        'iscale': iscale,
        'nhead': nhead,
        'ndes': ndes,
        'ldes': ldes,
        'nia': nia,
        'idatatype': idatatype,
        'iscale': iscale,
        'ixdeg_off': ixdeg_off,
        'iydeg_off': iydeg_off,
        'ideg': ideg_sc,
        'ia0_off': ia0_off,
        'ib0_off': ib0_off,
        'i0_sc': i0_sc,
        'datatype': datatype,
        'sensor': sensor,
        'grid_dx': iascale * 1000,
        'grid_dy': ibscale * 1000,
        'xgrid_min': xgrid_min,
        'xgrid_max': xgrid_max,
        'ygrid_min': ygrid_min,
        'ygrid_max': ygrid_max,
        'proj4_string': proj4_string}

    return header_dict
    
def main(sir_fname):
    # A test driver
    sir = loadsir(sir_fname)
    print("File name: %s" % sir_fname)
    printsirhead(sir[1])

if __name__ == "__main__":
    main(sys.argv[1])



