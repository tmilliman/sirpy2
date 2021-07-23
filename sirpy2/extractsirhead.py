#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Created on Dec 10, 2018

@author: TEM
"""
from __future__ import division
from .loadsir import loadsir
from .ease2helper import ease2_map_info
from .printsirhead import printsirhead
import numpy as np
import sys

REGIONS = {
    1: "Ama",
    2: "NAm",
    202: "Grn",
    203: "Ala",
    204: "CAm",
    205: "NAm",
    206: "SAm",
    207: "NAf",
    208: "SAf",
    209: "Sib",
    210: "Eur",
    211: "SAs",
    212: "ChJ",
    213: "Ind",
    214: "Aus",
    256: "Ber",
}


def sirstring(a, cnt):
    """sirstring(a,cnt)
    decodes and returns a string starting at a that is cnt bytes long
    """
    string = np.zeros(2 * cnt, dtype="i1")
    for i in range(cnt):
        j = (i - 1) * 2 + 1
        string[j - 1] = np.mod(a[i - 1], 256)
        string[j] = np.floor(a[i - 1] / 256)
    res = str(bytearray(string).decode("utf-8")).rstrip(" \n\0\1\r")

    return res


def extractsirhead(head):
    """
    printsirhead(head)

    extract information from the SIR file header info array from
    loadsir and return a dictionary of values

    INPUTS:
        head:       scaled header information block
    """

    region = int(head[17])
    region_type = int(head[18])
    nsx = int(head[0])
    nsy = int(head[1])
    iopt = int(head[16])
    xdeg = head[2]
    ydeg = head[3]
    ascale = head[5]
    bscale = head[6]
    a0 = head[7]
    b0 = head[8]

    ioff = int(head[9])
    iscale = int(head[10])
    year = int(head[11])
    doy_start = int(head[12])
    doy_end = int(head[13])
    nhead = int(head[40])
    nhtype = int(head[4])
    if nhtype == 1:
        nhead = 1

    ndes = int(head[41])
    ldes = int(head[42])
    nia = int(head[43])
    idatatype = int(head[47])
    nodata = head[48]

    iscale_sc = int(head[30])
    ixdeg_off = int(head[126])
    iydeg_off = int(head[127])
    ideg_sc = int(head[168])
    ia0_off = int(head[189])
    ib0_off = int(head[240])
    i0_sc = int(head[255])

    proj4_string = None
    if iopt == 2:
        # create proj4 string for projection
        semimajor_radius = 6378135.0
        f = 298.260
        semiminor_radius = semimajor_radius * (1.0 - 1.0 / f)
        latitude_of_projection_origin = ydeg
        longitude_of_projection_origin = xdeg
        # latitude_of_true_scale = ydeg

        # calculate local radius used for projection following
        # code used in latlon2pix.py
        eradearth = semimajor_radius
        radearth = 6378.135  # equitorial earth radius
        dtr = 3.141592654 / 180.0
        # orglon = xdeg
        # orglon1 = np.mod(orglon + 720.0, 360.0)
        orglat = ydeg
        era = 1.0 - 1.0 / f
        eradearth = (
            radearth
            * era
            / np.sqrt(era * era * np.cos(orglat * dtr) ** 2 + np.sin(orglat * dtr) ** 2)
        )

        iascale = 1.0 / ascale
        ibscale = 1.0 / bscale
        xgrid_min = a0 * 1000.0
        xgrid_max = (a0 + nsx / ascale) * 1000.0
        ygrid_min = b0 * 1000
        ygrid_max = (b0 + nsy / bscale) * 1000.0
        proj4_string_1 = "+proj=laea +lat_0={} "
        proj4_string_1 = proj4_string_1.format(latitude_of_projection_origin)
        proj4_string_2 = "+lon_0={} +k=1 +x_0=0 +y_0=0 "
        proj4_string_2 = proj4_string_2.format(longitude_of_projection_origin)
        # proj4_string_3 = "+a={:.3f} +rf={:.3f} +units=m +no_defs "
        # proj4_string_3 = proj4_string_3.format(semimajor_radius, 500000000.)
        proj4_string_3 = "+R={:.3f} +units=m +no_defs "
        proj4_string_3 = proj4_string_3.format(eradearth * 1000)
        proj4_string_4 = "+a_ullr {:.3f} {:.3f} {:.3f} {:.3f}"
        proj4_string_4 = proj4_string_4.format(
            xgrid_min, ygrid_max, xgrid_max, ygrid_min
        )
        proj4_string = proj4_string_1 + proj4_string_2 + proj4_string_3

        geoTransform = proj4_string_4

    sensor = sirstring(head[19:39], 20)
    datatype = sirstring(head[57:79], 22)

    # make a dictionary of extracted values
    header_dict = {
        "region": REGIONS[region],
        "region_type": region_type,
        "year": year,
        "doy_start": doy_start,
        "doy_end": doy_end,
        "nhtype": nhtype,
        "nsx": nsx,
        "nsy": nsy,
        "iopt": iopt,
        "xdeg": xdeg,
        "ydeg": ydeg,
        "ascale": ascale,
        "bscale": bscale,
        "iascale": iascale,
        "ibscale": ibscale,
        "a0": a0,
        "b0": b0,
        "semimajor_radius": semimajor_radius,
        "semiminor_radius": semiminor_radius,
        "ioff": ioff,
        "iscale": iscale,
        "nhead": nhead,
        "ndes": ndes,
        "ldes": ldes,
        "nia": nia,
        "nodata": nodata,
        "iscale_sc": iscale_sc,
        "idatatype": idatatype,
        "iscale": iscale,
        "ixdeg_off": ixdeg_off,
        "iydeg_off": iydeg_off,
        "ideg": ideg_sc,
        "ia0_off": ia0_off,
        "ib0_off": ib0_off,
        "i0_sc": i0_sc,
        "datatype": datatype,
        "sensor": sensor,
        "grid_dx": iascale * 1000,
        "grid_dy": ibscale * 1000,
        "xgrid_min": xgrid_min,
        "xgrid_max": xgrid_max,
        "ygrid_min": ygrid_min,
        "ygrid_max": ygrid_max,
        "proj4_string": proj4_string,
        "geotransform": geoTransform,
    }

    return header_dict


def main(sir_fname):
    # A test driver
    sir = loadsir(sir_fname)
    print("File name: %s" % sir_fname)
    printsirhead(sir[1])


if __name__ == "__main__":
    main(sys.argv[1])
