#!/usr/bin/env python
"""
Created on Apr 2, 2013
Revised on Apr 7, 2017 + include EASE2 support

@author: DGL
"""

from .loadsir import loadsir
from .ease2helper import ease2_map_info
import numpy as np
import sys
import argparse


def sirstring(a, cnt):
    """
    sirstring(a,cnt)
    decodes and returns a string starting at a that is cnt bytes long
    """
    string = np.zeros(2 * cnt, dtype="i1")
    for i in range(cnt):
        j = (i - 1) * 2 + 1
        string[j - 1] = np.mod(a[i - 1], 256)
        string[j] = np.floor(a[i - 1] / 256)
    res = str(bytearray(string).decode("utf-8")).rstrip(" \n\0\1\r")

    return res


def printsirhead(head):
    """
    printsirhead(head)

    prints information from the SIR file header info array from loadsir

    INPUTS:
        head:       scaled header information block
    """
    nhtype = int(head[4])
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
    nhead = int(head[40])
    if nhtype == 1:
        nhead = 1

    # ndes = int(head[41])
    # ldes = int(head[42])
    # nia = int(head[43])
    idatatype = int(head[47])

    iscale_sc = int(head[30])
    ixdeg_off = int(head[126])
    iydeg_off = int(head[127])
    ideg_sc = int(head[168])
    ia0_off = int(head[189])
    ib0_off = int(head[240])
    i0_sc = int(head[255])

    if iopt == -1:  # image only
        print("Image only form (%d)" % iopt)
        print("  x,y scale:   %f %f  %d" % (ascale, bscale, iscale_sc))
        print("  x,y offsets: %f %f  %d %d %d" % (a0, b0, ia0_off, ib0_off, i0_sc))
        print(
            "  x,y span:    %f %f  %d %d %d"
            % (xdeg, ydeg, ixdeg_off, iydeg_off, ideg_sc)
        )

    elif iopt == 0:  # rectalinear lat/lon
        print("Lat/Long Rectangular form (%d) ", iopt)
        print("  x,y scale: (pix/deg) %f %f  %d" % (ascale, bscale, iscale_sc))
        print(
            "  x,y offsets: (deg)   %f %f  %d %d %d" % (a0, b0, ia0_off, ib0_off, i0_sc)
        )
        print(
            "  x,y span: (deg)      %f %f  %d %d %d"
            % (xdeg, ydeg, ixdeg_off, iydeg_off, ideg_sc)
        )

    elif iopt == 1 or iopt == 2:  # lambert
        if iopt == 1:
            print("Lambert form (global radius) (%d)" % iopt)
        else:
            print("Lambert form (local radius) (%d)" % iopt)

        iascale = 1 / ascale
        ibscale = 1 / bscale
        print("  x,y scale: (km/pix) %f %f  %d" % (iascale, ibscale, iscale_sc))
        print(
            "  x,y corner: (km)    %f %f  %d %d %d" % (a0, b0, ia0_off, ib0_off, i0_sc)
        )
        print(
            "  x,y span: (deg)     %f %f  %d %d %d "
            % (xdeg, ydeg, ixdeg_off, iydeg_off, ideg_sc)
        )

        # create proj4 string for projection
        # semimajor_radius = 6378135.0
        f = 298.260
        # semiminor_radius = semimajor_radius * (1.0 - 1.0 / f)
        latitude_of_projection_origin = ydeg
        longitude_of_projection_origin = xdeg
        # latitude_of_true_scale = ydeg
        xgrid_min = a0 * 1000.0
        xgrid_max = (a0 + nsx / ascale) * 1000.0
        ygrid_min = b0 * 1000
        ygrid_max = (b0 + nsy / bscale) * 1000.0

        # calculate local radius used for projection following
        # that used in latlon2pix.py
        if iopt == 2:
            radearth = 6378.135  # equitorial earth radius
            dtr = 3.141592654 / 180.0
            # orglon = xdeg
            # orglon1 = np.mod(orglon + 720.0, 360.0)
            orglat = ydeg
            era = 1.0 - 1.0 / f
            eradearth = (
                radearth
                * era
                / np.sqrt(
                    era * era * np.cos(orglat * dtr) ** 2 + np.sin(orglat * dtr) ** 2
                )
            )

        proj4_string_1 = "+proj=laea +lat_0={} "
        proj4_string_1 = proj4_string_1.format(latitude_of_projection_origin)
        proj4_string_2 = "+lon_0={} +k=1 +x_0=0 +y_0=0 "
        proj4_string_2 = proj4_string_2.format(longitude_of_projection_origin)
        # proj4_string_3 = "+a={:.3f} +rf={:.3f} +units=m +no_defs "
        # proj4_string_3 = proj4_string_3.format(semimajor_radius, 50000000.)
        proj4_string_3 = "+R={:.3f} +units=m +no_defs "
        proj4_string_3 = proj4_string_3.format(eradearth * 1000)

        proj4_string_4 = "-a_ullr {:.3f} {:.3f} {:.3f} {:.3f}"
        proj4_string_4 = proj4_string_4.format(
            xgrid_min, ygrid_max, xgrid_max, ygrid_min
        )

        proj4_string = proj4_string_1 + proj4_string_2 + proj4_string_3

        print("  PROJ4 string:", proj4_string)
        print("  GeoTransform:", proj4_string_4)

    elif iopt == 5:  # polar stereographic
        print("Polar Stereographic form (%d)" % iopt)
        print("  x,y scale: (km/pix) %f %f  %d" % (ascale, bscale, iscale_sc))
        print(
            "  x,y offsets: (km)   %f %f  %d %d %d" % (a0, b0, ia0_off, ib0_off, i0_sc)
        )
        print(
            "  Center: (lon,lat)   %f %f  %d %d %d"
            % (xdeg, ydeg, ixdeg_off, iydeg_off, ideg_sc)
        )

    elif (iopt == 8) or (iopt == 9) or (iopt == 10):  # EASE2
        if iopt == 8:
            print("EASE2 form Northern Hemisphere(%d)" % iopt)
        elif iopt == 9:
            print("EASE2 form Southern Hemisphere (%d)" % iopt)
        else:
            print("EASE2 form Global (%d)" % iopt)
        (
            map_equatorial_radius_m,
            map_eccentricity,
            e2,
            map_reference_latitude,
            map_reference_longitude,
            map_second_reference_latitude,
            sin_phi1,
            cos_phi1,
            kz,
            map_scale,
            bcols,
            brows,
            r0,
            s0,
            epsilon,
        ) = ease2_map_info(iopt, int(ascale), int(bscale))
        print("  EASE2 scale: (%d %d) %f km/pix" % (ascale, bscale, map_scale * 0.001))
        print(
            "  Pixel origin: (c,r) %d %d   %d %d %d %d"
            % (a0, b0, ia0_off, ib0_off, i0_sc, iscale_sc)
        )
        print(
            "  Center pixel: (c,r) %d %d   %d %d %d"
            % (xdeg, ydeg, ixdeg_off, iydeg_off, ideg_sc)
        )

    elif (iopt == 11) or (iopt == 12) or (iopt == 13):  # EASE
        if iopt == 11:
            print("EASE form Northern Hemisphere(%d)" % iopt)
        elif iopt == 12:
            print("EASE form Southern Hemisphere (%d)" % iopt)
        else:
            print("EASE form Global (%d)" % iopt)
        print("  Scale: (c,r)  %d %d  %d" % (ascale, ascale, iscale_sc))
        print("  Origin: (c,r) %d %d  %d %d %d" % (a0, b0, ia0_off, ib0_off, i0_sc))
        print(
            "  Center: (c,r) %d %d  %d %d %d"
            % (xdeg, ydeg, ixdeg_off, iydeg_off, ideg_sc)
        )

    else:
        print("*** Unknown SIR transformation in file header: %d" % iopt)

    print("  Pixels:  %d by %d " % (nsx, nsy))
    print(
        "  Year: %d  Start Day: %d   Min:  %d  End Day: %d  Min: %d"
        % (
            int(head[11]),
            int(head[12]),
            int(head[13]),
            int(head[14]),
            int(head[16]),
        )
    )
    print(
        "  Region: %d  Type: %d  Htype: %d  Nhead: %d"
        % (int(head[17]), int(head[18]), int(head[4]), nhead)
    )
    print("  Nodata value: %f  Vmin: %f  Vmax: %f" % (head[48], head[49], head[50]))
    if idatatype == 4:
        print("  Datatype: float")
    elif idatatype == 1:
        print("  Datatype: byte  offset: %d  scale: %d" % (ioff, iscale))
    else:
        print("  Datatype: short  offset: %d  scale: %d" % (ioff, iscale))
    print("  Frequency: %.2f   Polarization: %d" % (head[45], int(head[44])))
    print("  Type: %s" % sirstring(head[57:79], 22))
    print("  Sensor: %s" % sirstring(head[19:39], 20))
    print("  Title: %s" % sirstring(head[128:168], 40))
    print("  Tag: %s" % sirstring(head[169:189], 20))
    print("  Creator: %s" % sirstring(head[190:240], 50))
    print("  Created: %s" % sirstring(head[241:255], 14))


def main():

    parser = argparse.ArgumentParser(description="Print summary of SIR header values")

    # positional arguments
    parser.add_argument("filename", help="SIR filename")

    # get arguments
    args = parser.parse_args()

    sir_fname = args.filename

    # A test driver
    sir = loadsir(sir_fname)
    print("File name: %s" % sir_fname)
    printsirhead(sir[1])


if __name__ == "__main__":

    main()
