#!/usr/bin/env python3
# -*- coding:utf-8 -*-

"""
script to convert BYU-MERS .sir file to geotiff.  Here we
preserve data values (no scaling to 0-255) like the c routine
supplied by the BYU-MERS sir2geotiff.c routine.
"""

import os
import sys
import argparse
import numpy as np
from osgeo import gdal, osr
from pathlib import Path

import sirpy2 as sp2


def main():

    parser = argparse.ArgumentParser(description="convert .sir file to geoTIFF file")

    # options
    parser.add_argument(
        "-v",
        "--verbose",
        help="increase output verbosity",
        action="store_true",
        default=False,
    )

    parser.add_argument(
        "-d",
        "--datadir",
        nargs="?",
        help=("data directory for output and finding landmask files"),
        const="./",
        default="./",
    )

    # positional args
    parser.add_argument("infile", help="Input .sir file")

    args = parser.parse_args()
    verbose = args.verbose
    infile = args.infile
    datadir = args.datadir

    if verbose:
        print("verbose: {}".format(verbose))
        print("input file: {}".format(infile))
        print("data directory: {}".format(datadir))

    # load sir file
    [image, head, descrip, iaopt] = sp2.loadsir(infile)
    if verbose:
        sp2.printsirhead(head)

    # extract the header info
    hd = sp2.extractsirhead(head)

    proj4_string = hd["proj4_string"]
    if verbose:
        print("proj4 string: {}".format(proj4_string))

    sig0_max = image.max()
    sig0_min = image.min()

    if verbose:
        print("Max Value: {}".format(sig0_max))
        print("Min Value: {}".format(sig0_min))

    # BUG fix: It looks like for the ERS-1,2 files the header
    # information is slightly different.  The result is the extracted
    # region for these files is always 100 (NAm).  So if sensor
    # is "ERS-1" or "ERS-2" then extract the region from the
    # filename

    # read land mask for this region
    sensor = hd["sensor"]
    datatype = hd["datatype"]
    if sensor.startswith("ERS"):
        fparts = sp2.parseFilename(os.path.basename(infile))
        region = fparts["region"]
    elif sensor == "QuikScat L1B":
        fparts = sp2.parseFilename(os.path.basename(infile))
        product = fparts["product"]
        region = fparts["region"]
    else:
        region = hd["region"]

    year = hd["year"]
    nodata_value = np.round(hd["nodata"])
    if verbose:
        print("Sensor: {}  Data Type: {}".format(sensor, datatype))
        print("Data type: {}".format(datatype))
        print("Region: {}".format(region))
        print("No Data Value: {}".format(nodata_value))

    # set name of landmask file based on sensor from file header
    if sensor.startswith("ASCAT"):
        landmask_file = "msf-{}.sir.lmask".format(region)
    elif sensor.startswith("OSCAT"):
        landmask_file = "oue-{}.sir.lmask".format(region)
    elif sensor.startswith("QuikScat L1B"):
        if product[2] == "e":
            landmask_file = "que-{}.sir.lmask".format(region)
        elif product[2] == "s":
            landmask_file = "qus-{}.sir.lmask".format(region)
    elif sensor.startswith("ERS-1/2"):
        landmask_file = "ers-{}.sir.lmask".format(region)
    elif sensor.startswith("NSCAT"):
        landmask_file = "nsc-{}.sir.lmask".format(region)
    elif sensor.startswith("SASS"):
        landmask_file = "sas-{}.sir.lmask".format(region)
    else:
        errmsg = "Unknown sensor"
        raise ValueError(errmsg)

    landmask_dir = os.path.join(datadir, "info", "landmasks")
    if verbose:
        print("Data directory: {}".format(datadir))
        print("Landmask directory: {}".format(landmask_dir))
        print("Landmask file: {}".format(landmask_file))

    landmask_path = os.path.join(landmask_dir, landmask_file)
    [mimage, mhead, mdescrip, miaopt] = sp2.loadsir(landmask_path)

    # # extract the header info
    # mhd = sp2.extractsirhead(mhead)
    # if verbose:
    #     sp2.printsirhead(mhd)

    # Okay start to convert data.  There are a few steps:
    # 1) change nodata values to NaNs
    # 2) change landmask 0's to NaNs
    # 3) multiply sig0 image by landmask
    # 4) convert NaNs to -9999

    # 1) change nodata values to NaNs
    image[image == nodata_value] = np.nan

    # 2) change landmask 0's to NaNs
    mimage = mimage.round()
    mimage[mimage == 0] = np.nan

    # 3) apply land mask
    masked_image = np.multiply(image, mimage)

    # 4) convert NaNs to -9999
    nodata_value = -9999.0
    masked_image[np.isnan(masked_image)] = nodata_value

    # finally, save data as a geotiff using GDAL library routines
    format = "GTiff"
    driver = gdal.GetDriverByName(format)
    src_filename = os.path.basename(infile)
    # src_dir = os.path.dirname(infile)
    dst_filename = os.path.splitext(os.path.basename(src_filename))[0] + ".tif"
    dst_dir = os.path.join(datadir, "geotiffs", region, str(year))
    dst_path = os.path.join(dst_dir, dst_filename)

    # create the destination directory if it doesn't already
    # exist
    Path(dst_dir).mkdir(parents=True, exist_ok=True)

    if verbose:
        print("Output file: {}".format(dst_path))

    nx = hd["nsx"]
    ny = hd["nsy"]
    dst_ds = driver.Create(dst_path, nx, ny, 1, gdal.GDT_Float32)
    dst_ds.GetRasterBand(1).WriteArray(masked_image)
    dst_ds.GetRasterBand(1).SetNoDataValue(nodata_value)

    dx = hd["grid_dx"]
    dy = hd["grid_dy"]
    xgrid_min = hd["xgrid_min"]
    ygrid_max = hd["ygrid_max"]
    geotransform = (xgrid_min, dx, 0.0, ygrid_max, 0.0, -1.0 * dy)

    csr = osr.SpatialReference()
    rc = csr.ImportFromProj4(proj4_string)

    if rc != 0:
        sys.stderr.write("Importing PROJ4 projection failed\n")
        sys.exit(1)
    else:
        projection = csr.ExportToWkt()
        if verbose:
            print("Output WKT Projection: {}".format(projection))

    dst_ds.SetGeoTransform(geotransform)
    dst_ds.SetProjection(projection)

    dst_ds = None
    sys.exit(0)


if __name__ == "__main__":

    main()
