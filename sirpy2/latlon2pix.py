"""
Created on Sep 21, 2011
Revised on Apr 7, 2017 + include EASE2 support

@author: Bradley, DGL
"""
import numpy as np
from numpy import cos, sin, tan, mod, sqrt

from .ease2helper import ease2_map_info, easeconv_normalize_degrees


def latlon2pix(alon, alat, head):
    """Latitude/longitude to pixels
    (x, y) = latlon2pix(lon,lat,head)

    Convert a lat,lon coordinate (lon,lat) to an image pixel location
    (x,y) (in floating point, matlab convention).
    To compute integer pixel indices (ix,iy): check to insure
    1 <= x < nsx+1 and 1 <= y < nsy+1 then ix=floor(x) iy=floor(y)

    INPUTS:
        lon,lat - longitude, latitude
        head - header array from load sir

    OUTPUTS:
        x,y - pixel location (matlab coordinates y_matlab=nxy-y_sir+1)
    """
    # nsx = head[0]
    nsy = head[1]
    iopt = head[16]
    xdeg = head[2]
    ydeg = head[3]
    ascale = head[5]
    bscale = head[6]
    a0 = head[7]
    b0 = head[8]

    if iopt == -1:  # image only (can't transform!)
        x = ascale * (alon - a0)
        y = bscale * (alat - b0)
    elif iopt == 0:  # rectalinear lat/lon
        thelon = alon
        thelat = alat
        x = ascale * (thelon - a0)
        y = bscale * (thelat - b0)
    elif (iopt == 1) or (iopt == 2):  # lambert
        thelon, thelat = lambert1(alat, alon, ydeg, xdeg, iopt)
        x = ascale * (thelon - a0)
        y = bscale * (thelat - b0)
    elif iopt == 5:  # polar stereographic
        thelon, thelat = polster(alon, alat, xdeg, ydeg)
        x = (thelon - a0) / ascale
        y = (thelat - b0) / bscale
    elif (iopt == 8) or (iopt == 9) or (iopt == 10):  # EASE2
        thelon, thelat = ease2grid(iopt, alat, alon, ascale, bscale)
        x = thelon + 1.0 - a0
        y = thelat + 1.0 + b0
    elif (iopt == 11) or (iopt == 12) or (iopt == 13):  # EASE
        thelon, thelat = easegrid(iopt, alat, alon, ascale)
        thelon = thelon + xdeg
        thelat = thelat + ydeg
        x = thelon - (xdeg + a0)
        y = thelat - (ydeg + b0)
    else:
        print("*** Unknown SIR transformation: %d" % iopt)

    y = nsy - y - 1.0  # convert from matlab coordinates to SIR coordinates
    return x, y


def lambert1(lat, lon, orglat, orglon, iopt):
    """Lambert azimuthal equal-area projection
    (x,y)=lambert1(lat,lon,orglat,orglon,iopt)

    Computes the transformation from lat/lon to x/y for the
    lambert azimuthal equal-area projection

    inputs:
        lat    (r): latitude +90 to -90 deg with north positive
        lon    (r): longitude 0 to +360 deg with east positive
            or -180 to +180 with east more positive
        orglat    (r): origin parallel +90 to -90 deg with north positive
        orglon    (r): central meridian (longitude) 0 to +360 deg
            or -180 to +180 with east more positive
        iopt    (i): earth radius option
            for iopt=1 a fixed, nominal earth radius is used.
            for iopt=2 the local radius of the earth is used.
    outputs:
        x,y    (r): rectangular coordinates in km

    see "map projections used by the u.s. geological survey"
    geological survey bulletin 1532, pgs 157-173
    for this routine, a spherical earth is assumed for the projection
    the 1972 wgs ellipsoid model (bulletin pg 15).
    the error will be small for small-scale maps.
    """

    radearth = 6378.135  # equitorial earth radius
    f = 298.26  # 1/f wgs 72 model values
    dtr = 3.141592654 / 180.0

    lon1 = mod(lon + 720.0, 360.0)
    orglon1 = mod(orglon + 720.0, 360.0)
    #
    # compute local radius of the earth at center of image
    #
    # eradearth = 6378.0         # use fixed nominal value
    eradearth = radearth
    if iopt == 2:  # local radius
        era = 1.0 - 1.0 / f
        eradearth = (
            radearth
            * era
            / sqrt(era * era * cos(orglat * dtr) ** 2 + sin(orglat * dtr) ** 2)
        )

    denom = (
        1.0
        + sin(orglat * dtr) * sin(lat * dtr)
        + cos(orglat * dtr) * cos(lat * dtr) * cos(dtr * (lon1 - orglon1))
    )
    if np.all(denom > 0.0):
        ak = sqrt(2.0 / denom)
    else:
        print("*** division error in lambert1 routine ***")
        ak = 1.0

    x = ak * cos(lat * dtr) * sin(dtr * (lon1 - orglon1))
    y = ak * (
        cos(dtr * orglat) * sin(dtr * lat)
        - sin(dtr * orglat) * cos(dtr * lat) * cos(dtr * (lon1 - orglon1))
    )
    x = x * eradearth
    y = y * eradearth
    return x, y


def polster(alon, alat, xlam, slat):
    """Polar stereographic trasnformation
    (x,y)=polster(lon,lat,xlam,slat)

    computes the polar sterographic transformation for a lon,lat
    input of (alon,alat) with reference origin  lon,lat=(xlam,slat).
    output is (x,y) in km


    algorithm is the same as used for processing ers-1 sar images
    as received from m. drinkwater (1994)
    """
    # ported from polster.m by JPB 21 Sept 2011
    e2 = 0.006693883
    re = 6378.273
    dtr = 3.141592654 / 180.0
    e = sqrt(e2)
    if slat < 0:
        sn = -1.0
        rlat = -alat
    else:
        sn = 1.0
        rlat = alat

    t = ((1.0 - e * sin(rlat * dtr)) / (1.0 + e * sin(rlat * dtr))) ** (e * 0.5)
    ty = tan(dtr * (45.0 - 0.5 * rlat)) / t
    if slat < 0:
        rlat = -slat
    else:
        rlat = slat

    t = ((1.0 - e * sin(dtr * rlat)) / (1.0 + e * sin(dtr * rlat))) ** (e * 0.5)
    tx = tan(dtr * (45.0 - 0.5 * rlat)) / t
    cm = cos(dtr * rlat) / sqrt(1.0 - e2 * sin(dtr * rlat) ** 2)
    rho = re * cm * ty / tx
    x = (sn * sin(dtr * (sn * alon - xlam))) * rho
    y = -(sn * cos(dtr * (sn * alon - xlam))) * rho
    return x, y


def easegrid(iopt, alat, alon, ascale):
    """EASE grid transformation
    (thelon thelat)=easegrid(iopt,lat,lon,ascale)

    computes the forward "ease" grid transform

    given a lat,lon (alat,alon) and the scale (ascale) the image
    transformation coordinates (thelon,thelat) are comuted
    using the "ease grid" (version 1.0) transformation given in fortran
    source code supplied by nsidc.

    the radius of the earth used in this projection is imbedded into
    ascale while the pixel dimension in km is imbedded in bscale
    the base values are: radius earth= 6371.228 km
                 pixel dimen =25.067525 km
    then, bscale = base_pixel_dimen
          ascale = radius_earth/base_pixel_dimen

    iopt is ease type: iopt=11=north, iopt=12=south, iopt=13=cylindrical
    """
    # ported from easegrid.m by JPB 21 Sept 2011
    pi2 = np.pi / 2.0
    dtr = pi2 / 90.0

    if iopt == 11:  # ease grid north
        thelon = ascale * sin(alon * dtr) * sin(dtr * (45.0 - 0.5 * alat))
        thelat = ascale * cos(alon * dtr) * sin(dtr * (45.0 - 0.5 * alat))
    elif iopt == 12:  # ease grid south
        thelon = ascale * sin(alon * dtr) * cos(dtr * (45.0 - 0.5 * alat))
        thelat = ascale * cos(alon * dtr) * cos(dtr * (45.0 - 0.5 * alat))
    elif iopt == 13:  # ease cylindrical
        thelon = ascale * pi2 * alon * cos(30.0 * dtr) / 90.0
        thelat = ascale * sin(alat * dtr) / cos(30.0 * dtr)

    return thelon, thelat


def ease2grid(iopt, alat, alon, ascale, bscale):
    """EASE2 grid transformation

    (thelon thelat)=ease2grid(iopt,lat,lon,ascale,bscale)

        given a lat,lon (alat,alon) and the scale (ascale) the image
        transformation coordinates (thelon,thelat) are comuted
        using the "ease2 grid" (version 2.0) transformation given in IDL
        source code supplied by MJ Brodzik

        RADIUS EARTH=6378.137 KM (WGS 84)
        MAP ECCENTRICITY=0.081819190843 (WGS84)

        inputs:
          iopt: projection type 8=EASE2 N, 9-EASE2 S, 10=EASE2 T/M
          alon, alat: lon, lat (deg) to convert (can be outside of image)
          ascale and bscale should be integer valued)
          ascale: grid scale factor (0..5)  pixel size is (bscale/2^ascale)
          bscale: base grid scale index (ind=int(bscale))

          NSIDC .grd file for isc=0
           project type    ind=0     ind=1         ind=3
               N         EASE2_N25km EASE2_N30km EASE2_N36km
              S         EASE2_S25km EASE2_S30km EASE2_S36km
              T/M       EASE2_T25km EASE2_M25km EASE2_M36km

          cell size (m) for isc=0 (scale is reduced by 2^isc)
           project type    ind=0     ind=1            ind=3
               N          25000.0     30000.0         36000.0
              S          25000.0     30000.0         36000.0
              T/M       T25025.26   M25025.2600081  M36032.220840584

          for a given base cell size isc is related to NSIDC .grd file names
             isc        N .grd name   S .grd name   T .grd name
              0       EASE2_N25km     EASE2_S25km     EASE2_T25km
              1       EASE2_N12.5km   EASE2_S12.5km   EASE2_T12.5km
              2       EASE2_N6.25km   EASE2_S6.25km   EASE2_T6.25km
              3       EASE2_N3.125km  EASE2_S3.125km  EASE2_T3.125km
              4       EASE2_N1.5625km EASE2_S1.5625km EASE2_T1.5625km

        outputs:
          thelon: X coordinate in pixels (can be outside of image)
          thelat: Y coordinate in pixels (can be outside of image)
    """

    DTR = 0.01745329241994
    ind = round(bscale)
    isc = round(ascale)
    dlon = alon
    phi = DTR * alat
    lam = dlon

    # get base EASE2 map projection parameters
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
    ) = ease2_map_info(iopt, isc, ind)

    dlon = dlon - map_reference_longitude
    dlon = easeconv_normalize_degrees(dlon)
    lam = DTR * dlon

    sin_phi = np.sin(phi)

    q = (1.0 - e2) * (
        (sin_phi / (1.0 - e2 * sin_phi * sin_phi))
        - (1.0 / (2.0 * map_eccentricity))
        * np.log(
            (1.0 - map_eccentricity * sin_phi) / (1.0 + map_eccentricity * sin_phi)
        )
    )

    if iopt == 8:  # EASE2 grid north
        qp = 1.0 - (
            (1.0 - e2)
            / (2.0 * map_eccentricity)
            * np.log((1.0 - map_eccentricity) / (1.0 + map_eccentricity))
        )
        rho = map_equatorial_radius_m * np.sqrt(qp - q)
        if np.size(rho) > 1:
            rho[np.abs(qp - q) < epsilon] = 0.0
        else:
            if np.abs(qp - q) < epsilon:
                rho = 0

        x = rho * np.sin(lam)
        y = -rho * np.cos(lam)

    elif iopt == 9:  # EASE2 grid south
        qp = 1.0 - (
            (1.0 - e2)
            / (2.0 * map_eccentricity)
            * np.log((1.0 - map_eccentricity) / (1.0 + map_eccentricity))
        )
        rho = map_equatorial_radius_m * np.sqrt(qp + q)
        if np.size(rho) > 1:
            rho[np.abs(qp - q) < epsilon] = 0.0
        else:
            if np.abs(qp - q) < epsilon:
                rho = 0

        x = rho * np.sin(lam)
        y = rho * np.cos(lam)

    elif iopt == 10:  # EASE2 cylindrical
        x = map_equatorial_radius_m * kz * lam
        y = (map_equatorial_radius_m * q) / (2.0 * kz)

    else:
        print("*** invalid EASE2 projection specificaion in ease2grid")

    thelon = r0 + (x / map_scale) + 0.5
    thelat = s0 + (y / map_scale) + 0.5

    return (thelon, thelat)
