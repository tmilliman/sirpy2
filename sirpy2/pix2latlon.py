"""
Created on Sep 21, 2011
Revised on Apr 7, 2017 + include EASE2 support

@author: Bradley, DGL
"""
from numpy import (
    cos,
    sin,
    tan,
    mod,
    sqrt,
    abs,
    arcsin,
    arccos,
    arctan,
    arctan2,
    array,
    all,
    sign,
    zeros,
    nonzero,
    log,
    size,
)
from latlon2pix import polster
from ease2helper import ease2_map_info, easeconv_normalize_degrees


def pix2latlon(x, y1, head):
    """Pixels to latitude/longitude
    function (lon, lat)=pix2latlon(x,y,head)

    Given an image pixel location (x,y) (1..nsx,1..nsy)
    computes the lat,lon coordinates (lon,lat).   The lat,lon returned
    corresponds to the lower-left corner of the pixel.  If lat,lon
    of pixel center is desired use (x+0.5,y+0.5) where x,y are integer
    valued pixels

    Note:  while routine will attempt to convert any (x,y)
    values, only (x,y) values with 1 <= x <= nsx+1 and 1 <= y <= nsy+1
    are contained within image.

    INPUTS:
        x,y - input pixel location (matlab coordinates y_matlab=nxy-y_sir+1)
        head - header array from loadsir

    OUTPUTS:
        lon,lat - longitude, latitude
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

    y = nsy - y1 - 1.0  # convert from MATLAB coord to SIR coordinates

    if iopt == -1:  # image only (can't transform!)
        thelon = x / ascale + a0
        thelat = y / bscale + b0
        alon = thelon
        alat = thelat
    elif iopt == 0:  # rectalinear lat/long
        thelon = x / ascale + a0
        thelat = y / bscale + b0
        alon = thelon
        alat = thelat
    elif (iopt == 1) or (iopt == 2):  # lambert
        thelon = x / ascale + a0
        thelat = y / bscale + b0
        alon, alat = ilambert1(thelon, thelat, ydeg, xdeg, iopt)
    elif iopt == 5:  # polar stereographic
        thelon = x * ascale + a0
        thelat = y * bscale + b0
        alon, alat = ipolster(thelon, thelat, xdeg, ydeg)
    elif (iopt == 8) or (iopt == 9) or (iopt == 10):  # EASE2
        thelon = x - 1.0 + a0
        thelat = y - 1.0 + b0
        alon, alat = iease2grid(iopt, thelon, thelat, ascale, bscale)
    elif (iopt == 11) or (iopt == 12) or (iopt == 13):  # EASE
        thelon = x - xdeg + (xdeg + a0)
        thelat = y - ydeg + (ydeg + b0)
        alon, alat = ieasegrid(iopt, thelon, thelat, ascale)
    else:
        print("*** unknown SIR transformation ***")

    return alon, alat


def ilambert1(x, y, orglat, orglon, iopt):
    """Inverse Lambert azimuthal equal-area projection
    (lon, lat)=ilambert1(x,y,orglat,orglon,iopt)

    computes the inverse transformation from lat/lon to x/y for the
    lambert azimuthal equal-area projection

    inputs:
        lat    (r): latitude +90 to -90 deg with north positive
        lon    (r): longitude 0 to +360 deg with east positive
            or -180 to +180 with east more positive
        orglat (r): origin parallel +90 to -90 deg with north positive
        orglon (r): central meridian (longitude) 0 to +360 deg
            or -180 to +180 with east more positive
        iopt   (i): earth radius option
            for iopt=1 a fixed, nominal earth radius is used.
            for iopt=2 the local radius of the earth is used
    outputs:
        x,y    (r): rectangular coordinates in km


    see "map projections used by the u.s. geological survey"
    geological survey bulletin 1532, pgs 157-173

    for this routine, a spherical earth is assumed for the projection.
    the error will be small for small-scale maps.
    for local radius the 1972 wgs ellipsoid model (bulletin pg 15).
    """
    dtr = 3.141592654 / 180.0
    radearth = 6378.135  # equitorial earth radius
    f = 298.26  # 1/f
    orglon1 = mod(orglon + 720.0, 360.0)
    #
    # compute local radius of the earth at center of image
    #
    eradearth = 6378.0  # use fixed nominal value
    if iopt == 2:  # use local radius
        era = 1.0 - 1.0 / f
        eradearth = (
            radearth
            * era
            / sqrt(era * era * cos(orglat * dtr) ** 2 + sin(orglat * dtr) ** 2)
        )

    x1 = x / eradearth
    y1 = y / eradearth
    rho = x1 * x1 + y1 * y1
    if all(rho > 0):
        rho = sqrt(rho)
        c = 2 * arcsin(rho * 0.5)
        lat = (
            arcsin(cos(c) * sin(orglat * dtr) + y1 * sin(c) * cos(orglat * dtr) / rho)
            / dtr
        )
    else:
        lat = orglat

    lon = 0
    if abs(orglat) != 90.0:
        if all(rho == 0):
            lon = orglon1
        else:
            t1 = x1 * sin(c)
            t2 = rho * cos(orglat * dtr) * cos(c) - y1 * sin(orglat * dtr) * sin(c)
            lon = orglon1 + arctan2(t1, t2) / dtr
    elif orglat == 90.0:
        lon = orglon1 + arctan2(x1, -y1) / dtr
    else:
        lon = orglon1 + arctan2(x1, y1) / dtr

    lon = mod(lon + 720.0, 360.0)
    if all(lon > 180):
        lon = lon - 360

    alon = lon
    return alon, lat


def ipolster(x, y, xlam, slat):
    """Inverse polar stereographic transformation
    (lon,lat)=ipolster(x,y,xlam,slat)

    computes the inverse polar sterographic transformation for (x,y)
    given in km with references lon,lat=(xlam,slat).
    output lon,lat=alon,alat

    algorithm is the same as used for processing ers-1 sar images
    as received from m. drinkwater (1994).  updated by d. long to
    improve accuracy using iteration with forward transform.
    """
    e2 = 0.006693883
    re = 6378.273
    pi2 = 1.570796327
    dtr = pi2 / 90.0
    #
    # first use approximate inverse calculation
    #
    e = sqrt(e2)
    e22 = e2 * e2
    e23 = e2 * e2 * e2
    x1 = array([x])
    y1 = array([y])
    rho = x1 * x1 + y1 * y1
    rho[rho > 0] = sqrt(rho[rho > 0])
    if all(rho < 0.05):
        alon = xlam
        alat = sign(90.0, slat)
    else:
        sn = 1.0
        slat1 = slat
        if slat < 0:
            sn = -1.0
            slat1 = -slat

        cm = cos(slat1 * dtr) / sqrt(1.0 - e2 * sin(slat1 * dtr) ** 2)
        t = tan(dtr * (45.0 - 0.5 * slat1)) / (
            (1.0 - e * sin(slat1 * dtr)) / (1.0 + e * sin(slat1 * dtr))
        ) ** (e * 0.5)
        t = rho * t / (re * cm)
        chi = pi2 - 2.0 * arctan(t)
        t = (
            chi
            + (0.5 * e2 + 5.0 * e22 / 24.0 + e23 / 12.0) * sin(2.0 * chi)
            + (7.0 * e22 / 48.0 + 29.0 * e23 / 240.0) * sin(4.0 * chi)
            + (7.0 * e23 / 120.0) * sin(6.0 * chi)
        )
        alat = sn * (t * 90.0 / pi2)
        alon = sn * arctan2(sn * x1, -sn * y1) / dtr + xlam
        alon[alon < -180] = alon[alon < -180] + 360.0
        alon[alon > 180] = alon[alon > 180] - 360.0
    #
    # using the approximate result as a starting point, iterate to
    # improve the accuracy of the inverse solution
    #
    sn1 = 1.0
    if slat < 0:
        sn1 = -1.0

    a = arctan2(y, x) / dtr
    r = sqrt(x * x + y * y)

    for icnt in range(1, 21):
        xx, yy = polster(alon, alat, xlam, slat)
        rr = sqrt(xx * xx + yy * yy)
        rerr = sn1 * (rr - r) / 180.0
        aa = arctan2(yy, xx) / dtr
        aerr = aa - a
        aerr[abs(aerr) > 180] = 360.0 - aerr[abs(aerr) > 180]
        #
        # check for convergence
        #
        if (
            (max(abs(rerr.flatten(1))) < 0.001) and (max(abs(aerr.flatten(1))) < 0.001)
        ) or (icnt > 9):
            if max(abs(alon.flatten(1))) > 360.0:
                alon = mod(alon, 360)
            return alon, alat
        #
        # constrain updates
        #
        alon = alon + aerr
        if max(abs(alon.flatten(1))) > 360:
            alon = mod(alon, 360)
        if alat * slat < 0:
            rerr = rerr * (1.0 - sin(dtr * abs(alat)))
            rerr[abs(rerr) > 2] = sign(2.0, rerr[abs(rerr) > 2]) / icnt
            alat = alat + rerr
            alat[abs(alat) > 90] = sign(90, alat[abs(alat) > 90])

    return alon, alat


def ieasegrid(iopt, thelon, thelat, ascale):
    """Inverse EASE grid transformation
    (lon, lat)=ieasegrid(iopt,thelon,thelat,ascale)

    computes the inverse "ease" grid transform

    given the image transformation coordinates (thelon,thelat) and
    the scale (ascale) the corresponding lon,lat (lon,lat) is computed
    using the "ease grid" (version 1.0) transformation given in fortran
    source code supplied by NSIDC
    iopt is ease type: iopt=11=north, iopt=12=south, iopt=13=cylindrical

    the radius of the earth used in this projection is imbedded into
    ascale while the pixel dimension in km is imbedded in bscale
    the base values are: radius earth= 6371.228 km
                 pixel dimen =25.067525 km
    then, bscale = base_pixel_dimen
          ascale = radius_earth/base_pixel_dimen
    """
    pi2 = 1.57079633  # pi/2 at standard precision
    dtr = pi2 / 90.0
    x1 = thelon
    y1 = thelat
    temp = zeros(thelon.ndim)
    if iopt == 11:  # ease grid north
        alon = arctan2(x1, -y1) / dtr
        # vectorized code
        alat = zeros(alon.shape)
        ind = nonzero(abs(sin(dtr * alon)) > abs(cos(alon * dtr)))
        temp[ind] = (x1[ind] / sin(alon[ind] * dtr)) / ascale
        ind = nonzero(abs(sin(dtr * alon)) <= abs(cos(alon * dtr)))
        temp[ind] = (-y1[ind] / cos(alon[ind] * dtr)) / ascale
        alat[abs(temp) <= 1] = 90.0 - 2.0 * arcsin(temp[abs(temp) <= 1]) / dtr
        alat[abs(temp) > 1] = 90.0 * sign(temp[abs(temp) > 1])
    elif iopt == 12:  # ease grid south
        alon = arctan2(x1, y1) / dtr
        # vectorized code
        alat = zeros(alon.shape)
        ind = abs(cos(alon * dtr)) > abs(sin(alon * dtr))
        temp[ind] = (y1[ind] / cos(alon[ind] * dtr)) / ascale
        ind = abs(cos(alon * dtr)) <= abs(sin(alon * dtr))
        temp[ind] = (x1[ind] / sin(alon[ind] * dtr)) / ascale
        alat[abs(temp) <= 1] = 90.0 - 2.0 * arccos(temp[abs(ind) <= 1]) / dtr
        alat[abs(temp) > 1] = 90.0 * sign(temp[abs(temp) > 1])
    elif iopt == 13:  # ease cylindrical
        alon = ((x1 / ascale) / cos(30.0 * dtr)) * 90.0 / pi2
        temp = (y1 * cos(30.0 * dtr)) / ascale
        # vectorized code
        alat = zeros(alon.shape)
        alat[abs(temp) <= 1] = arcsin(temp[abs(temp) <= 1]) / dtr
        alat[abs(temp) > 1] = 90.0 * sign(temp[abs(temp) > 1])

    return alon, alat


def iease2grid(iopt, thelon, thelat, ascale, bscale):
    """Inverse EASE2 grid transformation
    (lon, lat)=iease2grid(iopt,thelon,thelat,ascale,bscale)

    computes the inverse "ease2" grid transform

        RADIUS EARTH=6378.137 KM (WGS 84)
        MAP ECCENTRICITY=0.081819190843 (WGS84)

        inputs:
          iopt: projection type 8=EASE2 N, 9-EASE2 S, 10=EASE2 T/M
          thelon: X coordinate in pixels (can be outside of image)
          thelat: Y coordinate in pixels (can be outside of image)
          alon, alat: lon, lat (deg) to convert (can be outside of image)
          ascale and bscale should be integer valued)
          ascale: grid scale factor (0..5)  pixel size is (bscale/2^ascale)
          bscale: base grid scale index (ind=int(bscale))

          NSIDC .grd file for isc=0
           project type    ind=0     ind=1         ind=3
               N         EASE2_N25km EASE2_N30km EASE2_N36km
               S         EASE2_S25km EASE2_S30km EASE2_S36km
              T/M        EASE2_T25km EASE2_M25km EASE2_M36km

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
          alon, alat: lon, lat location in deg  (can be outside of image)

    """
    RTD = 57.29577951308232
    ind = round(bscale)
    isc = round(ascale)

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

    e4 = e2 * e2
    e6 = e4 * e2

    # qp is the function q evaluated at phi = 90.0 deg
    qp = (1.0 - e2) * (
        (1.0 / (1.0 - e2))
        - (1.0 / (2.0 * map_eccentricity))
        * log((1.0 - map_eccentricity) / (1.0 + map_eccentricity))
    )

    x = (thelon - r0 - 0.5) * map_scale
    y = (thelat - 0.5 - s0) * map_scale

    if iopt == 8:  # EASE2 grid north
        rho2 = (x * x) + (y * y)
        arg = 1.0 - (rho2 / (map_equatorial_radius_m * map_equatorial_radius_m * qp))
        if size(arg) > 1:
            arg[arg > 1.0] = 1.0
            arg[arg < -1.0] = -1.0
        else:
            if arg > 1.0:
                arg = 1.0
            if arg < -1.0:
                arg = -1.0

        beta = arcsin(arg)
        lam = arctan2(x, -y)

    elif iopt == 9:  # EASE2 grid south
        rho2 = (x * x) + (y * y)
        arg = 1.0 - (rho2 / (map_equatorial_radius_m * map_equatorial_radius_m * qp))
        if size(arg) > 1:
            arg[arg > 1.0] = 1.0
            arg[arg < -1.0] = -1.0
        else:
            if arg > 1.0:
                arg = 1.0
            if arg < -1.0:
                arg = -1.0
        beta = -arcsin(arg)
        lam = arctan2(x, y)

    elif iopt == 10:  # EASE2 cylindrical
        arg = 2.0 * y * kz / (map_equatorial_radius_m * qp)
        if size(arg) > 1:
            arg[arg > 1.0] = 1.0
            arg[arg < -1.0] = -1.0
        else:
            if arg > 1.0:
                arg = 1.0
            if arg < -1.0:
                arg = -1.0
        beta = arcsin(arg)
        lam = x / (map_equatorial_radius_m * kz)
    else:
        print("*** invalid EASE2 projection specificaion in iease2grid")
    phi = (
        beta
        + (
            ((e2 / 3.0) + ((31.0 / 180.0) * e4) + ((517.0 / 5040.0) * e6))
            * sin(2.0 * beta)
        )
        + ((((23.0 / 360.0) * e4) + ((251.0 / 3780.0) * e6)) * sin(4.0 * beta))
        + (((761.0 / 45360.0) * e6) * sin(6.0 * beta))
    )

    alat = RTD * phi
    alon = map_reference_longitude + (RTD * lam)
    alon = easeconv_normalize_degrees(alon)

    return (alon, alat)
