sirpy2 Repository
=================

Python routines for interacting with SIR formatted image files
from the BYU Scatterometer Climate Record Pathfinder project:

http://www.scp.byu.edu/

Many of the routines here are based on or directly copied from the `sirpy`
routines written by David G. Long and provided
[here](http://www.scp.byu.edu/software/sirpy/).

The motivation for starting this package was to be able to generate
GeoTIFF files with a float data type.  The `sir2geotiff.py` script
provided by the `sirpy` scripts above which (at the time) scaled the
data to byte values and generated a grayscale image with a 0-255 range
of values.  I assume this was done because many image display
utilities and software libraries could not handle the float data type.
To work with the images using utilities like GDAL or with GIS programs
like QGIS which do handle other data types I wanted to have the full
range data values in the GeoTIFF files.

---------------
