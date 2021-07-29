# sirpy2 Repository

Python routines for interacting with SIR formatted image files
from the BYU Scatterometer Climate Record Pathfinder (SCP) project:

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

## System Requirements

The sir2geotiff utility uses the GDAL python binding.  To build this
python binding you will need to have the GDAL library installed and
the gdal library header files available.  This is usually done at a
system level and the procedure depends on the operating system used.

## Python Virtual Environments

All development and use of these scripts was done using virtual
environments on computers running a unix-like (linux or OSX).  To set
up a local environment for using the script you can set up a virtual
environment either using python native virtual environments (or
virtualenv) or the equivalent conda/anaconda virtual environments.

### Using Python's venv module

To setup and activate the virtual environment in the directory where
you have extracted the git repository:

    python3 -m venv .venv
    . .venv/bin/activate
    python -m pip install git+https://github.com/tmilliman/sirpy2.git 


### Using conda based environments

You can find documentation on using conda environments at
[https://docs.conda.io/en/latest/](https://docs.conda.io/en/latest/)
Again in the directory where you have extracted the git repository:

    conda env create --file environment.yml
