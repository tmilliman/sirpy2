# sirpy2 Repository

Python routines for interacting with SIR formatted image files
from the BYU Scatterometer Climate Record Pathfinder (SCP) project:

http://www.scp.byu.edu/

Many of the routines here are based on or directly copied from the `sirpy`
routines written by David G. Long and provided
[here](http://www.scp.byu.edu/software/sirpy/).

The motivation for starting this package was to be able to generate
GeoTIFF files with a `float32` data type.  The `sir2geotiff.py` script
provided by the `sirpy` scripts above which (at the time) scaled the
data to byte values and generated a grayscale image with a 0-255 range
of values.  I assume this was done because many image display
utilities and software libraries could not handle the float data type.
To work with the images using utilities like GDAL or with GIS programs
like QGIS which do handle other data types I wanted to have the full
range data values in the GeoTIFF files.

# Installation

## Python Virtual Environments

All development and use of these scripts was done using virtual
environments on computers running a unix-like (linux or OSX).  To set
up a local environment for using the script you can set up a virtual
environment either using python native virtual environments (or
virtualenv) or the equivalent conda/anaconda virtual environments.

### Using Python3's venv module and pip

The sir2geotiff utility uses the GDAL python binding.  To build this
python binding in a virtual environment you will need to have the GDAL
library installed and the gdal library header files available.  This
is usually done at a system level and the procedure depends on the
operating system used.  Once the GDAL system libraries are installed
you can create a virtualenv, activate and install the sirpy2 package
by doing the following:

    python3 -m venv .venv
    . .venv/bin/activate
    pip install -U pip
    python -m pip install git+https://github.com/tmilliman/sirpy2.git 

To test the installation you can issue the following command:

   sir2geotiff --help
   
which will display the script's help message if the installation was
successful.

### Using conda/miniconda based environments

You can find documentation on using conda environments at
[https://docs.conda.io/en/latest/](https://docs.conda.io/en/latest/)
The `conda` package manager can install GDAL into the virtual
environment locally so there is no need to install it system wide.  In
the directory where you have extracted the sirpy2 git repository you
can create a conda environment, activate it and install the sirpy2
package by doing the following:

    conda env create --file environment.yml
    conda activate sirpy
    python -m pip install git+https://github.com/tmilliman/sirpy2.git

To test the installation you can issue the following command:

   sir2geotiff --help

which will display the script's help message if the installation was
successful.

## Usage Notes

After the `sirpy2` package is installed you will have 4 command
line scripts installed available from the command line:

     printsirhead
     showsir
     sir2png
     sir2geotiff

These scripts can be used to examine the `.sir` image data files from
the BYU SCP sites.  There is a sample image in the `tests` sub-directory
of the repository.

### `printsirhead`

      usage: printsirhead [-h] filename

      Print summary of SIR header values

      positional arguments:
        filename    SIR filename

      optional arguments:
        -h, --help  show this help message and exit

#### Example

     printsirhead msfa-a-NAm07-181-185.sir
     File name: msfa-a-NAm07-181-185.sir
     Lambert form (local radius) (2)
       x,y scale: (km/pix) 4.450000 4.450000  11632
       x,y corner: (km)    -4200.000000 -2300.000000  0 0 1
       x,y span: (deg)     -92.500000 45.000000  0 0 100
       PROJ4 string: +proj=laea +lat_0=45.0 +lon_0=-92.5 +k=1 +x_0=0 +y_0=0 +R=6367415.828 +units=m +no_defs
       GeoTransform: -a_ullr -4200000.000 2817500.000 4210500.000 -2300000.000
       Pixels:  1890 by 1150
       Year: 2007  Start Day: 181   Min:  0  End Day: 185  Min: 2
       Region: 205  Type: 1  Htype: 31  Nhead: 1
       Nodata value: -33.000000  Vmin: -32.000000  Vmax: 0.000000
       Datatype: short  offset: -33  scale: 1000
       Frequency: 5.30   Polarization: 2
       Type: A image  (msfa-a-NAm07-181-185-s102-02.sir)
       Sensor: ASCAT-A (ASCAT on MetOp-A)
       Title: SIR A image of north-ame
       Tag: (c) 2014 BYU MERS Laboratory
       Creator: BYU MERS:ascat_meta_sir3 version 0.12 Priors=0
       Created: 2015-08-05 07:10:42-0600


### `showsir`

    usage: showsir [-h] filename

    Display SIR images file.

    positional arguments:
      filename    SIR filename

    optional arguments:
      -h, --help  show this help message and exit

#### Example

     showsir msfa-a-NAm07-181-185.sir

     results in a pop window as shown below

     [showsir preview image](https://github.com/tmilliman/sirpy2/docs/showsir_figure.png)