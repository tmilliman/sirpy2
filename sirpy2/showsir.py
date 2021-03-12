#!/usr/bin/env python
"""
Created on Sep 21, 2011
Revised on Apr 7, 2017 + include EASE2 support

Plots a SIR image from a SIR file using header-default scaling

@author: DGL
"""

# Imports
import os
import sys
import numpy as np
from matplotlib import pyplot as plt
from loadsir import loadsir
from printsirhead import printsirhead


# Export function definition
def save_it(fname, img_dir="."):
    """ Helper function to export figures """
    if not os.access(img_dir, os.F_OK):
        os.mkdir(img_dir)
    plt.savefig("{0}/{1}.png".format(img_dir, fname), dpi=150)


# Display the contents of a sir image file
def showsir(sir_fname):
    """ Plots a sir image

    INPUTS:
      sir_fname: the SIR file name.

    OUTPUTS (on display):
       The SIR image and header

    """

    # Load in the SIR image
    sir = loadsir(sir_fname)
    # sir[0] is the image array; sir[1] is the header; sir[2]=iaopt; sir[3]=descript

    print("File name: ", sir_fname)
    printsirhead(sir[1])

    # Draw SIR image
    print(" ")
    print("Display SIR image")
    plt.figure()
    plt.imshow(sir[0], cmap=plt.cm.gray, interpolation="nearest", aspect="equal")

    # add axes and labels
    plt.axis("image")
    plt.title(sir_fname)
    plt.xlabel("x")
    plt.ylabel("y")


def main(sir_fname):
    sir = showsir(sir_fname)
    plt.show()
    return sir


if __name__ == "__main__":
    main(sys.argv[1])
