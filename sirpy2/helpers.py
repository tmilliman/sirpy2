#!/usr/bin/env python
""" Small helper utility functions for sirpy"""


# Imports
import os
import sys
import numpy as np
from matplotlib import pyplot as plt
from printsirhead import printsirhead


# The helper functions
def clearline():
    """ Clears the line right of the cursor position """
    sys.stdout.write(chr(27) + "[2K")  # ANSI escape sequence


def printit(msg):
    """ Print a line to stdout, reset cursor to beginning of line """
    print(msg + "\r", end="")
    sys.stdout.flush()
    clearline()


def _save_it(fname, img_dir="."):
    """ Helper function to export figures """
    if not os.access(img_dir, os.F_OK):
        os.mkdir(img_dir)
    plt.savefig(
        "{}.png".format(os.path.join(img_dir, fname)), bbox_inches="tight", dpi=150
    )
