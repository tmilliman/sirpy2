#!/usr/bin/env python
""" Small helper utility functions for sirpy"""

######
# Imports
######
from __future__ import print_function, division
from common import *

##############################
# The helper functions
##############################
def clearline():
    """ Clears the line right of the cursor position """
    sys.stdout.write(chr(27) + '[2K') # ANSI escape sequence

def printit(msg):
    """ Print a line to stdout, reset cursor to beginning of line """
    print(msg+'\r', end='')
    sys.stdout.flush()
    clearline()

def _save_it(fname, img_dir='.'):
    """ Helper function to export figures """
    if not os.access(img_dir, os.F_OK):
        os.mkdir(img_dir)
    plt.savefig('{}.png'.format(os.path.join(img_dir,fname)), 
            bbox_inches='tight', dpi=150)
    #plt.savefig('{}.pdf'.format(img_dir,fname))

