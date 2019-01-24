#!/usr/bin/env python
"""
Created Apr 7, 2017
This script creates a png image of the contents of a SIR file using
default vmin/vmax values in file header

@author: DGL at BYU
"""

######
# Imports
######
from common import *
from loadsir import loadsir
import png as png  # PyPng module

######
# Export function definition
######
def save_it(fname, img_dir='.'):
    """ Helper function to export figures """
    if not os.access(img_dir, os.F_OK):
        os.mkdir(img_dir)
    plt.savefig('{0}/{1}.png'.format(img_dir,fname), dpi=150)
    #plt.savefig('{0}/{1}.pdf'.format(img_dir,fname))


######
# Plot the lat/lon box over the sir
######
def sir2png(sir_fname,png_fname,vmin=0,vmax=0):
    """ 
    sir2png(sir_fname,png_fname,vmin=vmin,vmax=vmax)

    sir_fname: input SIR file name
    png_fname: output png file name
    vmin,vmax: min,max values to display.  uses SIR head default if vmin=vmax

    """
    # Load in the SIR image
    sir = loadsir(sir_fname)
    # sir[0] is the image array; sir[1] is the header; sir[2]=iaopt; sir[3]=descript
    if vmin == vmax: # use default
        vmin=sir[1][49]
        vmax=sir[1][50]

    # create a .png image from the file data
    img=np.round((sir[0]-vmin)*256.0/(vmax-vmin))
    img[img>255]=255
    img[img<0]=0
    img=np.uint8(img)
    nsx=np.int(sir[1][0])
    nsy=np.int(sir[1][1])

    # write out as greyscale png
    f=open(png_fname,'wb')
    w=png.Writer(width=nsx,height=nsy,greyscale=True)
    w.write(f,img);f.close()
    print('png file %s created from %s' % (png_fname,sir_fname))

def main(argv): # pass in a filenaem as argument
    vmin=0
    vmax=0
    if len(argv)>1:
        sir_fname = argv[1]
        if len(argv)>2:
            png_fname = argv[2]
            if len(argv)>4:
                vmin=float(argv[3])
                vmax=float(argv[4])
        else:
            png_fname = sir_fname + '.png'
    else: # no argument?, use a test image
        sir_fname = 'greeni.sir'
        #sir_fname = 'queh-a-E2N04-10-10.sir'
        png_fname = sir_fname + '.png'

    sir2png(sir_fname,png_fname,vmin=vmin,vmax=vmax)

if __name__ == '__main__':
    main(sys.argv)
