# -*- coding: utf-8 -*
# Script to create image cbar images from fits image file
# Author: Andres Gurpide Lasheras  andres.gurpide@gmail.com
# 25-03-2019
    #!/usr/bin/env python3

# imports
import argparse
import os
import sys
sys.path.insert(0, '/home/agurpide/scripts/pythonscripts')
import plot_utils.plot_functions as pf
import muse_utils as muse_utils


# read arguments
ap = argparse.ArgumentParser(description='Spectrum fits to be loaded')
ap.add_argument("input_images", nargs='+', help="List of images to be saved")
ap.add_argument("-z", "--zlabel", nargs='?', help="Label for the z axis", default=None)
ap.add_argument("-o", "--outdir", nargs='?', help="Output dir", default='images')
args = ap.parse_args()

images = args.input_images
zlabel = args.zlabel
outdir = args.outdir

if not os.path.isdir(outdir):
    os.mkdir(outdir)

for image in images:
    if os.path.isfile(image):
        print("Storing image %s" % image)
        outputfilename = image.replace(".fits", "")
        muse_utils.plot_image(image, outdir + "/" + outputfilename, zlabel)
