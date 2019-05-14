# -*- coding: utf-8 -*
# Script to create image cbar images from fits image file
# Author: Andres Gurpide Lasheras  andres.gurpide@gmail.com
# 25-03-2019
    #!/usr/bin/env python3

# imports
import argparse
import os
import sys
import logging
sys.path.insert(0, '/home/agurpide/scripts/pythonscripts')
import muse_utils as muse_utils

# read arguments
ap = argparse.ArgumentParser(description='Spectrum fits to be loaded')
ap.add_argument("input_images", nargs='+', help="List of images to be saved")
ap.add_argument("-o", "--outdir", nargs='?', help="Output dir", default='images')
ap.add_argument("-r", "--region", nargs='?', help="Region file to be overlaid on the image", default=None)
ap.add_argument("-c", "--cut", nargs='?', help="Region to cut the original image", default=None, type=str)
ap.add_argument("--vmin", nargs='?', help="Minimum value for the z scale of the image", default=None, type=float)
args = ap.parse_args()

images = args.input_images
outdir = args.outdir
region_file = args.region
spatial_cut = args.cut
vmin = args.vmin

scriptname = os.path.basename(__file__)
# logger

logging.basicConfig(level=logging.DEBUG, format='%(message)s')
logger = logging.getLogger(scriptname)

if not os.path.isdir(outdir):
    os.mkdir(outdir)


for image in images:
    if os.path.isfile(image):
        logging.info("Storing image %s" % image)
        outputfilename = image.replace(".fits", "")
        muse_utils.plot_image(image, outdir + "/" + outputfilename, region=region_file, cutting_region=spatial_cut, vmin=vmin)
    else:
        logging.warning("Image %s not found." % image)
