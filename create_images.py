# !/usr/bin/env python3
# -*- coding: utf-8 -*
# Script to create image cbar images from fits image file
# Author: Andres Gurpide Lasheras  andres.gurpide@gmail.com
# 25-03-2019

# imports
import argparse
import os
import sys
import logging
import muse_utils as muse_utils

# read arguments
ap = argparse.ArgumentParser(description='Spectrum fits to be loaded')
ap.add_argument("input_images", nargs='+', help="List of images to be saved")
ap.add_argument("-o", "--outdir", nargs='?', help="Output dir", default='images')
ap.add_argument("-r", "--region", nargs='?', help="Region file to be overlaid on the image", default=None)
ap.add_argument("-c", "--cut", nargs='?', help="Region to cut the original image", default=None, type=str)
ap.add_argument("--vmin", nargs='?', help="Minimum value for the z scale of the image", default=None, type=float)
ap.add_argument("--vmax", nargs='?', help="Maximum value for the z scale of the image", default=None, type=float)
ap.add_argument("--add", help="Switch to plot images in same figure, default is False", action='store_true')

args = ap.parse_args()

images = args.input_images
outdir = args.outdir
region_file = args.region
spatial_cut = args.cut
vmin = args.vmin
vmax = args.vmax

scriptname = os.path.basename(__file__)
# logger

logger = logging.getLogger(scriptname)
logger.setLevel(logging.DEBUG)
out_format = logging.Formatter('%(name)s - %(levelname)s - %(message)s')

# handler
stream_handler = logging.StreamHandler()
stream_handler.setLevel(logging.DEBUG)
stream_handler.setFormatter(out_format)
logger.addHandler(stream_handler)

if not os.path.isdir(outdir):
    os.mkdir(outdir)

if args.add:

    if len(images) == 1:
        nrow = 1
        ncols = 1
    elif len(images) == 4:
        nrows = 2
        ncols = 2
    elif len(images) == 2:
        nrows = 1
        ncols = 2
    else:
        logger.error("Add mode only works for 1, 2 or 4 images")
        sys.exit()
    logger.info("Creating plot grid with %i rows and %i columns" % (nrows, ncols))
    images_string = [image.replace(".fits", "") for image in images]
    outname = "".join(images_string)
    muse_utils.create_grid(outdir + "/" + outname, *images, region=region_file, cutting_region=spatial_cut, vmin=vmin, nrows=nrows, ncols=ncols, colorbar=None)

else:
    for image in images:
        if os.path.isfile(image):
            logging.info("Storing image %s" % image)
            outputfilename = os.path.basename(image).replace(".fits", "")
            muse_utils.plot_image(outdir + "/" + outputfilename, image, region=region_file, cutting_region=spatial_cut, vmin=vmin, vmax=vmax)
        else:
            logging.warning("Image %s not found." % image)
