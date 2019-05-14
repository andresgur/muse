# -*- coding: utf-8 -*
# Script to clean images using a signal-to-noise ratio
# Author: Andres Gurpide Lasheras  andres.gurpide@gmail.com
# 20-04-2019
    #!/usr/bin/env python3

# imports
import argparse
import os
import sys
from astropy.io import fits
import numpy as np

# read arguments
ap = argparse.ArgumentParser(description='Spectrum fits to be loaded')
ap.add_argument("input_images", nargs='+', help="List of images fits files to be cleaned")
ap.add_argument("-o", "--outdir", nargs='?', help="Output dir", default='cleaned_images')
ap.add_argument("-t", "--threshold", nargs='?', help="Threshold signal to noise ratio to clean the input images", default=5, type=float)
ap.add_argument("-s", "--signaltonoiseratiomap", nargs=1, help="Signal to noise ratio map to clean the input images")
args = ap.parse_args()

images_toclean = args.input_images
sthreshold = args.threshold
outdir = args.outdir
signal_map = args.signaltonoiseratiomap[0]

if not os.path.isdir(outdir):
    os.mkdir(outdir)

if os.path.isfile(signal_map):
    signal_map = fits.open(signal_map)

    print('Loaded signal to noise ratio map successfully')
else:
    print("Signal to noise ratio map %s does not exist!" % signal_map)
    sys.exit()
for image in images_toclean:
    if os.path.isfile(image):
        hdul = fits.open(image)
        data = hdul[0].data
        hdul[0].data[np.where(signal_map[0].data < sthreshold)] = None
        hdul[0].header['COMMENT'] = "Used a signal to noise ratio of %.2f to filter this image" % sthreshold
        hdul.writeto(outdir + "/" + "clean" + image, overwrite=True)
        print("Cleaned and stored %s image" % image)
        hdul.close()
    else:
        print("Image %s not found." % image)
signal_map.close()
