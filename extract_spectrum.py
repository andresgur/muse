
# -*- coding: utf-8 -*-
# Script to extract a spectrum from a certain region (circle) around a center of the cube. \
# The final spectrum is the mean spectrum of the region
# Author: Andres Gurpide Lasheras  andres.gurpide@gmail.com
# 25-03-2019
# !/usr/bin/env python3
# imports
from mpdaf.obj import Cube
import argparse
import os
import re
import numpy as np
import sys
import matplotlib.pyplot as plt
import logging
from regions import read_ds9

# read arguments
ap = argparse.ArgumentParser(description='Parameters for spectrum extraction of a subcube region of the input cube')
ap.add_argument("input_files", nargs='+', help="List of fits muse data cubes to be loaded")
ap.add_argument("-r", "--region", nargs=1, help="Region file with the regions to be extracted", type=str)
ap.add_argument("-o", "--outdir", nargs='?', help="Name of the output directory", default="subcubes", type=str)
ap.add_argument('-s', "--sum", help="Switch to indicate whether the output spectrum should be summed or averaged", action='store_true')

args = ap.parse_args()

muse_cubes = args.input_files
region_file = args.region[0]
outdir = args.outdir
sum = args.sum

scriptname = os.path.basename(__file__)

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(message)s', datefmt='%d-%b-%y %H:%M:%S')
logger = logging.getLogger(scriptname)

if os.path.isfile(region_file):
    regs = read_ds9(region_file)
    logger.info("Loaded %i region(s) from %s" % (len(regs), region_file))
    logger.info(regs)
else:
    logger.error("Region file %s not found" % region_file)
    sys.exit()

if not os.path.isdir(outdir):
    os.mkdir(outdir)


for cubefile in muse_cubes:

    if os.path.isfile(cubefile):
        logger.debug('Loading cube %s ...' % cubefile)
        cube = Cube(cubefile)
    else:
        logger.warning("Cube %s not found. Skipping..." % cubefile)
        continue

    for reg, reg_index in zip(regs, np.arange(0, len(regs))):
        plt.figure()
        if type(regs[0]).__name__ != 'CircleSkyRegion':
            logger.warning("Skipping region %s. Only circular regions are accepted for now.")
            continue
        else:
            centerra = reg.center.fk5.ra.value
            centerdec = reg.center.fk5.dec.value
            subcube_radius = reg.radius.value

            logger.info("Extracting spectrum number %i from cube %s around position (%s) with \
            a radius of %.2f arcsec" % (reg_index, cubefile, (centerra, centerdec), subcube_radius))
            try:
                subcube = cube.subcube_circle_aperture((centerdec, centerra), subcube_radius)
            except Exception as e:
                logging.error("Exception occurred", exc_info=True)
                continue
            if sum:
                print("Summing output spectra")
                spe = subcube.sum(axis=(1, 2))
            else:

                print("Averaging output spectra")
                spe = subcube.mean(axis=(1, 2))
            image = subcube.sum(axis=0)

            spe.info()
            spe.plot()
            plt.show()

            text = reg.meta['label']

            outname = re.sub("#.*text.*=\{|\}", "", text)
            if outname == "":
                outname = 'subcube%i' % reg_index

            outcubename = cubefile.replace(".fits", "")
            spe.write("%s/%s%s_spec.fits" % (outdir, outcubename, outname))
            image.write("%s/%s%s_img.fits" % (outdir, outcubename, outname))
            subcube.write("%s/%s%s_subcube.fits" % (outdir, outcubename, outname))
    cube = None
