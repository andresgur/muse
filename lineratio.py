# -*- coding: utf-8 -*
# Script to create line map ratios out of two line map images fits file
# Author: Andres Gurpide Lasheras  andres.gurpide@gmail.com
# 29-04-2019
# !/usr/bin/env python3

# imports
import argparse
import os
import sys
import numpy as np
import logging
from astropy.io import fits

# read arguments
ap = argparse.ArgumentParser(description='Spectrum fits to be loaded')
ap.add_argument("-n", "--numerators", nargs='+', help="The line maps to be added in the numerator of the line ratio")
ap.add_argument("-d", "--denominator", nargs=1, help="The line map for the denominator")
ap.add_argument("-o", "--outdir", nargs='?', help="Output dir", default='lineratios')
ap.add_argument("-f", "--file", nargs='?', help="Output file name (without fits ending)", default='lineratio.fits')
args = ap.parse_args()

linemaps = args.numerators
linemap2 = args.denominator[0]
outdir = args.outdir
outname = args.file

if not os.path.isdir(outdir):
    os.mkdir(outdir)

added_maps_log = ' '.join(linemaps)
numerator_data = []
# add line maps
logging.info('Adding %s maps together' % added_maps_log)
for linemap in linemaps:
    if os.path.isfile(linemap):
        linemapfits = fits.open(linemap)
        if len(numerator_data) == 0:
            numerator_data = linemapfits[0].data
        else:
            numerator_data += linemapfits[0].data
    else:
        logging.warning("Line map %s not found." % linemap)
        continue

if len(numerator_data) == 0:
    logging.error("No line maps were found")
    sys.exit()

if os.path.isfile(linemap2):
    linemap2fits = fits.open(linemap2)
    map2data = linemap2fits[0].data
else:
    logging.warning("Line map %s not found." % linemap2)
    sys.exit()

ratio_fits = fits.PrimaryHDU(data=numerator_data / map2data, header=linemap2fits[0].header)
ratio_fits.data[np.where(map2data is None)] = None

ratio_fits.header['COMMENT'] = "Ratio of %s/%s line maps" % (added_maps_log, linemap2)
if ratio_fits.header['WCSAXES'] == 3:
    ratio_fits.header['WCSAXES'] = 2  # set number of axes 2 for the image

ratio_fits.writeto(outdir + "/" + outname + ".fits", overwrite=True)
print('Line ratio %s/%s written to %s/%s.fits' % (added_maps_log, linemap2, outdir, outname))
