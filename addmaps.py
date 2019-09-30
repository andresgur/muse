## Script to add maps together, useful for line doublets maps


# imports
import argparse
import os
import sys
import numpy as np
import logging
from astropy.io import fits


# read arguments
ap = argparse.ArgumentParser(description='Spectrum fits to be loaded')
ap.add_argument("maps", nargs='+', help="The line maps to be added in the numerator of the line ratio")
ap.add_argument("-o", "--outdir", nargs='?', help="Output dir", default='addedmaps')
ap.add_argument("-f", "--file", nargs='?', help="Output file name (without fits ending)", default='addedmap.fits')
args = ap.parse_args()


linemaps = args.maps
fits_linemaps = [fits.open(linemap) for linemap in linemaps if os.path.isfile(linemap)]
added_data = np.sum([fitsmap["DATA"].data for fitsmap in fits_linemaps], axis=0)
added_map = fits.PrimaryHDU(data=added_data, header=fits_linemaps[0][0].header)
added_map.header['MAP_SUM'] = "Added %s line maps togeter" % (str(linemaps).strip('[]'))
if not os.path.isdir(args.outdir):
    os.mkdir(args.outdir)
added_map.writeto("%s/%s" % (args.outdir, args.file), overwrite=True)
