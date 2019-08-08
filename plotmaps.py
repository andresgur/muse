# !/usr/bin/env python3
# -*- coding: utf-8 -*
# Script to plot line maps
# Author: Andres Gurpide Lasheras  andres.gurpide@gmail.com
# 5-07-2019





# imports
import argparse
import os
import sys
import logging
import muse_utils as muse_utils
from scipy.ndimage.filters import gaussian_filter
from glob import glob
import pyregion
from mpdaf.obj import Image
import matplotlib.pyplot as plt
import logging
import muse_utils as mu


# read arguments
ap = argparse.ArgumentParser(description='To run in the directory containing the cleaned maps')
ap.add_argument("--hst", nargs='?', help="Path to HST image", default=None)
ap.add_argument("-r", "--region", nargs='?', help="Path to the region to overlay", default=None)

args = ap.parse_args()

region = args.region

# logger
scriptname = os.path.basename(__file__)

logger = logging.getLogger(scriptname)
logger.setLevel(logging.DEBUG)
out_format = logging.Formatter('%(name)s - %(levelname)s - %(message)s')

# line maps keywors
map_keywords = {'flux': '*[!e]flux*', 'disp': '*[!e]disp*', 'vel': '*[!e]vel*', 'cont': 'white*'}

lines_ha = ['HALPHA', 'NII6583', 'SII6716', 'SII6731']
directory = '.'
#lines_hb = ['HBETA', 'OIII5007']
#dir_hb = args.hbetapath

# line maps
flux_files = []

flux_files = [glob("%s/cleaned_images/%s%s*.fits" % (directory, map_keywords['flux'], line))[0] for line in lines_ha]

logger.info("Flux files:")
print(flux_files)

disp_files = [glob("%s/cleaned_images/%s%s*.fits" % (directory, map_keywords['flux'], line))[0] for line in lines_ha]
logger.info("Dispersion files:")
print(disp_files)


# velocity maps and continuum maps
vel_maps = glob(directory + "/cleaned_images/" + map_keywords['vel'])
cont_maps = glob(directory + "/cleaned_images/" + map_keywords['cont'])

nrows = len(flux_files)
ncols = 3

flux_images = [Image(flux_map).gaussian_filter(1, interp='no', inplace=True, truncate=2) for flux_map in flux_files]
disp_images = [Image(disp_map).gaussian_filter(1, interp='no', inplace=True, truncate=2) for disp_map in disp_files]
vel_images = [Image(vel_map).gaussian_filter(1, interp='no', inplace=True, truncate=2) for vel_map in vel_maps]
cont_images = [Image(cont_map).gaussian_filter(1, interp='no', inplace=True, truncate=2) for cont_map in cont_maps]

img_figure, axes = plt.subplots(nrows, ncols, figsize=(16, 10), sharey=True, subplot_kw={'projection': flux_images[0].wcs.wcs})
img_figure.tight_layout()
img_figure.subplots_adjust(hspace=0)
img_figure.subplots_adjust(vspace=0)

for flux_image, ax in zip(flux_images, axes.flat):
    flux_image.plot(ax=ax, scale='linear', colorbar='h')
    #if region is not None:
    #    mu.plot_regions(region, ax)
    im = ax.images
    cb = im[-1].colorbar
    cb.ax.set_ylabel(mu.add_zlabel(flux_image.filename, flux_image.unit), fontsize=14)

for flux_image, ax in zip(flux_images, axes.flat):
    flux_image.plot(ax=ax, scale='linear', colorbar='h')
    #if region is not None:
    #    mu.plot_regions(region, ax)
    im = ax.images
    cb = im[-1].colorbar
    cb.ax.set_ylabel(mu.add_zlabel(flux_image.filename, flux_image.unit), fontsize=14)
plt.show()
