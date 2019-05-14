# -*- coding: utf-8 -*
# Script to load, fit and plot spectra generated from MUSE cubes
# Author: Andres Gurpide Lasheras  andres.gurpide@gmail.com
# 09-05-2019
# !/usr/bin/env python3

# imports
from mpdaf.obj import Image
import matplotlib.pyplot as plt
import argparse
from regions import read_ds9
import os
from astropy.io import fits
import sys
import logging


scriptname = os.path.basename(__file__)
# logger

logging.basicConfig(level=logging.DEBUG, format='%(message)s')
logger = logging.getLogger(scriptname)

# read arguments
ap = argparse.ArgumentParser(description='Spectrum fits to be loaded')
ap.add_argument("input_maps", nargs='+', help="List of images to be fitted")
ap.add_argument("-e", "--error_maps", nargs='*', help="List error images for the fitting process")
ap.add_argument("-c", "--cut_region", nargs='?', help="Region for the spatial cut", default=None, type=str)

args = ap.parse_args()

input_maps = args.input_maps
cutting_region = args.cut_region
fit_dir = 'image_fit'

plot_out = 1


if not os.path.isdir(fit_dir):
    os.mkdir(fit_dir)


nmaps = len(input_maps)

if args.error_maps is not None:

    ne_maps = len(args.error_maps)
    weights = True
    if nmaps != ne_maps:
        logger.error("Number of maps and error maps has to be equal")
        sys.exit()
    e_maps = args.error_maps

    for fits_map, e_fits_map in zip(input_maps, e_maps):

        if os.path.isfile(fits_map):
            logger.info('Loading image %s ...' % fits_map)

            # read data and variance
            input_map = Image(fits_map)
            e_map = Image(e_fits_map)
            # set the variance of the input map to the values of the error map
            input_map.var = e_map.dat
else:
    weights = False

for processed_map in input_maps:
    if os.path.isfile(processed_map):

        logger.info("Loading %s map for Gaussian fitting" % processed_map)
        image = Image(processed_map)

        if 'flux' in processed_map:
            zlabel = "1e-20 erg/cm$^2$/s"
        else:
            zlabel = ""

        # cut image from region
        if cutting_region is not None:
            regs = read_ds9(cutting_region)
            if len(regs) > 1:
                logging.warning("Only one region will be used for the spatial cut")

            cut_region = regs[0]
            ra = cut_region.center.ra.value
            dec = cut_region.center.dec.value
            radius = cut_region.radius.value

            logging.debug("Cutting the image at (ra, dec) = (%.3f, %.3f) around %.2f asec" % (ra, dec, radius))
            image = image.subimage(center=(dec, ra), size=radius, minsize=0.1)

            center = (dec, ra)

        gauss = image.gauss_fit(center=center, cont=0, fit_back=False, weight=weights, full_output=True, plot=False)
        moffat = image.moffat_fit(center=center, cont=0, fit_back=False, weight=weights, full_output=True, plot=False)

        if plot_out:
            img_figure = plt.figure()
            ax = img_figure.add_subplot(111, projection=image.wcs.wcs)
            image.plot(ax=ax, scale='linear', colorbar='v', show_xlabel=True, show_ylabel=True)
            ax.set_xlabel('Ra')
            ax.set_ylabel('Dec')
            im = ax.images
            cb = im[-1].colorbar
            cb.ax.set_ylabel(zlabel, fontsize=14)
            gauss_figure = plt.figure()
            gauss_ax = gauss_figure.add_subplot(111, projection=image.wcs.wcs)
            gauss.ima.plot(ax=gauss_ax, scale='linear', colorbar='v', show_xlabel=True, show_ylabel=True)
            gauss_ax.set_xlabel('Ra')
            gauss_ax.set_ylabel('Dec')
            im = gauss_ax.images
            cb = im[-1].colorbar
            cb.ax.set_ylabel(zlabel, fontsize=14)
            moffat_figure = plt.figure()
            moffat_ax = moffat_figure.add_subplot(111, projection=image.wcs.wcs)
            moffat.ima.plot(ax=moffat_ax, scale='linear', colorbar='v', show_xlabel=True, show_ylabel=True)
            moffat_ax.set_xlabel('Ra')
            moffat_ax.set_ylabel('Dec')
            im = moffat_ax.images
            cb = im[-1].colorbar
            cb.ax.set_ylabel(zlabel, fontsize=14)

        # store gaussian fit.
        gauss_out_name = processed_map.replace(".fits", "_gauss_fit.fits")
        gauss.ima.write(fit_dir + "/" + gauss_out_name)
        # store moffat fit.
        moffat_out_name = processed_map.replace(".fits", "_moffat_fit.fits")
        moffat.ima.write(fit_dir + "/" + moffat_out_name)

        gauss_res = Image(data=image.data - gauss.ima.data)
        moffat_res = Image(data=image.data - moffat.ima.data)

        res_gauss_out_name = gauss_out_name.replace("fit.fits", "res.fits")
        gauss_res.write(fit_dir + "/" + res_gauss_out_name)

        res_moffat_out_name = moffat_out_name.replace("fit.fits", "res.fits")
        moffat_res.write(fit_dir + "/" + res_moffat_out_name)

        # compute energy enclosed as a function of radius
        radius, energy_values = image.eer_curve(cont=0, center=gauss.center)
        eer_circle_fig = plt.figure()
        ax_eer_circle = eer_circle_fig.add_subplot(111)
        ax_eer_circle.scatter(radius, energy_values, color='green')
        ax_eer_circle.minorticks_on()
        plt.axvline(x=gauss.fwhm[0], label='FWHM - y', ls='--', alpha=0.5)
        plt.axvline(x=gauss.fwhm[1], label='FWHM - x', ls='--', alpha=0.5)
        plt.legend(fontsize=14, loc='best')
        eer_circle_fig.savefig('image_fit/encircled_energy.pdf')
        ax_eer_circle.set_xlabel('Radius (asec)', fontsize=20)
        ax_eer_circle.set_ylabel('EEF', fontsize=20)
        plt.show()

    else:
        logger.warning("Map %s does not exist" % processed_map)
