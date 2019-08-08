# -*- coding: utf-8 -*-
# Script to find stars in an image
# Author: Andres Gurpide Lasheras  andres.gurpide@gmail.com
# 1-06-2019
# !/usr/bin/env python3
# imports

import logging

from photutils import make_source_mask
import argparse
from astropy.stats import sigma_clipped_stats
from photutils import DAOStarFinder
from mpdaf.obj import Image
import os
from astropy.table import Column
from photutils import Background2D, MedianBackground
from astropy.stats import SigmaClip
import astropy.units as u
from muse_utils import region_to_mask
import numpy as np
import matplotlib.pyplot as plt


# ---------------------------------------------------------------------------------
# main
# read arguments
parser = argparse.ArgumentParser(description='Parameters to find sources in the input image.')
parser.add_argument("images", help="Image files where to look for sources", nargs='+', type=str)
parser.add_argument("-s", "--sigma", type=float, help='Sigma level for the sigma clipping filter to estimate the background', default=3, nargs='?')
parser.add_argument("-t", "--threshold", help='Sigma threshold to consider a detection', type=float, nargs='?', default=5)
parser.add_argument("-f", "--fwhm", help='FWHM in arcseconds for the star detection', type=float, nargs='?', default=1.0)
parser.add_argument("-m", "--mask", help='Optional mask region to mask certain parts of the image (include/exclude). Alternatively a fits file with 0 (mask) and 1 not mask can be provided too.', type=str, nargs='?', default=None)

args = parser.parse_args()

show_mask = 0

sigma_clip = args.sigma

sigma_threshold = args.threshold

fwhm = args.fwhm * u.arcsec

masked_region = args.mask

brightest = 6

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

for image_file in args.images:
    if os.path.isfile(image_file):
        # load image
        image = Image(image_file)
        # MUSE image
        if 'HIERARCH ESO OCS IPS PIXSCALE' in image.primary_header:
            fwhm_pixels = fwhm.value / image.primary_header['HIERARCH ESO OCS IPS PIXSCALE']
        elif "D005SCAL" in image.primary_header:
            fwhm_pixels = fwhm.value / image.primary_header['D005SCAL']

        logger.info("FWHM = %.1f pixels " % fwhm_pixels)

        # mask most prominent sources for background computation
        mask = make_source_mask(image.data, snr=2, filter_fwhm=int(fwhm_pixels), filter_size=int(3 * fwhm_pixels), dilate_size=int(2 * fwhm_pixels), npixels=int(fwhm_pixels))

        nx = image.wcs.naxis1
        ny = image.wcs.naxis2

        sigma_clipping = SigmaClip(sigma=sigma_clip)

        bkg_box_size = (image.wcs.naxis1 // 6, image.wcs.naxis2 // 6)

        logger.info("Box size for the 2D background calculation pixels: \n")
        logger.info(bkg_box_size)
        bkg_estimator = MedianBackground()
        bkg = Background2D(image.data, bkg_box_size, sigma_clip=sigma_clipping)  # uses sextractor background estimator by default

        fig, axs = plt.subplots(nrows=2, ncols=2)
        titles = ["Mask", "Background", "Data", "Data - Background"]
        data = [mask, bkg.background, image.data, image.data - bkg.background]
        
        for title, ax, dataset in zip(titles, axs.flat, data):
            ax.imshow(dataset, origin="lower", cmap='inferno')
            ax.set_title(title)
        fig.savefig("daofind.pdf")

        if show_mask:
            plt.show()

        # compute background sigma clipping the masked image
        # mean, median, std = sigma_clipped_stats(image.data, sigma=sigma_clip, mask=mask)
        mean, median, std = sigma_clipped_stats(image.data, sigma=sigma_clip)
        # compute the standard deviation of the original image

        daofind_brightest = DAOStarFinder(fwhm=fwhm_pixels, threshold=sigma_threshold * std, brightest=brightest, exclude_border=True, sigma_radius=3)
        # apply mask if provided
        if masked_region is not None and 'reg' in masked_region:
            logger.info("Masking region inside %s" % masked_region)
            mask = region_to_mask(image, masked_region)
        elif masked_region is not None:
            mask = masked_region
        else:
            mask = None

        sources_brightest = daofind_brightest(image.data - bkg.background, mask=mask)

        daofind_all = DAOStarFinder(fwhm=fwhm_pixels, threshold=sigma_threshold * std, exclude_border=True, sigma_radius=3)

        sources_all = daofind_all(image.data - bkg.background, mask=mask)

        yx = np.stack((sources_all['ycentroid'].data, sources_all['xcentroid'].data)).T

        dec_ra_values = image.wcs.pix2sky(yx)

        sources_all.add_column(Column(name='ra', data=dec_ra_values.T[1], unit=u.degree), name='ra')
        sources_all.add_column(Column(name='dec', data=dec_ra_values.T[0], unit=u.degree), name='dec')
        sources_all.add_column(Column(name='err_ra', data=np.ones(dec_ra_values.T[0].shape) * fwhm), name='err_ra')
        sources_all.add_column(Column(name='err_dec', data=np.ones(dec_ra_values.T[0].shape) * fwhm), name='err_dec')

        yx = np.stack((sources_brightest['ycentroid'].data, sources_brightest['xcentroid'].data)).T

        dec_ra_values = image.wcs.pix2sky(yx)

        sources_brightest.add_column(Column(name='ra', data=dec_ra_values.T[1], unit=u.degree), name='ra')
        sources_brightest.add_column(Column(name='dec', data=dec_ra_values.T[0], unit=u.degree), name='dec')
        sources_brightest.add_column(Column(name='err_ra', data=np.ones(dec_ra_values.T[0].shape) * fwhm, unit=u.arcsec), name='err_ra')
        sources_brightest.add_column(Column(name='err_dec', data=np.ones(dec_ra_values.T[0].shape) * fwhm, unit=u.arcsec), name='err_dec')

        logger.info("Writing out results...")
        out_table_all = "allstars%s" % (os.path.basename(image_file).replace(".fits", ".csv"))
        out_table_brightest = "brighteststars%s" % (os.path.basename(image_file).replace(".fits", ".csv"))
        logger.info("All stars saved to %s" % out_table_all)
        sources_all.write(out_table_all, overwrite=True, format='csv', delimiter='\t')
        logger.info("Found %i sources" % len(sources_all))
        logger.info(sources_all)
        logger.info("%i brightest stars saved to %s" % (len(sources_brightest), out_table_brightest))
        sources_brightest.write(out_table_brightest, overwrite=True, format='csv', delimiter='\t')
