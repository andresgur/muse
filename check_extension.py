# -*- coding: utf-8 -*
# Script to check if a source is extended or not. Needs PSF file determined with determine_psf script.
# Author: Andres Gurpide Lasheras  andres.gurpide@gmail.com
# 25-06-2019
# !/usr/bin/env python3

# imports
    from mpdaf.obj import Image, Moffat2D
    import matplotlib.pyplot as plt
    import argparse
    import os
    import muse_utils as mu

    import logging
    import pyregion
    import numpy as np
    from astropy.io import fits
    import mpdaf.MUSE.PSF as psf
    import glob


# ------------------------------
# main
# read arguments
ap = argparse.ArgumentParser(description='Determines the FWHM of the cube and/or the energy encircled in an image')
ap.add_argument("-d", "--camel_dir", nargs=1, help="Directory where the images are stored", type=str)
ap.add_argument("-l", "--line", nargs=1, help="Line name to check for extension", type=str)
ap.add_argument('-i', '--intensity_image', nargs=1, type=str)
ap.add_argument("-psf", "--psf_fits", nargs=1, help="PSF fits file with PSF information")
ap.add_argument("-e", "--extension", nargs='?', help="Extension to read the PSF from the fits file", default=1, type=int)
ap.add_argument("-r", "--source_region", nargs=1, help="Source selecting the region to check for extended emission.", default=None, type=str)
args = ap.parse_args()

# logger
scriptname = os.path.basename(__file__)
logger = logging.getLogger(scriptname)
logger.setLevel(logging.DEBUG)
out_format = logging.Formatter('%(name)s - %(levelname)s - %(message)s')

# handler
stream_handler = logging.StreamHandler()
stream_handler.setLevel(logging.DEBUG)
stream_handler.setFormatter(out_format)
logger.addHandler(stream_handler)

camel_dir = args.camel_dir[0]
line = args.line[0]

walength_file = glob.glob('%s/*[!e]wave*%s.fits' % (camel_dir, line))[0]
ewalength_file = glob.glob('%s/*ewave*%s.fits' % (camel_dir, line))[0]
logger.info("Wavelength files found:")
print(walength_file)
print(ewalength_file)
walength_map = Image(walength_file)
ewalength_map = Image(walength_file)

intensity_file = glob.glob('%s/*[!e]int*%s.fits' % (camel_dir, line))[0]

eintensity_file = glob.glob('%s/*eint*%s.fits' % (camel_dir, line))[0]
logger.info("Intensity files found:")
print(intensity_file)
print(eintensity_file)

intensity_map = Image(intensity_file)
eintensity_map = Image(eintensity_file)


source_region = args.source_region[0]
out_dir = 'check_extension'

if source_region is not None:
    mask = mu.region_to_mask(intensity_map, source_region)
    walength_map.mask[np.where(mask)] = True
    intensity_map.mask[np.where(mask)] = True
    eintensity_map.mask[np.where(mask)] = True
    intensity_map.crop()
    eintensity_map.crop()
    walength_map.crop()
    region = pyregion.open(source_region)[0]
    center = (region.coord_list[1], region.coord_list[0])

# get mean observed wavelength
mean_wavelength = np.mean(walength_map.data)

logger.info("Mean wavelength: %.1f" % mean_wavelength)

psf_table, header, headername = mu.read_psf_fits(args.psf_fits[0], args.extension)

# create output dir
out_dir = "%s/%s" % (out_dir, headername)
if not os.path.isdir(out_dir):
    os.mkdir("%s" % out_dir)

# find nearest wavelength in the wavelength array
wa_idx = (np.abs(psf_table["Wavelength"] - mean_wavelength)).argmin()

if center is None:
    center = walength_map.wcs.pix2sky((walength_map.shape[0] // 2, walength_map.shape[1] // 2))[0]  # returns degrees

# determine peak
peak = intensity_map.peak(center=center)
logger.info("Image peak:")
print(peak)
if "Beta" in header:
    avg_beta = float(header["Beta"])
else:
    avg_beta = 2

logger.info("Average beta index used for the PSF creation: %.1f" % avg_beta)
psf_cube = psf.create_psf_cube((psf_table["Wavelength"].size, intensity_map.shape[0], intensity_map.shape[1]), psf_table["FWHM"], avg_beta, wcs=intensity_map.wcs, unit_fwhm=psf_table.columns["FWHM"].unit)

psf_moffat = Image(wcs=intensity_map.wcs, data=psf_cube[wa_idx], cmap='inferno')
psf_moffat.mask = walength_map.mask

# PSF vs real image
figure_comparison, contours_ax = plt.subplots(nrows=2, ncols=1, figsize=(16.0, 10.0), subplot_kw={'projection': intensity_map.wcs.wcs}, sharex=True)

intensity_map.plot(ax=contours_ax[0], zscale=True, colorbar='v', cmap='inferno')

im_0 = contours_ax[0].images
cb = im_0[-1].colorbar
cb.ax.set_ylabel(mu.add_zlabel(intensity_map.filename, intensity_map.unit), fontsize=14)
contours_ax[1].set_xlabel('Ra', fontsize=20)
psf_moffat.plot(ax=contours_ax[1], zscale=True, cmap='inferno', colorbar='v')
contours_ax[1].set_ylabel('Dec', fontsize=20)
contours_ax[0].set_ylabel('Dec', fontsize=20)
im = contours_ax[1].images
cb = im[-1].colorbar
cb.ax.set_ylabel(mu.add_zlabel("moffat_fit", intensity_map.unit), fontsize=14)

#figure_comparison.tight_layout()
figure_comparison.subplots_adjust(hspace=0)

# Prepare plot
figure_eer = plt.figure(figsize=(16.0, 10.0))
eer_ax = figure_eer.add_subplot(111)

eer_ax.set_xlabel("Radius (%s)" % psf_table.columns["FWHM"].unit, fontsize=20)
eer_ax.set_ylabel("Intensity fraction", fontsize=20)
eer_ax.ticklabel_format(axis='both', style='plain')
mu.format_axis(eer_ax)

radius_moffat, energy_values_moffat = psf_moffat.eer_curve(cont=0)
radius, energy_values = intensity_map.eer_curve(cont=0, center=(peak['y'], peak['x']), unit_center='deg', unit_radius='arcsec')

# errors
error_int_map = Image(data=intensity_map.data + eintensity_map.data, wcs=intensity_map.wcs)
radius_err, energy_values_err = eintensity_map.eer_curve(cont=0, center=(peak['y'], peak['x']), unit_center='deg', unit_radius='arcsec')

eer_ax.errorbar(x=radius_moffat, y=energy_values_moffat, color='green', label='PSF')
# yerr=energy_values - energy_values_err
eer_ax.errorbar(x=radius, y=energy_values, color='blue', label='source')
figure_eer.tight_layout()
eer_ax.legend(loc='best', prop={'size': 14}, frameon=False)
figure_eer.savefig("%s/eer_%s" % (out_dir, os.path.basename(source_region).replace("reg", "pdf")))

# subtract PSF to image
res = Image(data=intensity_map.data / intensity_map.data.sum() - psf_moffat.data, wcs=intensity_map.wcs)
res.primary_header['PSF'] = "PSF used %s" % args.psf_fits[0]
res.primary_header['REGION'] = "Region used %s" % source_region
res.write("%s/%sres%s" % (out_dir, line, os.path.basename(source_region).replace("reg", "fits")))
psf_moffat.write("%s/%spsf%s" % (out_dir, line, os.path.basename(source_region).replace("reg", "fits")))
logger.info("Output written in %s" % out_dir)
plt.show()
