# -*- coding: utf-8 -*
# Determines the FWHM polynomial parameters from a fits file resulting from determine_fwhm
# Author: Andres Gurpide Lasheras  andres.gurpide@gmail.com
# 28-05-2019
# !/usr/bin/env python3

import argparse
import os
import matplotlib.pyplot as plt
from lmfit.models import MoffatModel, GaussianModel
from astropy.io import fits
import logging


# read arguments
ap = argparse.ArgumentParser(description='Retrieve the FWHM as a function of wavelenght from a fits file output of determine_fwhm.py')
ap.add_argument("input_fits", nargs='+', help="List of fits files containing as header the wavelenght as two columns with radius and flux fractions.")
ap.add_argument('-p', '--profile', type=str, choices=["moffat", "gaussian"], default="moffat")

args = ap.parse_args()

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

peak_profile = 'fwhm_model_'

if args.profile == 'moffat':
    peak_model = MoffatModel(prefix=peak_profile)
elif args.profile == 'gaussian':
    peak_model = GaussianModel(prefix=peak_profile)

peak_model.make_params(amplitude=1, center=0, sigma=1.5, beta=1)
peak_model.set_param_hint("%scenter" % peak_profile, value=0, vary=False)

fwhms = []
fwhm_err = []
waves = []

for fits_file in args.input_fits:
    with fits.open(fits_file, mode='append') as hdul:
        # skip header file
        for hdu in hdul[1:]:

            waves.append(float(hdu.header['EXTNAME']))
            output = peak_model.fit(hdu.data.field(1), x=hdu.data.field(0))
            if not output.errorbars:
                logger.warning("Fit not successful for wavelenght %s" % hdu.header['EXTNAME'])
                continue

            best_fit_params = output.params

            try:
                hdu.columns.add_col(fits.Column(array=output.best_fit, name='%s fit' % args.profile, unit='', format='E'))
            except ValueError:
                hdu.columns['%s fit' % args.profile].array = output.best_fit
            fwhms.append(best_fit_params.get("%sfwhm" % peak_profile).value)
            fwhm_err.append(best_fit_params.get("%sfwhm" % peak_profile).stderr)
            plt.plot(hdu.data.field(1), output.best_fit)
            plt.plot(hdu.data.field(1), output.init_fit)

        hdul.flush()
    plt.show()
