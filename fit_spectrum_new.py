# -*- coding: utf-8 -*
# Script to load, fit and plot spectra generated from MUSE cubes
# Author: Andres Gurpide Lasheras  andres.gurpide@gmail.com
# 25-03-2019
# !/usr/bin/env python3

# imports
from mpdaf.obj import Spectrum
import matplotlib.pyplot as plt
import argparse
import os
import numpy as np
from astropy.stats import gaussian_fwhm_to_sigma
import logging
from specutils import Spectrum1D
from specutils.fitting import find_lines_threshold
from mpdaf.sdetect import matchlines, get_emlines
from astropy.stats import sigma_clipped_stats
from astropy.nddata.nduncertainty import StdDevUncertainty
from mpdaf.sdetect.source import emlines
from specutils.manipulation import noise_region_uncertainty
from astropy.table import Table
from astropy.table import Column
from math import pi, sqrt
from lmfit.models import LorentzianModel, PolynomialModel, GaussianModel
import sys


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

# read arguments
ap = argparse.ArgumentParser(description='Spectrum fits to be loaded')
ap.add_argument("input_files", nargs='+', help="List of spectra to be loaded")
ap.add_argument("-l", "--line_catalog", nargs='?', help="Line catalog of spectral lines in fits format",
                default=os.environ['HOME'] + "/scripts/pythonscripts/muse/optical_line_catalog/illss.fits")
ap.add_argument("-z", "--redshift", nargs='?', help="Redshift of the object", default=0, type=float)
ap.add_argument("-t", "--threshold", nargs='?', help="SNR to consider a line significant after fitting", default=5, type=float)
ap.add_argument("-k", "--kernel_size", nargs='?', help="Kernel size (A) of the median filter for the spike masking of the continuum fitting", default=20.0, type=float)
ap.add_argument("-d", "--degree", nargs='?', help="Polynomial degree for the continuum fitting", default=1, type=int)
ap.add_argument("-e", "--epsilon", nargs='?', help="Epsilon for the wavelet smoothing (0,1)", default=0.05, type=float)
ap.add_argument("--sigma_clip", nargs='?', help="Sigma clip for the sigma clipping noise calculation", default=3.0, type=float)

args = ap.parse_args()

spectra = args.input_files
redshift = args.redshift
line_threshold = args.threshold
kernel_size = args.kernel_size
line_catalog = args.line_catalog

# lines to be fitted
line_groups = np.array([np.array(['HBETA4861', 'Fe II', '[OIII]4959', '[OIII]5007']), np.array(['FEVI5722']), np.array(['[HeI]5877']), np.array(["[OI]6302", "[OI]6365"]), np.array(['[NII]6548', 'HALPHA6563', '[NII]6583']), np.array(['[SII]6716', '[SII]6731']), np.array(['[Arɪɪɪ]'])])
threshold = 50
# create color array

fig, axs = plt.subplots(nrows=2, ncols=1)

table_line = Table.read(line_catalog)  # or other format

filter = ((table_line['n_lambda'] != '?') & (table_line['Intens'] == 0) & (table_line['lambda'] < 9500) & (table_line['lambda'] > 4500))

table_line = table_line[filter]

lines = dict(zip(table_line['lambda'], table_line['Element']))

lines = emlines

for spectrum in spectra:
    if os.path.isfile(spectrum):
        logger.info('Loading spectrum %s ...' % spectrum)
        # read data and variance
        spe = Spectrum(filename=spectrum)
        spe.info()
        smoothed_spe = spe.median_filter(kernel_size, spline=True, unit=spe.wave.unit, inplace=False)

        spe.plot(ax=axs.flat[0])

        smoothed_spe.plot(ax=axs.flat[0], label='Median filtered', color='green')

        cont_poly = smoothed_spe.poly_fit(args.degree, weight=True, verbose=False)

        continuum_spe = smoothed_spe.clone()

        continuum_spe.poly_val(cont_poly)  # create spectrum with continuum coefficient of the poly fit

        continuum_spe.plot(ax=axs.flat[0], label='Cont fit', color='red')

        norm_spec = spe / continuum_spe

        norm_spec.plot(ax=axs.flat[1])

        clean_normed_spec = norm_spec.wavelet_filter(epsilon=args.epsilon, inplace=False)

        clean_normed_spec.plot(ax=axs.flat[1], label='Wavelet smoothed')

        cont_sub = spe - continuum_spe

        clean_sub = cont_sub.wavelet_filter(epsilon=args.epsilon, inplace=False)

        # use the standard deviation in the smoothed spectrum to find initial line guesses
        mean, median, stddev = sigma_clipped_stats(clean_sub.data, sigma=args.sigma_clip)

        uncertainty = StdDevUncertainty(np.ones(clean_sub.shape) * stddev * clean_sub.unit)

        spec1D = Spectrum1D(flux=clean_sub.data * clean_sub.unit, spectral_axis=clean_sub.wave.coord() * clean_sub.wave.unit,
                            uncertainty=uncertainty)

        lines_table = find_lines_threshold(spec1D, noise_factor=5)

        em_lines = lines_table[lines_table['line_type'] == 'emission']
        # returns overall error and indexes of the lines found in the line catalog
        error, indexes = matchlines(len(em_lines), em_lines['line_center'].value, redshift, lines)

        wavelengths_catalog = np.array(list(lines.keys()))[indexes]

        wa_diff = abs(wavelengths_catalog - em_lines['line_center'].value / (1 + redshift))

        indexes = indexes[np.where(wa_diff < threshold)]

        logger.info("List of identified lines (rest frame positions): \n")

        found_line_names = np.array(list(lines.values()))[indexes]

        print(np.array(list(lines.items()))[indexes])

        fwhm_est = []

        for row_index, line_index in enumerate(lines_table['line_center_index']):
            try:
                fwhm_est.append(spe.fwhm(spe.wave.coord()[line_index], cont=continuum_spe[line_index], unit=spe.wave.unit))
            except ValueError:
                logger.warning("FWHM could not be estimated. Removing inital guessed line at %.2f" % spe.wave.coord()[line_index])
                lines_table.remove_row(row_index)

        lines_table.add_column(Column(data=fwhm_est, format='E', name='FWHM'))

        # set continuum to model retrieved before

        complete_model = PolynomialModel(degree=args.degree, prefix='cont_')
        for paramname, c in zip(complete_model.param_names, cont_poly):

            complete_model.set_param_hint(paramname, value=c, vary=False)

        # set up parameter hints
        for row in lines_table:
            line_prefix = "%s%.0f_" % (row['line_type'], row['line_center'].value)
            line_model = GaussianModel(prefix=line_prefix)

            if row['line_type'] == 'emission':
                complete_model += line_model
            elif row['line_type'] == 'absorption':
                complete_model -= line_model
            # set initial guesses for Gaussian line
            for paramname in line_model.param_names:
                if 'center' in paramname:
                    value = float(row['line_center'].value)
                    min = value - 10
                    max = value + 10
                elif 'sigma' in paramname:
                    min = 1.0 * gaussian_fwhm_to_sigma
                    value = row['FWHM'] * gaussian_fwhm_to_sigma
                    max = 20.0 * gaussian_fwhm_to_sigma
                elif 'amplitude' in paramname:
                    index = spe.wave.pixel(row['line_center'].value, nearest=True, unit=spe.wave.unit)
                    peak = spe.data[index]
                    continuum = continuum_spe.data[index]
                    value = peak * sqrt(2 * pi) + continuum
                    min = 0
                    max = None
                else:
                    continue
                print(paramname)
                print(value)
                print(max)
                print(min)
                complete_model.set_param_hint(paramname, value=value, vary=True, max=max, min=min)

        # fit spectrum
        if spe.var is not None:
            weights = 1.0 / np.sqrt(spe.var)
            np.ma.fix_invalid(weights, copy=False, fill_value=1)

        else:
            logger.warning("Variances are null, using 1 as weights for the curve fitting")
            weights = np.ones(spe.shape)

        output = complete_model.fit(spe.data, x=spe.wave.coord(), weights=weights)
        best_fit_params = output.params
        mod_components = output.eval_components(x=spe.wave.coord())
        print(output.best_fit)

        logger.debug(output.fit_report())

        # iterate over the lines fitted
        real_lines = []
        plt.figure()

        for res_center, line_type in zip(lines_table['line_center'], lines_table["line_type"]):
            line_prefix = "%s%.0f_" % (line_type, res_center.value)

            line_center = best_fit_params.get("%scenter" % line_prefix)
            line_peak = best_fit_params.get("%sheight" % line_prefix)
            fwhm = best_fit_params.get("%sfwhm" % line_prefix)
            line_flux = best_fit_params.get("%samplitude" % line_prefix)
            # get pixels with 2sigma around the line
            line_sigma = fwhm * gaussian_fwhm_to_sigma
            line_low2sigma = spe.wave.pixel(line_center + 2 * (line_sigma), nearest=True, unit=spe.wave.unit)
            line_high2sigma = spe.wave.pixel(line_center - 2 * (line_sigma), nearest=True, unit=spe.wave.unit)
            line_spe = spe[line_low2sigma:line_high2sigma]
            print(line_spe)

            line_noise = np.sum(abs(line_spe.var)) / line_spe.shape

            snr = abs(line_flux) / abs(line_noise)

            if snr < line_threshold:
                logger.info("Line at %.1f has low SNR" % line_center)
                print(line_sigma)
                print(line_peak)
                print(fwhm)
                print(line_prefix)

                print(line_prefix)

            line_spe.plot()
            plt.plot(line_spe.wave.coord(), mod_components[line_prefix])

            plt.show()









        '''
        for line_group in line_groups:
            for line in line_group:
                indexes_found = np.where(found_lines == line)
                if indexes_found.shape == 0:
                    logger.info("Line %s not found" % line)
                    # line not found
                elif indexes_found > 1:
                    wa_diff = (spe[indexes_found] / (1 + redshift) - lbdas[np.newaxis, :])
                    jfound = np.argmin(a, axis=1)
        '''





        [axs.flat[1].text(center, 10, line_name.replace(" ", "")) for center, line_name in zip(np.array(list(lines.keys()))[indexes] / (1 + redshift), np.array(list(lines.values()))[indexes])]


for ax in axs.flat:
    ax.legend(loc='best')

plt.show()
