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
from astropy.constants import c
import numpy as np
from astropy.stats import gaussian_fwhm_to_sigma
import logging
from specutils import Spectrum1D
from specutils.fitting import find_lines_threshold
from mpdaf.sdetect import matchlines, get_emlines, crackz
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


# read arguments
ap = argparse.ArgumentParser(description='Spectrum fits to be loaded')
ap.add_argument("input_files", nargs='+', help="List of spectra to be loaded")
ap.add_argument("-l", "--line_catalog", nargs='?', help="Line catalog of spectral lines in fits format",
                default=os.environ['HOME'] + "/scripts/pythonscripts/muse/optical_line_catalog/illss.fits")
ap.add_argument("-z", "--redshift", nargs='?', help="Redshift of the object", default=0, type=float)
ap.add_argument("-n", "--noise_factor", nargs='?', help="Noise factor for initial line finding", default=3, type=float)
ap.add_argument("-t", "--threshold", nargs='?', help="SNR to consider a line significant after fitting", default=3, type=float)
ap.add_argument("-k", "--kernel_size", nargs='?', help="Kernel size (A) of the median filter for the spike masking of the continuum fitting", default=50.0, type=float)
ap.add_argument("-d", "--degree", nargs='?', help="Polynomial degree for the continuum fitting", default=2, type=int)
ap.add_argument("-e", "--epsilon", nargs='?', help="Epsilon for the wavelet smoothing (0,1)", default=0.05, type=float)
ap.add_argument("--lrange", nargs='?', help="Range of wavelenghts to fit (lmin:lmax) in angstrom", default=None, type=str)
ap.add_argument("--sigma_clip", nargs='?', help="Sigma clip for the sigma clipping noise calculation", default=3.0, type=float)

# line constraints from Osterbrock D.E. & Ferland G.J., 2006, Astrophysics of Gaseous Nebulae and Active Galactic Nuclei. University Science Books

args = ap.parse_args()

spectra = args.input_files
redshift = args.redshift
line_threshold = args.threshold
kernel_size = args.kernel_size
line_catalog = args.line_catalog
noise_factor = args.noise_factor
lrange = args.lrange
only_emlines = 1
# lines to be fitted
line_groups = np.array([np.array(['HBETA4861', 'Fe II', '[OIII]4959', '[OIII]5007']), np.array(['FEVI5722']), np.array(['[HeI]5877']), np.array(["[OI]6302", "[OI]6365"]), np.array(['[NII]6548', 'HALPHA6563', '[NII]6583']), np.array(['[SII]6716', '[SII]6731']), np.array(['[Arɪɪɪ]'])])
threshold = 50
# create color array

fig, cont_axs = plt.subplots(nrows=2, ncols=1)

fig, preliminar_fit_ax = plt.subplots()

fig, initial_fit_ax = plt.subplots()

fig, best_fit_ax = plt.subplots()

#fig, final_lines = plt.subplots()

# filter line catalog for the range of interest
table_line = Table.read(line_catalog)  # or other format

filter = ((table_line['n_lambda'] != '?') & (table_line['Intens'] == 0) & (table_line['lambda'] < 9500) & (table_line['lambda'] > 4500))

table_line = table_line[filter]

catalog_lines = dict(zip(table_line['lambda'], table_line['Element']))

catalog_lines = emlines

for spec_index, spectrum in enumerate(spectra):

    if os.path.isfile(spectrum):

        logger.info('Loading spectrum %s ...' % spectrum)
        # read data and variance
        spe = Spectrum(filename=spectrum)
        spe.info()

        if lrange is not None:
            logger.info("Splitting spectrum in wavelenght range %s" % lrange)
            lmin, lmax = lrange.split(":")
            pmin = spe.wave.pixel(float(lmin), nearest=True, unit=spe.wave.unit)
            pmax = spe.wave.pixel(float(lmax), nearest=True, unit=spe.wave.unit)

            spe = spe[pmin: pmax]

        logger.info("Median filtering spectrum for continuum fitting...")

        smoothed_spe = spe.median_filter(kernel_size, spline=True, unit=spe.wave.unit, inplace=False)

        spe.plot(ax=cont_axs.flat[0], label='_nolegend_')

        smoothed_spe.plot(ax=cont_axs.flat[0], label='Median filtered', color='green')

        cont_poly = smoothed_spe.poly_fit(args.degree, weight=True, verbose=False)

        logger.info("Fitting continuum...")

        continuum_spe = smoothed_spe.clone()

        cont_model = PolynomialModel(degree=args.degree, prefix='cont_')

        # fit spectrum
        if spe.var is not None:
            weights = 1.0 / np.sqrt(spe.var)
            np.ma.fix_invalid(weights, copy=False, fill_value=1)

        else:
            logger.warning("Variances are null, using 1 as weights for the curve fitting")
            weights = np.ones(spe.shape)

        for param in cont_model.param_names:
            cont_model.set_param_hint(value=10000, name=param)

        cont_fit = cont_model.fit(smoothed_spe.data, x=smoothed_spe.wave.coord(), cont_c0=10000, weights=weights)

        print(cont_fit.fit_report)
        # store errors
        cont_errors = cont_fit.eval_uncertainty()

        continuum_spe.poly_val(cont_poly)  # create spectrum with continuum coefficient of the poly fit

        continuum_spe.plot(ax=cont_axs.flat[0], label='Cont fit', color='red')

        cont_axs.flat[0].plot(continuum_spe.wave.coord(), cont_fit.best_fit, label='Best fit')

        norm_spec = spe / continuum_spe

        norm_spec.plot(ax=cont_axs.flat[1], label='_nolegend_')

        clean_normed_spec = norm_spec.wavelet_filter(epsilon=args.epsilon, inplace=False, levels=4)

        clean_normed_spec.plot(ax=cont_axs.flat[1], label='Wavelet smoothed')

        cont_sub = spe - continuum_spe

        clean_sub = cont_sub.wavelet_filter(epsilon=args.epsilon, inplace=False, levels=4)

        # use the standard deviation in the smoothed spectrum to find initial line guesses
        mean, median, stddev = sigma_clipped_stats(clean_sub.data, sigma=args.sigma_clip)

        uncertainty = StdDevUncertainty(np.ones(clean_sub.shape) * stddev * clean_sub.unit)

        spec1D = Spectrum1D(flux=clean_sub.data * clean_sub.unit, spectral_axis=clean_sub.wave.coord() * clean_sub.wave.unit,
                            uncertainty=uncertainty)

        lines_table = find_lines_threshold(spec1D, noise_factor=noise_factor)

        lines_table.rename_column('line_center', 'init_center')

        logger.info("Found %i lines" % len(lines_table))

        if only_emlines:
            lines_table = lines_table[np.where(lines_table['line_type'] == 'emission')]
            logger.info("Only %i emission lines" % len(lines_table))

        '''
        em_lines = lines_table[lines_table['line_type'] == 'emission']
        # returns overall error and indexes of the lines found in the line catalog
        error, indexes = matchlines(len(em_lines), em_lines['init_center'].value, redshift, catalog_lines)

        wavelengths_catalog = np.array(list(catalog_lines.keys()))[indexes]

        wa_diff = abs(wavelengths_catalog - em_lines['init_center'].value / (1 + redshift))

        indexes = indexes[np.where(wa_diff < threshold)]

        logger.info("List of identified lines (rest frame positions): \n")

        found_line_names = np.array(list(lines.values()))[indexes]

        print(np.array(list(catalog_lines.items()))[indexes])
        '''

        fwhm_est = []

        rows_to_remove = []

        for row_index, line_index in enumerate(lines_table['line_center_index']):

            try:
                fwhm_est.append(spe.fwhm(spe.wave.coord()[line_index], cont=continuum_spe[line_index], unit=spe.wave.unit))
            except ValueError:
                logger.warning("FWHM could not be estimated. Removing inital guessed line at %.2f" % spe.wave.coord()[line_index])
                rows_to_remove.append(row_index)

        lines_table.add_column(Column(data=fwhm_est, format='E', name='init_FWHM'))

        lines_table.remove_rows(rows_to_remove)

        # set continuum frozen
        for paramname, best_values in zip(cont_model.param_names, cont_fit.best_values):
            cont_model.set_param_hint(paramname, value=cont_fit.best_values[paramname], vary=False)

        complete_model = cont_model
        # set up parameter hints and add lines to the model
        for row in lines_table:

            line_prefix = "%s%.0f_" % (row['line_type'], row['init_center'].value)

            line_model = GaussianModel(prefix=line_prefix)

            # set initial guesses for Gaussian line
            for paramname in line_model.param_names:

                if "center" in paramname:
                    value = row['init_center'].value
                    min = value - 10
                    max = value + 10

                elif 'sigma' in paramname:
                    min = 1.0 * gaussian_fwhm_to_sigma
                    value = row['init_FWHM'] * gaussian_fwhm_to_sigma
                    max = 20.0 * gaussian_fwhm_to_sigma
                elif 'amplitude' in paramname:
                    index = spe.wave.pixel(row['init_center'].value, nearest=True, unit=spe.wave.unit)
                    peak = spe.data[index]
                    continuum = continuum_spe.data[index]
                    value = peak * sqrt(2 * pi) + continuum
                    min = 0
                    max = np.inf
                else:
                    continue

                line_model.set_param_hint(paramname, value=value, vary=True, max=max, min=min)
            if row['line_type'] == 'emission':
                complete_model += line_model
            elif row['line_type'] == 'absorption':
                complete_model -= line_model

        # add weights
        if spe.var is not None:
            weights = 1.0 / np.sqrt(spe.var)
            np.ma.fix_invalid(weights, copy=False, fill_value=1)

        else:
            logger.warning("Variances are null, using 1 as weights for the curve fitting")
            weights = np.ones(spe.shape)

        # perform fit and plot result, filtering low SNR lines
        output = complete_model.fit(spe.data, x=spe.wave.coord(), weights=weights)

        best_fit_params = output.params

        mod_components = output.eval_components(x=spe.wave.coord())

        #logger.debug(output.fit_report())

        # iterate over the lines fitted

        preliminar_fit_ax.set_title("SNR > %.1f" % line_threshold)

        spe.plot(ax=preliminar_fit_ax, label='_nolegend_')
        # model with only the lines above the SNR
        refined_model = cont_model

        logger.info("Filtering low SNR lines...")

        rows_to_remove = []

        for (index, res_center), line_type in zip(enumerate(lines_table['init_center']), lines_table["line_type"]):

            line_prefix = "%s%.0f_" % (line_type, res_center.value)

            line_center = best_fit_params.get("%scenter" % line_prefix)
            line_peak = best_fit_params.get("%sheight" % line_prefix)
            fwhm = best_fit_params.get("%sfwhm" % line_prefix)
            line_flux = best_fit_params.get("%samplitude" % line_prefix)

            # get pixels with 2sigma around the line
            line_sigma = fwhm * gaussian_fwhm_to_sigma
            line_highsigma = spe.wave.pixel(line_center + 3 * (line_sigma), nearest=True, unit=spe.wave.unit)
            line_lowsigma = spe.wave.pixel(line_center - 3 * (line_sigma), nearest=True, unit=spe.wave.unit)

            line_spe = spe[line_lowsigma:line_highsigma]

            line_center_px = spe.wave.pixel(line_center, nearest=True, unit=spe.wave.unit)

            line_cont = continuum_spe.data[line_center_px]

            line_noise = spe.var[line_center_px]

            snr = abs(line_peak - line_cont) / np.sqrt(line_noise + line_cont)

            if snr < line_threshold:
                logger.info("Line at %.1f has low SNR (%.1f)" % (line_center, snr))
                rows_to_remove.append(index)
                continue

            line_model = GaussianModel(prefix=line_prefix)
            # set initial values to best previous model
            for paramname in line_model.param_names:
                line_model.set_param_hint(paramname, value=best_fit_params.get(paramname))

            if line_type == 'absorption':
                preliminar_fit_ax.plot(spe.wave.coord(), - mod_components[line_prefix] + mod_components['cont_'], label="%.1f" % line_center)
                refined_model -= line_model
            else:
                preliminar_fit_ax.plot(spe.wave.coord(), mod_components[line_prefix] + mod_components['cont_'], label="%.1f" % line_center)
                refined_model += line_model

        lines_table.remove_rows(rows_to_remove)

        initial_fit_ax.set_title("Initial fit")
        spe.plot(ax=initial_fit_ax, label='_nolegend_')
        initial_fit_ax.plot(spe.wave.coord(), output.best_fit, color='blue', ls='--')

        final_output = refined_model.fit(spe.data, x=spe.wave.coord(), weights=weights)

        best_fit_params = final_output.params

        mod_components = final_output.eval_components(x=spe.wave.coord())

        best_fit_ax.set_title("Best fit")
        spe.plot(ax=best_fit_ax, label='_nolegend_')
        best_fit_ax.plot(spe.wave.coord(), final_output.best_fit, label="%d" % spec_index)

        #final_lines.set_title("SNR > %.1f" % line_threshold)

        #spe.plot(ax=final_lines)

        fluxes = []
        fit_centers = []
        fwhms = []
        snrs = []
        peaks = []

        for res_center, line_type in zip(lines_table['init_center'], lines_table["line_type"]):

            line_prefix = "%s%.0f_" % (line_type, res_center.value)

            line_center = best_fit_params.get("%scenter" % line_prefix)
            line_peak = best_fit_params.get("%sheight" % line_prefix)
            fwhm = best_fit_params.get("%sfwhm" % line_prefix)
            line_flux = best_fit_params.get("%samplitude" % line_prefix)

            # get pixels with 2sigma around the line
            line_sigma = fwhm * gaussian_fwhm_to_sigma
            line_highsigma = spe.wave.pixel(line_center + 3 * (line_sigma), nearest=True, unit=spe.wave.unit)
            line_lowsigma = spe.wave.pixel(line_center - 3 * (line_sigma), nearest=True, unit=spe.wave.unit)

            line_spe = spe[line_lowsigma:line_highsigma]

            line_center_px = spe.wave.pixel(line_center, nearest=True, unit=spe.wave.unit)

            line_cont = continuum_spe.data[line_center_px]

            line_noise = spe.var[line_center_px]

            snr = abs(line_peak + line_cont) / np.sqrt(line_noise + line_cont)

            if snr < line_threshold:
                logger.info("Line at %.1f has low SNR (%.1f)" % (line_center, snr))
                continue
            '''
            if line_type == 'absorption':
                final_lines.plot(spe.wave.coord(), - mod_components[line_prefix] + mod_components['cont_'], label="%.1f" % line_center)

            else:
                final_lines.plot(spe.wave.coord(), mod_components[line_prefix] + mod_components['cont_'], label="%.1f" % line_center)
            '''
            fluxes.append(line_flux)
            fit_centers.append(line_center)
            fwhms.append(fwhm)
            snrs.append(snr)
            peaks.append(line_peak)

        lines_table.add_column(Column(name='FWHM', format='E', data=fwhms, unit=spe.wave.unit))
        lines_table.add_column(Column(name='Center', format='E', data=fit_centers, unit=spe.wave.unit))
        lines_table.add_column(Column(name='Disp', format='E', data=lines_table['FWHM'] / lines_table["Center"] * c.to('km/s')))
        lines_table.add_column(Column(name='SNR', format='E', data=snrs))
        lines_table.add_column(Column(name='Flux', format='E', data=fluxes, unit=spe.unit))

        em_lines = lines_table[lines_table['line_type'] == 'emission']

        crack_redshift, err_crack_redshift, n_lines_used, wavelenghts, fluxes, identified_lines = crackz(len(em_lines), em_lines['Center'].value, em_lines['Flux'].value,
                                                                                                         catalog_lines, zguess=None)

        logger.info("Redshift from %d lines: %.5f \n Input redsfhit : %.5f \n Difference: %.5f \n Error: %.5f \n Consistent: %r" % (n_lines_used,
                    crack_redshift, redshift, crack_redshift - redshift, err_crack_redshift, crack_redshift - redshift < 2 * err_crack_redshift))
        logger.info("Final list of identified lines:")
        print(identified_lines)
        # dummy method just to obtain the profile names given by the initial line guessed centers
        dummy_z, z_dummy_err, dummy_n, prefix_was, fluxes, identified_lines = crackz(len(em_lines), em_lines['init_center'].value, em_lines['Flux'].value,
                                                                                              catalog_lines, zguess=None)


        # add constraints on the fit based on the identified lines for the final fit
        for paramname, best_values in zip(cont_model.param_names, cont_fit.best_values):
            cont_model.set_param_hint(paramname, value=cont_fit.best_values[paramname], vary=True)
        last_model = cont_model

        for flux, central_wa, line_name, prefix_wa in zip(fluxes, wavelenghts, identified_lines, prefix_was):
            prefix = line_name.replace("[", "").replace("]", "")
            line_model = GaussianModel(prefix=prefix)
            line_model.set_param_hint(name="%scenter" % prefix, value=central_wa, min=central_wa * (1 - err_crack_redshift),
                                      max=central_wa * (1 + err_crack_redshift))
            line_model.set_param_hint(name="%samplitude" % prefix, value=final_output.best_values["emission%.0f_amplitude" % prefix_wa], min=0)
            line_model.set_param_hint(name="%ssigma" % prefix, value=final_output.best_values["emission%.0f_sigma" % prefix_wa], min=1, max=20.0 * gaussian_fwhm_to_sigma)
            last_model += line_model

        oxygen_bigline = '[OIII]5007'.replace("[", "").replace(']', "")
        oxygen_smallline = '[OIII]4959'.replace("[", "").replace(']', "")
        oxygen_line_ratio = 2.99
        nitrogen_bigline = '[NII]6583'.replace("[", "").replace(']', "")
        nitrogen_smallline = '[NII]6548'.replace("[", "").replace(']', "")
        nitrogen_line_ratio = 2.99

        for comp in last_model.components:
            # both lines identified
            if comp.prefix == oxygen_bigline and oxygen_smallline in last_model.name:
                central_wa = 4959
                last_model.set_param_hint(name="%scenter" % oxygen_smallline, value=central_wa * (1 + crack_redshift), min=central_wa * (1 + (err_crack_redshift - redshift)),
                                          max=central_wa * (1 + (crack_redshift + err_crack_redshift)))
                last_model.set_param_hint(name="%samplitude" % oxygen_smallline,
                                          expr="%samplitude / %.2f" % (oxygen_bigline, oxygen_line_ratio), min=0)

            # oxygen line identified but the second one is not
            elif comp.prefix == oxygen_bigline:
                line_model = GaussianModel(prefix=oxygen_smallline)
                central_wa = 4959
                line_model.set_param_hint(name="%scenter" % oxygen_smallline, value=central_wa * (1 + crack_redshift), min=central_wa * (1 + (err_crack_redshift - redshift)),
                                          max=central_wa * (1 + (crack_redshift + err_crack_redshift)))
                line_model.set_param_hint(name="%samplitude" % oxygen_smallline,
                                          expr="%samplitude / %.2f" % (oxygen_bigline, oxygen_line_ratio), min=0)
                last_model += line_model
            # both lines identified
            if comp.prefix == nitrogen_bigline and nitrogen_smallline in last_model.name:
                central_wa = 6548
                last_model.set_param_hint(name="%scenter" % nitrogen_smallline, value=central_wa * (1 + crack_redshift), min=central_wa * (1 + (err_crack_redshift - redshift)),
                                          max=central_wa * (1 + (crack_redshift + err_crack_redshift)))
                last_model.set_param_hint(name="%samplitude" % nitrogen_smallline,
                                          expr="%samplitude / %.2f" % (nitrogen_bigline, nitrogen_line_ratio), min=0)
            # oxygen line identified but the second one is not
            elif comp.prefix == nitrogen_bigline:
                central_wa = 6548
                line_model = GaussianModel(prefix=nitrogen_smallline)
                line_model.set_param_hint(name="%scenter" % nitrogen_smallline, value=central_wa * (1 + crack_redshift), min=central_wa * (1 + (err_crack_redshift - redshift)),
                                          max=central_wa * (1 + (crack_redshift + err_crack_redshift)))
                line_model.set_param_hint(name="%samplitude" % nitrogen_smallline,
                                          expr="%samplitude / %.2f" % (nitrogen_bigline, nitrogen_line_ratio), min=0)
                last_model += line_model

        final_fit_out = last_model.fit(spe.data, x=spe.wave.coord(), weights=weights)

        print(final_fit_out.fit_report())

        fig, final_fitax = plt.subplots()
        spe.plot(ax=final_fitax)
        final_fitax.plot(spe.wave.coord(), final_fit_out.best_fit, label="%d" % spec_index)
        [final_fitax.text(center, 0.7 * flux, line_name) for center, line_name, flux in zip(wavelenghts, identified_lines, fluxes)]
        plt.show()

        if final_fit_out.errorbars:
            uncertainties = final_fit_out.eval_uncertainty()

        '''
        for line, center_value in zip(identified_lines, wavelenghts):
            if line == '[OIII]5007':
                biglineprefix = "%s%.0f_" % ('emission', center_value)
                # if the line was identified, it is in the model
                if '[OIII]4959' in mod_components:
                    center_value = lines_to_wa["[OIII]4959"]
                    smalllineprefix = "%s%.0f_" % ('emission', center_value)

                    refined_model.set_param_hint(name="%samplitude" % smalllineprefix, value=best_fit_params["%samplitude" % biglineprefix])

                else:
                    GaussianModel(prefix=)

        '''





        logger.info("Writing output with lines information to lines%s" % spectrum)
        lines_table.write("lines%s" % spectrum, overwrite=True, format='csv', delimiter='\t')

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


for ax in cont_axs.flat:
    ax.legend(loc='best')
best_fit_ax.legend(loc="best")
initial_fit_ax.legend(loc='best')
#final_lines.legend(loc='best')

plt.show()
