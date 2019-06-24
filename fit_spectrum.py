# -*- coding: utf-8 -*
# Script to load, fit and plot spectra generated from MUSE cubes
# Author: Andres Gurpide Lasheras  andres.gurpide@gmail.com
# 25-03-2019
# !/usr/bin/env python3

# imports
from mpdaf.obj import Spectrum
import matplotlib.pyplot as plt
import argparse
import numpy as np
import os
import astropy.units as u
from astropy.constants import c
import math as m
from astropy.stats import gaussian_fwhm_to_sigma
from line import Line
import logging
from lmfit.models import LorentzianModel, PolynomialModel, GaussianModel
import sys
sys.path.insert(0, '/home/agurpide/scripts/pythonscripts')
import plot_utils.plot_functions as pf


def line_type_toint(line_type):
    """Check whether the input line type is an integer i.e. multiple gaussian have to be fitted."""
    try:
        int(line_type)
        return True
    except ValueError:
        return False


def guess_init_line(data, wavelengths, initial_fwhm, linetype):
    """Guess initial line parameters and continuum from the subspectrum.

    Computes an initial guess for the continuum level, the width, position and peak of the line.
    Parameters
    ----------
    data : the flux values of the spectrum
    wavelenghts : the wavelength values
    initial_fwhm : whether the full width half maximum should be guessed or not

    """
    fmin = data[0]
    fmax = data[-1]
    lmin = wavelengths[0]
    lmax = wavelengths[-1]
    if '-' in linetype:
        lpeak = wavelengths[np.argmin(data)]
    else:
        lpeak = wavelengths[np.argmax(data)]

    continuum = ((fmax - fmin) * lpeak + lmax * fmin - lmin * fmax) / (lmax - lmin)

    if line.init_fwhm == 'None' and '-' not in linetype:
        # estimate fwhm
        try:
            initial_fwhm = subspec.fwhm(lpeak, continuum, False, unit=u.angstrom)
        except ValueError:
            logger.warning("Could not compute fwhm. Setting it to a constant initial value", exc_info=True)
            initial_fwhm = def_fwhm
    elif line.init_fwhm != 'None':
        initial_fwhm = float(line.init_fwhm)
    elif '-' in linetype:
        initial_fwhm = def_fwhm

    sigma = initial_fwhm * gaussian_fwhm_to_sigma

    pixel = subspec.wave.pixel(lpeak, nearest=True, unit=u.angstrom)
    peak = np.abs(data[pixel] - continuum)
    return continuum, sigma, lpeak, peak


def fit_continuum(spectrum, min_wa, max_wa, mask_rest):
    """Fits the continuum of the input spectrum masking the wavelenghts arrays provided.

    Parameters
    ----------
    spectrum : the spectrum whose continuum is to be fitted.
    min_wa : the minimum wavelength to consider in the fitting process (this can be an array too)
    max_wa : the maximum wavelength to consider in the fitting process (this can be an array too)
    mask_rest : boolean to indicate whether everything above the last maximum wavelenght should be mask (i.e. ignored)

    """
    for index in np.arange(len(max_wa)):
        print("Masking region from %s to %s" % (min_wa[index], max_wa[index]))
        # mask wavelength ranges between lines
        spectrum.mask_region(lmin=float(min_wa[index]), lmax=float(max_wa[index]), unit=u.angstrom, inside=True)

    # mask wavelenghts from the highest input wavelenght until the end of the spectrum
    if mask_rest:
        print("Masking everything above %s A in continuum fit" % max_wa[-1])
        spectrum.mask_region(lmin=float(max_wa[-1]), lmax=None, unit=u.angstrom, inside=True)

    continuum = spectrum.poly_spec(1)
    spectrum.unmask()
    continuum.plot(title="Fitted continuum", color='pink')
    plt.show()
    return continuum


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
ap.add_argument("-f", "--config_file", nargs='?', help="Config file with the lines to be fitted", default=os.environ['HOME'] + "/scripts/pythonscripts/muse/config_files/fit_spectrum.txt")
ap.add_argument("-z", "--red_shift", nargs='?', help="Redshift of the object", default=0, type=float)
ap.add_argument("-s", "--significance", nargs='?', help="Threshold significance to consider a line as present in EQW units", default=0.4, type=float)

args = ap.parse_args()

spectra = args.input_files
config_file = args.config_file
red_shift = args.red_shift
line_threshold = args.significance
# create color array
colors = pf.create_color_array(len(spectra))
# default_fwhm in A
def_fwhm = 1.8

# fwhm
max_fwhm = 20
min_fwhm = 1.0

# flags
plot_initfit = 0
plot_reference = 1
plot_residuals = 0
plot_result = 1
plot_lines = 0
plot_chisq = 0
plot_comp = 0

# text label params
textfontsize = 7
textabove = 15
gauss_color = ['magenta']
lorentzian_color = ['green']
multiple_line_color = ['cyan', 'orange', 'brown', 'black']
# linear model prefix
continuum_prefix = 'continuum_'
line_prefix = 'line_'
cont_color = 'yellow'

with open(config_file) as line_file:

    info_lines = np.genfromtxt(line_file, names=True, dtype=('U10', float, float, 'U10', float, 'U10', float, int), delimiter='\t\t')
    logger.debug("Fitting %i line(s)" % info_lines.size)
    lines = [Line(line[0], line[1], line[2], line[3], line[4], line[5], line[6], line[7]) for line in info_lines]

spectrum_figure = plt.figure(1, figsize=(16.0, 10.0))
spectrum_ax = spectrum_figure.add_subplot(1, 1, 1)

for spectrum, color in zip(spectra, colors):
    spectrum_ax.minorticks_on()
    if os.path.isfile(spectrum):
        logger.info('Loading spectrum %s ...' % spectrum)
        # read data and variance
        spe = Spectrum(filename=spectrum)
        spe.info()
        #spectrum_ax.errorbar(spe.wave.coord(unit=u.angstrom), spe.data, yerr=spe.var, color=color, ls="-",
        #            linewidth=0.5, elinewidth=1, markersize=4, errorevery=2, drawstyle='steps-mid')
        spe.plot(color=color, ax=spectrum_ax)

        for line in lines:
            subspec = spe.subspec(line.lmin, line.lmax, unit=u.angstrom)
            data = subspec._interp_data(True)
            wavelengths = subspec.wave.coord(unit=u.angstrom)
            continuum, sigma, lpeak, peak = guess_init_line(data, wavelengths, line.init_fwhm, line.type)

            if subspec.var is not None:
                weights = 1.0 / np.sqrt(subspec.var)
                np.ma.fix_invalid(weights, copy=False, fill_value=1)
            else:
                logger.warning("Variances are null, using 1 as weights for the curve fitting")
                weights = np.ones(len(wavelengths))

            if line.lmin >= line.lmax:
                logger.error("lmin (%.2f) is equal or higher than lmax (%.2f) for line %s" % (line.lmin, line.lmax, line))
                continue

            # fit multiple gaussians case
            if line_type_toint(line.type):

                # name of first line model
                line_prefixes = ["%s%i" % (line_prefix, 0)]

                linecolor = multiple_line_color

                linemodel = GaussianModel(prefix=line_prefixes[0])

                for gauss_index in np.arange(1, int(line.type)):
                    current_prefix = "%s%i" % (line_prefix, gauss_index)
                    current_line = GaussianModel(prefix=current_prefix)
                    linemodel += current_line
                    line_prefixes.append(current_prefix)

            # one gaussian case
            elif 'G' in line.type:

                linecolor = gauss_color
                line_prefixes = ["%s%i" % (line_prefix, 0)]
                linemodel = GaussianModel(prefix=line_prefixes[0])

            # lorentzian line case
            elif 'L' in line.type:

                linecolor = lorentzian_color
                line_prefixes = ["%s%i" % (line_prefix, 0)]
                linemodel = LorentzianModel(prefix=line_prefixes[0])

            logger.info("\n Fitting line %s with %s profile and with polynomial of degree %i for the continuum" % (line, linemodel.name, line.deg_cont))

            if line.deg_cont > 7:
                logger.warning("Maximum degree for the polynomial continuum fitting is 7. We will use this value instead of the provided: %i" % line.deg_cont)
                line.deg_cont = 7

            continuumModel = PolynomialModel(degree=int(line.deg_cont), prefix=continuum_prefix)

            if '-' in line.type:
                logger.info("Modelling absorption line")
                model = continuumModel - linemodel
            else:
                logger.info("Modelling emission line")
                model = continuumModel + linemodel

            for paramname in model.param_names:
                if 'center' in paramname:
                    max = line.lmax
                    min = line.lmin
                    value = lpeak
                    lpeak += line.linesep

                elif 'sigma' in paramname:
                    min = min_fwhm * gaussian_fwhm_to_sigma
                    max = max_fwhm * gaussian_fwhm_to_sigma

                    value = sigma

                elif 'amplitude' in paramname:
                    min = 0
                    max = 8 * peak * m.sqrt(2 * m.pi) + continuum
                    value = peak * m.sqrt(2 * m.pi) + continuum

                elif 'c0' in paramname:
                    min = 0
                    max = 2 * continuum
                    value = continuum
                else:
                    min = -2
                    max = 2
                    value = 0

                print("Value Max Min %s \n %.2f %.2f %.2f" % (paramname, value, max, min))
                model.set_param_hint(paramname, value=value, vary=True, min=min, max=max)

            # expected line center
            line.ex_center = (1 + red_shift) * line.lref
            logger.info("Expected line position: %.1f" % line.ex_center)

            # fit model
            output = model.fit(data, x=wavelengths, weights=weights)
            logger.debug(output.fit_report())
            best_fit_params = output.params
            line.chisq = output.chisqr
            line.dof = output.nfree
            mod_components = output.eval_components(x=wavelengths)
            # continuum at line peak
            cont_peak = mod_components[continuum_prefix][np.argmax(data)]

            total_line_flux = sum(best_fit_params.get("%samplitude" % prefix).value for prefix in line_prefixes)

            line_eq = total_line_flux / cont_peak

            if line_eq <= line_threshold or not output.errorbars:
                logger.warning("Line %s significance %i is below the threshold")
                line.present = False
            else:
                line.present = True

            # iterate each line model component
            for prefix, index in zip(line_prefixes, np.arange(0, line.line_components)):

                line_center = best_fit_params.get("%scenter" % prefix)
                line_peak = best_fit_params.get("%sheight" % prefix)
                fwhm = best_fit_params.get("%sfwhm" % prefix)
                line_flux = best_fit_params.get("%samplitude" % prefix)

                plt.text(line_center, line_peak.value + textabove, line, fontsize=textfontsize)

                logger.debug("Flux %.2f" % line_flux.value)
                # FWHM
                line_shift = line_center.value - line.ex_center
                line.center[index] = line_center.value
                line.flux[index] = line_flux.value
                line.fwhm[index] = fwhm.value / line_center * c.to('km/s').value
                line.shift[index] = c.to('km/s').value * (line_center.value - line.ex_center) / line.ex_center
                line.int[index] = line_peak.value
                # equivalent width line flux / cont flux
                line.eq[index] = line.flux[index] / cont_peak
                logger.debug("FWHM %.2f (%.2f) A (km/s)  \n" % (fwhm.value, line.fwhm[index]))
                logger.debug("Line shift %.2f (%.2f) A (km/s)" % (line_shift, line.shift[index]))

                if output.errorbars:
                    fwhmkm_std = c.to('km/s').value * (fwhm.stderr / line_center.value + fwhm.value * line_center.stderr / line_center.value**2)
                    # store errors if they were calculated
                    line.center_std[index] = line_center.stderr
                    line.flux_std[index] = line_flux.stderr
                    line.fwhm_std[index] = fwhmkm_std
                    line.shift_std[index] = c.to('km/s').value * (line_center.stderr) / line.ex_center
                    line.int_std[index] = line_peak.stderr

            if plot_residuals:
                # output.plot(spectrum_ax)
                output.plot_residuals()

            if plot_reference:
                spectrum_ax.axvline(x=line.ex_center, ls='--', color='gray', alpha=0.5)

            if line.present:
                spectrum_ax.plot(wavelengths, output.best_fit, 'red')

                if plot_comp:
                    spectrum_ax.plot(wavelengths, mod_components[continuum_prefix], cont_color)
                    for prefix, color in zip(line_prefixes, linecolor):
                        spectrum_ax.plot(wavelengths, mod_components[continuum_prefix] + mod_components[prefix], color)

        #    spectrum_ax.fill_between(wavelengths, output.best_fit - uncertainties, output.best_fit + uncertainties, color='yellow')
            if plot_initfit:
                spectrum_ax.plot(wavelengths, output.init_fit, linecolor, ls='--')

    # write fit spectrum figure
    outputfile_fig = spectrum.replace(".fits", "fit.pdf")
    spectrum_figure.savefig(outputfile_fig, bbox_inches='tight')

    # write fit lines output
    outputfile_dat = spectrum.replace(".fits", "fit.dat")
    out_string = ""
    string_array = ["%s\t\t%.3f\t\t%.1f\t\t%.1f\t\t%.2f\t\t%.0f\t\t%.0f\t\t%.0f\t\t%.0f\t\t%.0f\t\t%.0f\t\t%.1f\t\t%.1f\t\t%.1f\t\t%.0f\t\t%i\t\t%.1f\t\t%r\n" % (line.name,
                    line.lref, line.ex_center, line.center[index], line.center_std[index], line.shift[index], line.shift_std[index], line.fwhm[index],
                    line.fwhm_std[index], line.flux[index], line.flux_std[index], line.int[index], line.int_std[index], line.eq[index], line.chisq, line.dof,
                    line.chisq / line.dof, line.present) for line in lines for index in np.arange(0, line.line_components)]
    out_string = "".join(string_array)

    with open(outputfile_dat, "w") as out_file:
        out_file.write("#name\t\trest_wa\t\tex_center\t\tfit_center\t\tstd\t\tshift\t\tstd\t\tfwhm\t\tstd\t\tflux\t\tstd\t\tint\t\tstd\t\tEW\t\tchi\t\tdof\t\tchisqr\n")
        out_file.write(out_string)

if plot_lines:
    line_figure = plt.figure(2, figsize=(16.0, 10.0))

    plt.minorticks_on()
    plt.xlabel('FWHM (km/s)', fontsize=20)
    plt.ylabel('Line shift (km/s)', fontsize=20)
    logger.debug("Line FWHM(km/s) Flux Shift(km/s)")
    #[print("%s: %.2f %.2f %.2f" % (line.name, line.fwhm, line.flux, line.shift)) for line in lines]

    for line in lines:
        if line.type == 'L':
            marker = 'o'
        else:
            marker = '^'
        if 'Fe' in line.name:
            color = 'orange'
        elif 'He' in line.name:
            color = 'yellow'
        elif 'H' in line.name or 'Pa' in line.name:
            color = 'green'
        elif 'O' in line.name:
            color = 'red'
        elif 'N' in line.name:
            color = 'cyan'
        elif 'Cl' in line.name:
            color = 'blue'
        else:
            color = 'magenta'
        for index in np.arange(0, line.line_components):
            plt.errorbar(line.fwhm[index], line.shift[index], xerr=line.fwhm_std[index], yerr=line.shift_std[index], label=line.name,
                         fmt=marker, markersize=2, color=color)
    plt.legend(fontsize=14, loc='best', ncol=3)
    line_out = spectrum.replace(".fits", "lines.pdf")
    line_figure.savefig(line_out, bbox_inches='tight')

if plot_chisq:
    chisq_figure = plt.figure(3, figsize=(16.0, 10.0))

    plt.minorticks_on()
    plt.xlabel('$\chi$', fontsize=20)
    plt.ylabel('dof', fontsize=20)

    for line in lines:
        if line.type == 'L':
            marker = 'o'
        else:
            marker = '^'
        if 'Fe' in line.name:
            color = 'orange'
        elif 'He' in line.name:
            color = 'yellow'
        elif 'H' in line.name or 'Pa' in line.name:
            color = 'green'
        elif 'O' in line.name:
            color = 'red'
        elif 'N' in line.name:
            color = 'cyan'
        elif 'Cl' in line.name:
            color = 'blue'
        else:
            color = 'magenta'

        plt.scatter(line.chisq, line.dof, label=line.name, color=color)
        plt.scatter(line.dof, line.dof, color='gray')
    plt.legend(fontsize=14, loc='best', ncol=3)
    chisq_out = spectrum.replace(".fits", "chisq.pdf")
    chisq_figure.savefig(chisq_out, bbox_inches='tight')

if plot_result:
    plt.show()
