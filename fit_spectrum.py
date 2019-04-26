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
import sys
from astropy.stats import gaussian_fwhm_to_sigma
from line import Line
from lmfit.models import LorentzianModel, ConstantModel, LinearModel, GaussianModel
sys.path.insert(0, '/home/agurpide/scripts/pythonscripts')
import plot_utils.plot_functions as pf


def compute_FWHM(fwhm, linepeak):
    """Compute FWHM in km/s."""
    fwhm = fwhm / linepeak * c.to('km/s').value
    return fwhm


def fit_continuum(spectrum, min_wa, max_wa,  mask_rest):
    """Fits the continuum of the input spectrum masking the wavelenghts arrays provided.
    The boolean mask_rest can be used to mask the spectrum from the maxixum wavelength given in max_wa
    # provided until the last wavelenght in the spectrumself.
    Parameters
    ----------
    """
    for index in np.arange(len(max_wa)):
        print("Masking region from %s to %s" % (min_wa[index], max_wa[index]))
        # mask wavelength ranges between lines
        spectrum.mask_region(lmin=float(min_wa[index]), lmax=float(max_wa[index]), unit=u.angstrom, inside=True)

    # mask wavelenghts from the highest input wavelenght until the end of the spectrum
    if mask_rest:
        print("Masking everything above %s A in continuum fit"  % max_wa[-1])
        spectrum.mask_region(lmin=float(max_wa[-1]), lmax=None, unit=u.angstrom, inside=True)

    continuum = spectrum.poly_spec(1)
    spectrum.unmask()
    continuum.plot(title="Fitted continuum", color='pink')
    plt.show()
    return continuum


# read arguments
ap = argparse.ArgumentParser(description='Spectrum fits to be loaded')
ap.add_argument("input_files", nargs='+', help="List of spectra to be loaded")

args = ap.parse_args()

spectra = args.input_files

# create color array
colors = pf.create_color_array(len(spectra))

# text label params
textfontsize = 5
textabove = 15

lorentzian = LorentzianModel()

path_lines_file = os.environ['HOME'] + "/scripts/pythonscripts/muse/config_files/fit_spectrum.txt"

with open(path_lines_file) as line_file:
    info_lines = np.loadtxt(line_file, dtype='str', delimiter='\t\t')
    lines = [Line(line[0], line[1], line[2], line[3], line[4], line[5]) for line in info_lines]

for spectrum in spectra:
    plt.figure()
    if os.path.isfile(spectrum):
        print('Loading spectrum %s ...' % spectrum)
        # read data and variance
        spe = Spectrum(filename=spectrum)
        spe.info()
        spe.plot(color='blue')

        # continuum =fit_continuum(spe,min_wa,max_wa,True)

        for line in lines:

            if line.lmin == line.lmax:
                print("Error: lmin (%.2f) has same value as lmax (%.2f) for line %s" % (line.lmin, line.lmax, line))
                sys.exit

            if line.type == 'G':
                print("\nFitting line %s with Gaussian profile" % line)
                guessed_fwhm = line.fwhm
                if guessed_fwhm == 'None':
                    guessed_fwhm = None
                else:
                    guessed_fwhm = float(guessed_fwhm)
                line_gauss = spe.gauss_fit(lmin=line.lmin, lmax=line.lmax, unit=u.angstrom, plot=True, fwhm=guessed_fwhm)

                FWHM = compute_FWHM(line_gauss.fwhm, line_gauss.lpeak)
                print("FWHM (km/s) %.2f \n" % FWHM)
                plt.text(line_gauss.lpeak, line_gauss.peak * 1.4 + line_gauss.cont * 1.3, line, fontsize=textfontsize)
                line_gauss.print_param()
            elif line.type == 'A':
                print("\nFitting line %s with asymmetric Gaussian profile" % line)
                line_agauss = spe.gauss_asymfit(lmin=line.lmin, lmax=line.lmax, unit=u.angstrom, plot=True)

                plt.text(line_agauss[0].lpeak, line_agauss[0].peak, line + line_agauss[0].cont + textabove, fontsize=textfontsize)
                FWHM = compute_FWHM(line_agauss[0].fwhm, line_gauss[0].lpeak)

                print("FWHM (km/s) %.2f \n" % FWHM)
                line_agauss[0].print_param()

            elif line.type == 'L':
                print("\nFitting line %s with Lorentzian profile" % line)
                subspec = spe.subspec(line.lmin, line.lmax, unit=u.angstrom)
                data = subspec._interp_data(False)
                wavelengths = subspec.wave.coord(unit=u.angstrom)
                fmin = data[0]
                fmax = data[-1]
                lmin = wavelengths[0]
                lmax = wavelengths[-1]
                lpeak = wavelengths[np.argmax(data)]
                continuum = ((fmax - fmin) * lpeak + lmax * fmin - lmin * fmax) / (lmax - lmin)

                continuumModel = LinearModel()
                fwhm = subspec.fwhm(lpeak, continuum, False, unit=u.angstrom)
                sigma = fwhm * gaussian_fwhm_to_sigma
                pixel = subspec.wave.pixel(lpeak, nearest=True, unit=u.angstrom)
                peak = data[pixel] - continuum

                if subspec.var is not None:
                    weights = 1.0 / np.sqrt(np.abs(subspec.var))
                else:
                    print("Warning:Variances are null, using 1 as weights for the curve fitting")
                    weights = np.ones(len(wavelengths))
                model = lorentzian + continuumModel
                # model.guess(center=lpeak,sigma=sigma,amplitude=flux,c=continuum)
                output = model.fit(data, x=wavelengths, center=lpeak, sigma=sigma, amplitude=peak, slope=0, intercept=continuum, weights=weights)
                best_fit_params = output.params
                line_center = best_fit_params.get("center").value
                line_peak = best_fit_params.get("height").value
                fwhm = best_fit_params.get("fwhm").value

                FWHM = compute_FWHM(fwhm, line_center)

                plt.text(line_center, line_peak + textabove, line, fontsize=textfontsize)

                print("FWHM (km/s) %.2f \n" % FWHM)

                continuum = best_fit_params.get("intercept").value

                height = best_fit_params.get("height").value - continuum

                amplitude = best_fit_params.get("amplitude").value
                print("FWHM (A) %.2f" % fwhm)
                print("Flux %.2f" % amplitude)
                print("Line center %.2f vs Line reference wavelength %.2f (A)" % (line_center, line.lref))

                plt.plot(wavelengths, output.best_fit, 'g-')
            elif line.type == 'D':
                print("\nFitting line %s with Double Gaussian profile" % line)

                line_dgauss = spe.gauss_dfit(lmin=line.lmin, lmax=line.lmax, wratio=1.54, unit=u.angstrom, plot=True)

                # plt.text(line_dgauss[0].lpeak,line_dgauss[0].peak,line  + line_dgauss[0].cont  + textabove,fontsize=textfontsize)
                FWHM1 = compute_FWHM(line_dgauss[0].fwhm, line_dgauss[0].lpeak)
                FWHM2 = compute_FWHM(line_dgauss[1].fwhm, line_dgauss[1].lpeak)

                print("FWHM (km/s) of line 1 %.2f \n" % FWHM1)
                print("FWHM (km/s) %.2f \n" % FWHM2)
                line_dgauss[0].print_param()
                line_dgauss[1].print_param()

        plt.minorticks_on()
        outputfile = spectrum.replace(".fits", "fit.pdf")
        plt.savefig(outputfile, bbox_inches='tight')

plt.show()
