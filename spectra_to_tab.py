# -*- coding: utf-8 -*
# Script  to convert spectra from MUSE to tab separated value
# Author: Andres Gurpide Lasheras  andres.gurpide@gmail.com
# 05-08-2019
# !/usr/bin/env python3


from astropy.io import fits
import argparse
import numpy as np
import os
from math import sqrt


def convertospecviz(wave, data, variance, outname="specviz.fits"):
    """Convert the input data into tab separated values to be read by specviz software"""
    with open(outname, "w") as outfile:
        outfile.write("#Wave(angs)\tData\tVar\tStd\n")
        [outfile.write("%.2f\t%.2f\t%.2f\t%.2f\n" % (wave, flux, var, sqrt(var))) for wave, flux, var in zip(wave, data, variance)]
        print("File %s transformed successfully to %s specviz format" % (basename, outname))


def convertoxspec(wave, data, variance, delwave, outname="xspec.fits"):
    """Convert input data to xspec tab separated format to be used by ftflx2xsp command"""
    with open(outname, "w") as outfile:

        [outfile.write("%.2f\t%.2f\t%.2f\t%.2f\n" % (wave - delwave / 2, wave + delwave / 2, flux, sqrt(var))) for wave, flux, var in zip(wave, data, variance)]
        print("File %s transformed successfully to %s xspec format" % (basename, outname))


# read arguments
ap = argparse.ArgumentParser(description='Converts MUSE spectra to tab separated values')
ap.add_argument("inputspectra", nargs='+', help="Input spectra to be transformed")
args = ap.parse_args()

for spectrum in args.inputspectra:
    hdul = fits.open(spectrum)
    hdudata = hdul["DATA"]
    # convert 0-max to wave and add initial value to get wavelengths
    wave = hdudata.header['CDELT1'] * np.arange(hdudata.header['NAXIS1']) + hdudata.header['CRVAL1']
    variance = hdul["STAT"].data
    data = hdudata.data

    basename = os.path.basename(spectrum)
    outnamespecviz = basename.replace(".fits", "specviztab.dat")
    outnamexspec = basename.replace(".fits", "xspectab.dat")
    outnamestd = basename.replace(".fits", "withstd.fits")
    convertospecviz(wave, data, variance, outnamespecviz)
    convertoxspec(wave, data, variance, hdudata.header['CDELT1'], outnamexspec)

    # add std to the spectrum
    hdul.append(fits.ImageHDU(name="STD", data=np.sqrt(variance), header=hdudata.header.copy()))
    hdul.writeto(outnamestd, overwrite=True)
