
# -*- coding: utf-8 -*-
#Script to load, fit and plot spectra generated from MUSE cubes
# Author: Andres Gurpide Lasheras  andres.gurpide@gmail.com
#25-03-2019
    #!/usr/bin/env python3
    # coding=utf-8
#imports
from mpdaf.obj import Cube
from mpdaf.drs import PixTable
from mpdaf.obj import deg2sexa
import matplotlib.pyplot as plt
from mpdaf.obj import Spectrum, WaveCoord
import argparse
import numpy as np
import os
import astropy.units as u

#Compute FWHM in km/s
def compute_FWHM(line,c):
    fwhm =  line.fwhm/line.lpeak * c
    return fwhm


#read arguments
ap = argparse.ArgumentParser(description='Spectrum fits to be loaded')
ap.add_argument("input_files",nargs='+',help="List of spectra to be loaded")
#0:55:04.7485
#-37:41:43.420
#2.5 arcsec
args = ap.parse_args()

spectra=args.input_files

#speed of light in km/s
c=299792.485

#create color array
x = np.arange(len(spectra))
ys = [i+x+(i*x)**2 for i in range(len(spectra))]
setmap = plt.get_cmap('jet')
colors = setmap(np.linspace(0, 1, len(ys)))

plt.figure()
i=0

config_file=~/sc

with open(cstatlinescan_file) as f:
    print("Reading file: %s" %(cstatlinescan_file))
    table = np.loadtxt(cstatlinescan_file)

for spectrum in spectra:

    if os.path.isfile(spectrum):
        print('Loading spectrum %s ...' %spectrum)
        #read data and variance
        spe = Spectrum(filename=spectrum)
        spe.info()
        spe.plot(color=colors[i])
        i=i+1
        hbeta = spe.gauss_fit(lmin=4845, lmax=4886, unit=u.angstrom, plot=True)
        print("\n H$_\Beta$ line: ")
        hbetaFWHM = compute_FWHM(hbeta,c)
        print("FWHM (km/s) %.2f \n" %hbetaFWHM)

        hbeta.print_param()

        OIII4958 = spe.gauss_fit(lmin=4950, lmax=4970, unit=u.angstrom, plot=True)
        print("\n [OIII]4958$\lambda$ line: ")
        OIII4958FWHM = compute_FWHM(OIII4958,c)
        print("FWHM (km/s) %.2f \n" %OIII4958FWHM)
        OIII4958.print_param()

        OIII5006 = spe.gauss_fit(lmin=4995, lmax=5020, unit=u.angstrom, plot=True)
        print("\n [OIII]5006$\lambda$ line: ")
        OIII5006FWHM = compute_FWHM(OIII5006,c)
        print("FWHM (km/s) %.2f \n" %OIII5006FWHM)
        OIII5006.print_param()

        NII  = spe.gauss_asymfit(lmin=5576.83, lmax=5590, unit=u.angstrom, plot=True)
        print("\n [NII]5754$\lambda$ line: ")
        NIIFWHM = compute_FWHM(NII[0],c)
        print("FWHM (km/s) %.2f \n" %NIIFWHM)
        NII[0].print_param()

        FeVII5720 =spe.gauss_fit(lmin=5705, lmax=5740, unit=u.angstrom, plot=True)
        print("\n [FeVII]5720$\lambda$ line: ")
        FeVII5875FWHM = compute_FWHM(FeVII5720,c)
        print("FWHM (km/s) %.2f \n" %FeVII5875FWHM)
        FeVII5720.print_param()

        NeII5754 =spe.gauss_fit(lmin=5750, lmax=5765, unit=u.angstrom, plot=True)
        print("\n [FeVII]5720$\lambda$ line: ")
        FNeII5754FWHM = compute_FWHM(NeII5754,c)
        print("FWHM (km/s) %.2f \n" %FNeII5754FWHM)
        NeII5754.print_param()

        HeI5875 =spe.gauss_fit(lmin=5860, lmax=5900, unit=u.angstrom, plot=True)
        print("\n [HeI]5875$\lambda$ line: ")
        HeI5875FWHM = compute_FWHM(HeI5875,c)
        print("FWHM (km/s) %.2f \n" %HeI5875FWHM)
        HeI5875.print_param()

        FeVII =spe.gauss_fit(lmin=6080, lmax=6100, unit=u.angstrom, plot=True)
        print("\n [FeVIII]6087$\lambda$ line: ")
        FeVIIFWHM = compute_FWHM(FeVII,c)
        print("FWHM (km/s) %.2f \n" %FeVIIFWHM)
        FeVII.print_param()

        OI63000  = spe.gauss_asymfit(lmin=6300, lmax=6320, unit=u.angstrom, plot=True)
        print("\n [OI]6300$\lambda$ line: ")
        OI63000FWHM = compute_FWHM(OI63000[0],c)
        print("FWHM (km/s) %.2f \n" %OI63000FWHM)
        OI63000[0].print_param()

        halpha = spe.gauss_fit(lmin=6528, lmax=6608, unit=u.angstrom, plot=True)
        print("\n H$_\alpha$ line: ")
        halphaFWHM = compute_FWHM(halpha,c)
        print("FWHM (km/s) %.2f \n" %halphaFWHM)
        halpha.print_param()

        NII6583 = spe.gauss_fit(lmin=6675, lmax=6685, unit=u.angstrom, plot=True)
        print("\n [NII]6583$\lambda$ line: ")
        NII6583FWHM = compute_FWHM(NII6583,c)
        print("FWHM (km/s) %.2f \n" %NII6583FWHM)
        NII6583.print_param()

        HeI7065 = spe.gauss_fit(lmin=7060, lmax=7080, unit=u.angstrom, plot=True)
        print("\n [HeI]7065$\lambda$ line: ")
        HeI7065FWHM = compute_FWHM(HeI7065,c)
        print("FWHM (km/s) %.2f \n" %HeI7065FWHM)
        HeI7065.print_param()

        FeII = spe.gauss_fit(lmin=7145, lmax=7168, unit=u.angstrom, plot=True)
        print("\n [FeII]7155$\lambda$ line: ")
        FeIIFWHM = compute_FWHM(FeII,c)
        print("FWHM (km/s) %.2f \n" %FeIIFWHM)
        FeII.print_param()

        OI8446 = spe.gauss_fit(lmin=8440, lmax=8460, unit=u.angstrom, plot=True)
        print("\n [OI]8446$\lambda$ line: ")
        OI8446FWHM = compute_FWHM(OI8446,c)
        print("FWHM (km/s) %.2f \n" %OI8446FWHM)
        OI8446.print_param()










plt.show()
#0:55:04.7485
#-37:41:43.420
#2.5 arcsec
