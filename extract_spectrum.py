
# -*- coding: utf-8 -*-
# Script to extract a spectrum from a certain region (circle) around a center of the cube. \
# The final spectrum is the mean spectrum of the region
# Author: Andres Gurpide Lasheras  andres.gurpide@gmail.com
# 25-03-2019
# !/usr/bin/env python3
# imports
from mpdaf.obj import Cube
import argparse
import os
import matplotlib.pyplot as plt


# read arguments
ap = argparse.ArgumentParser(description='Muse data to be loaded')
ap.add_argument("input_files", nargs='+', help="List of fits muse data cubes to be loaded")
ap.add_argument("-c", "--center", nargs=1, help="Subcube center in degrees, separated by coma")
ap.add_argument("-r", "--radius", nargs=1, help="Radius of the circle of the subcube in arcsec", default="5", type=float)
ap.add_argument("-s", "--shape", nargs='?', help="Shape: box or circle, circle by default", default="circle", type=str)
ap.add_argument("-o", "--output", nargs='?', help="Name of the output extracted spectrum")

args = ap.parse_args()

muse_cubes = args.input_files
coordinates = args.center[0]
subcube_center = coordinates.split(",")
centerx = float(subcube_center[0])
centery = float(subcube_center[1])
subcube_radius = args.radius[0]
region_shape = args.shape
outputname = args.output
i = 0
plt.figure()

for cubefile in muse_cubes:

    if os.path.isfile(cubefile):
        print('Loading cube %s ...' % cubefile)
        cube = Cube(cubefile)

        print("Extracting spectrum number %i  from cube %s around position (%s) with \
        a radius of %.2f arcsec" % (i, cubefile, coordinates, subcube_radius))
        if region_shape == 'circle':
            subcube = cube.subcube_circle_aperture((centery, centerx), subcube_radius)
        elif region_shape == 'box':
            subcube = cube.subcube((centery, centerx), subcube_radius)
        # free memory immediately
        cube = None
        spe = subcube.mean(axis=(1, 2))
        image = subcube.sum(axis=0)

        spe.info()
        spe.plot()
        plt.show()
        if region_shape == 'box':
            extraction_text = "%s(%s,%.2f\",%.2f\",0)" % (region_shape, coordinates, subcube_radius, subcube_radius)
        elif region_shape == 'circle':
            extraction_text = "%s(%s,%.2f\")" % (region_shape, coordinates, subcube_radius)
        regtext = "# Region file format: DS9 version 4.1 \n global color=green dashlist=8 3 width=1 \
font=\"helvetica 10 normal roman\" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1 \n \
fk5 \n %s # color=cyan width=2 text={SRC} dash=1 \n" % extraction_text
        if outputname != "":
            regfile = open("%s%iextraction.reg" % (outputname, i), "w")
            regfile.write(regtext)
            spe.write("%s_spec%i.fits" % (outputname, i))
            image.write("%s_img%i.fits" % (outputname, i))
            subcube.write("%s_subcube%i.fits" % (outputname, i))
        else:
            regfile = open("%s%iextraction.reg" % (outputname, i), "w")
            regfile.write(regtext)
            spe.write("spectra%i.fits" % i)
            image.write("%s_img%i.fits" % (outputname, i))
            subcube.write("subcube%i.fits" % i)
        regfile.close()
