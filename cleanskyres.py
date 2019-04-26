# -*- coding: utf-8 -*
# Script to clean a cube from sky residuqls
# Author: Andres Gurpide Lasheras  andres.gurpide@gmail.com
# 12-04-2019
# !/usr/bin/env python3

# imports
import numpy as np
import matplotlib.pyplot as plt
from mpdaf.obj import Image
import zap
import argparse


def create_mask_region(image, percentile):
    """Create mask region from white light image and percentile flux threshold.

    Creates mask region out of a white ligh image.

    Parameters
    ----------
    image : white light image with same size as the cube to be cleaned
    percentile : threshold percentile of the histogram flux to set mask regions
    """
    print("Creating masked region")

    threshold_flux = np.percentile(image.data, percentile)
    print("Masking every pixel above %.2f %s" % (threshold_flux, image.unit.to_string()))
    image.info()
    mask_image = image.clone(data_init=np.zeros)

    mask_image[np.where(image.data > threshold_flux)] = 1

    return mask_image


# read arguments
ap = argparse.ArgumentParser(description='Clean input cubes from sky residuals using the ZAP tool')
ap.add_argument("input_cube", nargs=1, help="Muse data cube to be cleaned from sky residuals")
ap.add_argument("-i", "--image", nargs='?', help="White light image from the cube to compute the zero order sky residuals spectrum")
ap.add_argument("-p", "--percentatge", nargs='?', help="Percentatge threshold  of the white light image flux histogram to be masked \
 for the zero order sky calculation ", default="50", type=float)
ap.add_argument("-o", "--output", nargs='?', help="Name of the output cleaned cube", default="")

# parse args
args = ap.parse_args()

input_cube = args.input_cube[0]
image_file = args.image
percentile = args.percentatge
output_cube = args.output

white_light_image = Image(image_file)


if image_file is not None:
    mask_image = create_mask_region(white_light_image, percentile)
    print("Saving masked image")
    mask_file = image_file.replace(".fits", "masked.fits")
    mask_image.write(mask_file)
else:
    mask_file = None
    print("Warning: no image was provided to compute the masked regions")

if output_cube == "":
    output_cube = input_cube.replace(".fits", "cleaned.fits")

print("Name for the cleaned cube: %s" % output_cube)

# cfwidthSVD default is 300; SP filter for the calculation of the eigenvalues ~ 50
zap.process(input_cube, skycubefits='skycube.fits', mask=mask_file,
            outcubefits=output_cube, cfwidthSVD=300, cfwidthSP=50, nevals=15,
            varcurvefits='varcurve.fits', overwrite=True, clean=True)

'''
process.writecube(outcubefits='cube.fits', overwrite=True)
process.writeskycube(skycubefits='skycube.fits', overwrite=True)


process.plotvarcurve()
plt.figure()
plt.plot(process.cube[:, :, :].sum(axis=(1, 2)), 'b', alpha=0.5)
plt.plot(process.cleancube[:, :, :].sum(axis=(1, 2)), 'g')

plt.show()
'''
