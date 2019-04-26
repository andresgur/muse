# -*- coding: utf-8 -*
#Script to load, fit and plot spectra generated from MUSE cubes
# Author: Andres Gurpide Lasheras  andres.gurpide@gmail.com
#25-03-2019
    #!/usr/bin/env python3

#imports

import matplotlib.pyplot as plt
import argparse
import numpy as np
import os
import astropy.units as u
from math import *
import sys
sys.path.insert(0, '/home/agurpide/scripts/pythonscripts')
import plot_utils.plot_functions as pf
import muse_utils as muse_utils



#read arguments
ap = argparse.ArgumentParser(description='Spectrum fits to be loaded')
ap.add_argument("input_images",nargs='+',help="List of images to be saved")
args = ap.parse_args()

images=args.input_cubes

for image in images:
    if os.path.isfile(image):
        outputfilename = image.replace("fits","pdf")
        muse_utils.plot_image(image,outputfilename)
