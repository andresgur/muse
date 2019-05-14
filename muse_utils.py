# -*- coding: utf-8 -*-
# Methods to read config line files
# Author: Andres Gurpide Lasheras  andres.gurpide@gmail.com
# 04-04-2019
    #!/usr/bin/env python3
    # coding=utf-8
# imports
import matplotlib.pyplot as plt
import numpy as np
import sys
import re
from regions import read_ds9
import logging
from mpdaf.obj import Image
sys.path.insert(0, '/home/agurpide/scripts/pythonscripts')
import plot_utils.plot_functions as pf


def load_map_file(file):
    """Load config file for line map with line names and min and max wavelenghts to be masked.
    Parameters
    ----------
    file : the config file containing the lines to be processed."""
    with open(file) as line_file:
        info_lines = np.loadtxt(line_file, dtype='str', delimiter='\t\t')
        line_names = info_lines[:, 0]
        min_wa = np.array(info_lines[:, 1])
        max_wa = np.array(info_lines[:, 2])
        print("Lines found in %s \n" % file)
        print(line_names)
        return line_names, min_wa, max_wa


def load_spectrum_file(spectrum_file):
    """Load config file for spectrum fitting.
    Parameters
    ----------
    spectrum_file : the config file containing the lines to be fitted."""
    with open(spectrum_file) as config_file:
        info_lines = np.loadtxt(config_file, dtype='str', delimiter='\t\t')
        line_names = info_lines[:, 0]
        min_wa = np.array(info_lines[:, 1])
        max_wa = np.array(info_lines[:, 2])
        line_type = info_lines[:, 3]
        fwhms = info_lines[:, 5]
        return line_names, min_wa, max_wa, line_type, fwhms


def add_zlabel(image_file, units):
    '''Decide the z label for the input image based on the file name or the units.'''

    if units != "None":
        zlabel = units

    if 'flux' in image_file:
        zlabel = "1e-20 erg/cm$^2$/s"
    elif 'wave' in image_file:
        if 'disp' in image_file:
            zlabel = '$\Delta\lambda (\Delta\AA)$'
        else:
            zlabel = '$\lambda (\AA)$'

    elif 'disp' in image_file:
        zlabel = '$\Delta$v (km/s)'

    elif 'vel' in image_file:
        zlabel = 'v (km/s)'

    return zlabel


def plot_image(image_file, outputname, region=None, cutting_region=None, vmin=None):
    """Generate colorbar image from input string file.

    Parameters:
    image_file : The path to the fits image file to be plotted
    outputname : The output name to be given to the produced image
    region : a region file indicating the regions to be overlaid on the image
    cutting_region : a region file with a circular region to make a cut on the cube
    vmin : the minimum value for the colorbar scale

    """
    image = Image(image_file)

    if cutting_region is not None:
        regs = read_ds9(cutting_region)
        if len(regs) > 1:
            logging.warning("Only one region is allowed for the spatial cut")
        else:
            cut_region = regs[0]
            ra = cut_region.center.ra.value
            dec = cut_region.center.dec.value
            radius = cut_region.radius.value

            logging.debug("Cutting the image at (ra, dec) = (%.3f, %.3f) around %.2f asec" % (ra, dec, radius))
            image = image.subimage(center=(dec, ra), size=radius, minsize=0.1)

    img_figure = plt.figure()
    ax = plt.subplot(projection=image.wcs.wcs)
    image.plot(ax=ax, scale='linear', colorbar='v', show_xlabel=True, show_ylabel=True, zscale=False, vmin=vmin)
    ax.set_xlabel('Ra')
    ax.set_ylabel('Dec')
    if region is not None:
        regs = read_ds9(region)

        for reg in regs:
            print(reg)
            reg.plot(ax=ax)

    im = ax.images
    cb = im[-1].colorbar

    cb.ax.set_ylabel(add_zlabel(image_file, image.unit), fontsize=14)

    pf.set_plot_format(14)
    pf.save_plot(img_figure, outputname)
    plt.close(img_figure)


def read_regionfromfile(region_file):
    '''Read the region parameters from file. Units must be pixel detector.
    Parameters
    ----------
    region_file : The region file containing the desired regions to be plot.
    Returns all the parameters for every region.'''
    regions = []
    with open(region_file) as file:
        for line in file:
            if 'ellipse' in line.lower():
                region = 'ellipse'
                print("Region %s found " % region)
                cleaned = re.sub("\(|\).*\n", "", line)
                trimmed = cleaned.lower().replace(region, "")
                params = trimmed.split(',')
                print("%s parameters: " % region)
                print(params)
                regions.append([float(i) for i in params])
            elif 'circle' in line.lower():
                region = 'circle'
                print("Region %s found " % region)
                cleaned = re.sub("\(|\)\\n", "", line)
                print(cleaned)
                trimmed = cleaned.replace(region, "")
                params = trimmed.split(',')
                print("%s parameters: " % region)
                print(params)
                regions.append([float(i) for i in params])
    print("Found %i region(s)" % len(regions))

    return np.array(regions)
