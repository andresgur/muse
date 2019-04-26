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


def plot_image(image, outputname, zlabel=None):
    image = Image(image)
    img_figure = plt.figure()
    ax = plt.subplot(projection=image.wcs.wcs)
    image.plot(ax=ax, scale='linear', colorbar='v', show_xlabel=True, show_ylabel=True)
    ax.set_xlabel('Ra')
    ax.set_ylabel('Dec')
    im = ax.images
    cb = im[-1].colorbar
    if zlabel is not None:
        cb.ax.set_ylabel(zlabel, fontsize=20)
    pf.set_plot_format(14)
    pf.save_plot(img_figure, outputname)
    plt.close(img_figure)
