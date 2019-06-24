# -*- coding: utf-8 -*
# Class for line manipulation
# Author: Andres Gurpide Lasheras  andres.gurpide@gmail.com
# 24-04-2019
# !/usr/bin/env python3
import numpy as np

class Line:
    """Class that holds line parameters for easy manipulation."""

    def __init__(self, name, lmin, lmax, type, lref, init_fwhm=None, line_sep=0, deg_cont=1, fwhm=None, fwhm_std=None,
                 flux=None, flux_std=None, shift=None, shift_std=None, chisq=None, dof=None, present=None):
        """Shift and FWHM are in km/s."""
        self.name = name
        self.lmin = lmin
        self.lmax = lmax
        self.lref = lref
        self.type = type
        self.init_fwhm = init_fwhm
        self.linesep = line_sep
        self.deg_cont = deg_cont

        try:
            self.line_components = int(self.type)
        except Exception:
            self.line_components = 1
        self.ex_center = lref
        self.eq = np.empty(self.line_components)
        self.eq_std = np.empty(self.line_components)
        self.center = np.empty(self.line_components)
        self.center_std = np.empty(self.line_components)
        self.shift = np.empty(self.line_components)
        self.shift_std = np.empty(self.line_components)
        self.fwhm = np.empty(self.line_components)
        self.fwhm_std = np.empty(self.line_components)
        self.flux = np.empty(self.line_components)
        self.flux_std = np.empty(self.line_components)
        self.int = np.empty(self.line_components)
        self.int_std = np.empty(self.line_components)
        self.chisq = chisq
        self.dof = dof
        self.present = None
    def __str__(self):
        return self.name
