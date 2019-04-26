# -*- coding: utf-8 -*
# Class for line manipulation
# Author: Andres Gurpide Lasheras  andres.gurpide@gmail.com
# 24-04-2019
# !/usr/bin/env python3


class Line:
    '''Class that holds line parameters for easy manipulation.'''

    def __init__(self, name, lmin, lmax, type, lref, fwhm=None):
        self.name = name
        self.lmin = float(lmin)
        self.lmax = float(lmax)
        self.lref = float(lref)
        self.type = type
        self.fwhm = fwhm

    def __str__(self):
        return self.name
