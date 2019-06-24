# -*- coding: utf-8 -*
# Script to remove nan values from the spaxel pixels by interpolating from the neighbours. Uses ZAP method
# Author: Andres Gurpide Lasheras  andres.gurpide@gmail.com
# 24-05-2019
# !/usr/bin/env python3
import zap
import argparse
import os

ap = argparse.ArgumentParser(description='Removes NaN values from the input cube using ZAP tool')
ap.add_argument("input_cubes", nargs='+', help="Muse data cube to be cleaned from sky residuals")
ap.add_argument("-o", "--outdir", help="Output dir", nargs='?', type=str, default='nancleancubes')
ap.add_argument("-b", "--box", help="Defines the number of pixels around the offending NaN pixel.", default=1, type=int)
# parse args
args = ap.parse_args()

if not os.path.isdir(args.outdir):
    os.mkdir(args.outdir)

for cube in args.input_cubes:
    zap.nancleanfits(cube, "%s/nanclean%s" % (args.outdir, cube), overwrite=True, boxsz=args.box)
