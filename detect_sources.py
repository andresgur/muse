from mpdaf.sdetect import muselet
import argparse
import numpy as np
import os
import matplotlib.pyplot as plt
#read arguments
ap = argparse.ArgumentParser(description='Muse data to be loaded')
ap.add_argument("input_files",nargs='+',help="List of fits muse data cubes to be loaded")

args = ap.parse_args()

muse_cubes=args.input_files

for cubefile in muse_cubes:
    muselet(cubefile)
