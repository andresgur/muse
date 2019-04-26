# -*- coding: utf-8 -*
#Script to load, fit and plot spectra generated from MUSE cubes
# Author: Andres Gurpide Lasheras  andres.gurpide@gmail.com
#25-03-2019
    #!/usr/bin/env python3

#imports

import matplotlib.pyplot as plt
from mpdaf.obj import Spectrum, WaveCoord,iter_spe,Cube
import argparse
import numpy as np
import os
import astropy.units as u
from math import *
import sys
sys.path.insert(0, '/home/agurpide/scripts/pythonscripts')
import plot_utils.plot_functions as pf
import muse_utils as muse_utils
from mpdaf.obj import deg2sexa


#read arguments
ap = argparse.ArgumentParser(description='Spectrum fits to be loaded')
ap.add_argument("input_cubes",nargs='+',help="List of fits muse data cubes to be processed for the line map")
ap.add_argument("-l","--linefile",nargs='?',help="File with line limits to be masked for each linemap", default=os.environ['HOME'] +"/scripts/pythonscripts/muse_data/config_files/line_maps.txt",type=str)
ap.add_argument("-o","--outdir",nargs='?',help="Name of the output directory",default="linemaps")
args = ap.parse_args()

muse_cubes=args.input_cubes

line_file_path= args.linefile

outdir= args.outdir

line_names,min_wa,max_wa = muse_utils.load_map_file(line_file_path)

for cubefile in muse_cubes:

    if os.path.isfile(cubefile):
        print('Loading cube %s ...' %cubefile)
        original_cube = Cube(cubefile)

        continuum_cube = original_cube.clone(data_init=np.empty, var_init=np.zeros)
        #create image for line fluxes maps
        original_image = original_cube.sum(axis=0)
        #image_masked = original_image.clone(data_init= np.empty)

        for current_spec,(y,x) in iter_spe(original_cube,index=True):
            coord_sky = original_cube.wcs.pix2sky([y,x], unit=u.deg)
            dec, ra = deg2sexa(coord_sky)[0]
            print("Processing spectrum at position (x,y ) = (%i,%i)... (dec,ra) = (%s,%s)" %(x,y,dec,ra))
            #get continuum masking important wavelengths
            if not current_spec.data.any():
                #current_spec.plot()
                #plt.show()
                print("Skipping empty spectrum at position (x,y ) = (%i,%i)... (dec,ra) = (%s,%s)" %(x,y,dec,ra))
            #    print(current_spec.mask)
                #image_masked[y,x]=1
                #
                continue
            elif (y==0 and x==3) or (y==2 and x==4):
                    print("Storing spectrum at position (x,y ) = (%i,%i)... (dec,ra) = (%s,%s)" %(x,y,dec,ra))
                    plt.figure()
                    print(current_spec.data)
                    current_spec.plot()
                    plt.show()
                    current_spec.write("%i%ispectrum.fits" %(y,x))
        #    image_masked[y,x] =0
            #print("Processing spectrum at position %i,%i... (dec,ra) = (%s,%s)" %(x,y,ra,dec))

        #image_masked.write("masked.fits")
