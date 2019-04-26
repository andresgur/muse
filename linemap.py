# -*- coding: utf-8 -*
#Script to load, fit and plot spectra generated from MUSE cubes
# Author: Andres Gurpide Lasheras  andres.gurpide@gmail.com
#25-03-2019
    #!/usr/bin/env python3

#imports

import matplotlib.pyplot as plt
from mpdaf.obj import Spectrum, WaveCoord,iter_spe,Cube,deg2sexa
import argparse
import numpy as np
import os
import astropy.units as u
from math import *
import sys
sys.path.insert(0, '/home/agurpide/scripts/pythonscripts')
import plot_utils.plot_functions as pf
from muse_data import muse_utils


#fits the continuum of the input spectrum masking the wavelenghts arrays provided. The boolean mask_rest can be used to mask the spectrum from the maxixum wavelenght
#provided until the last wavelenght in the spectrum
def fit_continuum(spectrum, min_wa, max_wa, mask_rest):

    for index in np.arange(len(max_wa)):
        #mask wavelength ranges between lines
        #print("Masking region from %s to %s" %(min_wa[index],max_wa[index]))
        spectrum.mask_region(lmin=float(min_wa[index]),lmax=float(max_wa[index]),unit=u.angstrom,inside=True)

    #mask wavelenghts from the highest input wavelenght until the end of the spectrum
    if mask_rest:
        print("Masking everything above %s A in continuum fit"  %max_wa[-1])
        spectrum.mask_region(lmin=float(max_wa[-1]),lmax=None,unit=u.angstrom,inside=True)
    #print("Masked array")
    #print(spectrum.mask)
    print("Fitting continuum...")

    continuum = spectrum.poly_spec(1)

    spectrum.unmask()
    return continuum

#read arguments
ap = argparse.ArgumentParser(description='Spectrum fits to be loaded')
ap.add_argument("input_cubes",nargs='+',help="List of fits muse data cubes to be processed for the line map")
ap.add_argument("-l","--linefile",nargs='?',help="File with line limits to be masked for each linemap", default=os.environ['HOME'] +"/scripts/pythonscripts/muse_data/config_files/line_maps.txt",type=str)
ap.add_argument("-o","--outdir",nargs='?',help="Name of the output directory",default="linemaps")
args = ap.parse_args()

muse_cubes=args.input_cubes

line_file_path= args.linefile

outdir= args.outdir

plot_continuum = 0

line_names,min_wa,max_wa = muse_utils.load_map_file(line_file_path)

for cubefile in muse_cubes:

    if os.path.isfile(cubefile):
        print('Loading cube %s ...' %cubefile)
        original_cube = Cube(cubefile)
        print("Cube loaded successfully")

        continuum_cube = original_cube.clone(data_init=np.empty, var_init=np.zeros)

        continuum_cube.data.mask = original_cube.data.mask

        lines_cube = original_cube.clone(data_init=np.empty, var_init=np.zeros)

        lines_cube.data.mask = original_cube.data.mask

        #create continuum and lines cube
        for (current_spec,(y,x)),continuum_spec,line_spec in zip(iter_spe(original_cube,index=True),iter_spe(continuum_cube),iter_spe(lines_cube)):
            #skip masked or border pixels
            if not current_spec.data.any() or x==0 or y==0:

                #print("Skipping empty spectrum...")
                continue
            coord_sky = original_cube.wcs.pix2sky([y,x], unit=u.deg)
            dec, ra = deg2sexa(coord_sky)[0]
            print("Processing spectrum at position (x,y ) = (%i,%i)... (dec,ra) = (%s,%s)" %(x,y,dec,ra))
            continuum_spec[:] = fit_continuum(current_spec,min_wa,max_wa,True)

            if plot_continuum:
                plt.figure()
                current_spec.plot()
                continuum_spec.plot(color='red')
                plt.show()


            line_spec[:] = current_spec[:] - continuum_spec[:]

        if plot_continuum:
            plt.show()
        #correct for negative subtractions, where the continuum happen to be above the wavelenght pixel

        line_spec[np.where(line_spec<0)] = 0
        #create image for line fluxes maps
        original_image = original_cube.sum(axis=0)
        #free memory
        original_cube = None
        if not os.path.isdir(outdir):
            os.mkdir(outdir)
        #create continuum plot
        print("Storing white light continuum image...")
        white_light_image_cont = continuum_cube.sum(axis=0)
        white_light_image_cont.write(outdir + "/continuum_image.fits")#create line fluxes image

        print("Storing lines white light image...")
        lines_image = lines_cube.sum(axis=0)
        lines_image.write(outdir + "/line_image.fits")

        #create line maps for every line
        line_index = 0

        for wavelength in line_names:
            line_map = original_image.clone(data_init=np.empty)
            line_map.mask = original_image.mask
            print("\n Creating line map for %s A" %wavelength)
            print("Integrating from %s to %s A" %(min_wa[line_index],max_wa[line_index]))
            for line_spec,(p,q) in iter_spe(lines_cube,index=True):
                #skip masked pixels
                if line_spec.data is None or all(line_spec.var==0):
                    line_map.mask[p,q] = True
                    continue

                line_flux_quantity,error_quantity = line_spec.integrate(float(min_wa[line_index]),float(max_wa[line_index]),unit=u.angstrom)
                line_flux = line_flux_quantity.value
                units = line_flux_quantity.unit.to_string()
                line_flux_error = error_quantity.value
                line_flux_error_units = error_quantity.unit.to_string()

                print("Line flux %.2f (%s) with error %.2f (%s) for spectrum at pixel position %i,%i" %(line_flux,units,line_flux_error,line_flux_error_units,p,q))

                #sum,error = line_spec.sum(float(min_wa[line_index]),float(max_wa[line_index]),unit=u.angstrom)
                #step  = line_spec.get_step(unit=line_spec.wave.unit)
                #line_flux_num =  sum*step
                #print("Line flux by numeric integration: %.2f with error %.2f " %(line_flux_num,error*step))

                line_map[p,q] = line_flux
            line_index = line_index+1
            #store line map
            line_map.write(outdir + "/%s_map.fits" %wavelength)
