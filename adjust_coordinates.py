# Script to adjust the coordinates of a cube from an input image by cross-correlating it( Preferabley HST)
# Created the 10/06/2019 by Andres Gurpide
# imports
import argparse
from mpdaf.obj import Cube, Image
import os
import logging

# read arguments
ap = argparse.ArgumentParser(description='Adjust coordinates from HST image')
ap.add_argument("input_cube", nargs=1, help="Cube to be corrected", type=str)
ap.add_argument("--hst", nargs=1, help="HST image", type=str)
ap.add_argument("-o", "--outdir", nargs='?', type=str, help="Output directory", default="coordadjusted")

# parse args
args = ap.parse_args()
nsigma = 3

# create output dir
if not os.path.isdir(args.outdir):
    os.mkdir(args.outdir)

# logger
scriptname = os.path.basename(__file__)

logger = logging.getLogger(scriptname)
logger.setLevel(logging.DEBUG)
out_format = logging.Formatter('%(name)s - %(levelname)s - %(message)s')

# handler
stream_handler = logging.StreamHandler()
stream_handler.setLevel(logging.DEBUG)
stream_handler.setFormatter(out_format)
logger.addHandler(stream_handler)


input_cube = args.input_cube[0]
hst_image = args.hst[0]

outcube = '%s/corr%s' % (args.outdir, input_cube)
wlname = "%s/wlightcorr%s" % (args.outdir, input_cube)

if os.path.isfile(input_cube) and os.path.isfile(hst_image):
    logger.info("Loading cube %s" % input_cube)
    cube = Cube(input_cube)
    logger.info("Loading succesful")
    logger.info("Creating white light image...")
    white_image = cube.sum(axis=0)
    hst_im = Image(hst_image)
    logger.info("Estimating coordinate offset...")
    dy, dx = white_image.estimate_coordinate_offset(hst_im, nsigma=nsigma)
    white_image.adjust_coordinates(hst_im, nsigma=nsigma, inplace=True)
    white_image.write(wlname)
    # update cube WCS header
    cube.wcs.set_crpix1(white_image.wcs.get_crpix1())
    cube.wcs.set_crpix2(white_image.wcs.get_crpix2())
    cube.wcs.set_cd(white_image.wcs.get_cd())
    pixel_scale = cube.primary_header['HIERARCH ESO OCS IPS PIXSCALE']
    cube.primary_header["ASTRO_CORR_PX"] = ("%.2f, %.2f" % (dy, dx), 'Coordinate shift dy, dx in pixels')
    cube.primary_header["ASTRO_CORR_AS"] = ("%.2f, %.2f" % ((dy * pixel_scale), (dx * pixel_scale)), 'Coordinate shift dy, dx in arcseconds')

    cube.write(outcube)
else:
    logger.error("Input cube %s or HST image %s not found" % (input_cube, hst_image))
