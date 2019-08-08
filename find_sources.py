# -*- coding: utf-8 -*-
# Script to retrieve all sources from SIMBAD within a given region file
# Author: Andres Gurpide Lasheras  andres.gurpide@gmail.com
# 25-05-2019
# !/usr/bin/env python3
# imports
from astroquery.simbad import Simbad
from astropy.coordinates import SkyCoord, Distance
import astropy.units as u
from astropy.io import fits
from regions import read_ds9, write_ds9, DS9Parser
import logging
from math import sqrt
import os
from astropy.time import Time
from astropy.table import Column
from astroquery.gaia import Gaia
import numpy as np
import sys
import argparse


def search_sources(coords, radius, obs_epoch=Time('J2000'), equinox=2000, catalog='gaia'):
    if args.catalog == 'simbad':
        my_simbad = Simbad()
        # add the object type to the output results
        my_simbad.add_votable_fields('otype')
        result_table = my_simbad.query_region(coords,
                                              radius=radius, epoch=obs_epoch.jyear_str, equinox=equinox)
    elif args.catalog == 'gaia':
        j = Gaia.cone_search(coordinate=coords, radius=radius)
        result_table = j.get_results()
        pm_ra_cosdec = (result_table['pmra'] * np.cos(result_table['dec'].to("rad"))) * u.mas / u.yr

        # write coordinates of the sources propagated to the new epoch
        c = SkyCoord(ra=result_table['ra'], dec=result_table['dec'],
                     pm_ra_cosdec=pm_ra_cosdec, pm_dec=result_table['pmdec'],
                     radial_velocity=result_table["radial_velocity"],
                     obstime=Time(result_table['ref_epoch'], format='decimalyear'))

        # get position of the stars in the observation position
        c.apply_space_motion(new_obstime=obs_epoch)
        result_table.add_column(Column(data=c.ra.value, format='E', name='ra_%s' % obs_epoch.jyear_str))
        result_table.add_column(Column(data=c.dec.value, format='E', name='dec_%s' % obs_epoch.jyear_str))
    return result_table


def get_region_coords(region):
    # Circular region
    if type(region).__name__ == 'CircleSkyRegion':
        centerra = region.center.fk5.ra
        centerdec = region.center.fk5.dec
        radius = region.radius

    elif type(region).__name__ == 'EllipseSkyRegion':
        centerra = region.center.fk5.ra
        centerdec = region.center.fk5.dec
        radius = sqrt(region.height.to(u.arcsec) ** 2 + region.width.to(u.arcsec) ** 2)

    coords, radius = SkyCoord(ra=centerra.value, dec=centerdec.value, unit=(centerra.unit, centerdec.unit), frame='fk5')
    return coords, radius


# ---------------------------------------------------------------------------------
# main
# read arguments
parser = argparse.ArgumentParser(description='Regions must be in WCS in angles and radius in arcseconds')
parser.add_argument("--regions", help="Circular regions where to look for sources in the catalog. Defaults to image center if any", nargs='*', type=str)
parser.add_argument("fits", help="Fits image file with epoch, field of view size, coordinates, etc", nargs=1, type=str)
parser.add_argument("-c", "--catalog", type=str, choices=["simbad", "gaia"], default="gaia", help='Catalog search')
parser.add_argument("-r", "--radius", help='Radius search in arcseconds if no region is provided, default 70 arcsec (~ MUSE field of view in WFM)', type=float, nargs='?', default=70)

args = parser.parse_args()

# simbad radius in arseconds
simbad_radius = 5 * u.arcsec
simbad_color = 'blue'
input_fits = args.fits[0]
# default epoch for the search in case the fits file has not date-obs attribute
default_epoch = 'J2000'
default_equinox = "2000"

if args.radius is not None:
    radius = args.radius * u.arcsec

scriptname = os.path.basename(__file__)
# logger

logger = logging.getLogger(scriptname)
logger.setLevel(logging.DEBUG)
out_format = logging.Formatter('%(name)s - %(levelname)s - %(message)s')

# handler
stream_handler = logging.StreamHandler()
stream_handler.setLevel(logging.DEBUG)
stream_handler.setFormatter(out_format)
logger.addHandler(stream_handler)


if os.path.isfile(input_fits):
    # get primary header
    header = fits.getheader(input_fits, 0)
    # observation date in julian days
    try:
        obs_epoch = Time(header['DATE-OBS'])
        equinox = header['EQUINOX']

    except KeyError:
        logger.warning("Input %s fits file has not DATE-OBS keyword. Setting date-obs to default %s" % (input_fits, default_epoch))
        obs_epoch = Time(default_epoch)
        equinox = default_equinox
    try:
        radecsys = (header["RADECSYS"].lower())
        ra_center = header["RA"]
        dec_center = header["DEC"]
    except KeyError:
        logger.warning("No coordinate system found in %s fits file" % (input_fits))
        ra_center = None

else:
    logger.warning("Input %s file does not exist!" % input_fits)
    sys.exit()

if args.regions is not None:
    for region_file in args.regions:
        if not os.path.isfile(region_file):
            logger.warning("File %s not found" % region_file)
            continue
        out_region = "simbad_%s" % (region_file)
        read_regions = read_ds9(region_file)
        logger.info("Loaded %i region(s) from %s" % (len(read_regions), region_file))
        reg_string = 'global color=%s\nfk5\n' % simbad_color
        for index, region in enumerate(read_regions):

            coords, radius = get_region_coords(region, radecsys, obs_epoch.jyear_str, equinox, args.catalog)

            logger.debug("Radius search: %.3f" % radius.value)

            result_table = search_sources(coords, radius, obs_epoch, equinox, args.catalog)

            if result_table is None:
                logger.info("No objects found for region %s in file %s " % (region, region_file))
            else:
                logger.info("%i objects found for region %s in file %s" % (len(result_table), region, region_file))

                out_table = "%i%s" % (index, out_region.replace(".reg", ".csv"))
                result_table.write(out_table, overwrite=True, format='csv')
                
                for row in result_table:
                    ra = ("%s" % row['RA']).replace(" ", ":")
                    dec = ("%s" % row['DEC']).replace(" ", ":")
                    otype = ("%s" % row["OTYPE"]).replace("b\'", "").replace('\'', '')
                    reg_string += 'circle(%s,%s,%.2f\") # width=2 text={%s} \n' % (ra, dec, simbad_radius.value, otype)

                    if reg_string != "":
                        ds9_regions = DS9Parser(reg_string).shapes.to_regions()
                        print(ds9_regions)
                        write_ds9(ds9_regions, out_region)

elif ra_center is not None:
    coords = SkyCoord(ra=ra_center, dec=dec_center, unit=(u.deg, u.deg), frame=radecsys)
    logger.debug("Radius search: %.3f" % radius.value)

    result_table = search_sources(coords, radius, obs_epoch, equinox, args.catalog,)

    if result_table is None:
        logger.info("No objects found in field of view of %s with radius %.2f" % (input_fits, radius))
    else:
        logger.info("%i objects found in field of view of %s with radius %.2f" % (len(result_table), input_fits, radius.value))

        out_table = "%ssources%s" % (args.catalog, os.path.basename(input_fits).replace(".fits", ".csv"))
        result_table.write(out_table, overwrite=True, format='csv', delimiter='\t')
        print(result_table.colnames)
