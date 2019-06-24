# -*- coding: utf-8 -*
# Determines the FWHM of a cube or image by fitting a MOFFAT profile to the image
# Author: Andres Gurpide Lasheras  andres.gurpide@gmail.com
# 09-05-2019
# !/usr/bin/env python3

# imports
from mpdaf.obj import Cube, Image, iter_ima
import astropy.units as u
import matplotlib.pyplot as plt
import argparse
import os
import muse_utils
from math import sqrt
import sys
import logging
import pyregion
import numpy as np
from astropy.io import fits
from lmfit.models import PolynomialModel


def create_result_file(parentfile=''):
    hdr = fits.Header()
    hdr['AUTHOR'] = 'Andres'
    hdr['CREATOR'] = parentfile
    hdr['COMMENT'] = 'Fits file storing the values of the FWHM as a function of wavelength for each computation'
    primary_hdu = fits.PrimaryHDU(header=hdr)
    return fits.HDUList([primary_hdu])


def write_outputeer(wavelengths, radius_fractions, parentfile='', outname='fluxfraction.fits'):
    hdr = fits.Header()
    hdr['AUTHOR'] = 'Andres'
    hdr['CREATOR'] = parentfile
    hdr['COMMENT'] = 'Fits file storing the values of the energy flux contained within a certain radius as a function of wavelength'
    primary_hdu = fits.PrimaryHDU(header=hdr)
    hdul = fits.HDUList([primary_hdu])
    for index, l in enumerate(wavelengths):
        hdu = fits.BinTableHDU.from_columns([fits.Column(name='Radius', format='E', array=radius_fractions[index][0], unit="arcsec"),
                                            fits.Column(name='Flux fraction', format='E', array=radius_fractions[index][1], unit="")], name="%.2f" % l)

        hdul.append(hdu)

    hdul.writeto(outname, overwrite=True)


def parallel_eer(image, center=None, unit_center='deg'):
    """Computes the flux fraction contained within a radius from the given center or from the center of the image."""
    radius, flux_fraction = image.eer_curve(center=center, unit_center=unit_center)
    return radius, flux_fraction


def parallel_fit(image, verbose_level=False, circular=True, factor=2, profile='moffat', center=None, unit_center='deg', fwhm=None):
    if profile == 'moffat':
        return image.moffat_fit(cont=1, fit_back=True, weight=True, verbose=verbose_level, full_output=True,
                                circular=circular, plot=False, factor=factor, center=center, unit_center=unit_center, fwhm=fwhm, n=1)
    elif profile == 'gauss':
        return image.gauss_fit(cont=1, fit_back=True, weight=True, verbose=verbose_level, full_output=True,
                               circular=circular, plot=False, factor=factor, center=center, unit_center=unit_center, fwhm=fwhm)


def derive_fwhm(cube, center=None, threads=3):
    """Derive the FWHM PSF of the cube as a function of wavelength by computing the flux ratio around the center or the center of the image.
    Parameters:
    cube : the cube whose PSF is to be derived.
    center : the center for the calculation.
    """
    logger.debug("Deriving FWHM from raw image")
    radius_flux_fractions = cube.loop_ima_multiprocessing(f=parallel_eer, cpu=threads, center=center)
    return radius_flux_fractions


def estimate_fwhm(cube, center_guess=None, threads=3, fwhm_guess=None, circular=True):
    """Compute the FWHM PSF of the cube as a function of wavelength by fitting a 2D profile to each wavelength slice. A region file around a star can be provided to mask the rest of the cube.
    Parameters:
    cube : the cube whose PSF is to be estimated.
    center_guess : the initial guess for the center
    """
    # get fit ino
    logger.info("Fitting images...")
    fit_profiles = cube.loop_ima_multiprocessing(f=parallel_fit, cpu=threads, verbose=True, circular=circular, factor=1,
                                                 verbose_level=True, center=center_guess, fwhm=guessed_fwhm)
    out_shape = (cube.shape[0], fit_profiles[0].ima.shape[0], fit_profiles[0].ima.shape[1])
    fit_cube = Cube(wcs=cube.wcs, wave=cube.wave, data=np.zeros(out_shape), unit=fit_profiles[0].ima.unit,
                    data_header=cube.data_header.copy(), primary_header=cube.primary_header.copy())

    for l in range(fit_cube.shape[0]):
        fit_cube[l, :, :] = fit_profiles[l].ima
    return fit_profiles, fit_cube


# ------------------------------
# main
# read arguments
ap = argparse.ArgumentParser(description='Determines the FWHM of the cube or the energy encircled in an image')
ap.add_argument("input_file", nargs=1, help="Image or cube whose FWHM is to be determined")
ap.add_argument('-p', '--profile', type=str, choices=["moffat", "gaussian"], default="moffat")
ap.add_argument("-e", "--error_map", nargs=1, help="Error image for the fitting process if any", required=False)
ap.add_argument("-r", "--fitting_region", nargs='+', help="Regions for the image fitting. More than one to compare between them.", default=None, type=str)
ap.add_argument("-l", "--lbinning", nargs='?', help="Wavelength binning for the cube FWHM computation", type=int)
ap.add_argument("-t", "--threads", nargs='?', help="Number of threads to be used in the parallel image fits", default=2, type=int)
ap.add_argument("--elliptical", help="Use a circular profile or elliptical in the fitting. Default true (circular)", action='store_false')
ap.add_argument("--fit", help="Attempt to compute the FWHM as a function of wavelength also by fitting. Default false", action='store_true')

args = ap.parse_args()

circular_fit = args.elliptical

out_dir = 'fwhm_fit'
cube_file_type = 'cube'
image_file_type = 'image'
tolerance = 0.001
plot_out = 1
plot_errors = 0

# polynomial initial fit params
degree = 2
c0 = 0.65
c1 = -4.5 * 10 ** (-5)
c2 = 0
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


if not os.path.isdir(out_dir):
    os.mkdir(out_dir)

if args.fitting_region is not None:
    region_colors = muse_utils.create_color_array(len(args.fitting_region), cmap='PiYG')
    linestyles = muse_utils.get_linestyle_array(len(args.fitting_region))

input_file = args.input_file[0]
# skip files that are not found
if not os.path.isfile(input_file):
    logger.warn("File %s not found" % input_file)
    sys.exit()

cube_center_figure = None
eer_circle_fig = None

try:
    # is the input a cube?
    input_type = cube_file_type
    logger.info("Trying to load input %s as a %s..." % (input_file, input_type))
    input_mpdaf = Cube(input_file)
    logger.info("Loading successful")

    hdul_out = create_result_file(input_file)

    logger.info("Deriving theoretical model (Tokovinin et al 2002)")
    atm_seeing = muse_utils.get_atm_seeing(input_mpdaf)

    if atm_seeing is not None:

        c0 = atm_seeing.value
        psf_fwhm = muse_utils.tokovinin_model(atm_seeing, input_mpdaf.wave.coord() * input_mpdaf.wave.unit)
        model = PolynomialModel(degree=2)
        params = model.make_params(c0=c0, c1=c1, c2=c2)
        tokovinin_fit = model.fit(psf_fwhm, params, x=input_mpdaf.wave.coord())

        fwhm_wa = tokovinin_fit.eval(tokovinin_fit.params, x=input_mpdaf.wave.coord())

        avg_fwhm = np.mean(fwhm_wa)
        logger.info("MeanFWHM: %.3f arcsec" % avg_fwhm)

        tokovinin_results = fits.BinTableHDU.from_columns([fits.Column(name='FWHM', format='E', array=fwhm_wa, unit="arcsec"),
                                                          fits.Column(name='Wavelength', format='E', array=input_mpdaf.wave.coord(), unit=input_mpdaf.wave.unit.name)],
                                                          name="Tokovinin")

        tokovinin_results.header['COMMENT'] = 'Tokovinin PSF model Tokovining et al 2002'
        tokovinin_results.header['Model'] = (model.name)
        tokovinin_results.header['Degree'] = (model.poly_degree)
        tokovinin_results.header['c0'] = ("%.5f" % tokovinin_fit.params.get("c0"), 'c0 value of the polynomial fit')
        tokovinin_results.header['c1'] = ("%.5f" % tokovinin_fit.params.get("c1"), 'c1 value of the polynomial fit')
        tokovinin_results.header['c2'] = ("%.5f" % tokovinin_fit.params.get("c2"), 'c2 value of the polynomial fit')
        tokovinin_results.header['FWHMa'] = ("%.5f" % avg_fwhm, 'Mean FWHM value in arcsec')
        tokovinin_results.header['FWHMp'] = ("%.5f" % (avg_fwhm * input_mpdaf.primary_header["HIERARCH ESO OCS IPS PIXSCALE"]), 'Mean FWHM value in pixels')

        hdul_out.append(tokovinin_results)

        cube_psf_figure = plt.figure(figsize=(16.0, 10.0))
        cube_psf_ax = cube_psf_figure.add_subplot(111)

        cube_psf_ax.set_xlabel("$\lambda$ (%s)" % u.angstrom, fontsize=20)
        cube_psf_ax.set_ylabel("FWHM (%s)" % u.arcsec, fontsize=20)
        cube_psf_ax.ticklabel_format(axis='both', style='plain')
        muse_utils.format_axis(cube_psf_ax)

        cube_psf_ax.errorbar(input_mpdaf.wave.coord(), psf_fwhm, color='red', markersize=1, linewidth=1, label='PSF model')
        #cube_psf_ax.errorbar(input_mpdaf.wave.coord(), out_fit.best_fit, color='green', ls='-', linewidth=2)
        logger.info("Fit parameters for Tokovinin model: \n")
        logger.info(tokovinin_fit.fit_report())

    else:
        logger.warning("Seeing keyword not found in header. Will not derive theoretical model for the seeing")

    original_wavelengths = input_mpdaf.wave.coord()

    if args.lbinning is not None:
        rebinning = (input_mpdaf.wave.get_step(unit=u.angstrom) * args.lbinning)
        logger.info("Rebinning cube to %.3f A" % rebinning)
        input_mpdaf.rebin((args.lbinning, 1, 1), inplace=True, margin='origin')

    line_colors = muse_utils.create_color_array(input_mpdaf.shape[0], cmap='plasma')

    cube_center_figure = plt.figure(figsize=(16.0, 10.0))
    cube_center_ax = cube_center_figure.add_subplot(111)

    cube_center_ax.set_xlabel("Ra (deg)", fontsize=20)
    cube_center_ax.set_ylabel("Dec (deg)", fontsize=20)
    cube_center_ax.legend(loc='best', prop={'size': 14}, frameon=False)
    muse_utils.format_axis(cube_center_ax)

    eer_circle_fig = plt.figure(figsize=(16.0, 10.0))

    ax_eer_circle = eer_circle_fig.add_subplot(111)
    muse_utils.format_axis(ax_eer_circle)
    ax_eer_circle.set_xlabel('Radius (asec)', fontsize=20)
    ax_eer_circle.set_ylabel('Flux fraction', fontsize=20)

    sky_res, sky_res_std = muse_utils.get_sky_res(input_mpdaf)
    if sky_res is not None:
        # cube_psf_ax.axhline(y=sky_res, label='Avg FWHM', ls='--', color='black')
        guessed_fwhm = (sky_res, sky_res)
    else:
        guessed_fwhm = None


except ValueError:

    # it is an image
    input_type = image_file_type
    logger.info("Trying to load input %s as an %s..." % (input_file, input_type))
    input_mpdaf = Image(input_file)

    eer_circle_fig = plt.figure(figsize=(16.0, 10.0))

    ax_eer_circle = eer_circle_fig.add_subplot(111)
    muse_utils.format_axis(ax_eer_circle)
    ax_eer_circle.set_xlabel('Radius (asec)', fontsize=20)
    ax_eer_circle.set_ylabel('Flux fraction', fontsize=20)

if args.fitting_region is not None:

    for region, color, linestyle in zip(args.fitting_region, region_colors, linestyles):

        stem = os.path.basename(region)
        stem_out = stem.replace(".reg", ".fits")

        # cube
        if input_type == cube_file_type:

            cube = input_mpdaf.copy()

            if args.lbinning is not None:
                cube.primary_header["HISTORY"] = "Rebinned to %.2f" % rebinning

            if region is not None:
                logger.info("Masking images with region %s ..." % region)

                # mask every image
                img = cube[0, :, :]
                # create spatial mask out of region file
                mask = muse_utils.region_to_mask(img, region)
                for ima in iter_ima(cube):
                    ima.mask[np.where(mask)] = True
                # get rid of the masked values
                logger.info("Cropping cube")
                cube.crop()

            center_region = pyregion.open(region)[0]

            fit_center = (center_region.coord_list[1], center_region.coord_list[0])

            if args.fit:
                logger.info("Fitting region %s. Circular profile: %r" % (stem, circular_fit))
                fit_profiles, fit_cube = estimate_fwhm(cube, center_guess=fit_center, threads=args.threads, fwhm_guess=guessed_fwhm, circular=args.elliptical)
                fit_cube.write('%s/fitcube%s' % (out_dir, stem_out))
                # filter out bad fittings
                good_fit_profiles, good_wavelengts = zip(*((profile, wa) for wa, profile in zip(fit_cube.wave.coord(), fit_profiles) if profile.err_fwhm[0] < 10**3
                                                         and sqrt((profile.center[0] - fit_center[0])**2 + (profile.center[1] - fit_center[1])**2) < tolerance and np.isfinite(profile.flux) and profile.flux > 0))
                # plot psf FWHM
                cube_psf_ax.errorbar(x=good_wavelengts, y=[profile.fwhm[1] for profile in good_fit_profiles], yerr=[profile.err_fwhm[1] for profile in good_fit_profiles], ls="-", linewidth=0.5,
                                     markersize=4, color=color, label='%s' % stem, fmt='.', elinewidth=1, errorevery=4)
                if not args.elliptical:
                    cube_psf_ax.errorbar(x=good_wavelengts, y=[profile.fwhm[0] for profile in good_fit_profiles], yerr=[profile.err_fwhm[0] for profile in good_fit_profiles], ls="-", linewidth=0.5,
                                         markersize=4, color=color, label='%s' % stem, fmt='.', elinewidth=1, errorevery=4)
                # plot PSF centers
                cube_center_ax.errorbar(x=[profile.center[1] for profile in good_fit_profiles], y=[profile.center[0] for profile in good_fit_profiles],
                                        markersize=4, color=color, label='%s x' % stem, fmt='.', ls=None)

                fit_cube = None  # free memory
                model = PolynomialModel(degree=2)
                params = model.make_params(c0=c0, c1=c1, c2=c2)
                psf_fit = model.fit([profile.fwhm[0] for profile in good_fit_profiles], params, x=good_wavelengts, weights=[profile.err_fwhm[0] for profile in good_fit_profiles])
                logger.info("Fit parameters for region %s: \n" % region)
                logger.info(psf_fit.fit_report())
                avg_beta = np.mean([profile.n for profile in good_fit_profiles])
                avg_fwhm_y = np.mean([profile.fwhm[0] for profile in good_fit_profiles])
                avg_fwhm_x = np.mean([profile.fwhm[1] for profile in good_fit_profiles])
                logger.info("Averga beta index: %.2f" % avg_beta)
                logger.info("Averga FWHM y: %.3f" % avg_fwhm_y)
                logger.info("Averga FWHM x: %.3f" % avg_fwhm_x)

                fwhm_wa = psf_fit.eval(psf_fit.params, x=original_wavelengths)

                psf_fit_results = fits.BinTableHDU.from_columns([fits.Column(name='FWHM', format='E', array=fwhm_wa, unit="arcsec"),
                                                                fits.Column(name='Wavelength', format='E', array=original_wavelengths, unit=input_mpdaf.wave.unit.name)],
                                                                name="%s" % stem)
                psf_fit_results.header['COMMENT'] = 'Fit to region %s' % region
                psf_fit_results.header['Model'] = (model.name)
                psf_fit_results.header['Degree'] = (model.poly_degree)
                psf_fit_results.header['c0'] = ("%.5f" % psf_fit.params.get("c0"), 'c0 value of the polynomial fit')
                psf_fit_results.header['c1'] = ("%.5f" % psf_fit.params.get("c1"), 'c1 value of the polynomial fit')
                psf_fit_results.header['c2'] = ("%.5f" % psf_fit.params.get("c2"), 'c2 value of the polynomial fit')
                psf_fit_results.header['FWHMa'] = ("%.5f" % (avg_fwhm_x), 'Mean FWHM value in arcseconds')
                psf_fit_results.header['FWHMp'] = ("%.5f" % (avg_fwhm_x * input_mpdaf.primary_header["HIERARCH ESO OCS IPS PIXSCALE"]), 'Mean FWHM value in pixels')

                psf_fit_results.header['Beta'] = ("%.5f" % avg_beta, 'Mean Beta value from the Moffat fits')

                hdul_out.append(psf_fit_results)
                #cube_psf_ax.errorbar(x=good_wavelengts, y=out_fit.best_fit, ls='--', color=color)

            # the cube is cut around the region at this point so no need for cutting it
            radius_fractions = derive_fwhm(cube, center=fit_center, threads=args.threads)

            cube.write("%s/rebinned_%s" % (out_dir, stem_out))

            write_outputeer(cube.wave.coord(), radius_fractions, input_file, "%s/eer_%s" % (out_dir, stem_out))

            for l, line_color in enumerate(line_colors):
                ax_eer_circle.errorbar(x=radius_fractions[l][0], y=radius_fractions[l][1], color=color, linewidth=0.5, ls=linestyle, label="%.2f $\AA$" % cube.wave.coord()[l])

        # image
        elif input_type == image_file_type:

            if args.error_map is not None:
                errormap = args.error_map[0]
                if os.path.isfile(errormap):
                    logger.info('Loading error map %s ...' % errormap)
                    e_map = Image(errormap)
                    # set the variance of the input map to the values of the error map
                    input_mpdaf.var = e_map.data
                else:
                    logger.warning("Error map %s not found" % errormap)

            if 'flux' in input_file:
                zlabel = "1e-20 erg/cm  $^2$/s"
            else:
                zlabel = ""

            # cut image from region
            logger.info("Cropping %s with region %s" % (input_type, region))
            mask = muse_utils.region_to_mask(input_mpdaf, region)
            image = input_mpdaf.copy()
            image.mask[np.where(mask)] = True
            image.crop()

            # create figures
            img_figure = plt.figure(figsize=(16.0, 10.0))
            im_ax = img_figure.add_subplot(111, projection=input_mpdaf.wcs.wcs)
            im_ax.minorticks_on()
            image.plot(ax=im_ax, scale='linear', colorbar='v', show_xlabel=True, show_ylabel=True)
            im_ax.set_xlabel('Ra', fontsize=20)
            im_ax.set_ylabel('Dec', fontsize=20)

            im = im_ax.images
            cb = im[-1].colorbar
            cb.ax.set_ylabel(zlabel, fontsize=14)
            image.write("%s%s" % (os.path.basename(region).replace(".reg", ""), os.path.basename(input_file).replace(".fits", ".pdf")))

            profile_figure = plt.figure(figsize=(16.0, 10.0))
            profile_ax = profile_figure.add_subplot(111, projection=image.wcs.wcs)
            profile_ax.minorticks_on()

            # infer center, fwhm initial guess
            center_region = pyregion.open(region)[0]
            fit_center = (center_region.coord_list[1], center_region.coord_list[0])
            sky_res, sky_res_std = muse_utils.get_sky_res(input_mpdaf)
            if sky_res is not None:
                guessed_fwhm = (sky_res, sky_res)
                ax_eer_circle.axvline(x=sky_res, label='FWHM - cube', ls='--', alpha=0.5, color='blue')
            else:
                guessed_fwhm = None
            if args.profile == 'moffat':
                profile_fit = image.moffat_fit(center=fit_center, cont=0, fit_back=False, weight=True, full_output=True,
                                               plot=False, fwhm=guessed_fwhm, circular=args.elliptical, factor=2)
                if profile_fit.cont < 0:
                    logger.warning('Fit for region %s was not successful' % region)
                    continue
                if profile_fit.err_fwhm[0] > 10**6:
                    logger.warning("Fit parameters are unbound")

            profile_fit.ima.plot(ax=profile_ax, scale='linear', colorbar='v', show_xlabel=True, show_ylabel=True)

            profile_ax.set_xlabel('Ra')
            profile_ax.set_ylabel('Dec')
            im = profile_ax.images
            cb = im[-1].colorbar
            cb.ax.set_ylabel(zlabel, fontsize=14)

            # store gaussian fit.
            fit_out_name = "%s%s" % (os.path.basename(region).replace(".reg", ""), os.path.basename(input_file).replace(".fits", "%s_fit.fits" % args.profile))
            profile_fit.ima.write(out_dir + "/" + fit_out_name)

            fit_res = Image(data=image.data - profile_fit.ima.data)

            res_out_name = fit_out_name.replace("fit.fits", "res.fits")
            fit_res.write(out_dir + "/" + res_out_name)

            # compute energy enclosed as a function of radius
            radius, energy_values = image.eer_curve(cont=profile_fit.cont, center=profile_fit.center)
            radius_fit, energy_values_fit = profile_fit.ima.eer_curve(cont=profile_fit.cont, center=profile_fit.center)

            ax_eer_circle.scatter(radius, energy_values, color='green', label='%s' % os.path.basename(region).replace(".reg", ""))
            ax_eer_circle.scatter(radius_fit, energy_values_fit, color='blue', label='%s fit' % os.path.basename(region).replace(".reg", ""))

            ax_eer_circle.axvline(x=profile_fit.fwhm[0], ls='--', alpha=0.5, color='blue')
            ax_eer_circle.axvline(x=profile_fit.fwhm[1], ls='--', alpha=0.5, color='blue')

if cube_center_figure is not None:
    cube_center_figure.savefig('%s/psf_center%s' % (out_dir, os.path.basename(input_file).replace('.fits', ".pdf")))
    cube_psf_ax.legend(loc='best', prop={'size': 14}, frameon=False)
    cube_psf_figure.savefig('%s/fwhm_wa%s' % (out_dir, os.path.basename(input_file).replace('.fits', ".pdf")))

if eer_circle_fig is not None:

    ax_eer_circle.legend(loc='best', prop={'size': 14}, frameon=False)
    muse_utils.remove_legend_repetitions(ax_eer_circle)
    eer_circle_fig.savefig('%s/encircled_energy.pdf' % out_dir)

hdul_out.writeto("%s/fwhm%s" % (out_dir, input_file), overwrite=True)
logger.info("Fits file written to %s/fwhm%s" % (out_dir, input_file))

plt.show()
