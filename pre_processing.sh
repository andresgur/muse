# paths
hst_image="$HOME/optical_data/NGC0625/MAST_2019-06-12T1101/HLA/hst_08708_02_wfpc2_f555w_pc/hst_08708_02_wfpc2_f555w_pc_drz.fits" # path to HST image to adjust coordintes
input_cube="ADP.2016-06-20T16_43_08.136.fits" # path to the input cube
chandra_image="/home/agurpide/x_ray_data/NGC1672/chandradata/5932/detectedsources/sourceimage_broadband.fits"
ulx_region="/home/agurpide/x_ray_data/NGC1672/chandradata/5932/repro/ixo27.reg"
run_Zap=false
clean_nan=false
ds9_command="$DS9_PATH/ds9 -cmap Heat -mode region"
exec_python="python $HOME/scripts/pythonscripts/muse"
redshift=0.004464

# stop script if something goes wrong
set -e

# remove NaN values
if $clean_nan ; then
	$exec_python/removenan.py $input_cube -o nancleancubes
	cd nancleancubes
	input_cube=nanclean$input_cube
fi

# adjust coordinates; white light image will be created
if [ $hst_image != "" ]; then
	$exec_python/adjust_coordinates.py $input_cube --hst $hst_image -o coordadjusted
	input_cube=corr$input_cube
	white_light=wlight$input_cube
	cd coordadjusted
	
else
	# create white light image otherwise
	white_light_image="from mpdaf.obj import Cube;cubein = '${input_cube}'; cube = Cube(cubein);img = cube.sum(axis=0); print('Saving image to white%s' % cubein);  img.write('white%s' % cubein)"
	python -c "${white_light_image}"
	white_light=white$input_cube

fi 

$exec_python/find_sources.py -c gaia $white_light

echo -e "\e[36m Find a star and save the region file for PSF computation.\e[0m"

catalog_out=${white_light//fits/csv}

$ds9_command $input_cube -catalog import csv gaiasources$catalog_out $hst_image 

read -p "Did you find a star? y/n"

if [[ $REPLY =~ ^[Yy]$ ]]; then
	echo "Input region file"
	read star_reg
	# determine PSF by fitting a Moffat profile to it
	$exec_python/determine_psf.py $input_cube --fit --fitting_region $star_reg -l 100 -t 6
else
# use Tokovinin mode instead
	$exec_python/determine_psf.py $input_cube -t 6
fi


echo "Input average FWHM in arcsec:"

read fwhm

$exec_python/findstars.py $white_light --fwhm $fwhm -t 3 -s 2

echo "Copying catalog regions and saving them as zap_mask.reg"


zap_region_size=$(echo "scale=3; 5*$fwhm" | bc)

zap_regions_out="zap_mask"${white_light//fits/reg}


echo "Converting dao sources to ZAP regions with radius $zap_region_size arcsec"

rm zap_mask.sym

printf "condition\tshape\tcolor\twidth\tdash\tfont\tfontsize\tfontweight\tfontslant\ttext\tsize\tsize2\tunits\tangle\n
---------\t-----\t-----\t-----\t----\t----\t--------\t----------\t---------\t----\t----\t-----\t-----\t-----\n" >> zap_mask.sym

printf "\tcircle\tcyan\t2\t0\thelvetica\t10\tnormal\troman\t\t$zap_region_size\t\tarcsec" >> zap_mask.sym

$ds9_command $input_cube -catalog import csv allstars$catalog_out -catalog symbol load zap_mask.sym -catalog regions -regions save $zap_regions_out -exit

dao_find_reg="daofind"${white_light//fits/reg}

rm $dao_find_reg

printf "condition\tshape\tcolor\twidth\tdash\tfont\tfontsize\tfontweight\tfontslant\ttext\tsize\tsize2\tunits\tangle\n
---------\t-----\t-----\t-----\t----\t----\t--------\t----------\t---------\t----\t----\t-----\t-----\t-----\n" >> daofind.sym
printf "\tcircle\tcyan\t2\t0\thelvetica\t10\tnormal\troman\t\t$fwhm\t\tarcsec" >> daofind.sym

$ds9_command $input_cube -catalog import csv allstars$catalog_out -catalog symbol load daofind.sym -catalog regions -regions save $dao_find_reg -exit


if $run_Zap ; then
	
	$exec_python/cleanskyres.py $input_cube -r $zap_regions_out -i $white_light
	input_cube=zap_cleaned/${input_cube//.fits/zapcleaned.fits}
	
fi


catalog_out=${white_light//fits/csv}

echo "Chose cut extension and create catalog file with source coordinates"

echo "Redshift $redshift"

$ds9_command -tile $chandra_image -scale mode 100 $input_cube -scale mode 99.5 $white_light -scale mode 99.5 -regions load all $ulx_region -regions load all $dao_find_reg & 
camel_catalog=camel_cat.txt
gedit camel_cat.txt &

get_pixelsize="from mpdaf.obj import Image;image_file = '${white_light}'; image = Image(image_file);print('Pixel size %.1f' % image.primary_header['HIERARCH ESO OCS IPS PIXSCALE'])"
	python -c "${get_pixelsize}"

echo 'Input pixel size of the cube'

read fwhm_pixels

gauss_smooth=$(echo "scale=0; $fwhm/$fwhm_pixels/2" | bc)

echo 'Input cut around the source in pixels'
read dxy


echo "Input source ID"
read sourceId
lines=n2has2
mkdir _$sourceId\_$lines

create_config="import os; cwd = os.getcwd(); os.chdir(os.environ['HOME'] + '/scripts/pythonscripts/muse/camel/'); import create_config as cg; os.chdir(cwd); cg.create_config('./','${input_cube}', '${input_cube}', '${camel_catalog}', '${lines}', '', commonw=False, dz=(float('${redshift}')/5), dxy=int('${dxy}'), wmax=1000, degcont=1, ssmooth='${gauss_smooth}')"
python -c "${create_config}"








