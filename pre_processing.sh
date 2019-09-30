# paths
hst_image="$HOME/optical_data/NGC1672/MAST_2019-06-11T1148/HLA/hst_10354_01_acs_wfc_f550m/hst_10354_01_acs_wfc_f550m_drz.fits" # path to HST image to adjust coordintes
input_cube="ADP.2016-06-20T16_43_08.136.fits" # path to the input cube
chandra_image="/home/agurpide/x_ray_data/NGC1672/chandradata/5932/detectedsources/sourceimage_broadband.fits"
ulx_region="/home/agurpide/x_ray_data/NGC1672/chandradata/5932/repro/ixo27.reg"
run_Zap=false
clean_nan=false
ds9_command="$DS9_PATH/ds9 -cmap Heat -mode region"
exec_python="python3 $HOME/scripts/pythonscripts/muse"
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
	printf "Adjusting coordinates with Hubble image: $hst_image \n"
	$exec_python/adjust_coordinates.py $input_cube --hst $hst_image -o coordadjusted
	input_cube=corr$input_cube
	white_light=wlight$input_cube
	cd coordadjusted
	
else
	# create white light image otherwise
	white_light_image="from mpdaf.obj import Cube;cubein = '${input_cube}'; cube = Cube(cubein);img = cube.sum(axis=0); print('Saving image to white%s' % cubein);  img.write('white%s' % cubein)"
	python3 -c "${white_light_image}"
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
	$exec_python/determine_psf.py $input_cube --fit -r $star_reg -l 100 -t 6
else
# use Tokovinin mode instead
	$exec_python/determine_psf.py $input_cube -t 6
fi


echo "Input average FWHM in arcsec:"

read fwhm

$exec_python/findstars.py $white_light --fwhm $fwhm -t 4 -s 2

echo "Copying catalog regions and saving them as zap_mask.reg"

zap_region_size=$(echo "scale=3; 2*$fwhm" | bc)

zap_regions_out="zap_mask"${white_light//fits/reg}


echo "Converting dao sources to ZAP regions with radius $zap_region_size arcsec"

rm zap_mask.sym

printf "condition\tshape\tcolor\twidth\tdash\tfont\tfontsize\tfontweight\tfontslant\ttext\tsize\tsize2\tunits\tangle\n---------\t-----\t-----\t-----\t----\t----\t--------\t----------\t---------\t----\t----\t-----\t-----\t-----\n" >> zap_mask.sym

printf "\tcircle\tcyan\t2\t0\thelvetica\t10\tnormal\troman\t\t$zap_region_size\t\tarcsec" >> zap_mask.sym

$ds9_command $input_cube -catalog import csv allstars$catalog_out -catalog symbol load zap_mask.sym -catalog regions -regions save $zap_regions_out -exit

dao_find_reg="daofind"${white_light//fits/reg}

rm $dao_find_reg

printf "condition\tshape\tcolor\twidth\tdash\tfont\tfontsize\tfontweight\tfontslant\ttext\tsize\tsize2\tunits\tangle\n---------\t-----\t-----\t-----\t----\t----\t--------\t----------\t---------\t----\t----\t-----\t-----\t-----\n" >> daofind.sym
printf "\tcircle\tcyan\t2\t0\thelvetica\t10\tnormal\troman\t\t$fwhm\t\tarcsec" >> daofind.sym

$ds9_command $white_light -catalog import csv allstars$catalog_out -catalog symbol load daofind.sym -catalog regions -regions save $dao_find_reg -exit


if $run_Zap ; then
	
	$exec_python/cleanskyres.py $input_cube -r $zap_regions_out -i $white_light
	input_cube=zap_cleaned/${input_cube//.fits/zapcleaned.fits}
	
fi


catalog_out=${white_light//fits/csv}

echo "Chose cut extension and create catalog file with source coordinates"

echo "Redshift $redshift"

$ds9_command -tile $chandra_image -scale mode 100 $input_cube -scale mode 99.5 $white_light -scale mode 99.5 -regions load all $ulx_region -regions load all $dao_find_reg & 

# catalog for camel
camel_catalog=camel_cat.txt
gedit camel_cat.txt &

get_pixelsize="from mpdaf.obj import Image;image_file = '${white_light}'; image = Image(image_file);print('Pixel size %.1f' % image.primary_header['HIERARCH ESO OCS IPS PIXSCALE'])"
	python3 -c "${get_pixelsize}"

echo 'Input pixel size of the cube'

read fwhm_pixels

gauss_smooth=$(echo "scale=0; $fwhm/$fwhm_pixels/2" | bc)

echo "Gaussian smoothing in pixels: $gauss_smooth"

echo 'Input cut around the source in pixels'
read dxy

echo "Input source ID"
read sourceId

# Halpha nitrogen complex
lines=n2ha
mkdir "camel_"$sourceId\_$lines

create_config="import os; cwd = os.getcwd(); os.chdir(os.environ['HOME'] + '/scripts/pythonscripts/muse/camel/'); import create_config as cg; os.chdir(cwd); cg.create_config('./','${input_cube}', '${input_cube}', '${camel_catalog}', '${lines}', 'camel', commonw=False, dz=0.005, dxy=int('${dxy}'), wmax=1000, degcont=1, ssmooth='${gauss_smooth}', initw=100)"
python3 -c "${create_config}"

# sulfur lines
lines=s2
mkdir "camel_"$sourceId\_$lines

create_config="import os; cwd = os.getcwd(); os.chdir(os.environ['HOME'] + '/scripts/pythonscripts/muse/camel/'); import create_config as cg; os.chdir(cwd); cg.create_config('./','${input_cube}', '${input_cube}', '${camel_catalog}', '${lines}', 'camel', commonw=False, dz=0.005, dxy=int('${dxy}'), wmax=1000, degcont=1, ssmooth='${gauss_smooth}', initw=100)"
python3 -c "${create_config}"

# hbeta
lines=hb
mkdir "camel_"$sourceId\_$lines

create_config="import os; cwd = os.getcwd(); os.chdir(os.environ['HOME'] + '/scripts/pythonscripts/muse/camel/'); import create_config as cg; os.chdir(cwd); cg.create_config('./','${input_cube}', '${input_cube}', '${camel_catalog}', '${lines}', 'camel', commonw=False, dz=0.005, dxy=int('${dxy}'), wmax=1000, degcont=1, ssmooth='${gauss_smooth}', initw=100)"
python3 -c "${create_config}"

# oxygen III lines 
lines=o3
mkdir "camel_"$sourceId\_$lines

create_config="import os; cwd = os.getcwd(); os.chdir(os.environ['HOME'] + '/scripts/pythonscripts/muse/camel/'); import create_config as cg; os.chdir(cwd); cg.create_config('./','${input_cube}', '${input_cube}', '${camel_catalog}', '${lines}', 'camel', commonw=False, dz=0.005, dxy=int('${dxy}'), wmax=1000, degcont=1, ssmooth='${gauss_smooth}', initw=100)"
python3 -c "${create_config}"

# oxygen I line
lines=oI
mkdir "camel_"$sourceId\_$lines

create_config="import os; cwd = os.getcwd(); os.chdir(os.environ['HOME'] + '/scripts/pythonscripts/muse/camel/'); import create_config as cg; os.chdir(cwd); cg.create_config('./','${input_cube}', '${input_cube}', '${camel_catalog}', '${lines}', 'camel', commonw=False, dz=0.005, dxy=int('${dxy}'), wmax=1000, degcont=1, ssmooth='${gauss_smooth}', initw=100)"
python3 -c "${create_config}"

# He I line
lines=heI
mkdir "camel_"$sourceId\_$lines

create_config="import os; cwd = os.getcwd(); os.chdir(os.environ['HOME'] + '/scripts/pythonscripts/muse/camel/'); import create_config as cg; os.chdir(cwd); cg.create_config('./','${input_cube}', '${input_cube}', '${camel_catalog}', '${lines}', 'camel', commonw=False, dz=0.005, dxy=int('${dxy}'), wmax=1000, degcont=1, ssmooth='${gauss_smooth}', initw=100)"
python3 -c "${create_config}"








