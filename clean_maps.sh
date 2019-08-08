threshold_OIII=3
threshold_H=3
threshold_S=3
threshold_N=3
threshold_He=3
threshold_OI=3

region=../../ixo25.reg

if [ "$@" -eq 1 ]; then
	region=$1
fi

files_root=*_ssmooth

echo "File root: $files_root"

python ~/scripts/pythonscripts/muse/cleanlinemaps.py -r $files_root -l HALPHA -t $threshold_H
python ~/scripts/pythonscripts/muse/cleanlinemaps.py -r $files_root -l SII6716 -t $threshold_S
python ~/scripts/pythonscripts/muse/cleanlinemaps.py -r $files_root -l SII6731 -t $threshold_S
python ~/scripts/pythonscripts/muse/cleanlinemaps.py -r $files_root -l NII6548 -t $threshold_N
python ~/scripts/pythonscripts/muse/cleanlinemaps.py -r $files_root -l NII6583 -t $threshold_N
python ~/scripts/pythonscripts/muse/cleanlinemaps.py -r $files_root -l HBETA -t $threshold_H
python ~/scripts/pythonscripts/muse/cleanlinemaps.py -r $files_root -l OIII4959 -t $threshold_OIII
python ~/scripts/pythonscripts/muse/cleanlinemaps.py -r $files_root -l OIII5007 -t $threshold_OIII
python ~/scripts/pythonscripts/muse/cleanlinemaps.py -r $files_root -l HEI5876 -t $threshold_OIII
python ~/scripts/pythonscripts/muse/cleanlinemaps.py -r $files_root -l OI6300 -t $threshold_OIII

python ~/scripts/pythonscripts/muse/cleanlinemaps.py -i *vel* -s cleaned_images/clean$files_root\_flux_indep_SII6731.fits
python ~/scripts/pythonscripts/muse/cleanlinemaps.py -i *z* -s cleaned_images/clean$files_root\_flux_indep_SII6731.fits

python ~/scripts/pythonscripts/muse/cleanlinemaps.py -i *vel* -s cleaned_images/clean$files_root\_flux_indep_HALPHA.fits
python ~/scripts/pythonscripts/muse/cleanlinemaps.py -i *z* -s cleaned_images/clean$files_root\_flux_indep_HALPHA.fits

python ~/scripts/pythonscripts/muse/cleanlinemaps.py -i *vel* -s cleaned_images/clean$files_root\_flux_indep_HBETA.fits
python ~/scripts/pythonscripts/muse/cleanlinemaps.py -i *z* -s cleaned_images/clean$files_root\_flux_indep_HBETA.fits

python ~/scripts/pythonscripts/muse/cleanlinemaps.py -i *vel* -s cleaned_images/clean$files_root\_flux_indep_OIII5007.fits
python ~/scripts/pythonscripts/muse/cleanlinemaps.py -i *z* -s cleaned_images/clean$files_root\_flux_indep_OIII5007.fits

python ~/scripts/pythonscripts/muse/cleanlinemaps.py -i *vel* -s cleaned_images/clean$files_root\_flux_indep_HEI5876.fits
python ~/scripts/pythonscripts/muse/cleanlinemaps.py -i *z* -s cleaned_images/clean$files_root\_flux_indep_HEI5876.fits

python ~/scripts/pythonscripts/muse/cleanlinemaps.py -i *vel* -s cleaned_images/clean$files_root\_flux_indep_OI6300.fits
python ~/scripts/pythonscripts/muse/cleanlinemaps.py -i *z* -s cleaned_images/clean$files_root\_flux_indep_OI6300.fits

echo "Createing continuum image..."

input_cube=$(ls $files_root"_contcube.fits")

white_light_image="from mpdaf.obj import Cube;cubein = '${input_cube}'; cube = Cube(cubein);img = cube.sum(axis=0); print('Saving image to white%s' % cubein);  img.write('white%s' % cubein)"
python -c "${white_light_image}"
mv white$files_root"_contcube.fits" cleaned_images/

#cd cleaned_images
#python ~/scripts/pythonscripts/muse/create_images.py * -r $region

