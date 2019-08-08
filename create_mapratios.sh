
dir1=./camel_*_n2has2
dir2=./camel_*_hbo3

printf "Dir 1: $dir1 \n" 
printf "Dir 2: $dir2 \n"


outdir=lineratios

region=~/optical_data/NGC1672/cube3_ixo25/nancleancubes/coordadjusted/ixo25.reg

nII_line=NII6583
sII_lines=SII
ha_line=HALPHA
hb_line=HBETA
oIII_line=OIII5007

hb_path=$dir2/cleaned_images/*_flux*$hb_line.fits
oIII_path=$dir2/cleaned_images/*_flux*$oIII_line.fits
ha_path=$dir1/cleaned_images/*_flux*$ha_line.fits
nII_path=$dir1/cleaned_images/*_flux*$nII_line.fits
slines_path=$dir1/cleaned_images/*_flux*$sII_lines*.fits

hb_line_map=$(ls $hb_path)
ha_line_map=$(ls $ha_path)
oIII_line_map=$(ls $oIII_path)
nII_line_map=$(ls $nII_path)
sline_map=$(ls $slines_path)

#echo $hb_line_map
#echo $ha_line_map
#echo $oIII_line_map
#echo $nII_line_map
#echo $sline_map


python ~/scripts/pythonscripts/muse/lineratio.py -n $oIII_line_map -d $hb_line_map --outdir $outdir -f $oIII_line\_$hb_line"ratio"

python ~/scripts/pythonscripts/muse/lineratio.py -n $nII_line_map -d $ha_line_map --outdir $outdir -f $nII_line\_$ha_line"ratio"

python ~/scripts/pythonscripts/muse/lineratio.py -n $sline_map -d $ha_line_map --outdir $outdir -f $sII_lines\_$ha_line"ratio"

#cd $outdir
#python ~/scripts/pythonscripts/muse/create_images.py * -r $region
