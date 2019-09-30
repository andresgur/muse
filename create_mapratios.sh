# Create map line ratios given the directories for each of the lines
# a region can be added to be plotted

camel_line_dirs=./camel_*_*



outdir=lineratios

region=~/optical_data/NGC1672/cube3_ixo25/nancleancubes/coordadjusted/ixo25.reg

nII_line=NII6583
sII_lines=SII
ha_line=HALPHA
hb_line=HBETA
oIII_line=OIII5007
oI_line=OI6300

sII_line6716=SII6716
sII_line6731=SII6731

hb_path=$camel_line_dirs/cleaned_images/*_flux*$hb_line.fits
oIII_path=$camel_line_dirs/cleaned_images/*_flux*$oIII_line.fits
ha_path=$camel_line_dirs/cleaned_images/*_flux*$ha_line.fits
nII_path=$camel_line_dirs/cleaned_images/*_flux*$nII_line.fits
slines_path=$camel_line_dirs/cleaned_images/*_flux*$sII_lines*.fits
oIline_path=$camel_line_dirs/cleaned_images/*_flux*$oI_line*.fits
sII_line6716_map=$(ls $camel_line_dirs/cleaned_images/*_int*$sII_line6716.fits
sII_line6731_map=$(ls $camel_line_dirs/cleaned_images/*_int*$sII_line6731.fits)

hb_line_map=$(ls $hb_path)
ha_line_map=$(ls $ha_path)
oIII_line_map=$(ls $oIII_path)
nII_line_map=$(ls $nII_path)
sline_map=$(ls $slines_path)
oIline_map=$(ls $oIline_path)


#echo $hb_line_map
#echo $ha_line_map
#echo $oIII_line_map
#echo $nII_line_map
#echo $sline_map


python ~/scripts/pythonscripts/muse/lineratio.py -n $oIII_line_map -d $hb_line_map --outdir $outdir -f $oIII_line\_$hb_line"ratio"

python ~/scripts/pythonscripts/muse/lineratio.py -n $nII_line_map -d $ha_line_map --outdir $outdir -f $nII_line\_$ha_line"ratio"

python ~/scripts/pythonscripts/muse/lineratio.py -n $sline_map -d $ha_line_map --outdir $outdir -f $sII_lines\_$ha_line"ratio"

python ~/scripts/pythonscripts/muse/lineratio.py -n $oIline_map -d $ha_line_map --outdir $outdir -f $oI_line\_$ha_line"ratio"

# electon density map
python ~/scripts/pythonscripts/muse/lineratio.py -n $sII_line6716_map -d $sII_line6731_map --outdir $outdir -f electron_density


#cd $outdir
#python ~/scripts/pythonscripts/muse/create_images.py * -r $region
