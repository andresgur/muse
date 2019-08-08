hbo3dir=./camel_*_hbo3
n2hadir=./camel_*_n2has2
heIoI=''
lineratios=lineratios

# regions
# scales
scale_region=~/optical_data/line_scale.reg 
scale_region_hst=~/optical_data/line_scale_hst.reg 
# labels
flux_maps_label=~/optical_data/flux_maps.reg
disp_maps_label=~/optical_data/disp_maps.reg
cont_maps_label=~/optical_data/cont_maps.reg
vel_maps_label=~/optical_data/vel_maps.reg

# lineratios
n2_halpha=~/optical_data/n2_halpha.reg
o3_hbeta=~/optical_data/o3_hbeta.reg
s2_halpha=~/optical_data/s2_halpha.reg

# lines
# fluxes
flux_halpha=~/optical_data/halpha_text.reg
flux_s2=~/optical_data/s2_text.reg
flux_hbeta=~/optical_data/hbeta_text.reg
flux_o3=~/optical_data/o3_text.reg
flux_n2=~/optical_data/n2_text.reg
# dispersion
disp_halpha=~/optical_data/disp_halpha.reg
disp_hbeta=~/optical_data/disp_hbeta.reg
disp_s2=~/optical_data/disp_s2.reg
disp_n2=~/optical_data/disp_n2.reg
disp_o3=~/optical_data/disp_o3.reg

# velocity map
velocity_shift=~/optical_data/velocity_shift.reg

region='ixo25_fk5.reg'

dir=$n2hadir

fluxn2hamaps="$dir/cleaned_images/*[^e]flux*HALPHA* -region $scale_region -region $flux_halpha $dir/cleaned_images/*[^e]flux*NII6583* -region $flux_n2 $dir/cleaned_images/*[^e]flux*SII* -region $flux_s2"

dispn2hamaps="$dir/cleaned_images/*[^e]disp*HALPHA* -region $disp_halpha $dir/cleaned_images/*[^e]disp*NII6583* -region $disp_n2 $dir/cleaned_images/*[^e]disp*SII* -region $disp_s2"

velandcontha="$dir/cleaned_images/*[^e]vel* -region $vel_maps_label $dir/cleaned_images/white* -region $cont_maps_label"

dir=$hbo3dir

fluxhbo3maps="$dir/cleaned_images/*[^e]flux*HBETA* -region $flux_hbeta $dir/cleaned_images/*[^e]flux*OIII5007* -region $flux_o3" 

dispnhbo3maps="$dir/cleaned_images/*[^e]disp*HBETA* -region $disp_hbeta $dir/cleaned_images/*[^e]disp*OIII5007* -region $disp_o3"

velandconthb="$dir/cleaned_images/*[^e]vel* -region $vel_maps_label $dir/cleaned_images/white* -region $cont_maps_label"

hst_image=hst_image.fits

$HOME/ds9 -cmap b -smooth yes -smooth function gaussian -smooth sigma 1 -smooth radius 2 -scale mode 99 $fluxn2hamaps $dispn2hamaps $velandcontha $fluxhbo3maps $dispnhbo3maps $velandconthb $lineratios/NII*.fits -region $n2_halpha $lineratios/O*.fits -region $o3_hbeta $lineratios/S*.fits -region $s2_halpha $hst_image -region $scale_region_hst -region load all $region -wcs align yes -tile grid gap 0

# for heI oI
#ds9 -cmap b -smooth yes -smooth function gaussian -smooth sigma 1 -smooth radius 2 -scale mode 99 *[^e]flux* *[^e]disp* *[^e]vel* white* ../../hst_image.fits -region ../../ixo27_fk5.reg	 -wcs align yes -tile grid gap 0
g
