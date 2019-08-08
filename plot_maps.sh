region=../../ulx1.reg
string=""
for arg in "$@"; do
	string="$arg $string"
done
echo "ds9 command: $string"
$DS9_PATH/ds9 -wcs align yes -cmap b -scale mode 90 $string -region load all $region -tile grid gap 0
