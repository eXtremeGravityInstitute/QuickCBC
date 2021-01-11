#!/bin/sh

# Needs Tobs ttrig ifo

./ResFit $1 -1 $2 $3
cp Qdata.dat Qtransform.dat
gnuplot Qscan.gnu
coo=$(printf "Qwhite_%d_%d.png" $2 $3)
mv Qscan.png $coo
cp Qresidual.dat Qtransform.dat
gnuplot Qscan.gnu
goo=$(printf "MAP_res_%d_%d.png" $2 $3)
mv Qscan.png $goo

woo=$(printf "blink_%d_%d.gif" $2 $3)

convert -delay 100 $coo $goo -loop 0 $woo



