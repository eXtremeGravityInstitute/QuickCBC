#!/bin/sh

./skydensity $1 128

touch sky.fits

rm sky.fits

./write_fits sky.dat

python makesky.py

foo=$(printf "%d.fits" $2)

mv sky.fits $foo

foo=$(printf "sky_%d.png" $2)

cp sky.png $foo

open sky.png 

