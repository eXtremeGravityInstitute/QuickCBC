#!/bin/sh

# Needs Tobs ttrig ifo

coo=$(printf "Qdata_%d_%d.png" $2 $3)

for ((j=1;j<6;j++))
do
k=$((j*10))
echo $j $k
./Res $1 $k $2 $3
gnuplot Qscan.gnu
goo=$(printf "Res_%d_%d_%d.png" $2 $3 $j)
mv Qscan.png $goo
printf -v line ' %-15s \n' "$goo"
coo+=$line
done

woo=$(printf "%d_%d.gif" $2 $3)

convert -delay 100 $coo -loop 0 $woo



