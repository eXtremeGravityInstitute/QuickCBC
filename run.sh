#!/bin/sh

# Needs Tobs trig ifo list


na=$#

nifo=$((na-2))

echo "# Detectors =" $nifo

T=$1
N=$((T*2048))

echo "# Points =" $N

trig=$(echo $2 | awk '{print int($1)}')


foo=$(printf "%d_%d_%d" $trig $1 $nifo)

mkdir $foo

cp SpecFit.c ResFit.c QuickCBC.c waveplot.c waveband.c Utilities.c skydensity.c write_fits.c IMRPhenomD_internals.c IMRPhenomD.c IMRPhenomD_internals.h IMRPhenomD.h Utilities.h QuickCBC.h ConstCBC.h Constants.h chealpix.h fitsio.h Qscan.gnu Qscanwide.gnu Qscanwidehigh.gnu Qtrack.gnu Qtrackwide.gnu MAP_residual.sh residual.sh corner*py sky.sh  makesky.py corner_sky.py $foo

if [ "$nifo" -eq 1 ]; then
goo=$(printf "frame_%d_%d_%d.dat" $1 $trig $3)
cp $goo $foo
fi

if [ "$nifo" -eq 2 ]; then
goo=$(printf "frame_%d_%d_%d.dat" $1 $trig $3)
cp $goo $foo
goo=$(printf "frame_%d_%d_%d.dat" $1 $trig $4)
cp $goo $foo
fi

if [ "$nifo" -eq 3 ]; then
goo=$(printf "frame_%d_%d_%d.dat" $1 $trig $3)
cp $goo $foo
goo=$(printf "frame_%d_%d_%d.dat" $1 $trig $4)
cp $goo $foo
goo=$(printf "frame_%d_%d_%d.dat" $1 $trig $5)
cp $goo $foo
fi

cd $foo

clang -Xpreprocessor -fopenmp -lomp -w -o SpecFit SpecFit.c Utilities.c -lgsl  -lm
clang -Xpreprocessor -fopenmp -lomp -w -o ResFit ResFit.c Utilities.c -lgsl  -lm
gcc -w -o waveband waveband.c -lgsl -lm
gcc -o waveplot waveplot.c -lm
gcc -o write_fits write_fits.c -lm -lcfitsio
gcc -o skydensity skydensity.c -lgsl
clang -Xpreprocessor -fopenmp -lomp -w -o QuickCBC QuickCBC.c Utilities.c IMRPhenomD_internals.c IMRPhenomD.c -lgsl -lgslcblas  -lm

mkdir waves
mkdir specs

if [ "$nifo" -eq 1 ]; then
./SpecFit $1 $trig $3 1024
foo=$(printf "Qdata_%d.dat" $3)
cp Qtransform.dat $foo
gnuplot Qscan.gnu
foo=$(printf "Qdata_%d_%d.png" $trig $3)
mv Qscan.png $foo
if [ "$1" -gt  4 ]; then
gnuplot Qscanwide.gnu
foo=$(printf "Qdatawide_%d_%d.png" $trig $3)
mv Qscanwide.png $foo
gnuplot Qscanwidehigh.gnu
foo=$(printf "Qdatawidehigh_%d_%d.png" $trig $3)
mv Qscanwidehigh.png $foo
fi
fi


if [ "$nifo" -eq 2 ]; then
./SpecFit $1 $trig $3 1024
foo=$(printf "Qdata_%d.dat" $3)
cp Qtransform.dat $foo
gnuplot Qscan.gnu
foo=$(printf "Qdata_%d_%d.png" $trig $3)
mv Qscan.png $foo
if [ "$1" -gt  4 ]; then
gnuplot Qscanwide.gnu
foo=$(printf "Qdatawide_%d_%d.png" $trig $3)
mv Qscanwide.png $foo
gnuplot Qscanwidehigh.gnu
foo=$(printf "Qdatawidehigh_%d_%d.png" $trig $3)
mv Qscanwidehigh.png $foo
fi
./SpecFit $1 $trig $4 1024
foo=$(printf "Qdata_%d.dat" $4)
cp Qtransform.dat $foo
gnuplot Qscan.gnu
foo=$(printf "Qdata_%d_%d.png" $trig $4)
mv Qscan.png $foo
gnuplot Qscanwide.gnu
if [ "$1" -gt  4 ]; then
foo=$(printf "Qdatawide_%d_%d.png" $trig $4)
mv Qscanwide.png $foo
gnuplot Qscanwidehigh.gnu
foo=$(printf "Qdatawidehigh_%d_%d.png" $trig $4)
mv Qscanwidehigh.png $foo
fi
fi

if [ "$nifo" -eq 3 ]; then
./SpecFit $1 $trig $3 1024
foo=$(printf "Qdata_%d.dat" $3)
cp Qtransform.dat $foo
gnuplot Qscan.gnu
foo=$(printf "Qdata_%d_%d.png" $trig $3)
mv Qscan.png $foo
if [ "$1" -gt  4 ]; then
gnuplot Qscanwide.gnu
foo=$(printf "Qdatawide_%d_%d.png" $trig $3)
mv Qscanwide.png $foo
gnuplot Qscanwidehigh.gnu
foo=$(printf "Qdatawidehigh_%d_%d.png" $trig $3)
mv Qscanwidehigh.png $foo
fi
./SpecFit $1 $trig $4 1024
foo=$(printf "Qdata_%d.dat" $4)
cp Qtransform.dat $foo
gnuplot Qscan.gnu
foo=$(printf "Qdata_%d_%d.png" $trig $4)
mv Qscan.png $foo
if [ "$1" -gt  4 ]; then
gnuplot Qscanwide.gnu
foo=$(printf "Qdatawide_%d_%d.png" $trig $4)
mv Qscanwide.png $foo
gnuplot Qscanwidehigh.gnu
foo=$(printf "Qdatawidehigh_%d_%d.png" $trig $4)
mv Qscanwidehigh.png $foo
fi
./SpecFit $1 $trig $5 1024
foo=$(printf "Qdata_%d.dat" $5)
cp Qtransform.dat $foo
gnuplot Qscan.gnu
foo=$(printf "Qdata_%d_%d.png" $trig $5)
mv Qscan.png $foo
if [ "$1" -gt  4 ]; then
gnuplot Qscanwide.gnu
foo=$(printf "Qdatawide_%d_%d.png" $trig $5)
mv Qscanwide.png $foo
gnuplot Qscanwidehigh.gnu
foo=$(printf "Qdatawidehigh_%d_%d.png" $trig $5)
mv Qscanwidehigh.png $foo
fi
fi


echo "Starting MCMC"


if [ "$nifo" -eq 1 ]; then
./QuickCBC $1 $2 $3
foo=$(printf "Qdata_%d.dat" $3)
cp $foo Qtransform.dat
foo=$(printf "tftrack_%d.dat" $3)
cp $foo tftrack.dat
gnuplot Qtrack.gnu
foo=$(printf "Qtrack_%d_%d.png" $trig $3)
mv Qtrack.png $foo
if [ "$1" -gt  4 ]; then
gnuplot Qtrackwide.gnu
foo=$(printf "Qtrackwide_%d_%d.png" $trig $3)
mv Qtrackwide.png $foo
fi
source MAP_residual.sh $1 $trig $3
./waveband $1 $trig $N 100 $3
./waveplot $1 $trig 1024 $3
fi

if [ "$nifo" -eq 2 ]; then
./QuickCBC $1 $2 $3 $4
foo=$(printf "Qdata_%d.dat" $3)
cp $foo Qtransform.dat
foo=$(printf "tftrack_%d.dat" $3)
cp $foo tftrack.dat
gnuplot Qtrack.gnu
foo=$(printf "Qtrack_%d_%d.png" $trig $3)
mv Qtrack.png $foo
if [ "$1" -gt  4 ]; then
gnuplot Qtrackwide.gnu
foo=$(printf "Qtrackwide_%d_%d.png" $trig $3)
mv Qtrackwide.png $foo
fi
foo=$(printf "Qdata_%d.dat" $4)
cp $foo Qtransform.dat
foo=$(printf "tftrack_%d.dat" $4)
cp $foo tftrack.dat
gnuplot Qtrack.gnu
foo=$(printf "Qtrack_%d_%d.png" $trig $4)
mv Qtrack.png $foo
if [ "$1" -gt  4 ]; then
gnuplot Qtrackwide.gnu
foo=$(printf "Qtrackwide_%d_%d.png" $trig $4)
mv Qtrackwide.png $foo
fi
source MAP_residual.sh $1 $trig $3
source MAP_residual.sh $1 $trig $4
./waveband $1 $trig $N 100 $3 $4
./waveplot $1 $trig 1024 $3 $4
fi

if [ "$nifo" -eq 3 ]; then
./QuickCBC $1 $2 $3 $4 $5
foo=$(printf "Qdata_%d.dat" $3)
cp $foo Qtransform.dat
foo=$(printf "tftrack_%d.dat" $3)
cp $foo tftrack.dat
gnuplot Qtrack.gnu
foo=$(printf "Qtrack_%d_%d.png" $trig $3)
mv Qtrack.png $foo
if [ "$1" -gt  4 ]; then
gnuplot Qtrackwide.gnu
foo=$(printf "Qtrackwide_%d_%d.png" $trig $3)
mv Qtrackwide.png $foo
fi
foo=$(printf "Qdata_%d.dat" $4)
cp $foo Qtransform.dat
foo=$(printf "tftrack_%d.dat" $4)
cp $foo tftrack.dat
gnuplot Qtrack.gnu
foo=$(printf "Qtrack_%d_%d.png" $trig $4)
mv Qtrack.png $foo
if [ "$1" -gt  4 ]; then
gnuplot Qtrackwide.gnu
foo=$(printf "Qtrackwide_%d_%d.png" $trig $4)
mv Qtrackwide.png $foo
fi
foo=$(printf "Qdata_%d.dat" $5)
cp $foo Qtransform.dat
foo=$(printf "tftrack_%d.dat" $5)
cp $foo tftrack.dat
gnuplot Qtrack.gnu
foo=$(printf "Qtrack_%d_%d.png" $trig $5)
mv Qtrack.png $foo
if [ "$1" -gt  4 ]; then
gnuplot Qtrackwide.gnu
foo=$(printf "Qtrackwide_%d_%d.png" $trig $5)
mv Qtrackwide.png $foo
fi
source MAP_residual.sh $1 $trig $3
source MAP_residual.sh $1 $trig $4
source MAP_residual.sh $1 $trig $5
./waveband $1 $trig $N 100 $3 $4 $5
./waveplot $1 $trig 1024 $3 $4 $5
fi

source sky.sh skychain.dat 1
python corner_sky.py
python corner_dframe.py
python corner_mass.py
python corner_distance.py


