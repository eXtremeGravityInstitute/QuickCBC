# QuickCBC
Rapid and Reliable Inference for Binary Mergers

The QuickCBC code package implements the analysis pipeline described in arXiv:2101.01188

The package include the main analysis code, scripts to run the code and make plots, and helper codes for extracting data from the GWOSC website. The main analysis code uses openMP https://www.openmp.org. 

The first step in running the code is to download data from GWOSC. To select an event, go to the event list https://www.gw-openscience.org/eventapi/html/GWTC/. For example, to run on GW150914, follow the link to https://www.gw-openscience.org/eventapi/html/GWTC-1-confident/GW150914/v3/ and download the data for H1 and L1 under the tabs 32 sec, 16 KHz, TXT.

Next you need to compile and run the handy code gwoscdump.c

gcc -o gwoscdump gwoscdump.c -lm

This code reads in the data and extracts a snipper for analysis. The arguments are filename, observation time, integer GPS time of the trigger and observatory label. The codes use the labels 0 for H1, 1 for L1 and 2 for Virgo. The code has not yet been set up to run on GEO or Kagra data. The GPS trigger time can be found on the GWOSC site. For GW150914 you would type

./gwoscdump H-H1_GWOSC_16KHZ_R1-1126259447-32.txt 4 1126259462 0
./gwoscdump L-L1_GWOSC_16KHZ_R1-1126259447-32.txt 4 1126259462 1

Here we have chosen to use a 4 second observation time. The gwoscdump code extracts 4 second of data, with the end time of the segment placed 2 seconds after the reference GPS time. This is true no matter how large you set the observation time.

For high mass systems like GW150914, 4 seconds of data is usually enough. For BNS systems, such as GW170817, more data is needed. Good results for GW170817 can be found using 16 seconds of data. If you are prepared to wait a little longer for results then 32 or 64 seconds are better.

The data for GW170817 can be dumped using

./gwoscdump H-H1_GWOSC_16KHZ_R1-1187008867-32.txt 16 1187008882 0
./gwoscdump L-L1_GWOSC_16KHZ_R1-1187008867-32.txt 16 1187008882 1
./gwoscdump V-V1_GWOSC_16KHZ_R1-1187008867-32.txt 16 1187008882 2

To use longer observation times you will need to download the 4096 second long data files from the GWOSC website.

From here the rest of the analysis can be run by calling the script run.sh. This script compiles and runs all the codes, sets up an analysis directory, makes plots etc. The script is hard coded to run on the Mac OS X operating system. You will need to modify the compile commands for a Linux system. Examples of Linux compile lines are included for each of the codes, but they may need to be tweaked depending on your setup.


To run a full analysis on GW150914 you simply type

source run.sh 4 1126259462 0 1

For GW170817 you type

source run.sh 16 1187008882 0 1 2

The arguments are the observation time, integer GPS trigger time, then the list of detectors to be used. The plotting scripts use a combination of gnuplot and python. You will need to have corner.py and healpy installed.

The script produces many plots. The plots named Qdata*png and are time-frequency maps of the data. Qtrack*png superimpose the MAP time-frequency track of the signal. blink*gif are animated gifs comparing time-frequency maps of the data and residual. sky_0.png is the low latency sky map and sky_1.png is the final sky map. waveforms.png show Bayesorams of the whitened waveform reconstructions and original data. masses.png and distance.png are example corner plots of the posteriors. 

The parameters that control the analysis are found in the file ConstCBC.h. The purpose of each constant is listed the file. 



