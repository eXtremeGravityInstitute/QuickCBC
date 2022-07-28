set term pdfcairo font 'Helvetica,14'
set output "specband.pdf" 
set xrange [1:1023]
set logscale y
set yrange [1e-47:1e-42]
set ylabel 'S(f) [Hz^{-1}]'
set xlabel 'f [Hz]'
plot "specstat.dat" using 1:2:4 notitle with filledcurves lc "skyblue","specstat.dat" using 1:3 notitle with lines lw 1 lc rgb "blue"


