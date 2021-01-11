#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

// gcc -o waveplot waveplot.c -lm

int main(int argc, char *argv[])
{
    int i, j, k1, k2, M;
    int Nifo, Tobs, ttrig, fmax;
    int N, NS;
    double dt, x, y, z;
    double t1, t2;
    int flag;
    int *labels;
    double *wave;
    char command[1024];
    
    FILE *out;
    FILE *in;
    
    if(argc<5)
    {
        printf("./waveplot Tobs trig_time fmax detector-list\n");
         return 1;
    }
    
    Nifo = argc-4;
    labels = (int*)malloc(sizeof(int)* (Nifo));
    for (i = 0; i < Nifo; ++i) labels[i] = atoi(argv[4+i]);
    
    
    for (i = 0; i < Nifo; ++i) printf("%d ", labels[i]);
    printf("\n");
    

    Tobs = atoi(argv[1]);
    ttrig = atoi(argv[2]);
    fmax = atoi(argv[3]);
    
    N = 2*fmax*Tobs;
    
    dt = (double)(Tobs)/(double)(N);
    
    wave = (double*)malloc(sizeof(double)* (N));
    
    sprintf(command, "waves/waveenv_%d_%d_%d.dat", Tobs, ttrig, labels[0]);
    //printf("%s\n", command);
    in = fopen(command,"r");
    for (i = 0; i < N; ++i) fscanf(in,"%lf%lf%lf%lf", &x, &x, &wave[i], &x);
    fclose(in);
    
    x = -1.0;
    for (i = 0; i < N; ++i)
    {
     if(wave[i] > x)
     {
         x = wave[i];
         j = i;
     }
    }
    
    printf("Max at %f\n", (double)(j)*dt-Tobs+2.0);
    
    
    y = 0.1*x;
    flag = 0;
    for (i = 0; i < N; ++i)
    {
        if(flag == 0 && wave[i] > y)
        {
            k1 = i;
            flag = 1;
        }
    }
    
    flag = 0;
    for (i = N-1; i > -1; --i)
    {
        if(flag == 0 && wave[i] > y)
        {
            k2 = i;
            flag = 1;
        }
    }
    
    k1 -= (int)(2.0e-2/dt);
    k2 += (int)(2.0e-2/dt);
    
    
    z = 0.0;
    in = fopen("dataw.dat","r");
    for (i = 0; i < N; ++i)
    {
      fscanf(in,"%lf", &x);
       for (j = 0; j < Nifo; ++j)
       {
        fscanf(in,"%lf", &x);
        if(i > k1 && i < k2)
        {
            if(fabs(x) > z) z = fabs(x);
        }
       }
    }
    fclose(in);
    
    
    y = ceil(z);
    if(y < 3.0) y = 3.0;
    
    printf("max data %f %f\n", z, y);
    

    
    t1 = (double)(k1)*dt-Tobs+2.0;
    t2 = (double)(k2)*dt-Tobs+2.0;
    
    printf("Plot between %f and %f\n", t1, t2);
    
    if(Nifo==1)
    {
    out = fopen("wave.gnu","w");
    fprintf(out,"set term png enhanced truecolor crop font 'Helvetica,16' size 1200,800\n");
    fprintf(out,"set output 'waveforms.png'\n");
    fprintf(out,"set xrange [%f:%f]\n", t1, t2);
    fprintf(out,"set xtics\n");
    fprintf(out,"set yrange [%f:%f]\n", -y, y);
    fprintf(out,"set ylabel 'h(t)'\n");
    fprintf(out,"set xlabel 't (s)'\n");
    sprintf(command, "waves/wavestat_%d_%d_%d.dat", Tobs, ttrig, labels[0]);
    fprintf(out,"plot '%s' using 1:2:4 notitle with filledcurves lc 'skyblue' fs transparent solid 0.5, 'dataw.dat' using 1:2 title 'ifo %d' with lines lw 1 lc rgb 'grey', '%s' using 1:3 notitle with lines lc rgb 'blue'\n", command, labels[0], command);
    fclose(out);
    }
    
    if(Nifo==2)
    {
    out = fopen("wave.gnu","w");
    fprintf(out,"set term png enhanced truecolor crop font 'Helvetica,16' size 1200,800\n");
    fprintf(out,"set output 'waveforms.png'\n");
    fprintf(out,"set size 1.0,1.0\n");
    fprintf(out,"set origin 0,0\n");
    fprintf(out,"set multiplot\n");
    fprintf(out,"set size 1.0, 0.47\n");
    fprintf(out,"set origin 0.0, 0.51\n");
    fprintf(out,"set xrange [%f:%f]\n", t1, t2);
    fprintf(out,"unset xtics\n");
    fprintf(out,"set yrange [%f:%f]\n", -y, y);
    fprintf(out,"set ylabel 'h(t)'\n");
    sprintf(command, "waves/wavestat_%d_%d_%d.dat", Tobs, ttrig, labels[0]);
    fprintf(out,"plot '%s' using 1:2:4 notitle with filledcurves lc 'skyblue' fs transparent solid 0.5, 'dataw.dat' using 1:2 title 'ifo %d' with lines lw 1 lc rgb 'grey', '%s' using 1:3 notitle with lines lc rgb 'blue'\n", command, labels[0], command);
    fprintf(out,"set origin 0.0, 0.0\n");
    fprintf(out,"set size 1.0, 0.52\n");
    fprintf(out,"set xtics\n");
    fprintf(out,"set xlabel 't (s)'\n");
    sprintf(command, "waves/wavestat_%d_%d_%d.dat", Tobs, ttrig, labels[1]);
    fprintf(out,"plot '%s' using 1:2:4 notitle with filledcurves lc 'skyblue' fs transparent solid 0.5, 'dataw.dat' using 1:3 title 'ifo %d' with lines lw 1 lc rgb 'grey', '%s' using 1:3 notitle with lines lc rgb 'blue'\n", command, labels[1], command);
    fprintf(out,"unset multiplot\n");
    fclose(out);
    }
    
    if(Nifo==3)
    {
    out = fopen("wave.gnu","w");
    fprintf(out,"set term png enhanced truecolor crop font 'Helvetica,16' size 1200,1200\n");
    fprintf(out,"set output 'waveforms.png'\n");
    fprintf(out,"set size 1.0,1.0\n");
    fprintf(out,"set origin 0,0\n");
    fprintf(out,"set multiplot\n");
    fprintf(out,"set size 1.0, 0.33\n");
    fprintf(out,"set origin 0.0, 0.66\n");
    fprintf(out,"set xrange [%f:%f]\n", t1, t2);
    fprintf(out,"unset xtics\n");
    fprintf(out,"set yrange [%f:%f]\n", -y, y);
    fprintf(out,"set ylabel 'h(t)'\n");
    sprintf(command, "waves/wavestat_%d_%d_%d.dat", Tobs, ttrig, labels[0]);
    fprintf(out,"plot '%s' using 1:2:4 notitle with filledcurves lc 'skyblue' fs transparent solid 0.5, 'dataw.dat' using 1:2 title 'ifo %d' with lines lw 1 lc rgb 'grey', '%s' using 1:3 notitle with lines lc rgb 'blue'\n", command, labels[0], command);
    fprintf(out,"set origin 0.0, 0.345\n");
    sprintf(command, "waves/wavestat_%d_%d_%d.dat", Tobs, ttrig, labels[1]);
    fprintf(out,"plot '%s' using 1:2:4 notitle with filledcurves lc 'skyblue' fs transparent solid 0.5, 'dataw.dat' using 1:3 title 'ifo %d' with lines lw 1 lc rgb 'grey', '%s' using 1:3 notitle with lines lc rgb 'blue'\n", command, labels[1], command);
    fprintf(out,"set origin 0.0, 0.0\n");
    fprintf(out,"set size 1.0, 0.36\n");
    fprintf(out,"set xtics\n");
    fprintf(out,"set xlabel 't (s)'\n");
    sprintf(command, "waves/wavestat_%d_%d_%d.dat", Tobs, ttrig, labels[2]);
    fprintf(out,"plot '%s' using 1:2:4 notitle with filledcurves lc 'skyblue' fs transparent solid 0.5, 'dataw.dat' using 1:4 title 'ifo %d' with lines lw 1 lc rgb 'grey', '%s' using 1:3 notitle with lines lc rgb 'blue'\n", command, labels[2], command);
    fprintf(out,"unset multiplot\n");
    fclose(out);
    }
    
    sprintf(command,"gnuplot wave.gnu");
    system(command);
    
    
    free(wave);
    free(labels);
    
    
    
    return 0;

}
