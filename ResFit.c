/**************************************************************************
 
 Copyright (c) 2019 Neil Cornish
 
 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 
 ************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "Utilities.h"
#include "ConstCBC.h"
#include "Constants.h"

#ifndef _OPENMP
#define omp ignore
#endif

//OSX
// clang -Xpreprocessor -fopenmp -lomp -w -o Res Res.c Utilities.c -lgsl  -lm

// CIT
// gcc -std=gnu99 -fopenmp -w -o Res Res.c Utilities.c -lgsl -lm


int main(int argc, char *argv[])
{
    
    
    int i, ND, N, Ns, m, dec;
    double *timeF, *dataF;
    double *times, *data, *signal, *residual;
    double GMST;
    double x, y, z, dt, Tobs, ttrig, fny, fmn, fmx;
    double *SN, *SM;
    int Nsp, Nl, iter;
    char command[1024];
    
    // allocate some arrays
    double **tfDR, **tfDI;
    double **tfD;
    int j, subscale, octaves, Nf;
    double *freqs, *sqf;
    double dx, dlnf, t, f, Q, alpha, fac, t_rise;
    
    
    
    FILE *in;
    FILE *out;
    
    if(argc!=5)
    {
        printf("./Resf Tobs iter trig_time ifo\n");
        return 1;
    }
    
    
    Tobs = atof(argv[1]);
    iter = atoi(argv[2]);
    ttrig = atof(argv[3]);
    m = atoi(argv[4]);
    
    
    fmn = 8.0;
    
    sprintf(command, "framed_%d_%d_%d.dat", (int)(Tobs), (int)ttrig, m);
    
    in = fopen(command,"r");
    
    
    N = -1;
    while(!feof(in))
    {
        fscanf(in,"%lf%lf", &x, &y);
        N++;
    }
    rewind(in);
    printf("Number of points = %d\n", N);
    
    times = (double*)malloc(sizeof(double)* (N));
    data = (double*)malloc(sizeof(double)* (N));
    residual = (double*)malloc(sizeof(double)* (N));
    signal = (double*)malloc(sizeof(double)* (N));
    
    for (i = 0; i < N; ++i)
    {
        fscanf(in,"%lf%lf", &times[i], &data[i]);
    }
    fclose(in);
    
    dt = times[1]-times[0];
    Tobs = (double)(N)*dt;  // duration
    fny = 1.0/(2.0*dt);  // Nyquist
    
    fmx = fny;
    
    
    sprintf(command, "waves/wavef_%d_%d_%d_%d.dat", iter, (int)(Tobs), (int)ttrig, m);
    in = fopen(command,"r");
    
    for (i = 1; i < N/2; ++i)
    {
        fscanf(in,"%lf%lf%lf", &x, &y, &z);
        signal[i] = y;
        signal[N-i] = z;
    }
    fclose(in);
    
    signal[0] = 0.0;
    signal[N/2] = 0.0;
    
    fac = sqrt(Tobs)/(double)(N);
    
    for (i = 0; i < N; ++i) signal[i] /= fac;
    
    SM = (double*)malloc(sizeof(double)*(N/2));
    SN = (double*)malloc(sizeof(double)*(N/2));
    
    printf("Data volume %f seconds, %f Hz\n", Tobs, fny);
    
    fac = Tobs/((double)(N)*(double)(N));
    
    sprintf(command, "PSD_%d_%d_%d.dat", (int)(Tobs), (int)ttrig, m);
    in = fopen(command,"r");
    for (i = 0; i < N/2; ++i)
    {
        fscanf(in,"%lf%lf\n", &x, &SM[i]);
        SM[i] /= fac;
    }
    fclose(in);
    
    // Tukey window parameter. Flat for (1-alpha) of data
    t_rise = 0.4; // Standard LAL setting
    alpha = (2.0*t_rise/Tobs);
    
    
    // Tukey window
    tukey(data, alpha, N);
    // FFT
    gsl_fft_real_radix2_transform(data, 1, N);
    
    // remove signal from data
    for (i = 0; i < N; ++i) residual[i] = data[i] - signal[i];
    
    // whiten data, signal and residual
    whiten(data, SM, N);
    whiten(signal, SM, N);
    whiten(residual, SM, N);
    
    Q = Qs;
    
    // logarithmic frequency spacing
    subscale = 40;  // number of semi-tones per octave
    octaves = (int)(rint(log(fmx/fmn)/log(2.0))); // number of octaves
    Nf = subscale*octaves+1;
    freqs = (double*)malloc(sizeof(double)* (Nf));   // frequencies used in the analysis
    sqf = (double*)malloc(sizeof(double)* (Nf));
    dx = log(2.0)/(double)(subscale);
    dlnf = dx;
    x = log(fmn);
    for(i=0; i< Nf; i++)
    {
        freqs[i] = exp(x);
        sqf[i] = sqrt(freqs[i]);
        x += dx;
    }
    
    printf("%d freqs %d times\n", Nf, N);
    
    tfDR = double_matrix(Nf,N);
    tfDI = double_matrix(Nf,N);
    tfD = double_matrix(Nf,N);
    
    
    // Wavelet transform
    TransformC(data, freqs, tfD, tfDR, tfDI, Qprint, Tobs, N, Nf);
    
    
    out = fopen("Qdata.dat","w");
    for(j = 0; j < Nf; j++)
    {
        f = freqs[j];
        
        for(i = 0; i < N; i++)
        {
            t = (double)(i)*dt-Tobs+2.0;
            fprintf(out,"%e %e %e\n", t, f, tfD[j][i]);
        }
        
        fprintf(out,"\n");
    }
    fclose(out);
    
    // Wavelet transform
       TransformC(residual, freqs, tfD, tfDR, tfDI, Qprint, Tobs, N, Nf);
       
       
       out = fopen("Qresidual.dat","w");
       for(j = 0; j < Nf; j++)
       {
           f = freqs[j];
           
           for(i = 0; i < N; i++)
           {
               t = (double)(i)*dt-Tobs+2.0;
               fprintf(out,"%e %e %e\n", t, f, tfD[j][i]);
           }
           
           fprintf(out,"\n");
       }
       fclose(out);
    
    
    // Wavelet transform
    TransformC(signal, freqs, tfD, tfDR, tfDI, Qprint, Tobs, N, Nf);
    
    
    out = fopen("Qsignal.dat","w");
    for(j = 0; j < Nf; j++)
    {
        f = freqs[j];
        
        for(i = 0; i < N; i++)
        {
            t = (double)(i)*dt-Tobs+2.0;
            fprintf(out,"%e %e %e\n", t, f, tfD[j][i]);
        }
        
        fprintf(out,"\n");
    }
    fclose(out);
    
    
    
    
    free(times);
    free(data);
    free(signal);
    free(residual);
    
    return 0;
    
}



