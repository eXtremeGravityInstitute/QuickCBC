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
#include "Constants.h"

#ifndef _OPENMP
#define omp ignore
#endif

static const gsl_rng_type *rngtype;
static const gsl_rng *rng;

#define dfmin 4.0  // minimum spline spacing
#define dfmax 32.0  // maximum spline spacing
#define smooth 8.0   // moving average
#define tol 0.2     // tolerance for difference in averages
#define mulmax 10.0   // max of smoothness hyper prior
#define mulmin 0.0001     // min of smoothness hyper prior
#define MM 200000    // iterations of MCMC
#define ptype 0     // 0 for Gaussian, 1 for exponential
#define pfix 0      // set to 1 to turn off hyper prior scaling

//OSX
// clang -Xpreprocessor -fopenmp -lomp -w -o SpecFit SpecFit.c Utilities.c -lgsl  -lm

// Linux
// gcc -std=gnu99 -fopenmp -w -o SpecFit SpecFit.c Utilities.c -lgsl -lm



int main(int argc, char *argv[])
{


  int i, j, k, kk, ND, N, Ns, Nm, m, mc, dec, MK, MS, skip;
  int imin, imax, acS, acL, cS, cL;
  int typ, sm;
  double *timeF, *dataF;
  double *times, *data, *Hf;
  double x, y, z, dt, df, Tobs, ttrig, fny, fmn, fmx;
  double f, finc;
  double *SN, *SM, *SL, *PS, *S1, *S2, *Smean, **SR, *sp;
  int Nsp, Nl;
  char command[1024];
  double spread, *deltafmax, *linew;
  double *linef, *lineh, *lineQ;
  double f2, f4, ff4, ff2, df2, df4;
  double ylinew, ylinef, ylineh, ylineQ, ydeltafmax;
  int ii, flag, Nlines;
  double xold, xnext, Abar;
  double max;
  double pmx, pmy;
  int Nspline;
  double a1, a2, sdd;
    
  double alpha;
  double tuke = 0.4;    // Tukey window rise (s)
    
   const gsl_rng_type * T;
   gsl_rng * r;
    
   gsl_rng_env_setup();
    
   T = gsl_rng_default;
   r = gsl_rng_alloc (T);
    
  FILE *in;
  FILE *out;
  FILE *spec;
  FILE *lfile;
  FILE *sfile;
    
    if(argc!=5)
    {
        printf("./Spec Tobs trig_time obs fmax\n");
        return 1;
    }
    
    Tobs = atof(argv[1]);
    ttrig = atof(argv[2]);
    m = atoi(argv[3]);
    fmx = atof(argv[4]);
    
    df = 1.0/Tobs;
    
    MS = 1000;

    
    sprintf(command, "frame_%d_%d_%d.dat", (int)(Tobs), (int)ttrig, m);
    
    in = fopen(command,"r");
 
    
    ND = -1;
    while(!feof(in))
    {
        fscanf(in,"%lf%lf", &x, &y);
        ND++;
    }
    rewind(in);
    printf("Number of points = %d\n", ND);
    
    timeF = (double*)malloc(sizeof(double)* (ND));
    dataF = (double*)malloc(sizeof(double)* (ND));
    
    for (i = 0; i < ND; ++i)
    {
      fscanf(in,"%lf%lf", &timeF[i], &dataF[i]);
    }
    fclose(in);
    
    dt = timeF[1]-timeF[0];
    Tobs = (double)(ND)*dt;  // duration
    fny = 1.0/(2.0*dt);  // Nyquist
    
    printf("Nyquist %f  fmax %f\n", fny, fmx);
    
    // if fmax < fny we can downsample the data by decimation
    // first we bandpass then decimate
    
    dec = (int)(fny/fmx);
    
    if(dec > 8) dec = 8;
    
    
    printf("Down sample = %d\n", dec);

    fmn = 8.0;
    
    // apply 8th order zero phase bandpass filter
    bwbpf(dataF, dataF, 1, ND, 8, 1.0/dt, fmx, fmn);
    bwbpf(dataF, dataF, -1, ND, 8, 1.0/dt, fmx, fmn);
    
    N = ND/dec;
    
    times = (double*)malloc(sizeof(double)* (N));
    data = (double*)malloc(sizeof(double)* (N));
    Hf = (double*)malloc(sizeof(double)* (N));
    
    // decimate
    for (i = 0; i < N; ++i)
    {
        times[i] = timeF[i*dec];
        data[i] = dataF[i*dec];
        Hf[i] = 0.0;
    }
    
    free(timeF);
    free(dataF);

    
    // reset sample rate and Nyquist
    dt = times[1]-times[0];
    fny = 1.0/(2.0*dt);
    
    sprintf(command, "framed_%d_%d_%d.dat", (int)(Tobs), (int)ttrig, m);
    
    out = fopen(command,"w");
    for (i = 0; i < N; ++i) fprintf(out,"%.16e %.16e\n", times[i], data[i]);
    fclose(out);
    
    Ns = (int)(Tobs*fmx);
    
    S1 = (double*)malloc(sizeof(double)*(Ns));
    S2 = (double*)malloc(sizeof(double)*(Ns));
    SL = (double*)malloc(sizeof(double)*(Ns));
    SM = (double*)malloc(sizeof(double)*(Ns));
    SN = (double*)malloc(sizeof(double)*(Ns));
    Smean = (double*)malloc(sizeof(double)*(Ns));
    SR = double_matrix(MS,Ns);
    PS = (double*)malloc(sizeof(double)*(Ns));
    
    printf("Data volume %f seconds, %f Hz\n", Tobs, fny);
    
    specest(data, Hf, N, Ns, dt, fmx, SN, SM, PS);
    
    
    // moving average
    sm = (int)(smooth*Tobs);
    x = 0.0;
    for (i = 0; i < sm; ++i) x += SM[i];
    for (i = sm; i < Ns; ++i)
    {
        S1[i-sm/2] = x/(double)(sm);
        x += SM[i] - SM[i-sm];
    }
    
    // moving average with wider window
    sm *= 2;
    x = 0.0;
    for (i = 0; i < sm; ++i) x += SM[i];
    for (i = sm; i < Ns; ++i)
    {
        S2[i-sm/2] = x/(double)(sm);
        x += SM[i] - SM[i-sm];
    }
    
    // fill initial bins
    for (i = 0; i < sm/2; ++i)
    {
        S1[i] = SM[i];
        S2[i] = 2.0*S1[i];
    }

    /*
    out = fopen("smooth.dat","w");
    for (i = sm/2; i < Ns; ++i)
    {
        f = (double)(i)/Tobs;
        fprintf(out,"%.15e %.15e %.15e\n", f, S1[i], S2[i]);
    }
    fclose(out);
    */
    
   // count the number of spline knots
    Nspline = 1;
    k = (int)(dfmin*Tobs);
    kk = (int)(dfmax*Tobs);
    j = 0;
    flag = 0;
    max = 0.0;
    for (i = 1; i < Ns; ++i)
    {
        x = fabs(S2[i]/S1[i]-1.0);
        if(x > max) max = x;
        j++;
         if(i%k == 0)
          {
              if(max > tol || j == kk)
              {
                 // printf("%f %f %f\n", (double)(i-j/2)/Tobs, (double)(j)/Tobs, max);
                  max = 0.0;
                  j = 0;
                  Nspline++;
              }
          }
        
        
    }
    
     Nspline++;
    
     printf("There are %d spline knots\n", Nspline);
    
     double *ffit, *Xspline, *Yspline;
     
     ffit = (double*)malloc(sizeof(double)*(Nspline));
     Xspline = (double*)malloc(sizeof(double)*(Nspline));
     Yspline = (double*)malloc(sizeof(double)*(Nspline));
    
    Nspline = 1;
    j = 0;
    flag = 0;
    max = 0.0;
    for (i = 1; i < Ns; ++i)
    {
        x = fabs(S2[i]/S1[i]-1.0);
        if(x > max) max = x;
        j++;
         if(i%k == 0)
          {
              if(max > tol || j == kk)
              {
                  max = 0.0;
                  ffit[Nspline] = (double)(i-j/2)/Tobs;
                  Xspline[Nspline] = log(S1[i-j/2]);
                  j = 0;
                  Nspline++;
              }
          }
        
        
    }
    
    Nspline++;
    
    ffit[0] = 0.0;
    Xspline[0] = log(SM[0]);
    ffit[Nspline-1] = (double)(Ns-1)/Tobs;
    Xspline[Nspline-1] = log(SM[Ns-1]);
    
    
     for (i = 0; i < Nspline; ++i)
        {
            Yspline[i] = Xspline[i];
        }
    
    // Allocate spline
    gsl_spline   *cspline = gsl_spline_alloc(gsl_interp_cspline, Nspline);
    gsl_interp_accel *acc    = gsl_interp_accel_alloc();
    
    /* compute spline */
    gsl_spline_init(cspline,ffit,Xspline,Nspline);
    
    /*
    out = fopen("control.dat","w");
       for (i = 0; i < Nspline; ++i)
       {
           fprintf(out,"%e %e\n", ffit[i], exp(Xspline[i]));
       }
       fclose(out);
     */
     
      // count the number of lines
     j = 0;
     flag = 0;
     for (i = 0; i < Ns; ++i)
     {
         x = PS[i]/SM[i];
         // start of a line
         if(x > 9.0 && flag == 0)
         {
             k = 1;
             flag = 1;
             max = x;
             ii = i;
         }
         // in a line
         if(x > 9.0  && flag ==1)
         {
             k++;
             if(x > max)
             {
                 max = x;
                 ii = i;
             }
         }
         // have reached the end of a line
         if(flag == 1)
         {
             if(x < 9.0)
             {
                 flag = 0;
                 j++;
             }
         }
     }
    
     
     Nlines = j;
    


       linef = (double*)malloc(sizeof(double)*(Nlines));  // central frequency
       lineh = (double*)malloc(sizeof(double)*(Nlines));  // line height
       lineQ = (double*)malloc(sizeof(double)*(Nlines)); // line Q
       linew = (double*)malloc(sizeof(double)*(Nlines));  // line width
       deltafmax = (double*)malloc(sizeof(double)*(Nlines));  // cut-off
    
    j = -1;
    xold = 1.0;
    flag = 0;
    for (i = 0; i < Ns; ++i)
    {
        x = PS[i]/SM[i];
        // start of a line
        if(x > 9.0 && flag == 0)
        {
            k = 1;
            flag = 1;
            max = x;
            ii = i;
        }
        // in a line
        if((x > 9.0) && flag ==1)
        {
            k++;
            if(x > max)
            {
                max = x;
                ii = i;
            }
        }
        // have reached the end of a line
        if(flag == 1)
        {
            if(x < 9.0)
            {
                flag = 0;
                j++;
                linef[j] = (double)(ii)/Tobs;
                lineh[j] = max;
                Abar = 0.5*(SN[ii+1]/SM[ii+1]+SN[ii-1]/SM[ii-1]);
                //lineQ[j] = sqrt((max/Abar-1.0))*linef[j]*Tobs;
                 lineQ[j] = sqrt(max)*linef[j]*Tobs/(double)(k);
                
                  spread = (1.0e-2*lineQ[j]);
                  if(spread < 50.0) spread = 50.0;  // maximum half-width is f_resonance/50
                  deltafmax[j] = linef[j]/spread;
                  linew[j] = 8.0*deltafmax[j];
                
                //printf("%d %e %e %e %e %e\n", j, linef[j], lineh[j], linew, (double)(k)/Tobs, lineQ[j]);
               
            }
        }

    }
    
    //printf("There are %d lines\n", Nlines);
    
    // initialize smooth spline spectrum
    SM[0] = exp(Xspline[0]);
    SM[Ns-1] = exp(Xspline[Nspline-1]);
    for (i = 1; i < Ns-1; ++i)
    {
        f = (double)(i)/Tobs;
        SM[i] = exp(gsl_spline_eval(cspline,f,acc));
    }

    
    // initialize line spectrum
    //out = fopen("lines.dat","w");
    for (i = 0; i < Ns; ++i)
    {
        f = (double)(i)/Tobs;
        y = 0.0;
        for (j = 0; j < Nlines; ++j)
        {
           
            x = fabs(f - linef[j]);
            
            if(x < linew[j])
            {
                z = 1.0;
                if(x > deltafmax[j]) z = exp(-(x-deltafmax[j])/deltafmax[j]);
                f2 = linef[j]*linef[j];
                ff2 = f*f;
                df2 = lineQ[j]*(1.0-ff2/f2);
                df4 = df2*df2;
                y += z*lineh[j]/(ff2/f2+df4);
               // printf("%d %d %e %e %e\n", i, j, f2, ff2, df4);
            }
        }
        SL[i] = y;
        //fprintf(out,"%e %e %e\n", f, SM[i]*(1.0+SL[i]), SL[i]);
        
    }
   // fclose(out);
    
    printf("Initially, there are %d lines\n", Nlines);
    
    // changing a control point value impacts the spline in two segments
    // either side of the point changed.
    
    double *freqs, *lnLR, *lnpR;
    
    freqs = (double*)malloc(sizeof(double)*(Ns));
    lnLR = (double*)malloc(sizeof(double)*(Ns));
    lnpR = (double*)malloc(sizeof(double)*(Ns));
    
    /*
       out = fopen("specstart.dat","w");
       for (i = 0; i < Ns; ++i)
       {
           f = (double)(i)/Tobs;
           fprintf(out,"%e %e %e %e\n", f, SM[i], SL[i], SM[i]*(1.0+SL[i]));
       }
       fclose(out);
    */
   
    for (i = 0; i < Ns; ++i)
     {
         freqs[i] = (double)(i)/Tobs;
         x = SM[i]*(1.0+SL[i]);
         lnLR[i] = -(log(x) + PS[i]/x);
         lnpR[i] = 0.0;
     }
    
     for (i = 1; i < Ns-1; ++i)
        {
            sdd = gsl_spline_eval_deriv2(cspline, freqs[i], acc);
            if(ptype == 0) lnpR[i] = -sdd*sdd;
            if(ptype == 1) lnpR[i] = -fabs(sdd);
        }
    
     double *DlnLR, *DlnpR, *DS;
     
    // This holds the updated values in the region impacted by the change in the
    // spline point. Allocaing these to the size of the full spectrum
     DlnLR = (double*)malloc(sizeof(double)*(Ns));
     DlnpR = (double*)malloc(sizeof(double)*(Ns));
     DS = (double*)malloc(sizeof(double)*(Ns));
    
    double logLx, logLy, logpx, logpy;
    double H;
    
      logLx = 0.0;
      for (i = 0; i < Ns; ++i) logLx += lnLR[i];
    
      logpx = 0.0;
      for (i = 1; i < Ns-1; ++i) logpx += lnpR[i];
    
    sprintf(command, "ffile_%d.dat", m);
    sfile = fopen(command,"w");
    for (k = 0; k < Nspline; ++k) fprintf(sfile, "%e ", ffit[k]);
    fprintf(sfile, "\n");
    fclose(sfile);
    
    sprintf(command, "sfile_%d.dat", m);
    sfile = fopen(command,"w");
    
    sprintf(command, "lfile_%d.dat", m);
    lfile = fopen(command,"w");
    
    acS = 0;
    acL = 0;
    cS = 1;
    cL = 1;
    
    if(ptype == 0) pmx = 0.03;
    if(ptype == 1) pmx = 0.005;
    pmy = pmx;
    
    out = fopen("schain.dat","w");
    
    a1 = 0.5;
    a2 = 1.0;
    
    if(pfix == 0) //turn on dynamic scaling of smoothness prior
     {
      a1 = 0.2;
      a2 = 0.4;
     }
    
    ii = 0;
    
    skip = 500;
    
    MK = skip*MS;  //
    
     for (mc = 0; mc < MM+MK; ++mc)
     {
         
         flag = 0;
         pmy = pmx;
         logLy = logLx;
         logpy = logpx;
         
         if(mc >= MM && mc%skip == 0)
         {
           for (i = 0; i < Ns; ++i)
            {
               SR[ii][i] = SM[i]*(1.0+SL[i]);
            }
             ii++;
             
             for (k = 0; k < Nspline; ++k) fprintf(sfile, "%e ", Xspline[k]);
             fprintf(sfile, "\n");
             
             for (k = 0; k < Nlines; ++k) fprintf(lfile, "%e %e %e %e %e ", linef[k], lineh[k], lineQ[k], linew[k], deltafmax[k]);
             fprintf(lfile, "\n");
             
          }
         
         // prune any weak lines
         if(mc == MM/2)
         {
             j = 0;
             for (k = 0; k < Nlines; ++k)
             {
                 i = (int)(linef[k]*Tobs);  // bin where line peaks
                 x = lineh[k];     // height of line relative to smooth
                 if(x > 5.0)  // keep this line
                 {
                     linef[j] = linef[k];
                     lineh[j] = lineh[k];
                     lineQ[j] = lineQ[k];
                     linew[j] = linew[k];
                     deltafmax[j] = deltafmax[k];
                     j++;
                 }
             }
             
             Nlines = j;
             
             // reset line spectrum
             for (i = 0; i < Ns; ++i)
             {
                 f = (double)(i)/Tobs;
                 y = 0.0;
                 for (j = 0; j < Nlines; ++j)
                 {
                    
                     x = fabs(f - linef[j]);
                     
                     if(x < linew[j])
                     {
                          z = 1.0;
                         if(x > deltafmax[j]) z = exp(-(x-deltafmax[j])/deltafmax[j]);
                         f2 = linef[j]*linef[j];
                         ff2 = f*f;
                         df2 = lineQ[j]*(1.0-ff2/f2);
                         df4 = df2*df2;
                         y += z*lineh[j]/(ff2/f2+df4);
                     }
                 }
                 SL[i] = y;
             }
             
             // reset likelihood
             for (i = 0; i < Ns; ++i)
             {
                 freqs[i] = (double)(i)/Tobs;
                 x = SM[i]*(1.0+SL[i]);
                 lnLR[i] = -(log(x) + PS[i]/x);
             }
          }
         
         
         alpha = gsl_rng_uniform(r);
         
         if(alpha < a1) // update spline
         {
             
          typ = 0;
          cS++;
             
        // pick a knot to update
        k = (int)((double)(Nspline)*gsl_rng_uniform(r));
         
        alpha = gsl_rng_uniform(r);
        if(alpha > 0.5)
        {
         Yspline[k] =  Xspline[k] + gsl_ran_gaussian(r,0.05);
        }
        else if (alpha > 0.3)
        {
         Yspline[k] =  Xspline[k] + gsl_ran_gaussian(r,0.02);
        }
        else
        {
         Yspline[k] =  Xspline[k] + gsl_ran_gaussian(r,0.1);
        }
             
        gsl_spline_init(cspline,ffit,Yspline,Nspline);
         
         if(k > 1)
         {
             imin = (int)(ffit[k-2]*Tobs);
         }
         else
         {
             imin = (int)(ffit[0]*Tobs);
         }
         
         if(k < Nspline-2)
         {
             imax = (int)(ffit[k+2]*Tobs);
         }
         else
         {
             imax = (int)(ffit[Nspline-1]*Tobs);
         }
             
            // Delta the smooth spectrum and prior
             logpy = logpx;
             logLy = logLx;
            for (i = imin; i < imax; ++i)
               {
                          y = 0.0;
                           if(i > 0 && i < Ns-1)
                           {
                               x = gsl_spline_eval(cspline,freqs[i],acc);
                               sdd = gsl_spline_eval_deriv2(cspline, freqs[i], acc);
                               if(ptype == 0) y = -sdd*sdd;
                               if(ptype == 1) y = -fabs(sdd);
                           }
                           else
                           {
                               if(i==0) x = Yspline[0];
                               if(i==Ns-1) x = Yspline[Nspline-1];
                           }
                           DS[i-imin] = exp(x);
                           DlnpR[i-imin] = y;
                           x = DS[i-imin]*(1.0+SL[i]);
                           DlnLR[i-imin] = -(log(x) + PS[i]/x);
                           logLy += (DlnLR[i-imin] - lnLR[i]);
                           logpy += (DlnpR[i-imin] - lnpR[i]);
               }
         
         
             
         }
         else if(alpha < a2) // line update
         {
             
            typ = 1;
            cL++;
             
            // pick a line to update
            k = (int)((double)(Nlines)*gsl_rng_uniform(r));
             
             alpha = gsl_rng_uniform(r);
             if(alpha > 0.5)
             {
              ylinef = linef[k] + gsl_ran_gaussian(r,0.1);
              ylineh = lineh[k]*(1.0+gsl_ran_gaussian(r,0.05));
              ylineQ = lineQ[k] + gsl_ran_gaussian(r,1.0);
             }
             else if (alpha > 0.3)
             {
              ylinef = linef[k] + gsl_ran_gaussian(r,0.4);
              ylineh = lineh[k]*(1.0+gsl_ran_gaussian(r,0.2));
              ylineQ = lineQ[k] + gsl_ran_gaussian(r,4.0);
             }
             else
             {
              ylinef = linef[k] + gsl_ran_gaussian(r,0.01);
              ylineh = lineh[k]*(1.0+gsl_ran_gaussian(r,0.01));
              ylineQ = lineQ[k] + gsl_ran_gaussian(r,0.2);
             }
             
            
             if(ylineQ < 10.0) flag = 1;
             if(ylineh < 1.0) flag = 1;
             
           
             spread = (1.0e-2*ylineQ);
             if(spread < 50.0) spread = 50.0;  // maximum half-width is f_resonance/20
             ydeltafmax = ylinef/spread;
             ylinew = 8.0*ydeltafmax;
             
             
            // need to cover old and new line
             imin = (int)((linef[k]-linew[k])*Tobs);
             imax = (int)((linef[k]+linew[k])*Tobs);
             i = (int)((ylinef-ylinew)*Tobs);
             if(i < imin) imin = i;
             i = (int)((ylinef+ylinew)*Tobs);
             if(i > imax) imax = i;
             if(imin < 0) imin = 0;
             if(imax > Ns-1) imax = Ns-1;
             
            // new line contribution
            f2 = ylinef*ylinef;
            f4 = f2*f2;
            for (i = imin; i <= imax; ++i)
               {
                    f = (double)(i)/Tobs;
                    y = 0.0;
                    x = fabs(f - ylinef);
                    if(x < ylinew)
                    {
                    z = 1.0;
                    if(x > ydeltafmax) z = exp(-(x-ydeltafmax)/ydeltafmax);
                    ff2 = f*f;
                    df2 = ylineQ*(1.0-ff2/f2);
                    df4 = df2*df2;
                    y += z*ylineh/(ff2/f2+df4);
                    }
                    DS[i-imin] = y;
               }
            
             
            // have to recompute and remove current line since lines can overlap
            f2 = linef[k]*linef[k];
            for (i = imin; i <= imax; ++i)
               {
                    f = (double)(i)/Tobs;
                    y = 0.0;
                    x = fabs(f - linef[k]);
                    if(x < linew[k])
                    {
                    z = 1.0;
                    if(x > deltafmax[k]) z = exp(-(x-deltafmax[k])/deltafmax[k]);
                    ff2 = f*f;
                    df2 = lineQ[k]*(1.0-ff2/f2);
                    df4 = df2*df2;
                    y += z*lineh[k]/(ff2/f2+df4);
                    }
                    DS[i-imin] -= y;
               }
            
             logpy = logpx;
             logLy = logLx;
             for (i = imin; i <= imax; ++i)
             {
                 x = SM[i]*(1.0+DS[i-imin]+SL[i]);
                 DlnLR[i-imin] = -(log(x) + PS[i]/x);
                 logLy += (DlnLR[i-imin] - lnLR[i]);
             }

             
             
         }
         else  // update hyper prior
         {
              typ = 2;
             
            
             alpha = gsl_rng_uniform(r);
             if(alpha > 0.6)
             {
                 pmy = pmx + gsl_ran_gaussian(r,0.001);
             }
             else if (alpha > 0.3)
             {
                 pmy = pmx + gsl_ran_gaussian(r,0.01);
             }
             else
             {
                 pmy = pmx + gsl_ran_gaussian(r,0.1);
             }
            
             
             if(pmy > mulmax || pmy < mulmin) flag = 1;
            
             
         }
         
         H = -1.0e60;
         if(flag == 0)
         {
           if(ptype == 0) H = (logLy-logLx) + (logpy/(2.0*pmy*pmy) - logpx/(2.0*pmx*pmx)) + (double)(Ns-2)*(log(pmx/pmy));
           if(ptype == 1) H = (logLy-logLx) + (logpy/pmy - logpx/pmx) +(double)(Ns-2)*(log(pmx/pmy));
         }
          
         alpha = log(gsl_rng_uniform(r));
         
         if(H > alpha)
         {
             logLx = logLy;
             logpx = logpy;
             
             if(typ == 0)
             {
             acS++;
             Xspline[k] = Yspline[k];
             for (i = imin; i < imax; ++i)
                {
                    SM[i] = DS[i-imin];  // replace this segment
                    lnLR[i] = DlnLR[i-imin];
                    lnpR[i] = DlnpR[i-imin];
                }
             }
             
             if(typ == 1)
             {
             acL++;
               
                linef[k] = ylinef;
                lineh[k] = ylineh;
                lineQ[k] = ylineQ;
                linew[k] = ylinew;
                deltafmax[k] = ydeltafmax;
                
               for (i = imin; i <= imax; ++i)
                {
                    SL[i] = DS[i-imin]+SL[i]; // delta this segment
                    lnLR[i] = DlnLR[i-imin];
                }
                 
             }
             
             if(typ == 2) pmx = pmy;
               
           }
          else
          {
           // have to re-set if not accepted
           if(typ == 0) Yspline[k] = Xspline[k];
          }
         
         //if(mc%100 == 0) printf("%d %e %e %f %f\n", mc, logLx, logpx, (double)acS/(double)(cS), (double)acL/(double)(cL));
         
         if(mc%10 == 0) fprintf(out, "%d %e %e %e %f %f\n", mc, logLx, logpx, pmx, (double)acS/(double)(cS), (double)acL/(double)(cL));
    
     }
    fclose(out);
    fclose(sfile);
    fclose(lfile);
    
    sprintf(command, "summary_%d.dat", m);
    sfile = fopen(command,"w");
    fprintf(sfile,"%d %d %d\n", Nspline, Nlines, MS);
    fclose(sfile);
    
    printf("There are %d lines\n", Nlines);
    
    /*
    out = fopen("specfit.dat","w");
    for (i = 0; i < Ns; ++i)
    {
        f = (double)(i)/Tobs;
        fprintf(out,"%e %e %e %e\n", f, SM[i], SL[i], SM[i]*(1.0+SL[i]));
    }
    fclose(out);
    
    // final check
    gsl_spline_init(cspline,ffit,Xspline,Nspline);
    SM[0] = exp(Xspline[0]);
    SM[Ns-1] = exp(Xspline[Nspline-1]);
    for (i = 1; i < Ns-1; ++i)
    {
        f = (double)(i)/Tobs;
        SM[i] = exp(gsl_spline_eval(cspline,f,acc));
    }
    
    
    out = fopen("speccheck.dat","w");
    for (i = 0; i < Ns; ++i)
    {
        f = (double)(i)/Tobs;
        y = 0.0;
        for (j = 0; j < Nlines; ++j)
        {
            x = fabs(f - linef[j]);
            if(x < linew[j])
            {
            z = 1.0;
            if(x > deltafmax[j]) z = exp(-(x-deltafmax[j])/deltafmax[j]);
            f2 = linef[j]*linef[j];
            ff2 = f*f;
            df2 = lineQ[j]*(1.0-ff2/f2);
            df4 = df2*df2;
            y += z*lineh[j]/(ff2/f2+df4);
            }
        }
        SL[i] = y;
        fprintf(out,"%e %e %e %e\n", f, SM[i], SL[i], SM[i]*(1.0+SL[i]));
        
    }
    fclose(out);
    */
    
    sp = (double*)malloc(sizeof(double)* (MS));
    
    // get upper and lower 90th
    ii = MS/5;
    
    sprintf(command, "PSD_%d_%d_%d.dat", (int)(Tobs), (int)ttrig, m);
    out = fopen(command,"w");
    for (i = 0; i < Ns; ++i)
    {
      for (j = 0; j < MS; ++j) sp[j] = SR[j][i];
      gsl_sort(sp, 1, MS);
      fprintf(out,"%.15e %.15e %.15e %.15e\n", (double)(i)/Tobs, sp[MS/2], sp[ii-1], sp[MS-ii-1]);
    }
    fclose(out);
    
    free(sp);
    
    sprintf(command, "spec_%d_%d_%d.dat", (int)(Tobs), (int)ttrig, m);
    out = fopen(command,"w");
    for (i = 0; i < Ns; ++i)
    {
        f = (double)(i)/Tobs;
        fprintf(out,"%.15e %.15e %.15e %.15e\n", f, SM[i]*(1.0+SL[i]), SM[i], PS[i]);
    }
    fclose(out);
    
    
    free(lnLR);
    free(DlnLR);
    free(DS);
    free(lnpR);
    free(DlnpR);
    free(lineh);
    free(linef);
    free(lineQ);
    free(linew);
    free(deltafmax);
    free(SM);
    free(SN);
    free(Smean);
    free_double_matrix(SR,MS);
    free(PS);
    free(times);
    free(data);

    
    return 0;

}

