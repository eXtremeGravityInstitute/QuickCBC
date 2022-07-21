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


#define MM 200000    // iterations of MCMC

#define verbose 1   // set to 0 for quiet run, 1 for verbose
#define dfmin 4.0  // minimum spline spacing
#define dfmax 32.0  // maximum spline spacing
#define smooth 8.0   // moving average
#define tol 0.2     // tolerance for difference in averages
#define mulmax 10.0   // max of smoothness hyper prior
#define mulmin 0.0001     // min of smoothness hyper prior

//OSX
// clang -Xpreprocessor -fopenmp -lomp -w -o SpecFit SpecFit.c Utilities.c -lgsl  -lm

// Linux
// gcc -std=gnu99 -fopenmp -w -o SpecFit SpecFit.c Utilities.c -lgsl -lm

void setupsline(int Nsline, int Nsten, double *sline, double *ffit, double *slopes, double *bb, double *cc, double *dd);
void delta_setupsline(int iu, int Nsline, int Nsten, double *sline, double *ffit, double *slopes,  double *bb, double *cc, double *dd);
void sline(int istart, int iend, int Nsline, int Nsten, double *farray,  double *SM, double asmooth, double *ffit, double *slinep, double *bb, double *cc, double *dd, double xx, int Nlook, double *lookup);
double ltc(double x, double xx, int Nlook, double *lookup);
double line(double f, double linef, double lineh, double linew, double deltafmax, double lineQ);

int main(int argc, char *argv[])
{


  int i, j, k, kk, ND, N, Ns, Nm, m, mc, dec, MK, MS, skip;
  int imin, imax, acS, acL, acH, cS, cL, cH;
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
  int Nsline, Nsten;
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
    
    if(verbose == 1)
       {
       out = fopen("spec.dat","w");
       for (i = 0; i < Ns; ++i)
       {
           f = (double)(i)/Tobs;
           fprintf(out,"%.15e %.15e %.15e %.15e\n", f, SN[i], SM[i], PS[i]);
       }
       fclose(out);
       }
    
    
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
    Nsline = 1;
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
                  Nsline++;
              }
          }
        
        
    }
    
     Nsline++;
    
     printf("There are %d spline knots\n", Nsline);
    
     double *ffitx, *ffity, *Xsline, *Ysline, *freqs;
     
     ffitx = (double*)malloc(sizeof(double)*(Nsline));
     ffity = (double*)malloc(sizeof(double)*(Nsline));
     Xsline = (double*)malloc(sizeof(double)*(Nsline));
     Ysline = (double*)malloc(sizeof(double)*(Nsline));
    
    freqs = (double*)malloc(sizeof(double)*(Ns));
    for (i = 0; i < Ns; ++i)   freqs[i] = (double)(i+1)/Tobs;
    
    // asmooth min for Nsten = 3 is 15
    // asmooth min for Nsten = 5 is 8
    // recommended settings ar Nsten =5 and asmooth between 8 and 15
    Nsten = 5;  // use either 3 or 5. Must be odd. Could use 7, 9 etc, but then very long
    
    // Note that Nsline must be >= Nsten
    
    double asmooth, mx, xx;
    double *lookup;
    
    int Nlook = 2000;
      
    asmooth = 10.0;
      
      // generate lookup table
      mx = 6.0;
      lookup = (double*)malloc(sizeof(double)*(Nlook));
      xx = mx/(double)(Nlook);
      for(i=0; i< Nlook; i++)
      {
          x = (double)(i)*xx;
          lookup[i] = log(2.0*cosh(x));
      }
    
    
    double *bx, *cx, *delx;
    double *by, *cy, *dely;
    double *slopex, *slopey;
    
    slopex = (double*)malloc(sizeof(double)*(Nsline));
    cx = (double*)malloc(sizeof(double)*(Nsline));
    bx = (double*)malloc(sizeof(double)*(Nsline));
    delx = (double*)malloc(sizeof(double)*(Nsline));
    
    slopey = (double*)malloc(sizeof(double)*(Nsline));
    cy = (double*)malloc(sizeof(double)*(Nsline));
    by = (double*)malloc(sizeof(double)*(Nsline));
    dely = (double*)malloc(sizeof(double)*(Nsline));
    
    ffitx[0] = 1.0/Tobs;
    Xsline[0] = log(S1[0]);
    
    Nsline = 1;
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
                  ffitx[Nsline] = (double)(i-j/2)/Tobs;
                  Xsline[Nsline] = log(S1[i-j/2]);
                  j = 0;
                  Nsline++;
              }
          }
        
        
    }
    
    // set up sline
    setupsline(Nsline, Nsten, Xsline, ffitx, slopex, bx, cx, delx);
    
     for (i = 0; i < Nsline; ++i)
        {
            //printf("%d %e %e\n", i, ffitx[i], Xsline[i]);
            Ysline[i] = Xsline[i];
            ffity[i] = ffitx[i];  // These never change (we don't update the locations or number of control points in this code)
            by[i] = bx[i];
            cy[i] = cx[i];
            dely[i] = delx[i];
            slopey[i] = slopex[i];
           // printf("%d %f %e %e %e %e %e\n", i, ffitx[i], slopex[i], ax[i], bx[i], cx[i], delx[i]);
        }
    
      
       if(verbose == 1)
       {
        out = fopen("control.dat","w");
         for (i = 0; i < Nsline; ++i)
         {
             fprintf(out,"%e %e\n", ffitx[i], exp(Xsline[i]));
         }
         fclose(out);
       }

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
                  linef[j] = freqs[ii];
                  lineh[j] = (max-1.0)*SM[ii];
                  lineQ[j] = sqrt(max)*linef[j]*Tobs/(double)(k);
                  
                    spread = (1.0e-2*lineQ[j]);
                    if(spread < 50.0) spread = 50.0;  // maximum half-width is f_resonance/50
                    deltafmax[j] = linef[j]/spread;
                    linew[j] = 8.0*deltafmax[j];
                  
                  //printf("%d %e %e %e %e %e\n", j, linef[j], lineh[j], linew, freqs[ii], lineQ[j]);
                 
              }
          }

      }
      
      printf("There are %d lines\n", Nlines);
    

    
      // initialize smooth spectrum
      sline(0, Ns, Nsline, Nsten, freqs, SM, asmooth, ffitx, Xsline, bx, cx, delx, xx, Nlook, lookup);
    
      // initialize line spectrum
      //out = fopen("lines.dat","w");
      for (i = 0; i < Ns; ++i)
      {
          f = freqs[i];
          y = 0.0;
          for (j = 0; j < Nlines; ++j) y += line(f, linef[j], lineh[j], linew[j], deltafmax[j], lineQ[j]);
          SL[i] = y;
          //fprintf(out,"%e %e %e\n", f, SM[i] + SL[i], SL[i]);
      }
     //fclose(out);
      
      
       if(verbose == 1)
       {
         out = fopen("specstart.dat","w");
         for (i = 0; i < Ns; ++i)
         {
             fprintf(out,"%e %e %e %e\n", freqs[i], SM[i], SL[i], SM[i]+SL[i]);
         }
         fclose(out);
       }
    
   // return 1;

      
      // This holds the updated values
      
      double *DSM, *DSL;
      
       DSM = (double*)malloc(sizeof(double)*(Ns));
       DSL = (double*)malloc(sizeof(double)*(Ns));
     
      double logLx, logLy, logpx, logpy, pmulx, pmuly;
      double H;
      
      pmulx = 1.0;
      pmuly = 1.0;
      
      // initialise the log likelihood
      logLx = 0.0;
      for (i = 0; i < Ns; ++i)
       {
           x = SM[i] + SL[i];
           DSM[i] = SM[i];
           DSL[i] = 0.0;
           logLx += -(log(x) + PS[i]/x);
       }
      
      // initialise the log prior on slopes
         logpx = 0.0;
         x = 0.0;
         for (i = 0; i < Nsline-2; ++i)
          {
              logpx += -fabs(slopex[i+1]-slopex[i]);
              x += fabs(slopex[i]);
          }
         x /= (double)(Nsline-2); // average magnitude of slope
        logpx /= x;  // makes it dimensionless
      
       printf("%e %e\n", logLx, logpx);
      
      
      acS = 0;
      acL = 0;
      acH = 0;
      cS = 1;
      cL = 1;
      cH = 1;
      
     if(verbose==1) out = fopen("schain.dat","w");
      
       for (mc = 0; mc < MM; ++mc)
       {
           logpy = logpx;
           pmuly = pmulx;
           
           // prune any weak lines
           if(mc > 0 && mc == MM/4)
           {
               j = 0;
               for (k = 0; k < Nlines; ++k)
               {
                   i = (int)(linef[k]*Tobs)-1;  // bin where line peaks
                   x = lineh[k]/SM[i];     // height of line relative to smooth
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
                  f = freqs[i];
                  y = 0.0;
                  for (j = 0; j < Nlines; ++j) y += line(f, linef[j], lineh[j], linew[j], deltafmax[j], lineQ[j]);
                  SL[i] = y;
               }
               
               // reset likelihood
               logLx = 0.0;
               for (i = 0; i < Ns; ++i)
               {
                   x = SM[i] + SL[i];
                   logLx -= (log(x) + PS[i]/x);
               }

               
            }
           
           
           alpha = gsl_rng_uniform(r);
           
           if(alpha < 0.4) // update sline
           {
               
            typ = 0;
            cS++;
              
            // important to make sure all the sline points and slopes are properly initialized
            // the proposed updates reset parts of slopey array, and the impact on the aa, bb etc
            // coefficients extends over several points
            for (i = 0; i < Nsline; ++i)
            {
                Ysline[i] =  Xsline[i];
                ffity[i] = ffitx[i]; // the knot locations don't change in this code
                slopey[i] = slopex[i];
                by[i] = bx[i];
                cy[i] = cx[i];
                dely[i] = delx[i];
            }
               
             // pick a knot to update
             k = (int)((double)(Nsline)*gsl_rng_uniform(r));
               
               alpha = gsl_rng_uniform(r);
               if(alpha > 0.5)
               {
                Ysline[k] =  Xsline[k] + gsl_ran_gaussian(r,0.05);
               }
               else if (alpha > 0.3)
               {
                Ysline[k] =  Xsline[k] + gsl_ran_gaussian(r,0.02);
               }
               else
               {
                Ysline[k] =  Xsline[k] + gsl_ran_gaussian(r,0.1);
               }
               

              // if control point k is updated, need to do a delta update on the sline model
              // region between ffit[k-(Nsten-2)] and ffit[k] is impacted. Care needs to be taken at boundary points.

               
               //setupsline(Nsline, Nsten, Ysline, ffity, slopey, by, cy, dely);
              delta_setupsline(k, Nsline, Nsten, Ysline, ffity, slopey, by, cy, dely);
               
                   
              if(k >= (Nsten-2))
              {
              imin = (int)(ffity[k-(Nsten-2)]*Tobs);
              }
              else
              {
              imin = 0;
              }
              
              if(k < Nsline-1)
               {
               imax = (int)(ffity[k+1]*Tobs)+1;
               }
              else
               {
               imax = Ns;
               }
               
               // start is special
               if(k < Nsten) imax = (int)(ffity[Nsten]*Tobs)+1;
                  
              sline(imin, imax, Nsline, Nsten, freqs, DSM, asmooth, ffity, Ysline, by, cy, dely, xx, Nlook, lookup);
               
               // compute change in likelihood
              logLy = logLx;
               
               //printf("%d \n", mc);
              for (i = imin; i < imax; ++i)
                {
                  x = SM[i] + SL[i];
                  y = DSM[i] + SL[i];
                  logLy +=  -(log(y) + PS[i]/y)+(log(x) + PS[i]/x);
                   // printf("%d %e %e\n", i, SM[i], DSM[i]);
                }
               
               logpy = 0.0;
                x = 0.0;
                for (i = 0; i < Nsline-2; ++i)
                 {
                     logpy += -fabs(slopey[i+1]-slopey[i]);
                     x += fabs(slopey[i]);
                 }
                x /= (double)(Nsline-2); // average magnitude of slope
               logpy /= x;  // makes it dimensionless
               
           }
           else if(alpha < 0.8) // update sline // line update
           {
               
              typ = 1;
              cL++;
               
               // pick a line to update
                          k = (int)((double)(Nlines)*gsl_rng_uniform(r));
                           
                           alpha = gsl_rng_uniform(r);
                           if(alpha > 0.5)
                           {
                            ylinef = linef[k] + gsl_ran_gaussian(r,df);
                            ylineh = lineh[k]*(1.0+gsl_ran_gaussian(r,0.05));
                            ylineQ = lineQ[k] + gsl_ran_gaussian(r,1.0);
                           }
                           else if (alpha > 0.3)
                           {
                               ylinef = linef[k] + gsl_ran_gaussian(r,4.0*df);
                            ylineh = lineh[k]*(1.0+gsl_ran_gaussian(r,0.2));
                            ylineQ = lineQ[k] + gsl_ran_gaussian(r,4.0);
                           }
                           else
                           {
                            ylinef = linef[k] + gsl_ran_gaussian(r,0.2*df);
                            ylineh = lineh[k]*(1.0+gsl_ran_gaussian(r,0.01));
                            ylineQ = lineQ[k] + gsl_ran_gaussian(r,0.2);
                           }
               
              if(ylineQ < 10.0) ylineQ = 10.0;
              
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
               
              // proposed line contribution
               for (i = imin; i <= imax; ++i) DSL[i] = line(freqs[i], ylinef, ylineh, ylinew, ydeltafmax, ylineQ);
                 
              // have to recompute and remove current line since lines can overlap
               for (i = imin; i <= imax; ++i) DSL[i] -= line(freqs[i], linef[k], lineh[k], linew[k], deltafmax[k], lineQ[k]);
              
               logLy = logLx;
               for (i = imin; i <= imax; ++i)
               {
                  x = SM[i] + SL[i];
                  y = SM[i] + SL[i] + DSL[i];
                  logLy +=  -(log(y) + PS[i]/y)+(log(x) + PS[i]/x);
               }

           }
           else  // update smoothness hyper parameter
           {
               typ = 2;
               cH++;
               pmuly = pmulx + gsl_ran_gaussian(r,0.01);
           }
           
           
           H = (logLy-logLx) + pmuly*(logpy - logpx);
           alpha = log(gsl_rng_uniform(r));
           
           if(pmuly < mulmin || pmuly > mulmax) H = -1.0e60;
           
           //printf("%d %d %e %e\n", mc, typ, logLx, logLy);
           
           if(H > alpha)
           {
               logLx = logLy;
               logpx = logpy;
               
               if(typ == 0)
               {
                   
               acS++;
               Xsline[k] = Ysline[k];
               ffitx[k] = ffity[k];
                   
               for (i = imin; i < imax; ++i)
                  {
                      SM[i] = DSM[i]; // replace
                  }
                   
                   for(i=0; i< Nsline; i++)
                    {
                      slopex[i] = slopey[i];
                      cx[i] = cy[i];
                      bx[i] = by[i];
                      delx[i] = dely[i];
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
                      SL[i] += DSL[i];   // delta this segment
                  }
                   
               }
               
               if(typ == 2)
               {
               acH++;
                 
                pmulx = pmuly;
                  
               }
                 
           }
           
          if(verbose == 1)  if(mc%1000 == 0) printf("%d %e %f %f %f %f\n", mc, logLx, pmulx, (double)acS/(double)(cS), (double)acL/(double)(cL), (double)acH/(double)(cH) );
           
          if(verbose == 1)  if(mc%100 == 0) fprintf(out, "%d %e %f %f %f %f\n", mc, logLx, pmulx, (double)acS/(double)(cS), (double)acL/(double)(cL), (double)acH/(double)(cH));
      
       }
       if(verbose == 1) fclose(out);
      
      printf("There are %d lines\n", Nlines);
       
      out = fopen("specfit.dat","w");
      for (i = 0; i < Ns; ++i)
      {
          f = freqs[i];
          fprintf(out,"%e %e %e %e\n", f, SM[i], SL[i], SM[i]+SL[i]);
      }
      fclose(out);
      
      
       if(verbose == 1)
       {
           // final check
        sline(0, Ns, Nsline, Nsten, freqs, SM, asmooth, ffitx, Xsline, bx, cx, delx, xx, Nlook, lookup);
          
           
      out = fopen("speccheck.dat","w");
      for (i = 0; i < Ns; ++i)
      {
          f = freqs[i];
          y = 0.0;
          for (j = 0; j < Nlines; ++j) y += line(f, linef[j], lineh[j], linew[j], deltafmax[j], lineQ[j]);
          SL[i] = y;
          fprintf(out,"%e %e %e %e\n", f, SM[i], SL[i], SM[i]+SL[i]);
          
      }
      fclose(out);
       }
    
    sprintf(command, "spec_%d_%d_%d.dat", (int)(Tobs), (int)ttrig, m);
    out = fopen(command,"w");
    for (i = 0; i < Ns; ++i)
    {
        f = (double)(i)/Tobs;
        fprintf(out,"%.15e %.15e %.15e %.15e\n", f, SM[i]+SL[i], SM[i], PS[i]);
    }
    fclose(out);
    


      free(freqs);
      free(lookup);
      free(bx);
      free(cx);
      free(delx);
      free(slopex);
      free(by);
      free(cy);
      free(dely);
      free(slopey);
      free(DSM);
      free(DSL);
      free(lineh);
      free(linef);
      free(lineQ);
      free(linew);
      free(deltafmax);
      free(SM);
      free(SN);
      free(PS);
      free(times);
      free(data);
      
      
    
    return 0;

}

double line(double f, double linef, double lineh, double linew, double deltafmax, double lineQ)
{
    double x, y, z;
    double ff2, f2, df2, df4;
    
             f2 = linef*linef;
    
             y = 0.0;
             x = fabs(f - linef);
             if(x < linew)
             {
             z = 1.0;
             if(x > deltafmax) z = exp(-(x-deltafmax)/deltafmax);
             ff2 = f*f;
             df2 = lineQ*(1.0-ff2/f2);
             df4 = df2*df2;
             y = z*lineh/(ff2/f2+df4);
             }
    
    return y;
                 
}

double ltc(double x, double xx, int Nlook, double *lookup)
{
    double y;
    int i;
    
    i = (int)(fabs(x)/xx);
    
    if(i > Nlook-1)  // ouside of lookup
    {
        y = log(2.0*cosh(x));
    }
    else
    {
        y = lookup[i];
    }
    
    return y;
    
}

void setupsline(int Nsline, int Nsten, double *sline, double *ffit, double *slopes, double *bb, double *cc, double *dd)
{
    int i, k;
    
      for(i=0; i< Nsline-1; i++)
        {
         slopes[i] = (sline[i+1]-sline[i])/(ffit[i+1]-ffit[i]);
         dd[i] = 0.5/(ffit[i+1]-ffit[i]);
        }
    
      for(i=0; i< Nsline-2; i++)
       {
           cc[i+1] = 0.5*(slopes[i+1]-slopes[i]);
           k = i+(Nsten-2);
           if(k > Nsline-2) k = Nsline-2;
           bb[i] = 0.5*(slopes[k]+slopes[i]);
       }
        
}

// when point iu is updated, the fit between control points at iu-(Nsten-1)/2;  and iu+(Nsten-1)/2 is changed
void delta_setupsline(int iu, int Nsline, int Nsten, double *slinep, double *ffit, double *slopes, double *bb, double *cc, double *dd)
{
    int ii, i, k;

    
      for(i=-1; i< 1; i++)
       {
        ii = iu+i;
        if(ii > -1 && ii < Nsline-1)
        {
        slopes[ii] = (slinep[ii+1]-slinep[ii])/(ffit[ii+1]-ffit[ii]);
        dd[ii] = 0.5/(ffit[ii+1]-ffit[ii]);
        }
       }
    
    //for(i=-(Nsten-1); i< 1; i++)
       for(i=-(Nsten); i< 1; i++)
        {
         ii = iu+i;
         if(ii > -1 && ii < Nsline-1)
         {
             k = ii+(Nsten-2);
             if(k > Nsline-2) k = Nsline-2;
             bb[ii] = 0.5*(slopes[k]+slopes[ii]);
         }
        }
    
    for(i=-2; i< 1; i++)
     {
      ii = iu+i;
      if(ii > -1 && ii < Nsline-2)
      {
          cc[ii+1] = 0.5*(slopes[ii+1]-slopes[ii]);
      }
     }
        
}

void sline(int istart, int iend, int Nsline, int Nsten, double *farray,  double *SM, double asmooth, double *ffit, double *slinep, double *bb, double *cc, double *dd, double xx, int Nlook, double *lookup)
{
    double f;
    int i, j, k, jx, jj, kk, ic;
    double x, y, z;
    double ax, bx;
    
    double *cx, *dx, *xv;
    
    cx = double_vector(Nsten);
    dx = double_vector(Nsten);
    xv = double_vector(Nsten);
    
    ic = (Nsten-1)/2;  // central point of stencil
    if(ic == 1) ic = 2;  // Nsten = 3 is a special case
    
    j = -ic+1;
   
    for (i = istart; i < iend; ++i)
    {
        f = farray[i];
        
          j -= 1;
          do
          {
           j++;
          } while(j < Nsline-Nsten && f > 0.5*(ffit[j+ic]+ffit[j+ic-1]));
         
            
        jx = j;
        if(j < 0) jx = 0;
        if(j > Nsline-Nsten) jx = Nsline-Nsten;
        
        for(k=0; k< Nsten-1; k++)
        {
            jj = jx+k;
            dx[k] = dd[jj];
            cx[k] = cc[jj];
            xv[k] = ffit[jj];
        }
        
        ax = slinep[jx];
        for(k=1; k< Nsten-1; k++)
         {
            ax -= cx[k]*fabs((xv[k]-xv[0]));
        }
        bx = bb[jx];
                   
        y = ax + bx*(f-xv[0]);
        for(k=1; k< Nsten-1; k++)
        {
            x = fabs(asmooth*dx[k]*(f-xv[k]));
            if(x > 6.0)
            {
            y += cx[k]*fabs(f-xv[k]);
            }
            else
            {
             y += cx[k]/(asmooth*dx[k])*ltc(x,xx,Nlook,lookup);
            }
        }
        
        SM[i] = exp(y);
        
    }

     free_double_vector(cx);
     free_double_vector(dx);
     free_double_vector(xv);
    
    
}




