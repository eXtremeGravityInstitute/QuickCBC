/*******************************************************************************************

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

**********************************************************************************************/


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <time.h>

#include "QuickCBC.h"
#include "ConstCBC.h"
#include "Utilities.h"
#include "IMRPhenomD.h"

#ifndef _OPENMP
#define omp ignore
#endif

//OSX
// clang -Xpreprocessor -fopenmp -lomp -w -o QuickCBC QuickCBC.c Utilities.c IMRPhenomD_internals.c IMRPhenomD.c -lgsl -lgslcblas  -lm

// Linux
// gcc -std=gnu99 -fopenmp -w -o QuickCBC QuickCBC.c IMRPhenomD_internals.c IMRPhenomD.c -lgsl -lgslcblas  -lm

//##############################################
//MT modifications

gsl_rng **rvec;
//##############################################

int main(int argc, char *argv[])
{
  int i, j, k, M, N, Nf, Nstep, Nclean, m, rs, tsi, tti;
    int ii, jj, kk, Nlines, id;
  int *mxc;
  int nt, bn;
  int oflag, flag;
  int imin, imax;
  double SNR, max;
  double Fplus, Fcross;
  double ciota, Fscale, Ap, Ac;
  double psi, alpha, sindelta;
  double junk, Tobs, fix, f, t, t0, dt, dtm, df, x, y, z, dx, dtx;
  double dfx, Q, fny, scale, dlnf;
  double *freqs, *ref;
  double *inp, *oup, *slice;
  double *H1dat, *L1dat;
  double *waveH, *waveL;
  double **D, **DW, **Dtime;
  double *Dds;
  double **SN, **SM;
  double *sdata;
  double *intime, *sqf;
  double sigmean, sigmedian;
  int subscale, octaves;
    int mmax;
    double SNRsq, SNRold, pmax;
    double SNRH, SNRL, pw;
   double t_rise, s1, s2, ascale, fac, Tpad;
    double av, var, Tfull;
  double ttrig, tstart, tstart_clean, Tclean, starttime, endtime, Dfmax;
    double q, mc, mt, eta, m1, m2;
    
    int Oflag;
    
  int modelprint;
    
  double *linef, *linew, *lineh, *lineQ;
    
    double *time;
    double *DHfull, *DLfull;
    double *DHcopy, *DLcopy;
	
  char filename[1024];
  char command[1024];
    char Dname[1024];
   

  int n;
    
    const gsl_rng_type * P;
    gsl_rng * r;

    gsl_rng_env_setup();
    
    P = gsl_rng_default;
    r = gsl_rng_alloc (P);


  FILE *in;
  FILE *ifp;
  FILE *out;
    
    struct Net *net  = malloc(sizeof(struct Net));

    gsl_rng_set(r, 18346443564);
    
    //##############################################
    //open MP modifications
    omp_set_num_threads(NC);
    rvec = (gsl_rng **)malloc(sizeof(gsl_rng *) * (NC+1));
    for(i = 0 ; i<= NC; i++){
        rvec[i] = gsl_rng_alloc(P);
        gsl_rng_set(rvec[i] , i);
    }
    //##############################################
    

    if(argc<4)
    {
        printf("./QuickCBC Tobs trig_time detector-list\n");
        printf("The detector list can be just one detector, e.g. 0 for H, 1 for L, 2 for V\n");
        printf("Can also enter 1 2 or 0 1 2 etc. The order doesn't matter, except that the first on the list becomes the reference detector\n");
        return 1;
    }
    
    // Note that H is always called 0, L is 1 and V is 2 for data read and antenna response. But the ifo order can be anything. The labels array handles the mapping. e.g. if the command line has 2 0 1, then ifo order is 0=Virgo, 1=H and 2=L.
    
    Tobs = atof(argv[1]);
    ttrig = atof(argv[2]);
    
    // Hour angle
    net->GMST = gmst(ttrig);
    
    net->Tobs = Tobs;
    
    printf("GMST = %f\n", net->GMST);
    
    if(Tobs < 4.0)
    {
        printf("Observation time too short - need at least 4 seconds of data\n");
        return(-1);
    }
    
    net->Nifo = argc-3;
    
    net->labels = int_vector(net->Nifo);
    for (i = 0; i < net->Nifo; ++i) net->labels[i] = atoi(argv[3+i]);
    net->tds = double_vector(net->Nifo);
    
    net->delays = double_matrix(net->Nifo,net->Nifo);
    
    pairs(net);
    time_delays(net);
    
    for (i = 0; i < net->Nifo; ++i)
    {
        for (j = 0; j < net->Nifo; ++j)
        {
            printf("delay %d-%d = %f ", net->labels[i], net->labels[j], net->delays[i][j]);
        }
        printf("\n");
    }
     printf("\n");
                                                                                                    
    
    for (i = 1; i < net->Nifo; ++i) printf("delay 0-%d = %f ", net->labels[i], net->tds[i]);
    printf("\n");
    
    starttime = floor(ttrig)-Tobs+2.0;
    
    if(fabs(ttrig-floor(ttrig)) < 1.0e-8) // integer merger time
    {
       net->tmax = Tobs-0.5; // upper peak time  (Trigger time is Tobs-2)
       net->tmin =  Tobs-3.5; // lower peak time (Trigger time is Tobs-2)
       if(net->tmin < 1.5) net->tmin = 1.5; // need at least one second of good data
    }
    else
    {
        net->tmax = Tobs-2.0+(ttrig-floor(ttrig))+0.5*twidth;
        net->tmin = Tobs-2.0+(ttrig-floor(ttrig))-0.5*twidth;
    }
    
    printf("%f %f\n", net->tmin, net->tmax);

    
  /*********************    This section reads in the data and the PSD estimates *************************/
    

    // Read in the frame data for the first detector and use it to identify the data length
    
    sprintf(command, "framed_%d_%d_%d.dat", (int)(Tobs), (int)ttrig, net->labels[0]);
    in = fopen(command,"r");

    N = -1;
    while(!feof(in))
    {
        fscanf(in,"%lf%lf", &x, &y);
        N++;
    }
    rewind(in);
    
    dt = Tobs/(double)(N);
    
    printf("%d\n", N);

    
    time = (double*)malloc(sizeof(double)* (N));
    D= double_matrix(net->Nifo,N);
    Dtime = double_matrix(net->Nifo,N);
    DW = double_matrix(net->Nifo,N);
    
    for (i = 0; i < N; ++i) fscanf(in,"%lf%lf", &time[i], &D[0][i]);
    fclose(in);

    for (j = 1; j < net->Nifo; ++j)
    {
    sprintf(command, "framed_%d_%d_%d.dat", (int)(Tobs), (int)ttrig, net->labels[j]);
    in = fopen(command,"r");
    for (i = 0; i < N; ++i) fscanf(in,"%lf%lf", &x, &D[j][i]);
    fclose(in);
    }
    
      // keep a copy of the time domain data for later glitch removal
      for (j = 0; j < net->Nifo; ++j)
       {
       for (i = 0; i < N; ++i) Dtime[j][i] = D[j][i];
       }
    
    SN = double_matrix(net->Nifo,N/2);  // full PSD (smooth plus lines)
    SM = double_matrix(net->Nifo,N/2);  // smooth PSD
    
    // read in the PSDs
    for (j = 0; j < net->Nifo; ++j)
    {
    sprintf(command, "spec_%d_%d_%d.dat", (int)(Tobs), (int)ttrig, net->labels[j]);
    in = fopen(command,"r");
    for (i = 0; i < N/2; ++i) fscanf(in,"%lf%lf%lf%lf", &x, &SN[j][i], &SM[j][i], &x);
    fclose(in);
    }
    
    // Tukey window parameter. Flat for (1-alpha) of data
    t_rise = 0.4; // Standard LAL setting
    alpha = (2.0*t_rise/Tobs);
    
    tukey_scale(&s1, &s2, alpha, N);
    
    // Apply Tukey window to each detector time series
    for (i = 0; i < net->Nifo; ++i)
    {
     tukey(D[i], alpha, N); // Tukey window
     gsl_fft_real_radix2_transform(D[i], 1, N); // FFT
    for (j = 0; j < N; ++j) DW[i][j] = D[i][j];
    }
    
    // Apply scaling to FFTed data
    fac = sqrt((double)(N/2))/sqrt(Tobs);
    for (i = 0; i < net->Nifo; ++i)
    {
        whiten(DW[i], SN[i], N);  // whiten
        gsl_fft_halfcomplex_radix2_inverse(DW[i], 1, N); // iFFT
        for (j = 0; j < N; ++j) DW[i][j] /= fac;
    }
    

    // output whitened data
    out = fopen("dataw.dat","w");
    for (j = 0; j < N; ++j)
    {
        fprintf(out,"%e ", (double)(j)*dt-Tobs+2.0);
        for (i = 0; i < net->Nifo; ++i) fprintf(out,"%e ", DW[i][j]);
        fprintf(out,"\n");
    }
    fclose(out);
    

    fac = Tobs/((double)(N)*(double)(N));
    
    
    sprintf(command, "pspec_%d_%d.dat", (int)(Tobs), (int)ttrig);
    out = fopen(command,"w");
    for (i = 0; i < N/2; ++i)
    {
        fprintf(out,"%.15e ", (double)(i)/Tobs);
        for (j = 0; j < net->Nifo; ++j) fprintf(out,"%.15e ", fac*2.0*(D[j][i]*D[j][i]+D[j][N-i]*D[j][N-i]));
        fprintf(out,"\n");
    }
    fclose(out);
    
    sprintf(command, "wcheck_%d_%d.dat", (int)(Tobs), (int)ttrig);
    out = fopen(command,"w");
    for (i = 0; i < N/2; ++i)
    {
        fprintf(out,"%.15e ", (double)(i)/Tobs);
        for (j = 0; j < net->Nifo; ++j) fprintf(out,"%.15e ", fac*2.0*(D[j][i]*D[j][i]+D[j][N-i]*D[j][N-i])/SN[j][i]);
        fprintf(out,"\n");
    }
    fclose(out);
    

    fac = sqrt(Tobs)/(double)(N);
    for (i = 0; i < N; ++i)
    {
        for (k = 0; k < net->Nifo; ++k) D[k][i] *= fac;
    }
    
      for (i = 0; i < net->Nifo; ++i)
       {
       for (j = 0; j < N; ++j) DW[i][j] = D[i][j];
       }
    
    
     struct PSD *psd  = malloc(sizeof(struct PSD));
     
    
     psd->Nspline = int_vector(net->Nifo);
     psd->Nlines = int_vector(net->Nifo);
     
     i = 0; k = 0;
     for (j = 0; j < net->Nifo; ++j)
        {
        sprintf(command, "summary_%d.dat", net->labels[j]);
        in = fopen(command,"r");
        fscanf(in,"%d%d%d", &psd->Nspline[j], &psd->Nlines[j], &psd->Nsample);
            if(psd->Nspline[j] > i) i = psd->Nspline[j];
            if(psd->Nlines[j] > k) k = psd->Nlines[j];
        fclose(in);
        }
     
     psd->xspline = double_tensor(net->Nifo,psd->Nsample,i);
     psd->ffit = double_matrix(net->Nifo,i);
     psd->linef = double_tensor(net->Nifo,psd->Nsample,k);
     psd->lineh = double_tensor(net->Nifo,psd->Nsample,k);
     psd->lineQ = double_tensor(net->Nifo,psd->Nsample,k);
     psd->linew = double_tensor(net->Nifo,psd->Nsample,k);
     psd->deltafmax = double_tensor(net->Nifo,psd->Nsample,k);
     
        for (j = 0; j < net->Nifo; ++j)
           {
           sprintf(command, "sfile_%d.dat", net->labels[j]);
           in = fopen(command,"r");
           for (i = 0; i < psd->Nsample; ++i)
            {
                for (k = 0; k < psd->Nspline[j]; ++k) fscanf(in, "%lf", &psd->xspline[j][i][k]);
            }
           fclose(in);
               
             sprintf(command, "ffile_%d.dat", net->labels[j]);
             in = fopen(command,"r");
             for (k = 0; k < psd->Nspline[j]; ++k) fscanf(in, "%lf", &psd->ffit[j][k]);
             fclose(in);
               
           }
     
     for (j = 0; j < net->Nifo; ++j)
        {
         
            sprintf(command, "lfile_%d.dat", net->labels[j]);
             in = fopen(command,"r");
            
          for (i = 0; i < psd->Nsample; ++i)
            {
             for (k = 0; k < psd->Nlines[j]; ++k)
             {
                 fscanf(in, "%lf%lf%lf%lf%lf", &psd->linef[j][i][k], &psd->lineh[j][i][k], &psd->lineQ[j][i][k], &psd->linew[j][i][k], &psd->deltafmax[j][i][k]);
             }
            }
            
           fclose(in);
        }
     
    
    
    
  /*********************   With the data read in, we can now do the CBC PE *************************/
    
    
    // These need to be passed back and forth between BayesWave and the CBC PE code
    // The CBC PE code also needs the PSDs and glitch removed residuals D
   
    double ***global;
    double **skyx;
    int *who;
    double *heat;
    double **paramx;
    double **pallx;
    double ***history;
    double ***historyall;
    RealVector *freq;
    FILE *chainE;
    FILE *chainI;
    FILE *chainA;
    FILE *chainS;
    
    chainS = fopen("searchchain.dat","w");
    chainE = fopen("extrinsicchain.dat","w");
    chainI = fopen("intrinsicchain.dat","w");
    chainA = fopen("allchain.dat","w");
    
    global = double_tensor(NQ,NM,N);
    skyx = double_matrix(NC+1,NS);
    who = int_vector(NC+1);
    heat = double_vector(NC+1);
    paramx = double_matrix(NC+1,NX+3*net->Nifo);
    history = double_tensor(NC+1,NH,NX);
    freq = CreateRealVector((N/2));
    mxc = int_vector(3);
    for (i = 0; i < 3; i++) mxc[i] = 0;
    
    freq->data[0] = 1.0/Tobs;
    for (i=1; i< N/2; i++) freq->data[i] = (double)(i)/Tobs;
    
    pallx =  double_matrix(NC+1,NP);
    historyall = double_tensor(NC+1,NH,NP);
    
    // The CBC_start function does a search to find the signal and initializes the arrays
    CBC_start(net, mxc, chainS, paramx, skyx, pallx, who, heat, history, global, freq, D, Dtime, SN, N, Tobs, r);
    
    /*
    in = fopen("maplike_all.dat","r");
    fscanf(in,"%lf", &x);
    for (i=0; i< NP; i++)  fscanf(in,"%lf", &pallx[1][i]);
    fclose(in);
    
    for(i=1; i<=NC; i++) who[i] = i;
    for(i=1; i<=NCC; i++) heat[i] = 1.0;
    for(i=NCC+1; i<=NC; i++) heat[i] = heat[i-1]*ladder; */

    MCMC_all(net, psd, mxc, Nall, chainA, pallx, who, heat, historyall, global, freq, D, SN, SM, N, Tobs, ttrig, r);
    
    //###############################################
    //MT modification
   	for(i =0 ;i<= NC; i++){
        gsl_rng_free(rvec[i]);
    }
    free(rvec);
    //###############################################
    
    
    
    return 0;

}

double gmst(double ttrig)
{
    double GMST, x, toff, lp, dsec;
    
    /* Code based on Javascript online calculator https://celnav.de/longterm.htm */
    
    lp = leap(ttrig);
    time_t ttime = ttrig+EPOCH_UNIX_GPS-lp;
    struct tm *timeinfo;
    /* Get current time, print it and modify */
    timeinfo = gmtime(&ttime);
    
    printf("Trigger time: %04d-%02d-%02d %02d:%02d:%02d (%s)\n", timeinfo->tm_year + 1900, timeinfo->tm_mon + 1, timeinfo->tm_mday, timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec, timeinfo->tm_zone);
    
    double year, month, day, hour, minute, second, dayfraction;
    
    year = (double)timeinfo->tm_year + 1900.0;
    month = (double)timeinfo->tm_mon + 1.0;
    day = (double)timeinfo->tm_mday;
    hour = (double)timeinfo->tm_hour;
    minute = (double)timeinfo->tm_min;
    second = (double)timeinfo->tm_sec;
    
    dayfraction = (hour + minute/60.0 + second/3600.0)/24.0;
    
    if(month <= 2.0) {year -=1.0; month += 12;}
    double A = floor(year/100.0);
    double B = 2.0-A+floor(A/4.0);
    double JD0h = floor(365.25*(year+4716.0))+floor(30.6001*(month+1.0))+day+B-1524.5;
    double JD = JD0h+dayfraction;
    
    //Julian centuries (GMT) since 2000 January 0.5
    double T = (JD-2451545.0)/36525.0;
    double T2 = T*T;
    double T3 = T*T2;
    double T4 = T*T3;
    double T5 = T*T4;
    
    double GHAAmean = (280.46061837+ 360.98564736629*(JD-2451545.0)+0.000387933*T2-T3/38710000.0);
    
    while(GHAAmean < 360.0) GHAAmean += 360.0;
    
    x = floor(GHAAmean/360.0);
    GHAAmean = GHAAmean - x*360.0;
    
    double GMSTdecimal = GHAAmean/15.0;
    double GMSTh = floor(GMSTdecimal);
    double GMSTmdecimal = 60.0*(GMSTdecimal-GMSTh);
    double GMSTm = floor(GMSTmdecimal);
    double GMSTsdecimal = 60.0*(GMSTmdecimal-GMSTm);
    double GMSTs = rint(1000.0*GMSTsdecimal)/1000.0;
    
    printf("GMST h %f m %f s %f\n", GMSTh, GMSTm, GMSTs);
    
    // convert to radians
    GMST = ((GMSTs/60.0 +GMSTm)/60.0+ GMSTh)*PIn/12.0;
    
    return GMST;

}

double leap(double tgps)
{
    int i;
    double x;
    // 18 leap seconds since 1980
    double lp[] = {46828800.0, 78364801.0, 109900802.0, 173059203.0, 252028804.0, 315187205.0, 346723206.0, \
        393984007.0,425520008.0, 457056009.0, 504489610.0, 551750411.0, 599184012.0, 820108813.0, \
        914803214.0, 1025136015.0, 1119744016.0, 1167264017.0};
    
    x = 0.0;
    for (i = 0; i < 18; ++i)
    {
        if(tgps > lp[i]) x += 1.0;
    }
    
    return(x);
    
}



void CBC_start(struct Net *net, int *mxc, FILE *chainS, double **paramx, double **skyx, double **pallx, int *who, double *heat, double ***history, double ***global, RealVector *freq, double **D, double **Dtime, double **SN, int N, double Tobs, gsl_rng *r)
{
    int i, j, jj, k, kk, id;
    double *params;
    double *logL, *logLsky, *logLstart, *logLfull;
    double q, mc, m2, m1, mt;
    double x, y, qx, tx, fac, alpha;
    double lMc, lMcmin, lMcmax, dlMc, lmsun;
    double lMcx, lMtx;
    double **data, ***wave;
    double *Larray, **Tarray;
    double *DD;
    double **WW;
    double ***DHc, ***DHs, ***HH;
    int imin, imax;
    int nt, bn;
    double dtx, dt;
    double *SNRsq;
    double **rho;
    double Lmax;
    int jmax;
    char command[1024];
    
    FILE *chain;
    FILE *in;
    FILE *out;
    
    dt = Tobs/(double)(N);
    
    
    // intialize the who, heat, history and counter arrays
    
    for(i=1; i<=NC; i++) who[i] = i;
    
    // run cold to force ML
    for(i=1; i<=NCC; i++) heat[i] = 0.25;
    for(i=NCC+1; i<=NC; i++) heat[i] = heat[i-1]*ladder;  // spacing can be big here since we start at low temperature
    // with 10 hot chains this gets the hottest up to heat > 1.
    
    for(j = 1; j <= NC; j++)
    {
        for(k = 0; k <= NH; k++)
        {
            for(i = 0; i < NX; i++)
            {
                history[j][k][i] = 0.0;
            }
        }
    }
    
    mxc[0] = 0;
    mxc[1] = 0;
    
    logL = double_vector(NC+1);
    logLstart = double_vector(NC+1);
    logLsky = double_vector(NC+1);
    logLfull = double_vector(NC+1);
    
    rho = double_matrix(NC+1,net->Nifo);
 

    params = (double*)malloc(sizeof(double)* (NX+3*net->Nifo));
    // whitened signal in each detector
    wave = double_tensor(NC+1,net->Nifo,N);
    
    Larray = double_vector(N);
    Tarray = double_matrix(net->Nifo,N);
    
    
    // sky parameter order
    //[0] alpha, [1] sin(delta) [2] psi [3] ellipticity [4] scale [5] dphi [6] dt
    // The sky paramters [4], [5] and [6] hold the shift in the waveform in the reference detector
    // [4] Holds the amplitude rescaling
    // [5] Holds the time shift
    // [6] Holds the phase shift
    
    
    params[0] = log(1.0*MSUN_SI);  // log(Mc)
    
    q = 1.5; // q = m1/m2
    mc = exp(params[0]);
    m2 = mc*pow(1.0+q, 0.2)/pow(q,0.6);
    m1 = q*m2;
    mt = m1+m2;
    
    params[1] = log(mt);   // log(Mt)
    params[2] = 0.0;  // chi1
    params[3] = 0.0;  // chi2
    params[4] = 0.0;  // phi0
    params[5] = Tobs/2.0; // peak time (at reference detector, not geocenter)
    params[6] = log(1.0e8 * PC_SI); //distance (referenced to a F=1 antenna response in the reference detector, need antenna pattern factor to get physical distance).
    
    for (i = 1; i < net->Nifo; ++i)
    {
        params[NX+(i-1)*3] = 0.0; // phase offset
        params[NX+(i-1)*3+1] = 0.0; // time offset
        params[NX+(i-1)*3+2] = 1.0; // amplitude ratio
    }
    
    
    freq->data[0] = 1.0/Tobs;
    for (i=1; i< N/2; i++) freq->data[i] = (double)(i)/Tobs;
    
    lMcmax = log(mcmax*MSUN_SI);
    lMcmin = log(mcmin*MSUN_SI);
    
    dlMc = (lMcmax-lMcmin)/(double)(NM);
    
    for(k = 0; k < NQ; k++)
    {
        for(j = 0; j < NM; j++)
        {
            for(i = -N/2; i < N/2; i++)
            {
                global[k][j][i+N/2] = cap;
            }
        }
    }
    
    
    // Initial search to find signal
     MCMC_intrinsic(net, 1, mxc, Nsearch, chainS, paramx, skyx, pallx, who, heat, history, global, freq, D, SN, N, Tobs, r);

    // find current max likelihood waveform
    Lmax = 0.0;
    for(k = 1; k <= NC; k++)
    {
    logL[k] = log_likelihood_intrinsic(net, D, paramx[k], freq, SN, N, Tobs);
     if(logL[k] > Lmax)
      {
        Lmax = logL[k];
        kk = k;
     }
    }
    
    // params[NX+(k-1)*3+1] = delt[k]-delt[0]; // time offset
    
    
    // record ML f-tf track
    double *tm, *fr;
    tm = double_vector(N);
    fr = double_vector(N);
    for (i = 0; i < N; i++) tm[i] = Tobs/(double)(N)*(double)(i);
    f_of_t(paramx[kk], tm, fr, N);
    for(j = 0; j < net->Nifo; j++)
    {
      sprintf(command, "tftrack_search_%d.dat", net->labels[j]);
      out = fopen(command,"w");
      if(j == 0)
      {
       for (i = 0; i < N; i++) fprintf(out,"%e %e\n",  tm[i]-Tobs+2.0, fr[i]);
      }
      else  // put in time shift
      {
      for (i = 0; i < N; i++) fprintf(out,"%e %e\n",  tm[i]-Tobs+2.0+params[NX+(j-1)*3+1], fr[i]);
      }
      fclose(out);
    }
    free_double_vector(tm);
    free_double_vector(fr);
    
    //Here we remove any glitches and re-compute the SNR. Uses best-fit template
    
    // compute waveform in each detector using maxL chain parameters
    double **twave;
    double *SM, *SX;
    
    twave = double_matrix(net->Nifo,N);
    
    templates(net, twave, freq, paramx[kk], N);
    
    //Restore time domain data
    for (j = 0; j < net->Nifo; ++j)
    {
    for (i = 0; i < N; ++i) D[j][i] = Dtime[j][i];
    }
    
    SM = (double*)malloc(sizeof(double)*(N/2));
    SX = (double*)malloc(sizeof(double)*(N/2));
    
    fac = sqrt(Tobs)/(double)(N);

    // the spec code uses a different scaling
    for (j = 0; j < net->Nifo; ++j)
    {
        for (i = 0; i < N; ++i) twave[j][i] /= fac;
    }
    
    // make Qscan for data - signal
     if(printQ == 1)
     {
     for (j = 0; j < net->Nifo; ++j)
       {
           qscanres(D[j], twave[j], SN[j], Tobs, N);
           sprintf(command, "cp Qtranres.dat Qres_%d.dat", net->labels[j]);
           system(command);
          
           sprintf(command, "cp Qtranres.dat Qtransform.dat", net->labels[j]);
           system(command);
           sprintf(command, "gnuplot Qscan.gnu");
           system(command);
           sprintf(command, "cp Qscan.png Qres_%d.png", net->labels[j]);
           system(command);
       }
     }
    
    // make Qscan for signal
    if(printQ == 1)
    {
    for (j = 0; j < net->Nifo; ++j)
       {
           qscanf(twave[j], SN[j], Tobs, N);
           sprintf(command, "cp Qsig.dat Qsig_%d.dat", net->labels[j]);
           system(command);
           sprintf(command, "cp Qsig.dat Qtransform.dat", net->labels[j]);
           system(command);
           sprintf(command, "gnuplot Qscan.gnu");
           system(command);
           sprintf(command, "cp Qscan.png Qsig_%d.png", net->labels[j]);
           system(command);
       }
    }
    
    // recompute the PSD and remove any glitches. The signal is removed for PSD estimation and
    // glitch finding, but remains in the cleaned data that is returned
    for (j = 0; j < net->Nifo; ++j)
    {
       specest(D[j], twave[j], N, N/2, dt, 1.0/(2.0*dt), SN[j], SM, SX);
       sprintf(command, "cp wglitch.dat wglitch_%d.dat", net->labels[j]);
       system(command);
        
        sprintf(command, "clean_%d.dat", net->labels[j]);
        out = fopen(command,"w");
        for (i = 0; i < N; ++i) fprintf(out,"%e %e\n", dt*(double)(i), D[j][i]);
        fclose(out);
    }
    
    if(printQ == 1)
    {
    // This section is just for making Q-scans of the de-glitched residual
    for (j = 0; j < net->Nifo; ++j)
    {
       for (i = 0; i < N; ++i) Dtime[j][i] = D[j][i]; // use the cleaned data
       specest(Dtime[j], twave[j], N, N/2, dt, 1.0/(2.0*dt), SX, SM, SX);
       sprintf(command, "cp Qtransform.dat Qfullres_%d.dat", net->labels[j]);
       system(command);
       sprintf(command, "gnuplot Qscan.gnu");
       system(command);
       sprintf(command, "cp Qscan.png Qfullres_%d.png", net->labels[j]);
       system(command);
    }
    
    // This section is just for making Q-scans of the de-glitched data
    for (j = 0; j < net->Nifo; ++j)
    {
        for (i = 0; i < N; ++i) twave[j][i] = 0.0;  // don't subtract the signal
        for (i = 0; i < N; ++i) Dtime[j][i] = D[j][i]; // use the cleaned data
       specest(Dtime[j], twave[j], N, N/2, dt, 1.0/(2.0*dt), SX, SM, SX);
       sprintf(command, "cp Qtransform.dat Qclean_%d.dat", net->labels[j]);
       system(command);
       sprintf(command, "gnuplot Qscan.gnu");
       system(command);
       sprintf(command, "cp Qscan.png Qclean_%d.png", net->labels[j]);
       system(command);
    }
    // end Qscan making section
    }
    
    double t_rise = 0.4; // Standard LAL setting
    alpha = (2.0*t_rise/Tobs);
    
    // Apply Tukey window to each detector time series and FFT
    for (i = 0; i < net->Nifo; ++i)
    {
     tukey(D[i], alpha, N); // Tukey window
     gsl_fft_real_radix2_transform(D[i], 1, N); // FFT
    }
    
    for (i = 0; i < N; ++i)
    {
        for (k = 0; k < net->Nifo; ++k) D[k][i] *= fac;
    }
    
    // kill contribution below fmin Hz
    j = (int)(fmin*Tobs);
    for (i = 0; i < j; ++i)
    {
        for (k = 0; k < net->Nifo; ++k) SN[k][i] = 1.0;
    }
    
    free(SM);
    free(SX);
    free_double_matrix(twave,net->Nifo);
    
    // re-set with multiple cold chains for regular MCMC
    for(i=1; i<=NCC; i++) heat[i] = 1.0;
    for(i=NCC+1; i<=NC; i++) heat[i] = heat[i-1]*ladder;
    
    // refine the intrinsic parameters now any glitches have been removed
    MCMC_intrinsic(net, 0, mxc, Nintrinsic, chainS, paramx, skyx, pallx, who, heat, history, global, freq, D, SN, N, Tobs, r);
    
    imin = (int)((fmin*Tobs));
    imax = (int)((fmax*Tobs));
    
    upsample(N, Tobs, &nt, &bn);
    
    dtx = Tobs/(double)(bn);
    
    printf("Time steps = %d  time resolution = %f\n", nt, dtx);
  
    DD = double_vector(net->Nifo);
    WW = double_matrix(NC+1,net->Nifo);
    
    DHc = double_tensor(NC+1,net->Nifo,nt);
    DHs = double_tensor(NC+1,net->Nifo,nt);
    HH = double_tensor(NC+1,net->Nifo,3);
    
    // whitened data in each detector
    data = double_matrix(net->Nifo,N);
    
    for (id = 0; id < net->Nifo; ++id)
    {
        data[id][0] = 0.0;
        data[id][N/2] = 0.0;
    }
    
    for (i = 1; i < N/2; ++i)
    {
        for (id = 0; id < net->Nifo; ++id)
        {
            x = 1.0/sqrt(SN[id][i]);
            data[id][i] = D[id][i]*x;
            data[id][N-i] = D[id][N-i]*x;
        }
    }
    
 
    Lmax = 0.0;
    
    
    
    
    printf("Finding sky location\n");
    
    // Find sky location and set up whitend signal arrays
    #pragma omp parallel for
    for(j = 1; j <= NC; j++)
    {
        int mtid, mti;
        double mty, Scale;
        double *SNRsq;
        double **hwave;
        
        SNRsq = double_vector(net->Nifo);
        hwave = double_matrix(net->Nifo,N);
        
        templates(net, hwave, freq, paramx[j], N);
        
        // might want to put a catch here for any chains that didn't reach a decent logL and re-run the rearch phase
        logLstart[j] = log_likelihood_intrinsic(net, D, paramx[j], freq, SN, N, Tobs);
        
        //printf("%d %f\n", j, logLstart[j]);
        
        mty = 0.0;
        for (mtid = 0; mtid < net->Nifo; ++mtid)
        {
            SNRsq[mtid] = 0.0;
           for (mti = 1; mti < N/2; ++mti)
            {
            // The sky mapping code puts in the amplitude and phase shifts. The signals only differ due to the whitening
             SNRsq[mtid] += 4.0*(hwave[0][mti]*hwave[0][mti]+hwave[0][N-mti]*hwave[0][N-mti])/SN[mtid][mti];
            }
            Scale = 1.0;
            if(mtid > 0) Scale = paramx[j][(mtid-1)*3+NX+2]*paramx[j][(mtid-1)*3+NX+2];
            SNRsq[mtid] *= Scale;
            mty += SNRsq[mtid];
        }
        
       // printf("%f %f\n", SNRsq[0], SNRsq[1]);
        
        // find a sky location roughly consistent with the time delays
        skyring(net, paramx[j], skyx[j], pallx[j], SNRsq, freq, D, SN, N, Tobs, rvec[j]);
        
         // find a sky location roughly consistent with the time delay, amplitude ratio and phase difference
        //skystart(net, paramx[j], skyx[j], pallx[j], SNRsq, freq, D, SN, N, Tobs, rvec[j]);
        
        dshifts(net, skyx[j], paramx[j]);
        
        pmap(net, pallx[j], paramx[j], skyx[j]);
         
         // get the geocenter reference template. This does depend on the assumed sky location, hence follows skystart
         geotemplate(hwave[0], freq,  pallx[j], N);
         
        for (mtid = 0; mtid < net->Nifo; ++mtid)
        {
            wave[j][mtid][0] = 0.0;
            wave[j][mtid][N/2] = 0.0;
            for (mti = 1; mti < N/2; ++mti)
            {
                // The sky mapping code puts in the amplitude and phase shifts. The signals only differ due to the whitening
                mty = 1.0/sqrt(SN[mtid][mti]);
                wave[j][mtid][mti] = hwave[0][mti]*mty;
                wave[j][mtid][N-mti] = hwave[0][N-mti]*mty;
            }
        }
        
        skylikesetup(net, data, wave[j], DD, WW[j], DHc[j], DHs[j], Tobs, N, bn, nt, imin, imax);
        fisherskysetup(net, wave[j], HH[j], Tobs, N);
        
        logLsky[j] = skylike(net, skyx[j], DD, WW[j], DHc[j], DHs[j], dtx, nt, 0);
        
        logL[j] = log_likelihood_intrinsic(net, D, paramx[j], freq, SN, N, Tobs);
        
        logLfull[j] = log_likelihood_full(net, D, pallx[j], freq, SN, rho[j], N, Tobs);
        
        printf("%d logL initial %f intrinsic %f full %f sky %f\n", j, logLstart[j], logL[j], logLfull[j], logLsky[j]);
        
        free_double_vector(SNRsq);
        free_double_matrix(hwave,net->Nifo);
        
    }
    

    
    chain = fopen("searchsky.dat","w");
    skymcmc(net, Nsky, mxc, chain, paramx, skyx, pallx, who, heat, dtx, nt, DD, WW, DHc, DHs, HH, Tobs, r);
    fclose(chain);
    
    // make a map
    sprintf(command, "source sky.sh searchsky.dat 0");
    system(command);
    
    for(j = 1; j <= NC; j++)
    {
     pmap(net, pallx[j], paramx[j], skyx[j]);
     logL[j] = log_likelihood_intrinsic(net, D, paramx[j], freq, SN, N, Tobs);
     logLfull[j] = log_likelihood_full(net, D, pallx[j], freq, SN, rho[j], N, Tobs);
    }

    
    for(jj = 1; jj <= NC; jj++)
    {
        printf("\n");
        printf("%d logL initial %f intrinsic %f full %f\n", jj, logLstart[jj], logL[jj], logLfull[jj]);
        for (id = 0; id < net->Nifo; ++id)
        {
            printf("ifo %d %f ", net->labels[id], rho[jj][id]);
        }
         printf("\n");
    }
    
    
    x = 0.0;
    k = who[1];
    for(jj = 1; jj <= NC; jj++)
    {
        if(logLfull[jj] > x)
        {
            x = logLfull[jj];
            k = jj;
        }
    }
    
    printf("%d %f\n", k, x);
    
    // set threshold that all chains start withing 20% of highest likelihood
    x *= 0.8;
    
    for(jj = 1; jj <= NC; jj++)
    {
        
        printf("%d %f  ", jj, logLfull[jj]);
        
        if(logLfull[jj] < x)
        {
            for(i = 0; i < NP; i++)
            {
                pallx[jj][i] = pallx[k][i];
            }
        }
        
        logLfull[jj] = log_likelihood_full(net, D, pallx[jj], freq, SN, rho[jj], N, Tobs);
        
        printf("%f\n", logLfull[jj]);

    }
    
    
    
    // re-set with multiple cold chains for regular MCMC
     for(i=1; i<=NCC; i++) heat[i] = 1.0;
     for(i=NCC+1; i<=NC; i++) heat[i] = heat[i-1]*ladder;
 
    
    free_double_matrix(rho, NC+1);
    free(params);
    free_double_tensor(wave,NC+1,net->Nifo);
    free_double_matrix(Tarray,net->Nifo);
    
    free_double_vector(logL);
    free_double_vector(logLfull);
    free_double_vector(logLstart);
    free_double_vector(logLsky);
    
    free_double_vector(DD);
    free_double_matrix(WW,NC+1);
    free_double_tensor(DHc,NC+1,net->Nifo);
    free_double_tensor(DHs,NC+1,net->Nifo);
    free_double_tensor(HH,NC+1,net->Nifo);
    free_double_matrix(data,net->Nifo);

    
 return;
    
}

void CBC_update(struct Net *net, int *mxc, FILE *chainI, FILE *chainE, double **paramx, double **skyx, double **pallx, int *who, double *heat, \
                double ***history, double ***global, RealVector *freq, double **D, double **SN, int N, double Tobs, gsl_rng *r)
{
    int i, j, k, id;
    double *params;
    double logL;
    double q, mc, m2, m1, mt;
    double x, qx, tx;
    double lMc, lMcmin, lMcmax, dlMc, lmsun;
    double lMcx, lMtx;
    double **data, ***wave;
    double **hwave;
    double *DD;
    double **WW;
    double ***DHc, ***DHs, ***HH;
    double *rho;
    int imin, imax;
    int nt, bn;
    double dtx;
     double Fc, Fp, Fs, y;
    
    hwave = double_matrix(net->Nifo,N);
    // whitened signal in each detector
    wave = double_tensor(NC+1,net->Nifo,N);

    
    imin = (int)((fmin*Tobs));
    imax = (int)((fmax*Tobs));
    
    upsample(N, Tobs, &nt, &bn);
    
    dtx = Tobs/(double)(bn);
    
    rho = double_vector(net->Nifo);
    DD = double_vector(net->Nifo);
    WW = double_matrix(NC+1,net->Nifo);
    
    DHc = double_tensor(NC+1,net->Nifo,nt);
    DHs = double_tensor(NC+1,net->Nifo,nt);
    HH = double_tensor(NC+1,net->Nifo,3);
    
    // whitened data in each detector
    data = double_matrix(net->Nifo,N);

        for (id = 0; id < net->Nifo; ++id)
        {
            data[id][0] = 0.0;
            data[id][N/2] = 0.0;
            for (i = 1; i < N/2; ++i)
            {
            x = 1.0/sqrt(SN[id][i]);
            data[id][i] = D[id][i]*x;
            data[id][N-i] = D[id][N-i]*x;
            }
        }

    

    
    // Now alternate between bunches of intrinsic and extrinsic updates
    
    for(k = 0; k < 5; k++)
    {
        
        MCMC_intrinsic(net, 0, mxc, 200, chainI, paramx, skyx, pallx, who, heat, history, global, freq, D, SN, N, Tobs, r);
        
    
        // Set up whitend signal arrays
        for(j = 1; j <= NC; j++)
        {
            
            pmap(net, pallx[j], paramx[j], skyx[j]);
            
            // get the geocenter reference template. This does depend on the assumed sky location, hence follows skystart
            geotemplate(hwave[0], freq,  pallx[j], N);
            
            for (id = 0; id < net->Nifo; ++id)
            {
                wave[j][id][0] = 0.0;
                wave[j][id][N/2] = 0.0;
                for (i = 1; i < N/2; ++i)
                {
                    // The sky mapping code puts in the amplitude and phase shifts. The signals only differ due to the whitening
                    x = 1.0/sqrt(SN[id][i]);
                    wave[j][id][i] = hwave[0][i]*x;
                    wave[j][id][N-i] = hwave[0][N-i]*x;
                }
            }
            
            
            skylikesetup(net, data, wave[j], DD, WW[j], DHc[j], DHs[j], Tobs, N, bn, nt, imin, imax);
            fisherskysetup(net, wave[j], HH[j], Tobs, N);

            
            skyx[j][4] = 1.0;
            skyx[j][5] = 0.0;
            skyx[j][6] = 0.0;
            
            x = skylike(net, skyx[j], DD, WW[j], DHc[j], DHs[j], dtx, nt, 0);
            
            logL = log_likelihood_intrinsic(net, D, paramx[j], freq, SN, N, Tobs);
            
            y = log_likelihood_full(net, D, pallx[j], freq, SN, rho, N, Tobs);
            
           // printf("%d %f %f %f\n", j, x, logL, y);
            
           
            
        }
        
        skymcmc(net, 2000, mxc, chainE, paramx, skyx, pallx, who, heat, dtx, nt, DD, WW, DHc, DHs, HH, Tobs, r);
        
        
    }
    
    
    free_double_vector(rho);
    free_double_matrix(hwave,net->Nifo);
    free_double_tensor(wave,NC+1,net->Nifo);

    free_double_vector(DD);
    free_double_matrix(WW,NC+1);
    free_double_tensor(DHc,NC+1,net->Nifo);
    free_double_tensor(DHs,NC+1,net->Nifo);
    free_double_tensor(HH,NC+1,net->Nifo);
    free_double_matrix(data,net->Nifo);
    
    
    return;
    
}



void MCMC_intrinsic(struct Net *net, int lmax, int *mxc, int M, FILE *chain, double **paramx, double **skyx, double **pallx, int *who, double *heat, double ***history, double ***global, RealVector *freq, double **D, double **SN, int N, double Tobs, gsl_rng *r)
{
    int i, j, q, k, ip,  mc, flag;
    int *m;
    double *logLx, logLy, logL, x;
    double **paramy;
    double *pref;
    double *max, *min;
    double *href;
    double *pml;
    double alpha, beta, H;
    double Mchirp, Mtot, M1, M2, ciota;
    double eta, dm, m1, m2, chieff, DL, f;
    double a, b, c, d;
    double sqh;
    double Lmax, mch;
    double logpx, logpy;
    double pMcMx, pMcMy;
    double dA, dTD, dP, z, y;
    double leta;
    int bin;
    int **av, **cv;
    int typ, hold;
    int scount, sacc, mcount;
    double *Fscale;
    double pDx, pDy, qx, qy;
    double maxL;
    double *pmax;
    
    FILE *like;
    FILE *chainx;
    
    double ***fish, ***evec;
    double **ejump;
    
    av = int_matrix(5,NC+1);
    cv = int_matrix(5,NC+1);
    
    m = int_vector(NC+1);
    
    ejump = double_matrix(NC+1,NX);
    fish = double_tensor(NC+1,NX,NX);
    evec = double_tensor(NC+1,NX,NX);

    pmax = double_vector(NX+3*net->Nifo);
    paramy = double_matrix(NC+1,NX+3*net->Nifo);
    logLx = double_vector(NC+1);
    Fscale = double_vector(NC+1);
    
    // reference antenna pattern scaling used to convert reference detector frame amplitude to a physical distance
    if(lmax == 0)
    {
      for(k = 1; k <= NC; k++) Fscale[k] =  Fmag(skyx[k], net->GMST, net->labels[0]);
    }
    else
    {
        chainx = fopen("maxchain.dat","w");
        for(k = 1; k <= NC; k++)
        {
            Fscale[k] =  1.0;
        }
    }
    
    // prior boundaries
    max = (double*)malloc(sizeof(double)* (NX));
    min = (double*)malloc(sizeof(double)* (NX));
    
    max[0] = log(mcmax*MSUN_SI);  // Mc is at most 0.435 times Mt
    max[1] = log(mtmax*MSUN_SI);  // Mt
    max[2] = smax;  // chi1
    max[3] = smax;  // chi2
    max[4] = PI;  // phi0
    max[5] = net->tmax; // peak time  (Trigger time is Tobs-2)
    max[6] = log(DLmax * PC_SI); //distance
    
    
    min[0] = log(mcmin*MSUN_SI);  // Mc
    min[1] = log(mtmin*MSUN_SI);  // Mt
    min[2] = -smax;  // chi1
    min[3] = -smax;  // chi2
    min[4] = 0.0;  // phi0
    min[5] = net->tmin; // peak time (Trigger time is Tobs-2)
    min[6] = log(DLmin * PC_SI); //distance

    /*
    for(k = 0; k <= 6; k++) printf("%e ", max[k]);
    printf("\n");
    for(k = 0; k <= 6; k++) printf("%e ", min[k]);
    printf("\n");
     */
    
    /*
    chain = fopen("test.dat","w");
    for(i = 0; i < 10000; i++)
    {
        qx = globe(global, max, min, Tobs, paramx[1], N, r);
        Mchirp = exp(paramx[1][0])/MSUN_SI;
        Mtot = exp(paramx[1][1])/MSUN_SI;
        fprintf(chain, "%d %f %f %f %f\n", i, Mchirp, Mtot, Tobs-paramx[1][5], qx);
    }
    fclose(chain); */
    
    // random start drawn from the global
    if(lmax == 1)
    {
        for(k = 1; k <= NC; k++)
        {
             qx = globe(global, max, min, Tobs, paramx[k], N, r);
        }
    }
    
    // initialize the likelihoods
    for(k = 1; k <= NC; k++)
    {
    logLx[k] = 0.0;
    if(lhold == 0)
       {
        if(lmax == 0)
        {
            logLx[k] = log_likelihood_intrinsic(net, D, paramx[k], freq, SN, N, Tobs);
        }
        else
        {
            // GlitchSafe maximization
            logLx[k] = log_likelihood_max(net, D, paramx[k], freq, SN, N, Tobs, min[5], max[5], 0);
        }
       }
    Fisher_Full(net, fish[k], NX, paramx[k], freq, SN, N, Tobs);
    FisherEvec(fish[k], ejump[k], evec[k], NX);
    }
    
    
    printf("Intrinsic MCMC\n");

    for(j = 1; j <= 4; j++)
    {
     for(k=1; k <= NC; k++)
     {
         av[j][k] = 0;
         cv[j][k] = 0;
     }
    }
    
    scount = 1;
    sacc = 0;
    mcount = 1;
    
    for(k=1; k <= NC; k++) m[k] = 0;
    
   // q = who[1];
   // printf("%f %f %f %f %f %f\n", paramx[q][4], paramx[q][5], paramx[q][6], paramx[q][7], paramx[q][8], paramx[q][9]);
    
    maxL = -1.0e20;
           
    
    for(mc = 1; mc <= M; mc++)
    {
        if(lmax == 1 && mc%510==0)
        {
            x = -1.0;
            for(k = 1; k <= NC; k++)
            {
             if(logLx[k] > x)
             {
                 x = logLx[k];
                 j = k;
             }
            }
            if(x > 25.0)
            {
              for(k = 1; k <= NC; k++)
               {
                if(logLx[k] < 32.0)
                {
                    for(i = 0; i < NX+3*net->Nifo; i++)
                    {
                        paramx[k][i] = paramx[j][i];
                    }
                    logLx[k] = logLx[j];
                }
               }
            }
            
            // always map hottest chain to highest likelihood chain
            logLx[who[NC]] = logLx[j];
            for(i = 0; i < NX+3*net->Nifo; i++)
            {
                paramx[who[NC]][i] = paramx[j][i];
            }
            
        }
        
        if(mc%500==0 )
        {
            // update the Fisher matrices
            #pragma omp parallel for
            for(k = 1; k <= NC; k++)
            {
                Fisher_Full(net, fish[k], NX, paramx[k], freq, SN, N, Tobs);
                FisherEvec(fish[k], ejump[k], evec[k], NX);
            }
        }
        
        alpha = gsl_rng_uniform(r);
        
        if((NC > 1) && (alpha < 0.2))  // decide if we are doing a MCMC update of all the chains or a PT swap
        {
            
            // chain swap
            scount++;
            
            alpha = (double)(NC-1)*gsl_rng_uniform(r);
            j = (int)(alpha) + 1;
            beta = exp((logLx[who[j]]-logLx[who[j+1]])/heat[j+1] - (logLx[who[j]]-logLx[who[j+1]])/heat[j]);
            alpha = gsl_rng_uniform(r);
            if(beta > alpha)
            {
                hold = who[j];
                who[j] = who[j+1];
                who[j+1] = hold;
                sacc++;
            }
            
        }
        else      // MCMC update
        {
            
        mcount++;
            
            for(j = 1; j <= NC; j++)
            {
                for(i = 0; i < NX+3*net->Nifo; i++)
                {
                    paramy[j][i] = paramx[j][i];
                }
            }
        
    #pragma omp parallel for
    for(k=1; k <= NC; k++)
    {
     updatei(k, net, lmax, logLx, paramx, paramy, min, max, Fscale, who, heat, history, global, freq, D, SN, ejump, evec, N, Tobs, cv, av, rvec[k]);
    }
    
    // add to the history file
    if(mc%20 == 0)
    {
     for(k=1; k <= NC; k++)
     {
        q = who[k];
        i = m[k]%1000;
        // the history file is kept for each temperature
        for(j=0; j<NX; j++) history[k][i][j] = paramx[q][j];
        m[k]++;
     }
    }
            
        }
            
        // update maxL parameters
        if(logLx[who[1]] > maxL)
        {
            maxL = logLx[who[1]];
            for(i = 0; i < NX+3*net->Nifo; i++) pmax[i] = paramx[who[1]][i];
        }
        
        
        if(mc%10 == 0)
        {
            for(ip = 1; ip <= NCC; ip++)
            {
            q = who[ip];
            Mchirp = exp(paramx[q][0])/MSUN_SI;
            Mtot = exp(paramx[q][1])/MSUN_SI;
            eta = pow((Mchirp/Mtot), (5.0/3.0));
            dm = sqrt(1.0-4.0*eta);
            m1 = Mtot*(1.0+dm)/2.0;
            m2 = Mtot*(1.0-dm)/2.0;
            chieff = (m1*paramx[q][2]+m2*paramx[q][3])/Mtot;
            ciota = skyx[q][3];
            
            
            // pall  [0] log(Mc) [1] log(Mt) [2] chi1 [3] chi2 [4] phi0 [5] tp [6] log(DL) [7] alpha [8] sindelta [9] psi [10] ciota
            pmap(net, pallx[q], paramx[q], skyx[q]);
            
            DL = Fscale[q]*exp(paramx[q][6])/(1.0e6*PC_SI);  // includes conversion to Mpc
            z = z_DL(DL);
            
           // printf("%e %e\n", DL, exp(pallx[q][6])/(1.0e6*PC_SI));
            
            // counter, log likelihood, chirp mass, total mass, effective spin, geocenter GW phase, geocenter arrival time, distance, RA , sine of DEC,
            // polarization angle, cos inclination
            
            if(lmax == 1)
            {
                logL = log_likelihood_intrinsic(net, D, paramx[q], freq, SN, N, Tobs);
                fprintf(chainx,"%d %e %e %e %e %e %e %e ", mc/10, logLx[q], logL, Mchirp, Mtot, chieff, m1, m2);
                for(i = NX; i < NX+3*(net->Nifo-1); i++)  fprintf(chainx,"%e ", paramx[q][i]);
                fprintf(chainx,"\n");
            }
            else
            {
             fprintf(chain,"%d %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n", mxc[0], logLx[q], Mchirp, Mtot, chieff, pallx[q][4], \
                                     pallx[q][5]-(Tobs-2.0), DL, skyx[q][0], skyx[q][1], \
                                     skyx[q][2], ciota, z, Mchirp/(1.0+z), Mtot/(1.0+z), m1/(1.0+z), m2/(1.0+z), m1, m2);
             mxc[0] += 1;
            }
            }
        }
        
        
       if(mc%100 == 0)
       {
           q = who[1];
           Mchirp = exp(paramx[q][0])/MSUN_SI;
           Mtot = exp(paramx[q][1])/MSUN_SI;
           
           f = fbegin(paramx[q]);
           
           k = mxc[0];
           if(lmax == 1) k = mc;
           
           printf("%d %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", k, logLx[q], Mchirp, Mtot, f, paramx[q][5], paramx[q][6], paramx[q][7], paramx[q][8], paramx[q][9], \
                     (double)(sacc)/(double)(scount), (double)(av[1][1])/(double)(cv[1][1]), (double)(av[2][1])/(double)(cv[2][1]), (double)(av[3][1])/(double)(cv[3][1]), (double)(av[4][1])/(double)(cv[4][1]));
           //printf("%d %f %f %f %f %f\n", mc/100, paramx[q][5], paramx[q][6], paramx[q][7], paramx[q][8], paramx[q][9]);
           
           //logL = log_likelihood_print(net, D, paramx[q], freq, SN, N, Tobs);
           //x = log_likelihood_max(net, D, paramx[q], freq, SN, N, Tobs, min[5], max[5], 1);
           
               for(k=1; k <= NC; k++)
               {
                   q = who[k];
                   Mchirp = exp(paramx[q][0])/MSUN_SI;
                   printf("%f %f %d ", logLx[q], Mchirp, q);
               }
               printf("\n");
           
       }
        
        
        
    }
    
    if(lmax == 1) fclose(chainx);
    
    
    // record maxL parameters
    like = fopen("maxlike_intrinsic.dat", "w");
    fprintf(like,"%f ", maxL);
    for(i = 0; i < NX+3*net->Nifo; i++) fprintf(like,"%e ", pmax[i]);
    fprintf(like,"\n");
    fclose(like);
    
    free_int_matrix(av,5);
    free_int_matrix(cv,5);
    
    free_int_vector(m);
    
    free_double_matrix(ejump,NC+1);
    free_double_tensor(fish,NC+1,NX);
    free_double_tensor(evec,NC+1,NX);
    
    free_double_vector(pmax);
    free_double_matrix(paramy,NC+1);
    free_double_vector(logLx);
    
    
    free(max);
    free(min);
    
    
}
    
void updatei(int k, struct Net *net, int lmax, double *logLx, double **paramx, double **paramy, double *min, double *max, double *Fscale, int *who, double *heat, double ***history, double ***global, RealVector *freq, double **D, double **SN, double **ejump, double ***evec, int N, double Tobs, int **cv, int **av, gsl_rng *r)
{
        int q, i, j;
        double qx, qy, a, b, c, d;
        double alpha, beta, DL, pDx, pDy, H;
        double x, pAx, pAy, psx, psy;
        double logLy, eta, leta, pMcMx, pMcMy;
        int typ, flag;
        
        q = who[k];
        
        qx = qy = 0.0;    // log proposal densities
        
        alpha = gsl_rng_uniform(r);
    
        if(lmax == 0)
        {
        a = 0.5;
        b = 0.2;
        c = 0.1;
        d = 0.1;
        }
        else
        {
            a = 0.7;
            b = 0.5;
            c = 0.45;
            d = 0.1;
        }
        
        if(logLx[q] < 20.0)
        {
            a = 1.0;
            b = 1.0;
            c = 0.9;
            d = 1.0;
        }
        
        
        if(alpha > a) // fisher jump
        {
            typ = 1;
            
            // pick an eigendirection to jump in
            beta = gsl_rng_uniform(r);
            i = (int)(beta*NX);
            
            // draw the jump size
            beta = sqrt(heat[k])*ejump[q][i]*gsl_ran_gaussian(r,1.0);
            
            for(j = 0; j < NX; j++) paramy[q][j] = paramx[q][j]+beta*evec[q][i][j];
            
        }
        else if (alpha > b) // differential evolution
        {
            typ = 2;
            
            // the history file is kept for each temperature
            de_jump(paramx[q], paramy[q], history[k], NH, NX, r);
        }
        else if (alpha > c) // jiggle (most useful early when Fisher not effective)
        {
            typ = 3;

            beta = 0.01*pow(10.0, -floor(3.0*gsl_rng_uniform(r)))*sqrt(heat[k]);
            for(j = 0; j < NX; j++) paramy[q][j] = paramx[q][j]+beta*gsl_ran_gaussian(r,1.0);
        }
        else   // global
        {
            typ = 4;
            
            qx = qy = 1.0;
            
            do
            {
            for(j = 0; j < 2; j++) paramy[q][j] = min[j] + (max[j]-min[j])*gsl_rng_uniform(r);
            eta = exp((5.0/3.0)*(paramy[q][0]-paramy[q][1]));
            }while(eta > 0.25 || eta < etamin);
            
            for(j = 2; j < NX; j++) paramy[q][j] = min[j] + (max[j]-min[j])*gsl_rng_uniform(r);
            
            for (j = 1;  j < net->Nifo; ++j)
            {
                paramy[q][NX+(j-1)*3] = gsl_rng_uniform(r)*PI; // phase offset
                paramy[q][NX+(j-1)*3+1] = net->delays[0][j]*(-1.0+2.0*gsl_rng_uniform(r)); // time offset
                paramy[q][NX+(j-1)*3+2] = exp(-2.3+4.6*gsl_rng_uniform(r)); // amplitude ratio
            }
            

            //qy = globe(global, max, min, Tobs, paramy[q], N, r);
            //qx = globeden(global, max, min, Tobs, paramx[q], N);
        }
        
        cv[typ][k]++;
        
        flag = 0;
        
        // check intrinsic
        for(i = 0; i < 4; i++)
        {
            if(paramy[q][i] > max[i] || paramy[q][i] < min[i]) flag = 1;
        }
        
        // eta cannot exceed 0.25
        leta = (5.0/3.0)*(paramy[q][0]-paramy[q][1]);
        if(leta > log(0.25)) flag = 1;
        
        // Jacobian that makes the prior flat in m1, m2.
        if(flag == 0)
        {
            eta = exp(leta);
            pMcMy = 2.0*paramy[q][1]+leta-0.5*log(1.0-4.0*eta);
            
            leta = (5.0/3.0)*(paramx[q][0]-paramx[q][1]);
            eta = exp(leta);
            pMcMx = 2.0*paramx[q][1]+leta-0.5*log(1.0-4.0*eta);
            
            // cap on mass ratio - for PhenomD this is 1:18
            if(eta < etamin) flag = 1;
        }
        
        
        logLy = -1.0e20;
        if(flag == 0)
        {
            logLy = 0.0;
            if(lhold == 0)
            {
                if(lmax == 0)
                {
                    logLy = log_likelihood_intrinsic(net, D, paramy[q], freq, SN, N, Tobs);
                }
                else
                {
                    // GlitchSafe maximization
                    beta = gsl_rng_uniform(r);
                    
                    if(beta < d)
                    {
                    // time, amplitude and phase maximized (slow)
                    logLy = log_likelihood_max(net, D, paramy[q], freq, SN, N, Tobs, min[5], max[5], 0);
                    }
                    else
                    {
                    // amplitude and phase maximized (faster)
                    logLy = log_likelihood_APmax(net, D, paramy[q], freq, SN, N, Tobs);
                    }
                    
                    // maximization can throw distance out of allowed range
                    DL = Fscale[q]*exp(paramy[q][6])/PC_SI;
                    if(DL > DLmax) flag = 1;
                    if(DL < DLmin) flag = 1;
                }
            }
        }
        
        // now check the extrinsic parameters
        
        if(paramy[q][5] > max[5] || paramy[q][5] < min[5]) flag = 1;
        DL = Fscale[q]*exp(paramy[q][6])/PC_SI;
        if(DL > DLmax) flag = 1;
        if(DL < DLmin) flag = 1;
        
        if(flag == 1) logLy = 0.0;
        
        
        
        //  DLx = Fscale[q]*exp(paramx[q][6]);
        //  DLy = Fscale[q]*exp(paramy[q][6]);
        // since the Fscale and sky parameters are held fixed in the extrinsic MCMC can
        // use the expresion below for the DL^2 prior
        
        
        // variable in MCMC is x=logD, so p(x) = dD/dx p(D) = D p(D) = D^3
        
        pDx = 3.0*paramx[q][6];   // uniform in volume prior
        pDy = 3.0*paramy[q][6];   // uniform in volume prior
    
        /*
        x = log10(paramx[q][NX+2]); // log10 amplitude ratio
        pAx = -(x*x)/0.0648;   // amplitude ratio prior
    
        x = log10(paramy[q][NX+2]); // log10 amplitude ratio
        pAy = -(x*x)/0.0648;   // amplitude ratio prior
        */
    
        pAx = 0.0;
        pAy = 0.0;
    

        
        if(lmax == 1)  // no need to maintain balance if in search mode
        {
            pDx = pDy = 0.0;   // otherwise the distance prior fights with the likelihood during search
            qy = qx = 0.0;
        }
        
        
        H = (logLy-logLx[q])/heat[k] + pMcMy + pDy - qy - pDx - pMcMx + qx;
    
         // prior to mimic IMRPhenomP priors
        if(Pmimic == 1 && flag == 0)
         {
         psx = log(-log((fabs(paramx[q][2])+1.0e-8)/smax))+log(-log((fabs(paramx[q][3])+1.0e-8)/smax));
         psy = log(-log((fabs(paramy[q][2])+1.0e-8)/smax))+log(-log((fabs(paramy[q][3])+1.0e-8)/smax));
         H += (psy - psx);
         }
    
        //H = (logLy-logLx[q])/heat[k] + pMcMy + pDy + pAy - qy - pDx - pAx - pMcMx + qx;
        
        
        alpha = log(gsl_rng_uniform(r));
        
        if(H > alpha && flag == 0)
        {
            // copy over new state if accepted
            logLx[q] = logLy;
            for(i = 0; i < NX+3*net->Nifo; i++) paramx[q][i] = paramy[q][i];
            av[typ][k]++;
        }
        
}


void MCMC_all(struct Net *net, struct PSD *psd, int *mxc, int M, FILE *chain, double **paramx, int *who, double *heat, double ***history, double ***global, RealVector *freq, double **D, double **SN, double **SM, int N, double Tobs, double ttrig, gsl_rng *r)
{
    int i, j, q, k, mc, flag, wp;
    int id1, id2;
    int *m;
    double *logLx, *logPx, logLy, logL, x;
    double **paramy;
    double *pref;
    double *max, *min;
    double *href;
    double *pml;
    double alpha, beta, H;
    double Ap, Ac, psi, sindelta;
    double Fcross, Fplus, Fp;
    double Mchirp, Mtot, M1, M2, ciota, tc;
    double eta, dm, m1, m2, chieff, DL;
    double pMcMx, pMcMy;
    double a, b, c, d;
    double sqh;
    double Lmax, mch;
    double **tp;
    double logpx, logpy;
    double Jack;
    double dA, dTD, dP, z, y;
    double leta;
    int bin;
    int c1, c2, c3, c4;
    int a1, a2, a3, a4;
    int typ, hold;
    int scount, sacc, mcount;
    double *Fscale;
    double pDx, pDy, qx, qy;
    double *rtheta, *rphi;
    double **rhox;
    int MR;
    double mapL, lmapx, lmapy;
    double *pmax;
    double *lch, *lht;
 
    double *dtimes;
    
    FILE *like;
    FILE *ring;
    FILE *sky;
    
    double ***fish, ***evec;
    double **ejump;
    
    int **av, **cv;
    
    char filename[1024];
    
    clock_t start, end;
    double cpu_time_used;
    
    struct Het *het  = malloc(sizeof(struct Het));
    
    av = int_matrix(5,NC+1);
    cv = int_matrix(5,NC+1);
    
    // number of points printed for each sky ring
    MR = 1000;
    rtheta = double_vector(MR+1);
    rphi = double_vector(MR+1);
    
    m = int_vector(NC+1);
    
    tp = double_matrix(net->Nifo,2);
    
    dtimes = double_vector(6);
    rhox = double_matrix(NC+1,net->Nifo);
    ejump = double_matrix(NC+1,NP);
    fish = double_tensor(NC+1,NP,NP);
    evec = double_tensor(NC+1,NP,NP);
    
    lch = double_vector(NC+1);
    lht = double_vector(NC+1);
    pmax = double_vector(NP);
    paramy = double_matrix(NC+1,NP);
    logLx = double_vector(NC+1);
    logPx = double_vector(NC+1);
    
    // prior boundaries
    max = (double*)malloc(sizeof(double)* (NP));
    min = (double*)malloc(sizeof(double)* (NP));
    
    // [0] log(Mc) [1] log(Mt) [2] chi1 [3] chi2 [4] phi0 [5] tp [6] log(DL) [7] alpha [8] sindelta [9] psi [10] ciota

    max[0] = log(mcmax*MSUN_SI);  // Mc is at most 0.435 times Mt
    max[1] = log(mtmax*MSUN_SI);  // Mt
    max[2] = smax;  // chi1
    max[3] = smax;  // chi2
    max[4] = TPI;  // phi0
    max[5] = net->tmax; // peak time  (Trigger time is Tobs-2)
    max[6] = log(DLmax * PC_SI); //distance
    max[7] = TPI;  // alpha
    max[8] = 1.0;  // sindelta
    max[9] = PI; // psi
    max[10] = 1.0; // ciota
    
    
    min[0] = log(mcmin*MSUN_SI);  // Mc
    min[1] = log(mtmin*MSUN_SI);  // Mt
    min[2] = -smax;  // chi1
    min[3] = -smax;  // chi2
    min[4] = 0.0;  // phi0
    min[5] = net->tmin; // peak time (Trigger time is Tobs-2)
    min[6] = log(DLmin * PC_SI); //distance
    min[7] = 0.0;  // alpha
    min[8] = -1.0;  // sindelta
    min[9] = 0.0; // psi
    min[10] = -1.0; // ciota
    
    // initialize the likelihoods
    for(k = 1; k <= NC; k++)
    {
        logLx[k] = 0.0;
        if(lhold == 0)  logLx[k] = log_likelihood_full(net, D, paramx[k], freq, SN, rhox[k], N, Tobs);
        printf("%f ", logLx[k]);
        for(i = 0; i < net->Nifo; i++) printf("%f ", rhox[k][i]);
        printf("\n");
    }
    
    // find max likelihood waveform
    x = 0.0;
    for(k = 1; k <= NC; k++)
    {
      if(logLx[k] > x)
      {
          x = logLx[k];
          i = k;
      }
    }
    
    for(k = 0; k < NP; k++) pmax[k] = paramx[i][k];
    
    // set up heterodyne using max L waveform
     het_space(net, het, paramx[i], min, max, freq, SN, N, Tobs);
     heterodyne(net, het, D, paramx[i], freq, SN, N, Tobs);
    
    
     for(k = 1; k <= NC; k++)
     {
        //Fisher_All(net, fish[k], paramx[k], freq, SN, N, Tobs);
        Fisher_Het(net, het, fish[k], paramx[k], Tobs);
        FisherEvec(fish[k], ejump[k], evec[k], NP);
        efix(net, het, paramx[k], min, max, freq, SN, N, Tobs, ejump[k], evec[k], 1.0, 1);
     }
    
    
    start = clock();
    for(k = 0; k < 1000; k++) x = log_likelihood_full(net, D, paramx[1], freq, SN, rhox[1], N, Tobs);
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("Regular took %e s\n", cpu_time_used/1000.0);
    
    start = clock();
    for(k = 0; k < 1000; k++) x = log_likelihood_het(net, het, paramx[1], Tobs);
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("Heterodyned took %e s\n", cpu_time_used/1000.0);
    
    // reset likelihoods (slight differences in starting values due to Tukey window etc)
    for(k = 1; k <= NC; k++) logLx[k] = log_likelihood_het(net, het, paramx[k], Tobs);
    
    // re-initialize history
    for(i=0; i < NH; i++)
      {
        for(k=1; k <= NC; k++)
        {
            q = (int)(gsl_rng_uniform(r)*(double)(NC)) + 1; // use values from all temperatures
            for(j=0; j< NP; j++) history[k][i][j] = paramx[q][j];
        }
     }
    
    sky = fopen("skychain.dat","w");
    
    
    printf("Full MCMC\n");
    
    for(j = 1; j <= 4; j++)
    {
        for(k=1; k <= NC; k++)
        {
            av[j][k] = 0;
            cv[j][k] = 0;
        }
    }
    
    scount = 1;
    sacc = 0;
    mcount = 1;
    
    for(k=1; k <= NC; k++) m[k] = 0;
    
    mapL = -1.0e20;
    
    // print 100 waveforms
    wp = M/100;
    
    for(mc = -M/8; mc <= M; mc++)
    {
        // update the entire heterodyne
        if(mc%55000 == 0)
        {
          freehet(net,het);
          het_space(net, het, pmax, min, max, freq, SN, N, Tobs);
          heterodyne(net, het, D, pmax, freq, SN, N, Tobs);
          for(k = 1; k <= NC; k++) logLx[k] = log_likelihood_het(net, het, paramx[k], Tobs);
        }
        
    
        // (marginalize over PSD)
        if(noisemarg == 1 && mc%100==0)
        {
         
         // pick a PSD
         k = (int)(gsl_rng_uniform(r)*(double)(psd->Nsample));
        for (j = 0; j < net->Nifo; ++j)
        {
         makespec(psd->Nspline[j], psd->Nlines[j], psd->ffit[j], psd->xspline[j][k], psd->linef[j][k], psd->lineh[j][k],  psd->lineQ[j][k], psd->linew[j][k], psd->deltafmax[j][k], SN[j], Tobs, N);
        }
            
            // here I'm using the original parameters to define the heterodyne waveform
            //  we don't need to re-set the heterodyne if just updating the PSD
            heterodyne(net, het, D, het->pref, freq, SN, N, Tobs);
            
         for(k = 1; k <= NC; k++) logLx[k] = log_likelihood_het(net, het, paramx[k], Tobs);
            
        }
        
        if(mc > 0 && mc%10000==0 ) fflush(chain);
        
        if(mc > 0 && mc%wp ==0) printwaveall(net, N, freq, paramx[who[1]], SN, Tobs, ttrig, mc/wp);
        
        if(mc%50000==0 )
        {
            // update the Fisher matrices
            #pragma omp parallel for
            for(k = 1; k <= NC; k++)
            {
                //Fisher_All(net, fish[k], paramx[k], freq, SN, N, Tobs);
                Fisher_Het(net, het, fish[k], paramx[k], Tobs);
                FisherEvec(fish[k], ejump[k], evec[k], NP);
                efix(net, het, paramx[k], min, max, freq, SN, N, Tobs, ejump[k], evec[k], 1.0, 1);
            }
        }
        
        alpha = gsl_rng_uniform(r);
        
        if((NC > 1) && (alpha < 0.2))  // decide if we are doing a MCMC update of all the chains or a PT swap
        {
            
            // chain swap
            scount++;
            
            alpha = (double)(NC-1)*gsl_rng_uniform(r);
            j = (int)(alpha) + 1;
            beta = exp((logLx[who[j]]-logLx[who[j+1]])/heat[j+1] - (logLx[who[j]]-logLx[who[j+1]])/heat[j]);
            alpha = gsl_rng_uniform(r);
            if(beta > alpha)
            {
                hold = who[j];
                who[j] = who[j+1];
                who[j+1] = hold;
                sacc++;
            }
            
        }
        else      // MCMC update
        {
            
            mcount++;
            
            for(j = 1; j <= NC; j++)
            {
                for(i = 0; i < NP; i++)
                {
                    paramy[j][i] = paramx[j][i];
                }
            }
            
            #pragma omp parallel for
            for(k=1; k <= NC; k++)
            {
                update(k, net, het,  lch, lht, logLx, logPx, rhox, paramx, paramy, min, max, who, heat, history, global, freq, D, SN, SM, ejump, evec, N, Tobs, cv, av, rvec[k]);
            }
            
            // add to the history file
            if(mc%20 == 0)
            {
                for(k=1; k <= NC; k++)
                {
                    q = who[k];
                    i = m[k]%1000;
                    // the history file is kept for each temperature
                    for(j=0; j< NP; j++) history[k][i][j] = paramx[q][j];
                    m[k]++;
                }
            }

        // update MAP parameters
         for(j=1; j<=NCC; j++)
         {
          if(logPx[who[j]] > mapL)
          {
            mapL = logPx[who[j]];
            for(i = 0; i < NP; i++) pmax[i] = paramx[who[j]][i];
          }
         }
        
        }
    
        
        if(mc%10 == 0)
        {
            for(j=1; j<=NCC; j++)
            {
            q = who[j];
            Mchirp = exp(paramx[q][0])/MSUN_SI;
            Mtot = exp(paramx[q][1])/MSUN_SI;
            eta = pow((Mchirp/Mtot), (5.0/3.0));
            dm = sqrt(1.0-4.0*eta);
            m1 = Mtot*(1.0+dm)/2.0;
            m2 = Mtot*(1.0-dm)/2.0;
            chieff = (m1*paramx[q][2]+m2*paramx[q][3])/Mtot;
           
            DL = exp(paramx[q][6])/(1.0e6*PC_SI);  // includes conversion to Mpc
            z = z_DL(DL);
            
            // [0] log(Mc) [1] log(Mt) [2] chi1 [3] chi2 [4] phi0 [5] tp [6] log(DL) [7] alpha [8] sindelta [9] psi [10] ciota
            
            // counter, log likelihood, chirp mass (DF), total mass (DF), effective spin, geocenter orbital phase, geocenter arrival time, distance, RA , sine of DEC,
            // polarization angle, cos inclination, redshift, chirp mass (SF), total mass (SF), m1 (SF), m2 (SF), m1 (DF), m2 (DF), arrival time at each detector, SNR at each detector
            
            // Meger time is printed relative to trigger time, which is at Tobs-2.
            
            //wavemax(net, N, tp, freq, paramx, who, SN, Tobs);
            
            Times(paramx[q][7], paramx[q][8], net->GMST, dtimes);
            
            tc = paramx[q][5]-(Tobs-2.0);
                
               // [7] alpha [8] sindelta
            fprintf(sky,"%d %e %f %f %f\n", mxc[2], logLx[q], paramx[q][7], paramx[q][8], DL);
            
            fprintf(chain,"%d %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e ", mxc[2], logLx[q], Mchirp, Mtot, chieff, paramx[q][4], \
                    tc, DL, paramx[q][7], paramx[q][8], \
                    paramx[q][9], paramx[q][10], z, Mchirp/(1.0+z), Mtot/(1.0+z), m1/(1.0+z), m2/(1.0+z), m1, m2);
            for(i = 0; i < net->Nifo; i++) fprintf(chain, "%e ", tc+dtimes[net->labels[i]]);
            // for(i = 0; i < net->Nifo; i++) fprintf(chain, "%e ", rhox[q][i]);
           // for(i = 0; i < net->Nifo; i++) fprintf(chain, "%e %e ", tp[net->labels[i]][0], tp[net->labels[i]][1]);
            fprintf(chain,"\n");
 
            
            mxc[2] += 1;
            }
            
        }
        
        if(mc%100 == 0)
        {
            q = who[1];
            Mchirp = exp(paramx[q][0])/MSUN_SI;
            Mtot = exp(paramx[q][1])/MSUN_SI;
            
            printf("%d %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", mc, logLx[q], Mchirp, Mtot, paramx[q][4], paramx[q][5], paramx[q][6], paramx[q][7], paramx[q][8], paramx[q][9], \
                   (double)(sacc)/(double)(scount), (double)(av[1][1])/(double)(cv[1][1]), (double)(av[2][1])/(double)(cv[2][1]), (double)(av[3][1])/(double)(cv[3][1]), (double)(av[4][1])/(double)(cv[4][1]));
            //printf("%d %f %f %f %f %f\n", mc/100, paramx[q][5], paramx[q][6], paramx[q][7], paramx[q][8], paramx[q][9]);
            
            printf("%d ", mc);
            for(k=1; k <= NC; k++)
            {
                q = who[k];
                printf("%f %d ", logLx[q], q);
            }
            printf("\n");
            
        }
        
        
        
    }
    
    fclose(sky);
    
    x = log_likelihood_full(net, D, pmax, freq, SN, rhox[1], N, Tobs);
    
    // record MAP parameters
    like = fopen("maplike_all.dat", "w");
    fprintf(like,"%f ", x);
    for(i = 0; i < NP; i++) fprintf(like,"%e ", pmax[i]);
    fprintf(like,"\n");
    fclose(like);
    
    
    printf("likelihood at MAP %f\n", x);
    
    for (i=0; i< net->Nifo; i++)
    {
        for (j=i+1; j< net->Nifo; j++)
        {
            printf("detector %d-%d  SNR ratio %f\n", net->labels[i], net->labels[j], rhox[1][i]/rhox[1][j]);
        }
    }
    
    detector_shifts(net, pmax);
    
    // record MAP waveform
    printwaveall(net, N, freq, pmax, SN, Tobs, ttrig, -1);
    
    // record MAP f-tb track
    double *tm, *fr;
    tm = double_vector(N);
    fr = double_vector(N);
    for (i = 0; i < N; i++) tm[i] = Tobs/(double)(N)*(double)(i);
    f_of_t(pmax, tm, fr, N);
    Times(pmax[7], pmax[8], net->GMST, dtimes);
    for(j = 0; j < net->Nifo; j++)
    {
      sprintf(filename, "tftrack_%d.dat", net->labels[j]);
      like = fopen(filename,"w");
      for (i = 0; i < N; i++) fprintf(like,"%e %e\n",  tm[i]-Tobs+2.0+dtimes[net->labels[j]], fr[i]);
      fclose(like);
    }
    free_double_vector(tm);
    free_double_vector(fr);

    free_int_matrix(av,5);
    free_int_matrix(cv,5);
    
    free_int_vector(m);
    
    free_double_vector(rtheta);
    free_double_vector(rphi);
    
    free_double_matrix(tp,net->Nifo);
    
    free_double_matrix(rhox,NC+1);
    free_double_matrix(ejump,NC+1);
    free_double_tensor(fish,NC+1,NP);
    free_double_tensor(evec,NC+1,NP);
    
    free_double_matrix(paramy,NC+1);
    free_double_vector(logLx);
    free_double_vector(logPx);
    
    free_double_vector(lch);
    free_double_vector(lht);
    
    
    free(max);
    free(min);
    
    
}

void update(int k, struct Net *net, struct Het *het, double *lch, double *lht, double *logLx, double *logPx, double **rhox, double **paramx, double **paramy, double *min, double *max, int *who, double *heat, double ***history, double ***global, RealVector *freq, double **D, double **SN, double **SM, double **ejump, double ***evec, int N, double Tobs, int **cv, int **av, gsl_rng *r)
    {
        int q, i, j;
        double qx, qy, a, b, c;
        double alpha, beta, DL, pDx, pDy, H, Jack;
        double logLy, eta, leta, pMcMx, pMcMy;
        double psx, psy, x;
        double *rhoy;
        int typ, flag, id1, id2;
        
        rhoy = (double*)malloc(sizeof(double)* (net->Nifo));
        
        a = 0.5;
        b = 0.3;
        c = 0.2;
        
        /*
        a = 0.4;
        b = 0.1;
        c = 0.0;
        */
        
        q = who[k];
        
        Jack = 0.0;
        
        qx = qy = 0.0;    // log proposal densities
        
        alpha = gsl_rng_uniform(r);
        
        if(alpha > a) // fisher jump
        {
            typ = 1;
            // pick an eigendirection to jump in
            beta = gsl_rng_uniform(r);
            i = (int)(beta*NP);
            
            // draw the jump size
            beta = sqrt(heat[k])*ejump[q][i]*gsl_ran_gaussian(r,1.0);
            
            for(j = 0; j < NP; j++) paramy[q][j] = paramx[q][j]+beta*evec[q][i][j];
            
        }
        else if (alpha > b) // differential evolution
        {
            typ = 2;
            // the history file is kept for each temperature
            de_jump(paramx[q], paramy[q], history[k], NH, NP, r);
        }
        else if (alpha > c) // jiggle (most useful early when Fisher not effective)
        {
            typ = 3;
            beta = 0.01*pow(10.0, -floor(3.0*gsl_rng_uniform(r)))*sqrt(heat[k]);
            for(j = 0; j < NP; j++) paramy[q][j] = paramx[q][j]+beta*gsl_ran_gaussian(r,1.0);
        }
        else   // sky ring
        {
            typ = 4;
            
            id1 = 0;
            id2 = 1;
            
            // Pick a pair of interferometers to define sky ring
            if(net->Nifo > 2)
            {
                id1 = (int)((double)(net->Nifo)*gsl_rng_uniform(r));
                do
                {
                    id2 = (int)((double)(net->Nifo)*gsl_rng_uniform(r));
                }while(id1==id2);
            }
            
            id1 = net->labels[id1];
            id2 = net->labels[id2];
            
            Ring_all(paramx[q], paramy[q], id1, id2, net->GMST, r);
            
            skymap_all(paramx[q], paramy[q], net->GMST, id1, id2);
            Jack = log(skydensity_all(paramx[q], paramy[q], net->GMST, id1, id2));
            
            // The mapping only overs half the phi, psi space. Can cover it all by radomly shifting both by a half period
            beta = gsl_rng_uniform(r);
            if(beta > 0.5)
            {
                paramy[q][4] += PI;
                paramy[q][9] += PI/2.0;
            }
            
        }
        
        cv[typ][k]++;
        
        // [0] log(Mc) [1] log(Mt) [2] chi1 [3] chi2 [4] phi0 [5] tp [6] log(DL) [7] alpha [8] sindelta [9] psi [10] ciota
        
        // re-map angular parameters to their proper range
        while(paramy[q][4] > TPI)   paramy[q][4] -= TPI;
        while(paramy[q][4] < 0.0)  paramy[q][4] += TPI;
        while(paramy[q][7] > TPI)  paramy[q][7] -= TPI;
        while(paramy[q][7] < 0.0)  paramy[q][7] += TPI;
        while(paramy[q][9] > PI)   paramy[q][9] -= PI;
        while(paramy[q][9] < 0.0)  paramy[q][9] += PI;
        
        
        // check proposed values are in prior range
        flag = 0;
        for(i = 0; i < NP; i++)
        {
            if(paramy[q][i] > max[i] || paramy[q][i] < min[i]) flag = 1;
        }
        
        // eta cannot exceed 0.25
        leta = (5.0/3.0)*(paramy[q][0]-paramy[q][1]);
        if(leta > log(0.25)) flag = 1;
        
        // Jacobian that makes the prior flat in m1, m2.
        if(flag == 0)
        {
            eta = exp(leta);
            pMcMy = 2.0*paramy[q][1]+leta-0.5*log(1.0-4.0*eta);
            // cap on mass ratio - for PhenomD this is 1:18
            if(eta < etamin) flag = 1;
            
            leta = (5.0/3.0)*(paramx[q][0]-paramx[q][1]);
            eta = exp(leta);
            pMcMx = 2.0*paramx[q][1]+leta-0.5*log(1.0-4.0*eta);

        }
        
        logLy = -1.0e20;
        if(flag == 0)
        {
            logLy = 0.0;
            //if(lhold == 0) logLy = log_likelihood_full(net, D, paramy[q], freq, SN, rhoy, N, Tobs);
            //x = log_likelihood_test(net, het, D, paramy[q], freq, SN, SM, N, Tobs);
            if(lhold == 0) logLy = log_likelihood_het(net, het, paramy[q], Tobs);
            //printf("%f %f %f\n", x, logLy, x-logLy);
        }
        
        // variable in MCMC is x=logD, so p(x) = dD/dx p(D) = D p(D) = D^3
        
        pDx = 3.0*paramx[q][6];   // unifoform in volume prior
        pDy = 3.0*paramy[q][6];   // unifoform in volume prior
        
        H = Jack + (logLy-logLx[q])/heat[k] + pMcMy + pDy - qy - pDx - pMcMx + qx;
        
         // prior to mimic IMRPhenomP priors
        if(Pmimic == 1 && flag == 0)
         {
         psx = log(-log((fabs(paramx[q][2])+1.0e-8)/smax))+log(-log((fabs(paramx[q][3])+1.0e-8)/smax));
         psy = log(-log((fabs(paramy[q][2])+1.0e-8)/smax))+log(-log((fabs(paramy[q][3])+1.0e-8)/smax));
         H += (psy - psx);
         }
        
        alpha = log(gsl_rng_uniform(r));
        
        if(H > alpha && flag == 0)
        {
            // copy over new state if accepted
            logLx[q] = logLy;
            logPx[q] = logLy + pMcMy + pDy;  // posterior density
            for(i = 0; i < net->Nifo; i++) rhox[q][i] = rhoy[i];
            for(i = 0; i < NP; i++) paramx[q][i] = paramy[q][i];
            av[typ][k]++;
        }
        
        free(rhoy);
        
}
    
void de_jump(double *paramsx, double *paramsy, double **history, int m, int d, gsl_rng *r)
{
    int i, j, k;
    double alpha, beta;
    
    // pick two points from the history
    i = (int)((double)(m)*gsl_rng_uniform(r));
    do
    {
        j = (int)((double)(m)*gsl_rng_uniform(r));
    } while(i==j);
    
    alpha = 1.0;
    beta = gsl_rng_uniform(r);
    if(beta < 0.9) alpha = gsl_ran_gaussian(r,0.5);
    
    for(k=0; k< d; k++) paramsy[k] = paramsx[k]+alpha*(history[i][k]-history[j][k]);
    
}

void fisher_skyproposal(gsl_rng * r, double **skyvecs, double *skyevals, double *jump)
{
    int i, j, k;
    int a;
    double scale, x;

    
    // pick eigenvector to jump along
    a = (int)(gsl_rng_uniform(r)*(double)(NS));
    
    // set size of jump
    scale = gsl_ran_gaussian(r,1.0)*skyevals[a];
    
    
    //decompose eigenjumps into parameter directions
    
    for (j=0; j<NS; j++) jump[j] = scale*skyvecs[a][j];
    
}



void FisherEvec(double **fish, double *ej, double **ev, int d)
{
    int i, j, ec, sc;
    double x, maxc;
    
    // printf("Evec\n");
    
    /*
    printf("\n");
    for (i = 0 ; i < d ; i++)
    {
        for (j = 0 ; j < d ; j++)
        {
            printf("%e ", fish[i][j]);
        }
        printf("\n");
    } */
    
 
    
    ec = 0;
    for (i = 0 ; i < d ; i++) if(fabs(fish[i][i]) < 1.0e-16) ec = 1;
    
    if(ec == 0)
    {
        
       gsl_matrix *m = gsl_matrix_alloc (d, d);

    for (i = 0 ; i < d ; i++)
    {
        for (j = 0 ; j < d ; j++)
        {
            gsl_matrix_set(m, i, j, fish[i][j]);
        }
    }
    
    
    gsl_vector *eval = gsl_vector_alloc (d);
    gsl_matrix *evec = gsl_matrix_alloc (d, d);
    
    gsl_eigen_symmv_workspace * w =
    gsl_eigen_symmv_alloc (d);
    
    ec = gsl_eigen_symmv (m, eval, evec, w);
    
    gsl_eigen_symmv_free (w);
    
 
    sc = gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_ABS_ASC);
    
    for (i = 0; i < d; i++)
    {
        ej[i] = gsl_vector_get (eval, i);
        
        //printf("eigenvalue = %g\n", ej[i]);
        for (j = 0 ; j < d ; j++)
        {
            ev[i][j] = gsl_matrix_get(evec, j, i);
            //printf("%f ", ev[i][j]);
        }
        //printf("\n");
        
    }
    
    for (i = 0; i < d; i++)
    {
        // make sure no eigenvalue is too small
        if(ej[i] < 10.0) ej[i] = 10.0;
        // turn into 1-sigma jump amplitudes
        ej[i] = 1.0/sqrt(ej[i]);
        //printf("jump %d = %g\n", i, ej[i]);
    }
        
        gsl_matrix_free (m);
        gsl_vector_free (eval);
        gsl_matrix_free (evec);
        
    }
    else
    {
        for (i = 0; i < d; i++)
        {
            ej[i] = 10000.0;
            for (j = 0 ; j < d ; j++)
            {
                ev[i][j] = 0.0;
                if(i==j) ev[i][j] = 1.0;
            }
        }
        
    }
    
    

    
    return;
    
}

// Need to write a heterodyned Fisher
//void Fisher_Het(struct Net *net, struct Het *het, double **fish, double *params, double Tobs)

void Fisher_Full(struct Net *net, double **fish, int d, double *params, RealVector *freq, double **SN, int N, double Tobs)
{
    double *paramsP, *paramsM;
    double epsilon;
    double Scale;
    double *AH, *PH, **AHP, **AHM, **PHP, **PHM;
    int i, j, k, l;
    int imin, imax;
   
    //printf("Fisher ");
    
    for (i = 0; i < d; i++)
    {
        for (j = 0; j < d; j++) fish[j][i] = 0.0;
    }
    
    for (i = 0; i < d; i++) fish[i][i] = 1.0e4;
    
    
    if(fabs(params[2]) < 0.9999 &&  fabs(params[3]) < 0.9999)
    {
        
        
        imin = (int)(Tobs*fmin);
        imax = (int)(Tobs*fmax);
        if(imax > N/2) imax = N/2;

    
    epsilon = 1.0e-5;
    
    paramsP = (double*)malloc(sizeof(double)*(d));
    paramsM = (double*)malloc(sizeof(double)*(d));
    
    AH = (double*)malloc(sizeof(double)*(N/2));
    PH = (double*)malloc(sizeof(double)*(N/2));
    
    AHP = double_matrix(d,N/2);
    PHP = double_matrix(d,N/2);
    AHM = double_matrix(d,N/2);
    PHM = double_matrix(d,N/2);
 
    
    // Reference phase and amplitude
    AmpPhase(AH, PH, freq, params, N);
    // Note that since this Fisher matrix only varies intrinsic parameters, the different overall phasing and arrival time in each
    // detector is irrelvant. Only the overall amplitude differences matter, and that are taken care of by the Scale parameter below.
    
    for (i = 0 ; i < d ; i++)
    {
        
        for (k = 0 ; k < d ; k++)
        {
            paramsP[k] = params[k];
            paramsM[k] = params[k];
        }
        
        paramsP[i] += epsilon;
        paramsM[i] -= epsilon;
        
        AmpPhase(AHP[i], PHP[i], freq, paramsP, N);
        AmpPhase(AHM[i], PHM[i], freq, paramsM, N);
        
        // store the central difference in the Plus arrays
        for (k = 0 ; k < N/2 ; k++)
        {
            AHP[i][k] -= AHM[i][k];
            PHP[i][k] -= PHM[i][k];
            AHP[i][k] /= (2.0*epsilon);
            PHP[i][k] /= (2.0*epsilon);
        }
        
       // printf("%d ", i);
        
    }
    
    
     //printf("\n");
    
    for (i = 0; i < d; i++)
    {
        for (j = 0; j < d; j++) fish[j][i] = 0.0;
    }

    for (l = 0 ; l < net->Nifo ; l++)  // loop over detectors
    {
        if(l==0)
        {
            Scale = 4.0;
        }
        else
        {
            Scale = 4.0*params[(l-1)*3+NX+2]*params[(l-1)*3+NX+2];
        }
        
        for (i = 0 ; i < d ; i++)
        {
            for (j = i ; j < d ; j++)
            {
                for (k = imin ; k < imax ; k++) fish[i][j] += Scale*(AHP[i][k]*AHP[j][k]+AH[k]*AH[k]*PHP[i][k]*PHP[j][k])/SN[l][k];
            }
        }
        
    }
    
    
    /* fill in lower triangle */
    
    for (i = 0; i < d; i++)
    {
        for (j = i+1; j < d; j++)
            fish[j][i] = fish[i][j];
    }
    
    for (i = 0; i < d; i++) if(fish[i][i] < 1.0) fish[i][i] = 1.0;
    
    /*
    for (i = 0; i < d; i++)
    {
        for (j = 0; j < d; j++)
        {
            printf("%.3e ", fish[i][j]);
        }
        printf("\n");
    }
    */
    
    free(paramsP);
    free(paramsM);
    free(AH);
    free(PH);
    free_double_matrix(AHP,d);
    free_double_matrix(PHP,d);
    free_double_matrix(AHM,d);
    free_double_matrix(PHM,d);
    
   }

}






void Fisher_All(struct Net *net, double **fish, double *params, RealVector *freq, double **SN, int N, double Tobs)
{
    double *paramsP, *paramsM;
    double epsilon;
    double Scale;
    double **AH, **PH, ***AHP, ***AHM, ***PHP, ***PHM;
    int i, j, k, l, id;
    int imin, imax;
    
    //printf("Fisher ");
    
    for (i = 0; i < NP; i++)
    {
        for (j = 0; j < NP; j++) fish[j][i] = 0.0;
    }
    
    for (i = 0; i < NP; i++) fish[i][i] = 1.0e4;
    
    
    if(fabs(params[2]) < 0.9999 &&  fabs(params[3]) < 0.9999)
    {
        
        imin = (int)(Tobs*fmin);
        imax = (int)(Tobs*fmax);
        if(imax > N/2) imax = N/2;
        
        epsilon = 1.0e-5;
        
        paramsP = (double*)malloc(sizeof(double)*(NP));
        paramsM = (double*)malloc(sizeof(double)*(NP));
        
        AH = double_matrix(net->Nifo,N/2);
        PH = double_matrix(net->Nifo,N/2);
        AHP = double_tensor(NP,net->Nifo,N/2);
        PHP = double_tensor(NP,net->Nifo,N/2);
        AHM = double_tensor(NP,net->Nifo,N/2);
        PHM = double_tensor(NP,net->Nifo,N/2);
        
        // Reference phase and amplitude
        fullphaseamp(net, AH, PH, freq, params, N);
        
        for (i = 0 ; i < NP ; i++)
        {
            
            for (k = 0 ; k < NP ; k++)
            {
                paramsP[k] = params[k];
                paramsM[k] = params[k];
            }
            
            paramsP[i] += epsilon;
            paramsM[i] -= epsilon;
            
            fullphaseamp(net, AHP[i], PHP[i], freq, paramsP, N);
            fullphaseamp(net, AHM[i], PHM[i], freq, paramsM, N);
            
            // store the central difference in the Plus arrays
            for (id = 0 ; id < net->Nifo ; id++)  // loop over detectors
            {
            for (k = 1 ; k < N/2 ; k++)
            {
                AHP[i][id][k] -= AHM[i][id][k];
                PHP[i][id][k] -= PHM[i][id][k];
                AHP[i][id][k] /= (2.0*epsilon);
                PHP[i][id][k] /= (2.0*epsilon);
            }
            }
            
        }
        
        
        //printf("\n");
        
        for (i = 0; i < NP; i++)
        {
            for (j = 0; j < NP; j++) fish[j][i] = 0.0;
        }
        
        for (id = 0 ; id < net->Nifo ; id++)  // loop over detectors
        {
            
            for (i = 0 ; i < NP ; i++)
            {
                for (j = i ; j < NP ; j++)
                {
                    for (k = imin ; k < imax ; k++) fish[i][j] += 4.0*(AHP[i][id][k]*AHP[j][id][k]+AH[id][k]*AH[id][k]*PHP[i][id][k]*PHP[j][id][k])/SN[id][k];
                }
            }
            
        }
        
        
        /* fill in lower triangle */
        
        for (i = 0; i < NP; i++)
        {
            for (j = i+1; j < NP; j++)
                fish[j][i] = fish[i][j];
        }
        
        /*
         for (i = 0; i < NP; i++)
         {
         for (j = 0; j < NP; j++)
         {
         printf("%.3e ", fish[i][j]);
         }
         printf("\n");
         }
        printf("\n"); */
        
        free(paramsP);
        free(paramsM);
        free_double_matrix(AH,net->Nifo);
        free_double_matrix(PH,net->Nifo);
        free_double_tensor(AHP,NP,net->Nifo);
        free_double_tensor(PHP,NP,net->Nifo);
        free_double_tensor(AHM,NP,net->Nifo);
        free_double_tensor(PHM,NP,net->Nifo);
        
    }
    
}








void skyring(struct Net *net, double *params, double *sky, double *pall, double *SNRsq, RealVector *freq, double **D, double **SN, int N, double Tobs, gsl_rng * r)
{
    int i, j, k, id;
    int iref;
    double *delt;
    double Fp, Fc, Fs;
    double logL, logLmax;
    double *stry;
    
    
    iref = net->labels[0];
    stry = double_vector(NS);
    delt = double_vector(net->Nifo);
    
    
    for(id = 1; id< net->Nifo; id++)
    {
        delt[id] = -params[(id-1)*3+NX+1];
    }
    
    // sky parameter order
    //[0] alpha, [1] sin(delta) [2] psi [3] cosi [4] scale [5] dt [6] dphi
    
    // Scale changes the overall distance reference at geocenter
    // dt shifts the geocenter time
    // dphi shifts the geocenter phase
    
    stry[4] = 1.0;
    stry[5] = 0.0;
    stry[6] = 0.0;
    
        // find sky locations consistent with time delays
        
        ringfind(net, delt, stry, SNRsq, r);
    
        stry[2] = PI*gsl_rng_uniform(r);
        stry[3] = -1.0+2.0*gsl_rng_uniform(r);
            
            
        pmap(net, pall, params, stry);
    
    for(j = 0; j< NS; j++) sky[j] = stry[j];
    
    free_double_vector(delt);
    free_double_vector(stry);
    
}




void skystart(struct Net *net, double *params, double *sky, double *pall, double *SNRsq, RealVector *freq, double **D, double **SN, int N, double Tobs, gsl_rng * r)
{
    int i, j, k, id;
    int iref;
    double *delt;
    double Fp, Fc, Fs;
    double logL, logLmax;
    double *stry;
    double *rho;
    
    
    iref = net->labels[0];
    rho = double_vector(net->Nifo);
    stry = double_vector(NS);
    delt = double_vector(net->Nifo);

    
    for(id = 1; id< net->Nifo; id++)
    {
        delt[id] = -params[(id-1)*3+NX+1];
    }
    
    logLmax = -1.0e60;
    
    // sky parameter order
    //[0] alpha, [1] sin(delta) [2] psi [3] cosi [4] scale [5] dt [6] dphi
    
    // Scale changes the overall distance reference at geocenter
    // dt shifts the geocenter time
    // dphi shifts the geocenter phase
    
    stry[4] = 1.0;
    stry[5] = 0.0;
    stry[6] = 0.0;
    
    for(k = 0; k<100; k++)
    {
        // find sky locations consistent with time delays
       
        ringfind(net, delt, stry, SNRsq, r);
        
        
        for(i = 0; i<100; i++)
        {
            
            stry[2] = PI*gsl_rng_uniform(r);
            stry[3] = -1.0+2.0*gsl_rng_uniform(r);

            
            pmap(net, pall, params, stry);
           
            
            logL = log_likelihood_full(net, D, pall, freq, SN, rho, N, Tobs);
            
            if(logL > logLmax)
            {
                logLmax = logL;
                for(j = 0; j< NS; j++) sky[j] = stry[j];
            }
            
        }
    }
    
    free_double_vector(rho);
    free_double_vector(delt);
    free_double_vector(stry);
    
}




void skymcmc(struct Net *net, int MCX, int *mxc, FILE *chain, double **paramx, double **skyx, double **pallx, int *who, double *heat, double dtx, int nt, double *DD, double **WW, double ***DHc,  double ***DHs, double ***HH, double Tobs, gsl_rng * r)
{
    int i, j, k, q, ic, id1, id2;
    int scount, sacc, hold, mcount;
    int ac, rc, rca, clc, cla, PRcnt, POcnt;
    int sdx, Ax, mc;
    double alpha, beta;
    double Mchirp, Mtot, eta, dm, m1, m2, chieff, ciota;
    double qxy, qyx, Jack, phi;
    int rflag, cflag;
    double **sky, **skyy, **skyh;
    double x, y, z, DL, scale, logLy;
    double *logLx;
    double pAy, logH, pAx;
    double *param;
    double *dtimes, *dtimes2;
    double ***fishskyx, ***fishskyy;
    double ***skyvecsx, ***skyvecsy;
    double **skyevalsx, **skyevalsy;
    double Fp, Fc, Fs, ps;
    double *jump, *sqH;
    double ldetx, ldety;
    double scmax, scmin;
    double DLx, DLy;
    int fflag, fc, fac;
    int uflag, uc, uac;
    double Ap, Ac, Fcross, Fplus, lambda, lambda2, Fs2;
    double sindelta, psi;

    
    dtimes = (double*)malloc(sizeof(double)*5);
    dtimes2 = (double*)malloc(sizeof(double)*5);
    
    param = (double*)malloc(sizeof(double)*(NX+3*net->Nifo));
    
    // sky parameter order
    //[0] alpha, [1] sin(delta) [2] psi [3] ciota [4] scale [5] phi0 [6] dt
    

    
    // max and min of rescaling parameter
    scmin = 0.1;
    scmax = 10.0;
    
    sky = double_matrix(NC+1,NS);
    skyh = double_matrix(NC+1,NS);
    skyy = double_matrix(NC+1,NS);
    logLx = double_vector(NC+1);
    sqH = double_vector(NC+1);
    
    for(k = 1; k <= NC; k++)
    {
        for(i = 0; i < NS; i++) sky[k][i] = skyx[k][i];
        for(i = 0; i < NS; i++) skyh[k][i] = skyx[k][i];
    }
    
    for(k = 1; k <= NC; k++) sqH[k] = sqrt(heat[k]);
    
    ic = who[1];
    
    for(k = 1; k <= NC; k++) logLx[k] = skylike(net, skyx[k], DD, WW[k], DHc[k], DHs[k], dtx, nt, 0);
    
    
    fishskyx = double_tensor(NC+1,NS,NS);
    fishskyy = double_tensor(NC+1,NS,NS);
    skyvecsx = double_tensor(NC+1,NS,NS);
    skyvecsy = double_tensor(NC+1,NS,NS);
    skyevalsx = double_matrix(NC+1,NS);
    skyevalsy = double_matrix(NC+1,NS);
    jump = double_vector(NS);
    
    
    for(k = 1; k <= NC; k++)
    {
    fisher_matrix_fastsky(net, skyx[k], fishskyx[k], HH[k]);
    FisherEvec(fishskyx[k], skyevalsx[k], skyvecsx[k], NS);
    }
    
    
    printf("Extrinsic MCMC\n");
    
    
    ac = 0;
    rc = 1;
    rca = 0;
    fc = 1;
    uc = 1;
    
    clc = 0;
    cla = 0;
    fac = 0;
    uac = 0;
    
    PRcnt = 0;
    POcnt = 0;
    
    sdx = 0.0;
    Ax = 0.0;
    
    scount = 0;
    sacc = 0;
    mcount = 0;
    
    for(mc = 0; mc < MCX; mc++)
    {
        
        
        if(mc > 1 && mc%1000==0)
        {
            // update the Fisher matrices
            for(k = 1; k <= NC; k++)
            {
                fisher_matrix_fastsky(net, skyx[k], fishskyx[k], HH[k]);
                FisherEvec(fishskyx[k], skyevalsx[k], skyvecsx[k], NS);
            }
        }
        
        alpha = gsl_rng_uniform(r);
        
        if((NC > 1) && (alpha < 0.2))  // decide if we are doing a MCMC update of all the chains or a PT swap
        {
            
            // chain swap
            scount++;
            
            alpha = (double)(NC-1)*gsl_rng_uniform(r);
            j = (int)(alpha) + 1;
            beta = exp((logLx[who[j]]-logLx[who[j+1]])/heat[j+1] - (logLx[who[j]]-logLx[who[j+1]])/heat[j]);
            alpha = gsl_rng_uniform(r);
            if(beta > alpha)
            {
                hold = who[j];
                who[j] = who[j+1];
                who[j+1] = hold;
                sacc++;
            }
            
            
            
        }
        else      // MCMC update
        {
            
            mcount++;

            
            for(k = 1; k <= NC; k++)
            {
                for(i = 0; i < NS; i++) skyy[k][i] = skyx[k][i];
            }
            
            for(k=1; k <= NC; k++)
            {
                q = who[k];
                
                qxy = 0.0;
                qyx = 0.0;
                
                Jack = 0.0;
                
                rflag = 0;
                cflag = 0;
                fflag = 0;
                uflag = 0;
                
                
                alpha = gsl_rng_uniform(r);
                
                if(alpha > 0.8 && net->Nifo > 1)  // ring
                {
                    

                    id1 = 0;
                    id2 = 1;
                    
                    // Pick a pair of interferometers to define sky ring
                    if(net->Nifo > 2)
                    {
                      id1 = (int)((double)(net->Nifo)*gsl_rng_uniform(r));
                      do
                      {
                       id2 = (int)((double)(net->Nifo)*gsl_rng_uniform(r));
                      }while(id1==id2);
                    }
                    
                    // map these labels to actual detectors
                    id1 = net->labels[id1];
                    id2 = net->labels[id2];
                    
                    Ring(skyx[q], skyy[q], id1, id2, net->GMST, r);
                    qyx = 1.0;
                    qxy = 1.0;
                
                    
                    skymap(skyx[q], skyy[q], net->GMST, id1, id2, net->labels[0]);
                    
                    Jack = log(skydensity(skyx[q], skyy[q], net->GMST, id1, id2, net->labels[0]));
                    
                    // The mapping only overs half the phi, psi space. Can cover it all by radomly shifting both by a half period
                    x = gsl_rng_uniform(r);
                    if(x > 0.5)
                    {
                        skyy[q][5] += PI;
                        skyy[q][2] += PI/2.0;
                    }
                    
                    if(k==1) rc++;
                    
                    rflag = 1;
                    
                }
                else if (alpha > 0.2)  // Fisher matrix
                {
                    
                    if(k==1) fc++;
                    
                    fflag = 1;
                    
                    fisher_skyproposal(r, skyvecsx[q], skyevalsx[q], jump);
                    
                    for(i=0; i< NS; i++) skyy[q][i] = skyx[q][i]+sqH[k]*jump[i];
                    
                    // If the Fisher matrix was updated after each Fisher jump we would
                    // need these proposal densities. Since Fisher held fixed for blocks
                    // of iterations, we don't need the densities
                    // pfishxy = fisher_density(fishskyx, ldetx, skyx, skyy);
                    // fisher_matrix_fastsky(net, skyy, fishskyy, HH);
                    // fisher_skyvectors(fishskyy, skyvecsy, skyevalsy, &ldety);
                    //  pfishyx = fisher_density(fishskyy, ldety, skyy, skyx);
                }
                else  // jiggle (most useful early when Fisher not effective)
                {
                    
                    uflag = 1;
                    
                    if(k==1) uc++;
                    
                    beta = 0.01*pow(10.0, -floor(3.0*gsl_rng_uniform(r)))*sqH[k];
                    for(i = 0; i < NS-1; i++) skyy[q][i] = skyx[q][i]+beta*gsl_ran_gaussian(r,1.0);
                    skyy[q][6] = skyx[q][6]+0.01*beta*gsl_ran_gaussian(r,1.0);
                }
                
                
                if(skyy[q][0] > TPI) skyy[q][0] -= TPI;
                if(skyy[q][0] < 0.0) skyy[q][0] += TPI;
                if(skyy[q][2] > PI) skyy[q][2] -= PI;
                if(skyy[q][2] < 0.0) skyy[q][2] += PI;
                if(skyy[q][5] > TPI) skyy[q][5] -= TPI;
                if(skyy[q][5] < 0.0) skyy[q][5] += TPI;
                
                
                //[0] alpha, [1] sin(delta) [2] psi [3] cos(iota) [4] scale [5] phi0 [6] dt
                
 
                DLy = exp(pallx[q][6])/(skyy[q][4]*PC_SI);
                
                if(DLy < DLmin || DLy > DLmax || fabs(skyy[q][1]) > 1.0 || fabs(skyy[q][3]) > 1.0 || fabs(skyy[q][6]) > dtmax)
                {
                    logLy = -1.0e60;
                    
                    pAy = 0.0;
                    pAx = 0.0;
                    
                }
                else
                {
                    logLy = skylike(net, skyy[q], DD, WW[q], DHc[q], DHs[q], dtx, nt, 0);
                    
                    // Need a Jacobian a factor here since we sample uniformly in amplitude.
                    // Jacobian between D cos(theta) phi cos(iota) psi and A cos(theta) phi cos(iota) psi
                    // Since D = D_0/A, boils down to just D^2 |dD/dA| = D_0^3/A^4 = D^3/A
                
                    DLx = exp(pallx[q][6])/(skyx[q][4]);
                    pAx = 3.0*log(DLx)-log(skyx[q][4]);
  
                    DLy = exp(pallx[q][6])/(skyy[q][4]);
                    pAy = 3.0*log(DLy)- log(skyy[q][4]);
                    
                }
                
                
                logH = Jack + (logLy-logLx[q])/heat[k] +pAy-qyx-pAx+qxy;
                
                alpha = log(gsl_rng_uniform(r));
                
                
                if(logH > alpha)
                {
                    for(i=0; i< NS; i++) skyx[q][i] = skyy[q][i];
                    logLx[q] = logLy;
                    
                    if(k==1)
                    {
                        ac++;
                        if(rflag == 1) rca++;
                        if(fflag == 1) fac++;
                        if(uflag == 1) uac++;
                    }
                    
                }
                
            }  // ends loop over chains
            
        }  // ends choice of update
        
        /*
        if(mc%100 == 0)
        {
            ic = who[1];
            phi = skyx[ic][0];
            if(phi > TPI) phi -= TPI;
            if(phi < 0.0) phi += TPI;
            skyx[ic][0] = phi;
            Mchirp = exp(paramx[ic][0])/MSUN_SI;
            Mtot = exp(paramx[ic][1])/MSUN_SI;
            eta = pow((Mchirp/Mtot), (5.0/3.0));
            dm = sqrt(1.0-4.0*eta);
            m1 = Mtot*(1.0+dm)/2.0;
            m2 = Mtot*(1.0-dm)/2.0;
            chieff = (m1*paramx[ic][2]+m2*paramx[ic][3])/Mtot;

           
            // counter, log likelihood, chirp mass, total mass, effective spin,  phase shift , time shift, distance, RA , sine of DEC,
            // polarization angle, cos inclination
            
            
            DL = exp(pallx[ic][6])/(1.0e6*PC_SI*skyx[ic][4]);
            z = z_DL(DL);
            
            // Note that skyx[ic][5], skyx[ic][6], hold different quantities than what is printed by the other MCMCs
 
            fprintf(chain,"%d %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n", mxc[1], logLx[ic], Mchirp, Mtot, chieff, skyx[ic][5], \
                    skyx[ic][6], DL, skyx[ic][0], skyx[ic][1], \
                    skyx[ic][2], skyx[ic][3], z, Mchirp/(1.0+z), Mtot/(1.0+z), m1/(1.0+z), m2/(1.0+z), m1, m2);
            
        
                    
            
             mxc[1] += 1;
            
          } */
        
        if(mc%10000 == 0)
        {
         ic = who[1];
         DL = exp(pallx[ic][6])/(1.0e6*PC_SI*skyx[ic][4]);
         printf("%d %f %f %f %f %f %f\n", mc, logLx[ic], DL, skyx[ic][0], skyx[ic][1], skyx[ic][2], skyx[ic][3]);
        }
        
        if(mc > MCX/4 && mc%10 == 0)
        {
         ic = who[1];
         DL = exp(pallx[ic][6])/(1.0e6*PC_SI*skyx[ic][4]);
         fprintf(chain, "%d %f %f %f %f\n", mc, logLx[ic], skyx[ic][0], skyx[ic][1], DL);
        }

    }
    

    // update the amplitude, time and phase shifts between detectors in preparation for extrinsic updates
    for(k=1; k <= NC; k++) dshifts(net, skyx[k], paramx[k]);
    
    
    // sky  [0] alpha, [1] sin(delta) [2] psi [3] ciota [4] scale [5] dphi [6] dt
    // param [0] log(Mc) [1] log(Mt) [2] chi1 [3] chi2 [4] phi0  [5] tp0 [6] log(DL0) then relative amplitudes, time, phases

    
    // update the extrinsic parameters
    for(k = 1; k <= NC; k++)
    {
        
        // move reference point from geocenter to ref detector
        // Note that sky[4],sky[5], sky[6] hold shifts relative to the reference geocenter waveform
        // To map back to the reference detector
        
        ciota = skyh[k][3];
        Ap = (1.0+ciota*ciota)/2.0;
        Ac = -ciota;
        alpha = skyh[k][0];
        sindelta = skyh[k][1];
        psi = skyh[k][2];
        
        F_ant(psi, alpha, sindelta, net->GMST, &Fplus, &Fcross, net->labels[0]);
        Fs = sqrt(Ap*Ap*Fplus*Fplus+Ac*Ac*Fcross*Fcross);
        lambda = atan2(Ac*Fcross,Ap*Fplus);
        if(lambda < 0.0) lambda += TPI;
        
        Times(alpha, sindelta, net->GMST, dtimes);
        
        ciota = skyx[k][3];
        Ap = (1.0+ciota*ciota)/2.0;
        Ac = -ciota;
        alpha = skyx[k][0];
        sindelta = skyx[k][1];
        psi = skyx[k][2];
        
        F_ant(psi, alpha, sindelta, net->GMST, &Fplus, &Fcross, net->labels[0]);
        Fs2 = sqrt(Ap*Ap*Fplus*Fplus+Ac*Ac*Fcross*Fcross);
        lambda2 = atan2(Ac*Fcross,Ap*Fplus);
        if(lambda2 < 0.0) lambda2 += TPI;
        
        Times(alpha, sindelta, net->GMST, dtimes2);
        
        paramx[k][4] += 0.5*(skyx[k][5]+lambda2-lambda);
        while(paramx[k][4] > PI) paramx[k][4] -= PI;
        while(paramx[k][4] < 0.0) paramx[k][4] += PI;
        paramx[k][5] += skyx[k][6]+dtimes2[net->labels[0]]-dtimes[net->labels[0]];
        paramx[k][6] -= log(skyx[k][4]*Fs2/Fs);
        
        // sky will be re-aligned with geocenter so reset
        skyx[k][4] = 1.0;
        skyx[k][5] = 0.0;
        skyx[k][6] = 0.0;
       
    }
    
    
    printf("Swap Acceptance = %f\n", (double)sacc/(double)(scount));
    printf("MCMC Acceptance = %f\n", (double)ac/(double)(mcount));
    printf("Ring Acceptance = %f\n", (double)rca/(double)(rc));
    printf("Fisher Acceptance = %f\n", (double)fac/(double)(fc));
    printf("Jiggle Acceptance = %f\n", (double)uac/(double)(uc));
    
    
    
    free_double_matrix(skyy,NC+1);
    free_double_matrix(skyh,NC+1);
    free_double_matrix(sky,NC+1);
    free_double_vector(logLx);
    free_double_vector(sqH);
    free_double_tensor(fishskyx,NC+1,NS);
    free_double_tensor(fishskyy,NC+1,NS);
    free_double_tensor(skyvecsx,NC+1,NS);
    free_double_tensor(skyvecsy,NC+1,NS);
    free_double_matrix(skyevalsx,NC+1);
    free_double_matrix(skyevalsy,NC+1);
    free_double_vector(jump);
    
    free(dtimes);
    free(dtimes2);
    free(param);
    
    return;
    
}

void AmpPhase(double *Af, double *Pf, RealVector *freq, double *params, int N)
{
    
    AmpPhaseFDWaveform *ap = NULL;
    double phi0, fRef_in, mc, q, m1_SI, m2_SI, chi1, chi2, f_min, f_max, distance;
    int ret, flag, i;
    double p, cp, sp;
    double f, x, y, deltaF, ts, Tobs;
    
    double mt, eta, dm;
    
    Tobs = 1.0/freq->data[1];
    
    mc = exp(params[0]);
    mt = exp(params[1]);
    
    eta = pow((mc/mt), (5.0/3.0));
    if(eta > 0.25)
    {
        dm = 0.0;
    }
    else
    {
        dm = sqrt(1.0-4.0*eta);
    }
    m1_SI = mt*(1.0+dm)/2.0;
    m2_SI = mt*(1.0-dm)/2.0;
    
 
    chi1 = params[2];
    chi2 = params[3];
    
    phi0 = params[4];
    ts = Tobs-params[5];
    distance = exp(params[6]);
    fRef_in = fref;
    
    ret = IMRPhenomDGenerateh22FDAmpPhase(&ap, freq, phi0, fRef_in, m1_SI, m2_SI, chi1, chi2, distance);
    
    for (i=1; i< N/2; i++)
    {
       f = freq->data[i];
       Pf[i] = TPI*f*ts-ap->phase[i];
       Af[i] = h22fac*ap->amp[i];
    }
    
    
    DestroyAmpPhaseFDWaveform(ap);
    
}



void PDwave(double *wavef, RealVector *freq, double tcmin, double *params, int N)
{
    
    AmpPhaseFDWaveform *ap = NULL;
    double phi0, fRef_in, mc, q,  m1_SI, m2_SI, chi1, chi2, f_min, f_max, distance;
    int ret, flag, i;
    double p, cp, sp, Amp, fs;
    double f, x, y, deltaF, ts, Tobs, sqT;
    
    double mt, eta, dm;
    
    Tobs = 1.0/freq->data[1];
    sqT = sqrt(Tobs);
    
    // find safe start frequency for the current masses, spins
    // We do not want the waveform to have (much) support outside the time segment
    // Use the reference trigger time to find the frequency at the start of the segment
    x = params[5];
    params[5] = tcmin; // temporarily set merger time lowest allowed value
    fs = fbegin(params);
    params[5] = x;  // restore to assigned value
    
    mc = exp(params[0]);
    mt = exp(params[1]);
    
    eta = pow((mc/mt), (5.0/3.0));
    if(eta > 0.25)
    {
        dm = 0.0;
    }
    else
    {
    dm = sqrt(1.0-4.0*eta);
    }
    m1_SI = mt*(1.0+dm)/2.0;
    m2_SI = mt*(1.0-dm)/2.0;
    
    /*q = params[1];
    m2_SI = mc*pow(1.0+q, 0.2)/pow(q,0.6);
    m1_SI = q*m2_SI;*/
    
    chi1 = params[2];
    chi2 = params[3];
    
    phi0 = params[4];
    ts = Tobs-params[5];
    distance = exp(params[6]);
    
 
    fRef_in = fref;
    
    ret = IMRPhenomDGenerateh22FDAmpPhase(&ap, freq, phi0, fRef_in, m1_SI, m2_SI, chi1, chi2, distance);
    
    for (i=1; i< N/2; i++)
    {
        wavef[i] = 0.0;
        wavef[N-i] = 0.0;
        
        f = freq->data[i];
        if(f > fs)
        {
        p = TPI*f*ts;
        cp = cos(p);
        sp = sin(p);
        x = cos(ap->phase[i]);
        y = sin(ap->phase[i]);
        Amp = h22fac*ap->amp[i]/sqT;
        wavef[i] = Amp*(x*cp+y*sp);
        wavef[N-i] = Amp*(x*sp-y*cp);
        }
    }
    
    wavef[0] = 0.0;
    wavef[N/2] = 0.0;
    
    DestroyAmpPhaseFDWaveform(ap);
    
}

void templates(struct Net *net, double **hwave, RealVector *freq, double *params, int N)
{
    
    AmpPhaseFDWaveform *ap = NULL;
    double phi0, fRef_in, mc, q, m1_SI, m2_SI, chi1, chi2, f_min, f_max, distance;
    int ret, flag, i, j;
    double p, cp, sp;
    double f, x, y, deltaF, ts, Tobs, sqT;
    double pd, Ar, td, A, fs;
    
    double mt, eta, dm;
    
    Tobs = 1.0/freq->data[1];
    sqT = sqrt(Tobs);
    
    fs = fbegin(params);
    
    mc = exp(params[0]);
    mt = exp(params[1]);
    
    eta = pow((mc/mt), (5.0/3.0));
    if(eta > 0.25)
    {
        dm = 0.0;
    }
    else
    {
        dm = sqrt(1.0-4.0*eta);
    }
    m1_SI = mt*(1.0+dm)/2.0;
    m2_SI = mt*(1.0-dm)/2.0;
    
    
    chi1 = params[2];
    chi2 = params[3];
    
    phi0 = params[4];
    ts = Tobs-params[5];

    distance = exp(params[6]);
    
    fRef_in = fref;
    
    ret = IMRPhenomDGenerateh22FDAmpPhase(&ap, freq, phi0, fRef_in, m1_SI, m2_SI, chi1, chi2, distance);
        
    for (j=0; j< net->Nifo; j++)
    {
        pd = 0.0; // phase offset
        td = 0.0; // time offset
        Ar = 1.0; // amplitude ratio
            
        if(j > 0)
        {
         pd = params[(j-1)*3+NX]; // phase offset
         td = params[(j-1)*3+NX+1]; // time offset
         Ar = params[(j-1)*3+NX+2]; // amplitude
        }
        hwave[j][0] = 0.0;
        hwave[j][N/2] = 0.0;
            
         for (i=1; i< N/2; i++)
         {
            f = freq->data[i];
             
            if(f > fs)
            {
            A = Ar*h22fac*ap->amp[i]/sqT;
            p = ap->phase[i];
    
            x = TPI*f*(ts-td)+pd-p;
        
             hwave[j][i] = A*cos(x);
             hwave[j][N-i] = A*sin(x);
             }
             else
             {
                 hwave[j][i] = 0.0;
                 hwave[j][N-i] = 0.0;
             }
        }
    }
    
    DestroyAmpPhaseFDWaveform(ap);
    
}

void geotemplate(double *gwave, RealVector *freq, double *params, int N)
{
    
    AmpPhaseFDWaveform *ap = NULL;
    double phi0, fRef_in, mc, q, m1_SI, m2_SI, chi1, chi2, f_min, f_max, distance;
    int ret, flag, i, j, id;
    double p, cp, sp;
    double f, x, y, deltaF, ts, Tobs, sqT;
    double pd, Ar, td, A;
    double mt, eta, dm, fs;
    
    fs = fbegin(params);
 
    
    Tobs = 1.0/freq->data[1];
    sqT = sqrt(Tobs);
    
    mc = exp(params[0]);
    mt = exp(params[1]);
    
    eta = pow((mc/mt), (5.0/3.0));
    if(eta > 0.25)
    {
        dm = 0.0;
    }
    else
    {
        dm = sqrt(1.0-4.0*eta);
    }
    m1_SI = mt*(1.0+dm)/2.0;
    m2_SI = mt*(1.0-dm)/2.0;
    
    
    chi1 = params[2];
    chi2 = params[3];
    
    phi0 = 0.5*params[4];  // I'm holding the GW phase in [4], while PhenomD wants orbital
    ts = Tobs-params[5];
    
    distance = exp(params[6]);
    
    fRef_in = fref;
    
    ret = IMRPhenomDGenerateh22FDAmpPhase(&ap, freq, phi0, fRef_in, m1_SI, m2_SI, chi1, chi2, distance);

        gwave[0] = 0.0;
        gwave[N/2] = 0.0;
        
        for (i=1; i< N/2; i++)
        {
            f = freq->data[i];
            if(f > fs)
            {
            A = h22fac*ap->amp[i]/sqT;
            p = ap->phase[i];
            x = TPI*f*ts-p;
            gwave[i] = A*cos(x);
            gwave[N-i] = A*sin(x);
            }
            else
            {
                gwave[i] = 0.0;
                gwave[N-i] = 0.0;
            }
        }
    
    DestroyAmpPhaseFDWaveform(ap);
    
}

void fulltemplates(struct Net *net, double **hwave, RealVector *freq, double *params, int N)
{
    
    AmpPhaseFDWaveform *ap = NULL;
    double phi0, fRef_in, mc, q, m1_SI, m2_SI, chi1, chi2, f_min, f_max, distance;
    int ret, flag, i, j, id;
    double p, cp, sp;
    double f, x, y, deltaF, ts, Tobs, sqT;
    double pd, Ar, td, A;
    double alpha, sindelta, psi, ciota, lambda;
    double Ap, Ac;
    double Fplus, Fcross, Fs;
    double mt, eta, dm, fs;
    double *dtimes;
    
    dtimes = (double*)malloc(sizeof(double)* 5);
    
    fs = fbegin(params);
    
    Tobs = 1.0/freq->data[1];
    sqT = sqrt(Tobs);
    
    mc = exp(params[0]);
    mt = exp(params[1]);
    
    eta = pow((mc/mt), (5.0/3.0));
    if(eta > 0.25)
    {
        dm = 0.0;
    }
    else
    {
        dm = sqrt(1.0-4.0*eta);
    }
    m1_SI = mt*(1.0+dm)/2.0;
    m2_SI = mt*(1.0-dm)/2.0;
    
    
    chi1 = params[2];
    chi2 = params[3];
    
    phi0 = 0.5*params[4];  // [4] is the GW phase, while PhenomD uses orbital phase
    ts = Tobs-params[5];
    
    distance = exp(params[6]);
    
    fRef_in = fref;
    
    ret = IMRPhenomDGenerateh22FDAmpPhase(&ap, freq, phi0, fRef_in, m1_SI, m2_SI, chi1, chi2, distance);
    
    alpha = params[7];
    sindelta = params[8];
    psi = params[9];
    ciota = params[10];
    
    Ap = (1.0+ciota*ciota)/2.0;
    Ac = -ciota;
    
    Times(alpha, sindelta, net->GMST, dtimes);
    
    for (id=0; id< net->Nifo; id++)
    {
        
        F_ant(psi, alpha, sindelta, net->GMST, &Fplus, &Fcross, net->labels[id]);
        Fs = sqrt(Ap*Ap*Fplus*Fplus+Ac*Ac*Fcross*Fcross);  // magnitude of response
        lambda = atan2(Ac*Fcross,Ap*Fplus);
        if(lambda < 0.0) lambda += TPI;
        td = dtimes[net->labels[id]];
        
        hwave[id][0] = 0.0;
        hwave[id][N/2] = 0.0;
        
        for (i=1; i< N/2; i++)
        {
            f = freq->data[i];
            if(f > fs)
            {
            A = Fs*h22fac*ap->amp[i]/sqT;
            p = ap->phase[i];
            x = TPI*f*(ts-td)+lambda-p;
            hwave[id][i] = A*cos(x);
            hwave[id][N-i] = A*sin(x);
            }
            else
            {
                hwave[id][i] = 0.0;
                hwave[id][N-i] = 0.0;
            }
        }
    }
    
    free(dtimes);
    DestroyAmpPhaseFDWaveform(ap);
    
}


void fullphaseamp(struct Net *net, double **amp, double **phase, RealVector *freq, double *params, int N)
{
    
    AmpPhaseFDWaveform *ap = NULL;
    double phi0, fRef_in, mc, q, m1_SI, m2_SI, chi1, chi2, f_min, f_max, distance;
    int ret, flag, i, j, id;
    double p, cp, sp;
    double f, x, y, deltaF, ts, Tobs, sqT;
    double pd, Ar, td, A;
    double alpha, sindelta, psi, ciota, lambda;
    double Ap, Ac;
    double Fplus, Fcross, Fs;
    double mt, eta, dm, fs;
    double *dtimes;
    
    dtimes = (double*)malloc(sizeof(double)* 5);
    
    fs = fbegin(params);
    
    Tobs = net->Tobs;
    sqT = sqrt(Tobs);
    
    mc = exp(params[0]);
    mt = exp(params[1]);
    
    eta = pow((mc/mt), (5.0/3.0));
    if(eta > 0.25)
    {
        dm = 0.0;
    }
    else
    {
        dm = sqrt(1.0-4.0*eta);
    }
    m1_SI = mt*(1.0+dm)/2.0;
    m2_SI = mt*(1.0-dm)/2.0;
    
    
    chi1 = params[2];
    chi2 = params[3];
    
    phi0 = 0.5*params[4]; // I'm holding the GW phase in [4], while PhenomD wants orbital
    ts = Tobs-params[5];
    
    distance = exp(params[6]);
    
    fRef_in = fref;
    
    ret = IMRPhenomDGenerateh22FDAmpPhase(&ap, freq, phi0, fRef_in, m1_SI, m2_SI, chi1, chi2, distance);
    
    alpha = params[7];
    sindelta = params[8];
    psi = params[9];
    ciota = params[10];
    
    Ap = (1.0+ciota*ciota)/2.0;
    Ac = -ciota;
    
    Times(alpha, sindelta, net->GMST, dtimes);
    
    for (id=0; id< net->Nifo; id++)
    {
        
        F_ant(psi, alpha, sindelta, net->GMST, &Fplus, &Fcross, net->labels[id]);
        Fs = sqrt(Ap*Ap*Fplus*Fplus+Ac*Ac*Fcross*Fcross);  // magnitude of response
        lambda = atan2(Ac*Fcross,Ap*Fplus);
        if(lambda < 0.0) lambda += TPI;
        td = dtimes[net->labels[id]];
        
        amp[id][0] = 0.0;
        phase[id][0] = 0.0;
        
        for (i=0; i< N/2; i++)
        {
            f = freq->data[i];
            A = Fs*h22fac*ap->amp[i]/sqT;
            p = ap->phase[i];
            x = TPI*f*(ts-td)+lambda-p;
            y = 1.0;
            if(f < fs) y = exp(10.0*(f-fs)/fs);
            amp[id][i] = A*y;
            phase[id][i] = x;
        }
        
    }
    
    free(dtimes);
    DestroyAmpPhaseFDWaveform(ap);
    
}


void PDtemplates(double *waveH, double *waveL, RealVector *freq, double *params, int N)
{
    
    AmpPhaseFDWaveform *ap = NULL;
    double phi0, fRef_in, mc, q, m1_SI, m2_SI, chi1, chi2, f_min, f_max, distance;
    int ret, flag, i;
    double p, cp, sp;
    double f, x, y, deltaF, ts, Tobs, sqT;
    double pd, Ar, td, A, fs;
    
    double mt, eta, dm;
    
    Tobs = 1.0/freq->data[1];
    sqT = sqrt(Tobs);
    
    mc = exp(params[0]);
    mt = exp(params[1]);
    
    eta = pow((mc/mt), (5.0/3.0));
    if(eta > 0.25)
    {
        dm = 0.0;
    }
    else
    {
        dm = sqrt(1.0-4.0*eta);
    }
    m1_SI = mt*(1.0+dm)/2.0;
    m2_SI = mt*(1.0-dm)/2.0;
    
    
    chi1 = params[2];
    chi2 = params[3];
    
    phi0 = params[4];
    ts = Tobs-params[5];
    distance = exp(params[6]);
    
    pd = params[7]; // phase offset L from H
    td = params[8]; // time offset L from H
    Ar = params[9]; // amplitude ratio L/H
    
    fRef_in = fref;
    
    ret = IMRPhenomDGenerateh22FDAmpPhase(&ap, freq, phi0, fRef_in, m1_SI, m2_SI, chi1, chi2, distance);

    
    
    for (i=1; i< N/2; i++)
    {
        f = freq->data[i];
        A = h22fac*ap->amp[i]/sqT;
        
        
        p = TPI*f*ts-ap->phase[i];
        waveH[i] = A*cos(p);
        waveH[N-i] = A*sin(p);
        
        A *= Ar;
        p = TPI*f*(ts-td)+pd-ap->phase[i];
   
        waveL[i] = A*cos(p);
        waveL[N-i] = A*sin(p);
    }
    
    waveH[0] = 0.0;
    waveL[0] = 0.0;
    waveH[N/2] = 0.0;
    waveL[N/2] = 0.0;

    
    DestroyAmpPhaseFDWaveform(ap);
    
}

void log_likelihood_scan(struct Net *net, double **D, double *params, RealVector *freq, double **SN, double *Larray, double **Tarray, int N, double Tobs)
{
    int i, j, k, NMAX;
    int ii, jj, id;
    double sum;
    double *HH;
    double logL;
    double *lmax;
    double HD, LD, dt, x;
    double *norm, *delt, *pshift;
    double **HS, *HLS;
    int **shift;
    int **max, *tk;
    double **HC, **HF;
    int imin, imax;
    
    double *h;
    
    //Chi^2 = (d-h|d-h)/Sn
    
    h = double_vector(N);
    HS = double_matrix(net->Nifo,N);
    HC = double_matrix(net->Nifo,N);
    HF = double_matrix(net->Nifo,N);
    lmax = double_vector(N);
    tk = int_vector(net->Nifo);
    max = int_matrix(net->Nifo,N); // holds reference offset for each detector
    delt = double_vector(net->Nifo);
    norm = double_vector(net->Nifo);
    HH = double_vector(net->Nifo);
    pshift = double_vector(net->Nifo);
    
    imin = (int)(Tobs*fmin);
    imax = (int)(Tobs*fmax);
    if(imax > N/2) imax = N/2;
    
    dt = Tobs/(double)(N);
    
    // We always shift reference time to center
    params[5] = Tobs/2.0;
    
    // reference intrinsic waveform
    PDwave(h, freq, net->tmin, params, N);
    
    for(k = 0; k < net->Nifo; k++)
    {
        pbt_shift(HC[k], HF[k], D[k], h, SN[k], imin, imax, N);
        for(i = -N/2; i < N/2; i++) HS[k][i+N/2] = 0.0;
        for(i = 0; i < N/2; i++) HS[k][i+N/2] += sqrt(HC[k][i+N/2]*HC[k][i+N/2]+HF[k][i+N/2]*HF[k][i+N/2]);
        for(i = -N/2; i < 0; i++) HS[k][i+N/2] += sqrt(HC[k][N+i+N/2]*HC[k][N+i]+HF[k][N+i+N/2]*HF[k][N+i+N/2]);
    }
    
    // Allow time-travel time difference between detectors
    for(id = 1; id < net->Nifo; id++)  tk[id] = (int)((net->tds[id])/dt);

    
    // Note: there may not be a sky location that allows for the time delays between all detector pairs
    // the skyfind subroutine does its best to find something close to an unphysical maximum
    
    for (i = -N/2; i < N/2; ++i)
    {
        for(k = 1; k < net->Nifo; k++)
        {
            // find largest value within allowed time diference between each detector and reference
            x = 0.0;
            for (j = -tk[k]; j <= tk[k]; ++j)
            {
                ii = i+j;
                if(ii >= -N/2 && ii < N/2)
                {
                    if(HS[k][ii+N/2] > x)
                    {
                        x = HS[k][ii+N/2];
                        max[k][i+N/2] = ii;
                    }
                }
            }
        }
        
    }
    
    
    for(k = 1; k < net->Nifo; k++) HH[k] = fourier_nwip(h, h, SN[k], imin, imax, N);

    
    for (j = -N/2; j < N/2; ++j)
    {
    
    // start with reference detector
    HD = 2.0*(double)(N)*HS[0][j+N/2];
    delt[0] = dt*(double)(j);
    norm[0] = HD/HH[0];
        if(j >= 0)
        {
            pshift[0] = atan2(HF[0][j+N/2],HC[0][j+N/2]);
        }
        else
        {
            pshift[0] = atan2(HF[0][N+j+N/2],HC[0][N+j+N/2]);
        }

    
    logL = (HD*HD/HH[0])/2.0;
    
    for(k = 1; k < net->Nifo; k++)
    {
        HD = 2.0*(double)(N)*HS[k][max[k][j+N/2]];
        
        
        delt[k] = dt*(double)(max[k][j+N/2]);
        if(max[k][j+N/2] >= 0)
        {
            pshift[k]  = atan2(HF[k][max[k][j+N/2]],HC[k][max[k][j+N/2]]);
        }
        else
        {
            pshift[k]  = atan2(HF[k][N+max[k][j+N/2]],HC[k][N+max[k][j+N/2]]);
        }
        norm[k] = HD/HH[k];
        
        logL += (HD*HD/HH[k])/2.0;
        Tarray[k][j+N/2] = dt*(double)(j-max[k][j+N/2]);
        
    }
        
        Larray[j+N/2] = logL;
        
    }
    
    free_double_matrix(HS,net->Nifo);
    free_double_matrix(HC,net->Nifo);
    free_double_matrix(HF,net->Nifo);
    free_double_vector(lmax);
    free_int_vector(tk);
    free_int_matrix(max,net->Nifo);
    free_double_vector(delt);
    free_double_vector(norm);
    free_double_vector(HH);
    free_double_vector(pshift);
    free_double_vector(h);
    
    return;
}



double log_likelihood_intrinsic(struct Net *net, double **D, double *params, RealVector *freq, double **SN, int N, double Tobs)
{
    int i, j;
    double logL;
    double **twave;
    double *snrsq;
    double HH, HD, x;
    int imin, imax;
    
    logL = 0.0;
    

    if(lhold == 0)
    {
        
        imin = (int)(Tobs*fmin);
        imax = (int)(Tobs*fmax);
        if(imax > N/2) imax = N/2;
    
    twave = double_matrix(net->Nifo,N);
    snrsq = double_vector(net->Nifo);
        
    templates(net, twave, freq, params, N);
    
    logL = 0.0;
    for(i = 0; i < net->Nifo; i++)
    {
        HH = fourier_nwip(twave[i], twave[i], SN[i], imin, imax, N);
        HD = fourier_nwip(D[i], twave[i], SN[i], imin, imax, N);
        x = HD-0.5*HH;
        snrsq[i] = 2.0*x;
        logL += x;
    }
        
    if(net->Nifo > 1)
    {
        // count detectors where snr > 5
        j = 0;
        for(i = 0; i < net->Nifo; i++)
        {
            if(snrsq[i] > 25.0) j++;
        }
        // require a signal in at least two
        if(j < 2) logL = 0.0;
    }
    
    free_double_matrix(twave,net->Nifo);
    free_double_vector(snrsq);
    }
    
    return(logL);
    
}

double log_likelihood_print(struct Net *net, double **D, double *params, RealVector *freq, double **SN, int N, double Tobs)
{
    int i;
    double logL;
    double **twave;
    double HH, HD, x;
    int imin, imax;
    
    logL = 0.0;
    

    if(lhold == 0)
    {
        
        imin = (int)(Tobs*fmin);
        imax = (int)(Tobs*fmax);
        if(imax > N/2) imax = N/2;
    
    twave = double_matrix(net->Nifo,N);
    
    templates(net, twave, freq, params, N);
    
    logL = 0.0;
    for(i = 0; i < net->Nifo; i++)
    {
        HH = fourier_nwip(twave[i], twave[i], SN[i], imin, imax, N);
        HD = fourier_nwip(D[i], twave[i], SN[i], imin, imax, N);
        x = HD-0.5*HH;
        printf("%d %e\n", i, x);
        logL += x;
    }
    
    free_double_matrix(twave,net->Nifo);
    }
    
    return(logL);
    
}


void SNRvsf(struct Net *net, double **D, double *params, RealVector *freq, double **SN, int N, double Tobs)
{
    int i, j;
    double logL;
    double **twave;
    double HH, HD, x, y;
    int imin, imax;
    char filename[1024];
    FILE *out;
    
    
        
        twave = double_matrix(net->Nifo,N);
        
        fulltemplates(net, twave, freq, params, N);
        
    
        for(i = 0; i < net->Nifo; i++)
        {
            
            sprintf(filename, "SNRf_%d.dat", net->labels[i]);
            out = fopen(filename,"w");
            
            HH = 0.0;
            HD = 0.0;
            
            for(j = 1; j < N/2; j++)
            {
                
            x = 4.0*(twave[i][j]*twave[i][j]+twave[i][N-j]*twave[i][N-j])/SN[i][j];
            y = 4.0*(twave[i][j]*D[i][j]+twave[i][N-j]*D[i][N-j])/SN[i][j];
                
            HH += x;
            HD += y;
                
                fprintf(out, "%e %e %e %e %e %e\n", (double)(j)/Tobs, x, y, HH, HD, HD/sqrt(HH+1.0e-6));
                
            }
            
            fclose(out);
            
        }
        
        free_double_matrix(twave,net->Nifo);

    
}


double log_likelihood_full(struct Net *net, double **D, double *params, RealVector *freq, double **SN, double *rho, int N, double Tobs)
{
    int i;
    double logL;
    double **twave;
    double HH, HD, x;
    int imin, imax;
    
    logL = 0.0;
    
    
    if(lhold == 0)
    {
        
        imin = (int)(Tobs*fmin);
        imax = (int)(Tobs*fmax);
        if(imax > N/2) imax = N/2;
        
        twave = double_matrix(net->Nifo,N);
    
        fulltemplates(net, twave, freq, params, N);
        
        logL = 0.0;
        for(i = 0; i < net->Nifo; i++)
        {
            HH = fourier_nwip(twave[i], twave[i], SN[i], imin, imax, N);
            HD = fourier_nwip(D[i], twave[i], SN[i], imin, imax, N);
            rho[i] = HD/sqrt(HH);
            logL += HD-0.5*HH;
            //printf("%d %f %f\n", i, HH, HD);
        }
        
        free_double_matrix(twave,net->Nifo);
    }
    
    return(logL);
    
}

void limits(struct Net *net, double *params, RealVector *freq, double **SN, int *stst, double *ratio, int N, double Tobs)
{
    int i, id, j, M, jmin;
    int flag1, flag2;
    double **Ar, **Pr;
    double *SNRsq, *dSNRsq;
    double SNRsqT;
    double Re, Im;
    double fstop, dx, x, fs;
    
    FILE *out;
    
    fstop = 2.0*fringdown(params);  // twice ringdown frequency
    if(fstop < 64.0) fstop = 64.0; // just to be safe
    M = (int)(fstop*Tobs); // Covers the signal range
    if(M > N/2) M = N/2;  // can't exceed Nyqusit

    
    fs = fbegin(params);
    
    printf("%f %f\n", fs, fstop);
    
    if(fmin > fs)
    {
     jmin = (int)(Tobs*fmin);
     }
    else
     {
     jmin = (int)(Tobs*fs);
    }
    
    Ar = double_matrix(net->Nifo,M);
    Pr = double_matrix(net->Nifo,M);
    SNRsq = double_vector(M);
    dSNRsq = double_vector(M);
    
    fullphaseamp(net, Ar, Pr, freq, params, 2*M);
    
     for(j = 0; j < M; j++)
     {
         dSNRsq[j] = 0.0;
         SNRsq[j];
     }
    

    dx = 0.0;
    for(j = jmin; j < M; j++)
      {
        for(id = 0; id < net->Nifo; id++)
           {
               x = 4.0*(Ar[id][j]*Ar[id][j])/SN[id][j];
               dSNRsq[j] += x;
               dx += x;
            }
           SNRsq[j] = dx;
       }
    
    flag1 = 0;
    flag2 = 0;
    
    SNRsqT = SNRsq[M-1];
    
    printf("SNR = %f\n", sqrt(SNRsqT));
    
    
    stst[0] = jmin;
    stst[1] = M-1;
    
    for(j = jmin; j < M; j++)
    {
        if(flag1 == 0 && SNRsq[j] > 0.001)
        {
            flag1 = 1;
            stst[0] = j;
        }
        
        if(flag2 == 0 && SNRsqT-SNRsq[j] < 0.001)
        {
            flag2 = 1;
            stst[1] = j;
        }
        
        //printf("%d %d %f\n", j, stst[1], SNRsqT-SNRsq[j]);
    
    }
    
    if(flag2 == 0) stst[1] = M-1;
    if(stst[0] < jmin) stst[0] = jmin;
    
    printf("start frequency %f\n", (double)(stst[0])/Tobs);
    printf("end frequency %f\n", (double)(stst[1])/Tobs);
    
    out = fopen("ratio.dat","w");
    x = (double)(stst[1]-stst[0]);
    x = SNRsqT/x; // average contributon per significant bin
    
    for(j = 1; j < M; j++)
    {
        dx = dSNRsq[j]/x;
        ratio[j] = 1.0;
        if(dx > 0.0) ratio[j] = 1.0/(dx);
        fprintf(out,"%e %e %e\n", (double)(j)/Tobs, dx, SNRsq[j]);
    }
    fclose(out);

    free_double_matrix(Ar,net->Nifo);
    free_double_matrix(Pr,net->Nifo);
    free_double_vector(SNRsq);

    
}

void efix(struct Net *net, struct Het *het, double *params, double *min, double *max, RealVector *freq, double **SM, int N, double Tobs, double *eval, double **evec, double zs, int hflag)
{
    int i, j, k, flag;
    double alpha0, x, z, z0, alpha, alphanew;
    double dz, alphax;
    double leta, eta;
    double zmx, zmn;
    double dzmin, alpham, zm;
    double *px;
    
   
    zmx = zs*2.0;
    zmn = zs/2.0;
    
    px = double_vector(NP);

   for (i = 0; i < NP; ++i)
     {
         
        dzmin = 1.0e20;
        alpham = 1.0;
         
      alpha0 = 1.0;
         
    // [0] log(Mc) [1] log(Mt) [2] chi1 [3] chi2 [4] phi0 [5] tp [6] log(DL) [7] alpha [8] sindelta [9] psi [10] ciota
        
      k = 0;
      do
      {
      for (j = 0; j < NP; ++j) px[j] = params[j] + alpha0*eval[i]*evec[i][j];
      
      leta = (5.0/3.0)*(px[0]-px[1]);
      eta = exp(leta);
       if(eta > 0.25)
        {
         for (j = 0; j < NP; ++j) px[j] = params[j] - alpha0*eval[i]*evec[i][j];
       }
          
          // re-map angular parameters to their proper range
          
          x = px[4]/TPI;
          px[4] = (x-floor(x))*TPI;
          if(x < 0.0) px[4] += 1.0;
          x = px[7]/TPI;
          px[7] = (x-floor(x))*TPI;
          if(x < 0.0) px[7] += 1.0;
          x = px[9]/PI;
          px[9] = (x-floor(x))*PI;
          if(x < 0.0) px[9] += 1.0;
          
          while(px[4] > TPI)  px[4] -= TPI;
          while(px[4] < 0.0)  px[4] += TPI;
          while(px[7] > TPI)  px[7] -= TPI;
          while(px[7] < 0.0)  px[7] += TPI;
          while(px[9] > PI)   px[9] -= PI;
          while(px[9] < 0.0)  px[9] += PI;
          
          for (j = 0; j < NP; ++j) if(px[j] > max[j]) px[j] = max[j];
          for (j = 0; j < NP; ++j) if(px[j] < min[j]) px[j] = min[j];
      
          leta = (5.0/3.0)*(px[0]-px[1]);
          eta = exp(leta);
          if(eta > 0.25) px[0] = px[1] + 3.0/5.0*log(0.2499);
          if(eta < etamin) px[0] = px[1] + 3.0/5.0*log(etamin);
        
      if(hflag == 0)
      {
       z0 = diff(net, params, px, freq, SM,  N, Tobs);
      }
      else
      {
       z0 = het_diff(net, het, params, px);
      }
          
     //printf("B %f %f\n", alpha0, z0);
          
          dz = fabs(z0-zs);
          if(dz < dzmin)
          {
              dzmin = dz;
              alpham = alpha0;
              zm = z0;
          }
          
          if(z0 < zmn) alpha0 *= 2.0;
          if(z0 > zmx) alpha0 /= 1.9;
          
          k++;
          
      } while ( (z0 > zmx || z0 < zmn) && k < 15 && fabs(alpha0 < 1.0e4));
      
      
      alpha = alpha0*1.1;
          
        k = 0;
        do
        {
        for (j = 0; j < NP; ++j) px[j] = params[j] + alpha*eval[i]*evec[i][j];
            
        leta = (5.0/3.0)*(px[0]-px[1]);
        eta = exp(leta);
        if(eta > 0.25)
        {
          for (j = 0; j < NP; ++j) px[j] = params[j] - alpha*eval[i]*evec[i][j];
        }
            
        // re-map angular parameters to their proper range
            x = px[4]/TPI;
            px[4] = (x-floor(x))*TPI;
            if(x < 0.0) px[4] += 1.0;
            x = px[7]/TPI;
            px[7] = (x-floor(x))*TPI;
            if(x < 0.0) px[7] += 1.0;
            x = px[9]/PI;
            px[9] = (x-floor(x))*PI;
            if(x < 0.0) px[9] += 1.0;
            
        while(px[4] > TPI)  px[4] -= TPI;
        while(px[4] < 0.0)  px[4] += TPI;
        while(px[7] > TPI)  px[7] -= TPI;
        while(px[7] < 0.0)  px[7] += TPI;
        while(px[9] > PI)   px[9] -= PI;
        while(px[9] < 0.0)  px[9] += PI;
            
        for (j = 0; j < NP; ++j) if(px[j] > max[j]) px[j] = max[j];
        for (j = 0; j < NP; ++j) if(px[j] < min[j]) px[j] = min[j];
            
        
            leta = (5.0/3.0)*(px[0]-px[1]);
            eta = exp(leta);
            if(eta > 0.25) px[0] = px[1] + 3.0/5.0*log(0.2499);
            if(eta < etamin) px[0] = px[1] + 3.0/5.0*log(etamin);
            
            if(hflag == 0)
            {
             z = diff(net, params, px, freq, SM,  N, Tobs);
            }
            else
            {
             z = het_diff(net, het, params, px);
            }
            
           //printf("R %f %f\n", alpha, z);
            
            if(alpha > alpha0)
            {
            alphanew = alpha0 +(zs-z0)/(z-z0)*(alpha-alpha0);
            }
            else
            {
             alphanew = alpha +(zs-z)/(z0-z)*(alpha0-alpha);
            }
            
            dz = fabs(z-zs);
            if(dz < dzmin)
            {
              dzmin = dz;
              alpham = alpha;
              zm = z;
            }
            
            z0 = z;
            alpha0 = alpha;
            alpha = alphanew;
 
             k++;
             
        } while (fabs(z-zs) > 0.2 && k < 10 && fabs(alpha < 1.0e4));
         
          //printf("F %f %f\n", alpham, zm);
         
         // printf("\n");
         
         // kill any failed direction
         if(dzmin/zs > 1.0) alpham = 0.0;
         
         
         eval[i] *= alpham;
      
     }
}


void het_space(struct Net *net, struct Het *het, double *params, double *min, double *max, RealVector *freq, double **SN, int N, double Tobs)
{
    int i, id, j, k, ii, jj, kk, MN, MM, jmin, flag, djmax, M, J, NR;
    int jmid, boost;
    int pflag;
    double **Ar, ***Ap, **Pr, ***Pp;
    double ***Cs, ***Sn, ***DH;
    int *dfa;
    int *fgrid;
    double fstop;
    char filename[1024];
    double k0, k1, k2;
    double norm0, norm1;
    double x, y, z;
    double tol, tl, tola, tla;
    double cmid, smid, dmid;
    FILE *out;
    double *ratio;
    int *stst;
    double **fisher, **cov;
    double **evec, *eval;
    double **dparams, *px;
    double alpha0, z0, zs, alphanew;
    double leta, eta;
    double alphax, dz, LDLmin, Lmax;
    
    pflag = 1;  // print diagnostics if pflag = 1
   
    ratio = double_vector(N/2);
    stst = int_vector(2);
    
    limits(net, params, freq, SN, stst, ratio, N, Tobs);
    
    MN = stst[0];
    MM = stst[1];
    
    het->MN = MN;
    het->MM = MM;
    
    J = 2;
    tol = 0.003;
    
    het->J = J;
    
    dparams = double_matrix(NP,NP);
    fisher = double_matrix(NP,NP);
    evec = double_matrix(NP,NP);
    eval = double_vector(NP);
    px = double_vector(NP);
    
    dfa = int_vector(N/2);
    
    // [0] log(Mc) [1] log(Mt) [2] chi1 [3] chi2 [4] phi0 [5] tp [6] log(DL) [7] alpha [8] sindelta [9] psi [10] ciota

     Fisher_All(net, fisher, params, freq, SN, N, Tobs);
     FisherEvec(fisher, eval, evec, NP);
     efix(net, het, params, min, max, freq, SN, N, Tobs, eval, evec, 30.0, 0);
    

    for (i = 0; i < NP; ++i)
    {
        for (j = 0; j < NP; ++j) px[j] = params[j] + eval[i]*evec[i][j];
        
        leta = (5.0/3.0)*(px[0]-px[1]);
        eta = exp(leta);
        if(eta > 0.25)
        {
            for (j = 0; j < NP; ++j) px[j] = params[j] - eval[i]*evec[i][j];
        }
        
        // re-map angular parameters to their proper range
        x = px[4]/TPI;
        px[4] = (x-floor(x))*TPI;
        if(x < 0.0) px[4] += 1.0;
        x = px[7]/TPI;
        px[7] = (x-floor(x))*TPI;
        if(x < 0.0) px[7] += 1.0;
        x = px[9]/PI;
        px[9] = (x-floor(x))*PI;
        if(x < 0.0) px[9] += 1.0;
        while(px[4] > TPI)  px[4] -= TPI;
        while(px[4] < 0.0)  px[4] += TPI;
        while(px[7] > TPI)  px[7] -= TPI;
        while(px[7] < 0.0)  px[7] += TPI;
        while(px[9] > PI)   px[9] -= PI;
        while(px[9] < 0.0)  px[9] += PI;
        
        
        for (j = 0; j < NP; ++j) if(px[j] > max[j]) px[j] = max[j];
        for (j = 0; j < NP; ++j) if(px[j] < min[j]) px[j] = min[j];
            
        
            leta = (5.0/3.0)*(px[0]-px[1]);
            eta = exp(leta);
            if(eta > 0.25) px[0] = px[1] + 3.0/5.0*log(0.2499);
            if(eta < etamin) px[0] = px[1] + 3.0/5.0*log(etamin);
    

        
        for (j = 0; j < NP; ++j) dparams[i][j] = px[j];
                   
        }
    
    
   if(pflag == 1) out = fopen("fspace.dat","w");
    
    fgrid = int_vector(MM);
    
    Ar = double_matrix(net->Nifo,MM);
    Ap = double_tensor(NP,net->Nifo,MM);
    Pr = double_matrix(net->Nifo,MM);
    Pp = double_tensor(NP,net->Nifo,MM);
    Cs = double_tensor(NP,net->Nifo,MM);
    Sn = double_tensor(NP,net->Nifo,MM);
    DH = double_tensor(NP,net->Nifo,MM);
    
    fullphaseamp(net, Ar, Pr, freq, params, 2*MM);
    
    for(i = 0; i < NP; i++) fullphaseamp(net, Ap[i], Pp[i], freq, dparams[i], 2*MM);
    
    
    // phase difference vs frequency and ampltide ratio vs frequency
    for(i = 0; i < NP; i++)
    {
        for(id = 0; id < net->Nifo; id++)
        {
         for(j = MN; j < MM; j++)
          {
              Cs[i][id][j] = 1.0-Ap[i][id][j]/Ar[id][j]*cos(Pr[id][j]-Pp[i][id][j]);
              Sn[i][id][j] = Ap[i][id][j]/Ar[id][j]*sin(Pr[id][j]-Pp[i][id][j]);
              DH[i][id][j] = Cs[i][id][j]*Cs[i][id][j]+Sn[i][id][j]*Sn[i][id][j];
          }
        }
    }
    
   
    
    jmin = MN;
    djmax = (int)(100.0*Tobs);  // maximum spacing
    
    M = 0;
    fgrid[M] = jmin;
    
     if(pflag == 1)
     {
    fprintf(out, "%e ", (double)(jmin)/Tobs);
    for(i = 0; i < NP; i++) fprintf(out, "%e %e %e ", Cs[i][0][jmin], Sn[i][0][jmin], DH[i][0][jmin]);
    fprintf(out, "\n");
     }
    
        for(j = MN; j < MM; j++)
        {
            
            flag = 0;
            
            tl = tol*ratio[j];  // scale the tolerance by SNR contribution
            if(tl < tol) tl = tol;
            
            //tl = tol;
            
          // printf("%f %e\n", (double)(j)/Tobs, tl);
           
            for(id = 0; id < net->Nifo; id++)
              {
                  for(i = 0; i < NP; i++)
                  {
                    
                      
                      // linear interpolation at the mid-point
                      jmid = (j+jmin)/2;
                      
                      cmid = (Cs[i][id][jmin]*(double)(j-jmid)+ Cs[i][id][j]*(double)(jmid-jmin))/(double)(j-jmin);
                      smid = (Sn[i][id][jmin]*(double)(j-jmid)+ Sn[i][id][j]*(double)(jmid-jmin))/(double)(j-jmin);
                      dmid = (DH[i][id][jmin]*(double)(j-jmid)+ DH[i][id][j]*(double)(jmid-jmin))/(double)(j-jmin);
                      
                      if(fabs(Cs[i][id][jmid]-cmid) > tl || fabs(Sn[i][id][jmid]-smid) > tl || fabs(DH[i][id][jmid]-dmid) > tl)
                      {
                          flag = 1;
                          //printf("%d %d %d %f\n", j, id, i, tl);
                      }
                      
                  }
              }
            
            // max frequency spacing
            if((j-jmin) > djmax) flag = 1;
     
            
            if(flag == 1) // re-start
            {
                jmin = j;
                M++;
                fgrid[M] = jmin;
                 if(pflag == 1)
                 {
                fprintf(out, "%e ", (double)(j)/Tobs);
                for(i = 0; i < NP; i++) fprintf(out, "%e %e %e ", Cs[i][0][j], Sn[i][0][j], DH[i][0][j]);
                fprintf(out, "\n");
                 }
            }
            
        }
    
    if(pflag == 1)
       {
       out = fopen("df.dat","w");
       for (i = 0; i < M; ++i)
       {
           fprintf(out,"%e %e\n", (double)(fgrid[i])/Tobs, (double)(fgrid[i+1]-fgrid[i])/Tobs);
       }
       fclose(out);
       }
    
    for (i = 0; i < N/2; ++i) dfa[i] = 100;
    
    for(ii = 0;  ii < M; ii++)
    {
        j = fgrid[ii+1]-fgrid[ii];
        //printf("%d %d %d\n", ii, fgrid[ii], j);
        for (i = fgrid[ii]; i < fgrid[ii+1]; ++i)
        {
           dfa[i] = j;
          // printf("%d %d\n", i, dfa[i]);
        }
     }
    
      j = MM-fgrid[M];
       for (i = fgrid[M]; i < MM; ++i)
       {
           dfa[i] = j;
       }
    
    int *fgflat;
    
    fgflat = int_vector(N/2);
    
    
     j = MN;
     ii = 0;
     fgflat[ii] = j;
    
     do
     {
     jmin = dfa[j];
     //printf("%d %d %d\n", j, j+jmin*J, jmin);
     do
     {
     flag = 0;
     k = j+jmin*J;
     for (i = j; i <= k; ++i)
     {
         if(i < N/2)
         {
         if(dfa[i] < jmin)
         {
             jmin = dfa[i];
             flag = 1;
         }
         }
         else
         {
             flag = 0;
         }
     }
     //printf("%d %d %d\n", j, j+jmin*J, jmin);
     }while(flag == 1);
    
     j = k;
      
     for (i = 0; i < J; ++i)
     {
         ii++;
         fgflat[ii] = fgflat[ii-1]+jmin;
         //printf("%d %d\n", ii, fgflat[ii]);
     }
         
         
     }while(j < MM);
    
    // fix the last section so it ends at MM
    jmin = (MM-fgflat[ii-J])/J;
    
    if(jmin >= 1)
    {

    ii -= J;
    
   // printf("\n%d %d %d\n", jmin, fgflat[ii-1], fgflat[ii-1]+J*jmin);
    
    for (i = 0; i < J; ++i)
    {
        ii++;
        fgflat[ii] = fgflat[ii-1]+jmin;
    }
        
    }
    
    
   // printf("%d %d %d\n", fgflat[ii], fgflat[ii-J], MM);
    
    if(pflag == 1)
    {
     out = fopen("df_flat.dat","w");
     for (i = 0; i < ii; ++i)
     {
         fprintf(out,"%e %e\n", (double)(fgflat[i])/Tobs, (double)(fgflat[i+1]-fgflat[i])/Tobs);
     }
     fclose(out);
    }
    
    free_int_vector(fgflat);

    
       M = ii+1;
       NR = ii/J;
       het->NR = NR;
       het->M = M;
       het->freq = CreateRealVector(M);
       het->fgrid = int_vector(M);
       for (i = 0; i < M; ++i) het->fgrid[i] = fgflat[i];
       for (i = 0; i < M; ++i) het->freq->data[i] = ((double)fgflat[i])/Tobs;
        

    free_int_vector(dfa);
    free_double_vector(ratio);
    free_int_vector(stst);
    free_double_matrix(dparams,NP);
    free_double_matrix(fisher,NP);
    free_double_matrix(evec,NP);
    free_double_vector(eval);
    free_double_vector(px);
    free_int_vector(fgrid);
    free_double_matrix(Ar,net->Nifo);
    free_double_tensor(Ap,NP,net->Nifo);
    free_double_matrix(Pr,net->Nifo);
    free_double_tensor(Pp,NP,net->Nifo);
    free_double_tensor(Cs,NP,net->Nifo);
    free_double_tensor(Sn,NP,net->Nifo);
    free_double_tensor(DH,NP,net->Nifo);
    
}


double diff(struct Net *net, double *params, double *dparams, RealVector *freq, double **SN, int N, double Tobs)
{
    int i, id, j, M, imin;
    double **Ar, **Ap, **Pr, **Pp;
    double Re, Im;
    double fstop, dx, fs;
   
    
    fs = fbegin(params);
           
    if(fmin > fs)
     {
      imin = (int)(Tobs*fmin);
     }
    else
    {
    imin = (int)(Tobs*fs);
    }
    
    fstop = 2.0*fringdown(params);  // twice ringdown frequency
    if(fstop < 64.0) fstop = 64.0; // just to be safe
    M = 2*(int)(fstop*Tobs); // Covers the signal range
    if(M > N) M = N;  // can't exceed Nyqusit
    
    Ar = double_matrix(net->Nifo,M/2);
    Ap = double_matrix(net->Nifo,M/2);
    Pr = double_matrix(net->Nifo,M/2);
    Pp = double_matrix(net->Nifo,M/2);
    
    fullphaseamp(net, Ar, Pr, freq, params, M);
    fullphaseamp(net, Ap, Pp, freq, dparams, M);
    
    dx = 0.0;
    
    for(id = 0; id < net->Nifo; id++)
        {
         for(j = imin; j < M/2; j++)
            {
              Re = (Ar[id][j]-Ap[id][j]*cos(Pr[id][j]-Pp[id][j]));
              Im = Ap[id][j]*sin(Pr[id][j]-Pp[id][j]);
              dx += 4.0*(Re*Re+Im*Im)/SN[id][j];
            }
        }

    free_double_matrix(Ar,net->Nifo);
    free_double_matrix(Ap,net->Nifo);
    free_double_matrix(Pr,net->Nifo);
    free_double_matrix(Pp,net->Nifo);
    
    return(dx);
    
}

double het_diff(struct Net *net, struct Het *het, double *params, double *dparams)
{
    int i, id, j, M, J, NR, imin;
    int ii, jj, kk;
    double **Ar, **Ap, **Pr, **Pp;
    double Re, Im;
    double fstop, dx, fs;
    double *cc, *uu, *vv;
    double HH, x;

    M = het->M;
    J = het->J;
    NR = het->NR;

    Ar = double_matrix(net->Nifo,M);
    Ap = double_matrix(net->Nifo,M);
    Pr = double_matrix(net->Nifo,M);
    Pp = double_matrix(net->Nifo,M);
    
    cc = double_vector(M);
    uu = double_vector(J+1);
    vv = double_vector(J+1);
    
    fullphaseamp(net, Ar, Pr, het->freq, params, 2*M);
    fullphaseamp(net, Ap, Pp, het->freq, dparams, 2*M);
    
    dx = 0.0;
    
    for (id = 0 ; id < net->Nifo ; id++)  // loop over detectors
     {
                       
        HH = 0.0;
        x = 0.0;
                       
        // these are the slow terms
        for(jj = 0; jj < M; jj++)
        {
            Re = (Ar[id][jj]-Ap[id][jj]*cos(Pr[id][jj]-Pp[id][jj]));
            Im = Ap[id][jj]*sin(Pr[id][jj]-Pp[id][jj]);
            cc[jj] = 4.0*(Re*Re+Im*Im)/(het->amp[id][jj]*het->amp[id][jj]);
        }
                       
                       
          for(ii = 0; ii < NR; ii++) // loop over regions
            {
                                                            
                for(jj = 0; jj <= J; jj++)
                {
                 vv[jj] = cc[ii*J+jj];
                 uu[jj] = 0.0;
                }
                                                            
              // get Legendre coefficicients for slow term
              for(jj = 0; jj <= J; jj++)
               {
               for(kk = 0; kk <= J; kk++)
                {
                uu[jj] += het->IP[ii][jj][kk]*vv[kk];
                }
               }
                                                            
                for(jj = 0; jj <= J; jj++) HH += uu[jj]*het->SL[ii][id][jj];
                                                            
                if(ii < NR-2) x += het->aa[id][ii+1]*vv[J];  // correction for overcount
                                             
            }
                       
                // correction for overcount of edge values
            HH -= x;
         
         dx += HH;
         
     }
    
    free_double_vector(cc);
    free_double_vector(uu);
    free_double_vector(vv);

    free_double_matrix(Ar,net->Nifo);
    free_double_matrix(Ap,net->Nifo);
    free_double_matrix(Pr,net->Nifo);
    free_double_matrix(Pp,net->Nifo);
    
    return(dx);
    
}



void legendre_maker(int J, int U, double **P)
{
    double *gs, *a, *b;
    double c0, c1, c2, d0, d1, d2;
    int i, k;
    
    gs = double_vector(J+1);  // downsampled function
    a = double_vector(J+1);  // coefficients
    b = double_vector(J+1);  // data

    for (k = 0; k <= U; ++k)
    {
       
     c0 = 0.0;
     c1 = (double)(U-2*k);
     c2 = (double)(U);
     d0 = c2;
     d1 = c1+c1;
     d2 = d0;

     P[0][k] = 1.0;
     P[1][k] = c1/(double)(U);
 
    for (i = 1; i < J; ++i)
     {
        d0 += 2.0;
        d2 -= 2.0;
        c0 += d0;
        c1 += d1;
        c2 += d2;
        P[i+1][k] = (c1*P[i][k]-c0*P[i-1][k])/c2;
    }
        
    }
        
    free_double_vector(gs);
    free_double_vector(a);
    free_double_vector(b);

}

void freehet(struct Net *net, struct Het *het)
{
    free_double_tensor(het->IP,het->NR,(het->J)+1);
    free_double_tensor(het->SL,het->NR,net->Nifo);
    free_double_matrix(het->amp,net->Nifo);
    free_double_matrix(het->phase,net->Nifo);
    free_double_matrix(het->rc,net->Nifo);
    free_double_matrix(het->rs,net->Nifo);
    free_double_matrix(het->aa,net->Nifo);
    free_double_tensor(het->lc,het->NR,net->Nifo);
    free_double_tensor(het->ls,het->NR,net->Nifo);
    free_double_vector(het->pref);
    DestroyRealVector(het->freq);
    free_int_vector(het->fgrid);
}

void heterodyne(struct Net *net, struct Het *het, double **D, double *params, RealVector *freq, double **SN, int N, double Tobs)
{
    int i, id, j, k, n;
    double f, fstop, x, y, z;
    double **hb, **rb, **ampb, **phaseb;
    double cp, sp;
    double **RC, **RS, **AA;
    double logL, HH, HR;
    int M, MM, MN, NR, ds, J, U;
    
    M = het->M ;
    MM = het->MM;
    MN = het->MN;
    J = het->J;
    NR = het->NR;
    
    het->IP = double_tensor(NR,J+1,J+1);
    het->SL = double_tensor(NR,net->Nifo,J+1);
    het->amp = double_matrix(net->Nifo,M);
    het->phase = double_matrix(net->Nifo,M);
    het->rc = double_matrix(net->Nifo,NR);
    het->rs = double_matrix(net->Nifo,NR);
    het->aa = double_matrix(net->Nifo,NR);
    het->lc = double_tensor(NR,net->Nifo,J+1);
    het->ls = double_tensor(NR,net->Nifo,J+1);
    het->pref = double_vector(NP);  // keep a copy of the source parameters used for the heterodyne (used for testing)
    
    for(j = 0; j < NP; j++) het->pref[j] = params[j];
    
    hb = double_matrix(net->Nifo,N);
    rb = double_matrix(net->Nifo,N);
    RS = double_matrix(net->Nifo,N);
    RC = double_matrix(net->Nifo,N);
    AA = double_matrix(net->Nifo,N);
    ampb = double_matrix(net->Nifo,N/2);
    phaseb = double_matrix(net->Nifo,N/2);
    
    // full frequency template, amplitude and phase
    fulltemplates(net, hb, freq, params, N);
    fullphaseamp(net, ampb, phaseb, freq, params, N);
    
    // store the downsampled reference amplitude and phase
    for(id = 0; id < net->Nifo; id++)
       {
           for(j = 0; j < M; j++)
              {
                n = het->fgrid[j];
                het->amp[id][j] = ampb[id][n];
                het->phase[id][j] = phaseb[id][n];
              }
       }
    
    // form up resdiual
    for(id = 0; id < net->Nifo; id++)
    {
     for(j = 0; j < N; j++) rb[id][j] = D[id][j] - hb[id][j];
    }
    
    logL = 0.0;
    for(id = 0; id < net->Nifo; id++)
    {
        HH = fourier_nwip(hb[id], hb[id], SN[id], MN, MM, N);
        HR = fourier_nwip(rb[id], hb[id], SN[id], MN, MM, N);
        logL += HR+0.5*HH;
        //printf("HetRef %d %f %f\n", id, HH, HR);
    }

    het->logLb = logL;
    
    printf("ref logL = %e\n", logL);
    
    for(id = 0; id < net->Nifo; id++)
    {
        RC[id][0] = 0.0;
        RS[id][0] = 0.0;
      for(j = 1; j < N/2; j++)
      {
        cp = cos(phaseb[id][j]);
        sp = sin(phaseb[id][j]);
        x = ampb[id][j]/SN[id][j];
        RC[id][j] = 2.0*(rb[id][j]*cp+rb[id][N-j]*sp)*x;
        RS[id][j] = 2.0*(-rb[id][j]*sp+rb[id][N-j]*cp)*x;
        AA[id][j] = ampb[id][j]*ampb[id][j]/SN[id][j];
      }
    }
    
    // store boundary values for correcting overcount
    for(id = 0; id < net->Nifo; id++)
    {
        for (i = 0; i < NR; ++i)
          {
            k = het->fgrid[i*J];
            het->rc[id][i] = RC[id][k];
            het->rs[id][i] = RS[id][k];
            het->aa[id][i] = AA[id][k];
        }
    }
    

    
    double **P, **PP;
    
    PP = double_matrix(J+1,J+1); // downsampled Legendres
    
    for (i = 0; i < NR; ++i)
    {
      for(id = 0; id < net->Nifo; id++)
        {
         for(j = 0; j <= J; j++)
         {
             het->SL[i][id][j] = 0.0;
             het->lc[i][id][j] = 0.0;
             het->ls[i][id][j] = 0.0;
         }
        }
    }
    
    for (i = 0; i < NR; ++i)
    {
        //printf("%d %d %d\n", (i+1)*J, M, het->fgrid[(i+1)*J]);
        
        if((i+1)*J < M)
        {
            
         U = het->fgrid[(i+1)*J]-het->fgrid[i*J];
         
         P = double_matrix(J+1,U+1);
        
         legendre_maker(J, U, P);
        
        // downsampled Legendre
        for (j = 0; j <= J; ++j)
          {
            for (k = 0; k <= J; ++k)
              {
                  PP[j][k] = P[j][het->fgrid[i*J+k]-het->fgrid[i*J]];
              }
          }
        
        // store inverse of downsampled Legendre
        Inverse(PP, het->IP[i], J+1);
            
            /*
            printf("%d\n", i);
            for (j = 0; j <= J; ++j)
            {
              for (k = 0; k <= J; ++k)
                {
                    printf("%f ", het->IP[i][j][k]);
                }
                printf("\n");
            }
            printf("\n");
             */
        
        for(id = 0; id < net->Nifo; id++)
        {
          for(n = 0; n <= U; n++)
            {
                k = het->fgrid[i*J]+n;
               for(j = 0; j <= J; j++)
                {
                    het->SL[i][id][j] += AA[id][k]*P[j][n];
                    het->lc[i][id][j] += RC[id][k]*P[j][n];
                    het->ls[i][id][j] += RS[id][k]*P[j][n];
                }
            }
         }
        
        
        free_double_matrix(P,J+1);
            
        }
    }
    
    //printf("done with setup\n");
    
    free_double_matrix(PP,J+1);
    
    free_double_matrix(RC,net->Nifo);
    free_double_matrix(RS,net->Nifo);
    free_double_matrix(AA,net->Nifo);
    free_double_matrix(hb,net->Nifo);
    free_double_matrix(rb,net->Nifo);
    free_double_matrix(ampb,net->Nifo);
    free_double_matrix(phaseb,net->Nifo);
        
}



double log_likelihood_het(struct Net *net, struct Het *het, double *params, double Tobs)
{
    double **amp, **phase;
    double **hc, **hs;
    double HH, HR, logL;
    double x, y;
    double L0,L1;
    double **cc, **ss;
    double *uu, *vv;
    int i, j, k, id, M, J, NR;
    
    logL = 0.0;
    
    if(lhold == 0)
    {
    
    M = het->M;
    J = het->J;
    NR = het->NR;

    amp = double_matrix(net->Nifo,M);
    phase = double_matrix(net->Nifo,M);
    cc = double_matrix(net->Nifo,M);
    ss = double_matrix(net->Nifo,M);
    hc = double_matrix(net->Nifo,M);
    hs = double_matrix(net->Nifo,M);
    uu = double_vector(J+1);
    vv = double_vector(J+1);
    
    fullphaseamp(net, amp, phase, het->freq, params, 2*M);
    
    
    for(id = 0; id < net->Nifo; id++)
    {
        for(j = 0; j < M; j++)
        {
            cc[id][j] = 2.0*(het->amp[id][j]-amp[id][j]*cos(het->phase[id][j]-phase[id][j]));
            ss[id][j] = 2.0*amp[id][j]*sin(het->phase[id][j]-phase[id][j]);
            hc[id][j] = cc[id][j]/het->amp[id][j];
            hs[id][j] = ss[id][j]/het->amp[id][j];
            cc[id][j] = hc[id][j]*hc[id][j]+hs[id][j]*hs[id][j];
        }
    }
        
        logL = het->logLb;
        
        for(id = 0; id < net->Nifo; id++)
        {
          HH = 0.0;
          HR = 0.0;
          x = 0.0;
          y = 0.0;
       
            for(i = 0; i < NR; i++)
            {
                
                for(j = 0; j <= J; j++)
                {
                  vv[j] = cc[id][i*J+j];
                  uu[j] = 0.0;
                }
                
                // get Legendre coefficicients for slow term
                for(j = 0; j <= J; j++)
                    {
                     for(k = 0; k <= J; k++)
                        {
                            uu[j] += het->IP[i][j][k]*vv[k];
                        }
                    }
                
                for(j = 0; j <= J; j++) HH += uu[j]*het->SL[i][id][j];
                
                if(i < NR-2) x += het->aa[id][i+1]*vv[J];  // correction for overcount
                
                
                for(j = 0; j <= J; j++)
                {
                  vv[j] = hc[id][i*J+j];
                  uu[j] = 0.0;
                }
                
                for(j = 0; j <= J; j++)
                    {
                     for(k = 0; k <= J; k++)
                        {
                            uu[j] += het->IP[i][j][k]*vv[k];
                        }
                    }
                
                for(j = 0; j <= J; j++) HR += uu[j]*het->lc[i][id][j];
                
                if(i < NR-2) y += het->rc[id][i+1]*vv[J];  // correction for overcount
                
                
                    for(j = 0; j <= J; j++)
                           {
                             vv[j] = hs[id][i*J+j];
                             uu[j] = 0.0;
                           }
                           
                           for(j = 0; j <= J; j++)
                               {
                                for(k = 0; k <= J; k++)
                                   {
                                       uu[j] += het->IP[i][j][k]*vv[k];
                                   }
                               }
                           
                           for(j = 0; j <= J; j++) HR += uu[j]*het->ls[i][id][j];
                           
                           if(i < NR-2) y += het->rs[id][i+1]*vv[J];  // correction for overcount
          
        
        }
        
        // correction for overcount of edge values
        HH -= x;
        HR -= y;
    
        
        logL -= (HR+0.5*HH);
        
       // printf("Het %d %f %f %f %f\n", id, HH, HR, x, y);
        
        
    }
    
    free_double_vector(uu);
    free_double_vector(vv);
    free_double_matrix(cc,net->Nifo);
    free_double_matrix(ss,net->Nifo);
    free_double_matrix(hc,net->Nifo);
    free_double_matrix(hs,net->Nifo);
    free_double_matrix(amp,net->Nifo);
    free_double_matrix(phase,net->Nifo);
        
    }
    
    return(logL);
    
}

void Fisher_Het(struct Net *net, struct Het *het, double **fish, double *params, double Tobs)
{
    double *paramsP, *paramsM;
    double epsilon;
    double Scale, HH, x;
    double **amp, **phase;
    double ***ampP, ***phaseP;
    double ***ampM, ***phaseM;
    double **cc, *uu, *vv;
    int i, j, k, ii, jj, kk, l, id;
    int M, J, NR;
    int imin, imax;
    
    //printf("Fisher ");
    
    for (i = 0; i < NP; i++)
    {
        for (j = 0; j < NP; j++) fish[j][i] = 0.0;
    }
    
    for (i = 0; i < NP; i++) fish[i][i] = 1.0e4;
    
    
    if(fabs(params[2]) < 0.9999 &&  fabs(params[3]) < 0.9999)
    {
        
        M = het->M;
        J = het->J;
        NR = het->NR;

        amp = double_matrix(net->Nifo,M);
        phase = double_matrix(net->Nifo,M);
        ampP = double_tensor(NP,net->Nifo,M);
        phaseP = double_tensor(NP,net->Nifo,M);
        ampM = double_tensor(NP,net->Nifo,M);
        phaseM = double_tensor(NP,net->Nifo,M);
        cc = double_matrix(net->Nifo,M);
        uu = double_vector(J+1);
        vv = double_vector(J+1);
        
        epsilon = 1.0e-5;
        
        paramsP = (double*)malloc(sizeof(double)*(NP));
        paramsM = (double*)malloc(sizeof(double)*(NP));
        
        // Reference phase and amplitude
        fullphaseamp(net, amp, phase, het->freq, params, 2*M);
        
        for (i = 0 ; i < NP ; i++)
        {
            
            for (k = 0 ; k < NP ; k++)
            {
                paramsP[k] = params[k];
                paramsM[k] = params[k];
            }
            
            paramsP[i] += epsilon;
            paramsM[i] -= epsilon;
            
            fullphaseamp(net, ampP[i], phaseP[i], het->freq, paramsP, 2*M);
            fullphaseamp(net, ampM[i], phaseM[i], het->freq, paramsM, 2*M);
            
            // store the central difference in the Plus arrays
            for (id = 0 ; id < net->Nifo ; id++)  // loop over detectors
            {
            for (k = 0 ; k < M ; k++)
            {
                ampP[i][id][k] -= ampM[i][id][k];
                phaseP[i][id][k] -= phaseM[i][id][k];
                ampP[i][id][k] /= (2.0*epsilon);
                phaseP[i][id][k] /= (2.0*epsilon);
            }
            }
            
        }
        
        
        //printf("\n");
        
        for (i = 0; i < NP; i++)
        {
            for (j = 0; j < NP; j++) fish[j][i] = 0.0;
        }
        
            for (i = 0 ; i < NP ; i++)
            {
                for (j = i ; j < NP ; j++)
                {
                    
                    for (id = 0 ; id < net->Nifo ; id++)  // loop over detectors
                    {
                    
                    HH = 0.0;
                    x = 0.0;
                    
                    // these are the slow terms in the Fisher matrix sum
                    for(jj = 0; jj < M; jj++)
                    {
                     cc[id][jj] = 4.0*(ampP[i][id][jj]*ampP[j][id][jj]
                                      +amp[id][jj]*amp[id][jj]*phaseP[i][id][jj]*phaseP[j][id][jj]);
                     cc[id][jj] /= (het->amp[id][jj]*het->amp[id][jj]);
                    }
                    
                    
                for(ii = 0; ii < NR; ii++)
                 {
                                                         
                    for(jj = 0; jj <= J; jj++)
                    {
                      vv[jj] = cc[id][ii*J+jj];
                      uu[jj] = 0.0;
                    }
                                                         
                    // get Legendre coefficicients for slow term
                   for(jj = 0; jj <= J; jj++)
                    {
                    for(kk = 0; kk <= J; kk++)
                     {
                      uu[jj] += het->IP[ii][jj][kk]*vv[kk];
                     }
                    }
                                                         
                    for(jj = 0; jj <= J; jj++) HH += uu[jj]*het->SL[ii][id][jj];
                                                         
                    if(ii < NR-2) x += het->aa[id][ii+1]*vv[J];  // correction for overcount
                                          
                 }
                    
                    // correction for overcount of edge values
                    HH -= x;
                    fish[i][j] += HH;
                }
            }
            
        }
        
        
        /* fill in lower triangle */
        
        for (i = 0; i < NP; i++)
        {
            for (j = i+1; j < NP; j++)
                fish[j][i] = fish[i][j];
        }
        
        /*
         for (i = 0; i < NP; i++)
         {
         for (j = 0; j < NP; j++)
         {
         printf("%.3e ", fish[i][j]);
         }
         printf("\n");
         }
        printf("\n"); */
        
        free(paramsP);
        free(paramsM);
        free_double_matrix(amp,net->Nifo);
        free_double_matrix(phase,net->Nifo);
        free_double_tensor(ampP,NP,net->Nifo);
        free_double_tensor(phaseP,NP,net->Nifo);
        free_double_tensor(ampM,NP,net->Nifo);
        free_double_tensor(phaseM,NP,net->Nifo);
        free_double_matrix(cc,net->Nifo);
        free_double_vector(uu);
        free_double_vector(vv);
        
    }
    
}



double log_likelihood_test(struct Net *net, struct Het *het, double **D, double *params, RealVector *freq, double **SN, int N, double Tobs)
{
    int i, j, MM, MN;
    double logL;
    double **twave;
    double **res, **rwave, **dh;
    double *rc, *rs;
    double HH, HHS, HD, x;
    int imin, imax;
    
    MM = het->MM;
    MN = het->MN;
    
        twave = double_matrix(net->Nifo,N);
        rwave = double_matrix(net->Nifo,N);
        res = double_matrix(net->Nifo,N);
        dh = double_matrix(net->Nifo,N);
    
        fulltemplates(net, twave, freq, params, N);
        fulltemplates(net, rwave, freq, het->pref, N);
        
        for(i = 0; i < net->Nifo; i++)
        {
          for(j = 0; j < N; j++)
          {
               res[i][j] = D[i][j] - rwave[i][j];
               dh[i][j] =  rwave[i][j]-twave[i][j];
          }
        }
        
        logL = het->logLb;
        for(i = 0; i < net->Nifo; i++)
        {
            HH = fourier_nwip(dh[i], dh[i], SN[i], MN, MM, N);
            HD = fourier_nwip(res[i], dh[i], SN[i], MN, MM, N);
            logL -= (HD+0.5*HH);
            // printf("Test %d %f %f\n", i, HH, HD);
        }
        
        free_double_matrix(dh,net->Nifo);
        free_double_matrix(twave,net->Nifo);
        free_double_matrix(rwave,net->Nifo);
        free_double_matrix(res,net->Nifo);
    
    return(logL);
    
}



// Glitch-safe likelihood maximized wrt to amoplitude and phase
double log_likelihood_APmax(struct Net *net, double **D, double *params, RealVector *freq, double **SN, int N, double Tobs)
{
    
    double f, logL;
    int i, id, flag;
    int imin, imax;
    double *pref;
    double *A, *P;
    double **twave;
    
    logL = 0.0;
    
     if(lhold == 0)
     {

         twave = double_matrix(net->Nifo,N);
         templates(net, twave, freq, params, N);
         A = double_vector(net->Nifo);
         P = double_vector(net->Nifo);
         
         for(id = 0; id < net->Nifo; id++)
         {
             logL += log_likelihood_penalized(params, D[id], twave[id], SN[id], &A[id], &P[id], N, Tobs, 0);
            // printf("%d %f %f\n", id, A[id], P[id]);
         }
         
        // printf("%e %e\n", params[6], params[4]);
         
         params[6] -= log(A[0]);
         params[4] -= P[0]/2.0;  // PhenomD uses orbital phase, while here we have GW phase
         if(params[4] < 0.0) params[4] += PI;
         if(params[4] > PI) params[4] -= PI;
         
        // printf("%e %e\n", params[6], params[4]);
         
        // printf("%e %e\n", params[NX], params[NX+2]);
         
         for(id = 1; id < net->Nifo; id++)
            {
            params[NX+(id-1)*3] -= P[id]-P[0]; // GW phase offset
            params[NX+(id-1)*3+2] *= A[id]/A[0]; // amplitude ratio
            if(params[NX+(id-1)*3] < 0.0) params[NX+(id-1)*3] += TPI;
            if(params[NX+(id-1)*3] > TPI) params[NX+(id-1)*3] -= TPI;
            }
         
         flag = 0;
         for(id = 0; id < net->Nifo; id++)
         {
             if(A[id] > 1.1 || A[id] < 0.9) flag = 1;
             if(fabs(P[id]) > 0.1) flag = 1;
         }
        // printf("%e %e\n", params[NX], params[NX+2]);
         
         //printf("flag %d\n", flag);
         
         // one more rinse if the parameters shifted a lot. The glitch-safe banding can be sensitive to
         // the starting point
         if(flag == 1)
         {
             logL = 0.0;
             templates(net, twave, freq, params, N);
             
             for(id = 0; id < net->Nifo; id++)
              {
                  logL += log_likelihood_penalized(params, D[id], twave[id], SN[id], &A[id], &P[id], N, Tobs, 0);
                 // printf("%d %f %f\n", id, A[id], P[id]);
              }

              params[6] -= log(A[0]);
              params[4] -= P[0]/2.0;
              if(params[4] < 0.0) params[4] += PI;
              if(params[4] > PI) params[4] -= PI;

              for(id = 1; id < net->Nifo; id++)
                 {
                 params[NX+(id-1)*3] -= P[id]-P[0]; // GW phase offset
                 params[NX+(id-1)*3+2] *= A[id]/A[0]; // amplitude ratio
                 if(params[NX+(id-1)*3] < 0.0) params[NX+(id-1)*3] += TPI;
                 if(params[NX+(id-1)*3] > TPI) params[NX+(id-1)*3] -= TPI;
                 }
         }
         
         free_double_vector(A);
         free_double_vector(P);
         free_double_matrix(twave,net->Nifo);
         
     }
    
    return(logL);
    
}


double log_likelihood_max(struct Net *net, double **D, double *params, RealVector *freq, double **SN, int N, double Tobs, double tmin, double tmax, int pflag)
{
    int i, j, k, kk, NMAX, M;
    int ii, jj, id, flag;
    int imin, imax;
    int mn, mx;
    double sum;
    double *HH;
    double logL;
    double *lmax;
    double Ax, Px;
    double HD, LD, dt, x, f;
    double *norm, *delt, *pshift;
    double **HS, *HLS, **LX;
    double **LL, **TD;
    double **LLU, **TDU;
    double **PU, **AU;
    double **maxL, ***tx;
    double ***px, ***ax;
    int **shift;
    int **max, *tk;
    int *NU;
    double **HC, **HF;
    double *twave;
    
    double nm, psh, tsh, ph, cph, sph;
    
    double t0, t1, t2;
    
    double *h;
    
    logL = 0.0;
    
    
    if(lhold == 0)
    {
        
   
        
    f = fbegin(params);
        
    if(f < fmax)
    {
           
           if(fmin > f)
           {
           imin = (int)(Tobs*fmin);
           }
           else
           {
           imin = (int)(Tobs*f);
           }
        
    imin = (int)(Tobs*fmin);
    imax = (int)(Tobs*fmax);
    if(imin < 1) imin = 1;
    if(imax > N/2) imax = N/2;
        
    mn = (int)(tmin/Tobs*(double)(N))-1;
    mx = (int)(tmax/Tobs*(double)(N))+1;
    if(mn < 0) mn = 0;
    if(mx > N) mx = N;
        
       
        
        M = mx-mn;
        LL = double_matrix(net->Nifo,M);
        TD = double_matrix(net->Nifo,M);
        
        LLU = double_matrix(net->Nifo,M);
        TDU = double_matrix(net->Nifo,M);
        PU = double_matrix(net->Nifo,M);
        AU = double_matrix(net->Nifo,M);
        
        NU = int_vector(net->Nifo);

        maxL = double_matrix(net->Nifo,net->Nifo);
        
        tx = double_tensor(net->Nifo,net->Nifo,net->Nifo);
        px = double_tensor(net->Nifo,net->Nifo,net->Nifo);
        ax = double_tensor(net->Nifo,net->Nifo,net->Nifo);

    
    h = double_vector(N);
    twave = double_vector(N);
    HH = double_vector(net->Nifo);
    HS = double_matrix(net->Nifo,N);
    HC = double_matrix(net->Nifo,N);
    HF = double_matrix(net->Nifo,N);
    LX = double_matrix(net->Nifo,N);
    lmax = double_vector(N);
    tk = int_vector(net->Nifo);
    max = int_matrix(net->Nifo,N); // holds reference offset for each detector
    delt = double_vector(net->Nifo);
    norm = double_vector(net->Nifo);
    pshift = double_vector(net->Nifo);
    
    dt = Tobs/(double)(N);
 
    // We always shift reference time to center
    params[5] = Tobs/2.0;
    
    // reference intrinsic waveform
        PDwave(h, freq, net->tmax, params, N);
        
    for(k = 0; k < net->Nifo; k++) HH[k] = fourier_nwip(h, h, SN[k], imin, imax, N);
        
    for(k = 0; k < net->Nifo; k++)
    {
        pbt_shift(HC[k], HF[k], D[k], h, SN[k], imin, imax, N);
        for(i = 0; i < N; i++) HS[k][i] = 0.0;
        for(i = 0; i < N/2; i++) HS[k][i+N/2] += sqrt(HC[k][i]*HC[k][i]+HF[k][i]*HF[k][i]);
        for(i = -N/2; i < 0; i++) HS[k][i+N/2] += sqrt(HC[k][N+i]*HC[k][N+i]+HF[k][N+i]*HF[k][N+i]);
    }
        
    // Allow time-travel time difference between detectors
    for(id = 1; id < net->Nifo; id++)  tk[id] = (int)((net->tds[id])/dt);
        
    gsl_vector *v1 = gsl_vector_alloc (M);
    gsl_vector *v2 = gsl_vector_alloc (M);
        
    
    for(k = 0; k < net->Nifo; k++)
    {
        NU[k] = -1;
        
      j = 0;
      for (i = mn; i < mx; ++i)
       {
            TD[k][j] = dt*(double)(i-N/2);
            HD = 2.0*(double)(N)*HS[k][i];
            LL[k][j] = (HD*HD/HH[k])/2.0;
            j++;
        }
        
    for(j=0; j< M; j++)
    {
     gsl_vector_set(v1, j, LL[k][j]);
     gsl_vector_set(v2, j, TD[k][j]);
    }
        
    gsl_sort_vector2(v1,v2);
        
    for(j=0; j< M; j++)
    {
        LL[k][j] = gsl_vector_get(v1,M-1-j);
        TD[k][j] = gsl_vector_get(v2,M-1-j);
    }
        
    
    for(j=0; j< M; j++)
    {
        if(j == 0)
        {
        NU[k]++;
        LLU[k][NU[k]] = LL[k][j];
        TDU[k][NU[k]] = TD[k][j];
        }
        else
        {
          flag = 0;
            
          if(LL[k][j] >  Lpmin)  // check if loud enough
          {
            for (i = 0; i <= NU[k]; ++i) // check if a new peak
            {
                if(fabs(TD[k][j]-TDU[k][i]) < dtp) flag = 1;
            }
          }
          else
          {
              flag = 1;
          }
        
          if(flag == 0)
          {
              NU[k]++;
              LLU[k][NU[k]] = LL[k][j];
              TDU[k][NU[k]] = TD[k][j];
          }
            
        }
        
    }

    }
        
    gsl_vector_free(v1);
    gsl_vector_free(v2);
        
        
     for(k = 0; k < net->Nifo; k++)
     {
        
     for (i = 0; i <= NU[k]; ++i)
        {
            
         j = (int)(rint(TDU[k][i]/dt))+N/2;
         HD = 2.0*(double)(N)*HS[k][j];
         delt[k] = TDU[k][i];
            
         if(j >= N/2)
          {
         PU[k][i]  = atan2(HF[k][j-N/2],HC[k][j-N/2]);
         }
         else
         {
          PU[k][i]  = atan2(HF[k][j+N/2],HC[k][j+N/2]);
         }
         AU[k][i] = HD/HH[k];
            
          x = LLU[k][i];
            
            nm = AU[k][i];
            psh = PU[k][i];
            tsh = TDU[k][i];
            
            // apply the phase, time and amplitude shifts
            for(ii = imin; ii < imax; ii++)
            {
              f = (double)(ii)/Tobs;
              ph = psh-TPI*f*tsh;
              // could use recursion to make this faster
              cph = cos(ph);
              sph = sin(ph);
              twave[ii] = nm*(h[ii]*cph-h[N-ii]*sph);
              twave[N-ii] = nm*(h[N-ii]*cph+h[ii]*sph);
            }
            
          LLU[k][i] = log_likelihood_penalized(params, D[k], twave, SN[k], &Ax, &Px, N, Tobs, pflag);
            
            AU[k][i] *= Ax;
            PU[k][i] -= Px;
            
          if(pflag == 1) printf("%d %d %f %f %f\n", k, i,  x, LLU[k][i], TDU[k][i]+Tobs/2);
        
            
        }
        
     }
     
        
    // need to do a network maximization if more than 1 detector
    if(net->Nifo > 1)
    {
        
        
    for(i = 0; i < net->Nifo; i++)
    {
     for(j = i+1; j < net->Nifo; j++)
     {
         maxL[i][j] = 0.0;
         
         for(k=0; k<= NU[i]; k++)
         {
             if(LLU[i][k] > Lpmin)
             {
                 for(kk=0; kk <= NU[j]; kk++)
                   {
                      if(LLU[j][kk] > Lpmin)
                      {
                        if(fabs(TDU[i][k]-TDU[j][kk]) < net->delays[i][j])
                        {
                            x = LLU[i][k] + LLU[j][kk];
                            if(x > maxL[i][j])
                            {
                                maxL[i][j] = x;
                                tx[i][j][i] = TDU[i][k];
                                px[i][j][i] = PU[i][k];
                                ax[i][j][i] = AU[i][k];
                                tx[i][j][j] = TDU[j][kk];
                                px[i][j][j] = PU[j][kk];
                                ax[i][j][j] = AU[j][kk];
                            }
                        }
                      }
                   }
             }
         }
       }
    }
       

    // Note: there may not be a sky location that allows for the time delays between all detector pairs
    // the skyfind subroutine does its best to find something close to an unphysical maximum
        

        
        // if more than two detectors fold in the additional power
        if(net->Nifo > 2)
        {
            
         for(i = 0; i < net->Nifo; i++)
          {
          for(j = i+1; j < net->Nifo; j++)
           {
             if(maxL[i][j] > 0.0)
             {
                t0 = tx[i][j][i];
                t1 = tx[i][j][j];
                 
                 
                 for(k = 0; k < net->Nifo; k++)
                 {
                    x = 0.0;
                     
                    if(k != i && k != j)
                    {
                       
                       // initialize just in case no good value found
                        tx[i][j][k] = (t0+t1)/2.0;
                        px[i][j][k] = 0.0;
                        ax[i][j][k] = 1.0;
                        
                        for (ii = 0; ii < NU[k]; ++ii)
                        {
 
                            t2 = TDU[k][ii];
                            if(fabs(t0-t2) < net->delays[i][k] && fabs(t1-t2) < net->delays[j][k])
                              {
                                if(LLU[k][ii] > x)
                                {
                                    x = LLU[k][ii];
                                    tx[i][j][k] = TDU[k][ii];
                                    px[i][j][k] = PU[k][ii];
                                    ax[i][j][k] = AU[k][ii];
                                }
                             }
                        }
                        
                      }
                      maxL[i][j] += x;
                 }
                 
              }
            }
           }
          }
        
        
        // find the overall maximum
         x = 0.0;
         for(i = 0; i < net->Nifo; i++)
            {
            for(j = i+1; j < net->Nifo; j++)
                {
                    if(maxL[i][j] > x)
                    {
                        x = maxL[i][j];
                        ii = i;
                        jj = j;
                    }
                }
            }
        
        if(x > 0.0)
        {
          logL = maxL[ii][jj];
          for(k = 0; k < net->Nifo; k++)
          {
            delt[k] = tx[ii][jj][k];
            pshift[k] = px[ii][jj][k];
            norm[k] = ax[ii][jj][k];
          }
        }
        else
        {
            logL = 0.0;
            for(k = 0; k < net->Nifo; k++)
            {
                norm[k] = 1.0;
                delt[k] = 0.0;
                pshift[k] = 0.0;
            }
        }

    } // end multi-detector section
    else
    {
        // for a single detector
        
        x = 0.0;
        for(k=0; k<= NU[0]; k++)
        {
            if(LLU[0][k] > x)
            {
                x = LLU[0][k];
                jj = k;
            }
            
        }
        
        if(x > 0.0)
        {
        delt[0] = TDU[0][jj];
        norm[0] = AU[0][jj];
        pshift[0] = PU[0][jj];
        logL = LLU[0][jj];
        }
        else
        {
         norm[0] = 1.0;
         delt[0] = 0.0;
         pshift[0] = 0.0;
         logL = 0.0;
        }
        
            
    }

    
    // Shift the waverform to match in reference detector
    
    params[6] -= log(norm[0]);
    params[5] += delt[0];
    params[4] += pshift[0]/2.0;  // PhenomD uses orbital phase, while here we have GW phase
    
    if(params[4] < 0.0) params[4] += PI;
    if(params[4] > PI) params[4] -= PI;
    
    for(k = 1; k < net->Nifo; k++)
    {
    
    params[NX+(k-1)*3] = pshift[k]-pshift[0]; // GW phase offset
    params[NX+(k-1)*3+1] = delt[k]-delt[0]; // time offset
    params[NX+(k-1)*3+2] = norm[k]/norm[0]; // amplitude ratio
    if(params[NX+(k-1)*3] < 0.0) params[NX+(k-1)*3] += TPI;
    if(params[NX+(k-1)*3] > TPI) params[NX+(k-1)*3] -= TPI;
        
    }
        
        
       
    free_double_matrix(maxL,net->Nifo);
        
    free_double_tensor(tx,net->Nifo,net->Nifo);
    free_double_tensor(px,net->Nifo,net->Nifo);
    free_double_tensor(ax,net->Nifo,net->Nifo);
        
    free_double_matrix(LL,net->Nifo);
    free_double_matrix(TD,net->Nifo);
        
    free_double_matrix(LLU,net->Nifo);
    free_double_matrix(TDU,net->Nifo);
    free_double_matrix(PU,net->Nifo);
    free_double_matrix(AU,net->Nifo);
    free_int_vector(NU);
        
    free_double_matrix(LX,net->Nifo);
    free_double_matrix(HS,net->Nifo);
    free_double_matrix(HC,net->Nifo);
    free_double_matrix(HF,net->Nifo);
    free_double_vector(lmax);
    free_int_vector(tk);
    free_int_matrix(max,net->Nifo);
    free_double_vector(delt);
    free_double_vector(norm);
    free_double_vector(pshift);
    free_double_vector(HH);
    free_double_vector(h);
    free_double_vector(twave);
        
    }
    }
    
    return logL;
}

double log_likelihood_penalized(double *params, double *D, double *twave, double *SN, double *Ax, double *Px, int N, double Tobs, int pflag)
{
    int i, j, k;
    double logL, f, t;
    double HHX, HCX, HSX;
    double hh, dh, dhs;
    double *rho, *snr, *rhoc;
    double HH, HD, HC, HS;
    double ph, cph, sph, x, xs, y, fs;
    double max, min, alpha;
    double maxL, minL;
    int imin, imax;
    int *fbs, *bs;
    
    int bw, nb, ib, kend, flag, kstart;
    
    bw = (int)(fbw*Tobs);
    
    logL = 0.0;
    
    if(lhold == 0)
    {
        
        fs = fbegin(params);
               
               if(fmin > fs)
               {
               imin = (int)(Tobs*fmin);
               }
               else
               {
               imin = (int)(Tobs*fs);
               }
    
        imax = (int)(Tobs*fmax);
        if(imax > N/2) imax = N/2;
        nb = (imax-imin)/bw;
        if(nb < 1) nb = 1;
        
        rho = double_vector(nb);
        rhoc = double_vector(nb);
        snr = double_vector(nb);
        fbs = int_vector(nb);
        bs = int_vector(nb);
        
        for(i = 0; i < nb; i++) rho[i] = 0.0;
        
        flag = 0;
        j = 0;
        k = -1;
        x = 0.0;
        xs = 0.0;
        y = 0.0;
        
        for(i = imin; i < imax; i++)
        {
            if(j%bw==0) // moving to next block
            {
                if(y > 0.0)
                {
                x = sqrt(x*x+xs*xs); // phase maximized over frequency band
                rho[k] = (x-y)/sqrt(y);
                }
                
                if(k > -1) snr[k] = y;

                k++;
                x = 0.0;
                xs = 0.0;
                y = 0.0;
            }
            
            dh =  4.0*(D[i]*twave[i]+D[N-i]*twave[N-i])/SN[i];
            dhs = 4.0*(D[i]*twave[N-i]-D[N-i]*twave[i])/SN[i];
            hh =  4.0*(twave[i]*twave[i]+twave[N-i]*twave[N-i])/SN[i];
    
            HD += dh;
            HH += hh;
            x += dh;
            xs += dhs;
            y += hh;
            
            j++;
        }
        
        // find the region where the SNR is significant
        flag = 0;
        kstart = 0;
        for(k = 0; k < nb; k++)
        {
            if(snr[k] > 0.1 && flag == 0)
            {
                flag = 1;
                kstart = k;
            }
        }
        
        flag = 0;
        kend = nb;
        for(k = kstart; k < nb; k++)
        {
            if(snr[k] < 0.1 && flag == 0)
            {
                flag = 1;
                kend = k;
            }
        }
        
        // now cut out any low SNR bands
        j = 0;
        for(k = kstart; k < kend; k++)
        {
                fbs[j] = imin + k*bw;
                bs[j] = 0;
                rhoc[j] = rho[k];
                j++;
        }
        
           for(k = 0; k < j; k++)
            {
                if(rhoc[k] > rhocut) bs[k] = 1;  // exclude this band
            }
        
        // re-do amplitude and phase maximization with bad bands excluded
        HH = 0.0001;
        HC = 0.0;
        HS = 0.0;
       
        for(k = 0; k < j; k++)
         {
             
             if(pflag == 1)
             {
                 f = ((double)(fbs[k])+0.5*(double)(bw))/Tobs;
                 t = t_at_f(params, f);
                 HHX = 0.0;
                 HCX = 0.0;
                 HSX = 0.0;
                 for(i = fbs[k]; i < (fbs[k]+bw); i++)
                 {
                   HHX +=  4.0*(twave[i]*twave[i]+twave[N-i]*twave[N-i])/SN[i];
                   HCX +=  4.0*(D[i]*twave[i]+D[N-i]*twave[N-i])/SN[i];
                   HSX +=  4.0*(D[i]*twave[N-i]-D[N-i]*twave[i])/SN[i];
                 }
                 HD = sqrt(HCX*HCX+HSX*HSX);
                 printf("%f %f %f\n", f, rhoc[k], HD-0.5*HHX);
             }
             
             if(bs[k] == 0)
             {
                 for(i = fbs[k]; i < (fbs[k]+bw); i++)
                 {
                   HH +=  4.0*(twave[i]*twave[i]+twave[N-i]*twave[N-i])/SN[i];
                   HC +=  4.0*(D[i]*twave[i]+D[N-i]*twave[N-i])/SN[i];
                   HS +=  4.0*(D[i]*twave[N-i]-D[N-i]*twave[i])/SN[i];
                 }
             }
         }

        HD = sqrt(HC*HC+HS*HS);
        *Ax = HD/HH;
        *Px = atan2(HS,HC);
        
        logL = (HD*HD/HH)/2.0;
        
     free_int_vector(fbs);
     free_int_vector(bs);
     free_double_vector(rho);
     free_double_vector(rhoc);
     free_double_vector(snr);
        
    }
    
    return(logL);
    
}

void wavemax(struct Net *net, int N, double **tp, RealVector *freq, double **paramx, int *who, double **SN, double Tobs)
{
    
    double **twave, **tf;
    double **ww, **wf;
    double *A, *AW;
    double fr, f, x, y, dt, fac;
    int kx, ky;
    int i, j, id;
    char command[1024];
    FILE *out;
    
    dt = Tobs/(double)(N);
    
    A = double_vector(N);
    AW = double_vector(N);
    
    twave = double_matrix(net->Nifo,N);
    tf = double_matrix(net->Nifo,N);
    ww = double_matrix(net->Nifo,N);
    wf = double_matrix(net->Nifo,N);
    
    j = who[1];
    
    fulltemplates(net, twave, freq, paramx[j], N);
    
    fr = fmin/8.0;
    
    for (i = 1; i < N/2; ++i)
    {
        f = freq->data[i];
        x = 0.5*(1.0+tanh((f-(fmin+2.0*fr))/fr));
        for (id = 0; id < net->Nifo; ++id)
        {
            y = 1.0/sqrt(SN[id][i]);
            
            twave[id][i] *= x;
            twave[id][N-i] *= x;
            tf[id][i] = twave[id][N-i];
            tf[id][N-i] = -twave[id][i];
            
            ww[id][i] = y*twave[id][i];
            ww[id][N-i] = y*twave[id][N-i];
            wf[id][i] = ww[id][N-i];
            wf[id][N-i] = -ww[id][i];
            
        }
    }
    

    
    for (id = 0; id < net->Nifo; ++id)
    {
        twave[id][0] = 0.0;
        twave[id][N/2] = 0.0;
        tf[id][0] = 0.0;
        tf[id][N/2] = 0.0;
        ww[id][0] = 0.0;
        ww[id][N/2] = 0.0;
        wf[id][0] = 0.0;
        wf[id][N/2] = 0.0;
        
        gsl_fft_halfcomplex_radix2_inverse(twave[id], 1, N);
        gsl_fft_halfcomplex_radix2_inverse(tf[id], 1, N);
        gsl_fft_halfcomplex_radix2_inverse(ww[id], 1, N);
        gsl_fft_halfcomplex_radix2_inverse(wf[id], 1, N);
        
        x = 0.0;
        y = 0.0;
        for (i = 0; i < N; ++i)
        {
            A[i] = (twave[id][i]*twave[id][i]+tf[id][i]*tf[id][i]);
            AW[i] = ww[id][i]*ww[id][i]+wf[id][i]*wf[id][i];
            if(A[i] > x)
            {
                x = A[i];
                kx = i;
            }
            if(AW[i] > y)
            {
                y = AW[i];
                ky = i;
            }
        }
        
        // basic estimate of the peak time
        tp[id][0] = (double)(kx)*dt;
        tp[id][1] = (double)(ky)*dt;
        
        // refined estimate of the peak time
        if(kx > 1 && kx < N-1)  tp[id][0] = (double)(kx)*dt + 0.5*(A[kx+1]-A[kx-1])/(2.0*A[kx] -A[kx+1]-A[kx-1])*dt;
        if(ky > 1 && ky < N-1)  tp[id][1] = (double)(ky)*dt + 0.5*(AW[ky+1]-AW[ky-1])/(2.0*AW[ky] -AW[ky+1]-AW[ky-1])*dt;
 
    }

    free_double_vector(A);
    free_double_vector(AW);
    free_double_matrix(twave,net->Nifo);
    free_double_matrix(tf,net->Nifo);
    free_double_matrix(ww,net->Nifo);
    free_double_matrix(wf,net->Nifo);
    
}



void printwaveall(struct Net *net, int N, RealVector *freq, double *paramx, double **SN, double Tobs, double ttrig, int iter)
{
    
    double **twave, **tf, **thold;
    double **ww, **wf;
    double fr, f, x, y, dt, fac, fs;
    int i, j, id;
    char command[1024];
    FILE *out;
    
    dt = Tobs/(double)(N);
    
    thold = double_matrix(net->Nifo,N);
    twave = double_matrix(net->Nifo,N);
    tf = double_matrix(net->Nifo,N);
    ww = double_matrix(net->Nifo,N);
    wf = double_matrix(net->Nifo,N);
    
    fulltemplates(net, twave, freq, paramx, N);
    
    fs = fbegin(paramx);
    if(fs < fmin) fs = fmin;
    
    printf("%d %f %f\n", iter, exp(paramx[0])/MSUN_SI, fs);
    
    for (id = 0; id < net->Nifo; ++id)
    {
        thold[id][0] = 0.0;
        thold[id][N/2] = 0.0;
        for (i = 0; i < N; ++i) thold[id][i] = twave[id][i];
    }
    
    fr = fs/3.0;
    
    for (i = 1; i < N/2; ++i)
    {
        
        f = freq->data[i];
        x = 0.5*(1.0+tanh((f-(fs+fr))/fr));
        for (id = 0; id < net->Nifo; ++id)
        {
            y = 1.0/sqrt(SN[id][i]);
            
            twave[id][i] *= x;
            twave[id][N-i] *= x;
            tf[id][i] = twave[id][N-i];
            tf[id][N-i] = -twave[id][i];
            
            ww[id][i] = y*twave[id][i];
            ww[id][N-i] = y*twave[id][N-i];
            wf[id][i] = ww[id][N-i];
            wf[id][N-i] = -ww[id][i];
            
        }
    }
    
    x = sqrt(2.0)/dt;
    
    fac = sqrt(Tobs)*sqrt((double)(N/2));
    
    for (id = 0; id < net->Nifo; ++id)
    {
        twave[id][0] = 0.0;
        twave[id][N/2] = 0.0;
        tf[id][0] = 0.0;
        tf[id][N/2] = 0.0;
        ww[id][0] = 0.0;
        ww[id][N/2] = 0.0;
        wf[id][0] = 0.0;
        wf[id][N/2] = 0.0;
        
        gsl_fft_halfcomplex_radix2_inverse(twave[id], 1, N);
        gsl_fft_halfcomplex_radix2_inverse(tf[id], 1, N);
        gsl_fft_halfcomplex_radix2_inverse(ww[id], 1, N);
        gsl_fft_halfcomplex_radix2_inverse(wf[id], 1, N);
        
        for (i = 0; i < N; ++i)
        {
            twave[id][i] *= x;
            tf[id][i] *= x;
            ww[id][i] *= fac;
            wf[id][i] *= fac;
        }
    }
    
    for (id = 0; id < net->Nifo; ++id)
    {
        sprintf(command, "waves/wave_%d_%d_%d_%d.dat", iter, (int)(Tobs), (int)ttrig, net->labels[id]);
        out = fopen(command,"w");
        for (i = 0; i < N; ++i)
        {
            fprintf(out,"%e %.16e %e\n", (double)(i)*dt-Tobs+2.0, twave[id][i], sqrt(twave[id][i]*twave[id][i]+tf[id][i]*tf[id][i]));
        }
        fclose(out);
    }
    
    for (id = 0; id < net->Nifo; ++id)
    {
        sprintf(command, "waves/wavef_%d_%d_%d_%d.dat", iter, (int)(Tobs), (int)ttrig, net->labels[id]);
        out = fopen(command,"w");
        for (i = 1; i < N/2; ++i)
        {
            f = (double)(i)/Tobs;
            if(f > fs)
            {
            fprintf(out,"%e %.16e %.16e\n", f, thold[id][i], thold[id][N-i]);
            }
            else
            {
              fprintf(out,"%e %.16e %.16e\n", f, 0.0, 0.0);
            }
        }
        fclose(out);
    }
    
    
    for (id = 0; id < net->Nifo; ++id)
    {
        sprintf(command, "waves/wavewhite_%d_%d_%d_%d.dat", iter, (int)(Tobs), (int)ttrig, net->labels[id]);
        out = fopen(command,"w");
        for (i = 0; i < N; ++i)
        {
            fprintf(out,"%e %e %e\n", (double)(i)*dt-Tobs+2.0, ww[id][i], sqrt(ww[id][i]*ww[id][i]+wf[id][i]*wf[id][i]));
        }
        fclose(out);
    }
    
    free_double_matrix(thold,net->Nifo);
    free_double_matrix(twave,net->Nifo);
    free_double_matrix(tf,net->Nifo);
    free_double_matrix(ww,net->Nifo);
    free_double_matrix(wf,net->Nifo);
    
}


void printwave(struct Net *net, int N, RealVector *freq, double **paramx, int *who, double Tobs, double ttrig, int iter)
{
    
    double **twave;
    double fr, f, x, dt;
    int i, j, id;
    char command[1024];
    FILE *out;
    
    dt = Tobs/(double)(N);
    
    twave = double_matrix(net->Nifo,N);
    
    
    
    j = who[1];
    templates(net, twave, freq, paramx[j], N);
    
    printf("%d %f\n", iter, exp(paramx[j][0])/MSUN_SI);
    
    fr = fmin/8.0;
    
    for (i = 1; i < N/2; ++i)
    {
        f = freq->data[i];
        x = 0.5*(1.0+tanh((f-(fmin+2.0*fr))/fr));
        for (id = 0; id < net->Nifo; ++id)
        {
            twave[id][i] *= x;
            twave[id][N-i] *= x;
        }
    }
    
    x = sqrt(2.0)/dt;
    
    for (id = 0; id < net->Nifo; ++id)
    {
        gsl_fft_halfcomplex_radix2_inverse(twave[id], 1, N);
        for (i = 0; i < N; ++i)
        {
            twave[id][i] *= x;
        }
    }
    
    for (id = 0; id < net->Nifo; ++id)
    {
    sprintf(command, "wave_%d_%d_%d_%d.dat", iter, (int)(Tobs), (int)ttrig, net->labels[id]);
    out = fopen(command,"w");
    for (i = 0; i < N; ++i)
    {
        fprintf(out,"%e %e\n", (double)(i)*dt, twave[id][i]);
    }
    fclose(out);
    }
    
    free_double_matrix(twave,net->Nifo);
    
}

void fisherskysetup(struct Net *net, double **wave, double **HH, double Tobs, int n)
{
    double f;
    int i, j, k, l, bn, ii, id;
    double **wavef;
    
    wavef = double_matrix(net->Nifo,n);
    
    for(id = 0; id<net->Nifo; id++)
    {
        wavef[id][0] = 0.0;
        wavef[id][n/2] = 0.0;
        for(i = 1; i<n/2; i++)
        {
            f = (double)(i)/Tobs;
            j = i;
            k = n-i;
            wavef[id][j] = f*wave[id][j];
            wavef[id][k] = f*wave[id][k];
        }
    }
    
    for(id = 0; id<net->Nifo; id++)
    {
        HH[id][0] =  4.0*f_nwip(wave[id], wave[id], n);
        HH[id][1] =  4.0*f_nwip(wave[id], wavef[id], n);
        HH[id][2] =  4.0*f_nwip(wavef[id], wavef[id], n);
    }
    
    free_double_matrix(wavef,net->Nifo);
    
}





void skylikesetup(struct Net *net, double **data,  double **wave, double *D, double *H, double **DHc,  double **DHs, double Tobs, int n, int bn, int nt, int imin, int imax)
{
    double *corr, *corrf;
    double dt;
    int i, j, k, kk, l, ii, id;
    
    dt = Tobs/(double)(bn);
    
    corr = double_vector(bn);
    corrf = double_vector(bn);
    
    for(id = 0; id<net->Nifo; id++)
    {
        
        D[id] = 4.0*f_nwip(data[id], data[id], n);
        H[id] = 4.0*f_nwip(wave[id], wave[id], n);
        
        // The arrays are zero padded to increase the sample cadence
        
        for (i=0; i < bn; i++)
        {
            corr[i]    = 0.0;
            corrf[i] = 0.0;
        }
        
        // the non-zero part
        for (i=1; i < n/2; i++)
        {
            l=i;
            k=bn-i;
            kk = n-i;
            corr[l]    = ( data[id][l]*wave[id][l] + data[id][kk]*wave[id][kk]);
            corr[k]    = ( data[id][kk]*wave[id][l] - data[id][l]*wave[id][kk]);
            corrf[l] = corr[k];
            corrf[k] = -corr[l];
        }
        
        
        gsl_fft_halfcomplex_radix2_inverse(corr, 1, bn);
        gsl_fft_halfcomplex_radix2_inverse(corrf, 1, bn);
        
        
        k = (nt-1)/2;
        
        for(i=-k; i<k; i++)
        {
            j = i+k;
            
            if(i < 0)
            {
                DHc[id][j] = 2.0*(double)(bn)*corr[bn+i];
                DHs[id][j] = 2.0*(double)(bn)*corrf[bn+i];
            }
            else
            {
                DHc[id][j] = 2.0*(double)(bn)*corr[i];
                DHs[id][j] = 2.0*(double)(bn)*corrf[i];
            }
            
        }
        
        
    }
    
    free_double_vector(corr);
    free_double_vector(corrf);
    
    
    
}

double skylike(struct Net *net, double *params, double *D, double *H, double **DHc,  double **DHs, double dt, int nt, int flag)
{
    int i, j, k, l, id;
    double alpha, sindelta, dphi, t0, A, FA, ecc, ciota, psi;
    double *dtimes, *F, *lambda;
    double tdelay, tx, toff, dc, ds, DH;
    double dcc, dss;
    double Fplus, Fcross;
    double Ap, Ac;
    double clam, slam, x;
    double cphi, sphi;
    double logL;
    
    
    logL = 0.0;
    
    if(lhold == 0)
    {
        
        dtimes = (double*)malloc(sizeof(double)* 5);
        F = (double*)malloc(sizeof(double)* 5);
        lambda = (double*)malloc(sizeof(double)* 5);
        
        k = (nt-1)/2;
        
        alpha = params[0];
        sindelta = params[1];
        psi = params[2];
        ciota = params[3];
        
        Ap = (1.0+ciota*ciota)/2.0;
        Ac = -ciota;
        
        A = params[4];
        dphi = params[5];
        toff = params[6];
        
        cphi = cos(dphi);
        sphi = sin(dphi);
        
        Times(alpha, sindelta, net->GMST, dtimes);
        
        for(id = 0; id<net->Nifo; id++)
        {
            F_ant(psi, alpha, sindelta, net->GMST, &Fplus, &Fcross, net->labels[id]);
            F[id] = sqrt(Ap*Ap*Fplus*Fplus+Ac*Ac*Fcross*Fcross);  // magnitude of response
            lambda[id] = atan2(Ac*Fcross,Ap*Fplus);
            if(lambda[id] < 0.0) lambda[id] += TPI;
        }
        
        logL = 0.0;
        
        // Note that F[] and lambda[] have already used the label array to get the correct detector ordering while dtimes uses a fixed ordering
        // so needs the label array to get the correct delays
        
        for(id = 0; id<net->Nifo; id++)
        {
            // everything is reference to geocenter
            tdelay = toff+dtimes[net->labels[id]];
            
            i = (int)(floor(tdelay/dt));
            
            tx = (tdelay/dt - (double)(i));
            
           // if (flag == 1) printf("%d %d %d %f %f\n", i, i+k, nt, tx, tdelay/dt);
            
            // linear interpolation
            i += k;
            
            slam = sin(lambda[id]);
            clam = cos(lambda[id]);
            
            FA = A*F[id];
            
            if(i >= 0 && i < nt-2)
            {
                dc = DHc[id][i]*(1.0-tx)+DHc[id][i+1]*tx;
                ds = DHs[id][i]*(1.0-tx)+DHs[id][i+1]*tx;
                
                
                // put in phase rotation
                dcc = cphi*dc+sphi*ds;
                dss = -sphi*dc+cphi*ds;
                
                DH = clam*dcc+slam*dss;
                
                // relative likelihood
                x = -(FA*FA*H[id]-2.0*FA*DH)/2.0;
                
                // if (flag == 1) printf("%d %f %f\n", id, tdelay, x);
                
                logL += x;
                
            }
            else
            {
                logL -= 1.0e10;
            }
            
        }
        
        free(dtimes);
        free(F);
        free(lambda);
        
    }
    
    return(logL);
    
}

void fisher_matrix_fastsky(struct Net *net, double *params, double **fisher, double **HH)
{
    int i, j, k, l, id;
    double alpha, sindelta, phi0, t0, A, ecc, ciota, psi;
    double Ap, Ac;
    double dtimesP[4], dtimesM[4], dt0[4], dt1[4];
    double tdelay, tx, toff, dc, ds, DH;
    double dcc, dss;
    double Fplus, Fcross, epss;
    double FplusP, FcrossP, FplusM, FcrossM;
    double *dFplus, *dFcross;
    double *dlambda, *dtoff, *dF;
    double F, lambda, clam, slam, x, y, z;
    double cphi, sphi;
    double logL;
    
    dFplus = double_vector(NS);
    dFcross = double_vector(NS);
    dF = double_vector(NS);
    dlambda = double_vector(NS);
    dtoff = double_vector(NS);
    
    for(i = 3; i<NS; i++) dFplus[i] = 0.0;
    for(i = 3; i<NS; i++) dFcross[i] = 0.0;
    
    dF[5] = 0.0;
    dF[6] = 0.0;
    
    dlambda[4] = 0.0;
    dlambda[5] = 1.0;  // minus sign in phase definition
    dlambda[6] = 0.0;
    
    
    dtoff[2] = 0.0;
    dtoff[3] = 0.0;
    dtoff[4] = 0.0;
    dtoff[5] = 0.0;
    dtoff[6] = -1.0;
    
    alpha = params[0];
    sindelta = params[1];
    psi = params[2];
    ciota = params[3];
    A = params[4];
    phi0 = params[5];
    toff = params[6];
    
    Ap = 0.5*(1.0+ciota*ciota);
    Ac = -ciota;
    
    
    cphi = cos(phi0);
    sphi = sin(phi0);
    
    epss = 1.0e-6;
    
    Times(alpha+epss, sindelta, net->GMST, dtimesP);
    Times(alpha-epss, sindelta, net->GMST, dtimesM);
    for(id = 0; id<net->Nifo; id++)
    {
        dt0[id] = (dtimesP[net->labels[id]]-dtimesM[net->labels[id]])/(2.0*epss);
    }
    Times(alpha, sindelta+epss, net->GMST, dtimesP);
    Times(alpha, sindelta-epss, net->GMST, dtimesM);
    for(id = 0; id<net->Nifo; id++)
    {
        dt1[id] = (dtimesP[net->labels[id]]-dtimesM[net->labels[id]])/(2.0*epss);
    }
    
    for(i = 0; i<NS; i++)
    {
        for(j = 0; j<NS; j++)
        {
            fisher[i][j] = 0.0;
        }
    }
    
    for(id = 0; id<net->Nifo; id++)
    {
        
        dtoff[0] = -dt0[id];
        dtoff[1] = -dt1[id];
        
        //Calculate antenna pattern once for each detector
        //This is only OK because psi is constant over the orbits...not so if spin is added
        F_ant(psi, alpha, sindelta, net->GMST, &Fplus, &Fcross, net->labels[id]);
              
        F = A*sqrt(Ap*Ap*Fplus*Fplus+Ac*Ac*Fcross*Fcross);
        
        F_ant(psi, alpha+epss, sindelta, net->GMST, &FplusP, &FcrossP, net->labels[id]);
        F_ant(psi, alpha-epss, sindelta, net->GMST, &FplusM, &FcrossM, net->labels[id]);
        
        dFplus[0] = (FplusP-FplusM)/(2.0*epss);
        dFcross[0] = (FcrossP-FcrossM)/(2.0*epss);
        
        F_ant(psi, alpha, sindelta+epss, net->GMST, &FplusP, &FcrossP, net->labels[id]);
        F_ant(psi, alpha, sindelta-epss, net->GMST, &FplusM, &FcrossM, net->labels[id]);
        
        dFplus[1] = (FplusP-FplusM)/(2.0*epss);
        dFcross[1] = (FcrossP-FcrossM)/(2.0*epss);
        
        F_ant(psi+epss, alpha, sindelta, net->GMST, &FplusP, &FcrossP, net->labels[id]);
        F_ant(psi-epss, alpha, sindelta, net->GMST, &FplusM, &FcrossM, net->labels[id]);
        
        dFplus[2] = (FplusP-FplusM)/(2.0*epss);
        dFcross[2] = (FcrossP-FcrossM)/(2.0*epss);
        
        x = A*A/(F);
        dF[0] = x*(Ap*Ap*Fplus*dFplus[0]+Ac*Ac*Fcross*dFcross[0]);
        dF[1] = x*(Ap*Ap*Fplus*dFplus[1]+Ac*Ac*Fcross*dFcross[1]);
        dF[2] = x*(Ap*Ap*Fplus*dFplus[2]+Ac*Ac*Fcross*dFcross[2]);
        dF[3] = x*(Ap*ciota*Fplus*Fplus+Ac*Fcross*Fcross);
        dF[4] = F/A;
        
        x = A*A/(F*F);
        dlambda[0] = x*Ap*Ac*(Fplus*dFcross[0]-Fcross*dFplus[0]);
        dlambda[1] = x*Ap*Ac*(Fplus*dFcross[1]-Fcross*dFplus[1]);
        dlambda[2] = x*Ap*Ac*(Fplus*dFcross[2]-Fcross*dFplus[2]);
        dlambda[3] = x*Fcross*Fplus*(Ap-Ac*ciota);
        
        x = F*F;
        y = TPI*x;
        z = TPI*y;
        for(i = 0; i<NS; i++)
        {
            for(j = 0; j<NS; j++)
            {
                fisher[i][j] += HH[id][0]*(dF[i]*dF[j]+x*dlambda[i]*dlambda[j])+HH[id][1]*y*(dlambda[i]*dtoff[j]+dlambda[j]*dtoff[i])+HH[id][2]*z*dtoff[i]*dtoff[j];
            }
        }
        
        
    }
    
    /*
     printf("\n");
     for(i = 0; i<NS; i++)
     {
     for(j = 0; j<NS; j++)
     {
     printf("%e ", fisher[i][j]);
     }
     printf("\n");
     }*/
    
    free_double_vector(dFplus);
    free_double_vector(dFcross);
    free_double_vector(dF);
    free_double_vector(dlambda);
    free_double_vector(dtoff);
}



// finds the extrinsic parameters at the new sky location that preserve the waveforms
void skymap_all(double *paramsx, double *paramsy, double GMST, int ifo1, int ifo2)
{
    double phi, costheta, sintheta, cosdelta;
    double alpha, delta, sindelta, phi0, Phase, psi, psi0;
    double ecc, scale;
    double tl;
    double lAmin, lAmax, A;
    double Ap, Ac, ciota;
    double *fyplus, *fycross;
    double *fxplus, *fxcross;
    double ux, vx, wx, zx;
    double uy, vy, wy, zy;
    double x, y, xx, yy;
    double phiy, psiy, DLy, ciotay;
    int i, j, mc, ii, test;
    double *dtimesx, *dtimesy;
    
    int nfunk;
    
    fyplus = double_vector(5);
    fycross = double_vector(5);
    fxplus = double_vector(5);
    fxcross = double_vector(5);
    dtimesx = double_vector(5);
    dtimesy = double_vector(5);
    
    // [0] log(Mc) [1] log(Mt) [2] chi1 [3] chi2 [4] phi0 [5] tp [6] log(DL) [7] alpha [8] sindelta [9] psi [10] ciota
    
    // extract the primitives f+, fx for each detector at x
    alpha = paramsx[7];
    sindelta = paramsx[8];
    
    Times(alpha, sindelta, GMST, dtimesx);
    
    F_ant(0.0, alpha, sindelta, GMST, &fxplus[1], &fxcross[1], ifo1);
    F_ant(0.0, alpha, sindelta, GMST, &fxplus[2], &fxcross[2], ifo2);
    
    // extract the primitives f+, fx for each detector at y
    alpha = paramsy[7];
    sindelta = paramsy[8];
    
    Times(alpha, sindelta, GMST, dtimesy);
    
    F_ant(0.0, alpha, sindelta, GMST, &fyplus[1], &fycross[1], ifo1);
    F_ant(0.0, alpha, sindelta, GMST, &fyplus[2], &fycross[2], ifo2);
    
    uvwz_all(&ux, &vx, &wx, &zx, paramsx);
    
    uvwz_sol(&uy, &vy, &wy, &zy, ux, vx, wx, zx, fxplus[1], fyplus[1], fxcross[1], fycross[1], fxplus[2], fyplus[2], fxcross[2], fycross[2]);
    
    exsolve_all(&phiy, &psiy, &DLy, &ciotay, uy, vy, wy, zy);
    
    paramsy[9] = psiy;
    paramsy[10] = ciotay;
    paramsy[6] = log(DLy);
    paramsy[4] = phiy;
    
    paramsy[5] = paramsx[5]+dtimesx[ifo1]-dtimesy[ifo1];
    
    free_double_vector(dtimesx);
    free_double_vector(dtimesy);
    free_double_vector(fyplus);
    free_double_vector(fycross);
    free_double_vector(fxplus);
    free_double_vector(fxcross);
    
}


// finds the extrinsic parameters at the new sky location that preserve the waveforms
void skymap(double *paramsx, double *paramsy, double GMST, int ifo1, int ifo2, int iref)
{
    double phi, costheta, sintheta, cosdelta;
    double alpha, delta, sindelta, phi0, Phase, psi, psi0;
    double ecc, scale;
    double tl;
    double lAmin, lAmax, A;
    double *fyplus, *fycross;
    double *fxplus, *fxcross;
    double ux, vx, wx, zx;
    double uy, vy, wy, zy;
    double phiy, psiy, Ay, ciotay;
    double phix, psix, Ax, ciotax;
    int i, j, mc, ii, test;
    double Fp, Fc, Fs;
    double x, y, px, py;
    int nfunk;
    double *dtimesx, *dtimesy;
    
    fyplus = double_vector(5);
    fycross = double_vector(5);
    fxplus = double_vector(5);
    fxcross = double_vector(5);
    dtimesx = double_vector(5);
    dtimesy = double_vector(5);
    
    // extract the primitives f+, fx for each detector at x
    alpha = paramsx[0];
    sindelta = paramsx[1];
    F_ant(0.0, alpha, sindelta, GMST, &fxplus[1], &fxcross[1], ifo1);
    F_ant(0.0, alpha, sindelta, GMST, &fxplus[2], &fxcross[2], ifo2);
    
    Times(alpha, sindelta, GMST, dtimesx);
    
    // extract the primitives f+, fx for each detector at y
    alpha = paramsy[0];
    sindelta = paramsy[1];
    F_ant(0.0, alpha, sindelta, GMST, &fyplus[1], &fycross[1], ifo1);
    F_ant(0.0, alpha, sindelta, GMST, &fyplus[2], &fycross[2], ifo2);
    
    Times(alpha, sindelta, GMST, dtimesy);
    
    uvwz(&ux, &vx, &wx, &zx, paramsx);
    
    uvwz_sol(&uy, &vy, &wy, &zy, ux, vx, wx, zx, fxplus[1], fyplus[1], fxcross[1], fycross[1], fxplus[2], fyplus[2], fxcross[2], fycross[2]);
    
    exsolve(&phiy, &psiy, &Ay, &ciotay, uy, vy, wy, zy);
    
    paramsy[2] = psiy;
    paramsy[3] = ciotay;
    paramsy[4] = Ay;
    paramsy[5] = phiy;
    
    paramsy[6] = paramsx[6]+dtimesx[ifo1]-dtimesy[ifo1];
    
    //[0] alpha, [1] sin(delta) [2] psi [3] cos(iota) [4] scale [5] phi0 [6] dt
    
    
    free_double_vector(dtimesx);
    free_double_vector(dtimesy);
    free_double_vector(fyplus);
    free_double_vector(fycross);
    free_double_vector(fxplus);
    free_double_vector(fxcross);
    
}

double skydensity(double *paramsx, double *paramsy, double GMST, int ifo1, int ifo2, int iref)
{
    int i, ii, j;
    double *paramsxp, *paramsyp;
    double *Jmat;
    double x, Jack;
    double ep;
    
    ep = 1.0e-6;
    
    paramsxp = double_vector(NS);
    paramsyp = double_vector(NS);
    Jmat = double_vector(16);
    
    // Geocenter arrival time is also adjusted. Probably should be included in the Jacobian, currently isn't
    
    // At this stage all we have for the new location y is the sky location. So first we find what the rest
    // of the parameters (amplitude, polarization, ciota, phase) need to be to keep the signal the same
    
    skymap(paramsx, paramsy, GMST, ifo1, ifo2, iref);
    
    // The sky locations are chosen to be on the equal time delay sky ring. Presumably the density factors for that
    // mapping have been taken care of. What we need to know here is how the volume of a small volume
    // dAx /\ dpsix /\ dphix /\ dciotax onto the image volume dAy /\ dpsiy /\ dphiy /\ dciotay. This is
    // given by the Jacobian.
    
    // compute the Jacobian
    for(i=0; i< 4; i++)  // step through the parameters at x
    {
        
        for(j=0; j<= 5; j++) paramsxp[j] = paramsx[j];
        for(j=0; j<= 1; j++) paramsyp[j] = paramsy[j];
        paramsxp[i+2] += ep;
        
        skymap(paramsxp, paramsyp, GMST, ifo1, ifo2, iref);
        
        // can get thrown to the alternative solution of phi, psi. Map these back
        if( fabs(paramsyp[2]+PI/2.0-paramsy[2]) < fabs(paramsyp[2]-paramsy[2]) ) paramsyp[2] += PI/2.0;
        if( fabs(paramsyp[2]-PI/2.0-paramsy[2]) < fabs(paramsyp[2]-paramsy[2]) ) paramsyp[2] -= PI/2.0;
        if( fabs(paramsyp[5]+PI-paramsy[5]) < fabs(paramsyp[5]-paramsy[5]) ) paramsyp[5] += PI;
        if( fabs(paramsyp[5]-PI-paramsy[5]) < fabs(paramsyp[5]-paramsy[5]) ) paramsyp[5] -= PI;
        
        for(j=0; j< 4; j++) Jmat[i*4+j] = (paramsyp[j+2]-paramsy[j+2])/ep;
        
    }
    
    
    Jack = fabs(det(Jmat, 4));
    
    // printf("%f %f %f\n", Jmat[2*4+2], paramsx[4], paramsy[4]);
    
     /*      printf("\n");
     for(i=0; i< 4; i++)
     {
     for(j=0; j< 4; j++)
     {
     printf("%f ", Jmat[i*4+j]);
     }
     printf("\n");
     }
     printf("\n");
    
     */
    
    free_double_vector(paramsyp);
    free_double_vector(paramsxp);
    free_double_vector(Jmat);
    
    return Jack;
    
    
}

void uvwz_all(double *u, double *v, double *w, double *z, double *params)
{
    double psi, ecc, Amp, phi, Ap, Ac;
    double cphi, sphi, c2p, s2p;
    double ciota;
    
    // [0] log(Mc) [1] log(Mt) [2] chi1 [3] chi2 [4] phi0 [5] tp [6] log(DL) [7] alpha [8] sindelta [9] psi [10] ciota
    
    psi = params[9];
    ciota = params[10];
    Ap = 0.5*(1.0+ciota*ciota);
    Ac = -ciota;
    ecc = Ac/Ap;
    Amp = Ap/exp(params[6]);
    phi  = params[4];
    
    cphi = cos(phi);
    sphi = sin(phi);
    c2p = cos(2.0*psi);
    s2p = sin(2.0*psi);
    
    *u = Amp*(cphi*c2p+ecc*sphi*s2p);
    *v = Amp*(cphi*s2p-ecc*sphi*c2p);
    *w = Amp*(sphi*c2p-ecc*cphi*s2p);
    *z = Amp*(sphi*s2p+ecc*cphi*c2p);
    
    return;
    
}

void exsolve_all(double *phiy, double *psiy, double *DLy, double *ciotay, double uy, double vy, double wy, double zy)
{
    
    double q, rad;
    double x, p, ecc;
    int flag1, flag2, flag3;
    double Ay, ci;
    
    q = 2.0*(uy*wy+vy*zy)/((uy*uy+vy*vy) - (wy*wy+zy*zy));
    
    rad = sqrt(1.0+q*q);
    
    x = atan2(2.0*(uy*wy+vy*zy),  ((uy*uy+vy*vy) - (wy*wy+zy*zy)));
    if(x < 0.0) x += TPI;
    
    *phiy = 0.5*x;
    
    flag1 = 0;
    if(cos(x) > 0.0) flag1 = 1;
    
    flag2 = 0;
    if(cos(x/2.0) > 0.0) flag2 = 1;
    
    flag3 = 0;
    if(sin(x/2.0) > 0.0) flag3 = 1;
    
    p = ((rad+1.0)*(uy*uy+vy*vy)+(rad-1.0)*(wy*wy+zy*zy)+2.0*q*(uy*wy+vy*zy))/(2.0*rad);
    
    
    
    if (flag1 == 1)
    {
        ecc = (uy*zy-vy*wy)/p;
    }
    else
    {
        ecc = p/(uy*zy-vy*wy);
    }
    
    ci = -(1.0-sqrt(1.0-ecc*ecc))/ecc;
    
    *ciotay = ci;
    
    
    if(flag1 == 1 && flag2 == 1)
    {
        x = atan2((1.0+rad)*vy+q*zy,(1.0+rad)*uy+q*wy);
    }
    else if (flag1 == 1 && flag2 == 0)
    {
        x = atan2(-((1.0+rad)*vy+q*zy),-((1.0+rad)*uy+q*wy));
    }
    else if (flag1 == 0 && flag3==1)
    {
        x = atan2( ecc*((1.0+rad)*uy+q*wy), -ecc*((1.0+rad)*vy+q*zy));
    }
    else if (flag1 == 0 && flag3 == 0)
    {
        x = atan2(-ecc*((1.0+rad)*uy+q*wy), ecc*((1.0+rad)*vy+q*zy));
    }
    
    // printf("flag1 %d flag2 %d flag3 %d x %f\n", flag1, flag2, flag3, x);
    
    if(x < 0.0) x += TPI;
    
    
    
    *psiy = 0.5*x;
    
    if(flag1 == 1)
    {
        Ay = sqrt(p);
    }
    else
    {
        Ay = sqrt(p)/fabs(ecc);
    }
    
    *DLy = (1.0+ci*ci)/(2.0*Ay);
    
    return;
    
}



double skydensity_all(double *paramsx, double *paramsy, double GMST, int ifo1, int ifo2)
{
    int i, ii, j;
    double *paramsxp, *paramsyp;
    double *Jmat;
    double x, Jack;
    double ep;
    
    ep = 1.0e-6;
    
    // [0] log(Mc) [1] log(Mt) [2] chi1 [3] chi2 [4] phi0 [5] tp [6] log(DL) [7] alpha [8] sindelta [9] psi [10] ciota
    
    paramsxp = double_vector(NP);
    paramsyp = double_vector(NP);
    
    // 4 x 4 Jacobian
    Jmat = double_vector(16);
    
    
    // With the full signal (as opposed to the skymcmc) we are also adjusting the geocenter arrival time. Do we need this in the Jacobian?
    
    // At this stage all we have for the new location y is the sky location. So first we find what the rest
    // of the parameters (amplitude, polarization, ellipticity, phase) need to be to keep the signal the same
    
    skymap_all(paramsx, paramsy, GMST, ifo1, ifo2);
    
    // The sky locations are chosen to be on the equal time delay sky ring. Presumably the density factors for that
    // mapping have been taken care of. What we need to know here is how the volume of a small volume
    // dlnDLx /\ dpsix /\ dphix /\ dciotax onto the image volume dlnDLy /\ dpsiy /\ dphiy /\ dciotay. This is
    // given by the Jacobian.
    
    // compute the Jacobian
    for(i=0; i<= 3; i++)  // step through the parameters at x
    {
        
        for(j=0; j< NP; j++) paramsxp[j] = paramsx[j];
        for(j=0; j< NP; j++) paramsyp[j] = paramsy[j];
        
        if(i==0)  paramsxp[6] += ep;
        if(i==1)  paramsxp[9] += ep;
        if(i==2)  paramsxp[4] += ep;
        if(i==3)  paramsxp[10] += ep;
        
        skymap_all(paramsxp, paramsyp, GMST, ifo1, ifo2);
        
        
        // sky parameter order
        //[0] alpha, [1] sin(delta) [2] psi [3] ellipticity [4] scale [5] phi0 [6] dt
        // [0] log(Mc) [1] log(Mt) [2] chi1 [3] chi2 [4] phi0 [5] tp [6] log(DL) [7] alpha [8] sindelta [9] psi [10] ciota
        
        // can get thrown to the alternative solution of phi, psi. Map these back
        //psi
        if( fabs(paramsyp[9]+PI/2.0-paramsy[9]) < fabs(paramsyp[9]-paramsy[9]) ) paramsyp[9] += PI/2.0;
        if( fabs(paramsyp[9]-PI/2.0-paramsy[9]) < fabs(paramsyp[9]-paramsy[9]) ) paramsyp[9] -= PI/2.0;
        //gw phase
        if( fabs(paramsyp[4]+PI-paramsy[4]) < fabs(paramsyp[4]-paramsy[4]) ) paramsyp[4] += PI;
        if( fabs(paramsyp[4]-PI-paramsy[4]) < fabs(paramsyp[4]-paramsy[4]) ) paramsyp[4] -= PI;
        
        for(j=0; j<= 3; j++)
        {
            if(j==0) Jmat[i*4+j] = (paramsyp[6]-paramsy[6])/ep;
            if(j==1) Jmat[i*4+j] = (paramsyp[9]-paramsy[9])/ep;
            if(j==2) Jmat[i*4+j] = (paramsyp[4]-paramsy[4])/ep;
            if(j==3) Jmat[i*4+j] = (paramsyp[10]-paramsy[10])/ep;
        }
        
    }
    
    
    Jack = fabs(det(Jmat, 4));
    
    free_double_vector(paramsyp);
    free_double_vector(paramsxp);
    free_double_vector(Jmat);
    
    
    return Jack;
    
    
}




void uvwz(double *u, double *v, double *w, double *z, double *params)
{
    double psi, ciota, ecc, Amp, phi;
    double cphi, sphi, c2p, s2p;
    double Ap, Ac;
    
    psi = params[2];
    ciota = params[3];
    Amp = params[4];
    phi  = params[5];
    
    Ap = 0.5*(1.0+ciota*ciota);
    Ac = -ciota;
    ecc = Ac/Ap;
    
    Amp *= Ap;
    
    cphi = cos(phi);
    sphi = sin(phi);
    c2p = cos(2.0*psi);
    s2p = sin(2.0*psi);
    
    *u = Amp*(cphi*c2p+ecc*sphi*s2p);
    *v = Amp*(cphi*s2p-ecc*sphi*c2p);
    *w = Amp*(sphi*c2p-ecc*cphi*s2p);
    *z = Amp*(sphi*s2p+ecc*cphi*c2p);
    
    return;
    
}

void uvwz_sol(double *uy, double *vy, double *wy, double *zy, double ux, double vx, double wx, double zx, \
              double fp1x, double fp1y, double fc1x, double fc1y, double fp2x, double fp2y, double fc2x, double fc2y)
{
    
    double den;
    
    den = fc1y*fp2y-fc2y*fp1y;
    
    *uy = (vx*(fc1y*fc2x-fc1x*fc2y)+ux*(fc1y*fp2x - fc2y*fp1x))/den;
    
    *vy = (vx*(fc1x*fp2y-fc2x*fp1y)+ux*(fp2y*fp1x-fp1y*fp2x))/den;
    
    *wy = (zx*(fc1y*fc2x-fc1x*fc2y)+wx*(fc1y*fp2x - fc2y*fp1x))/den;
    
    *zy = (zx*(fc1x*fp2y-fc2x*fp1y)+wx*(fp2y*fp1x-fp1y*fp2x))/den;
    
    return;
    
}

void exsolve(double *phiy, double *psiy, double *Ay, double *ciotay, double uy, double vy, double wy, double zy)
{
    
    double q, rad;
    double x, p, ecc, ci;
    int flag1, flag2, flag3;
    
    q = 2.0*(uy*wy+vy*zy)/((uy*uy+vy*vy) - (wy*wy+zy*zy));
    
    rad = sqrt(1.0+q*q);
    
    x = atan2(2.0*(uy*wy+vy*zy),  ((uy*uy+vy*vy) - (wy*wy+zy*zy)));
    if(x < 0.0) x += TPI;
    
    *phiy = 0.5*x;
    
    flag1 = 0;
    if(cos(x) > 0.0) flag1 = 1;
    
    flag2 = 0;
    if(cos(x/2.0) > 0.0) flag2 = 1;
    
    flag3 = 0;
    if(sin(x/2.0) > 0.0) flag3 = 1;
    
    p = ((rad+1.0)*(uy*uy+vy*vy)+(rad-1.0)*(wy*wy+zy*zy)+2.0*q*(uy*wy+vy*zy))/(2.0*rad);
    
    
    
    if (flag1 == 1)
    {
        ecc = (uy*zy-vy*wy)/p;
    }
    else
    {
        ecc = p/(uy*zy-vy*wy);
    }
    
    ci = -(1.0-sqrt(1.0-ecc*ecc))/ecc;
    
    *ciotay = ci;
    
    
    if(flag1 == 1 && flag2 == 1)
    {
        x = atan2((1.0+rad)*vy+q*zy,(1.0+rad)*uy+q*wy);
    }
    else if (flag1 == 1 && flag2 == 0)
    {
        x = atan2(-((1.0+rad)*vy+q*zy),-((1.0+rad)*uy+q*wy));
    }
    else if (flag1 == 0 && flag3==1)
    {
        x = atan2( ecc*((1.0+rad)*uy+q*wy), -ecc*((1.0+rad)*vy+q*zy));
    }
    else if (flag1 == 0 && flag3 == 0)
    {
        x = atan2(-ecc*((1.0+rad)*uy+q*wy), ecc*((1.0+rad)*vy+q*zy));
    }
    
    // printf("flag1 %d flag2 %d flag3 %d x %f\n", flag1, flag2, flag3, x);
    
    if(x < 0.0) x += TPI;
    
    
    
    *psiy = 0.5*x;
    
    if(flag1 == 1)
    {
        x = sqrt(p);
    }
    else
    {
        x = sqrt(p)/fabs(ecc);
    }
    
    *Ay = 2.0*x/(1.0+ci*ci);
    
    return;
    
}

void pmap(struct Net *net, double *pall, double *param, double *sky)
{
    int j;
    double ciota, Ap, Ac, alpha, sindelta, psi, Fs;
    double *dtimes;
    double Fplus, Fcross, lambda;
    
    dtimes = (double*)malloc(sizeof(double)* 5);
    
        for (j = 0; j < NX; ++j)
        {
            pall[j] = param[j];
        }
        
        ciota = sky[3];
        Ap = (1.0+ciota*ciota)/2.0;
        Ac = -ciota;
        alpha = sky[0];
        sindelta = sky[1];
        psi = sky[2];
        
        F_ant(psi, alpha, sindelta, net->GMST, &Fplus, &Fcross, net->labels[0]);
        Fs = sqrt(Ap*Ap*Fplus*Fplus+Ac*Ac*Fcross*Fcross);
        lambda = atan2(Ac*Fcross,Ap*Fplus);
        if(lambda < 0.0) lambda += TPI;
    
        Times(alpha, sindelta, net->GMST, dtimes);
    
        // move reference point from ref detector to geocenter
        pall[4] = 2.0*param[4] - lambda;  // I'm holding the GW phase in pall[4]
        pall[5] = param[5] - dtimes[net->labels[0]];
        pall[6] = param[6] + log(Fs);
    
        pall[7] = alpha;
        pall[8] = sindelta;
        pall[9] = psi;
        pall[10] = ciota;
    
    free(dtimes);
}

void pmap_back(struct Net *net, double *pall, double *param, double *sky)
{
    int j;
    double ciota, Ap, Ac, alpha, sindelta, psi, Fs;
    double *dtimes;
    double Fplus, Fcross, lambda;
    
    dtimes = (double*)malloc(sizeof(double)* 5);
    
    // sky  [0] alpha, [1] sin(delta) [2] psi [3] ciota [4] scale [5] dphi [6] dt
    // param [0] log(Mc) [1] log(Mt) [2] chi1 [3] chi2 [4] phi0  [5] tp0 [6] log(DL0) then relative amplitudes, time, phases
    // pall  [0] log(Mc) [1] log(Mt) [2] chi1 [3] chi2 [4] 2phi0 [5] tp [6] log(DL) [7] alpha [8] sindelta [9] psi [10] ciota
    
    
    for (j = 0; j < NX; ++j)
    {
        param[j] = pall[j];
    }
    
    
    ciota = pall[10];
    Ap = (1.0+ciota*ciota)/2.0;
    Ac = -ciota;
    alpha = pall[7];
    sindelta = pall[8];
    psi = pall[9];
    
    sky[0] = alpha;
    sky[1] = sindelta;
    sky[2] = psi;
    sky[3] = ciota;
    sky[4] = 1.0;
    sky[5] = 0.0;
    sky[6] = 0.0;
    
    F_ant(psi, alpha, sindelta, net->GMST, &Fplus, &Fcross, net->labels[0]);
    Fs = sqrt(Ap*Ap*Fplus*Fplus+Ac*Ac*Fcross*Fcross);
    lambda = atan2(Ac*Fcross,Ap*Fplus);
    if(lambda < 0.0) lambda += TPI;
    
    Times(alpha, sindelta, net->GMST, dtimes);
    
    param[4] = 0.5*(pall[4]+lambda);
    param[5] = pall[5] + dtimes[net->labels[0]];
    param[6] = pall[6] - log(Fs);
    
    free(dtimes);
}

void pairs(struct Net *net)
{
    
  if(net->Nifo > 1)
   {
    //tds hold the maximum time delays between the reference detector and the other detectors
    if(net->labels[0] == 0)
    {
        if(net->labels[1] == 1)
        {
            net->tds[1] = HLdt;
          if(net->Nifo > 2)  net->tds[2] = HVdt;
        }
        if(net->labels[1] == 2)
        {
            net->tds[1] = HVdt;
           if(net->Nifo > 2) net->tds[2] = HLdt;
        }
        
    }
    if(net->labels[0] == 1)
    {
        if(net->labels[1] == 0)
        {
            net->tds[1] = HLdt;
            if(net->Nifo > 2) net->tds[2] = LVdt;
        }
        if(net->labels[1] == 2)
        {
            net->tds[1] = LVdt;
            if(net->Nifo > 2) net->tds[2] = HLdt;
        }
    }
    if(net->labels[0] == 2)
    {
        if(net->labels[1] == 1)
        {
            net->tds[1] = LVdt;
           if(net->Nifo > 2)  net->tds[2] = HVdt;
        }
        if(net->labels[1] == 0)
        {
            net->tds[1] = HVdt;
           if(net->Nifo > 2) net->tds[2] = LVdt;
        }
    }
        
    }
    
    return;
}

void time_delays(struct Net *net)
{
    
   int i, j;
    
   for(i = 0; i < net->Nifo; i++) net->delays[i][i] = 0.0;
   for(i = 0; i < net->Nifo; i++)
   {
    for(j = i+1; j < net->Nifo; j++)
     {
     if(net->labels[i]==0 && net->labels[j]==1) net->delays[i][j] = HLdt;
     if(net->labels[i]==1 && net->labels[j]==0) net->delays[i][j] = HLdt;
         
     if(net->labels[i]==0 && net->labels[j]==2) net->delays[i][j] = HVdt;
     if(net->labels[i]==2 && net->labels[j]==0) net->delays[i][j] = HVdt;
         
     if(net->labels[i]==1 && net->labels[j]==2) net->delays[i][j] = LVdt;
     if(net->labels[i]==2 && net->labels[j]==1) net->delays[i][j] = LVdt;

   }
  }
    
  for(i = 0; i < net->Nifo ; i++)
   {
   for(j = i+1; j < net->Nifo; j++)
    {
     net->delays[j][i] = net->delays[i][j];
   }
  }
    
}

void detector_shifts(struct Net *net, double *params)
{
    

    int i, j, id;
    double alpha, sindelta, psi, ciota;
    double dt, dp, dA;
    double Ap, Ac;
    double Fplus, Fcross;
    double *dtimes;
    double *lambda, *Fs;
    
    dtimes = (double*)malloc(sizeof(double)*5);
    lambda = (double*)malloc(sizeof(double)*net->Nifo);
    Fs = (double*)malloc(sizeof(double)*net->Nifo);

    alpha = params[7];
    sindelta = params[8];
    psi = params[9];
    ciota = params[10];
    
    Ap = (1.0+ciota*ciota)/2.0;
    Ac = -ciota;
    
    Times(alpha, sindelta, net->GMST, dtimes);
    
    for (id=0; id< net->Nifo; id++)
    {
        
        F_ant(psi, alpha, sindelta, net->GMST, &Fplus, &Fcross, net->labels[id]);
        Fs[id] = sqrt(Ap*Ap*Fplus*Fplus+Ac*Ac*Fcross*Fcross);  // magnitude of response
        lambda[id] = atan2(Ac*Fcross,Ap*Fplus);
        if(lambda[i] < 0.0) lambda[id] += TPI;
    }
    
    for (i=0; i< net->Nifo; i++)
    {
        for (j=i+1; j< net->Nifo; j++)
        {
            dt = (dtimes[net->labels[i]]-dtimes[net->labels[j]]);
            dp = lambda[i]-lambda[j];
            if(dp < 0.0) dp += TPI;
            if(dp > TPI) dp -= TPI;
            dA = Fs[i]/Fs[j];
            printf("detector %d-%d  delay %e  phase difference %f amp ratio %f\n", net->labels[i], net->labels[j], dt, dp, dA);
        }
    }
    
    free(lambda);
    free(Fs);
    free(dtimes);

}

void dshifts(struct Net *net, double *sky, double *params)
{
    int i, j, k, id;
    double dtimes[4];
    double Fplus[4], Fcross[4], F[4];
    double lambda[4];
    double logL, logLmax, ecc;
    
    // sky parameter order
    //[0] alpha, [1] sin(delta) [2] psi [3] cos(iota) [4] scale [5] dt [6] phi0
    
             Times(sky[0], sky[1], net->GMST, dtimes);
    
            for(id = 0; id<net->Nifo; id++)
            {
                ecc = -2.0*sky[3]/(1.0+sky[3]*sky[3]);
                F_ant(sky[2], sky[0], sky[1], net->GMST, &Fplus[id], &Fcross[id], net->labels[id]);
                F[id] = sqrt(Fplus[id]*Fplus[id]+ecc*ecc*Fcross[id]*Fcross[id]);
                lambda[id] = atan2(ecc*Fcross[id],Fplus[id]);
                if(lambda[id] < 0.0) lambda[id] += TPI;
            }
    
    // Note that F[] and lambda[] have already used the label array to get the correct detector ordering while dtimes uses a fixed detector ordering so needs the label array to get the correct delays
    
    for(id = 1; id<net->Nifo; id++)
    {
    params[(id-1)*3+NX] = lambda[id]-lambda[0]; // GW phase offset
    if(params[(id-1)*3+NX] < 0.0) params[(id-1)*3+NX] += TPI;
    if(params[(id-1)*3+NX] > TPI) params[(id-1)*3+NX] -= TPI;
    params[(id-1)*3+NX+1] = dtimes[net->labels[id]]-dtimes[net->labels[0]]; // time offset
    params[(id-1)*3+NX+2] = F[id]/F[0]; // amplitude ratio
    }

    
    
}

void ringfind(struct Net *net, double *tdelays, double *params, double *SNRsq, gsl_rng * r)
{
    double alpha, sindelta, beta;
    double alphax, sindeltax, lx;
    double *dtimes;
    double SNRT;
    double tolt, x, y, logp, dtd;
    int i, id, j, k;
    
    dtimes = (double*)malloc(sizeof(double)* 5);
    
    tolt = 2.0e-4;
    
    lx = -1.0e30;
    
    SNRT = 0.0;
    for(id = 1; id<net->Nifo; id++) SNRT += SNRsq[id];

    for(i = 0; i< 40000; i++)  // 40K trials is roughly a 1 degree search
    {

        alpha = TPI*gsl_rng_uniform(r);
        sindelta = -1.0+2.0*gsl_rng_uniform(r);
        Times(alpha, sindelta, net->GMST, dtimes);
        
        logp = 0.0;
        for(id = 1; id<net->Nifo; id++)
        {
            dtd = dtimes[net->labels[0]]-dtimes[net->labels[id]];
            x = (tdelays[id]-dtd)/tolt;
            y = x*x*SNRsq[id]/SNRT;
            logp -= y/2.0;
        }
 
        if(logp > lx)
        {
            lx = logp;
            alphax = alpha;
            sindeltax = sindelta;
        }
  
    }
    
    params[0] = alphax;
    params[1] = sindeltax;
    
    free(dtimes);
    
    return;
    
}




