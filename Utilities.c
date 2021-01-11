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

#include "Utilities.h"
#include "ConstCBC.h"
#include "Constants.h"

void qscan(double *data, double *Sn, double Tobs, int N)
{
    double t_rise, alpha;
    double dt, fmx, fmn;
    double *dcopy;
    double *freqs;
    double x, dx, dlnf;
    double t, f;
    int subscale, octaves, Nf, i, j;
    double **tfDR, **tfDI, **tfD;
    double fac;
    double *SX;
    
    FILE *out;
    
    dt = Tobs/(double)(N);
    fmx = 1.0/(2.0*dt);
    fmn = 8.0;
    
    dcopy = (double*)malloc(sizeof(double)* (N));
    SX = (double*)malloc(sizeof(double)* (N/2));
    
    // Qscan uses different scaling convention
    fac = Tobs/((double)(N)*(double)(N));
    for (i = 0; i < N/2; ++i) SX[i] = Sn[i]/fac;

    // Tukey window parameter. Flat for (1-alpha) of data
    t_rise = 0.4; // Standard LAL setting
    alpha = (2.0*t_rise/Tobs);
    
    // make copy
       for (i = 0; i < N; ++i) dcopy[i] = data[i];

    // Tukey window
    tukey(dcopy, alpha, N);
    // FFT
    gsl_fft_real_radix2_transform(dcopy, 1, N);
    // whiten data
    whiten(dcopy, SX, N);
    
    // logarithmic frequency spacing
    subscale = 40;  // number of semi-tones per octave
    octaves = (int)(rint(log(fmx/fmn)/log(2.0))); // number of octaves
    Nf = subscale*octaves+1;
    freqs = (double*)malloc(sizeof(double)* (Nf));   // frequencies used in the analysis
    dx = log(2.0)/(double)(subscale);
    dlnf = dx;
    x = log(fmn);
    for(i=0; i< Nf; i++)
    {
        freqs[i] = exp(x);
        x += dx;
    }
    
   // printf("%d freqs %d times\n", Nf, N);
    
    tfDR = double_matrix(Nf,N);
    tfDI = double_matrix(Nf,N);
    tfD = double_matrix(Nf,N);
    
    // Wavelet transform
    TransformC(dcopy, freqs, tfD, tfDR, tfDI, Qprint, Tobs, N, Nf);
    
    
    out = fopen("Qtran.dat","w");
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
    
    free(SX);
    free(dcopy);
    free(freqs);
    free_double_matrix(tfDR,Nf);
    free_double_matrix(tfDI,Nf);
    free_double_matrix(tfD,Nf);
    
}

//uses frequency domain data
void qscanf(double *data, double *Sn, double Tobs, int N)
{
    double t_rise, alpha;
    double dt, fmx, fmn;
    double *dcopy;
    double *freqs;
    double x, dx, dlnf;
    double t, f;
    int subscale, octaves, Nf, i, j;
    double **tfDR, **tfDI, **tfD;
    double fac;
    double *SX;
    
    FILE *out;
    
    dt = Tobs/(double)(N);
    fmx = 1.0/(2.0*dt);
    fmn = 8.0;
    
    dcopy = (double*)malloc(sizeof(double)* (N));
    SX = (double*)malloc(sizeof(double)* (N/2));
    
    // Qscan uses different scaling convention
    fac = Tobs/((double)(N)*(double)(N));
    for (i = 0; i < N/2; ++i) SX[i] = Sn[i]/fac;

    // Tukey window parameter. Flat for (1-alpha) of data
    t_rise = 0.4; // Standard LAL setting
    alpha = (2.0*t_rise/Tobs);
    
    // make copy
    for (i = 0; i < N; ++i) dcopy[i] = data[i];

    whiten(dcopy, SX, N);
    
    // logarithmic frequency spacing
    subscale = 40;  // number of semi-tones per octave
    octaves = (int)(rint(log(fmx/fmn)/log(2.0))); // number of octaves
    Nf = subscale*octaves+1;
    freqs = (double*)malloc(sizeof(double)* (Nf));   // frequencies used in the analysis
    dx = log(2.0)/(double)(subscale);
    dlnf = dx;
    x = log(fmn);
    for(i=0; i< Nf; i++)
    {
        freqs[i] = exp(x);
        x += dx;
    }
    
   // printf("%d freqs %d times\n", Nf, N);
    
    tfDR = double_matrix(Nf,N);
    tfDI = double_matrix(Nf,N);
    tfD = double_matrix(Nf,N);
    
    // Wavelet transform
    TransformC(dcopy, freqs, tfD, tfDR, tfDI, Qprint, Tobs, N, Nf);
    
    
    out = fopen("Qsig.dat","w");
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
    
    free(SX);
    free(dcopy);
    free(freqs);
    free_double_matrix(tfDR,Nf);
    free_double_matrix(tfDI,Nf);
    free_double_matrix(tfD,Nf);
    
}


void qscanres(double *data, double *signal, double *Sn, double Tobs, int N)
{
    double t_rise, alpha;
    double dt, fmx, fmn, fac;
    double *dcopy;
    double *freqs;
    double x, dx, dlnf;
    double t, f;
    int subscale, octaves, Nf, i, j;
    double **tfDR, **tfDI, **tfD;
    double *SX;
    
    FILE *out;
    
    dt = Tobs/(double)(N);
    fmx = 1.0/(2.0*dt);
    fmn = 8.0;
    
    dcopy = (double*)malloc(sizeof(double)* (N));
    SX = (double*)malloc(sizeof(double)* (N/2));

    // Tukey window parameter. Flat for (1-alpha) of data
    t_rise = 0.4; // Standard LAL setting
    alpha = (2.0*t_rise/Tobs);
    
    // make copy
    for (i = 0; i < N; ++i) dcopy[i] = data[i];
    
    // Qscan uses different scaling convention
    fac = Tobs/((double)(N)*(double)(N));
    for (i = 0; i < N/2; ++i) SX[i] = Sn[i]/fac;

    // Tukey window
    tukey(dcopy, alpha, N);
    // FFT
    gsl_fft_real_radix2_transform(dcopy, 1, N);
    // subtract the signal
    for (i = 0; i < N; ++i) dcopy[i] -= signal[i];
    // whiten data
    whiten(dcopy, SX, N);
    
    // logarithmic frequency spacing
    subscale = 40;  // number of semi-tones per octave
    octaves = (int)(rint(log(fmx/fmn)/log(2.0))); // number of octaves
    Nf = subscale*octaves+1;
    freqs = (double*)malloc(sizeof(double)* (Nf));   // frequencies used in the analysis
    dx = log(2.0)/(double)(subscale);
    dlnf = dx;
    x = log(fmn);
    for(i=0; i< Nf; i++)
    {
        freqs[i] = exp(x);
        x += dx;
    }
    
   // printf("%d freqs %d times\n", Nf, N);
    
    tfDR = double_matrix(Nf,N);
    tfDI = double_matrix(Nf,N);
    tfD = double_matrix(Nf,N);
    
    // Wavelet transform
    TransformC(dcopy, freqs, tfD, tfDR, tfDI, Qprint, Tobs, N, Nf);
    
    
    out = fopen("Qtranres.dat","w");
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
    
    free(SX);
    free(dcopy);
    free(freqs);
    free_double_matrix(tfDR,Nf);
    free_double_matrix(tfDI,Nf);
    free_double_matrix(tfD,Nf);
    
}


void whiten(double *data, double *Sn, int N)
{
    double f, x, y, fix;
    int i;
    
    data[0] = 0.0;
    data[N/2] = 0.0;
    
    for(i=1; i< N/2; i++)
    {
        x = 1.0/sqrt(Sn[i]);
        data[i] *= x;
        data[N-i] *= x;
    }
    
}


// QNM frequency
double fQNM(double m1, double m2, double chi1, double chi2)
{
 
    double q, eta, f;
    double Z1, Z2;
    double risco, Eisco, Lisco;
    double Mf, atot, aeff, af;
    double p0, p1, k00, k01, k02, k10, k11, k12, scri;
    
    p0 = 0.04827;   // fit parameters
    p1 = 0.01707;
    k00 = -3.82;
    k01 = -1.2019;
    k02 = -1.20764;
    k10 = 3.79245;
    k11 = 1.18385;
    k12 = 4.90494;
    scri = 0.41616;

    if(m1 < m2)
    {
        q = m1/m2;
    }
    else
    {
        q = m2/m1;
    }
    
    eta = q/((1.0+q)*(1.0+q));    // symmetric mass ratio
    atot = (chi1+q*q*chi2)/((1.0+q)*(1.0+q));
    aeff = atot+scri*eta*(chi1+chi2);   //technically, only 1605.01938 uses this, but the result for Mf is the same whether you use atot or aeff
    Z1 = 1.0 + pow((1.0 - aeff*aeff),(1./3.))*(pow((1. + aeff),(1./3.)) + pow((1. - aeff),(1./3.)));
    Z2 = sqrt(3.0*aeff*aeff + Z1*Z1);
    risco = 3.0 + Z2 - (aeff/(fabs(aeff)+1.0e-20))*sqrt((3.0 - Z1)*(3.0 + Z1 + 2.0*Z2));
    Eisco = sqrt(1.0-2./(3.0*risco));
    Lisco = (2./(3.0*sqrt(3.)))*(1.+2.0*sqrt(3.0*risco-2.0));
    
   // printf("%f %f %f %f %f %f\n", eta, atot, Z1, Z2, risco, Eisco);
    
    Mf = (1.0-eta*(1.0-Eisco)-4.0*eta*eta*(4.0*p0+16.0*p1*atot*(atot+1.0)+Eisco-1.0))*(m1+m2);
    
    af = atot + eta*(Lisco-2*atot*(Eisco-1.0)) + k00*eta*eta + k01*eta*eta*aeff \
    + k02*eta*eta*aeff*aeff + k10*eta*eta*eta + k11*eta*eta*eta*aeff + k12*eta*eta*eta*aeff*aeff;
    
    f = (1.5251-1.1568*pow((1.0-af),0.1292))/(TPI*Mf);  // Berti fit, from gr-qc/0512160
    
    return(f);
    
}

double t_at_f(double *param, double f)
{
    int i;
    double t, x;
    double dt, fring, a, Mchirp, M, eta, dm, m1, m2, chi, chi1, chi2, tc, theta;
    double gamma_E=0.5772156649; //Euler's Constant-- shows up in 3PN term
    
    Mchirp = exp(param[0])/MSUN_SI*TSUN;
    M = exp(param[1])/MSUN_SI*TSUN;
    eta = pow((Mchirp/M), (5.0/3.0));
    dm = sqrt(1.0-4.0*eta);
    m1 = M*(1.0+dm)/2.0;
    m2 = M*(1.0-dm)/2.0;
    chi1 = param[2];
    chi2 = param[3];

    chi = (m1*chi1+m2*chi2)/M;
    
    tc = param[5];
    
    x = pow(PI*M*f,2.0/3.0);
    
    t =  tc-pow(x,-4.0)*5.0*M/(256.0*eta)*(1.0      // 0 PN
                        +((1/252.0)*(743.0+924.0*eta)*x)  // 1 PN
                        +(2.0/15.0)*(-48.0*PI+113.0*chi-38.0*eta*(chi1+chi2))*pow(x,3.0/2.0)    // 1.5 PN
                        +(1.0/508032.0)*(3058673.0+5472432.0*eta+4353552.0*eta*eta-5080320.0*chi*chi+127008.0*eta*chi1*chi2)*pow(x,4.0/2.0) // 2 PN
                        +(1.0/756.0)*(3.0*PI*(-7729.0+1092.0*eta)+1512.0*chi*chi*chi-56.0*eta*(1285.0+153.0*eta)*(chi1+chi2)+chi*(147101.0+504.0*eta*(13.0-9.0*chi1*chi2)))*pow(x,5.0/2.0)// 2.5 PN
                        +((6848.0*gamma_E)/105.0+PI*PI/12.0*(512.0-451.0*eta)+(25565.0*eta*eta*eta)/1296.0+(-10052469856691.0+24236159077900.0*eta)/23471078400.0+(35.0*chi*chi)/3.0+eta*eta*(-(15211.0/1728.0)+19.0*chi1*chi1+(245.0*chi1*chi2)/6.0+19.0*chi2*chi2)+8.0*PI/3.0*(-73.0*chi+28.0*eta*(chi1+chi2))-eta/336.0*(26992.0*chi*chi-995.0*chi1*chi2+30688.0*chi*(chi1+chi2))+3424.0/105.0*log(16.0*x))*pow(x,6.0/2.0));  // 3PN
    
    return(t);
    
}


void t_of_f(double *param, double *times, double *freqs, int N)
{
    int i;
    double t, x;
    double f, dt, fring, a, Mchirp, M, eta, dm, m1, m2, chi, chi1, chi2, tc, theta;
    double gamma_E=0.5772156649; //Euler's Constant-- shows up in 3PN term
    
    Mchirp = exp(param[0])/MSUN_SI*TSUN;
    M = exp(param[1])/MSUN_SI*TSUN;
    eta = pow((Mchirp/M), (5.0/3.0));
    dm = sqrt(1.0-4.0*eta);
    m1 = M*(1.0+dm)/2.0;
    m2 = M*(1.0-dm)/2.0;
    chi1 = param[2];
    chi2 = param[3];

    chi = (m1*chi1+m2*chi2)/M;
    
    tc = param[5];
    
    for (i = 0; i < N; ++i)
    {
        
    f = freqs[i];
    
    x = pow(PI*M*f,2.0/3.0);
    
    times[i] =  tc-pow(x,-4.0)*5.0*M/(256.0*eta)*(1.0      // 0 PN
                                           +((1/252.0)*(743.0+924.0*eta)*x)  // 1 PN
                                           +(2.0/15.0)*(-48.0*PI+113.0*chi-38.0*eta*(chi1+chi2))*pow(x,3.0/2.0)    // 1.5 PN
                                           + (1.0/508032.0)*(3058673.0+5472432.0*eta+4353552.0*eta*eta-5080320.0*chi*chi+127008.0*eta*chi1*chi2)*pow(x,4.0/2.0) // 2 PN
                                           +(1.0/756.0)*(3.0*PI*(-7729.0+1092.0*eta)+1512.0*chi*chi*chi-56.0*eta*(1285.0+153.0*eta)*(chi1+chi2)+chi*(147101.0+504.0*eta*(13.0-9.0*chi1*chi2)))*pow(x,5.0/2.0)// 2.5 PN
                                                  +((6848.0*gamma_E)/105.0+PI*PI/12.0*(512.0-451.0*eta)+(25565.0*eta*eta*eta)/1296.0+(-10052469856691.0+24236159077900.0*eta)/23471078400.0+(35.0*chi*chi)/3.0+eta*eta*(-(15211.0/1728.0)+19.0*chi1*chi1+(245.0*chi1*chi2)/6.0+19.0*chi2*chi2)+8.0*PI/3.0*(-73.0*chi+28.0*eta*(chi1+chi2))-eta/336.0*(26992.0*chi*chi-995.0*chi1*chi2+30688.0*chi*(chi1+chi2))+3424.0/105.0*log(16.0*x))*pow(x,6.0/2.0));  // 3PN
 
    }
    
    
    
}

double f_at_t(double *param, double t)
{
    // [0] log(Mc) [1] log(Mt) [2] chi1 [3] chi2 [4] phi0 [5] tp [6] log(DL) [7] alpha [8] sindelta [9] psi [10] ciota
    
    int i;
    double f, dt, fring, a, Mchirp, M, eta, dm, m1, m2, chi, chi1, chi2, tc, theta;
    double gamma_E=0.5772156649; //Euler's Constant-- shows up in 3PN term
    double PN1, PN15, PN2, PN25, PN3, PN35;
    double theta2, theta3, theta4, theta5, theta6, theta7;
    
    Mchirp = exp(param[0])/MSUN_SI*TSUN;
    M = exp(param[1])/MSUN_SI*TSUN;
    eta = pow((Mchirp/M), (5.0/3.0));
    dm = sqrt(1.0-4.0*eta);
    m1 = M*(1.0+dm)/2.0;
    m2 = M*(1.0-dm)/2.0;
    chi1 = param[2];
    chi2 = param[3];

    chi = (m1*chi1+m2*chi2)/M;
    
    tc = param[5];
    
    fring = fQNM(m1,m2,chi1,chi2);
    
    //printf("fring %f\n", fring);
        
    if(tc > t)
    {
        
        theta = pow(eta*(tc-t)/(5.0*M),-1.0/8.0);
        theta2 = theta*theta;
        theta3 = theta2*theta;
        theta4 = theta2*theta2;
        theta5 = theta2*theta3;
        theta6 = theta3*theta3;
        
        
        PN1 = (11.0/32.0*eta+743.0/2688.0)*theta2;
        PN15 = -3.0*PI/10.0*theta3 + (1.0/160.0)*(113.0*chi-38.0*eta*(chi1+chi2))*theta3;
        PN2 = (1855099.0/14450688.0+56975.0/258048.0*eta+371.0/2048.0*eta*eta)*theta4 + (1.0/14450688.0)*(-3386880.0*chi*chi+1512.0*chi1*chi2)*theta4;
        PN25 = -(7729.0/21504.0-13.0/256.0*eta)*PI*theta5;
        PN3 = (-720817631400877.0/288412611379200.0+53.0/200.0*PI*PI+107.0/280.0*gamma_E+(25302017977.0/4161798144.0-451.0/2048.0*PI*PI)*eta-30913.0/1835008.0*eta*eta+235925.0/1769472.0*eta*eta*eta + 107.0/280.0*log(2.0*theta))*theta6;
        
        //printf("%f %f %e %e %e %e %e\n", t, theta3/(8.0*M*PI), PN1, PN15, PN2, PN25, PN3);
        
        f = theta3/(8.0*M*PI)*(1.0 + PN1 + PN15 + PN2 + PN25 + PN3);
        
        if(f > fring) f = fring;
    }
    else
    {
        f = fring;
    }
    
    return(f);
    
}



void f_of_t(double *param, double *times, double *freqs, int N)
{
    // [0] log(Mc) [1] log(Mt) [2] chi1 [3] chi2 [4] phi0 [5] tp [6] log(DL) [7] alpha [8] sindelta [9] psi [10] ciota
    
    int i;
    double f, t, dt, fring, a, Mchirp, M, eta, dm, m1, m2, chi, chi1, chi2, tc, theta;
    double gamma_E=0.5772156649; //Euler's Constant-- shows up in 3PN term
    double PN1, PN15, PN2, PN25, PN3, PN35;
    double theta2, theta3, theta4, theta5, theta6, theta7;
    
    Mchirp = exp(param[0])/MSUN_SI*TSUN;
    M = exp(param[1])/MSUN_SI*TSUN;
    eta = pow((Mchirp/M), (5.0/3.0));
    dm = sqrt(1.0-4.0*eta);
    m1 = M*(1.0+dm)/2.0;
    m2 = M*(1.0-dm)/2.0;
    chi1 = param[2];
    chi2 = param[3];

    chi = (m1*chi1+m2*chi2)/M;
    
    tc = param[5];
    
    fring = fQNM(m1,m2,chi1,chi2);
    
    //printf("fring %f\n", fring);
    
    for (i = 0; i < N; ++i)
    {
        
    t = times[i];
        
    freqs[i] = -1.0;
        
    if(tc > t)
    {
        
        theta = pow(eta*(tc-t)/(5.0*M),-1.0/8.0);
        theta2 = theta*theta;
        theta3 = theta2*theta;
        theta4 = theta2*theta2;
        theta5 = theta2*theta3;
        theta6 = theta3*theta3;
        
        
        PN1 = (11.0/32.0*eta+743.0/2688.0)*theta2;
        PN15 = -3.0*PI/10.0*theta3 + (1.0/160.0)*(113.0*chi-38.0*eta*(chi1+chi2))*theta3;
        PN2 = (1855099.0/14450688.0+56975.0/258048.0*eta+371.0/2048.0*eta*eta)*theta4 + (1.0/14450688.0)*(-3386880.0*chi*chi+1512.0*chi1*chi2)*theta4;
        PN25 = -(7729.0/21504.0-13.0/256.0*eta)*PI*theta5;
        PN3 = (-720817631400877.0/288412611379200.0+53.0/200.0*PI*PI+107.0/280.0*gamma_E+(25302017977.0/4161798144.0-451.0/2048.0*PI*PI)*eta-30913.0/1835008.0*eta*eta+235925.0/1769472.0*eta*eta*eta + 107.0/280.0*log(2.0*theta))*theta6;
        
        //printf("%f %f %e %e %e %e %e\n", t, theta3/(8.0*M*PI), PN1, PN15, PN2, PN25, PN3);
        
        freqs[i] = theta3/(8.0*M*PI)*(1.0 + PN1 + PN15 + PN2 + PN25 + PN3);
        
        
        if(freqs[i] > fring) freqs[i] = fring;
    }
    else
    {
        freqs[i] = fring;
    }
        
    }
    
    return;
    
}

double fringdown(double *param)
{
    double f, Mchirp, M, eta, dm, m1, m2, chi1, chi2;
    
    Mchirp = exp(param[0])/MSUN_SI*TSUN;
    M = exp(param[1])/MSUN_SI*TSUN;
    eta = pow((Mchirp/M), (5.0/3.0));
    if(eta > 0.25) eta = 0.25;
    dm = sqrt(1.0-4.0*eta);
    m1 = M*(1.0+dm)/2.0;
    m2 = M*(1.0-dm)/2.0;
    chi1 = param[2];
    chi2 = param[3];
    
    f = fQNM(m1, m2, chi1, chi2);
    
    return(f);

}

double fbegin(double *param)
{
    // [0] log(Mc) [1] log(Mt) [2] chi1 [3] chi2 [4] phi0 [5] tp [6] log(DL) [7] alpha [8] sindelta [9] psi [10] ciota
    
    double f, Mchirp, M, eta, dm, m1, m2, chi, chi1, chi2, tc, theta;
    double gamma_E=0.5772156649; //Euler's Constant-- shows up in 3PN term
    double PN1, PN15, PN2, PN25, PN3, PN35;
    double theta2, theta3, theta4, theta5, theta6, theta7;
    
    
    Mchirp = exp(param[0])/MSUN_SI*TSUN;
    M = exp(param[1])/MSUN_SI*TSUN;
    eta = pow((Mchirp/M), (5.0/3.0));
    if(eta > 0.25) eta = 0.25;
    dm = sqrt(1.0-4.0*eta);
    m1 = M*(1.0+dm)/2.0;
    m2 = M*(1.0-dm)/2.0;
    chi1 = param[2];
    chi2 = param[3];
    chi = (m1*chi1+m2*chi2)/M;
    tc = param[5];
    
    
    theta = pow(eta*tc/(5.0*M),-1.0/8.0);
    theta2 = theta*theta;
    theta3 = theta2*theta;
    theta4 = theta2*theta2;
    theta5 = theta2*theta3;
    theta6 = theta3*theta3;
    
    
    PN1 = (11.0/32.0*eta+743.0/2688.0)*theta2;
    PN15 = -3.0*PI/10.0*theta3 + (1.0/160.0)*(113.0*chi-38.0*eta*(chi1+chi2))*theta3;
    PN2 = (1855099.0/14450688.0+56975.0/258048.0*eta+371.0/2048.0*eta*eta)*theta4 + (1.0/14450688.0)*(-3386880.0*chi*chi+1512.0*chi1*chi2)*theta4;
    PN25 = -(7729.0/21504.0-13.0/256.0*eta)*PI*theta5;
    PN3 = (-720817631400877.0/288412611379200.0+53.0/200.0*PI*PI+107.0/280.0*gamma_E+(25302017977.0/4161798144.0-451.0/2048.0*PI*PI)*eta-30913.0/1835008.0*eta*eta+235925.0/1769472.0*eta*eta*eta + 107.0/280.0*log(2.0*theta))*theta6;

    
    f = theta3/(8.0*M*PI)*(1.0 + PN1 + PN15 + PN2 + PN25 + PN3);

    return(f);
    
}





void tukey(double *data, double alpha, int N)
{
  int i, imin, imax;
  double filter;
  
  imin = (int)(alpha*(double)(N-1)/2.0);
  imax = (int)((double)(N-1)*(1.0-alpha/2.0));
  
    int Nwin = N-imax;

 for(i=0; i< N; i++)
  {
    filter = 1.0;
    if(i<imin) filter = 0.5*(1.0+cos(M_PI*( (double)(i)/(double)(imin)-1.0 )));
    if(i>imax) filter = 0.5*(1.0+cos(M_PI*( (double)(i-imax)/(double)(Nwin))));
    data[i] *= filter;
  }
  
}


void tukey_scale(double *s1, double *s2, double alpha, int N)
{
    int i, imin, imax;
    double x1, x2;
    double filter;
    
    imin = (int)(alpha*(double)(N-1)/2.0);
    imax = (int)((double)(N-1)*(1.0-alpha/2.0));
    
    int Nwin = N-imax;
    
    x1 = 0.0;
    x2 = 0.0;
    for(i=0; i< N; i++)
    {
        filter = 1.0;
        if(i<imin) filter = 0.5*(1.0+cos(M_PI*( (double)(i)/(double)(imin)-1.0 )));
        if(i>imax) filter = 0.5*(1.0+cos(M_PI*( (double)(i-imax)/(double)(Nwin))));
        x1 += filter;
        x2 += filter*filter;
    }
    x1 /= (double)(N);
    x2 /= (double)(N);
    
    *s1 = x1;
    *s2 = sqrt(x2);
    
}



void pbt_shift(double *corr, double *corrf, double *data1, double *data2, double *Sn, int imin, int imax, int N)
{
    int nb2, i, l, k, j;
    
    for (i = 0; i < N; i++)
    {
        corr[i]	= 0.0;
        corrf[i] = 0.0;
    }
    
    for (i=imin; i < imax; i++)
    {
        l=i;
        k=N-i;
        
        corr[l]	= (data1[l]*data2[l] + data1[k]*data2[k])/Sn[i];
        corr[k]	= (data1[k]*data2[l] - data1[l]*data2[k])/Sn[i];
        corrf[l] = corr[k];
        corrf[k] = -corr[l];
        
    }
    
    gsl_fft_halfcomplex_radix2_inverse(corr, 1, N);
    gsl_fft_halfcomplex_radix2_inverse(corrf, 1, N);
    
    
}


double fourier_nwip(double *a, double *b, double *Sn, int imin, int imax, int N)
{
    int i, j, k;
    double arg, product;
    double test;
    double ReA, ReB, ImA, ImB;
    
    arg = 0.0;
    for(i=imin; i<imax; i++)
    {
        j = i;
        k = N-i;
        ReA = a[j]; ImA = a[k];
        ReB = b[j]; ImB = b[k];
        product = ReA*ReB + ImA*ImB;
        arg += product/(Sn[i]);
    }
    
    return(4.0*arg);
    
}

double globe(double ***global, double *max, double *min, double Tobs, double *params, int N, gsl_rng *r)
{
    double lMc, q, lq, tc;
    double lqmax, lqmin, dlq, dlMc, dt, alpha;
    double Mc, m1, m2, mt;
    double den;
    double leta;
    int flag;
    int k, i, j;
    
 
    
    for(k = 2; k < NX; k++) params[k] = min[k] + (max[k]-min[k])*gsl_rng_uniform(r);
    
    dlq = log(qfac);
    dlMc = (max[0]-min[0])/(double)(NM);
    dt = (max[5]-min[5])/(double)(N);
    
    lqmin = log(1.0);
    lqmax = (double)(NQ-1)*dlq;
    
    do
    {
    lMc = min[0] + (max[0]-min[0])*gsl_rng_uniform(r);
    lq = lqmin + (lqmax-lqmin)*gsl_rng_uniform(r);
    tc = min[5] + (max[5]-min[5])*gsl_rng_uniform(r);
    
    k = (int)((lq-lqmin)/dlq);
    i = (int)((lMc-min[0])/dlMc);
    j = (int)((tc-min[5])/dt);
        
    den = 0.0;
    if(k >= 0 && i >=0 && j >= 0 && k < NQ && i < NM && j < N) den = global[k][i][j];
    if(den > cap) den = cap;
        
    alpha = gsl_rng_uniform(r)*cap;
        
        
    } while(den < alpha);
    
    params[0] = lMc;
    q = exp(lq);
    Mc = exp(params[0]);
    m2 = Mc*pow(1.0+q, 0.2)/pow(q,0.6);
    m1 = q*m2;
    mt = m1+m2;
    //printf("%f \n", mt/MSUN_SI);
    params[1] = log(mt);   // Mt
    params[5] = Tobs-tc;   // PhenomD puts merger at Tobs, and params[5] tells us how far to pull it back
    

    
    return(den);
    
}


double globeden(double ***global, double *max, double *min, double Tobs, double *params, int N)
{
    double lMc, q, lq, tc;
    double lqmax, lqmin, dlq, dlMc, dt;
    double den, mc, mt, eta, dm;
    int k, i, j;
    
    den = 1.0;
    
    /*
    
    dlq = log(qfac);
    dlMc = (max[0]-min[0])/(double)(NM);
    dt = (max[5]-min[5])/(double)(N);
    
    lqmin = log(1.0);
    lqmax = (double)(NQ-1)*dlq;
    
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
    
    q = (1.0+dm)/(1.0-dm);
    
    lMc = params[0];
    lq = log(q);
    tc = Tobs - params[5]; // PhenomD puts merger at Tobs, and params[5] tells us how far to pull it back
    
    k = (int)((lq-lqmin)/dlq);
    i = (int)((lMc-min[0])/dlMc);
    j = (int)((tc-min[5])/dt);
    
    // printf("%d %d %d %f %f %e %e %e %f\n", k, i, j, q, exp(lMc)/MSUN_SI, min[0], lMc, dlMc, tc);
    
    den = 0.0;
    if(k >= 0 && i >=0 && j >= 0 && k < NQ && i < NM && j < N) den = global[k][i][j];
    
    if(den > cap) den = cap;
     
    */
    
    return(den);
}









void max_array_element(double *max, int *index, double *array, int n)
{
    int i;
    
    *max = array[0];
    *index = 0;
    
    for(i = 1; i <= n-1; i++)
    {
        if(array[i] > *max)
        {
            *max = array[i];
            *index = i;
        }
    }
}




void phase_blind_time_shift(double *corr, double *corrf, double *data1, double *data2, int n)
{
    int nb2, i, l, k, j;
    int imax, imin;
    
    nb2 = n / 2;
    
    corr[0] = 0.0;
    corrf[0] = 0.0;
    corr[nb2] = 0.0;
    corrf[nb2] = 0.0;
    
    for (i=1; i < nb2; i++)
    {
        l=i;
        k=n-i;
        
        corr[l]    = (data1[l]*data2[l] + data1[k]*data2[k]);
        corr[k]    = (data1[k]*data2[l] - data1[l]*data2[k]);
        corrf[l] = corr[k];
        corrf[k] = -corr[l];
    }
    
    gsl_fft_halfcomplex_radix2_inverse(corr, 1, n);
    gsl_fft_halfcomplex_radix2_inverse(corrf, 1, n);
    
    
}




double Fmag(double *sky, double GMST, int id)
{
    double Fplus, Fcross;
    double Ap, Ac;
    double F;
    
    Ap = 0.5*(1.0+sky[3]*sky[3]);
    Ac = -sky[3];
    
    F_ant(sky[2], sky[0], sky[1], GMST, &Fplus, &Fcross, id);
    F = sqrt(Ap*Ap*Fplus*Fplus+Ac*Ac*Fcross*Fcross);
    
    return(F);
    
}




void upsample(int n, double Tobs, int *nt, int *bn)
{
    double dt, dtres;
    int k;
    
    dtres = 1.0e-4;  // desired time resolution
    
    dt = Tobs/(double)(n);
    
    // figure out the upscaling, must be a power of 2
    k = (int)(pow(2.0,ceil(log(dt/dtres)/log(2.0))));
    
    // upsampled
    *bn = n*k;
    
    dt = Tobs/(double)(n*k);
    
    // allow for geocenter - surface light travel time
    // and for attempts to slide the waveform by the max delay time
    *nt = 2*(int)(dtmax/dt)+1;
    
}




double det(double *A, int N)
{
    int *IPIV;
    int LWORK = N*N;
    int INFO;
    int i, j;
    double dx, dy;
    int s;
    
    gsl_permutation *p = gsl_permutation_alloc(N);
    
    gsl_matrix *m = gsl_matrix_alloc (N, N);
    
    for (i = 0 ; i < N ; i++)
    {
        for (j = 0 ; j < N ; j++)
        {
            gsl_matrix_set(m, i, j, A[j*N+i]);
        }
    }
    
    // Compute the LU decomposition of this matrix
    gsl_linalg_LU_decomp(m, p, &s);
    
    dx = 1.0;
    for (i = 0; i < N; i++) dx *= gsl_matrix_get (m, i, i);
    dx = fabs(dx);
    
    //returns the absolute value of the determinant.
    
    gsl_permutation_free(p);
    gsl_matrix_free(m);
    
    
    return dx;
    
    
}




void F_ant(double psi, double alpha, double sindelta, double GMST, double *Fplus, double *Fcross, int id)
{
    if(id == 0) F_LHO(psi, alpha, sindelta, GMST, Fplus, Fcross);
    if(id == 1) F_LLO(psi, alpha, sindelta, GMST, Fplus, Fcross);
    if(id == 2) F_VIRGO(psi, alpha, sindelta, GMST, Fplus, Fcross);
}



void F_LHO(double psi, double alpha, double sindelta, double GMST, double *Fplus, double *Fcross)
{
    double cpsi, spsi, cpsi2;
    double cdelta, sdelta, sdelta2;
    double calpha, salpha, calpha2;
    double F_plus_LHO, F_plus_LLO, F_plus_VIRGO, F_plus_GEO;
    double F_cross_LHO, F_cross_LLO, F_cross_VIRGO, F_cross_GEO;
    
    cpsi = cos(psi);
    spsi = sin(psi);
    cpsi2 = cpsi*cpsi;
    sdelta = sindelta;
    cdelta = sqrt(1.0-sindelta*sindelta);
    sdelta2 = sdelta*sdelta;
    calpha = cos(alpha-GMST);
    calpha2 = calpha*calpha;
    salpha = sin(alpha-GMST);  //This is defined with a minus sign relative to the definitions in Times() - that is correct!
    
    
    *Fplus = .3104536562*cpsi2*salpha*calpha+.4947780910*cdelta*sdelta*calpha+1.424276331*cpsi2*sdelta2*calpha2-.4559956644*cdelta*salpha*sdelta-.1552268274*sdelta2*calpha*salpha+.3104536562*cpsi*spsi*sdelta-.1552268276*salpha*calpha-.4928680677*sdelta2*cpsi2+1.424276334*cpsi2*calpha2+.3104536558*sdelta2*calpha*salpha*cpsi2-.9119913289*cdelta*spsi*cpsi*calpha-.9314082637*cpsi2-.7121381660*calpha2+.9119913288*cdelta*cpsi2*salpha*sdelta-.6209073089*spsi*sdelta*calpha2*cpsi+2.848552667*cpsi*salpha*spsi*sdelta*calpha-.9895561804*cdelta*cpsi2*sdelta*calpha+.2464340337*sdelta2-.9895561811*cdelta*spsi*cpsi*salpha-.7121381653*sdelta2*calpha2+.4657041323;
    
    *Fcross = -1.424276332*spsi*sdelta2*calpha2*cpsi-.3104536551*cpsi*salpha*spsi*calpha-.6209073087*cpsi2*sdelta*calpha2-.9895561820*cdelta*cpsi2*salpha+.3104536561*sdelta*calpha2+.9895561812*cdelta*cpsi*spsi*sdelta*calpha-1.424276332*salpha*sdelta*calpha-.9119913298*cdelta*cpsi2*calpha+2.848552668*cpsi2*salpha*sdelta*calpha-.9119913290*cdelta*cpsi*spsi*salpha*sdelta+.9314082638*cpsi*spsi-1.424276333*cpsi*spsi*calpha2-.1552268276*sdelta+.3104536564*cpsi2*sdelta+.4947780906*cdelta*salpha-.3104536556*spsi*sdelta2*calpha*cpsi*salpha+.4928680679*spsi*sdelta2*cpsi+.4559956636*cdelta*calpha;
    
}

void F_LLO(double psi, double alpha, double sindelta, double GMST, double *Fplus, double *Fcross)
{
    double cpsi, spsi, cpsi2;
    double cdelta, sdelta, sdelta2;
    double calpha, salpha, calpha2;
    double F_plus_LHO, F_plus_LLO, F_plus_VIRGO, F_plus_GEO;
    double F_cross_LHO, F_cross_LLO, F_cross_VIRGO, F_cross_GEO;
    
    cpsi = cos(psi);
    spsi = sin(psi);
    cpsi2 = cpsi*cpsi;
    cdelta = sqrt(1.0-sindelta*sindelta);
    sdelta = sindelta;
    sdelta2 = sdelta*sdelta;
    calpha = cos(alpha-GMST);
    calpha2 = calpha*calpha;
    salpha = sin(alpha-GMST);
    
    
    *Fplus = -.5608410611*cpsi2*salpha*calpha-.4945891833*cdelta*sdelta*calpha-1.040573095*cpsi2*sdelta2*calpha2+.3632312610*cdelta*salpha*sdelta+.2804205304*sdelta2*calpha*salpha-.5608410611*cpsi*spsi*sdelta+.2804205305*salpha*calpha-.3865389475*sdelta2*cpsi2-1.040573095*cpsi2*calpha2-.5608410611*sdelta2*calpha*salpha*cpsi2+.7264625218*cdelta*spsi*cpsi*calpha-.7135560209+1.427112041*cpsi2+.5202865467*calpha2-.7264625224*cdelta*cpsi2*salpha*sdelta+1.121682123*spsi*sdelta*calpha2*cpsi-2.081146187*cpsi*salpha*spsi*sdelta*calpha+.9891783667*cdelta*cpsi2*sdelta*calpha+.1932694762*sdelta2+.9891783666*cdelta*spsi*cpsi*salpha+.5202865474*sdelta2*calpha2;
    
    *Fcross = 1.040573095*spsi*sdelta2*calpha2*cpsi+.5608410611*cpsi*salpha*spsi*calpha+1.121682122*cpsi2*sdelta*calpha2+.9891783668*cdelta*cpsi2*salpha-.5608410611*sdelta*calpha2-.9891783668*cdelta*cpsi*spsi*sdelta*calpha+1.040573095*salpha*sdelta*calpha+.7264625222*cdelta*cpsi2*calpha-2.081146187*cpsi2*salpha*sdelta*calpha+.7264625218*cdelta*cpsi*spsi*salpha*sdelta-1.427112042*cpsi*spsi+1.040573094*cpsi*spsi*calpha2+.2804205304*sdelta-.5608410610*cpsi2*sdelta-.4945891833*cdelta*salpha+.5608410611*spsi*sdelta2*calpha*cpsi*salpha+.3865389466*spsi*sdelta2*cpsi-.3632312614*cdelta*calpha;
    
}

void F_VIRGO(double psi, double alpha, double sindelta, double GMST, double *Fplus, double *Fcross)
{
    double cpsi, spsi, cpsi2;
    double cdelta, sdelta, sdelta2;
    double calpha, salpha, calpha2;
    double F_plus_LHO, F_plus_LLO, F_plus_VIRGO, F_plus_GEO;
    double F_cross_LHO, F_cross_LLO, F_cross_VIRGO, F_cross_GEO;
    
    cpsi = cos(psi);
    spsi = sin(psi);
    cpsi2 = cpsi*cpsi;
    sdelta = sindelta;
    cdelta = sqrt(1.0-sindelta*sindelta);
    sdelta2 = sdelta*sdelta;
    calpha = cos(alpha-GMST);
    calpha2 = calpha*calpha;
    salpha = sin(alpha-GMST);
    
    *Fplus = .3963351156*cpsi2*salpha*calpha+.4651524340*cdelta*sdelta*calpha-1.383399753*cpsi2*sdelta2*calpha2-.3756662059*cdelta*salpha*sdelta-.1981675579*sdelta2*calpha*salpha+.3963351153*cpsi*spsi*sdelta-.1981675578*salpha*calpha+1.303555289*sdelta2*cpsi2-1.383399752*cpsi2*calpha2-0.3992223200e-1+.3963351158*sdelta2*calpha*salpha*cpsi2-.7513324118*cdelta*spsi*cpsi*calpha+0.7984446420e-1*cpsi2+.6916998761*calpha2+.7513324118*cdelta*cpsi2*salpha*sdelta-.7926702330*spsi*sdelta*calpha2*cpsi-2.766799504*cpsi*salpha*spsi*sdelta*calpha-.9303048681*cdelta*cpsi2*sdelta*calpha-.6517776444*sdelta2-.9303048681*cdelta*spsi*cpsi*salpha+.6916998765*sdelta2*calpha2;
    
    *Fcross = 1.383399753*spsi*sdelta2*calpha2*cpsi-.3963351160*cpsi*salpha*spsi*calpha-.7926702326*cpsi2*sdelta*calpha2-.9303048681*cdelta*cpsi2*salpha+.3963351158*sdelta*calpha2+.9303048686*cdelta*cpsi*spsi*sdelta*calpha+1.383399752*salpha*sdelta*calpha-.7513324118*cdelta*cpsi2*calpha-2.766799505*cpsi2*salpha*sdelta*calpha-.7513324118*cdelta*cpsi*spsi*salpha*sdelta-0.7984446460e-1*cpsi*spsi+1.383399752*cpsi*spsi*calpha2-.1981675579*sdelta+.3963351158*cpsi2*sdelta+.4651524340*cdelta*salpha-.3963351158*spsi*sdelta2*calpha*cpsi*salpha-1.303555288*spsi*sdelta2*cpsi+.3756662059*cdelta*calpha;
    
}

void Times(double alpha, double sindelta, double GMST, double *dtimes)
{
    
    double sind, cosd;
    double calpha, salpha;
    double LHO, LLO, VIRGO, GEO;
    
    sind = sindelta;
    cosd = sqrt(1.0-sindelta*sindelta);
    calpha = cos(GMST-alpha);
    salpha = sin(GMST-alpha);
    
    // Note the minus sign flip on y-axis. Matches LAL and is du to sign flip on sin
    
    LHO = -2.161414928e6*cosd*calpha+3.834695183e6*cosd*salpha+4.600350224e6*sind;
    dtimes[0] = -LHO/CLIGHT;
    
    // From LIGO-P000006-C-E
    // H -2.161414928e6, -3.834695183e6, 4.600350224e6
    
    LLO = -74276.04192*cosd*calpha+5.496283721e6*cosd*salpha+3.224257016e6*sind;
    dtimes[1] = -LLO/CLIGHT;
    
    // From LIGO-P000006-C-E
    // L -74276.04192, -5.496283721e6, 3.224257016e6
    
    VIRGO = 4546374.098*cosd*calpha-842989.6972*cosd*salpha+4378576.963*sind;
    dtimes[2] = -VIRGO/CLIGHT;
    
    // V [4.54637409900e+06 8.42989697626e+05 4.37857696241e+06];
    
}



void Ring(double *skyx, double *skyy, int d1, int d2, double GMST, gsl_rng * r)
{
    
    int i;
    double *H1;
    double *L1;
    double *V1;
    double *kv, *kvr;
    double *zv, *uv, *vv;
    double sind, cosd, rot, alpha, x;
    double cr, sr;
    double calpha, salpha;
    double H1mag, L1mag, zmag;
    double cmu, smu;
    
    H1 = (double*)malloc(sizeof(double)* 3);
    L1 = (double*)malloc(sizeof(double)* 3);
    V1 = (double*)malloc(sizeof(double)* 3);
    zv = (double*)malloc(sizeof(double)* 3);
    uv = (double*)malloc(sizeof(double)* 3);
    vv = (double*)malloc(sizeof(double)* 3);
    kv = (double*)malloc(sizeof(double)* 3);
    kvr = (double*)malloc(sizeof(double)* 3);
    
    // From https://dcc.ligo.org/public/0072/P000006/000/P000006-C.pdf  (note sign flip on y axis)
    
    // Location of Hanford vertex from Geocenter in s
    H1[0] = -2.161414928e6/CLIGHT;
    H1[1] = 3.834695183e6/CLIGHT;
    H1[2] =  4.600350224e6/CLIGHT;
    
    //H1mag = sqrt(H1[0]*H1[0]+H1[1]*H1[1]+H1[2]*H1[2]);
    
    // Location of Livingston vertex from Geocenter in s
    L1[0] = -74276.04192/CLIGHT;
    L1[1] = 5.496283721e6/CLIGHT;
    L1[2] =  3.224257016e6/CLIGHT;
    
    //L1mag = sqrt(L1[0]*L1[0]+L1[1]*L1[1]+L1[2]*L1[2]);
    
    // Location of Virgo vertex from Geocenter in s
    V1[0] = 4546374.098/CLIGHT;
    V1[1] = -842989.6972/CLIGHT;
    V1[2] = 4378576.963/CLIGHT;
    
    if(d1 == 0)
    {
        if(d2 == 1)
        {
            for(i = 0; i <3; i++) zv[i] = H1[i]-L1[i];
        }
        if(d2 == 2)
        {
            for(i = 0; i <3; i++) zv[i] = H1[i]-V1[i];
        }
    }
    
    if(d1 == 1)
    {
        if(d2 == 0)
        {
            for(i = 0; i <3; i++) zv[i] = L1[i]-H1[i];
        }
        if(d2 == 2)
        {
            for(i = 0; i <3; i++) zv[i] = L1[i]-V1[i];
        }
    }
    
    if(d1 == 2)
    {
        if(d2 == 0)
        {
            for(i = 0; i <3; i++) zv[i] = V1[i]-H1[i];
        }
        if(d2 == 1)
        {
            for(i = 0; i <3; i++) zv[i] = V1[i]-L1[i];
        }
    }
    
    zmag = 0.0;
    for(i = 0; i <3; i++) zmag += zv[i]*zv[i];
    zmag = sqrt(zmag);
    
    // Unit vector between the sites
    for(i = 0; i <3; i++) zv[i] /= zmag;
    
    alpha = skyx[0];
    sind = skyx[1];
    cosd = sqrt(1.0-sind*sind);
    calpha = cos(alpha-GMST);
    salpha = -sin(alpha-GMST);
    
    // Unusual convention on spherical angle
    // propagation direction
    kv[0] = cosd*calpha;
    kv[1] = cosd*salpha;
    kv[2] = sind;
    
    
    // cosine of angle between propagation direction and vector connecting the sites
    cmu = 0.0;
    for(i = 0; i <3; i++) cmu += zv[i]*kv[i];
    
    if(cmu > 0.999999) cmu = 0.999999;
    if(cmu < -0.999999) cmu = -0.999999;
    
    smu = sqrt(1.0-cmu*cmu);
    
    // vector in plane of prop vector and line connecting sites
    for(i = 0; i <3; i++) uv[i] = (kv[i] - cmu*zv[i])/smu;
    
    //completing the triad
    vv[0] = (zv[1]*kv[2]-zv[2]*kv[1])/smu;
    vv[1] = (zv[2]*kv[0]-zv[0]*kv[2])/smu;
    vv[2] = (zv[0]*kv[1]-zv[1]*kv[0])/smu;
    
    
    // Rotating about zv axis will produce new sky locations with the same time delay between H1 and L1
    rot = TPI*gsl_rng_uniform(r);
    cr = cos(rot);
    sr = sin(rot);
    
    // find rotated k vector
    for(i = 0; i <3; i++) kvr[i] = cmu*zv[i] + smu*(cr*uv[i]+sr*vv[i]);
    
    skyy[1] = kvr[2];
    
    rot = atan2(kvr[1],kvr[0]);
    rot = GMST - rot;
    if(rot < 0.0) rot += TPI;
    
    skyy[0] = rot;
    
    free(H1);
    free(L1);
    free(V1);
    free(kv);
    free(kvr);
    free(zv);
    free(uv);
    free(vv);
    
    
}

void RingPlot(double *paramx, int d1, int d2, double GMST, double *theta, double *phi, int M)
{
    
    int i, j;
    double *H1;
    double *L1;
    double *V1;
    double *kv, *kvr;
    double *zv, *uv, *vv;
    double sind, cosd, rot, alpha, x;
    double cr, sr;
    double calpha, salpha;
    double H1mag, L1mag, zmag;
    double cmu, smu;
    double rd;
    
    H1 = (double*)malloc(sizeof(double)* 3);
    L1 = (double*)malloc(sizeof(double)* 3);
    V1 = (double*)malloc(sizeof(double)* 3);
    zv = (double*)malloc(sizeof(double)* 3);
    uv = (double*)malloc(sizeof(double)* 3);
    vv = (double*)malloc(sizeof(double)* 3);
    kv = (double*)malloc(sizeof(double)* 3);
    kvr = (double*)malloc(sizeof(double)* 3);
    
    rd = 180.0/PI;
    
    // From https://dcc.ligo.org/public/0072/P000006/000/P000006-C.pdf  (note sign flip on y axis)
    
    // Location of Hanford vertex from Geocenter in s
    H1[0] = -2.161414928e6/CLIGHT;
    H1[1] = 3.834695183e6/CLIGHT;
    H1[2] =  4.600350224e6/CLIGHT;
    
    //H1mag = sqrt(H1[0]*H1[0]+H1[1]*H1[1]+H1[2]*H1[2]);
    
    // Location of Livingston vertex from Geocenter in s
    L1[0] = -74276.04192/CLIGHT;
    L1[1] = 5.496283721e6/CLIGHT;
    L1[2] =  3.224257016e6/CLIGHT;
    
    //L1mag = sqrt(L1[0]*L1[0]+L1[1]*L1[1]+L1[2]*L1[2]);
    
    // Location of Virgo vertex from Geocenter in s
    V1[0] = 4546374.098/CLIGHT;
    V1[1] = -842989.6972/CLIGHT;
    V1[2] = 4378576.963/CLIGHT;
    

    
    if(d1 == 0)
    {
        if(d2 == 1)
        {
            for(i = 0; i <3; i++) zv[i] = H1[i]-L1[i];
        }
        if(d2 == 2)
        {
            for(i = 0; i <3; i++) zv[i] = H1[i]-V1[i];
        }
    }
    
    if(d1 == 1)
    {
        if(d2 == 0)
        {
            for(i = 0; i <3; i++) zv[i] = L1[i]-H1[i];
        }
        if(d2 == 2)
        {
            for(i = 0; i <3; i++) zv[i] = L1[i]-V1[i];
        }
    }
    
    if(d1 == 2)
    {
        if(d2 == 0)
        {
            for(i = 0; i <3; i++) zv[i] = V1[i]-H1[i];
        }
        if(d2 == 1)
        {
            for(i = 0; i <3; i++) zv[i] = V1[i]-L1[i];
        }
    }
    
    
    zmag = 0.0;
    for(i = 0; i <3; i++) zmag += zv[i]*zv[i];
    zmag = sqrt(zmag);
    
    // Unit vector between the sites
    for(i = 0; i <3; i++) zv[i] /= zmag;
    
    // [0] log(Mc) [1] log(Mt) [2] chi1 [3] chi2 [4] phi0 [5] tp [6] log(DL) [7] alpha [8] sindelta [9] psi [10] ciota
    
    alpha = paramx[7];
    sind = paramx[8];
    cosd = sqrt(1.0-sind*sind);
    calpha = cos(alpha-GMST);
    salpha = -sin(alpha-GMST);
    
    // Unusual convention on spherical angle
    // propagation direction
    kv[0] = cosd*calpha;
    kv[1] = cosd*salpha;
    kv[2] = sind;
    
    
    // cosine of angle between propagation direction and vector connecting the sites
    cmu = 0.0;
    for(i = 0; i <3; i++) cmu += zv[i]*kv[i];
    
    if(cmu > 0.999999) cmu = 0.999999;
    if(cmu < -0.999999) cmu = -0.999999;
    
    smu = sqrt(1.0-cmu*cmu);
    
    // vector in plane of prop vector and line connecting sites
    for(i = 0; i <3; i++) uv[i] = (kv[i] - cmu*zv[i])/smu;
    
    //completing the triad
    vv[0] = (zv[1]*kv[2]-zv[2]*kv[1])/smu;
    vv[1] = (zv[2]*kv[0]-zv[0]*kv[2])/smu;
    vv[2] = (zv[0]*kv[1]-zv[1]*kv[0])/smu;
    
    
    for(i = 0; i <= M; i++)
    {
    rot = (double)(i)*TPI/(double)(M);
    cr = cos(rot);
    sr = sin(rot);
    
    // find rotated k vector
    for(j = 0; j <3; j++) kvr[j] = cmu*zv[j] + smu*(cr*uv[j]+sr*vv[j]);
        
    
    theta[i] = acos(kvr[2]);
    
    rot = atan2(kvr[1],kvr[0]);
    rot = GMST - rot;
    if(rot < 0.0) rot += TPI;
    if(rot > TPI) rot -= TPI;
    
    phi[i] = rot;
  
    }
    
    free(H1);
    free(L1);
    free(V1);
    free(kv);
    free(kvr);
    free(zv);
    free(uv);
    free(vv);
    
    
}


double f_nwip(double *a, double *b, int n)
{
    int i, j, k;
    double arg, product;
    double test;
    double ReA, ReB, ImA, ImB;
    
    arg = 0.0;
    for(i=1; i<n/2; i++)
    {
        j = i;
        k = n-j;
        ReA = a[j]; ImA = a[k];
        ReB = b[j]; ImB = b[k];
        product = ReA*ReB + ImA*ImB;
        arg += product;
    }
    
    return(arg);
    
}



void Ring_all(double *paramx, double *paramy, int d1, int d2, double GMST, gsl_rng * r)
{
    
    int i;
    double *H1;
    double *L1;
    double *V1;
    double *kv, *kvr;
    double *zv, *uv, *vv;
    double sind, cosd, rot, alpha, x;
    double cr, sr;
    double calpha, salpha;
    double H1mag, L1mag, zmag;
    double cmu, smu;
    
    H1 = (double*)malloc(sizeof(double)* 3);
    L1 = (double*)malloc(sizeof(double)* 3);
    V1 = (double*)malloc(sizeof(double)* 3);
    zv = (double*)malloc(sizeof(double)* 3);
    uv = (double*)malloc(sizeof(double)* 3);
    vv = (double*)malloc(sizeof(double)* 3);
    kv = (double*)malloc(sizeof(double)* 3);
    kvr = (double*)malloc(sizeof(double)* 3);
    
    // From https://dcc.ligo.org/public/0072/P000006/000/P000006-C.pdf  (note sign flip on y axis)
    
    // Location of Hanford vertex from Geocenter in s
    H1[0] = -2.161414928e6/CLIGHT;
    H1[1] = 3.834695183e6/CLIGHT;
    H1[2] =  4.600350224e6/CLIGHT;
    
    //H1mag = sqrt(H1[0]*H1[0]+H1[1]*H1[1]+H1[2]*H1[2]);
    
    // Location of Livingston vertex from Geocenter in s
    L1[0] = -74276.04192/CLIGHT;
    L1[1] = 5.496283721e6/CLIGHT;
    L1[2] =  3.224257016e6/CLIGHT;
    
    //L1mag = sqrt(L1[0]*L1[0]+L1[1]*L1[1]+L1[2]*L1[2]);
    
    // Location of Virgo vertex from Geocenter in s
    V1[0] = 4546374.098/CLIGHT;
    V1[1] = -842989.6972/CLIGHT;
    V1[2] = 4378576.963/CLIGHT;
    
    if(d1 == 0)
    {
        if(d2 == 1)
        {
            for(i = 0; i <3; i++) zv[i] = H1[i]-L1[i];
        }
        if(d2 == 2)
        {
            for(i = 0; i <3; i++) zv[i] = H1[i]-V1[i];
        }
    }
    
    if(d1 == 1)
    {
        if(d2 == 0)
        {
            for(i = 0; i <3; i++) zv[i] = L1[i]-H1[i];
        }
        if(d2 == 2)
        {
            for(i = 0; i <3; i++) zv[i] = L1[i]-V1[i];
        }
    }
    
    if(d1 == 2)
    {
        if(d2 == 0)
        {
            for(i = 0; i <3; i++) zv[i] = V1[i]-H1[i];
        }
        if(d2 == 1)
        {
            for(i = 0; i <3; i++) zv[i] = V1[i]-L1[i];
        }
    }
    
    
    zmag = 0.0;
    for(i = 0; i <3; i++) zmag += zv[i]*zv[i];
    zmag = sqrt(zmag);
    
    // Unit vector between the sites
    for(i = 0; i <3; i++) zv[i] /= zmag;
    
    // [0] log(Mc) [1] log(Mt) [2] chi1 [3] chi2 [4] phi0 [5] tp [6] log(DL) [7] alpha [8] sindelta [9] psi [10] ciota
    
    alpha = paramx[7];
    sind = paramx[8];
    cosd = sqrt(1.0-sind*sind);
    calpha = cos(alpha-GMST);
    salpha = -sin(alpha-GMST);
    
    // Unusual convention on spherical angle
    // propagation direction
    kv[0] = cosd*calpha;
    kv[1] = cosd*salpha;
    kv[2] = sind;
    
    
    // cosine of angle between propagation direction and vector connecting the sites
    cmu = 0.0;
    for(i = 0; i <3; i++) cmu += zv[i]*kv[i];
    
    if(cmu > 0.999999) cmu = 0.999999;
    if(cmu < -0.999999) cmu = -0.999999;
    
    smu = sqrt(1.0-cmu*cmu);
    
    // vector in plane of prop vector and line connecting sites
    for(i = 0; i <3; i++) uv[i] = (kv[i] - cmu*zv[i])/smu;
    
    //completing the triad
    vv[0] = (zv[1]*kv[2]-zv[2]*kv[1])/smu;
    vv[1] = (zv[2]*kv[0]-zv[0]*kv[2])/smu;
    vv[2] = (zv[0]*kv[1]-zv[1]*kv[0])/smu;
    
    
    // Rotating about zv axis will produce new sky locations with the same time delay between H1 and L1
    rot = TPI*gsl_rng_uniform(r);
    cr = cos(rot);
    sr = sin(rot);
    
    // find rotated k vector
    for(i = 0; i <3; i++) kvr[i] = cmu*zv[i] + smu*(cr*uv[i]+sr*vv[i]);
    
    paramy[8] = kvr[2];
    
    rot = atan2(kvr[1],kvr[0]);
    rot = GMST - rot;
    if(rot < 0.0) rot += TPI;
    
    paramy[7] = rot;
    
    
    free(H1);
    free(L1);
    free(V1);
    free(kv);
    free(kvr);
    free(zv);
    free(uv);
    free(vv);
    
    
}


double z_DL(double DL)
{
    double zlow, zhigh;
    double zmid;
    double Dlow, Dhigh, Dmid;
    double u, v, w;
    int i;
    
    zlow = 1.0;
    zhigh = 1000.0;
    
    if(DL < 1.0)
    {
        zlow = 1.0e-10;
        zhigh = 3.0e-4;
    }
    else if (DL < 100.0)
    {
        zlow = 2.0e-4;
        zhigh = 5.0e-2;
    }
    else if (DL < 1000.0)
    {
        zlow = 1.0e-2;
        zhigh = 3.0e-1;
    }
    else if (DL < 10000.0)
    {
        zlow = 1.0e-1;
        zhigh = 2.0;
    }
    
    zmid = (zhigh+zlow)/2.0;
    
    Dlow = DL_fit(zlow);
    Dhigh = DL_fit(zhigh);
    Dmid = DL_fit(zmid);
    
    for (i=0; i< 30; i++)
    {
        
        u = Dlow-DL;
        v = Dmid-DL;
        
        if(u*v < 0.0)
        {
            zhigh = zmid;
            Dhigh = Dmid;
        }
        else
        {
            zlow = zmid;
            Dlow = Dmid;
        }
        zmid = (zhigh+zlow)/2.0;
        Dmid = DL_fit(zmid);
        
        
    }
    
    //  printf("%e %e\n", Dmid-DL, zmid);
    
    return(zmid);
    
    
}


double DL_fit(double z)
{
    double D;
    // Planck values https://arxiv.org/abs/1807.06209
    double Om = 0.315;
    double H0 = 67.4;
    double RH, x1, x2, Om1, Om2;
    
    // See http://arxiv.org/pdf/1111.6396v1.pdf
    
    RH = CLIGHT/H0*(1.0e-3*MPC);
    
    x1 = (1.0-Om)/(Om);
    
    x2 = (1.0-Om)/(Om*pow((1.0+z),3.0));
    
    Om1 = (1.0+1.32*x1+0.4415*x1*x1+0.02656*x1*x1*x1)/(1.0+1.392*x1+0.5121*x1*x1+0.03944*x1*x1*x1);
    
    Om2 = (1.0+1.32*x2+0.4415*x2*x2+0.02656*x2*x2*x2)/(1.0+1.392*x2+0.5121*x2*x2+0.03944*x2*x2*x2);
    
    D = 2.0*RH*(1.0+z)/sqrt(Om)*(Om1-Om2/sqrt(1.0+z));
    
    return(D/MPC);
}


void TransformC(double *a, double *freqs, double **tf, double **tfR, double **tfI, double Q, double Tobs, int n, int m)
{
    int j;
    double fix;
    
    // [0] t0 [1] f0 [2] Q [3] Amp [4] phi
    
    fix = sqrt((double)(n/2));
    
    #pragma omp parallel for
    for(j = 0; j < m; j++)
    {
        layerC(a, freqs[j], tf[j], tfR[j], tfI[j], Q, Tobs, fix, n);
    }
    
}


void layerC(double *a, double f, double *tf, double *tfR, double *tfI, double Q, double Tobs, double fix, int n)
{
    int i;
    double *AC, *AF;
    double *b;
    double *params;
    double bmag;
    
    params= double_vector(6);
    AC=double_vector(n);  AF=double_vector(n);
    b = double_vector(n);
    
    params[0] = 0.0;
    params[1] = f;
    params[2] = Q;
    params[3] = 1.0;
    params[4] = 0.0;
    
    SineGaussianC(b, params, Tobs, n);
    
    bmag = sqrt(f_nwip(b, b, n)/(double)n);
    
    bmag /= fix;
    
    phase_blind_time_shift(AC, AF, a, b, n);
    
    for(i = 0; i < n; i++)
    {
        tfR[i] = AC[i]/bmag;
        tfI[i] = AF[i]/bmag;
        tf[i] = tfR[i]*tfR[i]+tfI[i]*tfI[i];
    }
    
    free_double_vector(AC);  free_double_vector(AF);
    free_double_vector(b);  free_double_vector(params);
    
}

void SineGaussianC(double *hs, double *sigpar, double Tobs, int N)
{
    double f0, t0, Q, sf, sx, Amp;
    double fmx, fmn;//, fac;
    double phi, f;
    double tau;
    double re,im;
    double q, p, u;
    double A, B, C;
    double Am, Bm, Cm;
    double a, b, c;
    double dt, fac;
    double cosPhase_m, sinPhase_m, cosPhase_p, sinPhase_p;
    
    int i, imid, istart,istop,imin,imax,iend,even,odd;
    
    t0  = sigpar[0];
    f0  = sigpar[1];
    Q   = sigpar[2];
    Amp = sigpar[3];
    phi = sigpar[4];
    
    tau = Q/(TPI*f0);
    
    fmx = f0 + 3.0/tau;  // no point evaluating waveform past this time (many efolds down)
    fmn = f0 - 3.0/tau;  // no point evaluating waveform before this time (many efolds down)
    
    i = (int)(f0*Tobs);
    imin = (int)(fmn*Tobs);
    imax = (int)(fmx*Tobs);
    if(imax - imin < 10)
    {
        imin = i-5;
        imax = i+5;
    }
    
    if(imin < 0) imin = 1;
    if(imax > N/2) imax = N/2;
    
    hs[0] = 0.0;
    hs[N/2] = 0.0;
    
    for(i = 1; i < N/2; i++)
    {
        hs[i] = 0.0;
        hs[N-i] = 0.0;
    }
    
    dt = Tobs/(double)(N);
    fac = sqrt(sqrt(2.0)*PI*tau/dt);
    
    /* Use recursion relationship  */
    
    imid = (int)(f0*Tobs);
    
    p = PI*PI*tau*tau/(Tobs*Tobs);
    
    Bm = exp(-p*(((double)(imid)-f0*Tobs)*((double)(imid)-f0*Tobs)));
    Cm = 1.0;
    
    b = exp(-p*(1.0+2.0*((double)(imid)-f0*Tobs)));
    c = exp(-2.0*p);
    
    // start in the middle and work outwards
    
    B = Bm;
    C = Cm;

    
    for(i = imid; i < imax; i++)
    {
 
        f = (double)(i)/Tobs;
        
        sf = fac*B;
        
        //  printf("%e\n", exp(-pi2tau2*(f-f0)*(f-f0))-B);
        
        hs[i] = sf;
        hs[N-i] = 0.0;
        
        B *= (C*b);
        C *= c;
        
    }
    
    // reset to midpoint
    
    B = Bm;
    C = Cm;
    
    b = exp(p*(-1.0+2.0*((double)(imid)-f0*Tobs)));
    // c unchanged
    
    for(i = imid; i > imin; i--)
    {

        f = (double)(i)/Tobs;
        
        sf = fac*B;
        
        // printf("%e\n", exp(-pi2tau2*(f-f0)*(f-f0))-B);
        
        hs[i] = sf;
        hs[N-i] = 0.0;
        
        B *= (C*b);
        C *= c;
        

    }
    
    
}


void specest(double *data, double *Hf, int N, int Ns, double dt, double fmx, double *SN, double *SM, double *PS)
{
    int i, j, k, M, Nf, Nstep,  ii, m;
    int jj, kk, Nlines;
    int oflag, flag;
    int imin, imax;
    double SNR, max;
    double junk, Tobs, fix, f, t, t0, df, x, y, z, dx;
    double fmn, dfx, Q, fny, delt, scale, dlnf;
    double *freqs, *ref;
    double *Draw;
    double *D, *times;
    double *Sn;
    double *specD, *sspecD;
    double *sdata;
    double *intime, *sqf;
    double *pspline, *dspline;
    double sigmean, sigmedian;
    int subscale, octaves;
    int mmax;
    double SNRsq, SNRold, pH, pL, pmax;
    double SNRH, SNRL, pw, alpha;
    double t_rise, s1, s2, ascale, fac, Dfmax;
    double *linew;
    int pflag;
    
    int modelprint;

    char filename[1024];
    char command[1024];
    
    FILE *out;
    
    Q = Qs;  // Q of transform
    Dfmax = 8.0;  // width of smoothing window in Hz
    Tobs = (double)(N)*dt;  // duration
    
    // Tukey window parameter. Flat for (1-alpha) of data
    t_rise = 0.4; // Standard LAL setting
    alpha = (2.0*t_rise/Tobs);
    
    
    
    df = 1.0/Tobs;  // frequency resolution
    fny = 1.0/(2.0*dt);  // Nyquist
    
    // Choose the range of the spectrogram. fmin, fmax must be a power of 2
    fmn = 8.0;
    if(fmx > fny) fmx = fny;
    
    D = (double*)malloc(sizeof(double)* (N));
    Draw = (double*)malloc(sizeof(double)* (N));
    
    // Copy data over
    for(i=0; i< N; i++)
    {
     Draw[i] = data[i];
     D[i] = data[i];
    }
    
    // Normalization factor
    tukey_scale(&s1, &s2, alpha, N);
    
    printf("%f %f\n", s1, s2);
    
    // Time series data is corrupted by the Tukey window at either end
    // imin and imax define the "safe" region to use
    imin = (int)(2.0*t_rise/dt);
    imax = N-imin;
    
    
    // Prepare to make spectogram
    
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
    
    
    
    printf("frequency layers = %d\n", Nf);
    
    scale = Getscale(freqs, Q, Tobs, fmx, N, Nf);
    
    sspecD = (double*)malloc(sizeof(double)*(N/2));
    specD = (double*)malloc(sizeof(double)*(N/2));
    Sn = (double*)malloc(sizeof(double)*(N/2));
    

    SNRold = 0.0;
    clean(D, Draw, Hf, sqf, freqs, Sn, specD, sspecD, df, Q, Tobs, scale, alpha, Nf, N, imin, imax, &SNR, 0);
    
    
    // if big glitches are present we need to rinse and repeat
    i = 0;
    while(i < 10 && (SNR-SNRold) > 10.0)
    {
        SNRold = SNR;
        clean(D, Draw, Hf, sqf, freqs, Sn, specD, sspecD, df, Q, Tobs, scale, alpha, Nf, N, imin, imax, &SNR, 0);
        i++;
    }
    
    
    // pass back the cleaned data. This cleaned data is not passed forward
    // from the spectral estimation code. Only used inside the MCMC
    for(i=0; i< N; i++)
    {
     data[i] = D[i];
    }
    
    // re-compute the power spectrum using the cleaned data
    tukey(D, alpha, N);
    
    gsl_fft_real_radix2_transform(D, 1, N);
    
    // Form spectral model for whitening data (lines plus a smooth component)
    spectrum(D, Sn, specD, sspecD, df, N);
    
    fac = Tobs/((double)(N)*(double)(N));
    
    for (i = 0; i < Ns; ++i) SM[i] = sspecD[i]*fac;
    for (i = 0; i < Ns; ++i) SN[i] = specD[i]*fac;
    for (i = 0; i < Ns; ++i) PS[i] = Sn[i]*fac;
    
    // print the cleaned data
    for(i=0; i< N; i++)
    {
     D[i] = data[i];
    }
    
    // save Q-scan data to file
    clean(D, Draw, Hf, sqf, freqs, Sn, specD, sspecD, df, Qprint, Tobs, scale, alpha, Nf, N, imin, imax, &SNR, 1);
    
    
    free(D);
    free(Draw);
    free(sspecD);
    free(specD);
    free(Sn);
    free(freqs);
    free(sqf);

    
}


void CubicSplineGSL(int N, double *x, double *y, int Nint, double *xint, double *yint)
{
    int n;
    double tol=1.0e-6;
    
    
    /* set up GSL cubic spline */
   // gsl_spline       *cspline = gsl_spline_alloc(gsl_interp_cspline, N);
    gsl_spline *cspline = gsl_spline_alloc(gsl_interp_steffen, N);
    gsl_interp_accel *acc    = gsl_interp_accel_alloc();
    
    /* get derivatives */
    gsl_spline_init(cspline,x,y,N);
    
    /* interpolate */
    for(n=0; n<Nint; n++)
    {
        /*
         GSL cubic spline throws errors if
         interpolated points are at end of
         spline control points
         */
        if     (fabs(xint[n]-x[0])<tol)
            yint[n] = y[0];
        
        else if(fabs(xint[n]-x[N-1])<tol)
            yint[n] = y[N-1];
        
        else
            yint[n]=gsl_spline_eval(cspline,xint[n],acc);
    }
    
    gsl_spline_free (cspline);
    gsl_interp_accel_free (acc);
    
}

void clean(double *D, double *Draw, double *Hf, double *sqf, double *freqs, double *Sn, double *specD, double *sspecD, double df, double Q, double Tobs, double scale, double alpha, int Nf, int N, int imin, int imax, double *SNR, int pflag)
{
    
    int i, j, k;
    int flag;
    int ii, jj;
    double x, y, dt;
    double S;
    double fac;
    double t, f, fmn;
    
    // allocate some arrays
    double **tfDR, **tfDI;
    double **tfD;
    double **live;
    double **live2;
    
    double *Dtemp, *Drf;
    
    FILE *out;
    
    live= double_matrix(Nf,N);
    live2= double_matrix(Nf,N);
    tfDR = double_matrix(Nf,N);
    tfDI = double_matrix(Nf,N);
    tfD = double_matrix(Nf,N);
    
    dt = Tobs/(double)(N);
    
    
    Dtemp = (double*)malloc(sizeof(double)*(N));
    
    for (i = 0; i < N; i++) Dtemp[i] = Draw[i];
    
    // D holds the previously cleaned time-domain data
    // Draw holds the raw time-domain data
    
    // D is used to compute the spectrum. A copy of Draw, Dtemp, is then whitened using this spectrum and glitches are then identified
    // The glitches are then re-colored, subtracted from Draw to become the new D
 
 
    // Tukey window
    tukey(D, alpha, N);
    tukey(Dtemp, alpha, N);
   
    
    // FFT
    gsl_fft_real_radix2_transform(D, 1, N);
    gsl_fft_real_radix2_transform(Dtemp, 1, N);
    
    // remove the CBC signal if provided (in the spec code this is zero)
    for(i = 0; i < N; i++)
    {
        D[i] -= Hf[i];
        Dtemp[i] -= Hf[i];
    }
    
    // Form spectral model for whitening data (lines plus a smooth component)
    spectrum(D, Sn, specD, sspecD, df, N);
    
    // whiten data
    whiten(Dtemp, specD, N);
    
    // Wavelet transform
    TransformC(Dtemp, freqs, tfD, tfDR, tfDI, Q, Tobs, N, Nf);
    
    if(pflag == 1)
    {
        out = fopen("Qtransform.dat","w");
        for(j = 0; j < Nf; j++)
        {
            f = freqs[j];
            
            for(i = 0; i < N; i++)
            {
                t = (double)(i)*dt-Tobs+2.0;  // trigger time is set two seconds from end, set zero there
                fprintf(out,"%e %e %e\n", t, f, tfD[j][i]);
            }
            
            fprintf(out,"\n");
        }
        fclose(out);
        
    }
    
    k = 0;
    //  apply threshold
    for(j = 0; j < Nf; j++)
    {
        for (i = 0; i < N; i++)
        {
            live[j][i] = -1.0;
            if(tfD[j][i] > sthresh) live[j][i] = 1.0;
            if(tfD[j][i] > sthresh) k++;
            live2[j][i] = live[j][i];
        }
    }
   
    
    // dig deeper to extract clustered power
    for(j = 1; j < Nf-1; j++)
    {
        for (i = 1; i < N-1; i++)
        {
            
            flag = 0;
            for(jj = -1; jj <= 1; jj++)
            {
                for(ii = -1; ii <= 1; ii++)
                {
                    if(live[j+jj][i+ii] > 0.0) flag = 1;
                }
            }
            if(flag == 1 && tfD[j][i] > warm) live2[j][i] = 1.0;
        }
    }
    
    
    // build the excess power model
    for (i = 0; i < N; i++)
    {
        Dtemp[i] = 0.0;
    }
    
    k = 0;
    for(j = 0; j < Nf; j++)
    {
        for (i = imin; i < imax; i++)
        {
            if(live2[j][i] > 0.0) Dtemp[i] += scale*sqf[j]*tfDR[j][i];
        }
    }
    
    if(pflag == 1)
       {
        out = fopen("wglitch.dat", "w");
        for(i=0; i< N; i++)
        {
            fprintf(out,"%e %e\n", (double)(i)*dt-Tobs+2.0, Dtemp[i]);
         }
       fclose(out);
       }
    
    // Compute the excess power (relative to the current spectral model
    S = 0.0;
    for (i = imin; i < imax; ++i) S += Dtemp[i]*Dtemp[i];
    S = sqrt(S);
    
    printf("Excess SNR = %f\n", S);
    
    *SNR = S;
    

    //Unwhiten and subtract the excess power so we can compute a better spectral estimate
    // Back to frequency domain
    
    gsl_fft_real_radix2_transform(Dtemp, 1, N);
    
    
    fmn = freqs[0];
    
    // only use smooth spectrum in the un-whitening
    Dtemp[0] = 0.0;
    for(i=1; i< N/2; i++)
    {
        f = (double)(i)/Tobs;
        //y = 1.0;
        //if(f < fmn) y = 0.5*(1.0+tanh(8.0*(f-0.5*fmn)/fmn));
        //x = y*sqrt(sspecD[i]);
        x = sqrt(sspecD[i]);
        Dtemp[i] *= x;
        Dtemp[N-i] *= x;
    }
    Dtemp[N/2] = 0.0;
    
    gsl_fft_halfcomplex_radix2_inverse(Dtemp, 1, N);
    
    
    x = sqrt((double)(2*N));
    
    // note that this D[i] still contains the signal since it
    // is formed from the raw data minus the glitch model
    for(i=0; i< N; i++)
    {
        D[i] = Draw[i]-Dtemp[i]/x;
    }
    
    /*
    out = fopen("glitch.dat", "w");
    for(i=0; i< N; i++)
    {
        fprintf(out,"%d %e %e\n", i, Dtemp[i]/x, Draw[i]);
    }
    fclose(out);
    */
    
    free(Dtemp);
    free_double_matrix(live,Nf);
    free_double_matrix(live2,Nf);
    free_double_matrix(tfDR,Nf);
    free_double_matrix(tfDI,Nf);
    free_double_matrix(tfD,Nf);
    
    
    return;
    
    
}

void spectrum(double *data, double *S, double *Sn, double *Smooth, double df, int N)
{
    double Df, Dfmax, x, y;
    double Df1, Df2;
    int mw, k, i, j;
    int mm, kk;
    int end1, end2, end3;
    double med;
    double *chunk;
    
    // log(2) is median/2 of chi-squared with 2 dof
    
    
    for(i=1; i< N/2; i++) S[i] = 2.0*(data[i]*data[i]+data[N-i]*data[N-i]);
    S[0] = S[1];
    
    Dfmax = 16.0; // is the  width of smoothing window in Hz
    
    // Smaller windows used initially where the spectrum is steep
    Df2 = Dfmax/2.0;
    Df1 = Dfmax/4.0;
    
    // defines the ends of the segments where smaller windows are used
    end1 = (int)(32.0/df);
    end2 = 2*end1;
    
    mw = (int)(Dfmax/df)+1;  // size of median window
    //printf("numer of bins in smoothing window %d\n", mw);
    k = (mw+1)/2;
    chunk = double_vector(mw);
    
    end3 = N/2-k;  // end of final chunk
    
    // Fill the array so the ends are not empty - just to be safe
    for(i=0;i< N/2;i++)
    {
        Sn[i] = S[i];
        Smooth[i] = S[i];
    }
    

    mw = (int)(Df1/df)+1;  // size of median window
    k = (mw+1)/2;
    
    for(i=4; i< k; i++)
    {
        mm = i/2;
        kk = (mm+1)/2;
        
        for(j=0;j< mm;j++)
        {
            chunk[j] = S[i-kk+j];
        }
        
        gsl_sort(chunk, 1, mm);
        Sn[i] = gsl_stats_median_from_sorted_data(chunk, 1, mm)/LN2;
    
        Smooth[i] = Sn[i];
        
    }
    
    
    i = k;
    do
    {
        for(j=0;j< mw;j++)
        {
            chunk[j] = S[i-k+j];
        }
        
        gsl_sort(chunk, 1, mw);
        Sn[i] = gsl_stats_median_from_sorted_data(chunk, 1, mw)/LN2;
        
        Smooth[i] = Sn[i];
        
        i++;
        
    }while(i < end1);
    
    
    
    mw = (int)(Df2/df)+1;  // size of median window
    k = (mw+1)/2;
    
    
    do
    {
        for(j=0;j< mw;j++)
        {
            chunk[j] = S[i-k+j];
        }
        
        gsl_sort(chunk, 1, mw);
        Sn[i] = gsl_stats_median_from_sorted_data(chunk, 1, mw)/LN2;
        
        Smooth[i] = Sn[i];
        
        i++;
        
    }while(i < end2);
    
    mw = (int)(Dfmax/df)+1;  // size of median window
    k = (mw+1)/2;
    
    do
    {
        for(j=0;j< mw;j++)
        {
            chunk[j] = S[i-k+j];
        }
        
        gsl_sort(chunk, 1, mw);
        Sn[i] = gsl_stats_median_from_sorted_data(chunk, 1, mw)/LN2;
    
        Smooth[i] = Sn[i];
        
        i++;
        
    }while(i < end3);
    
    
    for(i=end3; i< N/2-4; i++)
    {
        mm = (N/2-i)/2;
        kk = (mm+1)/2;
        
        for(j=0;j< mm;j++)
        {
            chunk[j] = S[i-kk+j];
        }
        
        gsl_sort(chunk, 1, mw);
        Sn[i] = gsl_stats_median_from_sorted_data(chunk, 1, mm)/LN2;
        
        Smooth[i] = Sn[i];
        
    }
    
    
    free_double_vector(chunk);
    
    
    
    
    // zap the lines.
    for(i=1;i< N/2;i++)
    {
        x = S[i]/Sn[i];
        if(x > 10.0)
        {
            Sn[i] = S[i];
        }
    }
    
    
    
}



void bwbpf(double *in, double *out, int fwrv, int M, int n, double s, double f1, double f2)
{
    /* Butterworth bandpass filter
     n = filter order 4,8,12,...
     s = sampling frequency
     f1 = upper half power frequency
     f2 = lower half power frequency  */
    
    if(n % 4){ printf("Order must be 4,8,12,16,...\n"); return;}
    
    int i, j;
    double a = cos(PI*(f1+f2)/s)/cos(PI*(f1-f2)/s);
    double a2 = a*a;
    double b = tan(PI*(f1-f2)/s);
    double b2 = b*b;
    double r;
    
    n = n/4;
    double *A = (double *)malloc(n*sizeof(double));
    double *d1 = (double *)malloc(n*sizeof(double));
    double *d2 = (double *)malloc(n*sizeof(double));
    double *d3 = (double *)malloc(n*sizeof(double));
    double *d4 = (double *)malloc(n*sizeof(double));
    double *w0 = (double *)malloc(n*sizeof(double));
    double *w1 = (double *)malloc(n*sizeof(double));
    double *w2 = (double *)malloc(n*sizeof(double));
    double *w3 = (double *)malloc(n*sizeof(double));
    double *w4 = (double *)malloc(n*sizeof(double));
    double x;
    
    for(i=0; i<n; ++i)
    {
        r = sin(PI*(2.0*(double)i+1.0)/(4.0*(double)n));
        s = b2 + 2.0*b*r + 1.0;
        A[i] = b2/s;
        d1[i] = 4.0*a*(1.0+b*r)/s;
        d2[i] = 2.0*(b2-2.0*a2-1.0)/s;
        d3[i] = 4.0*a*(1.0-b*r)/s;
        d4[i] = -(b2 - 2.0*b*r + 1.0)/s;
        w0[i] = 0.0;
        w1[i] = 0.0;
        w2[i] = 0.0;
        w3[i] = 0.0;
        w4[i] = 0.0;
    }
    
    for(j=0; j< M; ++j)
    {
        if(fwrv == 1) x = in[j];
        if(fwrv == -1) x = in[M-j-1];
        for(i=0; i<n; ++i)
        {
            w0[i] = d1[i]*w1[i] + d2[i]*w2[i]+ d3[i]*w3[i]+ d4[i]*w4[i] + x;
            x = A[i]*(w0[i] - 2.0*w2[i] + w4[i]);
            w4[i] = w3[i];
            w3[i] = w2[i];
            w2[i] = w1[i];
            w1[i] = w0[i];
        }
        if(fwrv == 1) out[j] = x;
        if(fwrv == -1) out[M-j-1] = x;
    }
    
    free(A);
    free(d1);
    free(d2);
    free(d3);
    free(d4);
    free(w0);
    free(w1);
    free(w2);
    free(w3);
    free(w4);
    
    return;
}




double Getscale(double *freqs, double Q, double Tobs, double fmx, int n, int m)
{
    double *data, *intime, *ref, **tfR, **tfI, **tf;
    double f, t0, delt, t, x, fix, dt;
    double scale, sqf;
    int i, j;
    
    FILE *out;
    
    data = double_vector(n);
    ref = double_vector(n);
    intime = double_vector(n);
    
    tf = double_matrix(m,n);
    tfR = double_matrix(m,n);
    tfI = double_matrix(m,n);
    
    f = fmx/4.0;
    t0 = Tobs/2.0;
    delt = Tobs/8.0;
    dt = Tobs/(double)(n);
    
    //out = fopen("packet.dat","w");
    for(i=0; i< n; i++)
    {
        t = (double)(i)*dt;
        x = (t-t0)/delt;
        x = x*x/2.0;
        data[i] = cos(TPI*t*f)*exp(-x);
        ref[i] = data[i];
        //fprintf(out,"%e %e\n", t, data[i]);
    }
    // fclose(out);
    
    gsl_fft_real_radix2_transform(data, 1, n);

    TransformC(data, freqs, tf, tfR, tfI, Q, Tobs, n, m);
    
    for(i = 0; i < n; i++) intime[i] = 0.0;
    
    for(j = 0; j < m; j++)
    {
        
        f = freqs[j];
        
        
         sqf = sqrt(f);
        
        for(i = 0; i < n; i++)
        {
            intime[i] += sqf*tfR[j][i];
        }
        
    }
    
    x = 0.0;
    j = 0;
    // out = fopen("testtime.dat","w");
    for(i=0; i< n; i++)
    {
        // fprintf(out,"%e %e %e\n",times[i], intime[i], ref[i]);
        
        if(fabs(ref[i]) > 0.01)
        {
            j++;
            x += intime[i]/ref[i];
        }
    }
    //fclose(out);
    
    x /= sqrt((double)(2*n));
    
    scale = (double)j/x;
    
    // printf("scaling = %e %e\n", x/(double)j, (double)j/x);
    
    free_double_vector(data);
    free_double_vector(ref);
    free_double_vector(intime);
    free_double_matrix(tf,m);
    free_double_matrix(tfR,m);
    free_double_matrix(tfI,m);
    
    return scale;
    
}


void recursive_phase_evolution(double dre, double dim, double *cosPhase, double *sinPhase)
{
    /* Update re and im for the next iteration. */
    double cosphi = *cosPhase;
    double sinphi = *sinPhase;
    double x, y;
    
    x = (cosphi*dre + sinphi*dim);
    y = (sinphi*dre - cosphi*dim);
    
    double newRe = cosphi - x;
    double newIm = sinphi - y;
    
    *cosPhase = newRe;
    *sinPhase = newIm;
    
}

void SineGaussianF(double *hs, double *sigpar, double Tobs, int N)
{
    double f0, t0, Q, sf, sx, Amp;
    double fmx, fmn;//, fac;
    double phi, f;
    double tau;
    double re,im;
    double q, p, u;
    double A, B, C;
    double Am, Bm, Cm;
    double a, b, c;
    double cosPhase_m, sinPhase_m, cosPhase_p, sinPhase_p;
    
    int i, imid, istart,istop,imin,imax,iend,even,odd;
    
    t0  = sigpar[0];
    f0  = sigpar[1];
    Q   = sigpar[2];
    Amp = sigpar[3];
    phi = sigpar[4];
    
    tau = Q/(TPI*f0);
    
    fmx = f0 + 3.0/tau;  // no point evaluating waveform past this time (many efolds down)
    fmn = f0 - 3.0/tau;  // no point evaluating waveform before this time (many efolds down)
    
    i = (int)(f0*Tobs);
    imin = (int)(fmn*Tobs);
    imax = (int)(fmx*Tobs);
    if(imax - imin < 10)
    {
        imin = i-5;
        imax = i+5;
    }
    
    if(imin < 0) imin = 1;
    if(imax > N/2) imax = N/2;
    
    hs[0] = 0.0;
    hs[N/2] = 0.0;
    
    for(i = 1; i < N/2; i++)
    {
        hs[i] = 0.0;
        hs[N-i] = 0.0;
    }
    
    /* Use recursion relationship  */
    
    //incremental values of exp(iPhase)
    double dim = sin(TPI*t0/Tobs);
    double dre = sin(0.5*(TPI*t0/Tobs));
    dre = 2.0*dre*dre;
    
    double amplitude = 0.5*(Amp)*RTPI*tau;
    double pi2tau2   = PI*PI*tau*tau;
    double Q2        = Q*Q/f0;
    
    
    imid = (int)(f0*Tobs);
    
    
    q = Q*Q/(f0*Tobs);
    p = PI*PI*tau*tau/(Tobs*Tobs);
    u = PI*PI*tau*tau/Tobs*(2.0*f0*Tobs-1.0);
    
    Am = exp(-q*(double)(imid));
    Bm = exp(-p*(((double)(imid)-f0*Tobs)*((double)(imid)-f0*Tobs)));
    Cm = 1.0;
    
    a = exp(-q);
    b = exp(-p*(1.0+2.0*((double)(imid)-f0*Tobs)));
    c = exp(-2.0*p);
    
    // sine and cosine of phase at reference frequency
    f = (double)(imid)/Tobs;
    double phase = TPI*f*t0;
    double cosPhase_m0  = cos(phase-phi);
    double sinPhase_m0  = sin(phase-phi);
    double cosPhase_p0  = cos(phase+phi);
    double sinPhase_p0  = sin(phase+phi);
    
    // start in the middle and work outwards
    
    A = Am;
    B = Bm;
    C = Cm;
    
    cosPhase_m  = cosPhase_m0;
    sinPhase_m  = sinPhase_m0;
    cosPhase_p  = cosPhase_p0;
    sinPhase_p  = sinPhase_p0;
    
    for(i = imid; i < imax; i++)
    {
        even = 2*i;
        odd = even + 1;
        f = (double)(i)/Tobs;
        
        sf = amplitude*B;
        sx = A;
        re = sf*(cosPhase_m+sx*cosPhase_p);
        im = -sf*(sinPhase_m+sx*sinPhase_p);
        
        //  printf("%e %e\n", exp(-Q2*f)-A, exp(-pi2tau2*(f-f0)*(f-f0))-B);
        
        hs[i] = re;
        hs[N-i] = im;
        
        A *= a;
        B *= (C*b);
        C *= c;
        
        /* Now update re and im for the next iteration. */
        recursive_phase_evolution(dre, dim, &cosPhase_m, &sinPhase_m);
        recursive_phase_evolution(dre, dim, &cosPhase_p, &sinPhase_p);
    }
    
    // reset to midpoint
    
    A = Am;
    B = Bm;
    C = Cm;
    
    cosPhase_m  = cosPhase_m0;
    sinPhase_m  = sinPhase_m0;
    cosPhase_p  = cosPhase_p0;
    sinPhase_p  = sinPhase_p0;
    
    a = 1.0/a;
    b = exp(p*(-1.0+2.0*((double)(imid)-f0*Tobs)));
    // c unchanged
    
    // interate backwards in phase
    dim *= -1.0;
    
    for(i = imid; i > imin; i--)
    {
        even = 2*i;
        odd = even + 1;
        f = (double)(i)/Tobs;
        
        sf = amplitude*B;
        sx = A;
        re = sf*(cosPhase_m+sx*cosPhase_p);
        im = -sf*(sinPhase_m+sx*sinPhase_p);
        
        // printf("%e %e\n", exp(-Q2*f)-A, exp(-pi2tau2*(f-f0)*(f-f0))-B);
        
        hs[i] = re;
        hs[N-i] = im;
        
        A *= a;
        B *= (C*b);
        C *= c;
        
        /* Now update re and im for the next iteration. */
        recursive_phase_evolution(dre, dim, &cosPhase_m, &sinPhase_m);
        recursive_phase_evolution(dre, dim, &cosPhase_p, &sinPhase_p);
    }
    
    
}



void shift(double *a, double R, double delt, double pshift, int n, double Tobs)
{
    int i, j, k;
    double  f;
    double ReA, ImA;
    double cx, sx;
    
    for(i=0; i<n/2; i++)
    {
        f = (double)(i)/Tobs;
        j = i; // real
        k = n-i; // imaginary
        cx = cos(-pshift+TPI*delt*f);
        sx = sin(-pshift+TPI*delt*f);
        ReA = a[j]*cx-a[k]*sx;
        ImA = a[k]*cx+a[j]*sx;
        a[j] = R*ReA;
        a[k] = R*ImA;
    }
    
    return;
    
}


double fourier_nwip_shift(double *a, double *b, double delt, double pshift, int n, double Tobs, int imin, int imax)
{
    int i, j, k;
    double arg, product, f;
    double ReA, ReB, ImA, ImB;
    double cx, sx;
    
    // Does f_nwip with a given time and phase shift
    
    arg = 0.0;
    for(i=1; i<n/2; i++)
    {
        if(i > imin && i < imax)
        {
            f = (double)(i)/Tobs;
            j = i; // real
            k = n-i; // imaginary
            ReA = a[j];
            ImA = a[k];
            cx = cos(-pshift+TPI*delt*f);
            sx = sin(-pshift+TPI*delt*f);
            ReB = b[j]*cx-b[k]*sx;
            ImB = b[k]*cx+b[j]*sx;
            product = ReA*ReB + ImA*ImB;
            arg += product;
        }
    }
    
    return(arg);
    
}






int *int_vector(int N)
{
    return malloc( (N+1) * sizeof(int) );
}

void free_int_vector(int *v)
{
    free(v);
}

int **int_matrix(int N, int M)
{
    int i;
    int **m = malloc( (N+1) * sizeof(int *));
    
    for(i=0; i<N+1; i++)
    {
        m[i] = malloc( (M+1) * sizeof(int));
    }
    
    return m;
}

void free_int_matrix(int **m, int N)
{
    int i;
    for(i=0; i<N+1; i++) free_int_vector(m[i]);
    free(m);
}

double *double_vector(int N)
{
    return malloc( (N+1) * sizeof(double) );
}

void free_double_vector(double *v)
{
    free(v);
}

double **double_matrix(int N, int M)
{
    int i;
    double **m = malloc( (N+1) * sizeof(double *));
    
    for(i=0; i<N+1; i++)
    {
        m[i] = malloc( (M+1) * sizeof(double));
    }
    
    return m;
}

void free_double_matrix(double **m, int N)
{
    int i;
    for(i=0; i<N+1; i++) free_double_vector(m[i]);
    free(m);
}

double ***double_tensor(int N, int M, int L)
{
    int i,j;
    
    double ***t = malloc( (N+1) * sizeof(double **));
    for(i=0; i<N+1; i++)
    {
        t[i] = malloc( (M+1) * sizeof(double *));
        for(j=0; j<M+1; j++)
        {
            t[i][j] = malloc( (L+1) * sizeof(double));
        }
    }
    
    return t;
}

void free_double_tensor(double ***t, int N, int M)
{
    int i;
    
    for(i=0; i<N+1; i++) free_double_matrix(t[i],M);
    
    free(t);
}







