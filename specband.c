#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_sort_double.h>
#include <gsl/gsl_cdf.h>

#define TPI 6.2831853071795862319959269370884     // 2 Pi
#define SQPI 2.5066282746310002  // sqrt(TPI)

// gcc -w -o specband specband.c -lgsl -lm

void ADtest(int k, int Ntest, double *samples, double *AD);

int main(int argc, char *argv[])
{
    
    int i, j, id, Nifo, NS, M;
    int labels[4];
    double Tobs, ttrig;
    double x, y, u, v;
    double *freqs, *spec, *samp, *data, *wdata, *AD, *ng, *ngfs;
    char command[1024];
    
    FILE *in;
    FILE *out;
    FILE *white;
    
    M = 500;
    
    // how many samples in the data?
    NS=0;
    
    //get data file
    in = fopen("specs/spec_sample_0.dat","r");
    
    while(!feof(in))
    {
        fscanf(in,"%lg%lg", &x, &y);
        NS++;
    }
    NS--;
    
    fclose(in);
    
    
    printf("%d %d\n", NS, M);
    
    
    
    freqs = (double*)malloc(sizeof(double)* (NS));
    spec = (double*)malloc(sizeof(double)* (NS*M));
    samp = (double*)malloc(sizeof(double)* (M));
    data = (double*)malloc(sizeof(double)* (2*NS));
    wdata = (double*)malloc(sizeof(double)* (2*NS));
    ng = (double*)malloc(sizeof(double)* (M));
    
    in = fopen("specs/data.dat","r");
    for (i = 0; i < NS; ++i)
    {
        fscanf(in,"%lf%lf%lf", &freqs[i], &data[2*i], &data[2*i+1]);
    }
    fclose(in);
    
    double DF, df;
    
    DF = 8.0; // width of window for Gaussianity test
    df = freqs[1]-freqs[0];
    
    int Ntest, k;
    
    Ntest = (int)(2.0*DF/df);
    k = 2*NS/Ntest;
    
    printf("%d %d %f\n", Ntest, k, (double)(k)*0.05);
    
    AD = (double*)malloc(sizeof(double)* (k));
    ngfs = (double*)malloc(sizeof(double)* (k));
    
    for (j = 0; j < M; ++j) ng[j] = 0.0;
    for (j = 0; j < k; ++j) ngfs[j] = 0.0;
    
    out = fopen("ngcnt.dat","w");
    for (j = 0; j < M; ++j)
      {
       sprintf(command, "specs/spec_sample_%d.dat", j);
       in = fopen(command,"r");
       for (i = 0; i < NS; ++i)
        {
        fscanf(in,"%lf%lf", &freqs[i], &spec[j*NS+i]);
        x = sqrt(spec[j*NS+i]);
        wdata[2*i] = data[2*i]/x;
        wdata[2*i+1] = data[2*i+1]/x;
        //if(i== 10) printf("%e %e\n", wdata[2*i], wdata[2*i+1]);
       }
       fclose(in);
       ADtest(k, Ntest, wdata, AD);
        for (i = 0; i < k; ++i)
        {
            //printf("%f\n", AD[i]);
            if(AD[i] > 0.752)
            {
                ng[j] += 1.0;
                ngfs[i] += 1.0;
            }
                
        }
          fprintf(out, "%d %f\n", j, ng[j]);
      }
     fclose(out);
    
    for (i = 0; i < k; ++i) ngfs[i] /= (double)(M);
    
    out = fopen("ngf.dat","w");
    for (i = 0; i < k; ++i) fprintf(out,"%f %e\n", freqs[0]+((double)(i)+0.5)*DF, ngfs[i]);
    fclose(out);
      
      
      out = fopen("specstat.dat","w");
      white = fopen("medwhite.dat","w");
      for (i = 0; i < NS; ++i)
      {
          for (j = 0; j < M; ++j)
          {
              samp[j] = spec[j*NS+i];
          }
          gsl_sort(samp, 1, M);
          
          fprintf(out, "%e %e %e %e\n", freqs[i], samp[M/20], samp[M/2], samp[19*M/20]);
          
          x = sqrt(samp[M/2]);
          fprintf(white, "%e %e %e\n", freqs[i], data[2*i]/x, data[2*i+1]/x);
      }
      fclose(out);
      fclose(white);
    
     
    
      


}

void ADtest(int k, int Ntest, double *samples, double *AD)
{
    
    double *samps;
    int j, n;
    double mean, var, std, S, A, x, y, z, u, lx, ly;
    
    samps = (double*)malloc(sizeof(double)* (Ntest));
    
    for(j=0; j<k; j++)
        {
            
        for(n=0; n<Ntest; n++)
          {
              samps[n] = samples[n+j*Ntest];
          }
            
            mean = 0.0;
            var = 0.0;
            
            for(n=0; n<Ntest; n++)
               {
                   mean += samps[n];
                   var += samps[n]*samps[n];
               }
                   
               mean /= (double)(Ntest);
               var /= (double)(Ntest);
               var -= mean*mean;
               std = sqrt(var);
            
           // printf("%f %f\n", mean, std);
            
            for(n=0; n<Ntest; n++)
            {
                samps[n] -= mean;
                samps[n] /= std;
                //printf("%e\n", samps[n]);
            }
    
    
        gsl_sort(samps,1,Ntest);
    
       
        // Anderson-Darling statistics
        S=0.0;
        for(n=0; n<Ntest; n++)
        {
     
            x = gsl_cdf_ugaussian_P(samps[n]);
            
             //printf("%e\n", x);
            
            if(x < 0.999)
            {
                lx = log(x);
                ly = log(1.0-x);
            }
            else
            {
                u = samps[n];
                z = 1.0-1.0/(u*u)+3.0/(u*u*u*u)-15.0/(u*u*u*u*u*u)+105.0/(u*u*u*u*u*u*u*u);
                //printf("%e\n", z/u);
                y = z*exp(-u*u/2.0)/(u*SQPI);
                lx = -y;
                ly = -u*u/2.0 +log(z/(u*SQPI));
            }
            
            
             S += ((2.*(double)(n+1)-1.0)*lx + ((double)(2*Ntest)-2.*(double)(n+1)+1.0)*ly)/(double)(Ntest);
        }
        
        
        
        A = (double)(-1*Ntest) - S;
        
        A *= (1.0+0.75/(double)(Ntest)+2.25/(double)(Ntest*Ntest));
            
          //  printf("%e\n", A);
            
        AD[j] = A;
            
        }
            
        free(samps);
            
        return;
}
