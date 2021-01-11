#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_sort_double.h>


// gcc -w -o waveband waveband.c -lgsl -lm

int main(int argc, char *argv[])
{
    
    int i, j, id, Nifo, NS, M;
    int labels[4];
    double Tobs, ttrig;
    double x, y;
    double *times, *wave, *samp;
    char command[1024];
    
    FILE *in;
    FILE *out;
    
    if(argc<6)
    {
        printf("./waveband Tobs trig_time #samples #waveforms detector-list\n");
        return 1;
    }
    
    Nifo = argc-5;
    Tobs = atof(argv[1]);
    ttrig = atof(argv[2]);
    NS = atoi(argv[3]);
    M = atoi(argv[4]);
    
    printf("%d %f %f %d %d\n", Nifo, Tobs, ttrig, NS, M);
    
    times = (double*)malloc(sizeof(double)* (NS));
    wave = (double*)malloc(sizeof(double)* (NS*M));
    samp = (double*)malloc(sizeof(double)* (M));
    
    for (i = 0; i < Nifo; ++i) labels[i] = atoi(argv[5+i]);
    
  for (id = 0; id < Nifo; ++id)
  {
      
    for (j = 0; j < M; ++j)
      {
       sprintf(command, "waves/wavewhite_%d_%d_%d_%d.dat", j+1, (int)(Tobs), (int)ttrig, labels[id]);
       in = fopen(command,"r");
       for (i = 0; i < NS; ++i)
        {
        fscanf(in,"%lf%lf%lf", &times[i], &wave[j*NS+i], &y);
       }
       fclose(in);
      }
      
      sprintf(command, "waves/wavestat_%d_%d_%d.dat", (int)(Tobs), (int)ttrig, labels[id]);
      out = fopen(command,"w");
      
      for (i = 0; i < NS; ++i)
      {
          for (j = 0; j < M; ++j)
          {
              samp[j] = wave[j*NS+i];
          }
          gsl_sort(samp, 1, M);
          
          fprintf(out, "%f %f %f %f\n", times[i], samp[M/10], samp[M/2], samp[9*M/10]);
      }
      fclose(out);
      
  }
    
    for (id = 0; id < Nifo; ++id)
    {
        
        for (j = 0; j < M; ++j)
        {
            sprintf(command, "waves/wavewhite_%d_%d_%d_%d.dat", j+1, (int)(Tobs), (int)ttrig, labels[id]);
            in = fopen(command,"r");
            for (i = 0; i < NS; ++i)
            {
                fscanf(in,"%lf%lf%lf", &times[i], &y, &wave[j*NS+i]);
            }
            fclose(in);
        }
        
        sprintf(command, "waves/waveenv_%d_%d_%d.dat", (int)(Tobs), (int)ttrig, labels[id]);
        out = fopen(command,"w");
        
        for (i = 0; i < NS; ++i)
        {
            for (j = 0; j < M; ++j)
            {
                samp[j] = wave[j*NS+i];
            }
            gsl_sort(samp, 1, M);
            
            fprintf(out, "%f %f %f %f\n", times[i], samp[M/10], samp[M/2], samp[9*M/10]);
        }
        fclose(out);
        
    }



}
