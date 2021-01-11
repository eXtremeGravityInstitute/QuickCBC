#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>


#define NR_END 1
#define FREE_ARG char*

/* gcc -o gwoscdump gwoscdump.c -lm */


double *dvector(long nl, long nh);
void free_dvector(double *v, long nl, long nh);


int main(int argc, char *argv[])
{
  int i, j, k, M, N, Nf, Nstep, Nclean, ii, m, rs, tsi, tti;
    int jj, kk;
  int linlog;
  int imin, imax;
  double SNR;
  double junk, Tobs, fix, f, t, t0, dt, dtm, df, x, y, dx;
  double fmax, fmin, dfx, Q, fny, delt, scale, dlnf;
  double Hmax, Lmax;
  double pshift;
  double *freqs, *data, *ref;
  double *inp, *oup, *slice;
  double *tempH1;
  double *tempL1;
  double *H1dat, *L1dat;
  double *Hclean, *Lclean;
  double *Hraw, *Lraw, *traw;
  double *H, *L, *times;
  double *SH, *SL;
  double **geo;
  double *specH1, *specL1, *sspecH1, *sspecL1;
  double *sdata;
  double **tfmap, **tfavH1, **tfavL1;
  double *intime, *sqf;
  double sthresh, nthresh;
  double sigmean, sigmedian;
  int subscale, octaves;
    int mmax, flag;
    double SNRsq, pH, pL, pmax;
    double SNRH, SNRL, pw, alpha;
    int specprint, printglitch;
   double t_rise, s1, s2, ascale, warm;
  double ttrig, tstart, tstart_clean, Tclean, starttime, endtime;
	
  char filename[1024];
  char command[1024];
    char L1name[1024];
    char H1name[1024];
    char line1[100];
    char line2[100];
    char line3[100];
    
    char channel[1024];

  int n;


  FILE *in;
  FILE *ifp;
  FILE *out;
    
    
    // trigger time for BNS GW170817 = 1187008882.0;
    // trigger time for GW150914 = 1126259462.0;
    // trigger time for GW170814 = 1186741861.0;
   
    
    if(argc!=5)
    {
        printf("./gwoscdump filename Tobs trig_time obs\n");
        return 1;
    }
    
    // Designed to use data from https://www.gw-openscience.org/catalog/GWTC-1-confident/html/
    // The code is happy reading in the 4 kHz or 16 kHz data in .txt format. 32 seconds is usually enough
    
    in = fopen(argv[1],"r");
    if (in == NULL){
        printf("Could not open file %s", argv[1]);
        return 2;
    }
    
    
    
    Tobs = atof(argv[2]);
    ttrig = atof(argv[3]);
    m = atoi(argv[4]);
    tstart = ttrig + 2.0 - Tobs ;  // puts the trigger 2 seconds from the end
    //tstart = ttrig - 0.75*Tobs;  // puts the trigger 3/4 of the way into the data snippet

    starttime = tstart;
    endtime = tstart + Tobs;
    
    // Strip off text header
    
    fgets(line1, 100, in);
    fgets(line2, 100, in);
    fgets(line3, 100, in);
    
    /* printf("%s\n", line1);
    printf("%s\n", line2);
    printf("%s\n", line3); */
    
    long num[4];
    
    
    char *ptr = line2;
    
    i = 0;
    while (*ptr)
    {
        if  (isdigit (*ptr) ){
            long  val = strtol (ptr, &ptr, 10);
            num[i] = val;
            i++;
        }  else {
            ptr++;
        }
    }
    
    char *pt = line3;
    
    while (*pt)
    {
        if  (isdigit (*pt) ){
            long  val = strtol (pt, &pt, 10);
            num[i] = val;
            i++;
        }  else {
            pt++;
        }
    }
    
    printf("cadence %ld GPS start %ld duration %ld\n", num[0], num[1], num[2]);


    N = (int)num[0]*(int)num[2];
    t0 = (double)num[1];
    
    dt = (double)(num[2])/(double)(N);
    
    if(tstart < t0 || endtime > t0+(double)(num[2]))
    {
        printf("Requested data range outside of LOSC file\n");
    }
    
        sprintf(command, "frame_%d_%d_%d.dat", (int)(Tobs), (int)ttrig, m);
        

        out = fopen(command, "w");
        for(i=0; i< N; i++)
        {
            fscanf(in,"%lf", &x);
            t = t0 + (double)(i)*dt;
            if(t >= tstart && t < endtime)
            {
                fprintf(out,"%.16e %.16e\n",  t, x);
            }
        }
        fclose(out);
        fclose(in);
    
    
    return 0;

}

double *dvector(long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
    double *v=0;
    
    v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
    if (!v) fprintf(stderr,"allocation failure in dvector()");
    return v-nl+NR_END;
}

void free_dvector(double *v, long nl, long nh)
/* free a double vector allocated with dvector() */
{
    free((FREE_ARG) (v+nl-NR_END));
}





