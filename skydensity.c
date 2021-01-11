/* gcc -o skydensity skydensity.c -lgsl */

/* Standard Includes */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <gsl/gsl_sort_double.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_vector.h>

#define PIn 3.1415926535897931159979634685442

void ang2pix_ring( const long nside, double theta, double phi, long *ipix);
void pix2ang_ring( long nside, long ipix, double *theta, double *phi);

double *double_vector(int N);
void free_double_vector(double *v);
int *int_vector(int N);
void free_int_vector(int *v);

int main(int argc, char* argv[])
{

  long i, j, k, Nside, Npix, Nsample;
  double ttrig;
  double theta, phi;
  double pattern, GMT;
  double costh, x, y;
  double *map, *maps;
  long *indx;
    
  FILE *out;
    FILE *in;
    
 if(argc != 3) { printf("./skydensity file Nside\n"); return 0;}
    

    in = fopen(argv[1],"r");
    
    Nsample = -1;
    while (!feof(in))
    {
        fscanf(in,"%ld%lf%lf%lf%lf", &k, &x, &x, &x, &x);
        Nsample++;
    }
    rewind(in);
    
    
    // Nside must be a power of 2
    Nside = atoi(argv[2]);

    Npix = 12*Nside*Nside;
    
    map = double_vector(Npix);

    
    for(i=0; i< Npix; i++) map[i] = 0.0;
    
    out = fopen("skycheck.dat","w");
    for(i=0; i< Nsample; i++)
    {
        fscanf(in,"%ld%lf%lf%lf%lf", &k, &x, &phi, &costh, &x);
        theta = acos(costh);
        if(phi < 0.0) phi += 2.0*PIn;
        if(phi > 2.0*PIn) phi -= 2.0*PIn;
        fprintf(out,"%ld %f %f\n", i, phi, theta);
        ang2pix_ring(Nside, theta, phi, &k);
        //printf("%ld %ld %f %f\n", i, k, theta, phi);
        map[k] += 1.0;
    }
    fclose(in);
    fclose(out);
    
    for(i=0; i< Npix; i++) map[i] /= (double)(Nsample);
    
    gsl_vector *mv = gsl_vector_alloc (Npix);
    
    for(i=0; i< Npix; i++) gsl_vector_set(mv, i, map[i]);
    
    gsl_permutation *perm = gsl_permutation_alloc(Npix);
    
    gsl_sort_vector_index(perm, mv);
    
    k = 0;
    x = 0.0;
    do
    {
        k++;
        x += gsl_vector_get(mv,perm->data[Npix-k]);
    }while(x < 0.5);
    
    y = (double)(k)/(double)(Npix)*4.0*PIn*(180.0/PIn)*(180.0/PIn);
    printf("50 percent credible interval = %f  square degrees\n", y);
    
    do
    {
    k++;
    x += gsl_vector_get(mv,perm->data[Npix-k]);
    }while(x < 0.9);
    
    y = (double)(k)/(double)(Npix)*4.0*PIn*(180.0/PIn)*(180.0/PIn);
    printf("90 percent credible interval = %f  square degrees\n", y);
    
    pix2ang_ring(Nside, perm->data[Npix-1], &theta, &phi);
    
    printf("MAP location theta = %f phi = %f\n", theta*180.0/PIn, phi*180.0/PIn);
    
    gsl_permutation_free(perm);
    gsl_vector_free(mv);
    
  out = fopen("sky.dat", "w");
  for(i=0; i< Npix; i++) fprintf(out,"%f\n", map[i]);
  fclose(out);
    
  free_double_vector(map);

  return 1;

}


double *double_vector(int N)
{
    return malloc( (N+1) * sizeof(double) );
}

void free_double_vector(double *v)
{
    free(v);
}

int *int_vector(int N)
{
    return malloc( (N+1) * sizeof(int) );
}

void free_int_vector(int *v)
{
    free(v);
}



void ang2pix_ring( const long nside, double theta, double phi, long *ipix) {
  /*
    c=======================================================================
    c     gives the pixel number ipix (RING) 
    c     corresponding to angles theta and phi
    c=======================================================================
  */
  
  int nl2, nl4, ncap, npix, jp, jm, ipix1;
  double  z, za, tt, tp, tmp;
  int ir, ip, kshift;
  
  double piover2 = 0.5*M_PI;
  double PI=M_PI;
  double twopi=2.0*M_PI;
  double z0=2.0/3.0;
  long ns_max=8192;
  
  if( nside<1 || nside>ns_max ) {
    fprintf(stderr, "%s (%d): nside out of range: %ld\n", __FILE__, __LINE__, nside);
    exit(0);
  }
  
  if( theta<0. || theta>PI) {
    fprintf(stderr, "%s (%d): theta out of range: %f\n", __FILE__, __LINE__, theta);
    exit(0);
  }
  
  z = cos(theta);
  za = fabs(z);
  if( phi >= twopi)  phi = phi - twopi;
  if (phi < 0.)     phi = phi + twopi;
  tt = phi / piover2;//  ! in [0,4)
  
  nl2 = 2*nside;
  nl4 = 4*nside;
  ncap  = nl2*(nside-1);// ! number of pixels in the north polar cap
  npix  = 12*nside*nside;
  
  if( za <= z0 ) {
    
    jp = (int)floor(nside*(0.5 + tt - z*0.75)); /*index of ascending edge line*/
    jm = (int)floor(nside*(0.5 + tt + z*0.75)); /*index of descending edge line*/
    
    ir = nside + 1 + jp - jm;// ! in {1,2n+1} (ring number counted from z=2/3)
    kshift = 0;
    if (fmod(ir,2)==0.) kshift = 1;// ! kshift=1 if ir even, 0 otherwise
    
    ip = (int)floor( ( jp+jm - nside + kshift + 1 ) / 2 ) + 1;// ! in {1,4n}
    if( ip>nl4 ) ip = ip - nl4;
    
    ipix1 = ncap + nl4*(ir-1) + ip ;
  }
  else {
    
    tp = tt - floor(tt);//      !MOD(tt,1.d0)
    tmp = sqrt( 3.*(1. - za) );
    
    jp = (int)floor( nside * tp * tmp );// ! increasing edge line index
    jm = (int)floor( nside * (1. - tp) * tmp );// ! decreasing edge line index
    
    ir = jp + jm + 1;//        ! ring number counted from the closest pole
    ip = (int)floor( tt * ir ) + 1;// ! in {1,4*ir}
    if( ip>4*ir ) ip = ip - 4*ir;
    
    ipix1 = 2*ir*(ir-1) + ip;
    if( z<=0. ) {
      ipix1 = npix - 2*ir*(ir+1) + ip;
    }
  }
  *ipix = ipix1 - 1;// ! in {0, npix-1}
  
}



void pix2ang_ring( long nside, long ipix, double *theta, double *phi) {
  /*
    c=======================================================================
    c     gives theta and phi corresponding to pixel ipix (RING) 
    c     for a parameter nside
    c=======================================================================
  */
  
  int nl2, nl4, npix, ncap, iring, iphi, ip, ipix1;
  double  fact1, fact2, fodd, hip, fihip;
  double PI=M_PI;
  //      PARAMETER (pi     = 3.1415926535897932384626434d0)
  //      parameter (ns_max = 8192) ! 2^13 : largest nside available
  
  int ns_max=8192;
  
  if( nside<1 || nside>ns_max ) {
    fprintf(stderr, "%s (%d): nside out of range: %ld\n", __FILE__, __LINE__, nside);
    exit(0);
  }
  npix = 12*nside*nside;      // ! total number of points
  if( ipix<0 || ipix>npix-1 ) {
    fprintf(stderr, "%s (%d): ipix out of range: %ld\n", __FILE__, __LINE__, ipix);
    exit(0);
  }
  
  ipix1 = ipix + 1; // in {1, npix}
  nl2 = 2*nside;
  nl4 = 4*nside;
  ncap = 2*nside*(nside-1);// ! points in each polar cap, =0 for nside =1
  fact1 = 1.5*nside;
  fact2 = 3.0*nside*nside;
  
  if( ipix1 <= ncap ) {  //! North Polar cap -------------
    
    hip   = ipix1/2.;
    fihip = floor(hip);
    iring = (int)floor( sqrt( hip - sqrt(fihip) ) ) + 1;// ! counted from North pole
    iphi  = ipix1 - 2*iring*(iring - 1);
    
    *theta = acos( 1. - iring*iring / fact2 );
    *phi   = (1.*iphi - 0.5) * PI/(2.*iring);
  }
  else if( ipix1 <= nl2*(5*nside+1) ) {//then ! Equatorial region ------
    
    ip    = ipix1 - ncap - 1;
    iring = (int)floor( ip / nl4 ) + nside;// ! counted from North pole
    iphi  = (int)fmod(ip,nl4) + 1;
    
    fodd  = 0.5 * (1 + fmod((double)(iring+nside),2));//  ! 1 if iring+nside is odd, 1/2 otherwise
    *theta = acos( (nl2 - iring) / fact1 );
    *phi   = (1.*iphi - fodd) * PI /(2.*nside);
  }
  else {//! South Polar cap -----------------------------------
    
    ip    = npix - ipix1 + 1;
    hip   = ip/2.;
/* bug corrige floor instead of 1.* */
    fihip = floor(hip);
    iring = (int)floor( sqrt( hip - sqrt(fihip) ) ) + 1;//     ! counted from South pole
    iphi  = (int)(4.*iring + 1 - (ip - 2.*iring*(iring-1)));
    
    *theta = acos( -1. + iring*iring / fact2 );
    *phi   = (1.*iphi - 0.5) * PI/(2.*iring);
  }
}


