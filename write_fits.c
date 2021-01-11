#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "chealpix.h"
#include "fitsio.h"

// used by memory allocation
#define NR_END 1
#define FREE_ARG char*


// compile gcc -o write_fits write_fits.c -lm -lcfitsio
// usage:  ./write_fits sky.dat
// output will be called sky.fits

void setCoordSysHP(char *,char *);
float *vector(long nl, long nh);
void free_vector(float *v, long nl, long nh);
void nrerror(char error_text[]);

int main(int argc, char *argv[])
{
	int i, j, ND, nside;
	double x;
	float *data;
    FILE *in;	
	
	if (argc!=2) nrerror("usage: ./write_fits sky.dat");
	
 	ND = 0;
	
	in = fopen(argv[1],"r");
	
	while(!feof(in))
	{
		fscanf(in,"%lf", &x);
		
		ND++;
	}
	rewind(in);
	
	ND -= 1;
	
	nside = (int)(sqrt((double)(ND)/12.0));
	
	printf("NSIDE = %d\n", nside);
	
	data = vector(0,ND-1);
	
	// Read in  Data
    for(j=0; j< ND; j++)
    {
		fscanf(in,"%f", &data[j]);
    }
    fclose(in);

	write_healpix_map(data, nside, "sky.fits", 0, "C");

}



void nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
	int rx;
	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}

void free_vector(float *v, long nl, long nh)
/* free a double vector allocated with dvector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}


float *vector(long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
	float *v;
	
	v=(float *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(float)));
	if (!v) nrerror("allocation failure in dvector()");
	return v-nl+NR_END;
}



void write_healpix_map( float *signal, long nside, char *filename, char nest, char *coordsys) {

  /*******************************************************************/
  /* Create a binary table extension                                 */
  /*******************************************************************/
  fitsfile *fptr;       /* pointer to the FITS file, defined in fitsio.h */
  int status, hdutype;
  long firstrow, firstelem;
  
  int bitpix   =  SHORT_IMG;
  long naxis   =   0;
  long naxes[] = {0,0};
  
  int tfields   = 1;
  long nrows;
  
  char order[9];                 /* HEALPix ordering */
  char extname[] = "BINTABLE";   /* extension name */
  char *ttype[] = { "SIGNAL" };
  char *tform[] = { "1E" };
  char *tunit[] = { " " };
  char coordsys9[9];
  
  
  /* Calculate the number of pixels in the full map */
  nrows = 12L*nside*nside;
  
  /* initialize status before calling fitsio routines */
  status = 0;
  
  /* create new FITS file */
  if (fits_create_file(&fptr, filename, &status)) 
    fprintf(stderr, "%s (%d): Could not create new fits file.\n", 
	    __FILE__, __LINE__);
  
  if ( fits_create_img(fptr,  bitpix, naxis, naxes, &status) )
    fprintf(stderr, "%s (%d): Could not create new image file.\n", 
	    __FILE__, __LINE__);
  
  if ( fits_write_date(fptr, &status) )
    fprintf(stderr, "%s (%d): Could not add date.\n", 
	    __FILE__, __LINE__);
  
  /* move to 1nd HDU  */
  if ( fits_movabs_hdu(fptr, 1, &hdutype, &status) ) 
    fprintf(stderr, "%s (%d): Could not move to first HDU.\n", 
	    __FILE__, __LINE__);
  
  /* append a new empty binary table onto the FITS file */
  if ( fits_create_tbl( fptr, BINARY_TBL, nrows, tfields, ttype, tform,
			tunit, extname, &status) )
    fprintf(stderr, "%s (%d): Could not create new binary table.\n", 
	    __FILE__, __LINE__);
  
  if (fits_write_key(fptr, TSTRING, "PIXTYPE", "HEALPIX",
		     "HEALPIX Pixelisation", &status))
    fprintf(stderr, "%s (%d): Could not write PIXTYPE keyword.\n", 
	    __FILE__, __LINE__);
  
  if (nest) strcpy(order, "NESTED  ");
  else      strcpy(order, "RING    ");
  if (fits_write_key(fptr, TSTRING, "ORDERING", order, 
		     "Pixel ordering scheme, either RING or NESTED", &status))
    fprintf(stderr, "%s (%d): Could not write ORDERING keyword.\n", 
	    __FILE__, __LINE__);
  
  if (fits_write_key(fptr, TLONG, "NSIDE", &nside,
		     "Resolution parameter for HEALPIX", &status))  
    fprintf(stderr, "%s (%d): Could not write NSIDE keyword.\n", 
	    __FILE__, __LINE__);
 
  /*****NEW*************/
  setCoordSysHP(coordsys,coordsys9);
  if (fits_write_key(fptr, TSTRING, "COORDSYS", coordsys9,"Pixelisation coordinate system", &status))
    fprintf(stderr, "%s (%d): Could not write COORDSYS keyword.\n",__FILE__, __LINE__);

  /*****END**NEW*******/

  if (fits_write_comment(fptr,"           G = Galactic, E = ecliptic, C = celestial = equatorial  ", &status))
    fprintf(stderr, "%s (%d): Could not write COORDSYS explanation keyword.\n", 
	    __FILE__, __LINE__);
  
  firstrow  = 1;  /* first row in table to write   */
  firstelem = 1;  /* first element in row  (ignored in ASCII tables)  */
  
  if (fits_write_col(fptr, TFLOAT, 1, firstrow, firstelem, nrows, signal,
		     &status))
    fprintf(stderr, "%s (%d): Could not write signal.\n", __FILE__, __LINE__);

  /*fits_write_col(fptr, TLONG, 2, firstrow, firstelem, nrows, pixel,  &status);*/
  /*fits_write_col(fptr, TLONG, 3, firstrow, firstelem, nrows, n_obs,  &status);*/
  /*fits_write_col(fptr, TFLOAT, 4, firstrow, firstelem, nrows, serror,&status);*/
  
  if ( fits_close_file(fptr, &status) )       /* close the FITS file */
    fprintf(stderr, "%s (%d): Could not close file.\n", 
	    __FILE__, __LINE__);
  
  return;
}


void setCoordSysHP(char *coordsys,char *coordsys9){
  
  strcpy(coordsys9,"C       ");
  
  if(strncmp(coordsys,"G",1)!=0 &&  strncmp(coordsys,"E",1)!=0 &&  strncmp(coordsys,"C",1)!=0 && strncmp(coordsys,"Q",1)!=0)
    fprintf(stderr, "%s (%d): System Cordinates is not correct (Galactic,Ecliptic,Celestial=Equatorial). Celestial system was set.\n", __FILE__, __LINE__);
  
  
/*    if(strncmp(coordsys,"GAL",3)==0) */
/*      strcpy(coordsys9,"G       "); */
/*    else if(strncmp(coordsys,"ECL",3)==0) */
/*      strcpy(coordsys9,"E       "); */
/*    else if(strncmp(coordsys,"EQU",3)==0 || strncmp(coordsys,"CEL",3)==0) */
/*      strcpy(coordsys9,"C       "); */

  if(strncmp(coordsys,"G",1)==0)
    strcpy(coordsys9,"G       ");
  else if(strncmp(coordsys,"E",1)==0)
    strcpy(coordsys9,"E       ");
  else if(strncmp(coordsys,"C",1)==0 || strncmp(coordsys,"Q",1)==0)
    strcpy(coordsys9,"C       ");

  

}
