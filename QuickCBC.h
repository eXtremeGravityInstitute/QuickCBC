#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sort_double.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_eigen.h>
#include "IMRPhenomD.h"

struct Net
{
    int Nifo;
    double Tobs;
    double tmax;
    double tmin;
    double GMST;
    int *labels;
    double *tds;
    double **delays;
    int *lstart;
    int *lstop;
};

struct Het
{
    int M;
    int MM;
    int ell;
    double hh;
    double hr;
    double logLb;
    double **S0;
    double **S1;
    double **AS;
    double **SN;
    double **amp;
    double **phase;
    double **rc;
    double **rs;
    double ***lc;
    double ***ls;
    double *pref;
    RealVector *freq;
};


void detector_shifts(struct Net *net, double *params);
double gmst(double ttrig);
double leap(double tgps);

void pmap(struct Net *net, double *pallx, double *paramx, double *sky);
void pmap_back(struct Net *net, double *pall, double *param, double *sky);

void printwave(struct Net *net, int N, RealVector *freq, double **paramx, int *who, double Tobs, double ttrig, int iter);
void printwaveall(struct Net *net, int N, RealVector *freq, double *paramx, double **SN, double Tobs, double ttrig, int iter);

void wavemax(struct Net *net, int N, double **tp, RealVector *freq, double **paramx, int *who, double **SN, double Tobs);

void pairs(struct Net *net);
void time_delays(struct Net *net);

void pbt_shift_cut(double *corr, double *corrf, double *data1, double *data2, double *Sn, int imin, int imax, int N);

double fourier_nwip(double *a, double *b, double *Sn, int imin, int imax, int N);

double f_nwip(double *a, double *b, int n);
double fourier_nwip_shift(double *a, double *b, double delt, double pshift, int n, double Tobs, int imin, int imax);
double fb_nwip(double *a, double *b, int n, int imin, int imax);

void cossin(struct Net *net, double *hcos, double *hsin, RealVector *freq, double *params, int N);
void extrinsictemplates(struct Net *net, double **hwave, RealVector *freq, double *hcos, double *hsin, double *params, int N);

void templates(struct Net *net, double **hwave, RealVector *freq, double *params, int N);
void geotemplate(double *gwave, RealVector *freq, double *params, int N);
void fulltemplates(struct Net *net, double **hwave, RealVector *freq, double *params, int N);
void fullphaseamp(struct Net *net, double **amp, double **phase, RealVector *freq, double *params, int N);
void Fisher_All(struct Net *net, double **fish, double *params, RealVector *freq, double **SN, int N, double Tobs);
void freehet(struct Net *net, struct Het *het);
void heterodyne(struct Net *net, struct Het *het, double **D, double *params, RealVector *freq, double **SN, double **SM, int N, double Tobs);
double log_likelihood_het(struct Net *net, struct Het *het, double *params, double Tobs);
double log_likelihood_print(struct Net *net, double **D, double *params, RealVector *freq, double **SN, int N, double Tobs);
double log_likelihood(struct Net *net, double **D, double *params, RealVector *freq, double **SN, int N, double Tobs);
double log_likelihood_full(struct Net *net, double **D, double *params, RealVector *freq, double **SN, double *rho, int N, double Tobs);
double log_likelihood_ext(struct Net *net, double **D, double *params, RealVector *freq, double *hcos, double *hsin, double **SN, double *rho, int N, double Tobs);
double log_likelihood_max(struct Net *net, double **D, double *params, RealVector *freq, double **SN, int N, double Tobs, double tmin, double tmax, int pflag);
double log_likelihood_penalized(int ii, double *params, double *D, double *H, double *AU, double *PU, double *TDU, double *SN, int N, double Tobs, int pflag);
double log_likelihood_test(struct Net *net, struct Het *het, double **D, double *params, RealVector *freq, double **SN, double **SM, int N, double Tobs);

void jacobi(double **a, int n, double e[], double **v, int *nrot);

void fisherskysetup(struct Net *net, double **wave, double **HH, double Tobs, int n);

void ringfind(struct Net *net, double *tdelays, double *params, double *SNRsq, gsl_rng * r);
double skydensity(double *paramsx, double *paramsy, double GMST, int ifo1, int ifo2, int iref);
void skymap(double *paramsx, double *paramsy, double GMST, int ifo1, int ifo2, int iref);

void uvwz_all(double *u, double *v, double *w, double *z, double *params);
void exsolve_all(double *phiy, double *psiy, double *DLy, double *ciotay, double uy, double vy, double wy, double zy);
double skydensity_all(double *paramsx, double *paramsy, double GMST, int ifo1, int ifo2);
void skymap_all(double *paramsx, double *paramsy, double GMST, int ifo1, int ifo2);

void skyring(struct Net *net, double *params, double *sky, double *pall, double *SNRsq, RealVector *freq, double **D, double **SN, int N, double Tobs, gsl_rng * r);

void skystart(struct Net *net, double *params, double *sky, double *pall, double *SNRsq, RealVector *freq, double **D, double **SN, int N, double Tobs, gsl_rng * r);

double skylike(struct Net *net, double *params, double *D, double *H, double **DHc,  double **DHs, double dt, int nt, int flag);
void skylikesetup(struct Net *net, double **data,  double **wave, double *D, double *H, double **DHc,  double **DHs, double Tobs, int n, int bn, int nt, int imin, int imax);
void upsample(int n, double T, int *nt, int *bn);

void uvwz(double *u, double *v, double *w, double *z, double *params);
void exsolve(double *phiy, double *psiy, double *Ay, double *ey, double uy, double vy, double wy, double zy);
void uvwz_sol(double *uy, double *vy, double *wy, double *zy, double ux, double vx, double wx, double zx, \
              double fp1x, double fp1y, double fc1x, double fc1y, double fp2x, double fp2y, double fc2x, double fc2y);
void fisher_matrix_fastsky(struct Net *net, double *params, double **fisher, double **HH);
void fisher_skyproposal(gsl_rng * r, double **skyvecs, double *skyevals, double *jump);

static const gsl_rng_type *rngtype;
static const gsl_rng *rng;

void AmpPhase(double *Af, double *Pf, RealVector *freq, double *params, int N);

void Fisher_Full(struct Net *net, double **fish, int d, double *params, RealVector *freq, double **SN, int N, double Tobs);

void de_jump(double *paramsx, double *paramsy, double **history, int m, int d, gsl_rng *r);

void tukey(double *data, double alpha, int N);
void tukey_scale(double *s1, double *s2, double alpha, int N);

void log_likelihood_scan(struct Net *net, double **D, double *params, RealVector *freq, double **SN, double *Larray, double **Tarray, int N, double Tobs);

void PDwave(double *wavef, RealVector *freq, double tcmin, double *params, int N);

void PDtemplates(double *waveH, double *waveL, RealVector *freq, double *params, int N);

void FisherEvec(double **fish, double *eval, double **evec, int d);

void SNRvsf(struct Net *net, double **D, double *params, RealVector *freq, double **SN, int N, double Tobs);

double globeden(double ***global, double *max, double *min, double Tobs, double *params, int N);
double globe(double ***global, double *max, double *min, double Tobs, double *params, int N, gsl_rng *r);

void dshifts(struct Net *net, double *sky, double *params);

void CBC_start(struct Net *net, int *mxc, FILE *chainI, double **paramx, double **skyx, double **pallx, int *who, double *heat, double ***history, double ***global, RealVector *freq, double **D, double **Dtime, double **SN, int N, double Tobs, gsl_rng *r);

void CBC_update(struct Net *net, int *mxc, FILE *chainI, FILE *chainE, double **paramx, double **skyx, double **pallx, int *who, double *heat, double ***history, double ***global, RealVector *freq, double **D, double **SN, int N, double Tobs, gsl_rng *r);

void MCMC_intrinsic(struct Net *net, int lmax, int *mxc, int M, FILE *chain, double **paramx, double **skyx, double **pallx, int *who, double *heat, double ***history, double ***global, RealVector *freq, double **D, double **SN, int N, double Tobs, gsl_rng *r);

void MCMC_all(struct Net *net, int *mxc, int M, FILE *chain, double **paramx, int *who, double *heat, double ***history, double ***global, RealVector *freq, double **D, double **SN, double **SM, int N, double Tobs, double ttrig, gsl_rng *r);

void skymcmc(struct Net *net, int MCX, int *mxc, FILE *chain, double **paramx, double **skyx, double **pallx, int *who, double *heat, double dtx, int nt, double *DD, double **WW, double ***DHc,  double ***DHs, double ***HH, double Tobs, gsl_rng * r);

void updatei(int k, struct Net *net, int lmax, double *logLx, double **paramx, double **paramy, double *min, double *max, double *Fscale, int *who, double *heat, double ***history, double ***global, RealVector *freq, double **D, double **SN, double **ejump, double ***evec, int N, double Tobs, int **cv, int **av, gsl_rng *r);

void update(int k, struct Net *net, struct Het *het, double *logLx, double *logPx, double **rhox, double **paramx, double **paramy, double *min, double *max, int *who, double *heat, double ***history, double ***global, RealVector *freq, double **D, double **SN, double **SM, double **ejump, double ***evec, int N, double Tobs, int **cv, int **av, gsl_rng *r);

int *int_vector(int N);
void free_int_vector(int *v);
double *double_vector(int N);
void free_double_vector(double *v);
double **double_matrix(int N, int M);
void free_double_matrix(double **m, int N);
double ***double_tensor(int N, int M, int L);
void free_double_tensor(double ***t, int N, int M);
int **int_matrix(int N, int M);
void free_int_matrix(int **m, int N);


