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

#define sthresh 10.0
#define warm 6.0

typedef struct {
    int ncof,ioff,joff;
    double *cc,*cr;
} wavefilt;

static const gsl_rng_type *rngtype;
static const gsl_rng *rng;

void Inverse(double **M, double **IM, int d);
void whiten(double *data, double *Sn, int N);
void qscan(double *data, double *Sn, double Tobs, int N);
void qscanf(double *data, double *Sn, double Tobs, int N);
void qscanres(double *data, double *signal, double *Sn, double Tobs, int N);
double fQNM(double m1, double m2, double chi1, double chi2);
double t_at_f(double *param, double f);
void t_of_f(double *param, double *times, double *freqs, int N);
double f_at_t(double *param, double t);
void f_of_t(double *param, double *times, double *freqs, int N);
double fringdown(double *param);
double fbegin(double *param);
void tukey(double *data, double alpha, int N);
void tukey_scale(double *s1, double *s2, double alpha, int N);
void pbt_shift(double *corr, double *corrf, double *data1, double *data2, double *Sn, int imin, int imax, int N);
double fourier_nwip(double *a, double *b, double *Sn, int imin, int imax, int N);
double globe(double ***global, double *max, double *min, double Tobs, double *params, int N, gsl_rng *r);
double globeden(double ***global, double *max, double *min, double Tobs, double *params, int N);
void max_array_element(double *max, int *index, double *array, int n);
void phase_blind_time_shift(double *corr, double *corrf, double *data1, double *data2, int n);
double Fmag(double *sky, double GMST, int id);
void upsample(int n, double Tobs, int *nt, int *bn);
double det(double *A, int N);
void F_ant(double psi, double alpha, double sindelta, double GMST, double *Fplus, double *Fcross, int id);
void F_LHO(double psi, double alpha, double sindelta, double GMST, double *Fplus, double *Fcross);
void F_LLO(double psi, double alpha, double sindelta, double GMST, double *Fplus, double *Fcross);
void F_VIRGO(double psi, double alpha, double sindelta, double GMST, double *Fplus, double *Fcross);
void Times(double alpha, double sindelta, double GMST, double *dtimes);
void Ring(double *skyx, double *skyy, int d1, int d2, double GMST, gsl_rng * r);
void RingPlot(double *paramx, int d1, int d2, double GMST, double *theta, double *phi, int M);
double f_nwip(double *a, double *b, int n);
void Ring_all(double *paramx, double *paramy, int d1, int d2, double GMST, gsl_rng * r);
double z_DL(double DL);
double DL_fit(double z);
void TransformC(double *a, double *freqs, double **tf, double **tfR, double **tfI, double Q, double Tobs, int n, int m);
void layerC(double *a, double f, double *tf, double *tfR, double *tfI, double Q, double Tobs, double fix, int n);
void SineGaussianC(double *hs, double *sigpar, double Tobs, int N);
void specest(double *data, double *Hf, int N, int Ns, double dt, double fmx, double *SN, double *SM, double *PS);
void bwbpf(double *in, double *out, int fwrv, int M, int n, double s, double f1, double f2);
void CubicSplineGSL(int N, double *x, double *y, int Nint, double *xint, double *yint);
void clean(double *D, double *Draw, double *Hf, double *sqf, double *freqs, double *Sn, double *specD, double *sspecD, double df, double Q, double Tobs, double scale, double alpha, int Nf, int N, int imin, int imax, double *SNR, int pflag);
void spectrum(double *data, double *S, double *Sn, double *Smooth, double df, int N);
void recursive_phase_evolution(double dre, double dim, double *cosPhase, double *sinPhase);
void SineGaussianF(double *hs, double *sigpar, double Tobs, int NMAX);
double Getscale(double *freqs, double Q, double Tobs, double fmx, int n, int m);

void makespec(int Nspline, int Nlines, double *ffit, double *Sspline, double *linef, double *lineh,  double *lineQ, double *linew, double *deltafmax, double *SN, double Tobs,  int N);

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
double ****double_quad(int N, int M, int L, int K);
void free_double_quad(double ****t, int N, int M, int L);

