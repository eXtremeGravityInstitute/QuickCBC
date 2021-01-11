#define TSUN 4.92569043916e-6                                           // mass to seconds conversion
#define PI2 9.869604401089357992304940125905      // Pi^2
#define PIn 3.1415926535897931159979634685442     // Pi
#define TPI 6.2831853071795862319959269370884     // 2 Pi
#define RT2PI 2.5066282746310005024               // sqrt(2 Pi)
#define RTPI 1.772453850905516                    // sqrt(Pi)
#define RT2 1.414213562373095                      // sqrt(2)
#define IRT2 0.707106781186547549                      // 1/sqrt(2)
#define TRT2 2.8284271247461901                  // 2*sqrt(2)
#define LN2 0.6931471805599453                 // ln 2
#define CLIGHT 299792458.0                     // m/s
#define h22fac  0.31539156525252               //  2. * sqrt(5. / (64.*PI)) factor for h22 to h conversion
#define MPC 3.08568025e22       // Megaparsec in meters

#define EPOCH_J2000_0_JD 2451545.0         // < Julian Day (UTC) of the J2000.0 epoch (2000 JAN 1 12h UTC).
#define EPOCH_J2000_0_GPS 630763213.0       // < GPS seconds of the J2000.0 epoch (2000 JAN 1 12h UTC). */
#define EPOCH_UNIX_GPS 315964800.0

#define HLdt 0.0106                             // Maximum light travel time allowed H-L
#define HVdt 0.0278                             // Maximum light travel time allowed H-V
#define LVdt 0.0257                            // Maximum light travel time allowed L-V

#define dtp  0.05           // minimum time between likelihiood maximam: defualt 0.05
#define Lpmin 9.0           // minimum log likelihood to keep a peak: defualt 9.0
#define fbw 8.0             // width of frequency bands used in penalized likelihood: defualt 8.0
#define rhocut 4.0          // threshold for penalized likelihood: default 4.0

#define Edt  0.03             // Allow 26 ms either side of geocenter reference time (accounts for Earth radius and some extra)
#define dtmax 0.03            // Maximum time shift relative to geocenter waveform in sky likelihood

#define fref 100.0             // reference freqency for PhenomD waveform

#define fmax 1024.0            // Frequency range used in analysis
#define fmin 16.0
#define fres 4.0               // Frequency resolution used in heterodyne

#define DLmax 1.0e10          // Maximum  distance in pc
#define DLmin 1.0e6           // Minimum  distance in pc

#define cap 100.0     // cap on log likelihood proposal (sets SNR^2/2 threshold)

#define qfac 2.0    // geometric progression in mass ratio

#define ladder 1.15    // temperature ladder
#define twidth 0.1     // width of prior on merger time (used when given non-integer merger time)

// Mc =  Mt * eta^3/5. For eta=1/4 Mc = 0.43527528 Mt. No limit on how small Mc can be if we allow arbitrary q's
// Using q = 10 as the maximum mass ratio we have eta_min = 0.082644628 and Mc = 0.22404 Mt

#define mcmax 160.0   // maximum chirp mass in Msun.  Mc is at most 0.435 times Mt
#define mcmin 0.5    // minimum chirp mass in Msun
#define mtmax 300.0  // maximum total mass in Msun (should be not more than 2.3 times mcmax)
#define mtmin 0.5    // minimum total mass in Msun
#define smax  0.85      // max spin magnitude
#define etamin 0.05     // lower limit on the symmetric mass ratio (18:1 mass ratio)

#define printQ  0     // Set to 1 to print Qscans during cleaning phase
#define Qprint 8.0    // Q used to make output scans

#define Pmimic 1      // set to 1 to mimic precession priors on spin aligned magnitude, zero otherwise.
#define lhold 0       // set to 0 for normal running, set to 1 for prior recovery test


#define NX 7           // number of quasi-intrinsic parameters (Mc, Mt, chi1, chi2, phic, tc, DL)
#define NP 11          // total number of entries in the full parameter arrays NX+4
#define NS 7           // number of quasi-extrinsic parameters (alpha, sin(delta), psi, ellipticity, scale, phi0=2*phic, dt)
#define NC 12           // number of chains
#define NCC 4          // number of cold chains
#define NH 1000        // length of history
#define NQ 4           // number of mass ratios in global proposal
#define NM 1000        // number of chirp masses in global proposal
#define Nsearch 4000     // iterations for search
#define Nintrinsic 10000  // iterations for intrinsic refinement
#define Nsky 1000000      // iterations for sky mapping
#define Nall 400000       // iterations for full MCMC

