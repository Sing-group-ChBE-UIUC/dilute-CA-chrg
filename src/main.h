#ifndef Main_h
#define Main_h

// Header files for various libraries used in the code

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <mkl.h>
#include <string.h>
#include <time.h>
#include <getopt.h>

// Constants for the random number generator

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define EPS 1.2e-7
#define RNMX (1.0-EPS)
#define NDIV (1+IMM1/NTAB)

// User defined macros

#define max(a,b) ({ __typeof__ (a) _a = (a); __typeof__ (b) _b = (b); _a > _b ? _a : _b; })
#define min(a,b) ({ __typeof__ (a) _a = (a); __typeof__ (b) _b = (b); _a < _b ? _a : _b; })

// User defined functions
// Explanation of what each function does in main.c

///Charge Variables///
void coulombForce();
void chargeHop() ;
void chargeHop_new() ;
void calcDis() ;
void calcMSD() ;
void calcCH() ;
void printRestart() ; 
int *Charge;  //0 if uncharged //1 if charged
int *Charge_track ;  //tracks charge
int *Charge_old ; //tracks charge
int *Charge_indices ; //tracks charge
int Ncharges ; //Number of charged beads
int i,j,cnumber,counter1,counter2,counter3, counter4, counter5, counter6;
double lambda_d ; //debye screening length
double lambda_b ; //bjerrum length
double kappa_d ; //inverse debye screening length(1/lambda_d)
double txx; //average stress tensor tauxx
double txy; //average stress tensor tauxy
double txz; //average stress tensor tauxz
double tyy; //average stress tensor tauyy
double tyz; //average stress tensor tauyz
double tzz; //average stress tensor tauzz
int chargehop_type ; //0 - no charge hop, 1 - charge hop
int calcmsd ; //0 - don't calculate, 1 - calculate
int MSD_type ; //0 - displacement method | 1 - variable dt method | 2 - set dt method
int MStep ; // do MC every MSteps
int *i_values ;
double barrier ; //barrier for adjacent hopping
double barrier2 ; //"" non-adjacent hopping
char* str ;
char* str2 ;
char* outp ;
char* outp1 ;
char* dtfile ;
FILE *outputfile, *chfile;
FILE *datafile, *OAfile; //*tfile;
FILE *outputfile2;

double MSD_ch;
double *dxt, *dyt, *dzt;
double *rbi; //initial position of charge at timestep
double *prob; //probablility a monomer site is occupied by a charge
double A; //strength of fictitious potential U = Asin(...)
double *E_initial;
double *E_final;
int tcount;
//////////////////////

void readInput(int argc,const char *argv[]);
void allocate();
void printOutput();
void initChains();
void initHI();
void verletList();
void checkVerlet();
void blockNoise();
void bondForce();
void EVForce();
void updatePos();
void sampleHI();
void printHI();
void resetAverage();
void calcStress();
void calcConf();
void printTraj();
long initRan();
float ran1(long *idum); // Uniform RNG [0,1) from numerical recipes textbook
float gasdev();
double *rb; // particle positions
double *ro; // old particle positions at last Verlet list update for use in neighbor tables
double *drb; // deterministic change in particle positions (D.F term in Langevin eqn)
int *nlist,*list; // arrays for Verlet list
double *f; // particle forces
double *Dt; // running average diffusion tensor for the current iteration
double *D; // conformationally averaged diffusion tensor from previous iteration
double *chol; // conformationally averaged decomposed diffusion tensor, "B" in Langevin eqn
// double *R; // Gaussian random variables
double *Z; // Gaussian random variables
char *txt,*xyz,*rg,*ree,*cm,*strs,*ext,*rst,*Davg; // Strings for file names
int N; // number of particles
int tr; // trajectory number
int s; // number of time steps between updating RPY tensor. heuristic, s = 25-50
// need at least s*dt < 0.1 for accuracy. higher s is faster. speedup generally saturates for s > 50
int restart; // 0 -> new simulation. 1 -> restart existing simulation. 2 -> start from stretched conf (only available for lin. eqm.)
int over_write; // 0 -> don't write over previous files for this trajectory. 1 -> write over files for this trajectory if they exist
int printProps; // ts between printing polymer properties
int printxyz; // ts between printing trajectory
int flowType; // 0 -> none 1 -> planar extensional flow (PEF) 2 -> planar shear flow (PSF) 3-> planar rotational flow (PRF)
// HI_type is invalid for the CA method. if it = 0, then FD. if it > 0, then CA HI
// int HI_type; // How the diffusion tensor is calculated. 0 -> no HI. 1 -> RPY
int it; // iteration counter for CA method
int it_start; // starting iteration for CA method
int it_max; // max iteration for CA method
int count_HI; // counter for number of HI samples
int sample_HI; // number of ts between HI samples
int print_HI; // ts between printing average HI
unsigned long t; // current time step
unsigned long tmax; // max time step
unsigned long tstart; // starting ts for restart == 1
unsigned long equil_start; // time steps before startup of flow
unsigned long tcess; // time step for cessation of flow (not currently used)
long *idum; // RNG seed
double epsilon; // WCA and LJ potential prefactor
double kappas; // spring constant
double kappab; // bending constant
double dt; // time step size
double rc,rc2; // WCA/LJ neighbor list cutoff radius, radius squared
double rv,rnew; // Verlet list skin radius
double p,pc; // = sqrt(2.0*dt)
double wall_time; // for executing timing
double decomp_time; // decomposition execution timing
double update_time; // deterministic update execution timing
double edot; // strain rate
double bp; // TEA parameter
double strain; // accumulated strain = edot*t*dt
double strain_cess; // after this accumualted strain, the simulation will stop
double qmax,qmax2; // maximum length of the FENE spring
double ah; // hydrodynamic radius of a bead - leave at one for WCA/LJ potentials
double stress[6]; // stress tensor components
double tau; // user input polymer relaxation time. From simulations or theory (eg Zimm or Rouse)
double tdtau_max; // simulation length divided by polymer relaxation time
double Wi; // flow strength, strain rate made dimensionless by polymer relaxationtime, Wi = edot*t
char tplgy[100],dcmp_meth[100],EV[100],spring[100]; // input strings
FILE *txtfile,*xyzfile,*rgfile,*reefile,*cmfile,*stressfile,*extfile,*rstfile,*Davgfile; // file pointers

#endif // Main_h
