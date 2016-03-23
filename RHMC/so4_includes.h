// -----------------------------------------------------------------
// Include files for supersymmetric evolution
#include "../include/config.h"  // Keep this first
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>             // For print_var.c, setup.c, gauge_info.c
#include "../include/so4.h"
#include "../include/macros.h"
#include "lattice.h"
#include "../include/comdefs.h"
#include "../include/io_lat.h"
#include "../include/generic.h"
#include "../include/dirs.h"
#include "../include/field_alloc.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Prototypes for functions in high level code
int setup();
void setup_rhmc();
int readin(int prompt);
int update();
void update_h(Real eps);
void update_u(Real eps);

// Gaussian random source
int grsource(vector *src);

// Action routines
double action(vector **src, vector ***sol);
double scalar_action(double *plus_act, double *minus_act);
double fermion_action();

// Force routines
double scalar_force(Real eps);
double fermion_force(Real eps, vector *source, vector **psim);

// Fermion matrix--vector operators (D & D^2) and multi-mass CG
void fermion_op(vector *src, vector *dest, int sign);
void DSq(vector *src, vector *dest);
int congrad_multi(vector *src, vector **psim,
                  int MaxCG, Real RsdCG, Real *size_r);

// Epsilon tensor
Real order(int i, int j, int k, int l);
void epsilon();

// Utility to copy scalar field in site
void scalar_field_copy(field_offset src, field_offset dest);
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// More measurements
#ifdef CORR
// Four-fermion condensate and its susceptibility
int condensates();

// Two- and four-fermion correlator measurements
int correlators(int *pnt);
#endif
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Eigenvalue routines
#ifdef EIG
int make_evs(int Nvec, vector **eigVec, double *eigVal, int flag);
void check_Dmat(int Nvec, vector **eigVec);

// Use LAPACK to diagonalize <psi_j | D | psi_i>
// on the subspace of Ddag.D eigenvalues psi
// http://www.physics.orst.edu/~rubin/nacphy/lapack/routines/dgeev.html
// First two arguments turn off eigenvector computations
// Third and fifth arguments are the dimensions of the matrix
// Fourth argument is that matrix, which will be overwritten
// Sixth argument holds the computed eigenvalues
// Seventh argument holds their imaginary parts
// Eighth and tenth arguments are eigenvectors
// Ninth and eleventh arguments are the dimensions of the eigenvectors
// Twelfth argument is real workspace, of size given by the thirteenth argument
// Final argument reports success or information about failure
void dgeev_(char *doL, char *doR, int *N1, double *store, int *N2,
            double *eigs, double *imag,
            double *dumL, int *NL, double *dumR, int *NR,
            double *work, int *Nwork, int *stat);
#endif
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Pfaffian
#ifdef PFAFF
void pfaff();
#endif
// -----------------------------------------------------------------
