// -----------------------------------------------------------------
// Define global scalars and fields in the lattice
#ifndef _LATTICE_H
#define _LATTICE_H

#include "defines.h"
#include "../include/macros.h"    // For MAXFILENAME
#include "../include/io_lat.h"    // For gauge_file
#include "../include/so4.h"
#include "../include/random.h"    // For double_prn
#include "../include/dirs.h"      // For NDIMS
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// The lattice is an array of this site struct
typedef struct {
  short x, y, t;  // Coordinates of this site
  char parity;    // Is it even or odd?
  int index;      // Index in the array

#ifdef SITERAND
  // The state information for a random number generator
  double_prn site_prn;
#endif

  // Staggered phases, which have been absorbed into the matrices
  // Also includes the antiperiodic temporal boundary conditions
  Real phase[NDIMS];

  // Self-dual Hubbard--Stratonovich scalar field
  selfdual sigma;
#ifdef HMC_ALGORITHM
  selfdual old_sigma;  // For accept/reject
#endif
} site;
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Definition of global variables
#ifdef CONTROL
#define EXTERN
#else
#define EXTERN extern
#endif

// Global variables
EXTERN int nx, ny, nt;  // Lattice dimensions
EXTERN int PBC;         // Temporal fermion boundary condition
EXTERN int volume;      // Volume of lattice
EXTERN int iseed;       // Random number seed
EXTERN int warms, trajecs, niter, propinterval;
EXTERN Real traj_length;

// Translate (mu, nu) to linear index of anti-symmetric or self-dual matrix
EXTERN int as_index[DIMF][DIMF];
EXTERN int sd_index[DIMF][DIMF];

EXTERN Real rsqmin, rsqprop;
EXTERN Real G;
EXTERN double g_ssplaq, g_stplaq;
EXTERN double sigmasum;
EXTERN char startfile[MAXFILENAME], savefile[MAXFILENAME];
EXTERN int startflag; // Beginning lattice: CONTINUE, RELOAD, FRESH
EXTERN int fixflag;   // Gauge fixing: COULOMB_GAUGE_FIX, NO_GAUGE_FIX
EXTERN int saveflag;  // 1 if we will save the lattice;
EXTERN int total_iters;   // To be incremented by the multi-mass CG

// Some of these global variables are node dependent
// They are set in "make_lattice()"
EXTERN int sites_on_node;       // Number of sites on this node
EXTERN int even_sites_on_node;  // Number of even sites on this node
EXTERN int odd_sites_on_node;   // Number of odd sites on this node
EXTERN int number_of_nodes;     // Number of nodes in use
EXTERN int this_node;           // Node number of this node

// Stuff for multi-mass CG and RHMC
EXTERN int nsteps[2];           // Fermion and gauge steps
EXTERN Real ampdeg, *amp, *shift;
EXTERN Real ampdeg4, *amp4, *shift4;
EXTERN Real ampdeg8, *amp8, *shift8;
EXTERN int Nroot, Norder;
EXTERN Real snorm, *fnorm, max_sf, *max_ff;

// Momenta and forces for the scalars
selfdual *mom, *fullforce;

// Each node maintains a structure with the pseudorandom number
// generator state
EXTERN double_prn node_prn;

// Persistent fermions for matrix--vector operation
// Used in fermion_op and assemble_fermion_force
EXTERN vector *src;
EXTERN vector *dest;

// Temporary vectors, matrices and Twist_Fermion
EXTERN vector *tempvec;
EXTERN selfdual *tempsd;

EXTERN gauge_file *startlat_p;
EXTERN gauge_file *savelat_p;

// The lattice is a single global variable
// (actually this is the part of the lattice on this node)
EXTERN site *lattice;

// Vectors for addressing
// Generic pointers, for gather routines
// Probably need at least 6=2NDIMS
#define N_POINTERS 6
EXTERN char **gen_pt[N_POINTERS];

#ifdef EIG
// Eigenvalue stuff
EXTERN int Nvec;
EXTERN double *eigVal;
EXTERN vector **eigVec;
EXTERN Real eig_tol;          // Tolerance for the eigenvalue computation
EXTERN int maxIter;           // Maximum iterations
#endif

#ifdef PHASE
// Pfaffian phase stuff
EXTERN int Nmatvecs;                // For timing/counting
EXTERN int ckpt_load, ckpt_save;    // For checkpointing
#endif

#endif // _LATTICE_H
// -----------------------------------------------------------------
