// -----------------------------------------------------------------
// Three-dimensional four-fermion system setup
#include "params.h"
#include "so4_includes.h"
#define IF_OK if(status==0)

// Each node has a params structure for passing simulation parameters
params par_buf;
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// On node zero, read lattice size and seed, and send to others
int initial_set() {
  int prompt, status;
  if (mynode() == 0) {
    // Print banner
    // stringification kludge from GNU preprocessor manual
    // http://gcc.gnu.org/onlinedocs/cpp/Stringification.html
#define XSTR(s) STR(s)
#define STR(s) #s
    // end kludge
    printf("Three-dimensional four-fermion SO(%d) system\n", DIMF);
    printf("Microcanonical simulation with refreshing\n");
    printf("Machine = %s, with %d nodes\n", machine_type(), numnodes());
#ifdef HMC_ALGORITHM
    printf("Hybrid Monte Carlo algorithm\n");
#endif
#ifdef PHI_ALGORITHM
    printf("Phi algorithm\n");
#else   // Quit!
    printf("Only works for phi algorithm\n");
    exit(1);
#endif
    time_stamp("start");
    status = get_prompt(stdin,  &prompt);

    IF_OK status += get_i(stdin, prompt, "nx", &par_buf.nx);
    IF_OK status += get_i(stdin, prompt, "ny", &par_buf.ny);
    IF_OK status += get_i(stdin, prompt, "nt", &par_buf.nt);
    IF_OK status += get_i(stdin, prompt, "PBC", &par_buf.PBC);
    IF_OK status += get_i(stdin, prompt, "iseed", &par_buf.iseed);

    // Number of Nth roots to take (in addition to 1/4 power)
    IF_OK status += get_i(stdin, prompt, "Nroot", &par_buf.Nroot);

    // RHMC degree
    IF_OK status += get_i(stdin, prompt, "Norder", &par_buf.Norder);

    if (status > 0)
      par_buf.stopflag = 1;
    else
      par_buf.stopflag = 0;

    if (par_buf.PBC >= 0)
      printf("Periodic temporal boundary conditions\n");
    if (par_buf.PBC < 0)
      printf("Antiperiodic temporal boundary conditions\n");
  }

  // Broadcast parameter buffer from node 0 to all other nodes
  broadcast_bytes((char *)&par_buf, sizeof(par_buf));
  if (par_buf.stopflag != 0)
    normal_exit(0);

  nx = par_buf.nx;
  ny = par_buf.ny;
  nt = par_buf.nt;
  PBC = par_buf.PBC;
  iseed = par_buf.iseed;



  // Set up stuff for RHMC and multi-mass CG
  Nroot = par_buf.Nroot;
  fnorm = malloc(Nroot * sizeof(fnorm));
  max_ff = malloc(Nroot * sizeof(max_ff));

  Norder = par_buf.Norder;
  amp = malloc(Norder * sizeof(amp));
  amp4 = malloc(Norder * sizeof(amp4));
  amp8 = malloc(Norder * sizeof(amp8));
  shift = malloc(Norder * sizeof(shift));
  shift4 = malloc(Norder * sizeof(shift4));
  shift8 = malloc(Norder * sizeof(shift8));

  this_node = mynode();
  number_of_nodes = numnodes();
  volume = nx * ny * nt;
  total_iters = 0;
  return prompt;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Set up Kogut--Susskind phase factors, sum_{nu < mu} (-1)^{i[nu]}
void setup_phases() {
  register int i;
  register site *s;

  FORALLSITES(i, s) {
    s->phase[TUP] = 1.0;
    if ((s->t) % 2 == 1)
      s->phase[XUP] = -1.0;
    else
      s->phase[XUP] = 1.0;

    if ((s->x) % 2 == 1)
      s->phase[YUP] = -s->phase[XUP];
    else
      s->phase[YUP] = s->phase[XUP];
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Allocate space for fields
void make_fields() {
  double size = (double)(2.0 * sizeof(vector));
  FIELD_ALLOC(src, vector);
  FIELD_ALLOC(dest, vector);

  // Momenta and forces for the scalars
  FIELD_ALLOC(mom, antisym);
  FIELD_ALLOC(force, antisym);

  // Temporary vector and matrix
  size += (double)(sizeof(vector) + sizeof(antisym));
  FIELD_ALLOC(tempvec, vector);
  FIELD_ALLOC(tempas, antisym);

  size *= sites_on_node;
  node0_printf("Mallocing %.1f MBytes per core for fields\n", size / 1e6);
#ifdef PHASE
  // Total number of matvecs is (volume * DIMF)^2 / 4
  Nmatvecs = volume * volume * DIMF * DIMF / 4;

  // Total size of matrix is (volume * DIMF) x (sites_on_node * DIMF)
  size = (double)(volume * DIMF * sites_on_node * DIMF * sizeof(complex));
  node0_printf("Q has %d columns --> %d matvecs and %.1f MBytes per core...",
               volume * 16 * DIMF, Nmatvecs, size / 1e6);
#endif
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
int setup() {
  int prompt;

  // Print banner, get volume and seed
  prompt = initial_set();
  // Initialize the node random number generator
  initialize_prn(&node_prn, iseed, volume + mynode());
  // Initialize the layout functions, which decide where sites live
  setup_layout();
  // Allocate space for lattice, set up coordinate fields
  make_lattice();
  // Set up neighbor pointers and comlink structures
  make_nn_gathers();
  // Set up phases
  setup_phases();
  // Allocate space for fields
  make_fields();

  return prompt;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Read in parameters for Monte Carlo
int readin(int prompt) {
  // prompt=1 indicates prompts are to be given for input
  int status;
  Real x;

  // On node zero, read parameters and send to all other nodes
  if (this_node == 0) {
    printf("\n\n");
    status = 0;

    // Warms, trajecs
    IF_OK status += get_i(stdin, prompt, "warms", &par_buf.warms);
    IF_OK status += get_i(stdin, prompt, "trajecs", &par_buf.trajecs);
    IF_OK status += get_f(stdin, prompt, "traj_length", &par_buf.traj_length);

    // Number of fermion and gauge steps
    IF_OK status += get_i(stdin, prompt, "nstep", &par_buf.nsteps[0]);
    IF_OK status += get_i(stdin, prompt, "nstep_gauge", &par_buf.nsteps[1]);

    // Trajectories between propagator measurements
    IF_OK status += get_i(stdin, prompt, "traj_between_meas",
                          &par_buf.propinterval);

    // Four-fermion coupling
    IF_OK status += get_f(stdin, prompt, "G", &par_buf.G);

    // Maximum conjugate gradient iterations
    IF_OK status += get_i(stdin, prompt, "max_cg_iterations", &par_buf.niter);

    // Error per site for conjugate gradient
    IF_OK {
      status += get_f(stdin, prompt, "error_per_site", &x);
      par_buf.rsqmin = x;
    }

#ifdef EIG
    // Number of eigenvalues to calculate
    IF_OK status += get_i(stdin, prompt, "Nvec", &par_buf.Nvec);
    IF_OK status += get_f(stdin, prompt, "eig_tol", &par_buf.eig_tol);
    IF_OK status += get_i(stdin, prompt, "maxIter", &par_buf.maxIter);
#endif

#ifdef PHASE
    // Optional checkpointing for pfaffian computation
    IF_OK status += get_i(stdin, prompt, "ckpt_load", &par_buf.ckpt_load);
    IF_OK status += get_i(stdin, prompt, "ckpt_save", &par_buf.ckpt_save);
#endif

    // Find out what kind of starting lattice to use
    IF_OK status += ask_starting_lattice(stdin,  prompt, &par_buf.startflag,
                                         par_buf.startfile);

    // Find out what to do with lattice at end
    IF_OK status += ask_ending_lattice(stdin,  prompt, &(par_buf.saveflag),
                                       par_buf.savefile);

    if (status > 0)
      par_buf.stopflag = 1;
    else
      par_buf.stopflag = 0;
  }

  // Broadcast parameter buffer from node0 to all other nodes
  broadcast_bytes((char *)&par_buf, sizeof(par_buf));
  if (par_buf.stopflag != 0)
    normal_exit(0);

  warms = par_buf.warms;
  trajecs = par_buf.trajecs;
  traj_length = par_buf.traj_length;
  nsteps[0] = par_buf.nsteps[0];
  nsteps[1] = par_buf.nsteps[1];

  propinterval = par_buf.propinterval;
  fixflag = par_buf.fixflag;
  niter = par_buf.niter;
  rsqmin = par_buf.rsqmin;

  G = par_buf.G;

#ifdef EIG
  Nvec = par_buf.Nvec;
  eig_tol = par_buf.eig_tol;
  maxIter = par_buf.maxIter;
#endif
#ifdef PHASE
  ckpt_load = par_buf.ckpt_load;
  ckpt_save = par_buf.ckpt_save;
#endif

  startflag = par_buf.startflag;
  saveflag = par_buf.saveflag;
  strcpy(startfile, par_buf.startfile);
  strcpy(savefile, par_buf.savefile);

  // Do whatever is needed to get lattice
  startlat_p = reload_lattice(startflag, startfile);
  return 0;
}
// -----------------------------------------------------------------
