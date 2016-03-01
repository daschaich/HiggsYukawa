// -----------------------------------------------------------------
// Main procedure for four-fermion reversibility test
#define CONTROL
#include "so4_includes.h"

// For linking against update_o.c
void ranmom();
double update_scalar(Real eps);
int update_step(double *fnorm, double *snorm, vector **src, vector ***psim);
// -----------------------------------------------------------------



// -----------------------------------------------------------------
int main(int argc, char *argv[]){
#if (DIMF != 4)
  #error "Assuming DIMF=4!"
#endif
  register int i;
  register site *s;
  int prompt, s_iters = 0, j, n;
  Real final_rsq;
  double startAct, midAct, endAct, change, s_act, plus_act, minus_act, dtime;

  // Setup
  setlinebuf(stdout); // DEBUG
  initialize_machine(&argc, &argv);
  // Remap standard I/O
  if (remap_stdio_from_args(argc, argv) == 1)
    terminate(1);

  g_sync();
  prompt = setup();
  epsilon();
  setup_rhmc();

  // Always accept
#ifdef HMC_ALGORITHM
#undef HMC_ALGORITHM
#endif

  // Load input and run
  if (readin(prompt) != 0) {
    node0_printf("ERROR in readin, aborting\n");
    terminate(1);
  }
  dtime = -dclock();

  vector **src = malloc(Nroot * sizeof(**src));
  vector ***psim = malloc(Nroot * sizeof(***psim));
  for (n = 0; n < Nroot; n++) {
    src[n] = malloc(sites_on_node * sizeof(vector));
    psim[n] = malloc(Norder * sizeof(vector*));
    for (j = 0; j < Norder; j++)
      psim[n][j] = malloc(sites_on_node * sizeof(vector));
  }

  // Random momenta
  ranmom();

  // Set up the fermion variables
  // Compute g and src = (Mdag M)^(1 / 8) g
  for (n = 0; n < Nroot; n++)
    s_iters += grsource(src[n]);

  // Do a CG to get psim,
  // rational approximation to (Mdag M)^(-1 / 4) src = (Mdag M)^(-1 / 8) g
  for (j = 0; j < Norder; j++)
    shift[j] = shift4[j];

  // congrad_multi initializes psim
  for (n = 0; n < Nroot; n++)
    s_iters += congrad_multi(src[n], psim[n], niter, rsqmin, &final_rsq);

  // Observables at start of trajectory
  s_act = scalar_action(&plus_act, &minus_act) / (double)volume;
  plus_act /= (double)volume;
  minus_act /= (double)volume;
  startAct = action(src, psim);
  node0_printf("GMES %d %.8g %.8g %.8g\n",
               s_iters, plus_act, minus_act, s_act);
  node0_printf("ACT %.8g\n", startAct);

  // Evolve forward for one trajectory
  snorm = 0.0;
  max_sf = 0.0;
  for (n = 0; n < Nroot; n++) {
    fnorm[n] = 0.0;
    max_ff[n] = 0.0;
  }
  s_iters += update_step(fnorm, &snorm, src, psim);

  // Observables at end of trajectory
  // Scalar action and CG iterations at end of trajectory
  midAct = action(src, psim);
  change = midAct - startAct;
  node0_printf("delta S = %.4g start S = %.12g end S = %.12g\n",
               change, startAct, midAct);

  node0_printf("MONITOR_FORCE_SCALAR   %.4g %.4g\n",
               snorm / (double)(2 * nsteps[0]), max_sf);
  for (n = 0; n < Nroot; n++) {
    node0_printf("MONITOR_FORCE_FERMION%d %.4g %.4g\n",
                 n, fnorm[n] / (double)(2 * nsteps[0]), max_ff[n]);
  }
  s_act = scalar_action(&plus_act, &minus_act) / (double)volume;
  plus_act /= (double)volume;
  minus_act /= (double)volume;
  node0_printf("GMES %d %.8g %.8g %.8g\n",
               s_iters, plus_act, minus_act, s_act);
  node0_printf("ACT %.8g\n", midAct);

  // Reverse momenta and evolve backwards for one trajectory
  // Reverse the momenta, keeping all else the same
  FORALLSITES(i, s) {
    mom[i].e[0] *= -1.0;
    mom[i].e[1] *= -1.0;
    mom[i].e[2] *= -1.0;
    mom[i].e[3] *= -1.0;
    mom[i].e[4] *= -1.0;
    mom[i].e[5] *= -1.0;
  }
  snorm = 0.0;
  max_sf = 0.0;
  for (n = 0; n < Nroot; n++) {
    fnorm[n] = 0.0;
    max_ff[n] = 0.0;
  }
  s_iters = update_step(fnorm, &snorm, src, psim);

  // Observables hopefully back at the start of the trajectory
  endAct = action(src, psim);
  change = endAct - midAct;
  node0_printf("delta S = %.4g start S = %.12g end S = %.12g\n",
               change, midAct, endAct);

  node0_printf("MONITOR_FORCE_SCALAR   %.4g %.4g\n",
               snorm / (double)(2 * nsteps[0]), max_sf);
  for (n = 0; n < Nroot; n++) {
    node0_printf("MONITOR_FORCE_FERMION%d %.4g %.4g\n",
                 n, fnorm[n] / (double)(2 * nsteps[0]), max_ff[n]);
  }
  s_act = scalar_action(&plus_act, &minus_act) / (double)volume;
  plus_act /= (double)volume;
  minus_act /= (double)volume;
  node0_printf("GMES %d %.8g %.8g %.8g\n",
               s_iters, plus_act, minus_act, s_act);
  node0_printf("ACT %.8g\n", endAct);

  // Done
  node0_printf("RUNNING COMPLETED\n");
  node0_printf("Initial action: %.8g\n", startAct);
  node0_printf("Final action:   %.8g\n", endAct);
  node0_printf("Difference:     %.4g\n", fabs(endAct - startAct));
  dtime += dclock();
  node0_printf("\nTime = %.4g seconds\n", dtime);
  node0_printf("total_iters = %d\n", total_iters);
  fflush(stdout);
  for (n = 0; n < Nroot; n++) {
    free(src[n]);
    for (j = 0; j < Norder; j++)
      free(psim[n][j]);
    free(psim[n]);
  }
  free(src);
  free(psim);
  g_sync();         // Needed by at least some clusters
  return 0;
}
// -----------------------------------------------------------------
