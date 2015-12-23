// -----------------------------------------------------------------
// Main procedure for four-fermion evolution and measurements
#define CONTROL
#include "so4_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
int main(int argc, char *argv[]) {
  int prompt;
  int traj_done, s_iters, avs_iters = 0, avm_iters = 0, Nmeas = 0;
  Real f_eps, s_eps;
  double dtime, s_act;

#ifdef CORR
  int step;
  int *pnt = malloc(NDIMS * sizeof(*pnt));  // Measurement source point
#endif

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

  // Load input and run
  if (readin(prompt) != 0) {
    node0_printf("ERROR in readin, aborting\n");
    terminate(1);
  }
  dtime = -dclock();

  // Check: compute initial scalar action
  s_act = scalar_action();
  node0_printf("START %.8g\n", s_act / (double)volume);

  // Perform warmup trajectories
  f_eps = traj_length / (Real)nsteps[0];
  s_eps = f_eps / (Real)(2 * nsteps[1]);
  node0_printf("f_eps %.4g s_eps %.4g\n", f_eps, s_eps);
  for (traj_done = 0; traj_done < warms; traj_done++)
    update();
  node0_printf("WARMUPS COMPLETED\n");

  // Perform trajectories with measurements
  for (traj_done = 0; traj_done < trajecs; traj_done++) {
    s_iters = update();
    avs_iters += s_iters;

    // Do "local" measurements every trajectory!
    // Scalar action and CG iterations
    // Format: GMES cg_iters s_action
    s_act = scalar_action();
    node0_printf("GMES %d %.8g\n", s_iters, s_act / (double)volume);

    // Less frequent measurements every "propinterval" trajectories
    if ((traj_done % propinterval) == (propinterval - 1)) {
#ifdef CORR
      // Correlator measurements
      // TODO: Read in and loop over source points
      step = (int)(nt / 4);
      if (step < 1)   // Make sure we don't hit an infinite loop
        step = 1;
      for (pnt[2] = 0; pnt[2] < nt; pnt[2] += step) {
        pnt[0] = pnt[2] % nx;
        pnt[1] = pnt[2] % ny;
        node0_printf("Source point %d %d %d\n", pnt[0], pnt[1], pnt[2]);
        avm_iters += correlators(pnt);
      }
      Nmeas++;
#endif
    }
  }
  node0_printf("RUNNING COMPLETED\n");

  // Check: compute final scalar action
  s_act = scalar_action();
  node0_printf("STOP %.8g\n", s_act / (double)volume);
  node0_printf("Average CG iters for steps: %.4g\n",
               (double)avs_iters / trajecs);
  if (Nmeas > 0) {
    node0_printf("Average CG iters for measurements: %.4g\n",
                 (double)avm_iters / Nmeas);
  }
  dtime += dclock();
  node0_printf("\nTime = %.4g seconds\n", dtime);
  // total_iters is accumulated in the multiCG itself
  // Should equal total for steps plus measurements
  node0_printf("total_iters = %d\n\n", total_iters);
  fflush(stdout);

  // Save lattice if requested
  if (saveflag != FORGET)
    save_lattice(saveflag, savefile);
  g_sync();         // Needed by at least some clusters
  return 0;
}
// -----------------------------------------------------------------
