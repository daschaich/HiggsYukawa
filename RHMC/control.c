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
  double dtime, s_act, plus_act, minus_act;
#ifdef CORR
  int j;
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
  s_act = scalar_action(&plus_act, &minus_act) / (double)volume;
  plus_act /= (double)volume;
  minus_act /= (double)volume;
  node0_printf("START %.8g %.8g %.8g\n", plus_act, minus_act, s_act);

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
    // Format: GMES cg_iters plus_act minus_act s_action
    s_act = scalar_action(&plus_act, &minus_act) / (double)volume;
    plus_act /= (double)volume;
    minus_act /= (double)volume;
    node0_printf("GMES %d %.8g %.8g %.8g\n",
                 s_iters, plus_act, minus_act, s_act);

    // Less frequent measurements every "propinterval" trajectories
    if ((traj_done % propinterval) == (propinterval - 1)) {
#ifdef CORR
      // Correlator measurements
      for (j = 0; j < Nsrc; j++) {
        node0_printf("Source point %d %d %d %d\n",
                     pnts[j][0], pnts[j][1], pnts[j][2], pnts[j][3]);
        avm_iters += correlators(pnts[j]);
        avm_iters += condensates();
      }
      Nmeas++;
#endif
    }
  }
  node0_printf("RUNNING COMPLETED\n");

  // Check: compute final scalar action
  s_act = scalar_action(&plus_act, &minus_act) / (double)volume;
  plus_act /= (double)volume;
  minus_act /= (double)volume;
  node0_printf("STOP %.8g %.8g %.8g\n", plus_act, minus_act, s_act);
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
