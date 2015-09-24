// -----------------------------------------------------------------
// Main procedure for N=4 SYM evolution and measurements
#define CONTROL
#include "susy_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
int main(int argc, char *argv[]) {
  int prompt, dir;
  int traj_done, s_iters, avs_iters = 0, avm_iters = 0, Nmeas = 0;
  Real f_eps, s_eps;
  double dtime, s_act;

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
  s_act = d_gauge_action();
  node0_printf("START %.8g\n", s_act / (double)volume);

  // Perform warmup trajectories
  f_eps = traj_length / (Real)nsteps[0];
  s_eps = f_eps / (Real)(2 * nsteps[1]);
  node0_printf("f_eps %.4g s_eps %.4g\n", f_eps, s_eps);
  for (traj_done = 0; traj_done < warms; traj_done++)
    update();
  node0_printf("WARMUPS COMPLETED\n");

  // Perform trajectories with measurements
  // But WITHOUT reunitarizing!
  for (traj_done = 0; traj_done < trajecs; traj_done++) {
    s_iters = update();
    avs_iters += s_iters;

    // Do "local" measurements every trajectory!
    // Scalar action and CG iterations
    // Format: GMES cg_iters s_action
    s_act = d_scalar_action();
    node0_printf("GMES %d %.8g\n", s_iters, s_act / (double)volume);

    // Less frequent measurements every "propinterval" trajectories
    if ((traj_done % propinterval) == (propinterval - 1)) {
      // To be added...
    }
  }
  node0_printf("RUNNING COMPLETED\n");

  // Check: compute final scalar action
  s_act = d_gauge_action();
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
