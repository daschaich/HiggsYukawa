// -----------------------------------------------------------------
// Main procedure for four-fermion measurements
#define CONTROL
#include "so4_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
int main(int argc, char *argv[]) {
  int prompt;
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

  // Load input and run
  if (readin(prompt) != 0) {
    node0_printf("ERROR in readin, aborting\n");
    terminate(1);
  }
  dtime = -dclock();

  // Check: compute initial scalar action
  // This is the only local measurement for now
  s_act = scalar_action();
  node0_printf("START %.8g\n", s_act / (double)volume);

  // Main measurements
#ifdef CORR
  // Correlator measurements
  int j, avm_iters = 0;
  for (j = 0; j < Nsrc; j++) {
    node0_printf("Source point %d %d %d\n",
                 pnts[j][0], pnts[j][1], pnts[j][2]);
    avm_iters += correlators(pnts[j]);
  }
#endif

  node0_printf("RUNNING COMPLETED\n");
#ifdef CORR
  node0_printf("CG iters for measurements: %d\n", avm_iters);
  // total_iters is accumulated in the multiCG itself
  // Should equal total for correlator measurements
  node0_printf("total_iters = %d\n", total_iters);
#endif
  dtime += dclock();
  node0_printf("\nTime = %.4g seconds\n", dtime);
  fflush(stdout);
  g_sync();         // Needed by at least some clusters
  return 0;
}
// -----------------------------------------------------------------
