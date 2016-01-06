// -----------------------------------------------------------------
// Main procedure for four-fermion measurements
#define CONTROL
#include "so4_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
int main(int argc, char *argv[]) {
  int prompt;
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
  // TODO: Read in and loop over source points
  int avm_iters = 0;
  step = (int)(nt / 4);
  if (step < 1)   // Make sure we don't hit an infinite loop
    step = 1;
  for (pnt[0] = 0; pnt[0] < nt; pnt[0] += step) {
    for (pnt[1] = 0; pnt[1] < nt; pnt[1] += step) {
      for (pnt[2] = 0; pnt[2] < nt; pnt[2] += step) {
        node0_printf("Source point %d %d %d\n", pnt[0], pnt[1], pnt[2]);
        avm_iters += correlators(pnt);
      }
    }
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
