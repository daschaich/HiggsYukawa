// -----------------------------------------------------------------
// Main procedure for four-fermion eigenvalues
#define CONTROL
#include "so4_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
int main(int argc, char *argv[]) {
  int prompt;
  double dtime, s_act, plus_act, minus_act;
  int ivec, total_iters = 0;
#ifndef EIG
  node0_printf("Don't use control_eig unless compiling with -DEIG!\n");
  terminate(1);
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
  s_act = scalar_action(&plus_act, &minus_act) / (double)volume;
  plus_act /= (double)volume;
  minus_act /= (double)volume;
  node0_printf("START %.8g %.8g %.8g\n", plus_act, minus_act, s_act);

  // Main measurement: PRIMME eigenvalues
  // Allocate eigenvectors
  Nvec += 2;    // !!! Extra pair may help quartet formation...
  eigVal = malloc(Nvec * sizeof(*eigVal));
  eigVec = malloc(Nvec * sizeof(*eigVec));
  for (ivec = 0; ivec < Nvec; ivec++)
    eigVec[ivec] = malloc(sites_on_node * sizeof(vector));

  // Calculate and print smallest eigenvalues,
  // checking |D^dag D phi - lambda phi|^2
  total_iters = make_evs(Nvec, eigVec, eigVal, 1);

  // Check matrix elements of D with DDdag eigenmodes
  // The eigenvalues should be paired, with each pair producing
  // positive/negative matrix elements
  // In principle, one could tighten eig_tol until all pairs are found
  // For now we just print them all out to check offline
  // !!! Omit extra pair of eigenvalues
  check_Dmat(Nvec - 2, eigVec);

  // Calculate and print largest eigenvalues, for tuning RHMC
  // Don't need to compute many here...
  if (Nvec > 14) {    // !!! Extra pair of eigenvalues
    Nvec = 14;
    free(eigVal);
    eigVal = malloc(Nvec * sizeof(*eigVal));

    for (ivec = 0; ivec < Nvec; ivec++)
      free(eigVec[ivec]);
    free(eigVec);

    eigVec = malloc(Nvec * sizeof(*eigVec));
    for (ivec = 0; ivec < Nvec; ivec++)
      eigVec[ivec] = malloc(sites_on_node * sizeof(vector));
  }
  total_iters += make_evs(Nvec, eigVec, eigVal, -1);

  node0_printf("RUNNING COMPLETED\n");
  dtime += dclock();
  node0_printf("\nTime = %.4g seconds\n", dtime);
  node0_printf("total_iters = %d\n", total_iters);
  fflush(stdout);

  free(eigVal);
  for (ivec = 0; ivec < Nvec; ivec++)
    free(eigVec[ivec]);
  free(eigVec);
  g_sync();         // Needed by at least some clusters
  return 0;
}
// -----------------------------------------------------------------
