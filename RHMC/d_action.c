// -----------------------------------------------------------------
// Measure total action, as needed by the hybrid Monte Carlo algorithm
// When this routine is called the CG should already have been run,
// so that the vector **sol contains (M_adjoint*M+shift[n])^(-1) * src
#include "so4_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Scalar momenta contribution to the action
double d_mom_action() {
  register int i, index;
  register site *s;
  double sum = 0.0;

  FORALLSITES(i, s) {
    for (index = 0; index < NSD; index++)
      sum += (double)(mom[i].e[index] * mom[i].e[index]);
  }
  g_doublesum(&sum);
  return sum;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Scalar contribution to the action
// Identical terms in self-dual matrix cancel factor of 1/2
double d_scalar_action() {
  register int i;
  register site *s;
  int index;
  double s_action = 0.0;

  FORALLSITES(i, s) {
    for (index = 0; index < NSD; index++)
      s_action += (double)(s->sigma.e[index]) * (s->sigma.e[index]);
  }
  g_doublesum(&s_action);
  return s_action;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Fermion contribution to the action
// Include the ampdeg term to allow sanity check that the fermion action
// is 4*volume on average
// Since the pseudofermion src is fixed throughout the trajectory,
// ampdeg actually has no effect on Delta S (checked)
// sol, however, depends on the gauge fields through the CG
double d_fermion_action(vector *src, vector **sol) {
  register int i, j;
  register site *s;
  double sum = 0.0;
  Real tr;
#ifdef DEBUG_CHECK
  double im = 0.0;
#endif

  FORALLSITES(i, s) {
    sum += ampdeg4 * (double)magsq_vec(&(src[i]));
    for (j = 0; j < Norder; j++) {
      tr = dot(&(src[i]), &(sol[j][i]));   // src^dag.sol[j]
      sum += (double)(amp4[j] * tr);
    }
  }
  g_doublesum(&sum);
  return sum;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Print out total action and individual contributions
double d_action(vector **src, vector ***sol) {
  int n;
  double h_act, f_act, total;

  total = d_scalar_action();
  node0_printf("action: scalar %.8g ", total);

  for (n = 0; n < Nroot; n++) {
    f_act = d_fermion_action(src[n], sol[n]);
    node0_printf("fermion%d %.8g ", n, f_act);
    total += f_act;
  }

  h_act = d_mom_action();
  node0_printf("mom %.8g ", h_act);
  total += h_act;
  node0_printf("sum %.8g\n", total);
  return total;
}
// -----------------------------------------------------------------
