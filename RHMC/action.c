// -----------------------------------------------------------------
// Measure total action, as needed by the hybrid Monte Carlo algorithm
// When this routine is called the CG should already have been run,
// so that the vector **sol contains (M_adjoint*M+shift[n])^(-1) * src
#include "so4_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Scalar momenta contribution to the action
double mom_action() {
#if (DIMF != 4)
  #error "Assuming DIMF=4!"
#endif
  register int i;
  register site *s;
  double sum = 0.0;

  FORALLSITES(i, s) {
    sum += (double)(mom[i].e[0] * mom[i].e[0]);
    sum += (double)(mom[i].e[1] * mom[i].e[1]);
    sum += (double)(mom[i].e[2] * mom[i].e[2]);
    sum += (double)(mom[i].e[3] * mom[i].e[3]);
    sum += (double)(mom[i].e[4] * mom[i].e[4]);
    sum += (double)(mom[i].e[5] * mom[i].e[5]);
  }
  g_doublesum(&sum);
  return 0.5 * sum;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Scalar contribution to the action
// Separately compute sigma_+, sigma_- and the total
double scalar_action(double *plus_act, double *minus_act) {
#if (DIMF != 4)
  #error "Assuming DIMF=4!"
#endif
  register int i;
  register site *s;
  int a, b, c, d;
  double tr, s_action;
  antisym plus, minus;

  // Initialize plus_act and minus_act
  *plus_act = 0.0;
  *minus_act = 0.0;

  FORALLSITES(i, s) {
    // Clear +/- scalar fields
    for (a = 0; a < NAS; a++) {
      plus.e[a] = 0.0;
      minus.e[a] = 0.0;
    }

    // Compute +/- scalar fields
    for (a = 0; a < DIMF; a++) {
      for (b = a + 1; b < DIMF; b++) {
        plus.e[as_index[a][b]] += s->sigma.e[as_index[a][b]];
        minus.e[as_index[a][b]] += s->sigma.e[as_index[a][b]];
        for (c = 0; c < DIMF; c++) {
          for (d = c + 1; d < DIMF; d++) {
            tr = perm[a][b][c][d] * s->sigma.e[as_index[c][d]];
            plus.e[as_index[a][b]] += tr;
            minus.e[as_index[a][b]] -= tr;
          }
        }   // No half since not double-counting
      }
    }

    // Final factor of 1/2 and add square to plus_act, minus_act
    for (a = 0; a < NAS; a++) {
      *plus_act += 0.125 * plus.e[a] * plus.e[a];
      *minus_act += 0.125 * minus.e[a] * minus.e[a];
    }
  }
  g_doublesum(plus_act);
  g_doublesum(minus_act);

  // Total action is just sum
  s_action = *plus_act + *minus_act;
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
double fermion_action(vector *src, vector **sol) {
  register int i, j;
  register site *s;
  double sum = 0.0;
  Real tr;
#ifdef DEBUG_CHECK
  double im = 0.0;
#endif

  FORALLSITES(i, s) {
    sum += ampdeg2 * (double)magsq_vec(&(src[i]));
    for (j = 0; j < Norder; j++) {
      tr = dot(&(src[i]), &(sol[j][i]));   // src^dag.sol[j]
      sum += (double)(amp2[j] * tr);
    }
  }
  g_doublesum(&sum);
  return sum;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Print out total action and individual contributions
double action(vector **src, vector ***sol) {
  int n;
  double h_act, f_act, total, tr, tr2;

  total = scalar_action(&tr, &tr2);
  node0_printf("action: scalar %.8g ", total);

  for (n = 0; n < Nroot; n++) {
    f_act = fermion_action(src[n], sol[n]);
    node0_printf("fermion%d %.8g ", n, f_act);
    total += f_act;
  }

  h_act = mom_action();
  node0_printf("mom %.8g ", h_act);
  total += h_act;
  node0_printf("sum %.8g\n", total);
  return total;
}
// -----------------------------------------------------------------
