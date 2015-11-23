// -----------------------------------------------------------------
// Update the momentum matrices
#include "so4_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Update the momenta with the (fairly trivial) scalar force
double scalar_force(Real eps) {
  register int i;
  register site *s;
  Real tr = -1.0 * eps;   // To subtract
  double returnit = 0.0;

  // Just subtract sigma from the momenta and compute average force
  FORALLSITES(i, s) {
    scalar_mult_add_as(&(mom[i]), &(s->sigma), tr, &(mom[i]));
    returnit += magsq_as(&(s->sigma));
  }
  g_doublesum(&returnit);
  return (eps * sqrt(returnit) / volume);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Update the momenta with the fermion force
//   -F_{ab} = Ms^a s^b + 0.5 eps_{abcd} Ms^c s^d
// "s" is sol while "Ms" is tempvec
// Assume that the multiCG has been run, with the solution vectors in sol[j]
// Accumulate force for each pole, using tempvec for temporary storage
double fermion_force(Real eps, vector *src, vector **sol) {
#if (DIMF != 4)
  #error "Assuming DIMF=4!"
#endif
  register int i;
  register site *s;
  int n, a, b, c, d, p;
  Real tr, halfG = 0.5 * G;
  double returnit = 0.0;
  antisym tempas;

  // Clear the force accumulators
  FORALLSITES(i, s)
    clear_as(&(force[i]));

  // Accumulate vector product in force
  for (n = 0; n < Norder; n++) {
    fermion_op(sol[n], tempvec, PLUS);
    FORALLSITES(i, s) {
      // Makes sense to multiply here by amp4[n]...
      scalar_mult_vec(&(tempvec[i]), amp4[n], &(tempvec[i]));

      clear_as(&tempas);
      for (a = 0; a < DIMF; a++) {
        for (b = a + 1; b < DIMF; b++) {
          p = as_index[a][b];
          tempas.e[p] += tempvec[i].c[a] * sol[n][i].c[b];
          tempas.e[p] -= tempvec[i].c[b] * sol[n][i].c[a];
          for (c = 0; c < DIMF; c++) {
            for (d = c + 1; d < DIMF; d++) {
              tr = perm[a][b][c][d];    // No half since not double-counting
              tempas.e[p] += tr * tempvec[i].c[c] * sol[n][i].c[d];
              tempas.e[p] -= tr * tempvec[i].c[d] * sol[n][i].c[c];
            }
          }
        }
      }
      scalar_mult_add_as(&(force[i]), &tempas, halfG, &(force[i]));
    }
  }

  // Update the momenta from the fermion force
  // Note opposite sign compared to gauge force...
  FORALLSITES(i, s) {
    scalar_mult_add_as(&(mom[i]), &(force[i]), eps, &(mom[i]));
    returnit += magsq_as(&(force[i]));
  }
  g_doublesum(&returnit);
  return (eps * sqrt(returnit) / volume);
}
// -----------------------------------------------------------------
