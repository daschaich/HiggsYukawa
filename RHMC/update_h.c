// -----------------------------------------------------------------
// Update the momentum matrices
#include "so4_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Update momenta with the scalar force
double scalar_force(Real eps) {
  register int i;
  register site *s;
  double returnit = 0.0;

// !!!TODO

  return (eps * sqrt(returnit) / volume);  // Probably remove sqrt
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Assemble fermion contributions to scalar force,
//   F = !!!
//       Adj(Ms).D_U M(U, Ub).s - Adj[Adj(Ms).D_Ub M(U, Ub).s]
// "s" is sol while "Ms" is psol
void assemble_fermion_force(vector *sol, vector *psol) {
  register int i;
  register site *s;

  // Clear tempvec to store this force
  FORALLSITES(i, s)
    clear_sd(&(tempsd[i]));

// !!!TODO

  // Final minus sign
  FORALLSITES(i, s)
    scalar_mult_sd(&(tempsd[i]), -1.0, &(tempsd[i]));
}
// -----------------------------------------------------------------




// -----------------------------------------------------------------
// Update the momenta with the fermion force
// Assume that the multiCG has been run (updating the adjoint links),
// with the solution vectors in sol[j]
// Compute force for each pole in tempsd, accumulate into fullforce
// Also use tempvec for temporary storage
double fermion_force(Real eps, vector *src, vector **sol) {
#if (DIMF != 4)
  #error "Assuming DIMF=4!"
#endif
  register int i;
  register site *s;
  int n;
  double returnit = 0.0;

  // Clear the force accumulators
  FORALLSITES(i, s)
    clear_sd(&(fullforce[i]));

  for (n = 0; n < Norder; n++) {
    fermion_op(sol[n], tempvec, PLUS);
    // Makes sense to multiply here by amp4[n]...
    FORALLSITES(i, s)
      scalar_mult_vec(&(tempvec[i]), amp4[n], &(tempvec[i]));

    assemble_fermion_force(sol[n], tempvec);    // Overwrites tempsd
    FORALLSITES(i, s)
      add_sd(&(fullforce[i]), &(tempsd[i]), &(fullforce[i]));
  }

  // Update the momentum from the fermion force -- sum or eps
  // Opposite sign as to gauge force,
  // because dS_G / dU = 2F_g while ds_F / dU = -2F_f
  FORALLSITES(i, s) {
    scalar_mult_add_sd(&(mom[i]), &(fullforce[i]), eps, &(mom[i]));
    returnit += fullforce[i].e[0] + fullforce[i].e[1] + fullforce[i].e[2];
  }
  g_doublesum(&returnit);
  return (eps * returnit / volume);
}
// -----------------------------------------------------------------
