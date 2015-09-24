// -----------------------------------------------------------------
// Dirac operator and other helper functions for the action and force
#include "so4_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Copy self-dual scalar fields in the site struct
void scalar_field_copy(field_offset src, field_offset dest) {
  register int i;
  register site *s;

  FORALLSITES(i, s)
    sd_copy((selfdual *)F_PT(s, src), (selfdual *)F_PT(s, dest));
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Matrix--vector operation
// Applies either the operator (sign = 1) or its adjoint (sign = -1)
void fermion_op(vector *src, vector *dest, int sign) {
  register int i, j;
  register site *s;
  int mu;

  // Copy src TwistFermion into fieldwise site, link and plaq fermions
  // All of the latter are overwritten -- don't need to clear explicitly
  if (sign == 1) {
    FORALLSITES(i, s)
      vec_copy(&(src[i]), &(dest[i]));
  }
  else if (sign == -1) {
    FORALLSITES(i, s)
      vec_copy(&(src[i]), &(dest[i]));
  }
  else {
    node0_printf("Error: incorrect sign in fermion_op: %d\n", sign);
    terminate(1);
  }

  // Assemble separate routines for each term in the fermion operator
//  !!!TODO

  // Copy local plaquette, link and site fermions into dest TwistFermion
  if (sign == 1) {
    FORALLSITES(i, s)
      vec_copy(&(src[i]), &(dest[i]));
  }
  else if (sign == -1) {
    FORALLSITES(i, s)
      vec_copy(&(src[i]), &(dest[i]));
  }
}
// -----------------------------------------------------------------
