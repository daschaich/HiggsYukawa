// -----------------------------------------------------------------
// Scalar multiplication on vector
// c <-- s * a
#include "../include/config.h"
#include "../include/so4.h"

void scalar_mult_vec(vector *a, Real s, vector *c) {
#if (DIMF == 4)
  // May be marginally faster than loop
  c->c[0] = s * a->c[0];
  c->c[1] = s * a->c[1];
  c->c[2] = s * a->c[2];
  c->c[3] = s * a->c[3];
#else // DIMF != 4
  register int i;
  for (i = 0; i < DIMF; i++) {
    c->c[i] = s * a->c[i];
#endif
}
// -----------------------------------------------------------------
