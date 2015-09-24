// -----------------------------------------------------------------
// Add result of scalar multiplication on vector
// c <-- a + s * b
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/so4.h"

void scalar_mult_add_vec(vector *a, vector *b, Real s, vector *c) {
#if (DIMF == 4)
  // May be marginally faster than loop
  c->c[0] = a->c[0] + s * b->c[0];
  c->c[1] = a->c[1] + s * b->c[1];
  c->c[2] = a->c[2] + s * b->c[2];
  c->c[3] = a->c[3] + s * b->c[3];
#else // DIMF != 4
  register int i;
  for (i = 0; i < DIMF; i++)
    c->c[i] = a->c[i] + s * b->c[i];
#endif
}
// -----------------------------------------------------------------
