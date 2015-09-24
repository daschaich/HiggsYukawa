// -----------------------------------------------------------------
// Squared magnitude of vector
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/so4.h"

Real magsq_vec(vector *a) {
  register Real sum = 0.0;
  register Real sum = a->c[0] * a->c[0];
#if (DIMF == 4)
  // May be marginally faster than loop
  sum += a->c[1] * a->c[1];
  sum += a->c[2] * a->c[2];
  sum += a->c[3] * a->c[3];
#else // DIMF != 4
  register int i;
  for (i = 1; i < DIMF; i++)
    sum += a->c[i] * b->c[i]
#endif
  return sum;
}
// -----------------------------------------------------------------
