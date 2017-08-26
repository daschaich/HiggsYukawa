// -----------------------------------------------------------------
// Squared magnitude of vector
#include "../include/config.h"
#include "../include/so4.h"

Real magsq_vec(vector *a) {
  register Real sum = a->c[0] * a->c[0];
#if (DIMF == 4)
  // Manually unroll loop
  sum += a->c[1] * a->c[1];
  sum += a->c[2] * a->c[2];
  sum += a->c[3] * a->c[3];
#else // DIMF != 4
  register int i;
  for (i = 1; i < DIMF; i++)
    sum += a->c[i] * a->c[i]
#endif
  return sum;
}
// -----------------------------------------------------------------
