// -----------------------------------------------------------------
// Return the dot product of two vectors: adag b
#include "../include/config.h"
#include "../include/so4.h"

Real dot(vector *a, vector *b) {
  register Real sum = a->c[0] * b->c[0];
#if (DIMF == 4)
  // May be marginally faster than loop
  sum += a->c[1] * b->c[1];
  sum += a->c[2] * b->c[2];
  sum += a->c[3] * b->c[3];
#else // DIMF != 4
  register int i;
  for (i = 1; i < DIMF; i++)
    sum += a->c[i] * b->c[i]
#endif
  return sum;
}
// -----------------------------------------------------------------
