// -----------------------------------------------------------------
// Squared magnitude of vector
#include "../include/config.h"
#include "../include/so4.h"

Real magsq_as(antisym *m) {
  register Real sum = m->e[0] * m->e[0];
#if (DIMF == 4)
  // May be marginally faster than loop
  sum += m->e[1] * m->e[1];
  sum += m->e[2] * m->e[2];
  sum += m->e[3] * m->e[3];
  sum += m->e[4] * m->e[4];
  sum += m->e[5] * m->e[5];
#else // DIMF != 4
  register int i;
  for (i = 1; i < NAS; i++)
    sum += m->e[i] * m->e[i]
#endif
  return sum;
}

Real magsq_sd(selfdual *m) {
  register Real sum = m->e[0] * m->e[0];
#if (DIMF == 4)
  // May be marginally faster than loop
  sum += m->e[1] * m->e[1];
  sum += m->e[2] * m->e[2];
#else // DIMF != 4
  register int i;
  for (i = 1; i < NSD; i++)
    sum += m->e[i] * m->e[i]
#endif
  return sum;
}
// -----------------------------------------------------------------
