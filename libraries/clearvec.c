// -----------------------------------------------------------------
// Clear a vector
#include "../include/config.h"
#include "../include/so4.h"

void clearvec(vector *v) {
#if (DIMF == 4)
  // Manually unroll loop
  v->c[0] = 0.0;
  v->c[1] = 0.0;
  v->c[2] = 0.0;
  v->c[3] = 0.0;
#else  // DIMF != 4
  register int i;
  for (i = 0; i < DIMF; i++)
    v->c[i] = 0.0;
#endif
}
// -----------------------------------------------------------------
