// -----------------------------------------------------------------
// Print the given vector
#include "../include/config.h"
#include <stdio.h>
#include "../include/so4.h"

void dumpvec(vector *v) {
#if (DIMF == 4)
  // Manually unroll loop
  printf("%.4g  %.4g  %.4g  %.4g\n",
         v->c[0], v->c[1], v->c[2], v->c[3]);
#else  // DIMF != 4
  register int i;
  printf("%.4g", v->c[0]);
  for (i = 1; i < DIMF; i++)
    printf("  %.4g", v->c[i]);
  printf("\n");
#endif
}
// -----------------------------------------------------------------
