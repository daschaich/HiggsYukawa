// -----------------------------------------------------------------
// Print the given irrep matrix
#include "../include/config.h"
#include <stdio.h>
#include "../include/complex.h"
#include "../include/so4.h"

void dumpas(antisym *m) {
#if (DIMF != 4)
  #error "Assuming DIMF=4!"
#endif
  printf("%.4g  %.4g  %.4g  %.4g  %.4g  %.4g\n",
         v->e[0], v->e[1], v->e[2], v->e[3], v->e[4], v->e[5]);
}

void dumpsd(selfdual *m) {
#if (DIMF != 4)
  #error "Assuming DIMF=4!"
#endif
  printf("%.4g  %.4g  %.4g\n", v->e[0], v->e[1], m->e[2]);
}
// -----------------------------------------------------------------
