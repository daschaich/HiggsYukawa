// -----------------------------------------------------------------
// Print the given matrix
#include "../include/config.h"
#include <stdio.h>
#include "../include/so4.h"

void dumpas(antisym *m) {
#if (DIMF != 4)
  #error "Assuming DIMF=4!"
#endif
  printf("%.4g  %.4g  %.4g  %.4g  %.4g  %.4g\n",
         m->e[0], m->e[1], m->e[2], m->e[3], m->e[4], m->e[5]);
}

void dumpsd(selfdual *m) {
#if (DIMF != 4)
  #error "Assuming DIMF=4!"
#endif
  printf("%.4g  %.4g  %.4g\n", m->e[0], m->e[1], m->e[2]);
}
// -----------------------------------------------------------------
