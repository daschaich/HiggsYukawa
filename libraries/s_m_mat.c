// -----------------------------------------------------------------
// Scalar multiplication on matrices
// b <-- s * a
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/so4.h"

void scalar_mult_as(antisym *a, Real s, antisym *b) {
#if (DIMF != 4)
  #error "Assuming DIMF=4!"
#endif
  b->e[0] = s * a->e[0];
  b->e[1] = s * a->e[1];
  b->e[2] = s * a->e[2];
  b->e[3] = s * a->e[3];
  b->e[4] = s * a->e[4];
  b->e[5] = s * a->e[5];
}

void scalar_mult_sd(selfdual *a, Real s, selfdual *b) {
#if (DIMF != 4)
  #error "Assuming DIMF=4!"
#endif
  b->e[0] = s * a->e[0];
  b->e[1] = s * a->e[1];
  b->e[2] = s * a->e[2];
}
// -----------------------------------------------------------------
