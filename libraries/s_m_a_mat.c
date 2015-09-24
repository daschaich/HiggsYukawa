// -----------------------------------------------------------------
// Add result of scalar multiplication on matrices
// c <-- a + s * b
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/so4.h"

void scalar_mult_add_as(antisym *a, antisym *b, Real s, antisym *c) {
#if (DIMF != 4)
  #error "Assuming DIMF=4!"
#endif
  c->e[0] = a->e[0] + s * b->e[0];
  c->e[1] = a->e[1] + s * b->e[1];
  c->e[2] = a->e[2] + s * b->e[2];
  c->e[3] = a->e[3] + s * b->e[3];
  c->e[4] = a->e[4] + s * b->e[4];
  c->e[5] = a->e[5] + s * b->e[5];
}

void scalar_mult_add_sd(selfdual *a, selfdual *b, Real s, selfdual *c) {
#if (DIMF != 4)
  #error "Assuming DIMF=4!"
#endif
  c->e[0] = a->e[0] + s * b->e[0];
  c->e[1] = a->e[1] + s * b->e[1];
  c->e[2] = a->e[2] + s * b->e[2];
}
// -----------------------------------------------------------------
