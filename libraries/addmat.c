// -----------------------------------------------------------------
// Add two matrices
// c <-- a + b
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/so4.h"

void add_as(antisym *a, antisym *b, antisym *c) {
#if (DIMF != 4)
  #error "Assuming DIMF=4!"
#endif
  c->e[0] = a->e[0] + b->e[0];
  c->e[1] = a->e[1] + b->e[1];
  c->e[2] = a->e[2] + b->e[2];
  c->e[3] = a->e[3] + b->e[3];
  c->e[4] = a->e[4] + b->e[4];
  c->e[5] = a->e[5] + b->e[5];
}

void add_sd(selfdual *a, selfdual *b, selfdual *c) {
#if (DIMF != 4)
  #error "Assuming DIMF=4!"
#endif
  c->e[0] = a->e[0] + b->e[0];
  c->e[1] = a->e[1] + b->e[1];
  c->e[2] = a->e[2] + b->e[2];
}

// -----------------------------------------------------------------
