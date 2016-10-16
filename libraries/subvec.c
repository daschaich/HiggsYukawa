// -----------------------------------------------------------------
// Subtract two vectors
// c <-- a - b
#include "../include/config.h"
#include "../include/so4.h"

void sub_vec(vector *a, vector *b, vector *c) {
#if (DIMF != 4)
  #error "Assuming DIMF=4!"
#endif
  c->c[0] = a->c[0] - b->c[0];
  c->c[1] = a->c[1] - b->c[1];
  c->c[2] = a->c[2] - b->c[2];
  c->c[3] = a->c[3] - b->c[3];
}
// -----------------------------------------------------------------
