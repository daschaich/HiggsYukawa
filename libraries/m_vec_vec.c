// -----------------------------------------------------------------
// Vector--vector outer product for SO(4) fermion force
// See ../RHMC/setup_perm.c for as_index and sd_index mappings
// c <-- a * b
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/so4.h"

void mult_vec_vec(vector *a, vector *b, antisym *c) {
#if (DIMF != 4)
  #error "Assuming DIMF=4!"
#endif
  // First straightforward psol^a sol^b
  c->e[0] = a->c[0] * b->c[1] - a->c[1] * b->c[0];  // 01 --> 0
  c->e[1] = a->c[0] * b->c[2] - a->c[2] * b->c[0];  // 02 --> 1
  c->e[2] = a->c[0] * b->c[3] - a->c[3] * b->c[0];  // 03 --> 2
  c->e[3] = a->c[1] * b->c[2] - a->c[2] * b->c[1];  // 12 --> 3
  c->e[4] = a->c[1] * b->c[3] - a->c[3] * b->c[1];  // 13 --> 4
  c->e[5] = a->c[2] * b->c[3] - a->c[3] * b->c[2];  // 23 --> 5

  // Now 0.5 eps_{abcd} psol^c sol^d
  // Doubled thanks to antisymmetry of sigma and eps
  c->e[0] += a->c[2] * b->c[3] - a->c[3] * b->c[2]; // 0123 - 0132
  c->e[1] += a->c[3] * b->c[1] - a->c[1] * b->c[3]; // 0231 - 0213
  c->e[2] += a->c[1] * b->c[2] - a->c[2] * b->c[1]; // 0312 - 0321
  c->e[3] += a->c[0] * b->c[3] - a->c[3] * b->c[0]; // 1203 - 1230
  c->e[4] += a->c[2] * b->c[0] - a->c[0] * b->c[2]; // 1320 - 1302
  c->e[5] += a->c[0] * b->c[1] - a->c[1] * b->c[0]; // 2301 - 2310
}
// -----------------------------------------------------------------
