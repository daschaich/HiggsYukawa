// -----------------------------------------------------------------
// Matrices times vector operations
// See ../RHMC/setup_perm.c for as_index and sd_index mappings
// c <-- a * b
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/so4.h"

void mult_as_vec(antisym *a, vector *b, vector *c) {
#if (DIMF != 4)
  #error "Assuming DIMF=4!"
#endif
  c->c[0] = a->e[0] * b->c[1];          // 01 --> 0
  c->c[0] += a->e[1] * b->c[2];         // 02 --> 1
  c->c[0] += a->e[2] * b->c[3];         // 03 --> 2

  c->c[1] = -1.0 * a->e[0] * b->c[0];   // 10 --> -0
  c->c[1] += a->e[3] * b->c[2];         // 12 --> 3
  c->c[1] += a->e[4] * b->c[3];         // 13 --> 4

  c->c[2] = -1.0 * a->e[1] * b->c[0];   // 20 --> -1
  c->c[2] -= a->e[3] * b->c[1];         // 21 --> -3
  c->c[2] += a->e[5] * b->c[3];         // 23 --> 5

  c->c[3] = -1.0 * a->e[2] * b->c[0];   // 30 --> -2
  c->c[3] -= a->e[4] * b->c[1];         // 31 --> -4
  c->c[3] -= a->e[5] * b->c[2];         // 32 --> -5
}

void mult_sd_vec(selfdual *a, vector *b, vector *c) {
#if (DIMF != 4)
  #error "Assuming DIMF=4!"
#endif
  c->c[0] = a->e[0] * b->c[1];          // 01 --> 0
  c->c[0] += a->e[1] * b->c[2];         // 02 --> 1
  c->c[0] += a->e[2] * b->c[3];         // 03 --> 2

  c->c[1] = -1.0 * a->e[0] * b->c[0];   // 10 --> -0
  c->c[1] += a->e[2] * b->c[2];         // 12 --> 2
  c->c[1] += a->e[1] * b->c[3];         // 13 --> 1

  c->c[2] = -1.0 * a->e[1] * b->c[0];   // 20 --> -0
  c->c[2] -= a->e[2] * b->c[1];         // 21 --> -2
  c->c[2] += a->e[0] * b->c[3];         // 23 --> 0

  c->c[3] = -1.0 * a->e[2] * b->c[0];   // 30 --> -2
  c->c[3] -= a->e[1] * b->c[1];         // 31 --> -1
  c->c[3] -= a->e[0] * b->c[2];         // 32 --> -0
}
// -----------------------------------------------------------------
