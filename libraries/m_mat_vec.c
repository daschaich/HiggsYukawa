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

void mult_sd_vec(antisym *a, vector *b, vector *c) {
  // Start with psi^a sigma_{ab} psi^b terms (checks DIMF)
  mult_as_vec(a, b, c);

  // Now add 0.5 psi^a eps_{abcd} sigma_{cd} psi^b terms
  // For any given c & d we pick up 2*0.5 thanks to two negative signs,
  // one from eps_{abdc} = -eps_{abcd}, the other from sigma_{dc} = -sigma_{cd}
  // Just need sign of given eps_{abcd}
  c->c[0] += a->e[5] * b->c[1];         // 0123 23 --> 5
  c->c[0] -= a->e[4] * b->c[2];         // 0213 13 --> -4
  c->c[0] += a->e[2] * b->c[3];         // 0312 12 --> 3

  c->c[1] -= a->e[5] * b->c[0];         // 1023 23 --> -5
  c->c[1] += a->e[2] * b->c[2];         // 1203 03 --> 2
  c->c[1] -= a->e[1] * b->c[3];         // 1302 02 --> -1

  c->c[2] += a->e[4] * b->c[0];         // 2013 13 --> 4
  c->c[2] -= a->e[2] * b->c[1];         // 2103 03 --> -2
  c->c[2] += a->e[0] * b->c[3];         // 2301 01 --> 0

  c->c[3] -= a->e[3] * b->c[0];         // 3012 12 --> -3
  c->c[3] += a->e[1] * b->c[1];         // 3102 02 --> 1
  c->c[3] -= a->e[0] * b->c[2];         // 3201 01 --> -0
}
// -----------------------------------------------------------------
