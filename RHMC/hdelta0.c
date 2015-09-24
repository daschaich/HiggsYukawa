// -----------------------------------------------------------------
// Four-fermion matrix--vector operation, dest = D^2 src
// Use tempvec for temporary storage
#include "so4_includes.h"

void hdelta0(vector *src, vector *dest) {
  fermion_op(src, tempvec, PLUS);
  fermion_op(tempvec, dest, MINUS);
}
// -----------------------------------------------------------------
