// -----------------------------------------------------------------
// Copy anti-symmetric or self-dual matrices (hardly worth functions)
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/so4.h"

void as_copy(antisym *src, antisym *dest) {
  *dest = *src;
}

void sd_copy(selfdual *src, selfdual *dest) {
  *dest = *src;
}
// -----------------------------------------------------------------
