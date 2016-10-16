// -----------------------------------------------------------------
// Clear an irrep matrix
#include "../include/config.h"
#include "../include/so4.h"

void clear_as(antisym *m) {
#if (DIMF != 4)
  #error "Assuming DIMF=4!"
#endif
  m->e[0] = 0.0;
  m->e[1] = 0.0;
  m->e[2] = 0.0;
  m->e[3] = 0.0;
  m->e[4] = 0.0;
  m->e[5] = 0.0;
}

void clear_sd(selfdual *m) {
#if (DIMF != 4)
  #error "Assuming DIMF=4!"
#endif
  m->e[0] = 0.0;
  m->e[1] = 0.0;
  m->e[2] = 0.0;
}
// -----------------------------------------------------------------
