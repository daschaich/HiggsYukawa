// -----------------------------------------------------------------
// Set up generator matrices and epsilon^{ijklm}
#include "so4_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Set up four-index totally anti-symmetric tensor
// Initialize swap to avoid optimization dependence!!!
Real order(int i, int j, int k, int l) {
  int seq[4] = {i, j, k, l};
  int swap = 1, tmp, p, permutation = 1;
  while (swap > 0) {
    swap = 0;
    for (p = 0; p < 3; p++) {
      if (seq[p] > seq[p + 1]) {
        tmp = seq[p];
        seq[p] = seq[p + 1];
        seq[p + 1] = tmp;
        swap++;
        permutation *= -1;
      }
    }
  }
  return (Real)permutation;
}

// Set up translation of (mu, nu) to linear index of anti-symmetric matrix
void setup_as_index() {
#if (DIMF != 4)
  #error "Assuming DIMF=4!"
#endif
  int mu;
  for (mu = 0; mu < DIMF; mu++)
    as_index[mu][mu] = -1;
  for (mu = 1; mu < DIMF; mu++) {    // 0, 1 and 2
    as_index[0][mu] = mu - 1;
    as_index[mu][0] = mu - 1;
  }
  for (mu = 2; mu < DIMF; mu++) {    // 3 and 4
    as_index[1][mu] = mu + 1;
    as_index[mu][1] = mu + 1;
  }
  as_index[2][3] = 5;
  as_index[3][2] = 5;
}

// Set up translation of (mu, nu) to linear index of self-dual matrix
// Pair up {01, 23}, {02, 13} and {03, 12}
void setup_sd_index() {
#if (DIMF != 4)
  #error "Assuming DIMF=4!"
#endif
  int mu;
  for (mu = 0; mu < DIMF; mu++)
    sd_index[mu][mu] = -1;
  for (mu = 1; mu < DIMF; mu++) {    // 0, 1 and 2
    sd_index[0][mu] = mu - 1;
    sd_index[mu][0] = mu - 1;
  }
  for (mu = 2; mu < DIMF; mu++) {    // 2 and 1
    sd_index[1][mu] = 4 - mu;
    sd_index[mu][1] = 4 - mu;
  }
  sd_index[2][3] = 0;
  sd_index[3][2] = 0;
}

void epsilon() {
  int i, j, k, l;
  setup_as_index();
  setup_sd_index();
  for (i = 0; i < DIMF; i++) {
    for (j = 0; j < DIMF; j++) {
      for (k = 0; k < DIMF; k++) {
        for (l = 0; l < DIMF; l++)
          perm[i][j][k][l] = 0;
      }
    }
  }

  for (i = 0; i < DIMF; i++) {
    for (j = 0; j < DIMF; j++) {
      if (j == i)
        continue;
      for (k = 0; k < DIMF; k++) {
        if (k == j || k == i)
          continue;
        for (l = 0; l < DIMF; l++) {
          if (l == k || l == j || l == i)
            continue;
          perm[i][j][k][l] = order(i, j, k, l);
#ifdef DEBUG_CHECK
          if (perm[i][j][k][l] * perm[i][j][k][l] > 1.0e-4)
            node0_printf("PERM %d%d%d%d = %.4g\n",
                         i, j, k, l, perm[i][j][k][l]);
#endif
        }
      }
    }
  }
  return;
}
// -----------------------------------------------------------------
