// -----------------------------------------------------------------
// Dirac operator and other helper functions for the action and force
#include "so4_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Copy self-dual scalar fields in the site struct
void scalar_field_copy(field_offset src, field_offset dest) {
  register int i;
  register site *s;

  FORALLSITES(i, s)
    as_copy((antisym *)F_PT(s, src), (antisym *)F_PT(s, dest));
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Matrix--vector operation
// Applies either the operator (sign = 1) or its adjoint (sign = -1)
// Adjoint is simply overall negative sign...
void fermion_op(vector *src, vector *dest, int sign) {
  register int i;
  register site *s;
  int dir, a, b, c, d, p;
  Real tr, halfG = 0.5 * G;
  vector tvec, tvec_dir, tvec_opp;
  msg_tag *tag[6];    // 6 = 2NDIMS

  // Quick sanity check
  if (sign != 1 && sign != -1) {
    node0_printf("Error: incorrect sign in fermion_op: %d\n", sign);
    terminate(1);
  }

  // Start gathers for kinetic term
  for (dir = XUP; dir <= TUP; dir++) {
    tag[dir] = start_gather_field(src, sizeof(vector), dir,
                                  EVENANDODD, gen_pt[dir]);
    tag[OPP_DIR(dir)] = start_gather_field(src, sizeof(vector), OPP_DIR(dir),
                                           EVENANDODD, gen_pt[OPP_DIR(dir)]);
  }

  // Compute scalar term as gathers run
  // Initialize dest = 0.5G * sigma * src
  FORALLSITES(i, s) {
    clearvec(&(dest[i]));
    for (a = 0; a < DIMF; a++) {
      for (b = a + 1; b < DIMF; b++) {
        p = as_index[a][b];
        dest[i].c[a] += s->sigma.e[p] * src[i].c[b];
        dest[i].c[b] -= s->sigma.e[p] * src[i].c[a];
        for (c = 0; c < DIMF; c++) {
          for (d = c + 1; d < DIMF; d++) {
            p = as_index[c][d];
            tr = 0.5 * perm[a][b][c][d];
            dest[i].c[a] += tr * s->sigma.e[p] * src[i].c[b];
            dest[i].c[b] -= tr * s->sigma.e[p] * src[i].c[a];
          }
        }
      }
    }
    scalar_mult_vec(&(dest[i]), halfG, &(dest[i]));
  }

  // Accumulate kinetic term as gathers finish
  for (dir = XUP; dir <= TUP; dir++) {
    wait_gather(tag[dir]);
    wait_gather(tag[OPP_DIR(dir)]);
    FORALLSITES(i, s) {
      // Need to deal with BCs here since no gauge fields
      if (dir == TUP && PBC < 0 && s->t == nt - 1) {
        scalar_mult_vec((vector *)gen_pt[dir][i], -1.0, &tvec_dir);
        vec_copy((vector *)gen_pt[OPP_DIR(dir)][i], &tvec_opp);
      }
      else if (dir == TUP && PBC < 0 && s->t == 0) {
        vec_copy((vector *)gen_pt[dir][i], &tvec_dir);
        scalar_mult_vec((vector *)gen_pt[OPP_DIR(dir)][i], -1.0, &tvec_opp);
      }
      else {
        vec_copy((vector *)gen_pt[dir][i], &tvec_dir);
        vec_copy((vector *)gen_pt[OPP_DIR(dir)][i], &tvec_opp);
      }

      // Add 0.5 * phase(x)[dir] * [psi(x + dir) - psi(x - dir)] to dest
      sub_vec(&tvec_dir, &tvec_opp, &tvec);
      scalar_mult_add_vec(&(dest[i]), &tvec, 0.5 * s->phase[dir], &(dest[i]));
    }
    cleanup_gather(tag[dir]);
    cleanup_gather(tag[OPP_DIR(dir)]);
  }

  // Overall negative sign for adjoint
  if (sign == -1) {
    FORALLSITES(i, s)
      scalar_mult_vec(&(dest[i]), -1.0, &(dest[i]));
  }
}
// -----------------------------------------------------------------
