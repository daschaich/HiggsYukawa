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
  int dir, a, b, c, d;
  Real tr, halfG = 0.5 * G, m_ov_G, vev[DIMF][DIMF];
  vector tvec, tvec_dir, tvec_opp;
  msg_tag *tag[2 * NDIMS];

  // Quick sanity check
  if (sign != 1 && sign != -1) {
    node0_printf("Error: incorrect sign in fermion_op: %d\n", sign);
    terminate(1);
  }

  // Ignore site_mass if G = 0 to avoid dividing by zero
  // Could be made more robust, but unlikely to matter
  if (G == 0.0)
    m_ov_G = 0.0;
  else
    m_ov_G = 2.0 * site_mass / G;
   for(a=0;a<DIMF;a++){
	for(b=0;b<DIMF;b++){
         vev[a][b]=0.0;}}
    vev[0][1]=m_ov_G;
    vev[2][3]=m_ov_G;
    vev[1][0]=-m_ov_G;
    vev[3][2]=-m_ov_G;

  // Start gathers for kinetic term
  for (dir = XUP; dir <= TUP; dir++) {
    tag[dir] = start_gather_field(src, sizeof(vector), dir,
                                  EVENANDODD, gen_pt[dir]);
    tag[OPP_DIR(dir)] = start_gather_field(src, sizeof(vector), OPP_DIR(dir),
                                           EVENANDODD, gen_pt[OPP_DIR(dir)]);
  }

  // Compute scalar term as gathers run
  // Initialize dest = 0.5G * (sigma + 2m / G) * src
  // Add SO(4)-breaking 'site mass' term with same structure as sigma
  FORALLSITES(i, s) {
    clearvec(&(dest[i]));
    for (a = 0; a < DIMF; a++) {
      for (b = a + 1; b < DIMF; b++) {
        tr = s->sigma.e[as_index[a][b]] + vev[a][b];
        for (c = 0; c < DIMF; c++) {
          for (d = c + 1; d < DIMF; d++)
            tr += perm[a][b][c][d] * (s->sigma.e[as_index[c][d]] + vev[c][d]);
        }   // No half since not double-counting
        dest[i].c[a] += tr * src[i].c[b];
        dest[i].c[b] -= tr * src[i].c[a];
      }
    }
    scalar_mult_vec(&(dest[i]), halfG, &(dest[i]));
  }

  // Accumulate kinetic term as gathers finish
  for (dir = XUP; dir <= TUP; dir++) {
    wait_gather(tag[dir]);
    wait_gather(tag[OPP_DIR(dir)]);
    FORALLSITES(i, s) {
      // Deal with BCs here
      vec_copy((vector *)gen_pt[dir][i], &tvec_dir);
      vec_copy((vector *)gen_pt[OPP_DIR(dir)][i], &tvec_opp);
      if (dir == TUP && PBC < 0 && s->t == nt - 1)
        scalar_mult_vec(&tvec_dir, -1.0, &tvec_dir);
      else if (dir == TUP && PBC < 0 && s->t == 0)
        scalar_mult_vec(&tvec_opp, -1.0, &tvec_opp);

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



// -----------------------------------------------------------------
// Squared four-fermion matrix--vector operation
//   dest = D^2 src
// Use tempvec for temporary storage
void DSq(vector *src, vector *dest) {
  fermion_op(src, tempvec, PLUS);
  fermion_op(tempvec, dest, MINUS);
}
// -----------------------------------------------------------------
