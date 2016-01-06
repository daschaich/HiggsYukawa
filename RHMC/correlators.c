// -----------------------------------------------------------------
// Measure two- and four-fermion correlators
#include "so4_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Set dest to unit source at given point in given SO(4) index
// Then set src = Ddag dest
// so that (Ddag.D)^(-1).src will give D^(-1).pnt_src
// Return the number of iterations from the inversion
void pnt_src(int *pnt, int index) {
  register int i;
  register site *s;

  // Set dest to unit source at given point in given SO(4) index
  FORALLSITES(i, s)
    clearvec(&(dest[i]));
  if (node_number(pnt[0], pnt[1], pnt[2]) == mynode()) {
    i = node_index(pnt[0], pnt[1], pnt[2]);
    dest[i].c[index] = 1.0;
  }

  // Set src = Ddag dest
  fermion_op(dest, src, MINUS);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Set dest to random gaussian source at all sites
// Then set src = Ddag dest
// so that (Ddag.D)^(-1).src will give D^(-1).vol_src
void vol_src() {
  register int i;
  register site *s;

  // Set dest to random gaussian source at all sites
  FORALLSITES(i, s) {
#ifdef SITERAND
    dest[i].c[0] = gaussian_rand_no(&(s->site_prn));
    dest[i].c[1] = gaussian_rand_no(&(s->site_prn));
    dest[i].c[2] = gaussian_rand_no(&(s->site_prn));
    dest[i].c[3] = gaussian_rand_no(&(s->site_prn));
#else
    dest[i].c[0] = gaussian_rand_no(&node_prn);
    dest[i].c[1] = gaussian_rand_no(&node_prn);
    dest[i].c[2] = gaussian_rand_no(&node_prn);
    dest[i].c[3] = gaussian_rand_no(&node_prn);
#endif
  }

  // Set src = Ddag dest
  fermion_op(dest, src, MINUS);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Measure two- and four-fermion correlators
// Return total number of iterations
int correlators(int *pnt) {
  register int i, j, k, l;
  register site *s;
  int index, iters, tot_iters = 0, sav = Norder;
  Real size_r, four = 0.0, sus = 0.0, ***prop;
  double dtime;
  vector **psim;

  // Make sure pnt stays within lattice volume
  i = pnt[0] % nx;      pnt[0] = i;
  i = pnt[1] % ny;      pnt[1] = i;
  i = pnt[2] % nt;      pnt[2] = i;

  // Allocate structure to hold all DIMF propagators
  prop = malloc(DIMF * sizeof(***prop));
  for (i = 0; i < DIMF; i++) {
    prop[i] = malloc(DIMF * sizeof(vector*));
    for (j = 0; j < DIMF; j++)
      prop[i][j] = malloc(sites_on_node * sizeof(vector));
  }

  // Hack a basic CG out of the multi-mass CG
  Norder = 1;
  psim = malloc(sizeof(**psim));
  psim[0] = malloc(sites_on_node * sizeof(vector));
  shift[0] = 0;

  for (j = 0; j < DIMF; j++) {
    dtime = -dclock();
    pnt_src(pnt, j);
    iters = congrad_multi(src, psim, niter, rsqmin, &size_r);
    dtime += dclock();
    tot_iters += iters;
    node0_printf("Inversion %d of %d took %d iters and %.4g seconds\n",
                 j + 1, DIMF, iters, dtime);

    // Copy psim into f[j][k]
    FORALLSITES(i, s) {
      for (k = 0; k < DIMF; k++)
        prop[j][k][i] = psim[0][i].c[k];
    }

    // Now construct correlators
    // TODO: ...
  }

  // Compute four-fermion condensate
  if (node_number(pnt[0], pnt[1], pnt[2]) == mynode()) {
    index = node_index(pnt[0], pnt[1], pnt[2]);
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
          four += perm[i][j][k][l] * prop[i][j][index] * prop[k][l][index];
          }
        }
      }
    }
  }
  g_doublesum(&four);

  // Compute four-fermion susceptibility
  FORALLSITES(i, s) {
    for (j = 0; j < DIMF; j++) {
      for (k = 0; k < DIMF; k++)
        sus += prop[j][j][i] * prop[k][k][i] - prop[j][k][i] * prop[k][j][i];
    }
  }
  g_doublesum(&sus);

  // Print four-fermion condensate and susceptibility
  node0_printf("FOUR %d %d %d %.6g %.6g %d\n",
               pnt[0], pnt[1], pnt[2], four, sus, tot_iters);

  // Normalize correlators and print results
  // TODO: ...

  // Free structure to hold all DIMF propagators
  for (i = 0; i < DIMF; i++) {
    for (j = 0; j < DIMF; j++)
      free(prop[i][j]);
    free(prop[i]);
  }
  free(prop);

  // Reset multi-mass CG and clean up
  Norder = sav;
  free(psim[0]);
  free(psim);
  return tot_iters;
}
// -----------------------------------------------------------------
