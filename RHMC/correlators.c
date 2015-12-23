// -----------------------------------------------------------------
// Measure two- and four-fermion correlators
#include "so4_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Set dest to unit source at given point in given SO(4) index
// Then set src = Ddag dest
// so that (Ddag.D)^(-1).src will give D^(-1).pnt_src
// Make sure pnt ends up within lattice volume
// Return the number of iterations from the inversion
void pnt_src(int *pnt, int index) {
  register int i;
  register site *s;

  // Set dest to unit source at given point in given SO(4) index
  FORALLSITES(i, s)
    clearvec(&(dest[i]));
  if (node_number(pnt[0] % nx, pnt[1] % ny, pnt[2] % nt) == mynode()) {
    i = node_index(pnt[0] % nx, pnt[1] % ny, pnt[2] % nt);
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
  int i, iters, tot_iters = 0, sav = Norder;
  Real size_r;
  double dtime;
  vector **psim;

  // Hack a basic CG out of the multi-mass CG
  Norder = 1;
  psim = malloc(sizeof(**psim));
  psim[0] = malloc(sites_on_node * sizeof(vector));
  shift[0] = 0;

  for (i = 0; i < DIMF; i++) {
    dtime = -dclock();
    pnt_src(pnt, i);
    iters = congrad_multi(src, psim, niter, rsqmin, &size_r);
    dtime += dclock();
    tot_iters += iters;
    node0_printf("Inversion %d of %d took %d iters and %.4g secs\n",
                 i + 1, DIMF, iters, dtime);

    // Now construct correlators
    // TODO: ...
  }

  // Normalize correlators and print results
  // TODO: ...

  // Reset multi-mass CG and clean up
  Norder = sav;
  free(psim[0]);
  free(psim);
  return tot_iters;
}
// -----------------------------------------------------------------
