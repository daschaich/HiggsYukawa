// -----------------------------------------------------------------
// Construct a gaussian random vector R, return src = (Ddag.D)^{1 / 4} R
// Need to invert despite the positive power, since it is fractional
#include "so4_includes.h"

// Return the number of iterations from the inversion
int grsource(vector *src) {
#if (DIMF != 4)
  #error "Assuming DIMF=4!"
#endif
  register int i, j;
  register site *s;
  int avs_iters;
  Real size_r;
  vector **psim = malloc(Norder * sizeof(**psim));

  // Allocate psim (will be zeroed in congrad_multi)
  for (i = 0; i < Norder; i++)
    psim[i] = malloc(sites_on_node * sizeof(vector));

  // Begin with pure gaussian random numbers
  FORALLSITES(i, s) {
#ifdef SITERAND
    src[i].c[0] = gaussian_rand_no(&(s->site_prn));
    src[i].c[1] = gaussian_rand_no(&(s->site_prn));
    src[i].c[2] = gaussian_rand_no(&(s->site_prn));
    src[i].c[3] = gaussian_rand_no(&(s->site_prn));
#else
    src[i].c[0] = gaussian_rand_no(&node_prn);
    src[i].c[1] = gaussian_rand_no(&node_prn);
    src[i].c[2] = gaussian_rand_no(&node_prn);
    src[i].c[3] = gaussian_rand_no(&node_prn);
#endif
  }

#ifdef DEBUG_CHECK
  double source_norm = 0.0;
  FORALLSITES(i, s) {
//    if (i != 0)
//      clearvec(&(src[i]));
//    if (i == 0)
//      dumpvec(&src[i]);
    source_norm += (double)magsq_vec(&(src[i]));
  }
  g_doublesum(&source_norm);
  node0_printf("source_norm in grsource %.4g\n", source_norm);

//  fermion_op(src, psim[0], PLUS);
//  fermion_op(psim[0], psim[1], MINUS);
//  printf("\n\n TEST 1\n");
//  FORALLSITES(i, s) {
//    source_norm = (double)magsq_vec(&(psim[1][i]));
//    if (source_norm * source_norm > 0) {
//      printf("%d %d %d %.4g\n", s->x, s->y, s->t, source_norm);
//      dumpas(&(s->sigma));
//    }
//  }

//  fermion_op(src, psim[0], PLUS);
//  printf("\n\n TEST 2\n");
//  FORALLSITES(i, s) {
//    source_norm = (double)magsq_vec(&(psim[0][i]));
//    if (source_norm * source_norm > 0) {
//      printf("%d %d %d %.4g\n", s->x, s->y, s->t, source_norm);
//      dumpas(&(s->sigma));
//    }
//  }
#endif

  // We now compute (Mdag M)^{1 / 4}.src
  for (i = 0; i < Norder; i++)
    shift[i] = shift4[i];

  avs_iters = congrad_multi(src, psim, niter, rsqmin, &size_r);
#ifdef DEBUG_CHECK
  node0_printf("Iters for source %d\n", avs_iters);
#endif
  // Reconstruct (Mdag M)^{1 / 4}.src from multi-mass CG solution psim
  FORALLSITES(i, s) {
    scalar_mult_vec(&(src[i]), ampdeg4, &(src[i]));
    for (j = 0; j < Norder; j++)
      scalar_mult_sum_vec(&(psim[j][i]), amp4[j], &(src[i]));
  }

  // Correct normalization for real pseudofermions (factor of 1/sqrt(2))
  FORALLSITES(i, s)
    scalar_mult_vec(&(src[i]), 0.70710678118, &(src[i]));

  for (i = 0; i < Norder; i++)
    free(psim[i]);
  free(psim);
  return avs_iters;
}
// -----------------------------------------------------------------
