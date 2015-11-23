// -----------------------------------------------------------------
// Update lattice
// Omelyan integrator with Urbach, Jansen, Schindler, Wenger multiscale
// (CPC 174:87, 2006)

// Begin at "integral" time, with H and U evaluated at the same time
// For the final accept/reject, we already have a good solution to the CG
// The last update was of the momenta

// Uncomment to print out debugging messages
//#define UPDATE_DEBUG
#include "so4_includes.h"
#ifdef HAVE_IEEEFP_H
#include <ieeefp.h>         // For "finite"
#endif
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void ranmom() {
#if (DIMF != 4)
  #error "Assuming DIMF=4!"
#endif
  register int i;
  register site *s;

  FORALLSITES(i, s) {
#ifdef SITERAND
    mom[i].e[0] = gaussian_rand_no(&(s->site_prn));
    mom[i].e[1] = gaussian_rand_no(&(s->site_prn));
    mom[i].e[2] = gaussian_rand_no(&(s->site_prn));
    mom[i].e[3] = gaussian_rand_no(&(s->site_prn));
    mom[i].e[4] = gaussian_rand_no(&(s->site_prn));
    mom[i].e[5] = gaussian_rand_no(&(s->site_prn));
#else
    mom[i].e[0] = gaussian_rand_no(&(s->node_prn));
    mom[i].e[1] = gaussian_rand_no(&(s->node_prn));
    mom[i].e[2] = gaussian_rand_no(&(s->node_prn));
    mom[i].e[3] = gaussian_rand_no(&(s->node_prn));
    mom[i].e[4] = gaussian_rand_no(&(s->node_prn));
    mom[i].e[5] = gaussian_rand_no(&(s->node_prn));
#endif
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Omelyan version; ``dirty'' speeded-up version
double update_scalar(Real eps) {
  register int i;
  register site *s;
  int nsw = nsteps[1], isw;
  Real halfstep = 0.5 * eps;
  double norm;

#ifdef UPDATE_DEBUG
  node0_printf("scalar %d steps %.4g dt\n", nsw, eps);
#endif
  norm = scalar_force(eps * LAMBDA);
  for (isw = 1; isw <= nsw; isw++) {
    FORALLSITES(i, s)
      scalar_mult_add_as(&(s->sigma), &(mom[i]), halfstep, &(s->sigma));
    norm += scalar_force(eps * LAMBDA_MID);

    FORALLSITES(i, s)
      scalar_mult_add_as(&(s->sigma), &(mom[i]), halfstep, &(s->sigma));

    if (isw < nsw)
      norm += scalar_force(eps * TWO_LAMBDA);
    else
      norm += scalar_force(eps * LAMBDA);
  }
  return (norm / nsw);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
int update_step(double *fnorm, double *snorm, vector **src, vector ***psim) {
  int iters = 0, i_multi0, n;
  Real final_rsq, f_eps, s_eps, tr;

  f_eps = traj_length / (Real)nsteps[0];
  s_eps = f_eps / (Real)(2.0 * nsteps[1]);

  // Fermion update
  for (n = 0; n < Nroot; n++) {
    // CG called before update_step
    tr = fermion_force(f_eps * LAMBDA, src[n], psim[n]);
    fnorm[n] += tr;
    if (tr > max_ff[n])
      max_ff[n] = tr;
  }

  for (i_multi0 = 1; i_multi0 <= nsteps[0]; i_multi0++) {
    // Scalar update
    tr = update_scalar(s_eps);
    *snorm += tr;
    if (tr > max_sf)
      max_sf = tr;

    // Fermion update
    for (n = 0; n < Nroot; n++) {
      // Do conjugate gradient to get (Mdag M)^(-1 / 4) chi
      iters += congrad_multi(src[n], psim[n], niter, rsqmin, &final_rsq);
      tr = fermion_force(f_eps * LAMBDA_MID, src[n], psim[n]);
      fnorm[n] += tr;
      if (tr > max_ff[n])
        max_ff[n] = tr;
    }

    // Scalar update
    tr = update_scalar(s_eps);
    *snorm += tr;
    if (tr > max_sf)
      max_sf = tr;

    // Fermion update
    for (n = 0; n < Nroot; n++) {
      // Do conjugate gradient to get (Mdag M)^(-1 / 4) chi
      iters += congrad_multi(src[n], psim[n], niter, rsqmin, &final_rsq);

      if (i_multi0 < nsteps[0])
        tr = fermion_force(f_eps * TWO_LAMBDA, src[n], psim[n]);
      else
        tr = fermion_force(f_eps * LAMBDA, src[n], psim[n]);
      fnorm[n] += tr;
      if (tr > max_ff[n])
        max_ff[n] = tr;
    }
  }
  return iters;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
int update() {
  int i, n, iters = 0;
  Real final_rsq;
  double startaction, endaction, change;
  vector **src = malloc(Nroot * sizeof(**src));
  vector ***psim = malloc(Nroot * sizeof(***psim));

  for (n = 0; n < Nroot; n++) {
    src[n] = malloc(sites_on_node * sizeof(vector));
    psim[n] = malloc(Norder * sizeof(vector*));
    for (i = 0; i < Norder; i++)
      psim[n][i] = malloc(sites_on_node * sizeof(vector));
  }

  // Refresh the momenta
  ranmom();

  // Set up the fermion variables
  // Compute g and src = (Mdag M)^(1 / 8) g
  for (n = 0; n < Nroot; n++)
    iters += grsource(src[n]);

  // Do a CG to get psim,
  // rational approximation to (Mdag M)^(-1 / 4) src = (Mdag M)^(-1 / 8) g
  for (i = 0; i < Norder; i++)
    shift[i] = shift4[i];
#ifdef UPDATE_DEBUG
  node0_printf("Calling CG in update_o -- original action\n");
#endif
  // congrad_multi initializes psim
  for (n = 0; n < Nroot; n++)
    iters += congrad_multi(src[n], psim[n], niter, rsqmin, &final_rsq);

  // Find initial action
  startaction = action(src, psim);
  snorm = 0.0;
  max_sf = 0.0;
  for (n = 0; n < Nroot; n++) {
    fnorm[n] = 0.0;
    max_ff[n] = 0.0;
  }

#ifdef HMC_ALGORITHM
  Real xrandom;   // For accept/reject test
  // Copy scalar field to old_sigma
  scalar_field_copy(F_OFFSET(sigma), F_OFFSET(old_sigma));
#endif
  // Do microcanonical updating
  iters += update_step(fnorm, &snorm, src, psim);

  // Find ending action
  // Reuse data from update_step, don't need CG to get (Mdag M)^(-1) chi
  // If the final step were a gauge update, CG would be necessary
  endaction = action(src, psim);
  change = endaction - startaction;
#ifdef HMC_ALGORITHM
  // Reject configurations giving overflow
#ifndef HAVE_IEEEFP_H
  if (fabs((double)change) > 1e20) {
#else
  if (!finite((double)change)) {
#endif
    node0_printf("WARNING: Correcting Apparent Overflow: Delta S = %.4g\n",
                 change);
    change = 1.0e20;
  }

  // Decide whether to accept, if not, copy old link field back
  // Careful -- must generate only one random number for whole lattice
  if (this_node == 0)
    xrandom = myrand(&node_prn);
  broadcast_float(&xrandom);
  if (exp(-change) < (double)xrandom) {
    if (traj_length > 0.0) {
      scalar_field_copy(F_OFFSET(old_sigma), F_OFFSET(sigma));
    }
    node0_printf("REJECT: delta S = %.4g start S = %.12g end S = %.12g\n",
                 change, startaction, endaction);
  }
  else {
    node0_printf("ACCEPT: delta S = %.4g start S = %.12g end S = %.12g\n",
                 change, startaction, endaction);
  }
#else
  // Only print check if not doing HMC
  node0_printf("CHECK: delta S = %.4g\n", (double)(change));
#endif // ifdef HMC

  for (n = 0; n < Nroot; n++) {
    free(src[n]);
    for (i = 0; i < Norder; i++)
      free(psim[n][i]);
    free(psim[n]);
  }
  free(src);
  free(psim);

  if (traj_length > 0) {
    node0_printf("MONITOR_FORCE_SCALAR   %.4g %.4g\n",
                 snorm / (double)(2 * nsteps[0]), max_sf);
    for (n = 0; n < Nroot; n++) {
      node0_printf("MONITOR_FORCE_FERMION%d %.4g %.4g\n",
                   n, fnorm[n] / (double)(2 * nsteps[0]), max_ff[n]);
    }
    return iters;
  }
  else
    return -99;
}
// -----------------------------------------------------------------
