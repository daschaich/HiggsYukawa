// -----------------------------------------------------------------
// Measure four-fermion correlator, susceptibility and correlators
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
  if (node_number(pnt[0], pnt[1], pnt[2], pnt[3]) == mynode()) {
    i = node_index(pnt[0], pnt[1], pnt[2], pnt[3]);
    dest[i].c[index] = 1.0;
  }

  // Set src = Ddag dest
  fermion_op(dest, src, MINUS);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// For some reason this doesn't work in the libraries
Real Z2_rand_no(double_prn *prn_pt) {
  if (myrand(prn_pt) > 0.5)
    return 1.0;
  else
    return -1.0;
}

// Set dest to random Z2 source at all sites
// Then set src = Ddag dest
// so that (Ddag.D)^(-1).src will give D^(-1).vol_src
void vol_src() {
#if (DIMF != 4)
  #error "Assuming DIMF=4!"
#endif
  register int i;
  register site *s;

  FORALLSITES(i, s) {
#ifdef SITERAND
    dest[i].c[0] = Z2_rand_no(&(s->site_prn));
    dest[i].c[1] = Z2_rand_no(&(s->site_prn));
    dest[i].c[2] = Z2_rand_no(&(s->site_prn));
    dest[i].c[3] = Z2_rand_no(&(s->site_prn));
#else
    dest[i].c[0] = Z2_rand_no(&node_prn);
    dest[i].c[1] = Z2_rand_no(&node_prn);
    dest[i].c[2] = Z2_rand_no(&node_prn);
    dest[i].c[3] = Z2_rand_no(&node_prn);
#endif
  }
  fermion_op(dest, src, MINUS);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Measure four-fermion condensate and its susceptibility
// Use Z2 stochastic sources
// Return total number of iterations
int condensates() {
  register int i, ii;
  register site *s, *ss;
  int iters, tot_iters = 0, sav = Norder;
  int a, b, c, d, j, Nsrc = DIMF;
  Real size_r, norm = 1.0 / (Real)Nsrc, ***prop;
  double bilin = 0.0, four = 0.0, sus = 0.0, dtime;
  vector **psim;

  // Allocate and initialize stochastic propagator
  prop = malloc(DIMF * sizeof(***prop));
  for (a = 0; a < DIMF; a++) {
    prop[a] = malloc(DIMF * sizeof(Real*));
    for (b = 0; b < DIMF; b++) {
      prop[a][b] = malloc(sites_on_node * sizeof(Real));
      FORALLSITES(i, s) {
        prop[a][b][i] = 0.0;
      }
    }
  }

  // Hack a basic CG out of the multi-mass CG
  Norder = 1;
  psim = malloc(sizeof(**psim));
  psim[0] = malloc(sites_on_node * sizeof(vector));
  shift[0] = 0;

  // Construct stochastic propagator
  //   D_{ab}^{-1}(x) = (1/N) sum_N psi_a(x) dest_b(x)
  // where dest are random Z2 sources
  // Hit each dest with Mdag to get src_j, invert to get D_{kj}^{-1} dest_j
  for (j = 0; j < Nsrc; j++) {
    dtime = -dclock();
    vol_src();
    iters = congrad_multi(src, psim, niter, rsqmin, &size_r);
    dtime += dclock();
    tot_iters += iters;
    node0_printf("Inversion %d of %d took %d iters and %.4g seconds\n",
                 j + 1, DIMF, iters, dtime);

    // Copy psim into f[k][j]
    for (a = 0; a < DIMF; a++) {
      for (b = 0; b < DIMF; b++) {
        FORALLSITES(i, s) {
          prop[a][b][i] += psim[0][i].c[a] * dest[i].c[b];
//          FORALLSITES(ii, ss) {
//            sus += psim[0][i].c[a] * dest[ii].c[b] * psim[0][i].c[b] * dest[ii].c[a];
//            sus -= psim[0][i].c[a] * dest[ii].c[a] * psim[0][i].c[b] * dest[ii].c[b];
//          }
        }
      }
    }
  }

  // Normalize stochastic propagator by norm = 1 / Nsrc
  for (a = 0; a < DIMF; a++) {
    for (b = 0; b < DIMF; b++) {
      FORALLSITES(i, s)
        prop[a][b][i] *= norm;
    }
  }

  // Four-fermion condensate is just eps_{abcd} D_{ab}^{-1} D_{cd}^{-1}
  // Other permutations (-D_ac D_bd + D_ad D_bc) give identical contribution
  // TODO: Susceptibility is
  //   psi_a(x) psi_b(x) psi_a(y) psi_b(y)
  //     = D_{ab}^{-1} D_{ab}^{-1}
  //     - D_{aa}^{-1}(x - y) D_{bb}^{-1}(x - y)
  //     + D_{ab}^{-1}(x - y) D_{ba}^{-1}(x - y)
  for (a = 0; a < DIMF; a++) {
    for (b = 0; b < DIMF; b++) {
      if (b == a)
        continue;
      FORALLSITES(i, s) {
        bilin += prop[a][b][i];
//        sus_abab += prop[a][b][i] * prop[a][b][i];
//        FORALLSITES(ii, ss) {     // TODO: Need gathers...
//          sus_aabb -= prop[a][a][i] * prop[b][b][ii];
//          sus_abba += prop[a][b][i] * prop[b][a][ii];
//        }
      }
      for (c = 0; c < DIMF; c++) {
        if (c == b || c == a)
          continue;
        for (d = 0; d < DIMF; d++) {
          if (d == c || d == b || d == a)
            continue;
          FORALLSITES(i, s)
            four += perm[a][b][c][d] * prop[a][b][i] * prop[c][d][i];
        }
      }
    }
  }
  g_doublesum(&bilin);
//  g_doublesum(&sus_abab);
//  g_doublesum(&sus_aabb);
//  g_doublesum(&sus_abba);
  g_doublesum(&four);
  bilin /= (double)volume;
//  sus_abba /= (double)volume;
//  sus_aabb /= (double)volume;
//  sus_abab /= (double)volume;
  four /= (double)volume;

  // Print condensates and susceptibility
  node0_printf("STOCH BILIN %.6g %d\n", bilin, tot_iters);
  node0_printf("STOCH FOUR %.6g %d\n", four, tot_iters);
//  sus = sus_abab + sus_aabb + sus_abba;
//  node0_printf("STOCH SUS %.6g %.6g %.6g %.6g %d\n",
//               sus_abab, sus_aabb, sus_abba, sus, tot_iters);

  // Free structure to hold all DIMF propagators
  for (a = 0; a < DIMF; a++) {
    for (b = 0; b < DIMF; b++)
      free(prop[a][b]);
    free(prop[a]);
  }
  free(prop);

  // Reset multi-mass CG and clean up
  Norder = sav;
  free(psim[0]);
  free(psim);
  return tot_iters;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Measure two- and four-fermion correlators
// Return total number of iterations
int correlators(int *pnt) {
  register int i;
  register site *s;
  int a, b, c, d, iters, tot_iters = 0, sav = Norder;
  Real size_r, ***prop;
  double bilin = 0.0, four = 0.0, dtime;
  double sus_abba = 0.0, sus_aabb = 0.0, sus_abab = 0.0, sus;
  vector **psim;

  // Make sure pnt stays within lattice volume
  i = pnt[0] % nx;      pnt[0] = i;
  i = pnt[1] % ny;      pnt[1] = i;
  i = pnt[2] % nz;      pnt[2] = i;
  i = pnt[3] % nt;      pnt[3] = i;

  // Allocate structure to hold all DIMF propagators
  prop = malloc(DIMF * sizeof(***prop));
  for (a = 0; a < DIMF; a++) {
    prop[a] = malloc(DIMF * sizeof(Real*));
    for (b = 0; b < DIMF; b++)
      prop[a][b] = malloc(sites_on_node * sizeof(Real));
  }

  // Hack a basic CG out of the multi-mass CG
  Norder = 1;
  psim = malloc(sizeof(**psim));
  psim[0] = malloc(sites_on_node * sizeof(vector));
  shift[0] = 0;

  for (a = 0; a < DIMF; a++) {
    dtime = -dclock();
    pnt_src(pnt, a);
    iters = congrad_multi(src, psim, niter, rsqmin, &size_r);
    dtime += dclock();
    tot_iters += iters;
    node0_printf("Inversion %d of %d took %d iters and %.4g seconds\n",
                 a + 1, DIMF, iters, dtime);

    // Copy psim into f[j][k]
    FORALLSITES(i, s) {
      for (b = 0; b < DIMF; b++)
        prop[b][a][i] = psim[0][i].c[b];
    }

    // Now construct correlators
    // TODO: ...
  }

  // Compute four-fermion condensate
  if (node_number(pnt[0], pnt[1], pnt[2], pnt[3]) == mynode()) {
    i = node_index(pnt[0], pnt[1], pnt[2], pnt[3]);
    for (a = 0; a < DIMF; a++) {
      for (b = 0; b < DIMF; b++) {
        if (b == a)
          continue;
        bilin += prop[a][b][i];
        sus_abab += prop[a][b][i] * prop[a][b][i];
        for (c = 0; c < DIMF; c++) {
          if (c == b || c == a)
            continue;
          for (d = 0; d < DIMF; d++) {
            if (d == c || d == b || d == a)
              continue;
            four += perm[a][b][c][d] * prop[a][b][i] * prop[c][d][i];
          }
        }
      }
    }
  }
  g_doublesum(&bilin);
  g_doublesum(&sus_abab);
  g_doublesum(&four);

  // Compute four-fermion susceptibility
  for (a = 0; a < DIMF; a++) {
    for (b = 0; b < DIMF; b++) {
      FORALLSITES(i, s) {
        sus_aabb -= prop[a][a][i] * prop[b][b][i];
        sus_abba += prop[a][b][i] * prop[b][a][i];
      }
    }
  }
  g_doublesum(&sus_aabb);
  g_doublesum(&sus_abba);

  // Print condensates and susceptibility
  node0_printf("PNT BILIN %d %d %d %d %.6g %d\n",
               pnt[0], pnt[1], pnt[2], pnt[3], bilin, tot_iters);
  node0_printf("PNT FOUR %d %d %d %d %.6g %d\n",
               pnt[0], pnt[1], pnt[2], pnt[3], four, tot_iters);
  sus = sus_abab + sus_aabb + sus_abba;
  node0_printf("PNT SUS %d %d %d %d %.6g %.6g %.6g %.6g %d\n",
               pnt[0], pnt[1], pnt[2], pnt[3],
               sus_abab, sus_aabb, sus_abba, sus, tot_iters);

  // Normalize correlators and print results
  // TODO: ...

  // Free structure to hold all DIMF propagators
  for (a = 0; a < DIMF; a++) {
    for (b = 0; b < DIMF; b++)
      free(prop[a][b]);
    free(prop[a]);
  }
  free(prop);

  // Reset multi-mass CG and clean up
  Norder = sav;
  free(psim[0]);
  free(psim);
  return tot_iters;
}
// -----------------------------------------------------------------
