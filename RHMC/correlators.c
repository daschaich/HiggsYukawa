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
  register int i;//, ii;
  register site *s;//, *ss;
  int iters, tot_iters = 0, sav = Norder;
  int a, b, c, d, j;
  Real size_r, norm = 1.0 / (Real)DIMF, ***prop, ***prop2;
  double four = 0.0, sus = 0.0, dtime;
  vector **psim;
  double bilin[DIMF][DIMF],bilin2[DIMF][DIMF];
 
  for(a=0;a<DIMF;a++){
  for(b=0;b<DIMF;b++){
  bilin[a][b]=0.0;
  bilin2[a][b]=0.0;}}

  // Allocate and initialize stochastic propagator
  prop = malloc(DIMF * sizeof(***prop));
  prop2 = malloc(DIMF * sizeof(***prop2));
  for (a = 0; a < DIMF; a++) {
    prop[a] = malloc(DIMF * sizeof(Real*));
    prop2[a] = malloc(DIMF * sizeof(Real*));
    for (b = 0; b < DIMF; b++) {
      prop[a][b] = malloc(sites_on_node * sizeof(Real));
      prop2[a][b] = malloc(sites_on_node * sizeof(Real));
      FORALLSITES(i, s){
        prop[a][b][i] = 0.0;
        prop2[a][b][i] = 0.0;}
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
  for (j = 0; j < 32; j++) {
    dtime = -dclock();
    vol_src();
    iters = congrad_multi(src, psim, niter, rsqmin, &size_r);
    dtime += dclock();
    tot_iters += iters;
    node0_printf("1st inversion %d of %d took %d iters and %.4g seconds\n",
                 j + 1, DIMF, iters, dtime);

    // Copy psim into f[k][j]
    for (a = 0; a < DIMF; a++) {
      for (b = 0; b < DIMF; b++) {
        FORALLSITES(i, s)
          prop[a][b][i] += psim[0][i].c[a] * dest[i].c[b];
      }
    }

    dtime = -dclock();
    vol_src();
    iters = congrad_multi(src, psim, niter, rsqmin, &size_r);
    dtime += dclock();
    tot_iters += iters;
    node0_printf("2nd inversion %d of %d took %d iters and %.4g seconds\n",
                 j+1,Nsrc,iters,dtime);
   
    for (a = 0; a < DIMF; a++) {
      for (b = 0; b < DIMF; b++) {
        FORALLSITES(i, s)
          prop2[a][b][i] += psim[0][i].c[a] * dest[i].c[b];
      }
    }
  }

  // Normalize stochastic propagator by norm = 1 / DIMF
  for (a = 0; a < DIMF; a++) {
    for (b = 0; b < DIMF; b++) {
      FORALLSITES(i, s){
        prop[a][b][i] *= norm;
        prop2[a][b][i] *= norm;}
    }
  }

  // Four-fermion condensate is just eps_{abcd} D_{ab}^{-1} D_{cd}^{-1}
  // Other permutations (-D_ac D_bd + D_ad D_bc) give identical contribution
  for (a = 0; a < DIMF; a++) {
    for (b = 0; b < DIMF; b++) {
      if (b == a)
        continue;
        FORALLSITES(i, s){
          bilin[a][b] += prop[a][b][i];
          bilin2[a][b] += prop2[a][b][i];}

      for (c = 0; c < DIMF; c++) {
        if (c == b || c == a)
          continue;
        for (d = 0; d < DIMF; d++) {
          if (d == c || d == b || d == a)
            continue;
          FORALLSITES(i, s)
            four += perm[a][b][c][d] * prop[a][b][i] * prop2[c][d][i];
        }
      }
    }
  }
   for (a = 0; a < DIMF; a++){
  for (b = 0; b < DIMF; b++){
  g_doublesum(&bilin[a][b]);
  g_doublesum(&bilin2[a][b]);}}

  g_doublesum(&four);
  for (a = 0; a < DIMF; a++){
  for (b = 0; b < DIMF; b++){
  bilin[a][b] /= (double)volume;
  bilin2[a][b] /= (double)volume;}}
  
  sus=0.0;
  four /= (double)volume;
  for (a = 0; a < DIMF; a++){
  for (b = 0; b < DIMF; b++){
  sus += bilin[a][b]*bilin2[a][b];}}
  sus *=(double)volume;
  // Print condensates and susceptibility
  node0_printf("STOCH BILIN %.6g %d\n", 0.5*(bilin[0][1]+bilin2[0][1]), tot_iters);
  node0_printf("STOCH FOUR %.6g %d\n", four, tot_iters);

  // Free structure to hold all DIMF propagators
  for (a = 0; a < DIMF; a++) {
    for (b = 0; b < DIMF; b++){
      free(prop[a][b]);
      free(prop2[a][b]);}
    free(prop[a]);
    free(prop2[a]);
  }
  free(prop);
  free(prop2);

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
  int a, b, c, d, iters, tot_iters = 0, sav = Norder, dir;
  int L[NDIMS] = {nx, ny, nz, nt};
  Real size_r;
  double bilin = 0.0, four = 0.0, dtime;
  double sus_abba = 0.0, sus_aabb = 0.0, sus_abab = 0.0, sus;
  double one_link[NDIMS] = {0.0, 0.0, 0.0, 0.0};
  vector **psim;
  matrix *tm;
  msg_tag *tag[NDIMS];

  // Make sure pnt stays within lattice volume
  for (dir = XUP; dir <= TUP; dir++) {
    i = pnt[dir] % L[dir];
    pnt[dir] = i;
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
        prop[i].e[b][a] = psim[0][i].c[b];
    }
  }

  // Gather prop from -mu for one-link condensates
  for (dir = XUP; dir <= TUP; dir++) {
    if (L[dir] > 1) {
      tag[dir] = start_gather_field(prop, sizeof(matrix), OPP_DIR(dir),
                                    EVENANDODD, gen_pt[dir]);
    }
  }

  // Compute four-fermion condensate
  if (node_number(pnt[0], pnt[1], pnt[2], pnt[3]) == mynode()) {
    i = node_index(pnt[0], pnt[1], pnt[2], pnt[3]);
    for (a = 0; a < DIMF; a++) {
      for (b = 0; b < DIMF; b++) {
        if (b == a)
          continue;
        if ((a==0) && (b==1))
          bilin += prop[i].e[a][b];
        for (c = 0; c < DIMF; c++) {
          if (c == b || c == a)
            continue;
          for (d = 0; d < DIMF; d++) {
            if (d == c || d == b || d == a)
              continue;
            four += perm[a][b][c][d] * prop[i].e[a][b] * prop[i].e[c][d];
          }
        }
      }
    }
  }
  g_doublesum(&bilin);
  g_doublesum(&four);

  // Compute four-fermion susceptibility
  for (a = 0; a < DIMF; a++) {
    for (b = 0; b < DIMF; b++) {
      FORALLSITES(i, s) {
        sus_aabb -= prop[i].e[a][a] * prop[i].e[b][b];
        sus_abba += prop[i].e[a][b] * prop[i].e[b][a];
      }
    }
  }
  g_doublesum(&sus_aabb);
  g_doublesum(&sus_abba);

  // Compute one-link condensates
  for (dir = XUP; dir <= TUP; dir++) {
    if (L[dir] <= 1)          // Don't re-compute on-site bilinear
      continue;
    wait_gather(tag[dir]);
    if (node_number(pnt[0], pnt[1], pnt[2], pnt[3]) == mynode()) {
      i = node_index(pnt[0], pnt[1], pnt[2], pnt[3]);
      tm = (matrix *)(gen_pt[dir][i]);
      for (a = 0; a < DIMF; a++) {
        for (b = a + 1; b < DIMF; b++){
         if (pnt[dir] % 2 == 0)
            one_link[dir] += lattice[i].phase[dir] * tm->e[a][b];
          else
            one_link[dir] -= lattice[i].phase[dir] * tm->e[a][b];
      }
     }
    }
    g_doublesum(&(one_link[dir]));
    cleanup_gather(tag[dir]);
  }

  // Print condensates and susceptibility
  node0_printf("PNT BILIN %d %d %d %d %.6g %d\n",
               pnt[0], pnt[1], pnt[2], pnt[3], bilin, tot_iters);
  node0_printf("PNT FOUR %d %d %d %d %.6g %d\n",
               pnt[0], pnt[1], pnt[2], pnt[3], four, tot_iters);
  sus = sus_aabb + sus_abba;
  node0_printf("PNT SUS %d %d %d %d %.6g %.6g %.6g %d\n",
               pnt[0], pnt[1], pnt[2], pnt[3],
               sus_aabb, sus_abba, sus, tot_iters);
  node0_printf("PNT ONELINK %d %d %d %d %.6g %.6g %.6g %.6g %d\n",
               pnt[0], pnt[1], pnt[2], pnt[3],
               one_link[0], one_link[1], one_link[2], one_link[3], tot_iters);

  // Reset multi-mass CG and clean up
  Norder = sav;
  free(psim[0]);
  free(psim);
  return tot_iters;
}
// -----------------------------------------------------------------
