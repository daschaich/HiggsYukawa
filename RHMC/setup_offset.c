// -----------------------------------------------------------------
// label[k] labels the offset of the kth connection between psibar and psi
// by an integer 1 to 3
// The separation of the paths is in New York metric
// It corresponds to our labeling lambda_1 = lambda(1, 0, 0)
// offset[k] is the actual offset of the kth connection
// This program only generates the positive offsets
// The others can be found by using adjoints of these paths
// (e.g., (-1, -1, -1) is the adjoint of (1, 1, 1))

#include "susy_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void cubic_neighbor(int x, int y int t, int *arg, int forw_back,
                    int *xpt, int *ypt, int *tpt) {

  if (forw_back == FORWARDS) {
    *xpt = (x + nx + arg[0]) % nx;
    *ypt = (y + ny + arg[1]) % ny;
    *tpt = (t + nt + arg[3]) % nt;
  }
  else {
    *xpt = (x + nx - arg[0]) % nx;
    *ypt = (y + ny - arg[1]) % ny;
    *tpt = (t + nt - arg[3]) % nt;
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Set up Kogut--Susskind phase factors, sum(nu < mu) {-1^i[nu]}
// combined with boundary conditions
void setup_phases() {
  register int i;
  register site *s;

  FORALLSITES(i, s) {
    s->phase[TUP] = 1.0;
    if ((s->t) % 2 == 1)
      s->phase[XUP] = -1.0;
    else
      s->phase[XUP] = 1.0;

    if ((s->x) % 2 == 1)
      s->phase[YUP] = -s->phase[XUP];
    else
      s->phase[YUP] = s->phase[XUP];

    if (PBC < 0.0) {
      // Antiperiodic boundary conditions in time
      // All t phases for t = nt - 1 time slice get extra minus sign
      if (s->t == nt - 1)
        s->phase[TUP] = -s->phase[TUP];
    }
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void setup_offset() {
  int i, k;

  // Construct the link paths: one in each direction
  for (i = 0; i < NDIMS; i++) {
    for (k = 0; k < NDIMS; k++)
      offset[i][k] = 0;

    offset[i][i] = 1;
  }

#ifdef DEBUG_CHECK
  node0_printf("There are %d distinct paths:\n", NUMLINK);
#endif
  // goffset holds indices of gather_array in ../generic/com_mpi.c
  // The first six elements of gather_array are
  //   XUP, YUP, TUP, TDOWN, YDOWN, XDOWN
  // in that order!
  // In order to use XDOWN = XUP + 1, etc., we make the next six elements
  //   XUP, XDOWN, YUP, YDOWN, TUP, TDOWN
  // Then goffset[0]=6, goffset[1]=8 and goffset[2]=10
  // But we can't use these in EVEN or ODD gathers!
  for (i = 0; i < NUMLINK; i++) {
    goffset[i] = make_gather(cubic_neighbor, offset[i],
                             WANT_INVERSE, NO_EVEN_ODD, SCRAMBLE_PARITY);

#ifdef DEBUG_CHECK
    int dir;
    node0_printf("  %d ahead:", i);
    for (dir = 0; dir < NDIMS; dir++)
      node0_printf(" %d", offset[i][dir]);

    node0_printf(" (offset %d)\n", goffset[i]);
#endif
  }
}
// -----------------------------------------------------------------
