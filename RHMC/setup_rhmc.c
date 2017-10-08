// -----------------------------------------------------------------
// Coefficients for (Mdag M)^(-1 / 2) and (Mdag M)^(1 / 4)
// Note the relative sign between fractional powers!
// For now simply copy in remez output, switch using Norder
#include "so4_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// (17, 17) approximation for spectral range [1e-8, 500]
// Nroot 1 gives errors 1.344046e-05 and 9.402444e-06
void setup_rhmc17() {
  if (Nroot == 1) {
    // awk '/res_MD/{print("amp2["$2"] =",$3";")}' < remez/out.N1.D17
    ampdeg2 = 1.1147371798737100e-02;
    amp2[0] = 5.3722940572685301e-05;
    amp2[1] = 8.8385072480751058e-05;
    amp2[2] = 1.8007340231269867e-04;
    amp2[3] = 3.8794539542209139e-04;
    amp2[4] = 8.4612041500606707e-04;
    amp2[5] = 1.8502137829396362e-03;
    amp2[6] = 4.0480685084372699e-03;
    amp2[7] = 8.8577470397020776e-03;
    amp2[8] = 1.9382492160158894e-02;
    amp2[9] = 4.2413202346146879e-02;
    amp2[10] = 9.2812632978480553e-02;
    amp2[11] = 2.0313318738952396e-01;
    amp2[12] = 4.4491671801114813e-01;
    amp2[13] = 9.7797710754299438e-01;
    amp2[14] = 2.1868748327541017e+00;
    amp2[15] = 5.3125395060272789e+00;
    amp2[16] = 1.9677105020137233e+01;

    // awk '/pole_MD/{print("shift2["$2"] =",$3";")}' < remez/out.N1.D17
    shift2[0] = 1.6130042183256869e-09;
    shift2[1] = 2.1432786036861466e-08;
    shift2[2] = 1.2047405855410658e-07;
    shift2[3] = 5.9557262426526962e-07;
    shift2[4] = 2.8706425840136872e-06;
    shift2[5] = 1.3764277979797591e-05;
    shift2[6] = 6.5925709202639435e-05;
    shift2[7] = 3.1568761745969382e-04;
    shift2[8] = 1.5116111247316733e-03;
    shift2[9] = 7.2380282461060050e-03;
    shift2[10] = 3.4658438313205010e-02;
    shift2[11] = 1.6597498249559348e-01;
    shift2[12] = 7.9522916120383857e-01;
    shift2[13] = 3.8192462317609230e+00;
    shift2[14] = 1.8554468066394989e+01;
    shift2[15] = 9.5370123256473761e+01;
    shift2[16] = 6.6731352332296819e+02;

    // awk '/res_GR/{print("amp4["$2"] =",$3";")}' < remez/out.N1.D17
    ampdeg4 = 9.0494269232624820e+00;
    amp4[0] = -2.8414818252051702e-11;
    amp4[1] = -2.4019260853714093e-10;
    amp4[2] = -1.7395358164726721e-09;
    amp4[3] = -1.2355524981132256e-08;
    amp4[4] = -8.7448459325129812e-08;
    amp4[5] = -6.1848860184139651e-07;
    amp4[6] = -4.3736760924473751e-06;
    amp4[7] = -3.0927752903654763e-05;
    amp4[8] = -2.1870007726193438e-04;
    amp4[9] = -1.5465243973825584e-03;
    amp4[10] = -1.0937097367318850e-02;
    amp4[11] = -7.7379710345894914e-02;
    amp4[12] = -5.4854442303187756e-01;
    amp4[13] = -3.9256485440506355e+00;
    amp4[14] = -2.9402276895258002e+01;
    amp4[15] = -2.7524388588740709e+02;
    amp4[16] = -9.0024603556401817e+03;

    // awk '/pole_GR/{print("shift4["$2"] =",$3";")}' < remez/out.N1.D17
    shift4[0] = 5.6390697956428344e-09;
    shift4[1] = 4.2749125680623754e-08;
    shift4[2] = 2.2251161130347453e-07;
    shift4[3] = 1.0825970001930432e-06;
    shift4[4] = 5.1955188098676434e-06;
    shift4[5] = 2.4863017337093701e-05;
    shift4[6] = 1.1891054670472911e-04;
    shift4[7] = 5.6863429433976604e-04;
    shift4[8] = 2.7191624530224449e-03;
    shift4[9] = 1.3002845195175142e-02;
    shift4[10] = 6.2181019879431115e-02;
    shift4[11] = 2.9741103398251306e-01;
    shift4[12] = 1.4237655020084565e+00;
    shift4[13] = 6.8446425869024274e+00;
    shift4[14] = 3.3581277204051602e+01;
    shift4[15] = 1.8238276883896646e+02;
    shift4[16] = 1.8363089900772914e+03;
  }
  else {
    node0_printf("setup_rhmc: unrecognized Nroot %d\n", Nroot);
    terminate(1);
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Switch between different sets of coefficients
// for (Mdag M)^(-1 / 2) and (Mdag M)^(1 / 4),
// For given spectral range, run with the smallest Norder
// that keeps Remez errors < 2e-5
void setup_rhmc() {
  node0_printf("Using Nroot = %d\n", Nroot);
  node0_printf("RHMC Norder %d for spectral range ", Norder);

  switch(Norder) {
    case 17:
      node0_printf("[1e-8, 500]\n");
      setup_rhmc17();
      break;
    default:
      node0_printf("setup_rhmc: unrecognized Norder %d\n", Norder);
      terminate(1);
  }

  // Test by zeroing out all amp2 and amp4, optionally setting Norder to 1
//  int i;
//  Norder = 1;
//  for (i = 0; i < Norder; i++) {
//    amp2[i] = 0;
//    amp4[i] = 0;
//  }
//  node0_printf("TEST VERSION: internal Norder %d\n", Norder);
//  ampdeg2 = 0;
//  amp2[0] = 1;
//  shift2[0] = 1e-3;
//  ampdeg4 = 1;
//  amp4[0] = 0;
//  shift4[0] = 1;

#ifdef DEBUG_CHECK
  int i;
  node0_printf("RHMC ampdeg2 %e\n", ampdeg2);
  for (i = 0; i < Norder; i++)
    node0_printf("RHMC params %d amp2 %e shift2 %e\n", i, amp2[i], shift2[i]);

  node0_printf("RHMC ampdeg4 %e\n", ampdeg4);
  for (i = 0; i < Norder; i++)
    node0_printf("RHMC params %d amp4 %e shift4 %e\n", i, amp4[i], shift4[i]);
#endif
}
// -----------------------------------------------------------------
