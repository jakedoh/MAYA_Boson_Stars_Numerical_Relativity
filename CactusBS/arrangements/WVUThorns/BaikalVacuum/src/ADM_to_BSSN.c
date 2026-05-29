
#include <math.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

void BaikalVacuum_ADM_to_BSSN(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;

    CCTK_REAL *alphaSphorCartGF = alp;
    CCTK_REAL *BSphorCartU0GF = dtbetax;
    CCTK_REAL *BSphorCartU1GF = dtbetay;
    CCTK_REAL *BSphorCartU2GF = dtbetaz;
    CCTK_REAL *KSphorCartDD00GF = kxx;
    CCTK_REAL *KSphorCartDD01GF = kxy;
    CCTK_REAL *KSphorCartDD02GF = kxz;
    CCTK_REAL *KSphorCartDD11GF = kyy;
    CCTK_REAL *KSphorCartDD12GF = kyz;
    CCTK_REAL *KSphorCartDD22GF = kzz;
    CCTK_REAL *betaSphorCartU0GF = betax;
    CCTK_REAL *betaSphorCartU1GF = betay;
    CCTK_REAL *betaSphorCartU2GF = betaz;
    CCTK_REAL *gammaSphorCartDD00GF = gxx;
    CCTK_REAL *gammaSphorCartDD01GF = gxy;
    CCTK_REAL *gammaSphorCartDD02GF = gxz;
    CCTK_REAL *gammaSphorCartDD11GF = gyy;
    CCTK_REAL *gammaSphorCartDD12GF = gyz;
    CCTK_REAL *gammaSphorCartDD22GF = gzz;
    #pragma omp parallel for
    for (int i2 = 0; i2 < cctk_lsh[2]; i2++) {
        for (int i1 = 0; i1 < cctk_lsh[1]; i1++) {
            for (int i0 = 0; i0 < cctk_lsh[0]; i0++) {
                   /*
                    * NRPy+ Finite Difference Code Generation, Step 1 of 2: Read from main memory and compute finite difference stencils:
                    */
                   const double alphaSphorCart = alphaSphorCartGF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
                   const double betaSphorCartU0 = betaSphorCartU0GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
                   const double betaSphorCartU1 = betaSphorCartU1GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
                   const double betaSphorCartU2 = betaSphorCartU2GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
                   const double BSphorCartU0 = BSphorCartU0GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
                   const double BSphorCartU1 = BSphorCartU1GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
                   const double BSphorCartU2 = BSphorCartU2GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
                   const double gammaSphorCartDD00 = gammaSphorCartDD00GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
                   const double gammaSphorCartDD01 = gammaSphorCartDD01GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
                   const double gammaSphorCartDD02 = gammaSphorCartDD02GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
                   const double gammaSphorCartDD11 = gammaSphorCartDD11GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
                   const double gammaSphorCartDD12 = gammaSphorCartDD12GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
                   const double gammaSphorCartDD22 = gammaSphorCartDD22GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
                   const double KSphorCartDD00 = KSphorCartDD00GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
                   const double KSphorCartDD01 = KSphorCartDD01GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
                   const double KSphorCartDD02 = KSphorCartDD02GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
                   const double KSphorCartDD11 = KSphorCartDD11GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
                   const double KSphorCartDD12 = KSphorCartDD12GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
                   const double KSphorCartDD22 = KSphorCartDD22GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
                   /*
                    * NRPy+ Finite Difference Code Generation, Step 2 of 2: Evaluate SymPy expressions and write to main memory:
                    */
                   const double FDPart3_1 = gammaSphorCartDD00*gammaSphorCartDD11*gammaSphorCartDD22;
                   const double FDPart3_3 = 2*gammaSphorCartDD01*gammaSphorCartDD02*gammaSphorCartDD12;
                   const double FDPart3_5 = gammaSphorCartDD00*((gammaSphorCartDD12)*(gammaSphorCartDD12));
                   const double FDPart3_7 = ((gammaSphorCartDD01)*(gammaSphorCartDD01))*gammaSphorCartDD22;
                   const double FDPart3_9 = ((gammaSphorCartDD02)*(gammaSphorCartDD02))*gammaSphorCartDD11;
                   const double FDPart3_10 = FDPart3_1 + FDPart3_3 - FDPart3_5 - FDPart3_7 - FDPart3_9;
                   const double FDPart3_11 = (1.0/(FDPart3_10));
                   const double FDPart3_12 = cbrt(FDPart3_11);
                   const double FDPart3_13 = 2*FDPart3_11;
                   const double FDPart3_14 = FDPart3_11*KSphorCartDD00*(gammaSphorCartDD11*gammaSphorCartDD22 - ((gammaSphorCartDD12)*(gammaSphorCartDD12))) + FDPart3_11*KSphorCartDD11*(gammaSphorCartDD00*gammaSphorCartDD22 - ((gammaSphorCartDD02)*(gammaSphorCartDD02))) + FDPart3_11*KSphorCartDD22*(gammaSphorCartDD00*gammaSphorCartDD11 - ((gammaSphorCartDD01)*(gammaSphorCartDD01))) + FDPart3_13*KSphorCartDD01*(-gammaSphorCartDD01*gammaSphorCartDD22 + gammaSphorCartDD02*gammaSphorCartDD12) + FDPart3_13*KSphorCartDD02*(gammaSphorCartDD01*gammaSphorCartDD12 - gammaSphorCartDD02*gammaSphorCartDD11) + FDPart3_13*KSphorCartDD12*(-gammaSphorCartDD00*gammaSphorCartDD12 + gammaSphorCartDD01*gammaSphorCartDD02);
                   const double FDPart3_15 = (1.0/3.0)*FDPart3_14;
                   hDD00GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = FDPart3_12*gammaSphorCartDD00 - 1;
                   hDD01GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = FDPart3_12*gammaSphorCartDD01;
                   hDD02GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = FDPart3_12*gammaSphorCartDD02;
                   hDD11GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = FDPart3_12*gammaSphorCartDD11 - 1;
                   hDD12GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = FDPart3_12*gammaSphorCartDD12;
                   hDD22GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = FDPart3_12*gammaSphorCartDD22 - 1;
                   aDD00GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = FDPart3_12*(-FDPart3_15*gammaSphorCartDD00 + KSphorCartDD00);
                   aDD01GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = FDPart3_12*(-FDPart3_15*gammaSphorCartDD01 + KSphorCartDD01);
                   aDD02GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = FDPart3_12*(-FDPart3_15*gammaSphorCartDD02 + KSphorCartDD02);
                   aDD11GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = FDPart3_12*(-FDPart3_15*gammaSphorCartDD11 + KSphorCartDD11);
                   aDD12GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = FDPart3_12*(-FDPart3_15*gammaSphorCartDD12 + KSphorCartDD12);
                   aDD22GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = FDPart3_12*(-FDPart3_15*gammaSphorCartDD22 + KSphorCartDD22);
                   trKGF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = FDPart3_14;
                   vetU0GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = betaSphorCartU0;
                   vetU1GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = betaSphorCartU1;
                   vetU2GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = betaSphorCartU2;
                   betU0GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = BSphorCartU0;
                   betU1GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = BSphorCartU1;
                   betU2GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = BSphorCartU2;
                   alphaGF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = alphaSphorCart;
                   cfGF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = pow(FDPart3_10/(FDPart3_1*FDPart3_11 + FDPart3_11*FDPart3_3 - FDPart3_11*FDPart3_5 - FDPart3_11*FDPart3_7 - FDPart3_11*FDPart3_9), -1.0/6.0);
                
            } // END LOOP: for (int i0 = 0; i0 < cctk_lsh[0]; i0++)
        } // END LOOP: for (int i1 = 0; i1 < cctk_lsh[1]; i1++)
    } // END LOOP: for (int i2 = 0; i2 < cctk_lsh[2]; i2++)

    const CCTK_REAL invdx0 = 1.0/CCTK_DELTA_SPACE(0);
    const CCTK_REAL invdx1 = 1.0/CCTK_DELTA_SPACE(1);
    const CCTK_REAL invdx2 = 1.0/CCTK_DELTA_SPACE(2);
    if(FD_order == 4) {
        #include "ADM_to_BSSN__compute_lambdaU_FD_order_4.h"
    }
    if(FD_order == 6) {
        #include "ADM_to_BSSN__compute_lambdaU_FD_order_6.h"
    }
    if(FD_order == 8) {
        #include "ADM_to_BSSN__compute_lambdaU_FD_order_8.h"
    }

    ExtrapolateGammas(cctkGH,lambdaU0GF);
    ExtrapolateGammas(cctkGH,lambdaU1GF);
    ExtrapolateGammas(cctkGH,lambdaU2GF);
}
