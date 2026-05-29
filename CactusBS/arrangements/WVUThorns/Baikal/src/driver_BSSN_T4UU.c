
#include <math.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "SIMD/SIMD_intrinsics.h"

void Baikal_driver_BSSN_T4UU(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;

    const CCTK_REAL *restrict T4DD00GF = eTtt;
    const CCTK_REAL *restrict T4DD01GF = eTtx;
    const CCTK_REAL *restrict T4DD02GF = eTty;
    const CCTK_REAL *restrict T4DD03GF = eTtz;
    const CCTK_REAL *restrict T4DD11GF = eTxx;
    const CCTK_REAL *restrict T4DD12GF = eTxy;
    const CCTK_REAL *restrict T4DD13GF = eTxz;
    const CCTK_REAL *restrict T4DD22GF = eTyy;
    const CCTK_REAL *restrict T4DD23GF = eTyz;
    const CCTK_REAL *restrict T4DD33GF = eTzz;
#pragma omp parallel for
for (int i2 = 0; i2 < cctk_lsh[2]; i2++) {
    for (int i1 = 0; i1 < cctk_lsh[1]; i1++) {
        for (int i0 = 0; i0 < cctk_lsh[0]; i0 += SIMD_width) {
                  /*
                   * NRPy+ Finite Difference Code Generation, Step 1 of 2: Read from main memory and compute finite difference stencils:
                   */
                  const REAL_SIMD_ARRAY hDD00 = ReadSIMD(&hDD00GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)]);
                  const REAL_SIMD_ARRAY hDD01 = ReadSIMD(&hDD01GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)]);
                  const REAL_SIMD_ARRAY hDD02 = ReadSIMD(&hDD02GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)]);
                  const REAL_SIMD_ARRAY hDD11 = ReadSIMD(&hDD11GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)]);
                  const REAL_SIMD_ARRAY hDD12 = ReadSIMD(&hDD12GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)]);
                  const REAL_SIMD_ARRAY hDD22 = ReadSIMD(&hDD22GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)]);
                  const REAL_SIMD_ARRAY vetU0 = ReadSIMD(&vetU0GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)]);
                  const REAL_SIMD_ARRAY vetU1 = ReadSIMD(&vetU1GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)]);
                  const REAL_SIMD_ARRAY vetU2 = ReadSIMD(&vetU2GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)]);
                  const REAL_SIMD_ARRAY cf = ReadSIMD(&cfGF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)]);
                  const REAL_SIMD_ARRAY alpha = ReadSIMD(&alphaGF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)]);
                  const REAL_SIMD_ARRAY T4DD00 = ReadSIMD(&T4DD00GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)]);
                  const REAL_SIMD_ARRAY T4DD01 = ReadSIMD(&T4DD01GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)]);
                  const REAL_SIMD_ARRAY T4DD02 = ReadSIMD(&T4DD02GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)]);
                  const REAL_SIMD_ARRAY T4DD03 = ReadSIMD(&T4DD03GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)]);
                  const REAL_SIMD_ARRAY T4DD11 = ReadSIMD(&T4DD11GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)]);
                  const REAL_SIMD_ARRAY T4DD12 = ReadSIMD(&T4DD12GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)]);
                  const REAL_SIMD_ARRAY T4DD13 = ReadSIMD(&T4DD13GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)]);
                  const REAL_SIMD_ARRAY T4DD22 = ReadSIMD(&T4DD22GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)]);
                  const REAL_SIMD_ARRAY T4DD23 = ReadSIMD(&T4DD23GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)]);
                  const REAL_SIMD_ARRAY T4DD33 = ReadSIMD(&T4DD33GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)]);
                  /*
                   * NRPy+ Finite Difference Code Generation, Step 2 of 2: Evaluate SymPy expressions and write to main memory:
                   */
                  const double tmpFDPart3_Integer_1 = 1.0;
                  const REAL_SIMD_ARRAY FDPart3_Integer_1 = ConstSIMD(tmpFDPart3_Integer_1);
                  
                  const double tmpFDPart3_Integer_2 = 2.0;
                  const REAL_SIMD_ARRAY FDPart3_Integer_2 = ConstSIMD(tmpFDPart3_Integer_2);
                  
                  const double tmpFDPart3_NegativeOne_ = -1.0;
                  const REAL_SIMD_ARRAY FDPart3_NegativeOne_ = ConstSIMD(tmpFDPart3_NegativeOne_);
                  
                  const REAL_SIMD_ARRAY FDPart3_0 = DivSIMD(FDPart3_Integer_1, MulSIMD(MulSIMD(MulSIMD(alpha, alpha), alpha), alpha));
                  const REAL_SIMD_ARRAY FDPart3_1 = MulSIMD(FDPart3_0, T4DD00);
                  const REAL_SIMD_ARRAY FDPart3_2 = MulSIMD(vetU0, vetU0);
                  const REAL_SIMD_ARRAY FDPart3_4 = MulSIMD(vetU1, vetU1);
                  const REAL_SIMD_ARRAY FDPart3_6 = MulSIMD(vetU2, vetU2);
                  const REAL_SIMD_ARRAY FDPart3_8 = MulSIMD(FDPart3_0, vetU0);
                  const REAL_SIMD_ARRAY FDPart3_10 = MulSIMD(vetU1, vetU2);
                  const REAL_SIMD_ARRAY FDPart3_14 = MulSIMD(T4DD02, vetU1);
                  const REAL_SIMD_ARRAY FDPart3_15 = MulSIMD(T4DD03, vetU2);
                  const REAL_SIMD_ARRAY FDPart3_16 = DivSIMD(FDPart3_Integer_1, MulSIMD(alpha, alpha));
                  const REAL_SIMD_ARRAY FDPart3_18 = MulSIMD(FDPart3_16, MulSIMD(FDPart3_NegativeOne_, vetU0));
                  const REAL_SIMD_ARRAY FDPart3_19 = DivSIMD(FDPart3_Integer_1, MulSIMD(MulSIMD(MulSIMD(cf, cf), cf), cf));
                  const REAL_SIMD_ARRAY FDPart3_21 = AddSIMD(FDPart3_Integer_1, hDD22);
                  const REAL_SIMD_ARRAY FDPart3_23 = DivSIMD(FDPart3_Integer_1, MulSIMD(MulSIMD(MulSIMD(MulSIMD(MulSIMD(cf, cf), cf), cf), cf), cf));
                  const REAL_SIMD_ARRAY FDPart3_27 = AddSIMD(FDPart3_Integer_1, hDD11);
                  const REAL_SIMD_ARRAY FDPart3_30 = AddSIMD(FDPart3_Integer_1, hDD00);
                  const REAL_SIMD_ARRAY FDPart3_32 = DivSIMD(FDPart3_Integer_1, NegFusedMulAddSIMD(MulSIMD(FDPart3_23, FDPart3_27), MulSIMD(hDD02, hDD02), NegFusedMulAddSIMD(MulSIMD(FDPart3_23, FDPart3_30), MulSIMD(hDD12, hDD12), FusedMulAddSIMD(hDD01, MulSIMD(MulSIMD(FDPart3_23, FDPart3_Integer_2), MulSIMD(hDD02, hDD12)), FusedMulSubSIMD(MulSIMD(FDPart3_21, FDPart3_23), MulSIMD(FDPart3_27, FDPart3_30), MulSIMD(FDPart3_21, MulSIMD(FDPart3_23, MulSIMD(hDD01, hDD01))))))));
                  const REAL_SIMD_ARRAY FDPart3_33 = FusedMulAddSIMD(FDPart3_18, vetU1, MulSIMD(FDPart3_32, FusedMulSubSIMD(FDPart3_19, MulSIMD(hDD02, hDD12), MulSIMD(FDPart3_19, MulSIMD(FDPart3_21, hDD01)))));
                  const REAL_SIMD_ARRAY FDPart3_34 = MulSIMD(FDPart3_16, MulSIMD(FDPart3_NegativeOne_, T4DD02));
                  const REAL_SIMD_ARRAY FDPart3_35 = FusedMulAddSIMD(FDPart3_18, vetU2, MulSIMD(FDPart3_32, FusedMulSubSIMD(FDPart3_19, MulSIMD(hDD01, hDD12), MulSIMD(FDPart3_19, MulSIMD(FDPart3_27, hDD02)))));
                  const REAL_SIMD_ARRAY FDPart3_36 = MulSIMD(FDPart3_16, MulSIMD(FDPart3_NegativeOne_, T4DD03));
                  const REAL_SIMD_ARRAY FDPart3_37 = MulSIMD(FDPart3_33, T4DD12);
                  const REAL_SIMD_ARRAY FDPart3_38 = MulSIMD(FDPart3_16, vetU0);
                  const REAL_SIMD_ARRAY FDPart3_39 = MulSIMD(FDPart3_35, T4DD13);
                  const REAL_SIMD_ARRAY FDPart3_40 = MulSIMD(FDPart3_33, T4DD22);
                  const REAL_SIMD_ARRAY FDPart3_41 = MulSIMD(FDPart3_16, vetU1);
                  const REAL_SIMD_ARRAY FDPart3_42 = MulSIMD(FDPart3_35, T4DD23);
                  const REAL_SIMD_ARRAY FDPart3_43 = MulSIMD(FDPart3_33, T4DD23);
                  const REAL_SIMD_ARRAY FDPart3_44 = MulSIMD(FDPart3_16, vetU2);
                  const REAL_SIMD_ARRAY FDPart3_45 = MulSIMD(FDPart3_35, T4DD33);
                  const REAL_SIMD_ARRAY FDPart3_47 = FusedMulSubSIMD(FDPart3_32, FusedMulSubSIMD(FDPart3_19, MulSIMD(FDPart3_21, FDPart3_27), MulSIMD(FDPart3_19, MulSIMD(hDD12, hDD12))), MulSIMD(FDPart3_16, FDPart3_2));
                  const REAL_SIMD_ARRAY FDPart3_48 = MulSIMD(FDPart3_16, MulSIMD(FDPart3_NegativeOne_, T4DD01));
                  const REAL_SIMD_ARRAY FDPart3_49 = MulSIMD(FDPart3_47, T4DD11);
                  const REAL_SIMD_ARRAY FDPart3_50 = MulSIMD(FDPart3_47, T4DD12);
                  const REAL_SIMD_ARRAY FDPart3_51 = MulSIMD(FDPart3_47, T4DD13);
                  const REAL_SIMD_ARRAY FDPart3_53 = FusedMulSubSIMD(FDPart3_32, FusedMulSubSIMD(FDPart3_19, MulSIMD(hDD01, hDD02), MulSIMD(FDPart3_19, MulSIMD(FDPart3_30, hDD12))), MulSIMD(FDPart3_10, FDPart3_16));
                  const REAL_SIMD_ARRAY FDPart3_55 = MulSIMD(FDPart3_38, FDPart3_53);
                  const REAL_SIMD_ARRAY FDPart3_56 = MulSIMD(FDPart3_33, T4DD13);
                  const REAL_SIMD_ARRAY FDPart3_57 = MulSIMD(FDPart3_41, FDPart3_53);
                  const REAL_SIMD_ARRAY FDPart3_59 = FusedMulSubSIMD(FDPart3_32, FusedMulSubSIMD(FDPart3_19, MulSIMD(FDPart3_21, FDPart3_30), MulSIMD(FDPart3_19, MulSIMD(hDD02, hDD02))), MulSIMD(FDPart3_16, FDPart3_4));
                  const REAL_SIMD_ARRAY FDPart3_60 = MulSIMD(FDPart3_59, T4DD12);
                  const REAL_SIMD_ARRAY FDPart3_62 = MulSIMD(FDPart3_59, T4DD23);
                  const REAL_SIMD_ARRAY FDPart3_63 = MulSIMD(FDPart3_35, FDPart3_41);
                  const REAL_SIMD_ARRAY FDPart3_65 = FusedMulSubSIMD(FDPart3_32, FusedMulSubSIMD(FDPart3_19, MulSIMD(FDPart3_27, FDPart3_30), MulSIMD(FDPart3_19, MulSIMD(hDD01, hDD01))), MulSIMD(FDPart3_16, FDPart3_6));
                  const REAL_SIMD_ARRAY FDPart3_68 = MulSIMD(FDPart3_33, FDPart3_33);
                  const REAL_SIMD_ARRAY FDPart3_69 = MulSIMD(FDPart3_35, FDPart3_35);
                  const REAL_SIMD_ARRAY FDPart3_70 = MulSIMD(FDPart3_33, FDPart3_Integer_2);
                  const REAL_SIMD_ARRAY FDPart3_71 = MulSIMD(FDPart3_35, FDPart3_Integer_2);
                  const REAL_SIMD_ARRAY FDPart3_72 = MulSIMD(FDPart3_38, T4DD01);
                  const REAL_SIMD_ARRAY FDPart3_73 = MulSIMD(FDPart3_33, FDPart3_35);
                  const REAL_SIMD_ARRAY FDPart3_75 = MulSIMD(FDPart3_14, FDPart3_16);
                  const REAL_SIMD_ARRAY FDPart3_79 = MulSIMD(FDPart3_15, FDPart3_16);
                  const REAL_SIMD_ARRAY FDPart3_80 = MulSIMD(FDPart3_53, FDPart3_53);
                  const REAL_SIMD_ARRAY FDPart3_81 = MulSIMD(FDPart3_53, FDPart3_Integer_2);
                  const REAL_SIMD_ARRAY FDPart3_82 = MulSIMD(FDPart3_65, FDPart3_Integer_2);
                  const REAL_SIMD_ARRAY __RHS_exp_0 = SubSIMD(FusedMulAddSIMD(MulSIMD(FDPart3_8, FDPart3_Integer_2), MulSIMD(T4DD12, vetU1), FusedMulAddSIMD(MulSIMD(FDPart3_8, FDPart3_Integer_2), MulSIMD(T4DD13, vetU2), FusedMulAddSIMD(FDPart3_0, MulSIMD(FDPart3_6, T4DD33), FusedMulAddSIMD(MulSIMD(FDPart3_0, FDPart3_10), MulSIMD(FDPart3_Integer_2, T4DD23), FusedMulAddSIMD(FDPart3_0, MulSIMD(FDPart3_2, T4DD11), FusedMulAddSIMD(FDPart3_0, MulSIMD(FDPart3_4, T4DD22), SubSIMD(NegFusedMulAddSIMD(MulSIMD(FDPart3_0, FDPart3_Integer_2), MulSIMD(T4DD01, vetU0), FDPart3_1), MulSIMD(MulSIMD(FDPart3_0, FDPart3_Integer_2), MulSIMD(T4DD03, vetU2))))))))), MulSIMD(MulSIMD(FDPart3_0, FDPart3_Integer_2), MulSIMD(T4DD02, vetU1)));
                  const REAL_SIMD_ARRAY __RHS_exp_1 = FusedMulAddSIMD(FDPart3_47, FDPart3_48, FusedMulAddSIMD(FDPart3_44, FDPart3_45, FusedMulAddSIMD(FDPart3_44, FDPart3_51, FusedMulAddSIMD(FDPart3_41, FDPart3_50, FusedMulAddSIMD(FDPart3_43, FDPart3_44, FusedMulAddSIMD(FDPart3_40, FDPart3_41, FusedMulAddSIMD(FDPart3_41, FDPart3_42, FusedMulAddSIMD(FDPart3_38, FDPart3_39, FusedMulAddSIMD(FDPart3_38, FDPart3_49, FusedMulAddSIMD(FDPart3_35, FDPart3_36, FusedMulAddSIMD(FDPart3_37, FDPart3_38, FusedMulAddSIMD(FDPart3_15, FDPart3_8, FusedMulAddSIMD(FDPart3_33, FDPart3_34, FusedMulAddSIMD(FDPart3_0, MulSIMD(FDPart3_2, T4DD01), FusedMulSubSIMD(FDPart3_14, FDPart3_8, MulSIMD(FDPart3_1, vetU0))))))))))))))));
                  const REAL_SIMD_ARRAY __RHS_exp_2 = FusedMulAddSIMD(FDPart3_44, MulSIMD(FDPart3_53, T4DD33), FusedMulAddSIMD(FDPart3_33, MulSIMD(FDPart3_38, T4DD11), FusedMulAddSIMD(FDPart3_41, MulSIMD(FDPart3_59, T4DD22), FusedMulAddSIMD(FDPart3_0, MulSIMD(FDPart3_15, vetU1), FusedMulAddSIMD(FDPart3_0, MulSIMD(FDPart3_4, T4DD02), FusedMulAddSIMD(FDPart3_55, T4DD13, FusedMulAddSIMD(FDPart3_57, T4DD23, FusedMulAddSIMD(FDPart3_44, FDPart3_56, FusedMulAddSIMD(FDPart3_44, FDPart3_62, FusedMulAddSIMD(FDPart3_37, FDPart3_41, FusedMulAddSIMD(FDPart3_38, FDPart3_60, FusedMulAddSIMD(FDPart3_34, FDPart3_59, FusedMulAddSIMD(FDPart3_36, FDPart3_53, FusedMulAddSIMD(FDPart3_8, MulSIMD(T4DD01, vetU1), FusedMulSubSIMD(FDPart3_33, FDPart3_48, MulSIMD(FDPart3_1, vetU1))))))))))))))));
                  const REAL_SIMD_ARRAY __RHS_exp_3 = FusedMulAddSIMD(FDPart3_44, MulSIMD(FDPart3_65, T4DD33), FusedMulAddSIMD(FDPart3_41, MulSIMD(FDPart3_65, T4DD23), FusedMulAddSIMD(FDPart3_44, MulSIMD(FDPart3_53, T4DD23), FusedMulAddSIMD(FDPart3_35, MulSIMD(FDPart3_38, T4DD11), FusedMulAddSIMD(FDPart3_38, MulSIMD(FDPart3_65, T4DD13), FusedMulAddSIMD(FDPart3_0, MulSIMD(FDPart3_14, vetU2), FusedMulAddSIMD(FDPart3_0, MulSIMD(FDPart3_6, T4DD03), FusedMulAddSIMD(FDPart3_57, T4DD22, FusedMulAddSIMD(FDPart3_63, T4DD12, FusedMulAddSIMD(FDPart3_39, FDPart3_44, FusedMulAddSIMD(FDPart3_55, T4DD12, FusedMulAddSIMD(FDPart3_35, FDPart3_48, FusedMulAddSIMD(FDPart3_36, FDPart3_65, FusedMulAddSIMD(FDPart3_8, MulSIMD(T4DD01, vetU2), FusedMulSubSIMD(FDPart3_34, FDPart3_53, MulSIMD(FDPart3_1, vetU2))))))))))))))));
                  const REAL_SIMD_ARRAY __RHS_exp_4 = FusedMulAddSIMD(FDPart3_38, MulSIMD(FDPart3_70, T4DD02), FusedMulAddSIMD(FDPart3_38, MulSIMD(FDPart3_71, T4DD03), FusedMulAddSIMD(FDPart3_69, T4DD33, FusedMulAddSIMD(T4DD11, MulSIMD(FDPart3_47, FDPart3_47), FusedMulAddSIMD(FDPart3_51, FDPart3_71, FusedMulAddSIMD(FDPart3_68, T4DD22, FusedMulAddSIMD(FDPart3_47, MulSIMD(FDPart3_72, FDPart3_Integer_2), FusedMulAddSIMD(FDPart3_73, MulSIMD(FDPart3_Integer_2, T4DD23), FusedMulAddSIMD(FDPart3_1, FDPart3_2, MulSIMD(FDPart3_50, FDPart3_70))))))))));
                  const REAL_SIMD_ARRAY __RHS_exp_5 = FusedMulAddSIMD(FDPart3_73, T4DD13, FusedMulAddSIMD(FDPart3_1, MulSIMD(vetU0, vetU1), FusedMulAddSIMD(FDPart3_63, T4DD03, FusedMulAddSIMD(FDPart3_68, T4DD12, FusedMulAddSIMD(FDPart3_51, FDPart3_53, FusedMulAddSIMD(FDPart3_55, T4DD03, FusedMulAddSIMD(FDPart3_45, FDPart3_53, FusedMulAddSIMD(FDPart3_50, FDPart3_59, FusedMulAddSIMD(FDPart3_42, FDPart3_59, FusedMulAddSIMD(FDPart3_43, FDPart3_53, FusedMulAddSIMD(FDPart3_33, FDPart3_75, FusedMulAddSIMD(FDPart3_40, FDPart3_59, FusedMulAddSIMD(FDPart3_38, MulSIMD(FDPart3_59, T4DD02), FusedMulAddSIMD(FDPart3_41, MulSIMD(FDPart3_47, T4DD01), FusedMulAddSIMD(FDPart3_33, FDPart3_49, MulSIMD(FDPart3_33, FDPart3_72))))))))))))))));
                  const REAL_SIMD_ARRAY __RHS_exp_6 = FusedMulAddSIMD(FDPart3_1, MulSIMD(vetU0, vetU2), FusedMulAddSIMD(FDPart3_33, MulSIMD(FDPart3_44, T4DD02), FusedMulAddSIMD(FDPart3_69, T4DD13, FusedMulAddSIMD(FDPart3_73, T4DD12, FusedMulAddSIMD(FDPart3_51, FDPart3_65, FusedMulAddSIMD(FDPart3_55, T4DD02, FusedMulAddSIMD(FDPart3_45, FDPart3_65, FusedMulAddSIMD(FDPart3_50, FDPart3_53, FusedMulAddSIMD(FDPart3_42, FDPart3_53, FusedMulAddSIMD(FDPart3_43, FDPart3_65, FusedMulAddSIMD(FDPart3_35, FDPart3_79, FusedMulAddSIMD(FDPart3_40, FDPart3_53, FusedMulAddSIMD(FDPart3_38, MulSIMD(FDPart3_65, T4DD03), FusedMulAddSIMD(FDPart3_44, MulSIMD(FDPart3_47, T4DD01), FusedMulAddSIMD(FDPart3_35, FDPart3_49, MulSIMD(FDPart3_35, FDPart3_72))))))))))))))));
                  const REAL_SIMD_ARRAY __RHS_exp_7 = FusedMulAddSIMD(T4DD22, MulSIMD(FDPart3_59, FDPart3_59), FusedMulAddSIMD(FDPart3_41, MulSIMD(FDPart3_70, T4DD01), FusedMulAddSIMD(FDPart3_68, T4DD11, FusedMulAddSIMD(FDPart3_80, T4DD33, FusedMulAddSIMD(FDPart3_60, FDPart3_70, FusedMulAddSIMD(FDPart3_62, FDPart3_81, FusedMulAddSIMD(FDPart3_57, MulSIMD(FDPart3_Integer_2, T4DD03), FusedMulAddSIMD(FDPart3_59, MulSIMD(FDPart3_75, FDPart3_Integer_2), FusedMulAddSIMD(FDPart3_1, FDPart3_4, MulSIMD(FDPart3_56, FDPart3_81))))))))));
                  const REAL_SIMD_ARRAY __RHS_exp_8 = FusedMulAddSIMD(FDPart3_41, MulSIMD(FDPart3_65, T4DD03), FusedMulAddSIMD(FDPart3_44, MulSIMD(FDPart3_59, T4DD02), FusedMulAddSIMD(FDPart3_80, T4DD23, FusedMulAddSIMD(FDPart3_33, MulSIMD(FDPart3_44, T4DD01), FusedMulAddSIMD(FDPart3_63, T4DD01, FusedMulAddSIMD(FDPart3_73, T4DD11, FusedMulAddSIMD(FDPart3_56, FDPart3_65, FusedMulAddSIMD(FDPart3_62, FDPart3_65, FusedMulAddSIMD(FDPart3_53, FDPart3_75, FusedMulAddSIMD(FDPart3_53, FDPart3_79, FusedMulAddSIMD(FDPart3_37, FDPart3_53, FusedMulAddSIMD(FDPart3_39, FDPart3_53, FusedMulAddSIMD(FDPart3_53, MulSIMD(FDPart3_59, T4DD22), FusedMulAddSIMD(FDPart3_53, MulSIMD(FDPart3_65, T4DD33), FusedMulAddSIMD(FDPart3_1, FDPart3_10, MulSIMD(FDPart3_35, FDPart3_60))))))))))))))));
                  const REAL_SIMD_ARRAY __RHS_exp_9 = FusedMulAddSIMD(FDPart3_44, MulSIMD(FDPart3_71, T4DD01), FusedMulAddSIMD(FDPart3_44, MulSIMD(FDPart3_81, T4DD02), FusedMulAddSIMD(FDPart3_80, T4DD22, FusedMulAddSIMD(T4DD33, MulSIMD(FDPart3_65, FDPart3_65), FusedMulAddSIMD(FDPart3_69, T4DD11, FusedMulAddSIMD(FDPart3_79, FDPart3_82, FusedMulAddSIMD(FDPart3_53, MulSIMD(FDPart3_71, T4DD12), FusedMulAddSIMD(FDPart3_53, MulSIMD(FDPart3_82, T4DD23), FusedMulAddSIMD(FDPart3_1, FDPart3_6, MulSIMD(FDPart3_39, FDPart3_82))))))))));
                  WriteSIMD(&T4UU00GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)], __RHS_exp_0);
                  WriteSIMD(&T4UU01GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)], __RHS_exp_1);
                  WriteSIMD(&T4UU02GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)], __RHS_exp_2);
                  WriteSIMD(&T4UU03GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)], __RHS_exp_3);
                  WriteSIMD(&T4UU11GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)], __RHS_exp_4);
                  WriteSIMD(&T4UU12GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)], __RHS_exp_5);
                  WriteSIMD(&T4UU13GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)], __RHS_exp_6);
                  WriteSIMD(&T4UU22GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)], __RHS_exp_7);
                  WriteSIMD(&T4UU23GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)], __RHS_exp_8);
                  WriteSIMD(&T4UU33GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)], __RHS_exp_9);
            
        } // END LOOP: for (int i0 = 0; i0 < cctk_lsh[0]; i0 += SIMD_width)
    } // END LOOP: for (int i1 = 0; i1 < cctk_lsh[1]; i1++)
} // END LOOP: for (int i2 = 0; i2 < cctk_lsh[2]; i2++)

}
