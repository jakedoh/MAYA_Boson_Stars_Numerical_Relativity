#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

#include "mpi.h"

#include <iostream>
#include <vector>

using namespace std;

inline CCTK_REAL interpolant(const CCTK_REAL t, const CCTK_REAL f0, const CCTK_REAL f1, const CCTK_REAL f2, const CCTK_REAL f3, const CCTK_REAL f4, 
                                                const CCTK_REAL t0, const CCTK_REAL t1, const CCTK_REAL t2, const CCTK_REAL t3, const CCTK_REAL t4)
{
   CCTK_REAL res;

   res = f0 * (t-t1)/(t0-t1) * (t-t2)/(t0-t2) * (t-t3)/(t0-t3) * (t-t4)/(t0-t4)
       + f1 * (t-t0)/(t1-t0) * (t-t2)/(t1-t2) * (t-t3)/(t1-t3) * (t-t4)/(t1-t4)
       + f2 * (t-t0)/(t2-t0) * (t-t1)/(t2-t1) * (t-t3)/(t2-t3) * (t-t4)/(t2-t4)
       + f3 * (t-t0)/(t3-t0) * (t-t1)/(t3-t1) * (t-t2)/(t3-t2) * (t-t4)/(t3-t4)
       + f4 * (t-t0)/(t4-t0) * (t-t1)/(t4-t1) * (t-t2)/(t4-t2) * (t-t3)/(t4-t3);

   return res;
}


#define SINDEX2D(p, i, j) (p*null_lsh[0]*null_lsh[1] + j*null_lsh[0] + i)


extern "C" void NullNews_Interp_to_constant_uBondi(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  
  // define and point pointer arrays to the past timelevels
#include "set_pointers.h"

  // find min value on 1st, 2nd and 3rd timelevel
  double min, min_p, min_p_p;
  min = min_p = min_p_p = 1.0e20;
  for (int p=0; p < 2; ++p)
     for (int j=0; j < null_lsh[1]; ++j)
        for (int i=0; i < null_lsh[0]; ++i)
        {
           const int ij = SINDEX2D(p, i, j);
           // only compute on nominal grid!
           if (EV_mask[j*null_lsh[0] + i] > 0) {
              if (uBondi[ij] < min) min = uBondi[ij];
              if (uBondiP[1][ij] < min_p) min_p = uBondiP[1][ij];
              if (uBondiP[2][ij] < min_p_p) min_p_p = uBondiP[2][ij];
           }
        }
  
  double min_send[3] = { min, min_p, min_p_p };
  double min_receive[3];
  MPI_Allreduce(&min_send, &min_receive, 3, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

  min = min_receive[0];
  min_p = min_receive[1];
  min_p_p = min_receive[2];

  if (!(min > min_p && min_p > min_p_p))
     CCTK_WARN(0, "Minimum Bondi times are not strictly monotonic!");

  // next, we have to find the interpolation target Bondi-time
  // This will be the average time on the 3rd timelevel
  
  CCTK_REAL avg_uBondi = 0;
  CCTK_REAL npoints = 0;
  for (int p=0; p < 2; ++p)
     for (int j=0; j < null_lsh[1]; ++j)
        for (int i=0; i < null_lsh[0]; ++i)
        {
           const int ij = SINDEX2D(p, i, j);
           // only compute on nominal grid!
           if (EV_mask[j*null_lsh[0] + i] > 0) {
              avg_uBondi += uBondiP[2][ij];
              ++npoints;
           }
        }
   
   double sum_send[2] = { avg_uBondi, npoints };
   double sum_receive[2];
   MPI_Allreduce(&sum_send, &sum_receive, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

   *constant_uBondi = sum_receive[0] / sum_receive[1];

   // Our new average time is not allowed to be bigger than minimum uBondi on 2nd timelevel!

   if (*constant_uBondi > 0.5*(min_p+min_p_p)) *constant_uBondi = 0.5*(min_p+min_p_p);

   cout << "constant_uBondi = " << *constant_uBondi << endl;
   /*cout << "min = " << min << endl;
   cout << "min_p = " << min_p << endl;
   cout << "min_p_p = " << min_p_p << endl;*/

   // interpolate
   for (int p=0; p < 2; ++p)
     for (int j=0; j < null_lsh[1]; ++j)
        for (int i=0; i < null_lsh[0]; ++i)
        {
           const int ij = SINDEX2D(p, i, j);
           CCTK_REAL real, imag;
           
           // only compute on nominal grid!
           if (EV_mask[j*null_lsh[0] + i] == 0) continue;
           
           if (*constant_uBondi > uBondi[ij])
              CCTK_WARN(0, "Severe: constant_uBondi more advanced than latest uBondi!");
           
           if (*constant_uBondi < uBondiP[max_timelevels-1][ij])
              CCTK_WARN(0, "constant_uBondi lagging behind earliest uBondi! Increase your number of max_timelevels.");
           
           // identify past timelevel which is smaller than our current target constant Bondi time
           // for the current point
           int tl = max_timelevels-3;
           for (int t=1; t < max_timelevels-2; ++t)
              if (0.5*(uBondiP[t][ij]+uBondiP[t-1][ij]) >= *constant_uBondi && 
                  0.5*(uBondiP[t+1][ij]+uBondiP[t][ij]) < *constant_uBondi) {
                 tl = t;
                 break;
              }
           
           //cout << tl << "  ";
           
           if (*constant_uBondi > uBondiP[tl-2][ij])
           {
              cout << "uBondi[ij] = " << uBondi[ij] << endl;
              cout << "uBondiP[" << tl-2 << "][ij] = " << uBondiP[tl-2][ij] << endl;
              CCTK_WARN(0, "Severe: constant_uBondi more advanced than interpolation intervall!");
           }
           
           if (*constant_uBondi < uBondiP[tl+2][ij])
              CCTK_WARN(0, "Severe: constant_uBondi lagging behind interpolation intervall!");
           
           
           
           real = interpolant(*constant_uBondi,
                              Psi4P[tl-2][ij].real(), Psi4P[tl-1][ij].real(), Psi4P[tl][ij].real(), Psi4P[tl+1][ij].real(), Psi4P[tl+2][ij].real(),
                              uBondiP[tl-2][ij], uBondiP[tl-1][ij], uBondiP[tl][ij], uBondiP[tl+1][ij], uBondiP[tl+2][ij]);
           imag = interpolant(*constant_uBondi,
                              Psi4P[tl-2][ij].imag(), Psi4P[tl-1][ij].imag(), Psi4P[tl][ij].imag(), Psi4P[tl+1][ij].imag(), Psi4P[tl+2][ij].imag(),
                              uBondiP[tl-2][ij], uBondiP[tl-1][ij], uBondiP[tl][ij], uBondiP[tl+1][ij], uBondiP[tl+2][ij]);
           Psi4_uBondi[ij] = CCTK_Cmplx(real, imag);
           
           real = interpolant(*constant_uBondi,
                              NewsP[tl-2][ij].real(), NewsP[tl-1][ij].real(), NewsP[tl][ij].real(), NewsP[tl+1][ij].real(), NewsP[tl+2][ij].real(),
                              uBondiP[tl-2][ij], uBondiP[tl-1][ij], uBondiP[tl][ij], uBondiP[tl+1][ij], uBondiP[tl+2][ij]);
           imag = interpolant(*constant_uBondi,
                              NewsP[tl-2][ij].imag(), NewsP[tl-1][ij].imag(), NewsP[tl][ij].imag(), NewsP[tl+1][ij].imag(), NewsP[tl+2][ij].imag(),
                              uBondiP[tl-2][ij], uBondiP[tl-1][ij], uBondiP[tl][ij], uBondiP[tl+1][ij], uBondiP[tl+2][ij]);
           News_uBondi[ij] = CCTK_Cmplx(real, imag);
           
           real = interpolant(*constant_uBondi,
                              NewsBP[tl-2][ij].real(), NewsBP[tl-1][ij].real(), NewsBP[tl][ij].real(), NewsBP[tl+1][ij].real(), NewsBP[tl+2][ij].real(),
                              uBondiP[tl-2][ij], uBondiP[tl-1][ij], uBondiP[tl][ij], uBondiP[tl+1][ij], uBondiP[tl+2][ij]);
           imag = interpolant(*constant_uBondi,
                              NewsBP[tl-2][ij].imag(), NewsBP[tl-1][ij].imag(), NewsBP[tl][ij].imag(), NewsBP[tl+1][ij].imag(), NewsBP[tl+2][ij].imag(),
                              uBondiP[tl-2][ij], uBondiP[tl-1][ij], uBondiP[tl][ij], uBondiP[tl+1][ij], uBondiP[tl+2][ij]);
           NewsB_uBondi[ij] = CCTK_Cmplx(real, imag);
           
           if (compute_lin_strain) {
              real = interpolant(*constant_uBondi,
                                 linStrainP[tl-2][ij].real(), linStrainP[tl-1][ij].real(), linStrainP[tl][ij].real(), linStrainP[tl+1][ij].real(), linStrainP[tl+2][ij].real(),
                                 uBondiP[tl-2][ij], uBondiP[tl-1][ij], uBondiP[tl][ij], uBondiP[tl+1][ij], uBondiP[tl+2][ij]);
              imag = interpolant(*constant_uBondi,
                                 linStrainP[tl-2][ij].imag(), linStrainP[tl-1][ij].imag(), linStrainP[tl][ij].imag(), linStrainP[tl+1][ij].imag(), linStrainP[tl+2][ij].imag(),
                                 uBondiP[tl-2][ij], uBondiP[tl-1][ij], uBondiP[tl][ij], uBondiP[tl+1][ij], uBondiP[tl+2][ij]);
              linStrain_uBondi[ij] = CCTK_Cmplx(real, imag);
           }
        }
   
}


extern "C" void NullNews_InterpInit(CCTK_ARGUMENTS)
{
   DECLARE_CCTK_ARGUMENTS
   DECLARE_CCTK_PARAMETERS

   // define and point pointer arrays to the past timelevels
#include "set_pointers.h"

   // initialize previous timelevels to coincide with coordinate time (offset by 0.5*delta_t because
   // News is computed at staggered timelevels.
   for (int tl=0; tl < max_timelevels; ++tl)
     for (int p=0; p < 2; ++p)
        for (int j=0; j < null_lsh[1]; ++j)
           for (int i=0; i < null_lsh[0]; ++i)
           {
              const int ij = SINDEX2D(p, i, j);
              uBondiP[tl][ij]  = -(tl-0.5)*cctk_delta_time;
              
              Psi4P[tl][ij]    = CCTK_Cmplx(0,0);
              
              NewsP[tl][ij]    = CCTK_Cmplx(0,0);
              
              NewsBP[tl][ij]   = CCTK_Cmplx(0,0);
              
              if (compute_lin_strain)
                 linStrainP[tl][ij] = CCTK_Cmplx(0,0);
           }

   *constant_uBondi = uBondiP[2][0];
   cout << "Initial constant_uBondi = " << *constant_uBondi << endl;
}



extern "C" void NullNews_InterpCycleTimelevels(CCTK_ARGUMENTS)
{
   DECLARE_CCTK_ARGUMENTS
   DECLARE_CCTK_PARAMETERS

   // define and point pointer arrays to the past timelevels
#include "set_pointers.h"

   // copy current timelevels to previous ones.
   for (int tl=max_timelevels-2; tl >= 0; --tl)
     for (int p=0; p < 2; ++p)
        for (int j=0; j < null_lsh[1]; ++j)
           for (int i=0; i < null_lsh[0]; ++i)
           {
              const int ij = SINDEX2D(p, i, j);
              
              uBondiP[tl+1][ij] = uBondiP[tl][ij];
              Psi4P[tl+1][ij]   = Psi4P[tl][ij];
              NewsP[tl+1][ij]   = NewsP[tl][ij];
              NewsBP[tl+1][ij]  = NewsBP[tl][ij];
              if (compute_lin_strain)
                 linStrainP[tl+1][ij] = linStrainP[tl][ij];
           }
}


