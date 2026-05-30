/*
 * This program sets all the HydroBase primitive variables untouched by FFE to zero

 Please consult GiRaFFEfood/README for the list of authors and current maintainers
 * */

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <ctime>
#include <cstdlib>

using namespace std;

extern "C" void GiRaFFEfood_Set_HydroBase_to_Zero(CCTK_ARGUMENTS) {

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_INT imin = 0, jmin = 0, kmin = 0;
  CCTK_INT imax = cctk_lsh[0], jmax = cctk_lsh[1], kmax = cctk_lsh[2];

#pragma omp parallel for
  for(int k=kmin;k<kmax;k++)
    for(int j=jmin;j<jmax;j++)
      for(int i=imin;i<imax;i++) {
        int index = CCTK_GFINDEX3D(cctkGH,i,j,k);

        rho[index] = 0.0;
        press[index] = 0.0;
        eps[index] = 0.0;
        // FIXME: Set velocities to zero (conceivable we might want to add a case in which v's are not zero!
        vel[CCTK_GFINDEX4D(cctkGH,i,j,k,0)] = 0.0;
        vel[CCTK_GFINDEX4D(cctkGH,i,j,k,1)] = 0.0;
        vel[CCTK_GFINDEX4D(cctkGH,i,j,k,2)] = 0.0;

        //Set the magnetic potential to zero as well, before filling in with particular initial data
        Aphi[index] = 0.0;
        Avec[CCTK_GFINDEX4D(cctkGH,i,j,k,0)] = 0.0;
        Avec[CCTK_GFINDEX4D(cctkGH,i,j,k,1)] = 0.0;
        Avec[CCTK_GFINDEX4D(cctkGH,i,j,k,2)] = 0.0;

      }
}

