#include <math.h>
#include "Geodesics.hh"
#include "MersenneTwister.h" // For high quality/fast uniform rand()

extern "C" void Geodesics_Inject(CCTK_ARGUMENTS);

extern "C" void Geodesics_Inject(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if(start_t > cctk_time) return; // Save some cpu time if we aren't evolving yet...

  int lsh[1];
  int ierr = CCTK_GrouplshGN (cctkGH, 1, lsh, "Geodesics::geo_pos");
  int lnum_geo = lsh[0];

  MTRand mtrand;

  for( int i=0; i < lnum_geo; i++ ) {  
    if((active_geo[i] == -1 || geo_z[i] > 20.0 || geo_z[i] < -20.0) 
       && (cctk_iteration%output_every == 1 || output_every == 1)) {

      active_geo[i] = 0;

      if(CCTK_Equals(dist,"sphere")) {
        double randx, randy, randz, randr;
        do {
          randx = 25.0 - mtrand.rand(50.0);
          randy = 25.0 - mtrand.rand(50.0);
          randz = 25.0 - mtrand.rand(50.0);
          randr = sqrt(randx*randx+randy*randy+randz*randz);
        } while(randr < sphere_innerr || randr > sphere_outerr );

        geo_x[i] = randx;
        geo_y[i] = randy;
        geo_z[i] = randz;

      }

      if(CCTK_Equals(pdist,"zero")) {
        geo_ux[i] = 0.0;
        geo_uy[i] = 0.0;
        geo_uz[i] = 0.0;
      }

      //      CCTK_REAL rand_phi = 2.0*3.1415*mtrand.rand();

      //      CCTK_REAL rand_z = 1.0 - 2.0*mtrand.rand();
      //      CCTK_REAL geo_r = 10.0; //for testing...

      //      geo_x[i] = geo_r*cos(rand_phi);
      //      geo_y[i] = geo_r*sin(rand_phi);
      //      geo_z[i] = rand_z;

      //      CCTK_REAL vel_fac = 10000000.0;

      //      geo_ux[i] = -vel_fac*cos(rand_phi);
      //      geo_uy[i] = -vel_fac*sin(rand_phi);
      //      geo_uz[i] = -vel_fac*rand_z;

    }
  }
}
