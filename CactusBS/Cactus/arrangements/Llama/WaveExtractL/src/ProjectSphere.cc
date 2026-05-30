
/* Copyright 2013 Peter Diener, Nils Dorband, Roland Haas, Ian Hinder,
Christian Ott, Denis Pollney, Thomas Radke, Christian Reisswig, Erik
Schnetter, Barry Wardell and Burkhard Zink

This file is part of Llama.

Llama is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free
Software Foundation, either version 2 of the License, or (at your
option) any later version.

Llama is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with Llama.  If not, see <http://www.gnu.org/licenses/>. */

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

#include <vector>
#include "carpetinterp2.hh"

using namespace std;
using namespace CarpetInterp2;


struct interp_setup_t {
  fasterp_setup_t * fasterp_setup;
  int npoints;

  interp_setup_t (int npoints_) : fasterp_setup (NULL), npoints (npoints_) {}
  ~interp_setup_t () { if (fasterp_setup) delete fasterp_setup; }
};
vector<interp_setup_t *> interp_setups;   // has to be the same number as there are detectors!


#define SINDEX(lsh, i, j) j*lsh[0]+i

extern "C" void WavExtrL_ProjectSphere(CCTK_ARGUMENTS)
{
   DECLARE_CCTK_ARGUMENTS
   DECLARE_CCTK_PARAMETERS
   
   int ierr = 0;
   int di;
   int lsh[2], lbnd[2];
   
   if (verbose>2) CCTK_INFO("Interpolating metric and derivatives onto sphere");

   if (*do_nothing == 1) return;

   if (cctk_iteration != 0)
   {
       if (cctk_iteration % my_out_every_det[*current_detector-1] != 0)
       {
         if (verbose>2) CCTK_INFO("No time for this detector");
         return;
      }
   }

   *current_detector_radius = detector_radius[*current_detector-1];
   di=*current_detector-1;

   // Sanity Check
   if (*current_detector_radius < 1.e-10)
   {
      *do_nothing = 1;
      CCTK_WARN(1,"This should never happen: The detector radius is 0!");
      return;
   }

   if (calc_when_necessary == 1)
   {
      if (cctk_time < *current_detector_radius-50)
      {
         if (verbose>2) CCTK_INFO("No time for this detector");
         return;
      }
      int istat = CCTK_IsFunctionAliased("MergerHandler_WeHaveMerger");
      if (istat == 1)
      {
         if (MergerHandler_WeHaveMerger() == 1)
         {
            if (cctk_time > MergerHandler_MergerTime()+*current_detector_radius+ringdown_margin)
            {
               if (verbose>2) CCTK_INFO("No time for this detector");
               return;
            }
         }
      }
   }
    

   if (verbose>1)
      cout << "Analysing Detector No.: " << *current_detector << " Radius " << *current_detector_radius << endl;



   // local shape of the 2D grid arrays
   ierr = CCTK_GrouplshGN(cctkGH, 2, lsh, "WaveExtractL::surface_arrays");
   if ( ierr < 0 )
      CCTK_WARN(0, "cannot get local size for surface arrays");

   
   int number_of_vars = 12;

   vector<string> varnames(number_of_vars);
   
   varnames[0] = "admbase::gxx";
   varnames[1] = "admbase::gxy";
   varnames[2] = "admbase::gxz";
   varnames[3] = "admbase::gyy";
   varnames[4] = "admbase::gyz";
   varnames[5] = "admbase::gzz";
   varnames[6] = "admderivatives::gxx_dr";
   varnames[7] = "admderivatives::gxy_dr";
   varnames[8] = "admderivatives::gxz_dr";
   varnames[9] = "admderivatives::gyy_dr";
   varnames[10] = "admderivatives::gyz_dr";
   varnames[11] = "admderivatives::gzz_dr";

   // get the input variable's index
   vector<int> varindices (number_of_vars, -1);
   
   for (int i=0; i < number_of_vars; ++i)
   {
      varindices[i] = CCTK_VarIndex(varnames[i].c_str());
      if (varindices[i] < 0)
         CCTK_VWarn(CCTK_WARN_ABORT, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "couldn't get index of slice variable '%s'", varnames[i].c_str());
   }
   
   vector<CCTK_REAL*> values (number_of_vars, static_cast<CCTK_REAL*>(NULL));
   
   // set output pointers
   values[0] = gxxi;
   values[1] = gxyi;
   values[2] = gxzi;
   values[3] = gyyi;
   values[4] = gyzi;
   values[5] = gzzi;
   values[6] = dr_gxxi;
   values[7] = dr_gxyi;
   values[8] = dr_gxzi;
   values[9] = dr_gyyi;
   values[10] = dr_gyzi;
   values[11] = dr_gzzi;
   
   int interpolator_order = interpolation_stencil;
   
   // set up the interpolation
   interp_setups.resize(maximum_detector_number, NULL);
   assert(di > 0);
   assert(interp_setups.size() > (unsigned)di);
   interp_setup_t* &interp_setup = interp_setups[di];
   if (not interp_setup or interp_setup->fasterp_setup->outofdate()) {
      if (interp_setup)
        delete interp_setup;
      interp_setup = new interp_setup_t(lsh[0]*lsh[1]);

      // allocate storage for coordinates
      fasterp_glocs_t locations (interp_setup->npoints);

      // get Cartesian coordinate values of gridpoints on spherical surface
      CCTK_REAL rad = *current_detector_radius;
      for (int ip=0; ip < lsh[1]; ++ip) {
         for (int it=0; it < lsh[0]; ++it) {
            const int l = SINDEX(lsh, it, ip);
            locations.coords[0][l] = origin_x + rad * sintheta[l] * cosphi[l];
            locations.coords[1][l] = origin_y + rad * sintheta[l] * sinphi[l];
            locations.coords[2][l] = origin_z + rad * costheta[l];
         }
      }

      interp_setup->fasterp_setup =
         new fasterp_setup_t(cctkGH, locations, interpolator_order);
   }

   // do the interpolation
   assert(interp_setup->fasterp_setup);
   assert(interp_setup->npoints == lsh[0]*lsh[1]);
   interp_setup->fasterp_setup->interpolate (cctkGH, varindices, values);

   // Find out the lower bounds of the distributed integration grid arrays.
   ierr = CCTK_GrouplbndGN(cctkGH,2,lbnd,"WaveExtractL::surface_integrands");
   if ( ierr < 0 )
      CCTK_WARN(0, "cannot get lower bounds for surface integrands");

   // Convert to spherical coord. system
   // Note that these equations do not take the conformal factor into
   // account. That has to be done afterwards !!
   // FIXME : WARN, THESE EQUATIONS PROBABLY CANT WORK FOR PHYSICAL METRIC
   CCTK_REAL rad = *current_detector_radius;
   CCTK_REAL r2=rad*rad;
   for (int ip = 0; ip < lsh[1]; ++ip)
   {
      for (int it = 0; it < lsh[0]; ++it)
      {
         const int l = SINDEX(lsh, it, ip);
         
         CCTK_REAL ct = costheta[l]; 
         CCTK_REAL ct2 = ct*ct;
         CCTK_REAL st = sintheta[l]; 
         CCTK_REAL st2 = st*st;
         CCTK_REAL cp = cosphi[l]; 
         CCTK_REAL cp2 = cp*cp;
         CCTK_REAL sp = sinphi[l]; 
         CCTK_REAL sp2 = sp*sp;
         
         CCTK_REAL two = 2.0;
         
         CCTK_REAL tgxx = gxxi[l]; 
         CCTK_REAL dgxx = dr_gxxi[l];
         CCTK_REAL tgxy = gxyi[l]; 
         CCTK_REAL dgxy = dr_gxyi[l];
         CCTK_REAL tgxz = gxzi[l];
         CCTK_REAL dgxz = dr_gxzi[l];
         CCTK_REAL tgyy = gyyi[l];
         CCTK_REAL dgyy = dr_gyyi[l];
         CCTK_REAL tgyz = gyzi[l];
         CCTK_REAL dgyz = dr_gyzi[l];
         CCTK_REAL tgzz = gzzi[l]; 
         CCTK_REAL dgzz = dr_gzzi[l];
   
         grr[l] = st2*cp2*tgxx +st2*sp2*tgyy +ct2*tgzz 
                  +two*( st2*cp*sp*tgxy +st*cp*ct*tgxz +st*ct*sp*tgyz);
   
         grt[l] = rad*(st*cp2*ct*tgxx +two*st*ct*sp*cp*tgxy 
                  +cp*(ct2-st2)*tgxz +st*sp2*ct*tgyy 
                  +sp*(ct2-st2)*tgyz-ct*st*tgzz);
   
         grp[l] = rad*st*(-st*sp*cp*tgxx -st*(sp2-cp2)*tgxy 
                  -sp*ct*tgxz +st*sp*cp*tgyy +ct*cp*tgyz);
   
         gtt[l] = r2*(ct2*cp2*tgxx +two*ct2*sp*cp*tgxy 
                  -two*st*ct*cp*tgxz +ct2*sp2*tgyy 
                  -two*st*sp*ct*tgyz +st2*tgzz);
   
         gtp[l] = r2*st*(-cp*sp*ct*tgxx -ct*(sp2-cp2)*tgxy 
                  +st*sp*tgxz +cp*sp*ct*tgyy -st*cp*tgyz);
   
         gpp[l] = r2*st2*(sp2*tgxx -two*cp*sp*tgxy +cp2*tgyy);
   
   
         dr_gtt[l] = two/rad*gtt[l] +r2*(ct2*cp2*dgxx 
                        +two*ct2*sp*cp*dgxy -two*st*ct*cp*dgxz 
                        +ct2*sp2*dgyy -two*st*sp*ct*dgyz +st2*dgzz);
   
         dr_gtp[l] = two/rad*gtp[l] +r2*st*(-cp*sp*ct*dgxx 
                        -ct*(sp2-cp2)*dgxy +st*sp*dgxz 
                        +cp*sp*ct*dgyy -st*cp*dgyz);
   
         dr_gpp[l] = two/rad*gpp[l] +r2*st2*(sp2*dgxx 
                        -two*cp*sp*dgxy+cp2*dgyy);
                        
         if (it+lbnd[0]>=int_ntheta[*current_detector-1] ||
            ip+lbnd[1]>=int_nphi[*current_detector-1] )
         {
            grr[l]=0;
            grt[l]=0;
            grp[l]=0;
            gtt[l]=0;
            gtp[l]=0;
            gpp[l]=0;
            dr_gtt[l]=0;
            dr_gtp[l]=0;
            dr_gpp[l]=0;
         }
      }
   }

}
