#include <stdio.h>

#include "assert.h"
#include "math.h"
#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "EOS_Base.h"
#include "matterlibs.h"
#include "steffenson.h" // for gsl_test_interval
#include "HydroBase.h"

#define DEBUG
#define MINB2 1e-24

static inline double pow2 (const double x)
{
  return x*x;
}
static inline double pow3 (const double x)
{
  return x*x*x;
}
static inline double pow4 (const double x)
{
  const double y = x*x;
  return y*y;
}

void MatterLibs_SetPrims(CCTK_ARGUMENTS)
{
		  DECLARE_CCTK_ARGUMENTS;
		  DECLARE_CCTK_PARAMETERS;

	if(CCTK_EQUALS(combination_method, "rescale or background only") && rescale_hydro && provide_momentum)
	{
			  MatterLibs_CalculateRho(CCTK_PASS_CTOC);
	}
	else if((CCTK_EQUALS(combination_method, "rescale or background only") && rescale_hydro) || CCTK_EQUALS(combination_method, "rescale"))
	{
			  MatterLibs_FeedRho(CCTK_PASS_CTOC);
	}
   else if((CCTK_EQUALS(combination_method, "rescale or background only") && !rescale_hydro) || CCTK_EQUALS(combination_method, "background only"))
	{
			  MatterLibs_SetRho(CCTK_PASS_CTOC);
	}
   else if(CCTK_EQUALS(combination_method, "sum metric"))
	{
		     MatterLibs_SetRhoMetricSummed(CCTK_PASS_CTOC);
	}
}

void MatterLibs_CalculateRho(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  enum eos_type_enum {eos_invalid = 0, eos_polytype, eos_general};
  enum initial_guess_type_enum {guess_invalid = 0, guess_rescaled_E, guess_rho, guess_rescaled_rho};
  enum eos_type_enum eos_type = eos_invalid;
  enum initial_guess_type_enum initial_guess_type = guess_invalid;
  const char * eos_type_name;
  double dpsi, drho, deps, dpress, dvelx, dvely, dvelz; // primitive initial data
  double dE, dD, dtau, dJx, dJy, dJz;                   // conservative initial data
  double dEbar, dJbarx, dJbary, dJbarz;// conformal initial data
  double dlorentz2, dlorentz, dh, B2, vB;		// auxiliary variables
  double xx, yy, zz;
  double psiL, chiL, rootdetg;
  double twelfth, third;
  double g11, g21, g31, g22, g32, g33, detg;
  double uxx, uxy, uxz, uyy, uyz, uzz;
  static double linmom[3];                              // total linear momentum of matter
  static double angmom[3];                              // total angular momentum of matter 
  FILE *testsuitefh=NULL;
  static int firsttime = 1;				// chain output together 

  // check whether we have a B-field and set pointers accordingly
  int have_uniform_bfield=0;
  double dbvecx, dbvecy, dbvecz; // primitive B-field initial data (conservative is trivial)
  if ( b0_uniform[0] != 0 || b0_uniform[1] != 0 || b0_uniform[2] != 0 ) {
     have_uniform_bfield=1;
  }
     

  // open file to dump internal state to
  if(testsuite)
  {
    char fn[1024];

    snprintf(fn, sizeof(fn)/sizeof(fn[0]), "%s/matterlibs_internaldata.asc", out_dir);
    if((testsuitefh = fopen(fn, firsttime ? "w" : "a")) == NULL)
      CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,  "Cannot open file %s to dump internal data", fn);

    if(firsttime) {
       if ( have_uniform_bfield ) {
          fprintf(testsuitefh, "# 1:x 2:y 3:z 4:psi 5:rho 6:eps 7:press 8:velx 9:vely 10:velz 11:bx 12:by 13:bz 14:psi_final\n");
       } else {
          fprintf(testsuitefh, "# 1:x 2:y 3:z 4:psi 5:rho 6:eps 7:press 8:velx 9:vely 10:velz 11:psi_final\n");
       }
    }

    firsttime = 0;
  }

  // find out what kind of EOS we are using
  if (CCTK_EQUALS(matterlibs_eos_type, ""))
      eos_type_name = whisky_eos_type;
  else
      eos_type_name = matterlibs_eos_type;
  if (CCTK_EQUALS(eos_type_name, "Polytype")) {
      eos_type = eos_polytype;
      assert(*whisky_polytrope_handle>=0);
  } else if (CCTK_EQUALS(eos_type_name, "General")) {
      CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,  "MatterLibs_CalculateRho can only be used with EOS type 'Polytype'");
  } else {
      CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,  "Unknown EOS type %s is not 'Polytype' or 'General'.", eos_type_name);
  }
  assert(eos_type);

  // decide how to provide initial guesses for the E_to_rho solver
  if (CCTK_EQUALS(initial_guess, "rescaled E")) {
      initial_guess_type = guess_rescaled_E;
  } else if (CCTK_EQUALS(initial_guess, "rho")) {
      initial_guess_type = guess_rho;
  } else if (CCTK_EQUALS(initial_guess, "rescaled rho")) {
      initial_guess_type = guess_rescaled_rho;
  } else {
      CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,  "Unknown initial guess type '%s'.", initial_guess);
  }
  assert(initial_guess_type);

  // Whisky kindly sets the initial rho_min in Whisky_Rho_Minima_Setup

  twelfth = 1.0/12.0;
  third = 1.0/3.0;

  CCTK_REAL *CarpetWeights = (CCTK_REAL *)(CCTK_VarDataPtr(cctkGH, 0, "CarpetReduce::weight"));
  assert(CarpetWeights);
  for ( int kk = 0 ; kk < cctk_lsh[2] ; kk++) {
    for ( int j = 0 ; j < cctk_lsh[1] ; j++) {
      for (int i = 0 ; i < cctk_lsh[0] ; i++) {

        const int ind = CCTK_GFINDEX3D (cctkGH, i, j, kk);

        xx = x[ind];
        yy = y[ind];
        zz = z[ind];

        g11 = gxx[ind];
        g22 = gyy[ind];
        g33 = gzz[ind];
        g21 = gxy[ind];
        g31 = gxz[ind];
        g32 = gyz[ind];
        detg  =  2*g21*g31*g32 + g33*(g11*g22 - pow2(g21)) - g22*pow2(g31) - g11*pow2(g32);

	chiL = fmax(1/pow(detg,third),chiEps);
	// psiL is used in the rescaling from Ebar to E, it should not be
	// capped from above since that makes the initial data potentially
	// large. rootdetg on the other hand is only used to densitize the
	// conservative varibles at the end, here we want to cap it from above
	// to limit the magnitude of the ID
        // psiL = 1/pow(chiL,1.0/4.0); 
	rootdetg = 1/pow(chiL,3.0/2.0);
        psiL = pow(detg,twelfth);
	/*rootdetg = sqrt(detg);*/

        // obtain guesses for the primitive variables
        calcConformalRho( xx, yy, zz, &dpsi, &drho, &dEbar, &dJbarx, &dJbary, &dJbarz, &dbvecx, &dbvecy, &dbvecz, &vB, &dlorentz );
        dpsi = combineInitialGuesses(guessPuncturePsi(xx, yy, zz), dpsi);
        if ( have_uniform_bfield > 0 && CCTK_EQUALS(initial_Bvec,"uniformBn") ) {
           /* divide by new psi since b0_uniform should be the *conserved* bfield */
           CCTK_REAL inv_rdetg=pow(dpsi,-6);
           dbvecx = dbvecx*inv_rdetg;
           dbvecy = dbvecy*inv_rdetg;
           dbvecz = dbvecz*inv_rdetg;
        }

		/*
	if(testsuite) // dump internal data
	{
          fprintf(testsuitefh, "%g %g %g %g %g %g %g %g %g %g %g %g %g %g\n", 
                  xx, yy, zz, dpsi, drho, deps, dpress, dvelx,
		  dvely, dvelz, dbvecx, dbvecy, dbvecz, psiL);
	}
	*/

	if( drho < *whisky_rho_min*(1.0+atmo_tolerance) ) // enforce at least atmospheric density
	{
	      drho = *whisky_rho_min;
	      dpress = EOS_Pressure(*whisky_polytrope_handle, drho, matterlibs_poison);
	      deps = EOS_SpecificIntEnergy(*whisky_polytrope_handle, drho, matterlibs_poison);
	      dvelx = dvely = dvelz = 0.;
	}

        // calculate auxiliary variables
        dlorentz2 = pow2(dlorentz);
//        dh = 1 + deps + dpress/drho;
        B2 = pow4(dpsi) * ( pow2(dbvecx) + pow2(dbvecy) + pow2(dbvecz) );
//        vB = pow4(dpsi) * ( dvelx*dbvecx + dvely*dbvecy + dvelz*dbvecz );

        // Calculate linear and angular momentum through raw primitives
        if (CarpetWeights[ind] == 1.0 &&
             i >= cctk_nghostzones[0] &&  i < cctk_lsh[0]-cctk_nghostzones[0] &&
             j >= cctk_nghostzones[1] &&  j < cctk_lsh[1]-cctk_nghostzones[1] &&
            kk >= cctk_nghostzones[2] && kk < cctk_lsh[2]-cctk_nghostzones[2]) {

           // Calculate appropriate values for E and S_j after baring and unbaring 
           // as in gr-qc/0606104.

           linmom[0] += dJbarx*CCTK_DELTA_SPACE(0)*CCTK_DELTA_SPACE(1)*CCTK_DELTA_SPACE(2); 
           linmom[1] += dJbary*CCTK_DELTA_SPACE(0)*CCTK_DELTA_SPACE(1)*CCTK_DELTA_SPACE(2); 
           linmom[2] += dJbarz*CCTK_DELTA_SPACE(0)*CCTK_DELTA_SPACE(1)*CCTK_DELTA_SPACE(2); 
           angmom[0] += (yy*dJbarz-zz*dJbary)*CCTK_DELTA_SPACE(0)*CCTK_DELTA_SPACE(1)*CCTK_DELTA_SPACE(2);
           angmom[1] += (zz*dJbarx-xx*dJbarz)*CCTK_DELTA_SPACE(0)*CCTK_DELTA_SPACE(1)*CCTK_DELTA_SPACE(2);
           angmom[2] += (xx*dJbary-yy*dJbarx)*CCTK_DELTA_SPACE(0)*CCTK_DELTA_SPACE(1)*CCTK_DELTA_SPACE(2);

        }
        // DEBUG

        // Calculate appropriate values for E and S_j after baring and unbaring 
        // as in gr-qc/0606104.
        dE = dEbar/pow(psiL, conformal_density_power);
		dJx = dJbarx/pow(psiL, 10.);
		dJy = dJbary/pow(psiL, 10.);
		dJz = dJbarz/pow(psiL, 10.);

          // for a polytropic EOS D is linked to E and S^j

          // recover primitive values "by hand"
          switch(initial_guess_type) // provide initial guess to solver
          {
            case guess_rescaled_E:
              drho = dE;
              break;
            case guess_rho:
              drho = drho;
              break;
            case guess_rescaled_rho:
              drho = drho * pow(dpsi/psiL, conformal_density_power);
              break;
            case guess_invalid:
            default:
              assert(0);
              break;
          }
          int Etorho_err;
          if ( have_uniform_bfield == 0 || B2 < MINB2 ) {
               Etorho_err = matterlibs_E_to_rho(*whisky_polytrope_handle, dE, dJx, dJy, dJz, psiL,
                                      &drho, &dlorentz);
          } else {
               Etorho_err = matterlibs_E_to_rho_MHD(*whisky_polytrope_handle, dE, dJx, dJy, dJz, 
                              dbvecx, dbvecy, dbvecz, psiL, &drho, &dlorentz);
          }  

	  if( Etorho_err != 0 )
	  {
	      CCTK_VWarn(0,__LINE__,__FILE__,CCTK_THORNSTRING,
		      "initial recovery of primitive variables failed at "
		      "(%g,%g,%g), where (E,J^j)=[%g,%g,%g,%g], "
		      "psiL = %g, dpsi = %g => (rho,w_lorentz) = (%g,%g)",
		      xx,yy,zz, dE, dJx,dJy,dJz, psiL,dpsi, drho,dlorentz);
	      return; /* NOTREACHED */
	  }

	  if( drho < *whisky_rho_min*(1.0+atmo_tolerance) ) // enforce at least atmospheric density
	  {
	    drho = *whisky_rho_min;
	    dpress = EOS_Pressure(*whisky_polytrope_handle, drho, matterlibs_poison); // atmosphere always used a polytrope
	    deps = EOS_SpecificIntEnergy(*whisky_polytrope_handle, drho, matterlibs_poison);
	    dvelx = dvely = dvelz = 0.;
	    dlorentz2 = dlorentz = 1.;
            vB = 0.;

	    dD = drho * dlorentz;
	    dJx = dJy = dJz = 0.;
	    // hand crafted to guarantee positivity
	    dtau = drho*(dlorentz2-dlorentz) + drho*deps*dlorentz2 + dpress*(dlorentz2-1.) + (1 - 0.5/dlorentz2)*B2 - 0.5*pow2(vB);
	  }
	  else
	  {
	    dpress = EOS_Pressure(*whisky_polytrope_handle, drho, matterlibs_poison);
	    deps = EOS_SpecificIntEnergy(*whisky_polytrope_handle, drho, matterlibs_poison);
	    dlorentz2 = pow2(dlorentz);
        dh = 1 + deps + dpress/drho;
	    dvelx = ( dJx + vB*pow4(dpsi)*dbvecx ) / (drho * dh * dlorentz2 + B2);
	    dvely = ( dJy + vB*pow4(dpsi)*dbvecy ) / (drho * dh * dlorentz2 + B2);
	    dvelz = ( dJz + vB*pow4(dpsi)*dbvecz ) / (drho * dh * dlorentz2 + B2);

	    // calculate a suitable matter density (this is the reason we go 
	    // through the whole rigmalore...)
	    dD = drho * dlorentz;

	    // hand crafted to guarantee positivity
	    dtau = drho*(dlorentz2-dlorentz) + drho*deps*dlorentz2 + dpress*(dlorentz2-1.) + (1 - 0.5/dlorentz2)*B2 - 0.5*pow2(vB);

            // check result of conversion for consistency
            if(verbose >= 1 &&
              GSL_SUCCESS != matterlibs_gsl_root_test_interval(dE, drho*dh*dlorentz2-dpress+(1-0.5/dlorentz2)*B2 - 0.5*pow2(vB), 
                                                    abs_err, rel_err)) 
            {
              CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,  
                      "Back and forth conversion of relativistic energy density E failed: "
                      "E(rho) != E_orig; E(rho) = %g, rho = %g, E_orig = %g, diff = %e", 
                      drho*dh*dlorentz2-dpress +(1-0.5/dlorentz2)*B2 - 0.5*pow2(vB), 
                      drho, dE, dE - (drho*dh*dlorentz2-dpress+(1-0.5/dlorentz2)*B2-0.5*pow2(vB)));
            }
	  }
           
          // factor in the determinant that Whisky uses
          dJx *= rootdetg;
          dJy *= rootdetg;
          dJz *= rootdetg;
          dD *= rootdetg;
          dtau *= rootdetg;


        // fill in grid functions (Si = gij S^j!)
        dens[ind] = dD;
        sx[ind] = g11*dJx+g21*dJy+g31*dJz;
        sy[ind] = g21*dJx+g22*dJy+g32*dJz;
        sz[ind] = g31*dJx+g32*dJy+g33*dJz;
        tau[ind] = dtau;
        rho[ind] = drho;
#ifdef HAVE_MAYA_WHISKY
//        velx[ind] = dvelx;
        velx[ind] = dvelx;
        vely[ind] = dvely;
        velz[ind] = dvelz;
#else
        const int xidx = CCTK_VECTGFINDEX3D (cctkGH, i, j, kk, 0);
        const int yidx = CCTK_VECTGFINDEX3D (cctkGH, i, j, kk, 1);
        const int zidx = CCTK_VECTGFINDEX3D (cctkGH, i, j, kk, 2);
       // Migu vel[xidx] = dvelx;
        vel[xidx] = dvelx;
        vel[yidx] = dvely;
        vel[zidx] = dvelz;
#endif
        eps[ind] = deps;
        press[ind] = dpress;
        w_lorentz[ind] = dlorentz;

        /* HACK! Setting *conservative B* to constant, not primitive 
           Same as using B^i = B_0 psi^{-6} */
        if ( Bvec != NULL ) {
           const int indx = CCTK_VECTGFINDEX3D (cctkGH, i, j, kk, 0);
           const int indy = CCTK_VECTGFINDEX3D (cctkGH, i, j, kk, 1);
           const int indz = CCTK_VECTGFINDEX3D (cctkGH, i, j, kk, 2);
           Bvec[indx] = dbvecx;
           Bvec[indy] = dbvecy;
           Bvec[indz] = dbvecz;
           Bnx[ind] = dbvecx * rootdetg; 
           Bny[ind] = dbvecy * rootdetg; 
           Bnz[ind] = dbvecz * rootdetg;
        }
        if ( divclean_psi != NULL )
           divclean_psi[ind] = 0.; 

	if ( dtau <= 0. ) {
	   CCTK_VInfo(CCTK_THORNSTRING,"dtau error at (x,y,z)=(%g,%g,%g). (drho,dlorentz,deps,dpress)=(%g,%g,%g,%g)",
		xx,yy,zz, drho, dlorentz, deps, dpress);
	}
	assert(dtau > 0.);
	assert(dD > 0.);
	assert(dtau*(dtau+2*dD) > ( g11*pow2(dJx)+g22*pow2(dJy)+g33*pow2(dJz) +
                          2.*g21*dJx*dJy + 2.*g31*dJx*dJz + 
			  2.*g32*dJy*dJz ) );
	assert(drho > 0.);
	assert(deps > 0.);
	assert(dpress >= 0.);
	assert(dlorentz >= 1.);
	assert(( g11*pow2(dvelx)+g22*pow2(dvely)+g33*pow2(dvelz) +
                          2.*g21*dvelx*dvely + 2.*g31*dvelx*dvelz + 
			  2.*g32*dvely*dvelz ) <= 1. );
	if(testsuite) // dump internal data
	{
		fprintf(testsuitefh, "%g %g %g %g %g %g %g %g %g %g %g %g %g %g\n", 
				xx, yy, zz, dpsi, drho, deps, dpress, dvelx,
				dvely, dvelz, dbvecx, dbvecy, dbvecz, psiL);
	}

	  }
    }
  }

  if(verbose >= 2) {

    CCTK_VInfo (CCTK_THORNSTRING, 
        "Total matter linear momentum P^i = (%g, %g, %g)",
        linmom[0],linmom[1],linmom[2]);
    CCTK_VInfo (CCTK_THORNSTRING, 
        "Total matter angular momentum J^i = (%g, %g, %g)", 
        angmom[0],angmom[1],angmom[2]);

  }


  if(testsuitefh)
    fclose(testsuitefh);

}
void MatterLibs_FeedRho ( CCTK_ARGUMENTS ) {

  // Given the data fed to TwoPunctures and the resulting conformal factor, 
  // calculate the actual data.  Then fill in Whisky's variables.

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  enum eos_type_enum {eos_invalid = 0, eos_polytype, eos_general};
  enum initial_guess_type_enum {guess_invalid = 0, guess_rescaled_E, guess_rho, guess_rescaled_rho};
  enum eos_type_enum eos_type = eos_invalid;
  enum initial_guess_type_enum initial_guess_type = guess_invalid;
  const char * eos_type_name;
  double dpsi, drho, deps, dpress, dvelx, dvely, dvelz; // primitive initial data
  double dE, dD, dtau, dJx, dJy, dJz;                   // conservative initial data
  double dlorentz2, dlorentz, dh, B2, vB;		// auxiliary variables
  double xx, yy, zz;
  double psiL, chiL, rootdetg;
  double twelfth, third;
  double g11, g21, g31, g22, g32, g33, detg;
  double uxx, uxy, uxz, uyy, uyz, uzz;
  static double linmom[3];                              // total linear momentum of matter
  static double angmom[3];                              // total angular momentum of matter 
  FILE *testsuitefh=NULL;
  static int firsttime = 1;				// chain output together 

  // check whether we have a B-field and set pointers accordingly
  int have_uniform_bfield=0;
  double dbvecx, dbvecy, dbvecz; // primitive B-field initial data (conservative is trivial)
  if ( b0_uniform[0] != 0 || b0_uniform[1] != 0 || b0_uniform[2] != 0 ) {
     have_uniform_bfield=1;
  }
     

  // open file to dump internal state to
  if(testsuite)
  {
    char fn[1024];

    snprintf(fn, sizeof(fn)/sizeof(fn[0]), "%s/matterlibs_internaldata.asc", out_dir);
    if((testsuitefh = fopen(fn, firsttime ? "w" : "a")) == NULL)
      CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,  "Cannot open file %s to dump internal data", fn);

    if(firsttime) {
       if ( have_uniform_bfield ) {
          fprintf(testsuitefh, "# 1:x 2:y 3:z 4:psi 5:rho 6:eps 7:press 8:velx 9:vely 10:velz 11:bx 12:by 13:bz 14:psi_final\n");
       } else {
          fprintf(testsuitefh, "# 1:x 2:y 3:z 4:psi 5:rho 6:eps 7:press 8:velx 9:vely 10:velz 11:psi_final\n");
       }
    }

    firsttime = 0;
  }

  // find out what kind of EOS we are using
  if (CCTK_EQUALS(matterlibs_eos_type, ""))
      eos_type_name = whisky_eos_type;
  else
      eos_type_name = matterlibs_eos_type;
  if (CCTK_EQUALS(eos_type_name, "Polytype")) {
      eos_type = eos_polytype;
      assert(*whisky_polytrope_handle>=0);
  } else if (CCTK_EQUALS(eos_type_name, "General")) {
      eos_type = eos_general;
      assert(*whisky_eos_handle>=0);
  } else {
      CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,  "Unknown EOS type %s is not 'Polytype' or 'General'.", eos_type_name);
  }
  assert(eos_type);

  // decide how to provide initial guesses for the E_to_rho solver
  if (CCTK_EQUALS(initial_guess, "rescaled E")) {
      initial_guess_type = guess_rescaled_E;
  } else if (CCTK_EQUALS(initial_guess, "rho")) {
      initial_guess_type = guess_rho;
  } else if (CCTK_EQUALS(initial_guess, "rescaled rho")) {
      initial_guess_type = guess_rescaled_rho;
  } else {
      CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,  "Unknown initial guess type '%s'.", initial_guess);
  }
  assert(initial_guess_type);

  // Whisky kindly sets the initial rho_min in Whisky_Rho_Minima_Setup

  twelfth = 1.0/12.0;
  third = 1.0/3.0;

  CCTK_REAL *CarpetWeights = (CCTK_REAL *)(CCTK_VarDataPtr(cctkGH, 0, "CarpetReduce::weight"));
  assert(CarpetWeights);
  for ( int kk = 0 ; kk < cctk_lsh[2] ; kk++) {
    for ( int j = 0 ; j < cctk_lsh[1] ; j++) {
      for (int i = 0 ; i < cctk_lsh[0] ; i++) {

        const int ind = CCTK_GFINDEX3D (cctkGH, i, j, kk);

        xx = x[ind];
        yy = y[ind];
        zz = z[ind];

        g11 = gxx[ind];
        g22 = gyy[ind];
        g33 = gzz[ind];
        g21 = gxy[ind];
        g31 = gxz[ind];
        g32 = gyz[ind];
        detg  =  2*g21*g31*g32 + g33*(g11*g22 - pow2(g21)) - g22*pow2(g31) - g11*pow2(g32);

	chiL = fmax(1/pow(detg,third),chiEps);
	// psiL is used in the rescaling from Ebar to E, it should not be
	// capped from above since that makes the initial data potentially
	// large. rootdetg on the other hand is only used to densitize the
	// conservative varibles at the end, here we want to cap it from above
	// to limit the magnitude of the ID
        // psiL = 1/pow(chiL,1.0/4.0); 
	rootdetg = 1/pow(chiL,3.0/2.0);
        psiL = pow(detg,twelfth);
	/*rootdetg = sqrt(detg);*/

        // obtain guesses for the primitive variables
        calcMatterRho( xx, yy, zz, &dpsi, &drho, &deps, &dpress, 
                        &dvelx, &dvely, &dvelz, &dbvecx, &dbvecy, &dbvecz );
        dpsi = combineInitialGuesses(guessPuncturePsi(xx, yy, zz), dpsi);
        if ( have_uniform_bfield > 0 && CCTK_EQUALS(initial_Bvec,"uniformBn") ) {
           /* divide by new psi since b0_uniform should be the *conserved* bfield */
           CCTK_REAL inv_rdetg=pow(dpsi,-6);
           dbvecx = dbvecx*inv_rdetg;
           dbvecy = dbvecy*inv_rdetg;
           dbvecz = dbvecz*inv_rdetg;
        }

	if(testsuite) // dump internal data
	{
          fprintf(testsuitefh, "%g %g %g %g %g %g %g %g %g %g %g %g %g %g\n", 
                  xx, yy, zz, dpsi, drho, deps, dpress, dvelx,
		  dvely, dvelz, dbvecx, dbvecy, dbvecz, psiL);
	}

	if( drho < *whisky_rho_min*(1.0+atmo_tolerance) ) // enforce at least atmospheric density
	{
	      drho = *whisky_rho_min;
	      dpress = EOS_Pressure(*whisky_polytrope_handle, drho, matterlibs_poison);
	      deps = EOS_SpecificIntEnergy(*whisky_polytrope_handle, drho, matterlibs_poison);
	      dvelx = dvely = dvelz = 0.;
	}

        // calculate auxiliary variables
        dlorentz2 = 1 / ( 1 - pow4(dpsi) * ( pow2(dvelx)+pow2(dvely)+pow2(dvelz) ) );
        dh = 1 + deps + dpress/drho;
        B2 = pow4(dpsi) * ( pow2(dbvecx) + pow2(dbvecy) + pow2(dbvecz) );
        vB = pow4(dpsi) * ( dvelx*dbvecx + dvely*dbvecy + dvelz*dbvecz );

        // Calculate linear and angular momentum through raw primitives
        if (CarpetWeights[ind] == 1.0 &&
             i >= cctk_nghostzones[0] &&  i < cctk_lsh[0]-cctk_nghostzones[0] &&
             j >= cctk_nghostzones[1] &&  j < cctk_lsh[1]-cctk_nghostzones[1] &&
            kk >= cctk_nghostzones[2] && kk < cctk_lsh[2]-cctk_nghostzones[2]) {

           // Calculate appropriate values for E and S_j after baring and unbaring 
           // as in gr-qc/0606104.
           dJx = ( (drho * dh * dlorentz2 + B2 ) * dvelx - vB*pow4(dpsi)*dbvecx ) * pow(dpsi/psiL, 10.);
           dJy = ( (drho * dh * dlorentz2 + B2 ) * dvely - vB*pow4(dpsi)*dbvecy ) * pow(dpsi/psiL, 10.);
           dJz = ( (drho * dh * dlorentz2 + B2 ) * dvelz - vB*pow4(dpsi)*dbvecz ) * pow(dpsi/psiL, 10.);

           linmom[0] += dJx*pow(psiL,10.)*CCTK_DELTA_SPACE(0)*CCTK_DELTA_SPACE(1)*CCTK_DELTA_SPACE(2); 
           linmom[1] += dJy*pow(psiL,10.)*CCTK_DELTA_SPACE(0)*CCTK_DELTA_SPACE(1)*CCTK_DELTA_SPACE(2); 
           linmom[2] += dJz*pow(psiL,10.)*CCTK_DELTA_SPACE(0)*CCTK_DELTA_SPACE(1)*CCTK_DELTA_SPACE(2); 
           angmom[0] += (yy*dJz-zz*dJy)*pow(psiL,10.)*CCTK_DELTA_SPACE(0)*CCTK_DELTA_SPACE(1)*CCTK_DELTA_SPACE(2);
           angmom[1] += (zz*dJx-xx*dJz)*pow(psiL,10.)*CCTK_DELTA_SPACE(0)*CCTK_DELTA_SPACE(1)*CCTK_DELTA_SPACE(2);
           angmom[2] += (xx*dJy-yy*dJx)*pow(psiL,10.)*CCTK_DELTA_SPACE(0)*CCTK_DELTA_SPACE(1)*CCTK_DELTA_SPACE(2);

        }
        // DEBUG

        // Calculate appropriate values for E and S_j after baring and unbaring 
        // as in gr-qc/0606104.
        dE = ( drho * dh * dlorentz2 - dpress + (1 - 0.5/dlorentz2)*B2 - 0.5*pow2(vB) ) * pow(dpsi/psiL, conformal_density_power);
        dJx = ( (drho * dh * dlorentz2 + B2) * dvelx - vB*pow4(dpsi)*dbvecx ) * pow(dpsi/psiL, 10.);
        dJy = ( (drho * dh * dlorentz2 + B2) * dvely - vB*pow4(dpsi)*dbvecy ) * pow(dpsi/psiL, 10.);
        dJz = ( (drho * dh * dlorentz2 + B2) * dvelz - vB*pow4(dpsi)*dbvecz ) * pow(dpsi/psiL, 10.);

        if ( ( eos_type == eos_polytype ) ) {

          // for a polytropic EOS D is linked to E and S^j

          // recover primitive values "by hand"
          switch(initial_guess_type) // provide initial guess to solver
          {
            case guess_rescaled_E:
              drho = dE;
              break;
            case guess_rho:
              drho = drho;
              break;
            case guess_rescaled_rho:
              drho = drho * pow(dpsi/psiL, conformal_density_power);
              break;
            case guess_invalid:
            default:
              assert(0);
              break;
          }
          int Etorho_err;
          if ( have_uniform_bfield == 0 || B2 < MINB2 ) {
               Etorho_err = matterlibs_E_to_rho(*whisky_polytrope_handle, dE, dJx, dJy, dJz, psiL,
                                      &drho, &dlorentz);
          } else {
               Etorho_err = matterlibs_E_to_rho_MHD(*whisky_polytrope_handle, dE, dJx, dJy, dJz, 
                              dbvecx, dbvecy, dbvecz, psiL, &drho, &dlorentz);
          }  

	  if( Etorho_err != 0 )
	  {
	      CCTK_VWarn(0,__LINE__,__FILE__,CCTK_THORNSTRING,
		      "initial recovery of primitive variables failed at "
		      "(%g,%g,%g), where (E,J^j)=[%g,%g,%g,%g], "
		      "psiL = %g, dpsi = %g => (rho,w_lorentz) = (%g,%g)",
		      xx,yy,zz, dE, dJx,dJy,dJz, psiL,dpsi, drho,dlorentz);
	      return; /* NOTREACHED */
	  }

	  if( drho < *whisky_rho_min*(1.0+atmo_tolerance) ) // enforce at least atmospheric density
	  {
	    drho = *whisky_rho_min;
	    dpress = EOS_Pressure(*whisky_polytrope_handle, drho, matterlibs_poison); // atmosphere always used a polytrope
	    deps = EOS_SpecificIntEnergy(*whisky_polytrope_handle, drho, matterlibs_poison);
	    dvelx = dvely = dvelz = 0.;
	    dlorentz2 = dlorentz = 1.;
            vB = pow4(dpsi) * ( dvelx*dbvecx + dvely*dbvecy + dvelz*dbvecz );

	    dD = drho * dlorentz;
	    dJx = dJy = dJz = 0.;
	    // hand crafted to guarantee positivity
	    dtau = drho*(dlorentz2-dlorentz) + drho*deps*dlorentz2 + dpress*(dlorentz2-1.) + (1 - 0.5/dlorentz2)*B2 - 0.5*pow2(vB);
	  }
	  else
	  {
	    dpress = EOS_Pressure(*whisky_polytrope_handle, drho, matterlibs_poison);
	    deps = EOS_SpecificIntEnergy(*whisky_polytrope_handle, drho, matterlibs_poison);
	    dlorentz2 = pow2(dlorentz);
            dh = 1 + deps + dpress/drho;
	    dvelx = ( dJx + vB*pow4(dpsi)*dbvecx ) / (drho * dh * dlorentz2 + B2);
	    dvely = ( dJy + vB*pow4(dpsi)*dbvecy ) / (drho * dh * dlorentz2 + B2);
	    dvelz = ( dJz + vB*pow4(dpsi)*dbvecz ) / (drho * dh * dlorentz2 + B2);

	    // calculate a suitable matter density (this is the reason we go 
	    // through the whole rigmalore...)
	    dD = drho * dlorentz;

	    // hand crafted to guarantee positivity
	    dtau = drho*(dlorentz2-dlorentz) + drho*deps*dlorentz2 + dpress*(dlorentz2-1.) + (1 - 0.5/dlorentz2)*B2 - 0.5*pow2(vB);

            // check result of conversion for consistency
            if(verbose >= 1 &&
              GSL_SUCCESS != matterlibs_gsl_root_test_interval(dE, drho*dh*dlorentz2-dpress+(1-0.5/dlorentz2)*B2 - 0.5*pow2(vB), 
                                                    abs_err, rel_err)) 
            {
              CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,  
                      "Back and forth conversion of relativistic energy density E failed: "
                      "E(rho) != E_orig; E(rho) = %g, rho = %g, E_orig = %g, diff = %e", 
                      drho*dh*dlorentz2-dpress +(1-0.5/dlorentz2)*B2 - 0.5*pow2(vB), 
                      drho, dE, dE - (drho*dh*dlorentz2-dpress+(1-0.5/dlorentz2)*B2-0.5*pow2(vB)));
            }
	  }
           
          // factor in the determinant that Whisky uses
          dJx *= rootdetg;
          dJy *= rootdetg;
          dJz *= rootdetg;
          dD *= rootdetg;
          dtau *= rootdetg;

        } else if ( (eos_type == eos_general) ) {
          // for general EOS D is independend of E and S^j, and since Roland 
          // does not better we choose a convenient D.

          dlorentz = sqrt(dlorentz2);
          dD = drho * dlorentz * pow(dpsi/psiL, conformal_density_power);
          dtau = dE - dD;

          // factor in the determinant that Whisky uses

          dJx *= rootdetg;
          dJy *= rootdetg;
          dJz *= rootdetg;
          dD *= rootdetg;
          dtau *= rootdetg;

          // Create primitives from conservatives
          UpperMet(&uxx, &uxy, &uxz, &uyy, &uyz, &uzz, 
                 detg, g11, g21, g31, g22, g32, g33);

          assert(0 == "Con2PrimGen is not publicly accessibly in whisky");
#if 0 
          Con2PrimGen(*whisky_eos_handle, &dD, 
                             &dJx, &dJy, &dJz, 
                             &dtau, 
                             &drho, &dvelx, &dvely, &dvelz, 
                             &deps, &dpress, 
                             &dlorentz, 
                             uxx, uxy, uxz, uyy, uyz, uzz, detg, 
                             xx, yy, zz, r[ind],
                             *matterlibs_rho_min);
#endif
              
        } else {
            CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,  "Invalid eos_type: %d", eos_type);
            return; /* NOTREACHED */

        }

        // fill in grid functions (Si = gij S^j!)
        dens[ind] = dD;
        sx[ind] = g11*dJx+g21*dJy+g31*dJz;
        sy[ind] = g21*dJx+g22*dJy+g32*dJz;
        sz[ind] = g31*dJx+g32*dJy+g33*dJz;
        tau[ind] = dtau;
        rho[ind] = drho;
#ifdef HAVE_MAYA_WHISKY
//miguekl        velx[ind] = dvelx;
        velx[ind] = dvelx;
        vely[ind] = dvely;
        velz[ind] = dvelz;
#else
        const int xidx = CCTK_VECTGFINDEX3D (cctkGH, i, j, kk, 0);
        const int yidx = CCTK_VECTGFINDEX3D (cctkGH, i, j, kk, 1);
        const int zidx = CCTK_VECTGFINDEX3D (cctkGH, i, j, kk, 2);
//miguel        vel[xidx] = dvelx;
        vel[xidx] = dvelx;
        vel[yidx] = dvely;
        vel[zidx] = dvelz;
#endif
        eps[ind] = deps;
        press[ind] = dpress;
        w_lorentz[ind] = dlorentz;

        /* HACK! Setting *conservative B* to constant, not primitive 
           Same as using B^i = B_0 psi^{-6} */
        if ( Bvec != NULL ) {
           const int indx = CCTK_VECTGFINDEX3D (cctkGH, i, j, kk, 0);
           const int indy = CCTK_VECTGFINDEX3D (cctkGH, i, j, kk, 1);
           const int indz = CCTK_VECTGFINDEX3D (cctkGH, i, j, kk, 2);
           Bvec[indx] = dbvecx;
           Bvec[indy] = dbvecy;
           Bvec[indz] = dbvecz;
           Bnx[ind] = dbvecx * rootdetg; 
           Bny[ind] = dbvecy * rootdetg; 
           Bnz[ind] = dbvecz * rootdetg;
        }
        if ( divclean_psi != NULL )
           divclean_psi[ind] = 0.; 

	if ( dtau <= 0. ) {
	   CCTK_VInfo(CCTK_THORNSTRING,"dtau error at (x,y,z)=(%g,%g,%g). (drho,dlorentz,deps,dpress)=(%g,%g,%g,%g)",
		xx,yy,zz, drho, dlorentz, deps, dpress);
	}
	assert(dtau > 0.);
	assert(dD > 0.);
	assert(dtau*(dtau+2*dD) > ( g11*pow2(dJx)+g22*pow2(dJy)+g33*pow2(dJz) +
                          2.*g21*dJx*dJy + 2.*g31*dJx*dJz + 
			  2.*g32*dJy*dJz ) );
	assert(drho > 0.);
	assert(deps > 0.);
	assert(dpress >= 0.);
	assert(dlorentz >= 1.);
	assert(( g11*pow2(dvelx)+g22*pow2(dvely)+g33*pow2(dvelz) +
                          2.*g21*dvelx*dvely + 2.*g31*dvelx*dvelz + 
			  2.*g32*dvely*dvelz ) <= 1. );

      }
    }
  }

  if(verbose >= 2) {

    CCTK_VInfo (CCTK_THORNSTRING, 
        "Total matter linear momentum P^i = (%g, %g, %g)",
        linmom[0],linmom[1],linmom[2]);
    CCTK_VInfo (CCTK_THORNSTRING, 
        "Total matter angular momentum J^i = (%g, %g, %g)", 
        angmom[0],angmom[1],angmom[2]);

  }


  if(testsuitefh)
    fclose(testsuitefh);
}

void MatterLibs_SetRho ( CCTK_ARGUMENTS ) {

  // Given the data fed to TwoPunctures and the resulting conformal factor, 
  // calculate the actual data.  Then fill in Whisky's variables.

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  double dpsi, drho, deps, dpress, dvelx, dvely, dvelz; // primitive initial data
  double dE, dD, dtau, dJx, dJy, dJz;                   // conservative initial data
  double dlorentz2, dlorentz, dh;        		// auxiliary variables
  double B2, vB, Blowx, Blowy, Blowz;                   // Bfield shorthands
  double xx, yy, zz;
  double chiL, rootdetg;
  double twelfth, third;
  double g11, g21, g31, g22, g32, g33, detg;

  // check whether we have a B-field
  int have_uniform_bfield=0;
  double dbvecx, dbvecy, dbvecz; // primitive B-field initial data (conservative is trivial)
  if ( b0_uniform[0] != 0 || b0_uniform[1] != 0 || b0_uniform[2] != 0 ) {
     have_uniform_bfield=1;
  }

  // Whisky kindly sets the initial rho_min in Whisky_Rho_Minima_Setup

  twelfth = 1.0/12.0;
  third = 1.0/3.0;

  for ( int kk = 0 ; kk < cctk_lsh[2] ; kk++) {
    for ( int j = 0 ; j < cctk_lsh[1] ; j++) {
      for (int i = 0 ; i < cctk_lsh[0] ; i++) {

        const int ind = CCTK_GFINDEX3D (cctkGH, i, j, kk);

        xx = x[ind];
        yy = y[ind];
        zz = z[ind];

        g11 = gxx[ind];
        g22 = gyy[ind];
        g33 = gzz[ind];
        g21 = gxy[ind];
        g31 = gxz[ind];
        g32 = gyz[ind];
        detg  =  2*g21*g31*g32 + g33*(g11*g22 - pow2(g21)) - g22*pow2(g31) - g11*pow2(g32);

	chiL = fmax(1/pow(detg,third),chiEps);
	rootdetg = 1/pow(chiL,3.0/2.0);

        // obtain guesses for the primitive variables
        calcMatterRho( xx, yy, zz, &dpsi, &drho, &deps, &dpress, 
                        &dvelx, &dvely, &dvelz, &dbvecx, &dbvecy, &dbvecz );
        if ( have_uniform_bfield > 0 && CCTK_EQUALS(initial_Bvec,"uniformBn") ) {
           /* divide by new psi since b0_uniform should be the *conserved* bfield */
           CCTK_REAL inv_rootdetg=pow(chiL,3.0/2.0);
           dbvecx = dbvecx*inv_rootdetg;
           dbvecy = dbvecy*inv_rootdetg;
           dbvecz = dbvecz*inv_rootdetg;
        }
        // NOTE: dpsi is never used
        dpsi = -1.0; // trash variable

	if( drho < *whisky_rho_min*(1.0+atmo_tolerance) ) // enforce at least atmospheric density
	{
	      drho = *whisky_rho_min;
	      dpress = EOS_Pressure(*whisky_polytrope_handle, drho, matterlibs_poison);
	      deps = EOS_SpecificIntEnergy(*whisky_polytrope_handle, drho, matterlibs_poison);
	      dvelx = dvely = dvelz = 0.;
	}
//        if (xx*xx+yy*yy+zz*zz <= 0.4) {
//        dvelx= 0.0;
//	dvely=dvelz=0.0;
//        }
//        else {
        dvelx= -0.3/sqrt(g11);
	dvely=dvelz=0;
//        }
        // calculate auxiliary variables
        dlorentz2 = 1 / ( 1 - ( g11*pow2(dvelx) + g22*pow2(dvely) + g33*pow2(dvelz)
                                + 2*g21*dvelx*dvely + 2*g31*dvelx*dvelz + 2*g32*dvely*dvelz ) );
        dlorentz = sqrt(dlorentz2);
        dh = 1 + deps + dpress/drho;
        Blowx = g11*dbvecx + g21*dbvecy + g31*dbvecz;
        Blowy = g21*dbvecx + g22*dbvecy + g32*dbvecz;
        Blowz = g31*dbvecx + g32*dbvecy + g33*dbvecz;
        B2 = Blowx*dbvecx + Blowy*dbvecy + Blowz*dbvecz;
        vB = g11*dbvecx*dvelx + g22*dbvecy*dvely + g33*dbvecz*dvelz
              + g21*(dbvecx*dvely+dvelx*dbvecy) + g31*(dbvecx*dvelz+dbvecz*dvelx)
              + g32*(dbvecy*dvelz+dbvecz*dvely);

        // Calculate appropriate values for E and S^j after baring and unbaring 
        // as in gr-qc/0606104.
        dE = ( drho * dh * dlorentz2 - dpress + (1. - 0.5/dlorentz2)*B2 - 0.5*pow2(vB) );
        dJx = ( (drho * dh * dlorentz2 + B2 ) * dvelx - vB*Blowx );
        dJy = ( (drho * dh * dlorentz2 + B2 ) * dvely - vB*Blowy );
        dJz = ( (drho * dh * dlorentz2 + B2 ) * dvelz - vB*Blowz );
        dD = drho * dlorentz;
        dtau = drho*(dlorentz2-dlorentz) + drho*deps*dlorentz2 + dpress*(dlorentz2-1.) + (1.-0.5/dlorentz2)*B2 - 0.5*pow2(vB);

        // factor in the determinant that Whisky uses
        dJx *= rootdetg;
        dJy *= rootdetg;
        dJz *= rootdetg;
        dD *= rootdetg;
        dtau *= rootdetg;

        // fill in grid functions (Si = gij S^j!)
        dens[ind] = dD;
        sx[ind] = g11*dJx+g21*dJy+g31*dJz;
        sy[ind] = g21*dJx+g22*dJy+g32*dJz;
        sz[ind] = g31*dJx+g32*dJy+g33*dJz;
        tau[ind] = dtau;
        rho[ind] = drho;
#ifdef HAVE_MAYA_WHISKY
//Miguel        velx[ind] = dvelx;
	velx[ind] = dvelx;
        vely[ind] = dvely;
        velz[ind] = dvelz;
        
#else   
        const int xidx = CCTK_VECTGFINDEX3D (cctkGH, i, j, kk, 0);
        const int yidx = CCTK_VECTGFINDEX3D (cctkGH, i, j, kk, 1);
        const int zidx = CCTK_VECTGFINDEX3D (cctkGH, i, j, kk, 2);
// miguel       vel[xidx] = dvelx;
        vel[xidx] = dvelx;
        vel[yidx] = dvely;
        vel[zidx] = dvelz;
#endif
        eps[ind] = deps;
        press[ind] = dpress;
        w_lorentz[ind] = dlorentz;

        if ( Bvec != NULL ) {
           const int indx = CCTK_VECTGFINDEX3D (cctkGH, i, j, kk, 0);
           const int indy = CCTK_VECTGFINDEX3D (cctkGH, i, j, kk, 1);
           const int indz = CCTK_VECTGFINDEX3D (cctkGH, i, j, kk, 2);
           Bvec[indx] = dbvecx;
           Bvec[indy] = dbvecy;
           Bvec[indz] = dbvecz;
        }
        if ( Bnx != NULL ) {
           Bnx[ind] = dbvecx * rootdetg; 
           Bny[ind] = dbvecy * rootdetg; 
           Bnz[ind] = dbvecz * rootdetg;
        }
        if ( divclean_psi != NULL ) 
           divclean_psi[ind] = 0.; 

	if ( dtau <= 0. ) {
	   CCTK_VInfo(CCTK_THORNSTRING,"dtau error at (x,y,z)=(%g,%g,%g). (drho,dlorentz,deps,dpress)=(%g,%g,%g,%g)",
		xx,yy,zz, drho, dlorentz, deps, dpress);
	}
	if ( isnan(dtau) ) {
	   CCTK_VInfo(CCTK_THORNSTRING,"dtau is NaN at (x,y,z)=(%g,%g,%g). (drho,dlorentz,deps,dpress)=(%g,%g,%g,%g), dv=(%g,%g,%g), w2=%g, g=(%g,%g,%g,%g,%g,%g)",
		xx,yy,zz, drho, dlorentz, deps, dpress,dvelx,dvely,dvelz, dlorentz2, g11, g22, g33, g21,g31,g32);
	}
	assert(dtau > 0.);
	assert(dD > 0.);
	assert(dtau*(dtau+2*dD) > ( g11*pow2(dJx)+g22*pow2(dJy)+g33*pow2(dJz) +
                          2.*g21*dJx*dJy + 2.*g31*dJx*dJz + 
			  2.*g32*dJy*dJz ) );
	assert(drho > 0.);
	assert(deps > 0.);
	assert(dpress >= 0.);
	assert(dlorentz >= 1.);
	assert(( g11*pow2(dvelx)+g22*pow2(dvely)+g33*pow2(dvelz) +
                          2.*g21*dvelx*dvely + 2.*g31*dvelx*dvelz + 
			  2.*g32*dvely*dvelz ) <= 1. );

      }
    }
  }

}






static void Internal_SetRhoMetricSummed ( 
    CCTK_REAL * CCTK_RESTRICT gxx,  CCTK_REAL * CCTK_RESTRICT gxy,  CCTK_REAL * CCTK_RESTRICT gxz,  
    CCTK_REAL * CCTK_RESTRICT gyy,  CCTK_REAL * CCTK_RESTRICT gyz,  CCTK_REAL * CCTK_RESTRICT gzz,  
    CCTK_REAL * CCTK_RESTRICT betax,  CCTK_REAL * CCTK_RESTRICT betay,  CCTK_REAL * CCTK_RESTRICT betaz,  
    CCTK_REAL * CCTK_RESTRICT alp,  
    CCTK_REAL * CCTK_RESTRICT rho,  CCTK_REAL * CCTK_RESTRICT eps,  CCTK_REAL * CCTK_RESTRICT press,  
    CCTK_REAL * CCTK_RESTRICT w_lorentz,  
    CCTK_REAL * CCTK_RESTRICT velx,  CCTK_REAL * CCTK_RESTRICT vely,  CCTK_REAL * CCTK_RESTRICT velz, 
    CCTK_REAL * CCTK_RESTRICT Bvec, 
    CCTK_REAL * CCTK_RESTRICT dens,  CCTK_REAL * CCTK_RESTRICT tau,
    CCTK_REAL * CCTK_RESTRICT sx,  CCTK_REAL * CCTK_RESTRICT sy,  CCTK_REAL * CCTK_RESTRICT sz, 
    CCTK_REAL * CCTK_RESTRICT Bnx,  CCTK_REAL * CCTK_RESTRICT Bny,  CCTK_REAL * CCTK_RESTRICT Bnz, CCTK_REAL * CCTK_RESTRICT divpsi, 
    CCTK_REAL * CCTK_RESTRICT x,  CCTK_REAL * CCTK_RESTRICT y,  CCTK_REAL * CCTK_RESTRICT z, 
    const CCTK_REAL t0, 
    const CCTK_REAL boost_x, const CCTK_REAL boost_y, const CCTK_REAL boost_z, 
    const CCTK_REAL offset_x, const CCTK_REAL offset_y, const CCTK_REAL offset_z, 
    const int set_metric, const int set_lapse, const int set_shift, const int set_matter,
    const CCTK_REAL rho_min, const CCTK_REAL atmo_tolerance, const CCTK_INT eos_handle,
    const cGH *cctkGH, const CCTK_INT cctk_lsh[]);
void MatterLibs_SetRhoMetricSummed (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTS;

  // Set extra pointers for B-field vars 
  CCTK_REAL * CCTK_RESTRICT bvec=Bvec;
  CCTK_REAL * CCTK_RESTRICT bnx=Bnx;
  CCTK_REAL * CCTK_RESTRICT bny=Bny;
  CCTK_REAL * CCTK_RESTRICT bnz=Bnz;
  CCTK_REAL * CCTK_RESTRICT divpsi=divclean_psi;

  if (set_curvature) {
    CCTK_REAL dt = CCTK_DELTA_TIME;
    CCTK_REAL *gxx_next, *gxy_next, *gxz_next, *gyy_next, *gyz_next, *gzz_next;
    CCTK_REAL *gxx_prev, *gxy_prev, *gxz_prev, *gyy_prev, *gyz_prev, *gzz_prev;

    // get variable pointers via VarDataPtr to avoid inheritance issues for now.
    assert(gxx_next = CCTK_VarDataPtr(cctkGH, 0, "CalcK::gxx_next"));
    assert(gxy_next = CCTK_VarDataPtr(cctkGH, 0, "CalcK::gxy_next"));
    assert(gxz_next = CCTK_VarDataPtr(cctkGH, 0, "CalcK::gxz_next"));
    assert(gyy_next = CCTK_VarDataPtr(cctkGH, 0, "CalcK::gyy_next"));
    assert(gyz_next = CCTK_VarDataPtr(cctkGH, 0, "CalcK::gyz_next"));
    assert(gzz_next = CCTK_VarDataPtr(cctkGH, 0, "CalcK::gzz_next"));

    assert(gxx_prev = CCTK_VarDataPtr(cctkGH, 0, "CalcK::gxx_prev"));
    assert(gxy_prev = CCTK_VarDataPtr(cctkGH, 0, "CalcK::gxy_prev"));
    assert(gxz_prev = CCTK_VarDataPtr(cctkGH, 0, "CalcK::gxz_prev"));
    assert(gyy_prev = CCTK_VarDataPtr(cctkGH, 0, "CalcK::gyy_prev"));
    assert(gyz_prev = CCTK_VarDataPtr(cctkGH, 0, "CalcK::gyz_prev"));
    assert(gzz_prev = CCTK_VarDataPtr(cctkGH, 0, "CalcK::gzz_prev"));

    // copy metric in ADMBase into the various CalcK copies
    // this is not quite correct, in theory I should construct g_ADM(Delta_t)
    // from the definition of the extrinsic curvature and beta
    for (int ind = 0 ; ind < cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2] ; ind++) {
      gxx_prev[ind] = gxx_next[ind] = gxx[ind];
      gxy_prev[ind] = gxy_next[ind] = gxy[ind];
      gxz_prev[ind] = gxz_next[ind] = gxz[ind];
      gyy_prev[ind] = gyy_next[ind] = gyy[ind];
      gyz_prev[ind] = gyz_next[ind] = gyz[ind];
      gzz_prev[ind] = gzz_next[ind] = gzz[ind];
    }


    // FIXME: terrible hack we set gxx_prev and gxx_next first and ASSUME that
    // beta^i,alp is the same for all three times 
    // FIXME: it might be better to use first use CalcK to compute the extrinsic
    // curvature of the boosted metric only, then merge boosted metric, extrinsic
    // curvature, lapse and shift with those already on the grid. This most
    // likely means that I have to copy CalcK's code into a new routine.
    Internal_SetRhoMetricSummed(gxx_prev, gxy_prev, gxz_prev, gyy_prev, gyz_prev, gzz_prev, 
        betax, betay, betaz, alp,
        NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
        x, y, z,
        -dt, boost_x, boost_y, boost_z, offset_x, offset_y, offset_z,
        1, 0, 0, 0, 
        *whisky_rho_min, atmo_tolerance, *whisky_polytrope_handle,
        cctkGH, cctk_lsh);
    Internal_SetRhoMetricSummed(gxx_next, gxy_next, gxz_next, gyy_next, gyz_next, gzz_next, 
        betax, betay, betaz, alp,
        NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
        x, y, z,
        +dt, boost_x, boost_y, boost_z, offset_x, offset_y, offset_z,
        1, 0, 0, 0,
        *whisky_rho_min, atmo_tolerance, *whisky_polytrope_handle,
        cctkGH, cctk_lsh);
  }
  Internal_SetRhoMetricSummed(gxx, gxy, gxz, gyy, gyz, gzz, betax, betay, betaz,  
      alp, rho, eps, press, w_lorentz, velx, vely, velz, bvec, dens, tau, sx, sy, sz, bnx, bny, bnz, divpsi,
      x, y, z,
      0.0, boost_x, boost_y, boost_z, offset_x, offset_y, offset_z,
      set_metric, set_lapse, set_shift, 1,
      *whisky_rho_min, atmo_tolerance, *whisky_polytrope_handle,
      cctkGH, cctk_lsh);
}

static void Internal_SetRhoMetricSummed ( 
    CCTK_REAL * CCTK_RESTRICT gxx,  CCTK_REAL * CCTK_RESTRICT gxy,  CCTK_REAL * CCTK_RESTRICT gxz,  
    CCTK_REAL * CCTK_RESTRICT gyy,  CCTK_REAL * CCTK_RESTRICT gyz,  CCTK_REAL * CCTK_RESTRICT gzz,  
    CCTK_REAL * CCTK_RESTRICT betax,  CCTK_REAL * CCTK_RESTRICT betay,  CCTK_REAL * CCTK_RESTRICT betaz,  
    CCTK_REAL * CCTK_RESTRICT alp,  
    CCTK_REAL * CCTK_RESTRICT rho,  CCTK_REAL * CCTK_RESTRICT eps,  CCTK_REAL * CCTK_RESTRICT press,  
    CCTK_REAL * CCTK_RESTRICT w_lorentz,  
    CCTK_REAL * CCTK_RESTRICT velx,  CCTK_REAL * CCTK_RESTRICT vely,  CCTK_REAL * CCTK_RESTRICT velz,
    CCTK_REAL * CCTK_RESTRICT Bvec,
    CCTK_REAL * CCTK_RESTRICT dens,  CCTK_REAL * CCTK_RESTRICT tau,
    CCTK_REAL * CCTK_RESTRICT sx,  CCTK_REAL * CCTK_RESTRICT sy,  CCTK_REAL * CCTK_RESTRICT sz,
    CCTK_REAL * CCTK_RESTRICT Bnx,  CCTK_REAL * CCTK_RESTRICT Bny,  CCTK_REAL * CCTK_RESTRICT Bnz, CCTK_REAL * CCTK_RESTRICT divclean_psi,
    CCTK_REAL * CCTK_RESTRICT x,  CCTK_REAL * CCTK_RESTRICT y,  CCTK_REAL * CCTK_RESTRICT z, 
    const CCTK_REAL t0, 
    const CCTK_REAL boost_x, const CCTK_REAL boost_y, const CCTK_REAL boost_z, 
    const CCTK_REAL offset_x, const CCTK_REAL offset_y, const CCTK_REAL offset_z, 
    const int set_metric, const int set_lapse, const int set_shift, const int set_matter,
    const CCTK_REAL rho_min, const CCTK_REAL atmo_tolerance, const CCTK_INT eos_handle,
    const cGH *cctkGH, const CCTK_INT cctk_lsh[])
{
  // fill in gravity and matter variables by averaging the metric

  double dpsi, drho, deps, dpress, dvelx, dvely, dvelz; // primitive initial data
  double dalp;
  double dE, dD, dtau, dJx, dJy, dJz;                   // conservative initial data
  double dlorentz2, dlorentz, dh;			// auxiliary variables
  double dbvecx, dbvecy, dbvecz, B2, vB;                // Bfield shorthands
  double xx, yy, zz;
  double rootdetg;
  double g[4][4], g_boosted[4][4], detg;
  double u[4], u_boosted[4]; // covariant four-velocity

  CCTK_REAL betax_low, betay_low, betaz_low;
  double u11, u12, u13, u22, u23, u33;
  double g11, g12, g13, g22, g23, g33;
  double lapse, shiftx, shifty, shiftz;
  double velx_low, vely_low, velz_low;
  double bvecx_low, bvecy_low, bvecz_low;

  // unit vector, boost factor, lorentz factor
  double n[4], boost_mag, gamma_boost;
  // transformation matrix (2.44)
  double lambda[4][4];
  
  // Lorentz boost equations from MTW 2.44ff
  boost_mag = sqrt(pow2(boost_x) + pow2(boost_y) + pow2(boost_z));
  gamma_boost = 1/sqrt(1-pow2(boost_mag));
  if(boost_mag > 0)
  {
    n[0] = 0.0; // n is a purely spatial vector, but I use four components for convenience of coding further down
    n[1] = boost_x/boost_mag;
    n[2] = boost_y/boost_mag;
    n[3] = boost_z/boost_mag;
  }
  else
    n[0] = n[1] = n[2] = n[3] = 1.0; // avoid 0*nan situations

  lambda[0][0] = gamma_boost;
  for ( int i = 1 ; i < 4 ; i++ ) {
    lambda[0][i] = lambda[i][0] = -boost_mag*gamma_boost*n[i];
    for ( int j = i ; j < 4 ; j++ ) {
      lambda[i][j] = lambda[j][i] = (gamma_boost-1)*n[i]*n[j] + (i==j ? 1. : 0.);
    }
  }

  for ( int kk = 0 ; kk < cctk_lsh[2] ; kk++) {
    for ( int j = 0 ; j < cctk_lsh[1] ; j++) {
      for (int i = 0 ; i < cctk_lsh[0] ; i++) {

        const int ind = CCTK_GFINDEX3D (cctkGH, i, j, kk);

        // boost and offset coordinates. (inverse trafo)
        xx = t0*lambda[0][1] + (x[ind]-offset_x) * lambda[1][1] + (y[ind]-offset_y) * lambda[2][1] + (z[ind]-offset_z) * lambda[3][1];
        yy = t0*lambda[0][2] + (x[ind]-offset_x) * lambda[1][2] + (y[ind]-offset_y) * lambda[2][2] + (z[ind]-offset_z) * lambda[3][2];
        zz = t0*lambda[0][3] + (x[ind]-offset_x) * lambda[1][3] + (y[ind]-offset_y) * lambda[2][3] + (z[ind]-offset_z) * lambda[3][3];

        // obtain guesses for the primitive variables
        calcMatterRho( xx, yy, zz, &dpsi, &drho, &deps, &dpress, 
                        &dvelx, &dvely, &dvelz, &dbvecx, &dbvecy, &dbvecz );
        dalp = calcTOVLapse(xx, yy, zz);
        // Whisky kindly sets the initial rho_min in Whisky_Rho_Minima_Setup
	if( drho < rho_min*(1.0+atmo_tolerance) ) // enforce at least atmospheric density
	{
	      drho = rho_min;
	      dpress = EOS_Pressure(eos_handle, drho, matterlibs_poison);
	      deps = EOS_SpecificIntEnergy(eos_handle, drho, matterlibs_poison);
	      dvelx = dvely = dvelz = 0.;
	}

        // compute unboosted metric diag(-alp^2, psi^4, psi^4, psi^4)
        g[0][0] = -pow2(dalp);
        for ( int ii = 1 ; ii < 4 ; ii++ ) {
          g[ii][0] = g[0][ii] = 0.;
          for ( int jj = ii ; jj < 4 ; jj++ ) {
            g[ii][jj] = g[jj][ii] = (ii==jj ? pow4(dpsi) : 0.);
          }
        }

        // boost metric (could be made faster since g is diagonal...)
        // forward trafo for covariant objects
        for ( int ii = 0 ; ii < 4 ; ii++ ) {
          for ( int jj = 0 ; jj < 4 ; jj++ ) {
            g_boosted[ii][jj] = 0.0;
            for ( int ii2 = 0 ; ii2 < 4 ; ii2++ ) {
              for ( int jj2 = 0 ; jj2 < 4 ; jj2++ ) {
                g_boosted[ii][jj] += g[ii2][jj2] * lambda[ii2][ii] * lambda[jj2][jj];
              }
            }
          }
        }

        // rho, press and eps (and Bvec) do not change since they are defined in the rest
        // frame of the matter fluid

        // u_i = W*v_i
        dlorentz = 1 / sqrt(1 - pow4(dpsi) * (pow2(dvelx) + pow2(dvely) + pow2(dvelz)));
        u[1] = dlorentz * pow4(dpsi) * dvelx;
        u[2] = dlorentz * pow4(dpsi) * dvely;
        u[3] = dlorentz * pow4(dpsi) * dvelz;
        u[0] = -dalp * sqrt(1 + (pow2(u[1]) + pow2(u[2]) + pow2(u[3]))/pow4(dpsi)); // from normalization and diagonal metric

        // boost four-velocity
        // forward trafo for covariant objects
        for ( int ii = 0 ; ii < 4 ; ii++ ) {
          u_boosted[ii] = 0.0;
          for ( int ii2 = 0 ; ii2 < 4 ; ii2++ ) {
            u_boosted[ii] += u[ii2] * lambda[ii2][ii];
          }
        }

        // merge metrics (re-use g array)
        betax_low = gxx[ind] * betax[ind] + gxy[ind] * betay[ind] + gxz[ind] * betaz[ind]; 
        betay_low = gxy[ind] * betax[ind] + gyy[ind] * betay[ind] + gyz[ind] * betaz[ind]; 
        betaz_low = gxz[ind] * betax[ind] + gyz[ind] * betay[ind] + gzz[ind] * betaz[ind]; 
        g[0][0] = (-pow2(alp[ind]) + betax_low*betax[ind] + betay_low*betay[ind] + betaz_low*betaz[ind]) + g_boosted[0][0] + 1.0;
        g[0][1] = g[1][0] = betax_low + g_boosted[0][1];
        g[0][2] = g[2][0] = betay_low + g_boosted[0][2];
        g[0][3] = g[3][0] = betaz_low + g_boosted[0][3];
        g[1][1] = gxx[ind] + g_boosted[1][1] - 1.0;
        g[1][2] = g[2][1] = gxy[ind] + g_boosted[1][2];
        g[1][3] = g[3][1] = gxz[ind] + g_boosted[1][3];
        g[2][2] = gyy[ind] + g_boosted[2][2] - 1.0;
        g[2][3] = g[3][2] = gyz[ind] + g_boosted[2][3];
        g[3][3] = gzz[ind] + g_boosted[3][3] - 1.0;

        // 3+1 decompose merged metric
        g11 = g[1][1];
        g12 = g[1][2];
        g13 = g[1][3];
        g22 = g[2][2];
        g23 = g[2][3];
        g33 = g[3][3];

        detg  =  2*g12*g13*g23 + g33*(g11*g22 - pow2(g12)) - g22*pow2(g13) - g11*pow2(g23);

        // inverse three-metric
        u11=(g22*g33-g23*g23)/detg;
        u12=(g13*g23-g33*g12)/detg;
        u13=(g12*g23-g22*g13)/detg;
        u22=(g11*g33-g13*g13)/detg;
        u23=(g12*g13-g11*g23)/detg;
        u33=(g11*g22-g12*g12)/detg;

        shiftx = u11*g[0][1] + u12*g[0][2] + u13*g[0][3];
        shifty = u12*g[0][1] + u22*g[0][2] + u23*g[0][3];
        shiftz = u13*g[0][1] + u23*g[0][2] + u33*g[0][3];
        
        lapse = sqrt(fabs((shiftx*g[0][1]+shifty*g[0][2]+shiftz*g[0][3])-g[0][0]));

        dlorentz2 = 1 + u11*pow2(u_boosted[1]) + 2*u12*u_boosted[1]*u_boosted[2] + 
                  2*u13*u_boosted[1]*u_boosted[3] + u22*pow2(u_boosted[2]) +
                  2*u23*u_boosted[2]*u_boosted[3] + u33*pow2(u_boosted[3]);
        dlorentz = sqrt(dlorentz2);
        velx_low = u_boosted[1]/dlorentz;
        vely_low = u_boosted[2]/dlorentz;
        velz_low = u_boosted[3]/dlorentz;


        // calculate auxiliary variables
	rootdetg = sqrt(detg);
        dh = 1 + deps + dpress/drho;
        bvecx_low = g11*dbvecx + g12*dbvecy + g13*dbvecz;
        bvecy_low = g12*dbvecx + g22*dbvecy + g23*dbvecz;
        bvecz_low = g13*dbvecx + g23*dbvecy + g33*dbvecz;
        B2 = bvecx_low*dbvecx + bvecy_low*dbvecy + bvecz_low*dbvecz;
        vB = velx_low*dbvecx + vely_low*dbvecy + velz_low*dbvecz;

        // Calculate appropriate values for E and S_j
        dE = ( drho * dh * dlorentz2 - dpress + (1-0.5/dlorentz2)*B2 - 0.5*pow2(vB) );
        dJx = ( (drho * dh * dlorentz2 + B2) * velx_low - vB*bvecx_low);
        dJy = ( (drho * dh * dlorentz2 + B2) * vely_low - vB*bvecy_low);
        dJz = ( (drho * dh * dlorentz2 + B2) * velz_low - vB*bvecz_low);
        dD = drho * dlorentz;
        dtau = drho*(dlorentz2-dlorentz) + drho*deps*dlorentz2 + dpress*(dlorentz2-1.) + (1-0.5/dlorentz2)*B2 - 0.5*pow2(vB);

        // factor in the determinant that Whisky uses
        dJx *= rootdetg;
        dJy *= rootdetg;
        dJz *= rootdetg;
        dD *= rootdetg;
        dtau *= rootdetg;

        // fill in grid functions (vel[xyz] = v^i)
        if ( set_matter ) {
          dens[ind] = dD;
          sx[ind] = dJx;
          sy[ind] = dJy;
          sz[ind] = dJz;
          tau[ind] = dtau;
          rho[ind] = drho;
          CCTK_REAL dvelx = u11*velx_low + u12*vely_low + u13*velz_low;
          CCTK_REAL dvely = u12*velx_low + u22*vely_low + u23*velz_low;
          CCTK_REAL dvelz = u13*velx_low + u23*vely_low + u33*velz_low;
#ifdef HAVE_MAYA_WHISKY
//          velx[ind] = dvelx;
          velx[ind] = dvelx;
 
          vely[ind] = dvely;
          velz[ind] = dvelz;
#else
          const int xidx = CCTK_VECTGFINDEX3D (cctkGH, i, j, kk, 0);
          const int yidx = CCTK_VECTGFINDEX3D (cctkGH, i, j, kk, 1);
          const int zidx = CCTK_VECTGFINDEX3D (cctkGH, i, j, kk, 2);
//miguel          vel[xidx] = dvelx;
          vel[xidx] = dvelx;
          vel[yidx] = dvely;
          vel[zidx] = dvelz;
#endif
          eps[ind] = deps;
          press[ind] = dpress;
          w_lorentz[ind] = dlorentz;
          if ( Bvec != NULL ) {
             const int indx = CCTK_VECTGFINDEX3D (cctkGH, i, j, kk, 0);
             const int indy = CCTK_VECTGFINDEX3D (cctkGH, i, j, kk, 1);
             const int indz = CCTK_VECTGFINDEX3D (cctkGH, i, j, kk, 2);
             Bvec[indx] = dbvecx;
             Bvec[indy] = dbvecy;
             Bvec[indz] = dbvecz;
             Bnx[ind] = dbvecx * rootdetg; 
             Bny[ind] = dbvecy * rootdetg; 
             Bnz[ind] = dbvecz * rootdetg; 
          }
          if ( divclean_psi != NULL ) 
             divclean_psi[ind] = 0.; 

          assert(tau[ind] > 0.);
          assert(dens[ind] > 0.);
          assert(tau[ind]*(tau[ind]+2*dens[ind]) > ( u11*pow2(sx[ind])+u22*pow2(sy[ind])+u33*pow2(sz[ind]) +
                            2.*u12*sx[ind]*sy[ind] + 2.*u13*sx[ind]*sz[ind] + 
                            2.*u23*sy[ind]*sz[ind] ) );
          assert(rho[ind] > 0.);
          assert(eps[ind] > 0.);
          assert(press[ind] >= 0.);
          assert(w_lorentz[ind] >= 1.);
          assert(( gxx[ind]*pow2(dvelx)+gyy[ind]*pow2(dvely)+gzz[ind]*pow2(dvelz) +
                            2.*gxy[ind]*dvelx*dvely + 2.*gxz[ind]*dvelx*dvelz + 
                            2.*gyz[ind]*dvely*dvelz ) <= 1. );
        }

        if ( set_metric ) {
          gxx[ind] = g11;
          gyy[ind] = g22;
          gzz[ind] = g33;
          gxy[ind] = g12;
          gxz[ind] = g13;
          gyz[ind] = g23;
        }
        if ( set_shift ) {
          betax[ind] = shiftx;
          betay[ind] = shifty;
          betaz[ind] = shiftz;
        }
        if ( set_lapse ) {
          alp[ind] = lapse;
        }

      }
    }
  }

}

void MatterLibs_SetLapse ( CCTK_ARGUMENTS ) {

  // Provides a new option for initial_lapse "tov"

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_REAL xx,yy,zz;

  for ( int kk=0; kk < cctk_lsh[2]; kk++) {
    for ( int j=0; j < cctk_lsh[1]; j++) {
      for (int i=0; i < cctk_lsh[0]; i++) {

        const int ind = CCTK_GFINDEX3D (cctkGH, i, j, kk);

        xx = x[ind];
        yy = y[ind];
        zz = z[ind];

        // lapse for a stationary TOV star
        alp[ind] = calcTOVLapse(xx, yy, zz);
      }
    }
  }
}

void MatterLibs_SetMetric ( CCTK_ARGUMENTS ) {

  // Provides a new option for initial_data "tov"

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_REAL xx,yy,zz;
  CCTK_REAL dpsi, drho, deps, dpress, dvelx, dvely, dvelz;
  CCTK_REAL g_diag;

  // check whether we have a B-field and set pointers accordingly
  int have_uniform_bfield=0;
  CCTK_REAL dbvecx, dbvecy, dbvecz;

  for ( int kk=0; kk < cctk_lsh[2]; kk++) {
    for ( int j=0; j < cctk_lsh[1]; j++) {
      for (int i=0; i < cctk_lsh[0]; i++) {

        const int ind = CCTK_GFINDEX3D (cctkGH, i, j, kk);

        xx = x[ind];
        yy = y[ind];
        zz = z[ind];

        // obtain guesses for the primitive variables, just to get psi
        calcMatterRho( xx, yy, zz, &dpsi, &drho, &deps, &dpress, 
                        &dvelx, &dvely, &dvelz, &dbvecx, &dbvecy, &dbvecz );
        g_diag = pow4(dpsi);

        // metric and extrinsice curvature (if beta=0) for a stationary TOV star
        gxx[ind] = gyy[ind] = gzz[ind] = g_diag;
        gxy[ind] = gxz[ind] = gyz[ind] = 0.0;
        kxx[ind] = kyy[ind] = kzz[ind] = 0.0;
        kxy[ind] = kxz[ind] = kyz[ind] = 0.0;
      }
    }
  }
}
