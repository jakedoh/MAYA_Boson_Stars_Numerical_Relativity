/* Main functions for converting conservatives
   to primitives */

#include <assert.h>
#include <cstdlib>
#include <stdio.h>
#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "EOS_Base.h"
#include "RecoverMHD.h"

// In Whisky (probably has a aliased function as well...)
extern "C" {
void CCTK_FCALL CCTK_FNAME(Prim2ConMHD)(int *handle, 
 double* gxx, double* gxy, double* gxz, double* gyy, double* gyz, double* gzz, double* det, double* ddens, 
 double* dsx, double* dsy, double* dsz, double* dtau , double* drho, double* dvelx, double* dvely, double* dvelz, double* deps, double* dpress, double* w, 
 double* bnx, double* bny, double* bnz, double* bvecx, double* bvecy, double* bvecz);
}

// usual max function
static inline double max(double a, double b) {
  return a>b ? a : b;
}

// how many tests to run?
#define NUM_TESTS 1000

// external routine, test con2prim
extern "C"
void RecoverMHD_Test(CCTK_ARGUMENTS);

// test con2prim with random data
extern "C"
void RecoverMHD_Test(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  // type of EOS to use
  const int polytype = CCTK_Equals(whisky_eos_type,"Polytype");

  // noise for initial guess (0.2 corresponds to correct value +/-0.1)
  const double noise = 0.2;

  // how good are we?
  double maxerr;

  int ierr; // found an error in con2prim?

  // send in

  // returned 

  // temps for vel (no-constrained)
  double uxL, uyL, uzL, w;
  double uxR, uyR, uzR;

  /* Create handle for solver type */
  typedef enum GSL_SOLVER_TYPE gsolve_type;
  gsolve_type stype;
  if ( CCTK_EQUALS(solver_type,"Newton") )
	stype=GST_Newton;
  else if ( CCTK_EQUALS(solver_type,"gNewton") )
	stype=GST_gNewton;
  else if ( CCTK_EQUALS(solver_type,"HybridSJ") )
	stype=GST_HybridSJ;
  else if ( CCTK_EQUALS(solver_type,"HybridJ") )
	stype=GST_HybridJ;
  else if ( CCTK_EQUALS(solver_type,"HybridS") )
	stype=GST_HybridS;
  else if ( CCTK_EQUALS(solver_type,"Hybrid") )
	stype=GST_Hybrid;
  else if ( CCTK_EQUALS(solver_type,"dNewton") )
	stype=GST_dNewton;
  else if ( CCTK_EQUALS(solver_type,"Broyden") )
	stype=GST_Broyden;
  else {
	CCTK_VWarn(0,__LINE__,__FILE__,CCTK_THORNSTRING,"Unknown solver type '%s'",solver_type);
	return; /* NOTREACHED */
  }

  /* Handle choice of system of eqns for recovery here to avoid */
  enum SYSTEM_EQN_TYPE{ HD_ONLY=0, MHD_5D_UiEps=1, MHD_2D_WP=2, MHD_1D_hW2=3, MHD_2D_hv2=4 } solve_system_method;
  if ( CCTK_Equals(system_type, "u_i eps") ) {
     solve_system_method = MHD_5D_UiEps;
  } else if (CCTK_Equals(system_type, "W press") ) {
     solve_system_method = MHD_2D_WP;
  } else if (CCTK_Equals(system_type, "h*W^2") ) {
     solve_system_method = MHD_1D_hW2;
  } else if (CCTK_Equals(system_type, "h v^2") ) {
     solve_system_method = MHD_2D_hv2;
  } else {
     CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING, "unknown system type '%s'", system_type);
     return; /* NOTREACHED */
  }

  // initialize random number generator, error accumulator
  srand48(42);
  maxerr = 0.;

  for(int testno = 0 ; testno < NUM_TESTS ; testno++)
  {

    // Randomly generated input
    double gxxL, gxyL, gxzL, gyyL, gyzL, gzzL, detg;
    double rhoL, velxL, velyL, velzL, pressL, w_lorentzL, epsL,
           BvecxL, BvecyL, BveczL;

    // Generated conservatives
    double densL, sxL, syL, szL, s2L, BnxL, BnyL, BnzL, tauL;
    double v2L, rho_enthalpy_W2L;

    // Returned primitives
    double rhoR, velxR, velyR, velzR, pressR, w_lorentzR, epsR,
          BvecxR, BvecyR, BveczR;

    // try flat space only for now
    gxxL = gyyL= gzzL = 1.0;
    gxyL = gxzL = gyzL = 0.0;
    detg = 1.0;

    // we choose u^i instead of vel^i since it is unconstrained
    uxL = 10*(drand48()-0.5);
    uyL = 15*(drand48()-0.5);
    uzL =  5*(drand48()-0.5);
    w = sqrt(1+gxxL*uxL*uxL+gyyL*uyL*uyL+gzzL*uzL*uzL); // flat

    rhoL = drand48();
    velxL = uxL/w; // assumes lapse=1,shift=0
    velyL = uyL/w;
    velzL = uzL/w;
    epsL = 0.3*drand48(); // internal energy, if too large then recovery fails

    BvecxL = drand48();
    BvecyL = drand48();
    BveczL = drand48();

    // if we really care we could output to some file instead...
    printf("Test no. %d:\n"
           "  rho  = %g\n"
           "  velx = %g\n"
           "  vely = %g\n"
           "  velz = %g\n"
           "  ux = %g\n"
           "  uy = %g\n"
           "  uz = %g\n"
           "  eps  = %g\n"
           "  Bvecx = %g\n"
           "  Bvecy = %g\n"
           "  Bvecz = %g\n",
           testno, rhoL, velxL, velyL, velzL, uxL, uyL, uzL, epsL, BvecxL, BvecyL, BveczL);

    CCTK_FNAME(Prim2ConMHD)(whisky_eos_handle, 
        &gxxL, &gxyL, &gxzL, &gyyL, &gyzL, &gzzL, &detg,
        &densL, &sxL, &syL, &szL, &tauL, 
        &rhoL, &velxL, &velyL, &velzL, &epsL, &pressL, &w_lorentzL, 
        &BnxL, &BnyL, &BnzL,
        &BvecxL, &BvecyL, &BveczL);
    s2L = gxxL*sxL*sxL+gyyL*syL*syL+gzzL*szL*szL; // flat
    v2L = gxxL*velxL*velxL+gyyL*velyL*velyL+gzzL*velzL*velzL; // flat
    rho_enthalpy_W2L = (rhoL*(1.+epsL)+pressL)*SQR(w_lorentzL);

    printf("Conservatives:\n"
           "  dens      = %g\n"
           "  sx        = %g\n"
           "  sy        = %g\n"
           "  sz        = %g\n"
           "  s2        = %g\n"
           "  tau       = %g\n"
           "  press     = %g\n"
           "  w_lorentz = %g\n"
           "  v2        = %g\n"
           "  rho h W^2 = %g\n"
           "  Bnx       = %g\n"
           "  Bny       = %g\n"
           "  Bnz       = %g\n",
           densL, sxL, syL, szL, s2L, tauL, pressL, w_lorentzL, v2L, rho_enthalpy_W2L, BnxL, BnyL, BnzL);

    printf("Solving...");
    // make up some random initial guesses :-)
    uxR = (1.0+noise*(drand48()-0.5))*uxL;
    uyR = (1.0+noise*(drand48()-0.5))*uyL;
    uzR = (1.0+noise*(drand48()-0.5))*uzL;
    w = sqrt(1+gxxL*uxR*uxR+gyyL*uyR*uyR+gzzL*uzR*uzR); // flat

    rhoR = (1.0+noise*(drand48()-0.5))*rhoL;
    velxR = uxR/w; // assumes lapse=1,shift=0
    velyR = uyR/w;
    velzR = uzR/w;
    pressR = (1.0+noise*(drand48()-0.5))*pressL;
    w_lorentzR = (1.0+noise*(drand48()-0.5))*w_lorentzL;
    epsR = (1.0+noise*(drand48()-0.5))*epsL;
    BvecxR = (1.0+noise*(drand48()-0.5))*BvecxL;
    BvecyR = (1.0+noise*(drand48()-0.5))*BvecyL;
    BveczR = (1.0+noise*(drand48()-0.5))*BveczL;

#if 0
    // hack: test computation of primitives after solve
    rhoR = rhoL;
    velxR = velxL;
    velyR = velyL;
    velzR = velzL;
    pressR = pressL;
    w_lorentzR = w_lorentzL;
    epsR = epsL;
    BvecxR = BvecxL;
    BvecyR = BvecyL;
    BveczR = BveczL;
#endif

    ierr = RecoverMHD_Con2Prim_pointwise( *whisky_eos_handle, *whisky_polytrope_handle, polytype, cctk_iteration, stype,
                                (int) solve_system_method, &densL, &tauL, &sxL, &syL, &szL, &BnxL, &BnyL, &BnzL,
                                &rhoR, &epsR, &velxR, &velyR, &velzR, &BvecxR, &BvecyR, &BveczR,
                                &w_lorentzR, &pressR, 0.0, 0.0, 0.0,
                                gxxL, gxyL, gxzL, gyyL, gyzL, gzzL,
                                *whisky_reflevel, 1.0, 
                                false, *whisky_rho_min, *whisky_tau_min);
    assert(!ierr);

    w = 1/sqrt(1-(gxxL*SQR(velxR)+gyyL*SQR(velyR)+gzzL*SQR(velzR)));

    printf("Found results:\n"
           "  rho  = %g [%g err]\n"
           "  velx = %g [%g err]\n"
           "  vely = %g [%g err]\n"
           "  velz = %g [%g err]\n"
           "  ux   = %g [%g err]\n"
           "  uy   = %g [%g err]\n"
           "  uz   = %g [%g err]\n"
           "  eps  = %g [%g err]\n"
           "  Bvecx = %g [%g err]\n"
           "  Bvecy = %g [%g err]\n"
           "  Bvecz = %g [%g err]\n",
           rhoR, fabs(rhoR-rhoL), 
           velxR,fabs(velxR-velxL),  
           velyR,fabs(velyR-velyL),  
           velzR,fabs(velzR-velzL),  
           velxR*w,fabs(velxR*w-uxL),  
           velyR*w,fabs(velyR*w-uyL),  
           velzR*w,fabs(velzR*w-uzL),  
           epsR, fabs(epsR-epsL),
           BvecxR,fabs(BvecxR-BvecxL),  
           BvecyR,fabs(BvecyR-BvecyL),  
           BveczR,fabs(BveczR-BveczL)  
          );

    // keep track of largest error
    maxerr = max(maxerr, fabs(rhoR-rhoL));
    maxerr = max(maxerr, fabs(velxR-velxL));
    maxerr = max(maxerr, fabs(velyR-velyL));
    maxerr = max(maxerr, fabs(velzR-velzL));
    maxerr = max(maxerr, fabs(epsR-epsL));
    maxerr = max(maxerr, fabs(BvecxR-BvecxL));
    maxerr = max(maxerr, fabs(BvecyR-BvecyL));
    maxerr = max(maxerr, fabs(BveczR-BveczL));
  }
  printf("Ran %d tests, max error: %g\n", NUM_TESTS, maxerr);
  exit(0);
}


