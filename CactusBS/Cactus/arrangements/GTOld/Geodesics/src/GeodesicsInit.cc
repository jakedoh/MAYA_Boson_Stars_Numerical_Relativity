#include "Geodesics.hh"
#include "MersenneTwister.h" // For high quality/fast uniform rand()
#include <math.h>

//double mb_dist(double temp);

extern "C" void Geodesics_Init(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  int lsh[1];
  int const ierr = CCTK_GrouplshGN (cctkGH, 1, lsh, "Geodesics::geo_pos");
  int lnum_geo = lsh[0];

  if( verbose > 2)
    CCTK_VInfo(CCTK_THORNSTRING,
        "Geodesics_Init(): Geodesics initializing");

  MTRand mtrand;

  for( int i=0; i < lnum_geo; i++ )  {

    // Initialize positions (TODO: read in csv's and old snapshots)
    if(CCTK_Equals(dist,"custom")){
      geo_x[i] = geo_x0[i];
      geo_y[i] = geo_y0[i];
      geo_z[i] = geo_z0[i];
    }

    if(CCTK_Equals(dist,"sphere")) { // The default
      double randx, randy, randz, randr;
      do {
        randx = sphere_outerr - mtrand.rand(2*sphere_outerr);
        randy = sphere_outerr - mtrand.rand(2*sphere_outerr);
        randz = sphere_outerr - mtrand.rand(2*sphere_outerr);
        randr = sqrt(randx*randx+randy*randy+randz*randz);
      } while(randr < sphere_innerr || randr > sphere_outerr );

      geo_x[i] = randx;
      geo_y[i] = randy;
      geo_z[i] = randz;
    }

    if(CCTK_Equals(dist,"annulus")) {
      double randx, randy, randz, randr;
      do {
        randx = 25.0 - mtrand.rand(50.0);
        randy = 25.0 - mtrand.rand(50.0);
        randz = 10.0 - mtrand.rand(20);
        randr = sqrt(randx*randx+randy*randy);//+randz*randz); //(Uncomment for sph r instead of polar)
      } while(randr < 8.0 || randr > 25.0
          || randz > 10.0|| randz < -10.0);

      geo_x[i] = randx;
      geo_y[i] = randy;
      geo_z[i] = randz;
    }

    if(CCTK_Equals(dist,"raytrace")) {

      CCTK_REAL sqrtnumgeo = (int)sqrt(num_geo);

      CCTK_INT gi = i + lnum_geo * CCTK_MyProc(cctkGH);

      geo_x[i] = 10.0 - 20.0 * fmod(gi,sqrtnumgeo)/sqrtnumgeo;
      geo_y[i] = 10.0;// - 20.0 * (int)(gi/sqrtnumgeo)/sqrtnumgeo; //-10.0;
      geo_z[i] = 10.0 - 20.0 * (int)(gi/sqrtnumgeo)/sqrtnumgeo;
    }

    // Initialize velocities
    if(CCTK_Equals(pdist,"custom")){
      geo_ux[i] = geo_px0[i];
      geo_uy[i] = geo_py0[i];
      geo_uz[i] = geo_pz0[i];
    }

    if(CCTK_Equals(pdist,"zero")) { //the default
      geo_ux[i] = 0.0;
      geo_uy[i] = 0.0;
      geo_uz[i] = 0.0;
    }

    if(CCTK_Equals(pdist,"raytrace")) {
      geo_ux[i] = 0.0;
      geo_uy[i] = -1.0;
      geo_uz[i] = 0.0; //-1.0;
    }

    if(CCTK_Equals(pdist,"kepler")) {
      CCTK_REAL geo_v = sqrt(2.0/sqrt(geo_x[i]*geo_x[i]+geo_y[i]*geo_y[i]+0.000001));
      CCTK_REAL geo_phi = atan2((double)geo_y[i], (double)geo_x[i]);

      geo_ux[i] = geo_v*cos(geo_phi);
      geo_uy[i] = geo_v*sin(geo_phi);
      geo_uz[i] = 0.0;
    }

    if(CCTK_Equals(pdist,"rand")) {
      double randpx, randpy, randpz, magp;
      do {
        randpx = 1.0 - mtrand.rand(2.0);
        randpy = 1.0 - mtrand.rand(2.0);
        randpz = 1.0 - mtrand.rand(2.0);
        magp = sqrt(randpx*randpx + randpy*randpy + randpz*randpz);
      } while(magp > 1.0);

      geo_ux[i] = randpx;
      geo_uy[i] = randpy;
      geo_uz[i] = randpz;
    }

    // Initialize charge and mass in cgs planck units 
    if(CCTK_Equals(geo_charge,"zero")){  
      geo_q[i] = 0.0;
      geo_m[i] = 1.0;
    }

    if(CCTK_Equals(geo_charge,"e")) {
      geo_q[i] = -1.38047 * pow(10.0,-34);
      geo_m[i] = 6.76320 * pow(10.0,-56);
    }

    if(CCTK_Equals(geo_charge,"eplus")) {
      geo_q[i] = 1.38047 * pow(10.0,-34);
      geo_m[i] = 6.76320 * pow(10.0,-56);
    }

    if(CCTK_Equals(geo_charge,"p")){
      geo_q[i] = 1.38047 * pow(10.0,-34);
      geo_m[i] = 1.24183 * pow(10.0,-52);
    }

    if(CCTK_Equals(geo_charge,"leptonic")){
      geo_q[i] = -1.38047 * pow(10.0,-34);
      if(i%2==0) geo_q[i] *= -1.0;
      geo_m[i] = 6.76320 * pow(10.0,-56);
    }

    if(CCTK_Equals(geo_charge,"hadronic")){
      if(i%2==0) {
        geo_q[i] = -1.38047 * pow(10.0,-34);
        geo_m[i] = 6.76320 * pow(10.0,-56);
      } else {
        geo_q[i] = 1.38047 * pow(10.0,-34);
        geo_m[i] = 1.24183 * pow(10.0,-52);
      }
    }

    if(CCTK_Equals(geo_charge,"mixed")){
      switch (i%3) {
        case 0:
          geo_q[i] = -1.38047 * pow(10.0,-34);
          geo_m[i] = 6.76320 * pow(10.0,-56);
          break;
        case 1:
          geo_q[i] = 1.38047 * pow(10.0,-34);
          geo_m[i] = 6.76320 * pow(10.0,-56);
          break;
        case 2:
          geo_q[i] = 1.38047 * pow(10.0,-34);
          geo_m[i] = 1.24183 * pow(10.0,-52);
          break;
      } // if we need a default case % has a problem (a big one)
    }

    // just some stuff so the first iterations output isnt poisoned values
    active_geo[i] = 0;
    geo_iut[i] = 1.0;
    if(CCTK_Equals(calc_sync,"yes")) geo_sync[i] = 1.0;
    if(CCTK_Equals(geo_moltest,"yes")) geo_t[i]=0.0;

/*    if(CCTK_Equals(pdist,"mb")) {
      double randgam, randpx, randpy, randpz, magp, pxhat, pyhat, pzhat, randmagp;

      randgam = mb_dist(mb_temp);

      randpx = 1.0 - mtrand.rand(2.0);
      randpy = 1.0 - mtrand.rand(2.0);
      randpz = 1.0 - mtrand.rand(2.0);
      magp = sqrt(randpx*randpx + randpy*randpy + randpz*randpz);

      pxhat = randpx/magp;
      pyhat = randpy/magp;
      pzhat = randpz/magp;
      randmagp = sqrt(1.0-1.0/(randgam*randgam));

      geo_ux[i] = randmagp*pxhat;
      geo_uy[i] = randmagp*pyhat;
      geo_uz[i] = randmagp*pzhat;
    } */

    if( verbose > 0 ) {
      CCTK_VInfo(CCTK_THORNSTRING,
          "Geodesics_Init(): Init particle %d @ x=%f,y=%f,z=%f",
          i,
          (double)geo_x[i],
          (double)geo_y[i],
          (double)geo_z[i]);
    }
  }

  Geodesics_Flag(CCTK_PASS_CTOC); // run through Flag() just in case
  //Geodesics_Analysis(CCTK_PASS_CTOC);
}

/*
double mb_dist(double temp) { //this is a matlab script I converted by hand (so it might be broken, it's untested)

  MTRand mt;
  double aa = 6.0*pow(10.0,9)/temp;               // temperature parameter
  double gend = 10.0/aa + 2.0;                   // end value of gamma
  double bess_k2 = 1.0;//gsl_sf_bessel_Kn(2,aa);//besselk(2,aa);        // TODO:modified bessel function of the second kind
  // see GSL: double gsl_sf_bessel_Kn (int n, double x)
  double nfact = aa/bess_k2;             // normalization factor

  // Compute distribution function

  int ng = 1000;                   // No. of grid points 
  double dg = (gend-1.0)/(ng-1);      // Grid spacing
  double mbmx = 0.0;                    // Maximum value of distribution
  double gam[ng]; //= zeros(1,ng);
  double mbdist[ng]; //= zeros(1,ng);
  double gmx = 0.0;

  for (int j = 0; j<ng; j++){
    gam[j] = (j-1)*dg+1.0;                                             // gam factor
    mbdist[j] = nfact*gam[j]*sqrt(gam[j]*gam[j]-1.0)*exp(-aa*gam[j]);   // M-B distribution
    if (mbmx < mbdist[j]){
      mbmx = mbdist[j];
      gmx = gam[j];
    } // Maximum of M_B distribution
  }

  // Generate random numbers with M-B distribution
  while (1) { //TODO this is sketchy, while(1) is prolly a really bad idea... 
    double yy = mbmx*(gend-1)*mt.rand(1.0);       // random number between [0,mbmx*(gend-1)]
    double rr = mbmx*mt.rand(1.0);                // random number between [0,mbmx]
    double gg = 1+yy/mbmx;                // compute gam from yy
    double mb = nfact*gg*sqrt(gg*gg-1.0)*exp(-aa*gg);   // evaluate M-B dist at gg 
    if (rr <  mb){                // accept point

      return gg;

      break;
    }                         // reject point
  }

}*/
