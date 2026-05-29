#include "Geodesics.hh"
#include "carpet.hh"

extern "C" void Geodesics_Flag(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if(start_t > cctk_time && cctk_iteration != 0) return; // Save some cpu time if we aren't evolving yet...

  int lsh[1];
  int ierr = CCTK_GrouplshGN (cctkGH, 1, lsh, "Geodesics::geo_pos");
  int lnum_geo = lsh[0];

  CCTK_REAL xmn, xmx, ymn, ymx, zmn, zmx;
  CCTK_REAL rA;

  // TODO: rename active_geo to geo_ref or something like that
  // (it stores each geo's ref level or -1 if outside of the domain)
  // ^^that isn't the case (but it could need to be) 
  for( int i=0; i < lnum_geo; i++ )  {
    if( active_geo[i] != -1 ) {
      active_geo[i]=0;
    }
  }

  ierr = CCTK_CoordRange(cctkGH, &xmn, &xmx, -1, "x", "cart3d");
  ierr = CCTK_CoordRange(cctkGH, &ymn, &ymx, -1, "y", "cart3d");
  ierr = CCTK_CoordRange(cctkGH, &zmn, &zmx, -1, "z", "cart3d");

  if( verbose > 2 ) {
    CCTK_VInfo(CCTK_THORNSTRING,
        "Geodesics_Flag(): %f:%f, %f:%f %f:%f",
        (float)xmn, (float)xmx,
        (float)ymn, (float)ymx,
        (float)zmn, (float)zmx); }

  CCTK_REAL ahradius[num_bh][lnum_geo];
  CCTK_REAL ahcenter_x[num_bh],ahcenter_y[num_bh],ahcenter_z[num_bh];

  ierr = HorizonLocalCoordinateOrigin(1, &ahcenter_x[0], &ahcenter_y[0], &ahcenter_z[0]);
  ierr = HorizonRadiusInDirection(1, lnum_geo, geo_x, geo_y, geo_z, ahradius[0]);

  if((int)num_bh==2) {
    ierr = HorizonLocalCoordinateOrigin(2, &ahcenter_x[1], &ahcenter_y[1], &ahcenter_z[1]);
    ierr = HorizonRadiusInDirection(2, lnum_geo, geo_x, geo_y, geo_z, ahradius[1]);
  }

  
  for( int i=0; i < lnum_geo; i++ )  {
    if( active_geo[i] == 0 ) { 
      if(CCTK_IsThornActive("AHFinderDirect")) {
         rA = sqrt(pow(geo_x[i]-ahcenter_x[0],2)
                 + pow(geo_y[i]-ahcenter_y[0],2)
                 + pow(geo_z[i]-ahcenter_z[0],2));
    
        if((int)num_bh==1) {       

          if( rA < ahradius[0][i] )
          {
            active_geo[i] = -1;
//            in_horizon[i] = 0;

            geo_x[i] = 0.0;
            geo_y[i] = 0.0;
            geo_z[i] = 0.0;
            geo_ux[i] = 0.0;
            geo_uy[i] = 0.0;
            geo_uz[i] = 0.0;
            geo_iut[i] = 0.0;
          }
        }

        if((int)num_bh==2) {
          CCTK_REAL rB = sqrt(pow(geo_x[i]-ahcenter_x[1],2)
                   + pow(geo_y[i]-ahcenter_y[1],2)
                   + pow(geo_z[i]-ahcenter_z[1],2));
          if( rA < ahradius[0][i] || rB < ahradius[1][i] )
          {
            active_geo[i] = -1;
//            in_horizon[i] = 0;

            geo_x[i] = 0.0;
            geo_y[i] = 0.0;
            geo_z[i] = 0.0;
            geo_ux[i] = 0.0;
            geo_uy[i] = 0.0;
            geo_uz[i] = 0.0;
            geo_iut[i] = 0.0;
          }
        }
      } else {

        // If were too close to the origin turn the geodesic off
        if((int)num_bh==1)
        {
          rA = sqrt(pow(geo_x[i],2)+pow(geo_y[i],2)+pow(geo_z[i],2));
          if( rA < bh_rad )
          {
            active_geo[i] = -1;
//            in_horizon[i] = 0;

            geo_x[i] = 0.0;
            geo_y[i] = 0.0;
            geo_z[i] = 0.0;
            geo_ux[i] = 0.0;
            geo_uy[i] = 0.0;
            geo_uz[i] = 0.0;
            geo_iut[i] = 0.0;
          }
        }

        if((int)num_bh==2) // if we have 2 BH, get positions from shifttracker
        {
          CCTK_REAL rB, dxA,dyA,dzA,dxB,dyB,dzB;

          dxA = geo_x[i] - st_x[0];
          dyA = geo_y[i] - st_y[0];
          dzA = geo_z[i] - st_z[0];
          dxB = geo_x[i] - st_x[1];
          dyB = geo_y[i] - st_y[1];
          dzB = geo_z[i] - st_z[1];

          rA = sqrt(pow(dxA,2)+pow(dyA,2)+pow(dzA,2));
          rB = sqrt(pow(dxB,2)+pow(dyB,2)+pow(dzB,2));

          if( rA < bh_rad || rB < bh_rad )
          {
            active_geo[i] = -1;
//            in_horizon[i] = 0;

            geo_x[i] = 0.0;
            geo_y[i] = 0.0;
            geo_z[i] = 0.0;
            geo_ux[i] = 0.0;
            geo_uy[i] = 0.0;
            geo_uz[i] = 0.0;
            geo_iut[i] = 0.0;
          }
        }
      }
      // if we're outside the mesh, turn off geodesic
      // this assumes y_min will be the correct min if the run has sym
      // this is the case for my runs but it doesnt have to be (aka this isnt general enough)
      if( geo_x[i] < ymn || geo_x[i] > xmx
          || geo_y[i] < ymn || geo_y[i] > ymx
          || geo_z[i] < ymn || geo_z[i] > zmx 
          || // I don't always get nan's, but when I do, I turn that geodesic off and have a dos equis
             isnan(geo_x[i]) 
          || isnan(geo_y[i]) 
          || isnan(geo_z[i])
        )
      { 
        active_geo[i] = -1;//rlev[i];
//        in_horizon[i] = -1;

        geo_x[i] = 0.0;
        geo_y[i] = 0.0;
        geo_z[i] = 0.0;
        geo_ux[i] = 0.0;
        geo_uy[i] = 0.0;
        geo_uz[i] = 0.0;
        geo_iut[i] = 0.0;
      }

    }

    // Flood stdout (only for debugging)
    if( verbose > 2 ) {
      CCTK_VInfo(CCTK_THORNSTRING,
          "Geodesics_Flag(): Particle %d @ active_geo=%d",
          i,
          (int)active_geo[i]); }

  }
}
