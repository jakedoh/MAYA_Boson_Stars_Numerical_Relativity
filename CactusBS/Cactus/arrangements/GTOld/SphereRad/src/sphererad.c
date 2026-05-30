#include <stdio.h> 

#include "cctk.h" 
#include "cctk_Arguments.h" 
#include "cctk_Parameters.h" 

#define MAXSURFACES 10

int set_spherical_surfaces(CCTK_ARGUMENTS,
			   CCTK_REAL xah[MAXSURFACES],
			   CCTK_REAL yah[MAXSURFACES],
			   CCTK_REAL zah[MAXSURFACES],
			   CCTK_REAL rah[MAXSURFACES])
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  //--------------------------------------------

  int sind=-1;
  for (int i=0;i<max_nsurface;i++)  
  { 
    sind=surface_index[i];
    if (sind<0) {
      CCTK_WARN(1,"surface not selected");
      return -1;
    }
    if (sind>sphericalsurfaces_nsurfaces) {
      CCTK_WARN(1,"increase SphericalSurfaces::nsurfaces parameter!");
      return -1;
    }
    sf_valid[sind]=1;
    sf_active[sind]=1;
    // populate center first
    sf_centroid_x[sind]=xah[i];
    sf_centroid_y[sind]=yah[i];
    sf_centroid_z[sind]=zah[i];
    if (verbose>1) {
      CCTK_VInfo(CCTK_THORNSTRING,"surface %d %d: (%g,%g,%g)",
		 i,sind,xah[i],yah[i],zah[i]);
    }

    if (verbose>1) {
      CCTK_VInfo(CCTK_THORNSTRING,"  r=%g",rah[i]);
    }

    // all points have constant radius
    CCTK_INT ntheta=sf_ntheta[sind];
    CCTK_INT nphi=sf_nphi[sind];
    CCTK_INT ind=0;
    for (int a=0;a<ntheta;a++) {
      for (int b=0;b<nphi;b++) {
	ind=a+maxntheta*(b+maxnphi*sind);
	sf_radius[ind]=rah[i];
      }
    }

    sf_mean_radius[sind] = sf_max_radius[sind] = sf_min_radius[sind] = rah[i];
  }
  return 1;
}

int get_positions(CCTK_ARGUMENTS,
		  CCTK_REAL xah[MAXSURFACES],
		  CCTK_REAL yah[MAXSURFACES],
		  CCTK_REAL zah[MAXSURFACES])
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_INT sn;

  for (int i=0;i<max_nsurface;i++) {
    if (force_position[i]>0) {
      xah[i]=x_fixed[i];
      yah[i]=y_fixed[i];
      zah[i]=z_fixed[i];
    }
    else if (spherical_surface[i]>0) {
      sn=shifttracker_index[i];
      if (sn>=sphericalsurfaces_nsurfaces) {
        CCTK_VInfo(CCTK_THORNSTRING,
                   "surface number sn=%d too large, increase SphericalSurface::nsurfaces from its current value %d",
                   sn,sphericalsurfaces_nsurfaces);
        xah[i]=0;
        yah[i]=0;
        zah[i]=0;
        continue;
      }
      if (sf_valid[sn]<=0) {
        CCTK_VInfo(CCTK_THORNSTRING,
                   "didn't find valid horizon surface for sn=%d, index=%d",sn,i);
        continue;
      }
      xah[i]=sf_centroid_x[sn];
      yah[i]=sf_centroid_y[sn];
      zah[i]=sf_centroid_z[sn];
      if (verbose>1) {
        CCTK_VInfo(CCTK_THORNSTRING,"setting center from sphericalsurface %d to (%.19g,%.19g,%.19g)",
                   sn,xah[i],yah[i],zah[i]);
      }
    }
    else {
      CCTK_INT ind=shifttracker_index[i];
      if (ind<0) {
	if (verbose>1) {
	  CCTK_VInfo(CCTK_THORNSTRING,
		     "warning: ind=%d<0!",ind);
	}
	continue;
      }
      xah[i]=st_x[ind];
      yah[i]=st_y[ind];
      zah[i]=st_z[ind];
    }
  }
  return 1;
}  

int update_center_of_mass(CCTK_ARGUMENTS,
			  CCTK_REAL xah[MAXSURFACES],
			  CCTK_REAL yah[MAXSURFACES],
			  CCTK_REAL zah[MAXSURFACES])
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  const CCTK_REAL total_mass=mass[0]+mass[1];
  const CCTK_REAL xcm=(xah[0]*mass[0]+xah[1]*mass[1]) /total_mass;
  const CCTK_REAL ycm=(yah[0]*mass[0]+yah[1]*mass[1]) /total_mass;
  const CCTK_REAL zcm=(zah[0]*mass[0]+zah[1]*mass[1]) /total_mass;

  if (verbose>1) {
    CCTK_VInfo(CCTK_THORNSTRING,
	       "updating center of mass to (%g,%g,%g)",
	       xcm,ycm,zcm);
  }

  for (int i=0;i<max_nsurface;i++) {
    if (center_of_mass[i]) {
      if (verbose>1) {
	CCTK_VInfo(CCTK_THORNSTRING,"  for i=%d",i);
      }
      xah[i]=xcm;
      yah[i]=ycm;
      zah[i]=zcm;
    }
  }
  return 1;
}

void update_sphere_const_rad (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  CCTK_INT ierr;

  if (verbose>0) {
    CCTK_VInfo(CCTK_THORNSTRING,
	       "setting spherical surfaces: T=%f",cctk_time);
  }

  CCTK_REAL xah[MAXSURFACES], yah[MAXSURFACES],
            zah[MAXSURFACES], rah[MAXSURFACES];

  for (int i=0;i<MAXSURFACES;i++) {
    xah[i]=0;yah[i]=0;zah[i]=0;rah[i]=0;
  }

  // set xah, yah, zah the central positions
  ierr=get_positions(CCTK_PASS_CTOC, xah,yah,zah);

  for (int i=0;i<max_nsurface;i++) {
    rah[i]=radius_surface[i];
  }

  ierr=update_center_of_mass(CCTK_PASS_CTOC, xah,yah,zah);

  ierr=set_spherical_surfaces(CCTK_PASS_CTOC, xah,yah,zah,rah);

  if (verbose>2) {
    CCTK_INFO("done in sphererad");
  }
}
