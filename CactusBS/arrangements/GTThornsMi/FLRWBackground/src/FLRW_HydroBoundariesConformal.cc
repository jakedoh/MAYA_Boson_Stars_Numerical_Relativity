#include <algorithm>
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "vectors.h"
#include "FLRWBackground.h"



#define IDX(i,j,k) CCTK_GFINDEX3D(cctkGH,(i),(j),(k))
#define IDXX(i,j,k) CCTK_VECTGFINDEX3D(cctkGH, (i),(j),(k) , 0)
#define IDXY(i,j,k) CCTK_VECTGFINDEX3D(cctkGH, (i),(j),(k) , 1)
#define IDXZ(i,j,k) CCTK_VECTGFINDEX3D(cctkGH, (i),(j),(k) , 2)

#define XMAX_OB(FUNC,imax,val) for(int k=0;k<cctk_lsh[2];k++) for(int j=0;j<cctk_lsh[1];j++) FUNC[IDX(imax,j,k)] = val;
#define YMAX_OB(FUNC,jmax,val) for(int k=0;k<cctk_lsh[2];k++) for(int i=0;i<cctk_lsh[0];i++) FUNC[IDX(i,jmax,k)] = val;
#define ZMAX_OB(FUNC,kmax,val) for(int j=0;j<cctk_lsh[1];j++) for(int i=0;i<cctk_lsh[0];i++) FUNC[IDX(i,j,kmax)] = val;

#define XMIN_OB(FUNC,imin,val) for(int k=0;k<cctk_lsh[2];k++) for(int j=0;j<cctk_lsh[1];j++) FUNC[IDX(imin,j,k)] = val;
#define YMIN_OB(FUNC,jmin,val) for(int k=0;k<cctk_lsh[2];k++) for(int i=0;i<cctk_lsh[0];i++) FUNC[IDX(i,jmin,k)] = val;
#define ZMIN_OB(FUNC,kmin,val) for(int j=0;j<cctk_lsh[1];j++) for(int i=0;i<cctk_lsh[0];i++) FUNC[IDX(i,j,kmin)] = val;


#define XMAX_OB_VECT(FUNC,imax,val) for(int k=0;k<cctk_lsh[2];k++) for(int j=0;j<cctk_lsh[1];j++) { FUNC[IDXX(imax,j,k)] = val; FUNC[IDXY(imax,j,k)] = val; FUNC[IDXZ(imax,j,k)] = val;}
#define YMAX_OB_VECT(FUNC,jmax,val) for(int k=0;k<cctk_lsh[2];k++) for(int i=0;i<cctk_lsh[0];i++) {FUNC[IDXX(i,jmax,k)] = val;FUNC[IDXY(i,jmax,k)] = val;FUNC[IDXZ(i,jmax,k)] = val;}
#define ZMAX_OB_VECT(FUNC,kmax,val) for(int j=0;j<cctk_lsh[1];j++) for(int i=0;i<cctk_lsh[0];i++) {FUNC[IDXX(i,j,kmax)] = val;FUNC[IDXY(i,j,kmax)] = val;FUNC[IDXZ(i,j,kmax)] = val;}

#define XMIN_OB_VECT(FUNC,imin,val) for(int k=0;k<cctk_lsh[2];k++) for(int j=0;j<cctk_lsh[1];j++) {FUNC[IDXX(imin,j,k)] = val;FUNC[IDXY(imin,j,k)] = val;FUNC[IDXZ(imin,j,k)] = val;}
#define YMIN_OB_VECT(FUNC,jmin,val) for(int k=0;k<cctk_lsh[2];k++) for(int i=0;i<cctk_lsh[0];i++) {FUNC[IDXX(i,jmin,k)] = val;FUNC[IDXY(i,jmin,k)] = val;FUNC[IDXZ(i,jmin,k)] = val;}
#define ZMIN_OB_VECT(FUNC,kmin,val) for(int j=0;j<cctk_lsh[1];j++) for(int i=0;i<cctk_lsh[0];i++) {FUNC[IDXX(i,j,kmin)] = val;FUNC[IDXY(i,j,kmin)] = val;FUNC[IDXZ(i,j,kmin)] = val;}

extern "C" void FLRWBackground_HydroBoundaryConformal(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  bool NotSymx=false; if(CCTK_EQUALS(Symmetry_x,"no")) NotSymx=true;
  bool NotSymy=false; if(CCTK_EQUALS(Symmetry_y,"no")) NotSymy=true;
  bool NotSymz=false; if(CCTK_EQUALS(Symmetry_z,"no")) NotSymz=true;

  int levelnumber = GetRefinementLevel(cctkGH);
  if(cctk_iteration==0 || levelnumber!=0) return;

  const CCTK_REAL t = cctk_time;
  double expansion = a0*t/t0;
  double Hubble = H0*pow(expansion , -2.0);

  double rhocrit = 3.0*H0*H0/(8.0*PI);

  double epst = getEps( rhocrit , K_intrinsic , expansion , tol);
  double rhot = rhocrit*pow(expansion , -4)/(1+epst);
  double presst = rhot*epst/3.0;
  double velt = 0.;
  double w_lorentzt = 1.0;

  for(int which_bdry_pt=0;which_bdry_pt<cctk_nghostzones[0];which_bdry_pt++) {
    int imax=cctk_lsh[0]-cctk_nghostzones[0]+which_bdry_pt; // for cctk_nghostzones==3, this goes {cctk_lsh-3,cctk_lsh-2,cctk_lsh-1}; outer bdry pt is at cctk_lsh-1
    int jmax=cctk_lsh[1]-cctk_nghostzones[1]+which_bdry_pt;
    int kmax=cctk_lsh[2]-cctk_nghostzones[2]+which_bdry_pt;

    int imin=cctk_nghostzones[0]-which_bdry_pt-1; // for cctk_nghostzones==3, this goes {2,1,0}
    int jmin=cctk_nghostzones[1]-which_bdry_pt-1;
    int kmin=cctk_nghostzones[2]-which_bdry_pt-1;

    if(cctk_bbox[1]) { XMAX_OB(eps,imax,epst); XMAX_OB(rho,imax,rhot); XMAX_OB(press,imax,presst); XMAX_OB_VECT(vel,imax,velt); XMAX_OB(w_lorentz,imax,w_lorentzt);}
    if(cctk_bbox[3]) { YMAX_OB(eps,jmax,epst); YMAX_OB(rho,jmax,rhot); YMAX_OB(press,jmax,presst); YMAX_OB_VECT(vel,jmax,velt); YMAX_OB(w_lorentz,jmax,w_lorentzt);}
    if(cctk_bbox[5]) { ZMAX_OB(eps,kmax,epst); ZMAX_OB(rho,kmax,rhot); ZMAX_OB(press,kmax,presst); ZMAX_OB_VECT(vel,kmax,velt); ZMAX_OB(w_lorentz,kmax,w_lorentzt);}

    if((cctk_bbox[0]) && NotSymx) { XMIN_OB(eps,imin,epst); XMIN_OB(rho,imin,rhot); XMIN_OB(press,imin,presst); XMIN_OB_VECT(vel,imin,velt); XMIN_OB(w_lorentz,imin,w_lorentzt);}
    if((cctk_bbox[2]) && NotSymy) { YMIN_OB(eps,jmin,epst); YMIN_OB(rho,jmin,rhot); YMIN_OB(press,jmin,presst); YMIN_OB_VECT(vel,jmin,velt); YMIN_OB(w_lorentz,jmin,w_lorentzt);}
    if((cctk_bbox[4]) && NotSymz) { ZMIN_OB(eps,kmin,epst); ZMIN_OB(rho,kmin,rhot); ZMIN_OB(press,kmin,presst); ZMIN_OB_VECT(vel,kmin,velt); ZMIN_OB(w_lorentz,kmin,w_lorentzt);}
  }
}
