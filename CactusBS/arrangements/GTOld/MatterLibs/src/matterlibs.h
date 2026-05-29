// MatterLibs: matterlibs.h

#include "cctk.h"
#include "cctk_Arguments.h"

#ifdef __cplusplus
extern "C" {
#endif

// Routines in make_rho.c
void matterlibs_calcMatterRho( double x, double y, double z, double *psi, double *rho, double *eps, double *P, double *vel_x, double *vel_y, double *vel_z, double *bvec_x, double *bvec_y, double *bvec_z ); 
void matterlibs_calcConformalRho( double x, double y, double z, double *psi, double *rho, double *Ebar, double *Jbarx, double *Jbary, double *Jbarz, double *bvec_x, double *bvec_y, double *bvec_z, double *vB, double *w_lorentz );

void calcGaussRho( double x, double y, double z, double *rho, double *vel_x, double *vel_y, double *vel_z );
void calcCBDiskRho( double x, double y, double z, double *rho, double *press, double *vel_x, double *vel_y, double *vel_z );
double calcCloudRho( double x, double y, double z ); 

// Routines in make_psi.c
double matterlibs_calcMatterPsi( double x, double y, double z );

// TOV routines (tov.cc)
void calcTOVRho( double x, double y, double z, double *psi, double *rho, double *eps, double *P, double *v_x, double *v_y, double *v_z ); 
void calcConformalTOVRho( double x, double y, double z, double *dpsi, double *drho0, double *dEbar, double *dJbarx, double *dJbary, double *dJbarz, double *dlorentz ) ;
double calcTOVLapse( double x, double y, double z ); 
extern void matterlibs_preSolveTOVs( CCTK_ARGUMENTS );
void matterlibs_calcTOVQCNRhoPsi(int star, double distance, double *qq, double *cc, double *nn, double *rho, double *psi_nonflat);
void matterlibs_updateTOVs();
double matterlibs_getTOVMass(int star);
double matterlibs_getTOVRadius(int star);

// Bondi Routines (bondi.c)
void calcBondiRho( double x, double y, double z, double *psi, double *rho, double *eps, double *P, double *vel_x, 
                   double *vel_y, double *vel_z, double *bx, double *by, double *bz ); 

// this one is actually kind of private (in E_to_rho.c)
int matterlibs_E_to_rho(int eos_handle, double dE, double dJx, double dJy,
    double dJz, double psi, double *drho, double *dlorentz);
int matterlibs_E_to_rho_MHD(int eos_handle, double dE, double dJx, double dJy, 
    double dJz, double dBx, double dBy, double dBz, double psi, double *drho, double *dlorentz);

extern CCTK_REAL matterlibs_poison;

#ifdef __cplusplus
}
#endif
