/* Includes */

#include "cctk.h"

/* Definitions */
#define SQR(x) ((x)*(x))
#define SPATIAL_DET(gxx_,gxy_,gxz_,gyy_,gyz_,gzz_) \
                   (-(gxz_)**2*(gyy_) + 2*(gxy_)*(gxz_)*(gyz_) - (gxx_)*(gyz_)**2 - (gxy_)**2*(gzz_) \
                   + (gxx_)*(gyy_)*(gzz_))

//#ifndef
//#define HAVE_INLINE
//#endif

/************************/
/* Types, Ennumerations */
/************************/
enum GSL_SOLVER_TYPE { GST_Newton, GST_gNewton, GST_HybridSJ, GST_HybridJ, GST_HybridS, GST_Hybrid, GST_dNewton, GST_Broyden };

/*************************/
/** Function Prototypes **/
/*************************/

/* Pointwise.c */
extern "C" 
CCTK_INT RecoverMHD_Con2Prim_pointwise( CCTK_INT handle, CCTK_INT phandle, CCTK_INT polytype, CCTK_INT iteration, CCTK_INT stype,
	CCTK_INT sys_eqns, CCTK_REAL *dens, CCTK_REAL *tau, CCTK_REAL *sx, CCTK_REAL *sy, CCTK_REAL *sz, CCTK_REAL *Bnx, CCTK_REAL *Bny, CCTK_REAL *Bnz,
	CCTK_REAL *rho, CCTK_REAL *eps, CCTK_REAL *velx, CCTK_REAL *vely, CCTK_REAL *velz, CCTK_REAL *Bvecx, CCTK_REAL *Bvecy, CCTK_REAL *Bvecz,
	CCTK_REAL *w_lorentz, CCTK_REAL *press, CCTK_REAL x, CCTK_REAL y, CCTK_REAL z,
	CCTK_REAL gxx, CCTK_REAL gxy, CCTK_REAL gxz, CCTK_REAL gyy, CCTK_REAL gyz, CCTK_REAL gzz,
	CCTK_INT whisky_reflevel, CCTK_REAL carpetweight, 
	CCTK_INT do_impose_s2_limit, CCTK_REAL min_rho, CCTK_REAL min_tau );

/* RecoverMHD_Utils.c */
double spatialdeterminant( double gxx, double gxy, double gxz, double gyy, double gyz, double gzz );
void invertspatialdet( double gxx, double gxy, double gxz, double gyy, double gyz, double zz, 
	double det, double& ugxx, double& ugxy, double& ugxz, double& ugyy, double& ugyz, double& ugzz );
bool check_interval( double x_lower, double x_upper, double epsrel, double epsabs);

/* Solve_SysEqns.c */
int EOS_solvesys (int eos_handle, enum GSL_SOLVER_TYPE stype, bool is_polytype,
    double udens, double utau, double usx, double usy, double usz,
    double bnx, double bny, double bnz,
    double gxx, double gxy, double gxz, double gyy, double gyz, double gzz,
    double guxx, double guxy, double guxz, double guyy, double guyz, double guzz,
    double *ux, double *uy, double *uz, double *eps, double soln_quality[3]);

int EOS_solveWPresssys (int eos_handle, enum GSL_SOLVER_TYPE stype,
    bool is_polytype,
    double udens, double utau, double usx, double usy, double usz,
    double bnx, double bny, double bnz,
    double gxx, double gxy, double gxz, double gyy, double gyz, double gzz,
    double guxx, double guxy, double guxz, double guyy, double guyz, double guzz,
    double *ux, double *uy, double *uz, double *eps, double soln_quality[3]);

int EOS_solvehW2sys (int eos_handle, enum GSL_SOLVER_TYPE stype,
    bool is_polytype,
    double udens, double utau, double usx, double usy, double usz,
    double bnx, double bny, double bnz,
    double gxx, double gxy, double gxz, double gyy, double gyz, double gzz,
    double guxx, double guxy, double guxz, double guyy, double guyz, double guzz,
    double *ux, double *uy, double *uz, double *eps, double soln_quality[3]);

int EOS_solvehv2sys (int eos_handle, enum GSL_SOLVER_TYPE stype,
    bool is_polytype,
    double udens, double utau, double usx, double usy, double usz,
    double bnx, double bny, double bnz,
    double gxx, double gxy, double gxz, double gyy, double gyz, double gzz,
    double guxx, double guxy, double guxz, double guyy, double guyz, double guzz,
    double *ux, double *uy, double *uz, double *eps, double soln_quality[3]);

int EOS_solveHDsys (int eos_handle,
    bool is_polytype, double rho_min, double pmin,
    double udens, double utau, double usx, double usy, double usz,
    double gxx, double gxy, double gxz, double gyy, double gyz, double gzz,
    double guxx, double guxy, double guxz, double guyy, double guyz, double guzz,
    double *rho, double *ux, double *uy, double *uz, double *eps, double soln_quality[3]);

CCTK_INT EOS_solvev2Z (CCTK_INT eos_handle, CCTK_REAL rho_min, CCTK_REAL atmo_tol,
    CCTK_REAL udens, CCTK_REAL utau, CCTK_REAL usx, CCTK_REAL usy, CCTK_REAL usz,
    CCTK_REAL bnx, CCTK_REAL bny, CCTK_REAL bnz,
    CCTK_REAL gxx, CCTK_REAL gxy, CCTK_REAL gxz, CCTK_REAL gyy, CCTK_REAL gyz, CCTK_REAL gzz,
    CCTK_REAL guxx, CCTK_REAL guxy, CCTK_REAL guxz, CCTK_REAL guyy, CCTK_REAL guyz, CCTK_REAL guzz,
    CCTK_REAL *rho, CCTK_REAL *eps, CCTK_REAL *ux, CCTK_REAL *uy, CCTK_REAL *uz, CCTK_REAL soln_quality[3]);
