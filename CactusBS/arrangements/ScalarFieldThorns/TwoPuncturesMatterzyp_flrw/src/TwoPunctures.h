// TwoPunctures:  File  "TwoPunctures.h"

#define StencilSize 19
#define N_PlaneRelax 1
#define NRELAX 200
#define Step_Relax 1

typedef struct DERIVS
{
  double *d0, *d1, *d2, *d3, *d11, *d12, *d13, *d22, *d23, *d33;
} derivs;

/*
Files of "TwoPunctures":
	TwoPunctures.c
	FuncAndJacobian.c
	CoordTransf.c
	Equations.c
	Newton.c
	utilities.c (see utilities.h)
**************************
*/

// Routines in  "TwoPunctures.c"
double TestSolution (double A, double B, double X, double R, double phi);
void TestVector_w (double *par, int nvar, int n1, int n2, int n3, double *w);

// Routines in  "FuncAndJacobian.c"
int Index (int ivar, int i, int j, int k, int nvar, int n1, int n2, int n3);
void allocate_derivs (derivs * v, int n);
void free_derivs (derivs * v, int n);
void Derivatives_AB3 (int nvar, int n1, int n2, int n3, derivs v);
void F_of_v (int nvar, int n1, int n2, int n3, derivs v,
	     double *F, derivs u, double time);
void J_times_dv (int nvar, int n1, int n2, int n3, derivs dv,
		 double *Jdv, derivs u, double time);
void JFD_times_dv (int i, int j, int k, int nvar, int n1,
		   int n2, int n3, derivs dv, derivs u, double *values, double time);
void SetMatrix_JFD (int nvar, int n1, int n2, int n3,
		    derivs u, int *ncols, int **cols, double **Matrix, double time);
double PunctEvalAtArbitPosition (double *v, int ivar, double A, double B, double phi,
				 int nvar, int n1, int n2, int n3);
double PunctEvalCoeffsAtArbitPosition (double *v, int ivar, double A, double B, double phi,
                                 int nvar, int n1, int n2, int n3);
void calculate_derivs (int i, int j, int k, int ivar, int nvar, int n1,
		       int n2, int n3, derivs v, derivs vv);
double interpol (double a, double b, double c, derivs v);
double PunctTaylorExpandAtArbitPosition (int ivar, int nvar, int n1,
                                         int n2, int n3, derivs v, double x,
                                         double y, double z);
double PunctIntPolAtArbitPosition (int ivar, int nvar, int n1,
				   int n2, int n3, derivs v, double x,
				   double y, double z);
double PunctIntPolCoeffsAtArbitPosition (int ivar, int nvar, int n1,
				   int n2, int n3, derivs v, double x,
				   double y, double z);

// Routines in  "CoordTransf.c"
void AB_To_XR (int nvar, double A, double B, double *X,
	       double *R, derivs U);
void C_To_c (int nvar, double X, double R, double *x,
	     double *r, derivs U);
void rx3_To_xyz (int nvar, double x, double r, double phi, double *y,
		 double *z, derivs U);
void AB3_To_xyz(double A, double B, double phi, double *x, double *y, 
                 double *z);

// Routines in  "Equations.c"
double BY_KKofxyz (double x, double y, double z, double time);
void BY_Aijofxyz (double x, double y, double z, double time, double Aij[3][3]);
void Calc_AijofUP (double x, double y, double z, derivs UV, double Aij[3][3]);
void SetAijAij(int nvar, int n1, int n2, int n3, derivs v, double time, double *AijAij);
void GetAijMatter(int nvar, int n1, int n2, int n3, derivs v, derivs Aij_matter);
void FD_AijofUP(const derivs w, double (*interpolate)(int c, int nvar, int n1, int n2, int n3, derivs v, double x, double y, double z), double x, double y, double z, const double *h, CCTK_REAL Aij[6], int n1, int n2, int n3);
void NonLinEquations_Psi (int i, int j, int k, int n1, int n2, int n3,
                      double A, double B, double X, double R,
		      double x, double r, double phi,
		      double y, double z, derivs U, double *values, double time);
void LinEquations_Psi (int i, int j, int k, int n1, int n2, int n3, 
	           double A, double B, double X, double R,
		   double x, double r, double phi,
		   double y, double z, derivs dU, derivs U, double *values, double time);
void getRho( double x, double y, double z, double r_plus, double r_minus, double *rhohat, double jhat[3] );
void makeInitialGuess_Psi(double x, double y, double z, double *u, double time);
void NonLinEquations_Aij (int i, int j, int k, int n1, int n2, int n3,
                      double A, double B, double X, double R,
		      double x, double r, double phi,
		      double y, double z, derivs U, double *values, double time);
void LinEquations_Aij (int i, int j, int k, int n1, int n2, int n3, 
	           double A, double B, double X, double R,
		   double x, double r, double phi,
		   double y, double z, derivs dU, derivs U, double *values, double time);
void makeInitialGuess_Aij(double x, double y, double z, double *u, double time);

// Routine in "Teukolsky.c"
void addteuk (double x, double y, double z, double t, int w_ind, double Aij[3][3]);
void Origin_Eppley(double x, double y, double z, double acttime, int w_ind, double kktw[3][3]);
void ABC_Eppley(double x, double y, double z, double acttime, int w_ind, double cofs[5]);
void Texp_ABC_Eppley(double x, double y, double z, double acttime, int w_ind, double cofs[5]);
void ExtCurv( double x, double y, double z, int w_ind, double cofs[5], double kktw[3][3]);
void Shape_Eppley( double dff[5], double u, double amp, double sigma, int kij );
void Shape_SMEppley( double dff[5], double u, double amp, double sigma, double k );

// Routines in  "Newton.c"
double
norm_inf (double *F, int ntotal);
void resid (double *res, int ntotal, double *dv, double *rhs,
	    int *ncols, int **cols, double **JFD);
void LineRelax_al (double *dv, int j, int k, int nvar, int n1, int n2, int n3,
		   double *rhs, int *ncols, int **cols, double **JFD);
void LineRelax_be (double *dv, int i, int k, int nvar, int n1, int n2, int n3,
		   double *rhs, int *ncols, int **cols, double **JFD);
void TestRelax (int nvar, int n1, int n2, int n3, derivs v,
		double *dv, double time);
int bicgstab (int nvar, int n1, int n2, int n3, derivs v,
	      derivs dv, int output, int itmax, double tol, double *normres, double time);
int Newton (int nvar, int n1, int n2, int n3, derivs v,
	     double tol, int itmax, double time);
void getInitialGuess(int nvar, int n1, int n2, int n3, derivs v, double time);

// function pointers describing the phyical system to Newton.c
extern void (*NonLinEquations) (int i, int j, int k, int n1, int n2, int n3, 
	              double A, double B, double X, double R,
		      double x, double r, double phi,
		      double y, double z, derivs U, double *values, double time);
extern void (*LinEquations) (int i, int j, int k, int n1, int n2, int n3, 
	           double A, double B, double X, double R,
		   double x, double r, double phi,
		   double y, double z, derivs dU, derivs U, double *values, double time);
extern void (*makeInitialGuess)(double x, double y, double z, double *u, double time);
extern void *NewtonParam; // used to pass parameters to (Non)LinEquations and makeInitialGuess


/* 
 27: -1.325691774825335e-03
 37: -1.325691778944117e-03
 47: -1.325691778942711e-03
 
 17: -1.510625972641537e-03
 21: -1.511443006977708e-03
 27: -1.511440785153687e-03
 37: -1.511440809549005e-03
 39: -1.511440809597588e-03
 */
