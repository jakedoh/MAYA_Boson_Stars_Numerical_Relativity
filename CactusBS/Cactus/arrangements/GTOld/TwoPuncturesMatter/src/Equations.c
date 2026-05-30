// TwoPunctures:  File  "Equations.c"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <assert.h>
#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"
#include "TP_utilities.h"
#include "TwoPunctures.h"

#define DEBUG

// utility macros
#define SQR(x) ((x)*(x))
#define QUAD(x) (SQR(x)*SQR(x))
#define MAX(x, y) ((x) > (y) ? (x) : (y))

// U.d0[ivar]   = U[ivar];  (ivar = 0..nvar-1) 
// U.d1[ivar]   = U[ivar]_x;  
// U.d2[ivar]   = U[ivar]_y;  
// U.d3[ivar]   = U[ivar]_z;  
// U.d11[ivar]  = U[ivar]_xx; 
// U.d12[ivar]  = U[ivar]_xy; 
// U.d13[ivar]  = U[ivar]_xz;
// U.d22[ivar]  = U[ivar]_yy;
// U.d23[ivar]  = U[ivar]_yz;
// U.d33[ivar]  = U[ivar]_zz;

double
BY_KKofxyz (double x, double y, double z, double time)
{
  DECLARE_CCTK_PARAMETERS;
  int i, j;
  double r_plus, r2_plus, r3_plus, r_minus, r2_minus, r3_minus, np_Pp, nm_Pm,
    Aij, AijAij, n_plus[3], n_minus[3], np_Sp[3], nm_Sm[3], TWAij[3][3];

  r_plus = sqrt( (x - par_b) * (x - par_b) + y * y + z * z );
  r_minus = sqrt( (x + par_b) * (x + par_b) + y * y + z * z );

  r_plus = pow(pow(r_plus,4) + pow(epsilon,4),0.25);
  r_minus = pow(pow(r_minus,4) + pow(epsilon,4),0.25);

  r2_plus = r_plus*r_plus;
  r2_minus = r_minus*r_minus;
  r3_plus = r_plus * r2_plus;
  r3_minus = r_minus * r2_minus;

  n_plus[0] = (x - par_b) / r_plus;
  n_minus[0] = (x + par_b) / r_minus;
  n_plus[1] = y / r_plus;
  n_minus[1] = y / r_minus;
  n_plus[2] = z / r_plus;
  n_minus[2] = z / r_minus;

  // dot product: np_Pp = (n_+).(P_+); nm_Pm = (n_-).(P_-) 
  np_Pp = 0;
  nm_Pm = 0;
  for (i = 0; i < 3; i++)
  {
    np_Pp += n_plus[i] * par_P_plus[i];
    nm_Pm += n_minus[i] * par_P_minus[i];
  }
  // cross product: np_Sp[i] = [(n_+) x (S_+)]_i; nm_Sm[i] = [(n_-) x (S_-)]_i
  np_Sp[0] = n_plus[1] * par_S_plus[2] - n_plus[2] * par_S_plus[1];
  np_Sp[1] = n_plus[2] * par_S_plus[0] - n_plus[0] * par_S_plus[2];
  np_Sp[2] = n_plus[0] * par_S_plus[1] - n_plus[1] * par_S_plus[0];
  nm_Sm[0] = n_minus[1] * par_S_minus[2] - n_minus[2] * par_S_minus[1];
  nm_Sm[1] = n_minus[2] * par_S_minus[0] - n_minus[0] * par_S_minus[2];
  nm_Sm[2] = n_minus[0] * par_S_minus[1] - n_minus[1] * par_S_minus[0];
  AijAij = 0;

  for (i = 0; i < 3; i++)
  {
      for (j = 0; j < 3; j++)
      {
          TWAij[i][j] = 0;
      }
  }

  // Get the TW addition to Aij
  for (i=0; i < nwaves; i++)
  {
  if (add_teuk_wave)
     addteuk(x, y, z, time, i, TWAij);
  }

  for (i = 0; i < 3; i++)
  {
    for (j = 0; j < 3; j++)
    {				// Bowen-York-Curvature :
      Aij =
	+ 1.5 * (par_P_plus[i] * n_plus[j] + par_P_plus[j] * n_plus[i]
                 + np_Pp * n_plus[i] * n_plus[j]) / r2_plus
	+ 1.5 * (par_P_minus[i] * n_minus[j] + par_P_minus[j] * n_minus[i]
		 + nm_Pm * n_minus[i] * n_minus[j]) / r2_minus
	- 3.0 * (np_Sp[i] * n_plus[j] + np_Sp[j] * n_plus[i]) / r3_plus
	- 3.0 * (nm_Sm[i] * n_minus[j] + nm_Sm[j] * n_minus[i]) / r3_minus
	+ TWAij[i][j];
      if (i == j)
	Aij -= +1.5 * (np_Pp / r2_plus + nm_Pm / r2_minus);

      AijAij += Aij * Aij;
    }
  }

  return AijAij;
}

void
BY_Aijofxyz (double x, double y, double z, double time, double Aij[3][3])
{
  DECLARE_CCTK_PARAMETERS;
  int i, j;
  double r_plus, r2_plus, r3_plus, r_minus, r2_minus, r3_minus, np_Pp, nm_Pm,
    n_plus[3], n_minus[3], np_Sp[3], nm_Sm[3];

  r_plus = sqrt( (x - par_b) * (x - par_b) + y * y + z * z );
  r_minus = sqrt( (x + par_b) * (x + par_b) + y * y + z * z );

  r_plus = pow(pow(r_plus,4) + pow(epsilon,4),0.25);
  r_minus = pow(pow(r_minus,4) + pow(epsilon,4),0.25);

  r2_plus = r_plus*r_plus;
  r2_minus = r_minus*r_minus;
  r3_plus = r_plus * r2_plus;
  r3_minus = r_minus * r2_minus;

  n_plus[0] = (x - par_b) / r_plus;
  n_minus[0] = (x + par_b) / r_minus;
  n_plus[1] = y / r_plus;
  n_minus[1] = y / r_minus;
  n_plus[2] = z / r_plus;
  n_minus[2] = z / r_minus;

  // dot product: np_Pp = (n_+).(P_+); nm_Pm = (n_-).(P_-) 
  np_Pp = 0;
  nm_Pm = 0;
  for (i = 0; i < 3; i++)
  {
    np_Pp += n_plus[i] * par_P_plus[i];
    nm_Pm += n_minus[i] * par_P_minus[i];
  }
  // cross product: np_Sp[i] = [(n_+) x (S_+)]_i; nm_Sm[i] = [(n_-) x (S_-)]_i
  np_Sp[0] = n_plus[1] * par_S_plus[2] - n_plus[2] * par_S_plus[1];
  np_Sp[1] = n_plus[2] * par_S_plus[0] - n_plus[0] * par_S_plus[2];
  np_Sp[2] = n_plus[0] * par_S_plus[1] - n_plus[1] * par_S_plus[0];
  nm_Sm[0] = n_minus[1] * par_S_minus[2] - n_minus[2] * par_S_minus[1];
  nm_Sm[1] = n_minus[2] * par_S_minus[0] - n_minus[0] * par_S_minus[2];
  nm_Sm[2] = n_minus[0] * par_S_minus[1] - n_minus[1] * par_S_minus[0];
  for (i = 0; i < 3; i++)
  {
    for (j = 0; j < 3; j++)
    {				// Bowen-York-Curvature :
      Aij[i][j] =
        + 1.5 * (par_P_plus[i] * n_plus[j] + par_P_plus[j] * n_plus[i]
		 + np_Pp * n_plus[i] * n_plus[j]) / r2_plus
	+ 1.5 * (par_P_minus[i] * n_minus[j] + par_P_minus[j] * n_minus[i]
		 + nm_Pm * n_minus[i] * n_minus[j]) / r2_minus
	- 3.0 * (np_Sp[i] * n_plus[j] + np_Sp[j] * n_plus[i]) / r3_plus
	- 3.0 * (nm_Sm[i] * n_minus[j] + nm_Sm[j] * n_minus[i]) / r3_minus;
      if (i == j)
	Aij[i][j] -= +1.5 * (np_Pp / r2_plus + nm_Pm / r2_minus);
    }
  }

  // Adding the Teukolsky Wave
  for (i=0; i<nwaves; i++)
  {
      if (add_teuk_wave) {
	 addteuk(x, y, z, time, i, Aij);
      }
  }

}

void Calc_AijofUP(double x, double y, double z, derivs W, double Aij[3][3])
{
  // calculate Aij from Wi which in turn is given by U,Pi (see arXiv:gr-qc/0211028v1)
  double divP;
  #define U 3
  #define P1 0
  #define P2 1
  #define P3 2
  divP = W.d1[P1] + W.d2[P2] + W.d3[P3];
  Aij[0][0] = 1.5*W.d1[P1] - 0.25*W.d11[U] - 0.25*(x*W.d11[P1]+y*W.d11[P2]+z*W.d11[P3]) - 0.5*divP;
  Aij[0][1] = 0.75*(W.d1[P2]+ W.d2[P1]) - 0.25*W.d12[U] - 0.25*(x*W.d12[P1]+y*W.d12[P2]+z*W.d12[P3]);
  Aij[0][2] = 0.75*(W.d1[P3]+ W.d3[P1]) - 0.25*W.d13[U] - 0.25*(x*W.d13[P1]+y*W.d13[P2]+z*W.d13[P3]);

  Aij[1][0] = 0.75*(W.d2[P1]+ W.d1[P2]) - 0.25*W.d12[U] - 0.25*(x*W.d12[P1]+y*W.d12[P2]+z*W.d12[P3]);
  Aij[1][1] = 1.5*W.d2[P2] - 0.25*W.d22[U] - 0.25*(x*W.d22[P1]+y*W.d22[P2]+z*W.d22[P3]) - 0.5*divP;
  Aij[1][2] = 0.75*(W.d2[P3]+ W.d3[P2]) - 0.25*W.d23[U] - 0.25*(x*W.d23[P1]+y*W.d23[P2]+z*W.d23[P3]);

  Aij[2][0] = 0.75*(W.d3[P1]+ W.d1[P3]) - 0.25*W.d13[U] - 0.25*(x*W.d13[P1]+y*W.d13[P2]+z*W.d13[P3]);
  Aij[2][1] = 0.75*(W.d3[P2]+ W.d2[P3]) - 0.25*W.d23[U] - 0.25*(x*W.d23[P1]+y*W.d23[P2]+z*W.d23[P3]);
  Aij[2][2] = 1.5*W.d3[P3] - 0.25*W.d33[U] - 0.25*(x*W.d33[P1]+y*W.d33[P2]+z*W.d33[P3]) - 0.5*divP;
  #undef U
  #undef P1
  #undef P2
  #undef P3
}

void SetAijAij(int nvar, int n1, int n2, int n3, derivs v, double time, double *AijAij)
{
  // calculate the square AijAij and store it in a 1D array for later use
  DECLARE_CCTK_PARAMETERS;

  if(!zero_velocity)
    Derivatives_AB3 (nvar, n1, n2, n3, v);

  // loop over all grid points and calculate BY+Matter Aij, then square and add
#pragma omp parallel for
  for (int i = 0; i < n1; i++)
  {
    double al, be, A, B, X, R, x, r, phi, y, z, Am1;
    derivs U;
    double Aij_matter[3][3], Aij[3][3];

    allocate_derivs (&U, nvar);

    // if zero_velocity then matter does not contribute below
    for(int ii = 0 ; ii < 3 ; ii++)
      for(int jj = 0 ; jj < 3 ; jj++)
        Aij_matter[ii][jj] = 0.;

    for (int j = 0; j < n2; j++)
    {
      for (int k = 0; k < n3; k++)
      {
	// calculate the x,y,z coordinates of the grid point
	al = Pih * (2 * i + 1) / n1;
	A = -cos (al);
	be = Pih * (2 * j + 1) / n2;
	B = -cos (be);
	phi = 2. * Pi * k / n3;

        if (!zero_velocity)
        {
          Am1 = A - 1;
          for (int ivar = 0; ivar < nvar; ivar++)
          {
            int indx = Index (ivar, i, j, k, nvar, n1, n2, n3);
            U.d0[ivar] = Am1 * v.d0[indx];	// U
            U.d1[ivar] = v.d0[indx] + Am1 * v.d1[indx];	// U_A
            U.d2[ivar] = Am1 * v.d2[indx];	// U_B
            U.d3[ivar] = Am1 * v.d3[indx];	// U_3
            U.d11[ivar] = 2 * v.d1[indx] + Am1 * v.d11[indx];	// U_AA
            U.d12[ivar] = v.d2[indx] + Am1 * v.d12[indx];	// U_AB
            U.d13[ivar] = v.d3[indx] + Am1 * v.d13[indx];	// U_AB
            U.d22[ivar] = Am1 * v.d22[indx];	// U_BB
            U.d23[ivar] = Am1 * v.d23[indx];	// U_B3
            U.d33[ivar] = Am1 * v.d33[indx];	// U_33
          }
        }
	// Calculation of (X,R) and
	// (U_X, U_R, U_3, U_XX, U_XR, U_X3, U_RR, U_R3, U_33)
	AB_To_XR (nvar, A, B, &X, &R, U);
	// Calculation of (x,r) and
	// (U, U_x, U_r, U_3, U_xx, U_xr, U_x3, U_rr, U_r3, U_33)
	C_To_c (nvar, X, R, &x, &r, U);
	// Calculation of (y,z) and
	// (U, U_x, U_y, U_z, U_xx, U_xy, U_xz, U_yy, U_yz, U_zz)
	rx3_To_xyz (nvar, x, r, phi, &y, &z, U);
	

	// BY curvature
	BY_Aijofxyz(x,y,z, time, Aij);

	// Matter curvature (if needed)
	if(!zero_velocity)
	  Calc_AijofUP(x,y,z, U, Aij_matter);

	int indx = Index (0, i, j, k, 1, n1, n2, n3);
	AijAij[indx] = 0.;
	for(int ii = 0 ; ii < 3 ; ii++)
	{
	  for(int jj = 0 ; jj < 3 ; jj++)
	  {
	    double comp;	// one component of the total Aij
	    
	    comp = Aij[ii][jj]+Aij_matter[ii][jj];
	    AijAij[indx] += SQR(comp);
	  }
	}
      }
    }
    free_derivs (&U, nvar);
  }
}

void GetAijMatter(int nvar, int n1, int n2, int n3, derivs v, derivs Aij_matter)
{
  /* calculate the spectral decomposition of Aij */
  DECLARE_CCTK_PARAMETERS;
  const int ncomp = 6; /* six components to a symmetric tensor */


  Derivatives_AB3 (nvar, n1, n2, n3, v);

  /* loop over all grid points and calculate Matter Aij */
#pragma omp parallel for
  for (int k = 0; k < n3; k++)
  {
    double al, be, A, B, X, R, x, r, phi, y, z, Am1;
    derivs U;
    double Aij[3][3];

    allocate_derivs (&U, nvar);

    for (int j = 0; j < n2; j++)
    {
      for (int i = 0; i < n1; i++)
      {

	/* calculate the x,y,z coordinates of the grid point */
	al = Pih * (2 * i + 1) / n1;
	A = -cos (al);
	be = Pih * (2 * j + 1) / n2;
	B = -cos (be);
	phi = 2. * Pi * k / n3;

	Am1 = A - 1;
	for (int ivar = 0; ivar < nvar; ivar++)
	{
	  int indx = Index (ivar, i, j, k, nvar, n1, n2, n3);
	  U.d0[ivar] = Am1 * v.d0[indx];	/* U */
	  U.d1[ivar] = v.d0[indx] + Am1 * v.d1[indx];	/* U_A */
	  U.d2[ivar] = Am1 * v.d2[indx];	/* U_B */
	  U.d3[ivar] = Am1 * v.d3[indx];	/* U_3 */
	  U.d11[ivar] = 2 * v.d1[indx] + Am1 * v.d11[indx];	/* U_AA */
	  U.d12[ivar] = v.d2[indx] + Am1 * v.d12[indx];	/* U_AB */
	  U.d13[ivar] = v.d3[indx] + Am1 * v.d13[indx];	/* U_AB */
	  U.d22[ivar] = Am1 * v.d22[indx];	/* U_BB */
	  U.d23[ivar] = Am1 * v.d23[indx];	/* U_B3 */
	  U.d33[ivar] = Am1 * v.d33[indx];	/* U_33 */
	}
	/* Calculation of (X,R) and
	 * (U_X, U_R, U_3, U_XX, U_XR, U_X3, U_RR, U_R3, U_33) */
	AB_To_XR (nvar, A, B, &X, &R, U);
	/* Calculation of (x,r) and
	 * (U, U_x, U_r, U_3, U_xx, U_xr, U_x3, U_rr, U_r3, U_33) */
	C_To_c (nvar, X, R, &x, &r, U);
	/* Calculation of (y,z) and
	 * (U, U_x, U_y, U_z, U_xx, U_xy, U_xz, U_yy, U_yz, U_zz) */
	rx3_To_xyz (nvar, x, r, phi, &y, &z, U);
	
	/* Matter curvature */
        Calc_AijofUP(x,y,z, U, Aij);

        Aij_matter.d0[Index(0, i, j, k, ncomp, n1, n2, n3)] = Aij[0][0]/Am1;
        Aij_matter.d0[Index(1, i, j, k, ncomp, n1, n2, n3)] = Aij[0][1]/Am1;
        Aij_matter.d0[Index(2, i, j, k, ncomp, n1, n2, n3)] = Aij[0][2]/Am1;
        Aij_matter.d0[Index(3, i, j, k, ncomp, n1, n2, n3)] = Aij[1][1]/Am1;
        Aij_matter.d0[Index(4, i, j, k, ncomp, n1, n2, n3)] = Aij[1][2]/Am1;
        Aij_matter.d0[Index(5, i, j, k, ncomp, n1, n2, n3)] = Aij[2][2]/Am1;

      }
    }
    
    free_derivs (&U, nvar);
  }
  
  /* compute spectral derivatives */
  Derivatives_AB3 (ncomp, n1, n2, n3, Aij_matter);
}

void FD_AijofUP(const derivs w, double (*interpolate)(int c, int nvar, int n1, int n2, int n3, derivs v, double x, double y, double z), double x, double y, double z, const double *h, double Aij[6], int n1, int n2, int n3)
{
  double values[4][27];
  struct { /* [P1,P2,P3,U] */
    double d0[4], d1[4], d2[4], d3[4], d11[4], d12[4], d13[4], d22[4], d23[4], d33[4];
  } W;

  /* interpolate onto cube surrounding x,y,z */
  for (int c = 0; c < 4; c++) {
    for (int i = -1 ; i <= 1 ; i++) {
      for (int j = -1 ; j <= 1 ; j++) {
        for (int k = -1 ; k <= 1 ; k++) {
          values[c][9*(1+k)+3*(1+j)+(1+i)] = interpolate (c, 4, n1, n2, n3, w, x+i*h[0], y+j*h[1], z+k*h[2]);
        }
      }
    }
  }

  /* compute FD derivatives */
  for (int c = 0; c < 4; c++) {
    W.d0[c] = values[c][1*9+1*3+1];

    W.d1[c] = (values[c][1*9+1*3+2] - values[c][1*9+1*3+0])/(2*h[0]);
    W.d2[c] = (values[c][1*9+2*3+1] - values[c][1*9+0*3+1])/(2*h[1]);
    W.d3[c] = (values[c][2*9+1*3+1] - values[c][0*9+1*3+1])/(2*h[2]);

    // mixed derivatives from A&S 25.3.26 (p.885)
    W.d11[c] = (values[c][1*9+1*3+2] - 2*values[c][1*9+1*3+1] + values[c][1*9+1*3+0])/(h[0]*h[0]);
    W.d12[c] = (values[c][1*9+2*3+2] - values[c][1*9+2*3+0] - values[c][1*9+0*3+2] + values[c][1*9+0*3+0])/(4*h[0]*h[1]);
    W.d13[c] = (values[c][2*9+1*3+2] - values[c][2*9+1*3+0] - values[c][0*9+1*3+2] + values[c][0*9+1*3+0])/(4*h[0]*h[2]);
    W.d22[c] = (values[c][1*9+2*3+1] - 2*values[c][1*9+1*3+1] + values[c][1*9+0*3+1])/(h[1]*h[1]);
    W.d23[c] = (values[c][2*9+2*3+1] - values[c][2*9+0*3+1] - values[c][0*9+2*3+1] + values[c][0*9+0*3+1])/(4*h[1]*h[2]);
    W.d33[c] = (values[c][2*9+1*3+1] - 2*values[c][1*9+1*3+1] + values[c][0*9+1*3+1])/(h[2]*h[2]);
  }

  /* calculate Aij from Wi which in turn is given by U,Pi (see arXiv:gr-qc/0211028v1) */
  double divP;
  #define U 3
  #define P1 0
  #define P2 1
  #define P3 2
  divP = W.d1[P1] + W.d2[P2] + W.d3[P3];
  Aij[0] = 1.5*W.d1[P1] - 0.25*W.d11[U] - 0.25*(x*W.d11[P1]+y*W.d11[P2]+z*W.d11[P3]) - 0.5*divP;
  Aij[1] = 0.75*(W.d1[P2]+ W.d2[P1]) - 0.25*W.d12[U] - 0.25*(x*W.d12[P1]+y*W.d12[P2]+z*W.d12[P3]);
  Aij[2] = 0.75*(W.d1[P3]+ W.d3[P1]) - 0.25*W.d13[U] - 0.25*(x*W.d13[P1]+y*W.d13[P2]+z*W.d13[P3]);

  Aij[3] = 1.5*W.d2[P2] - 0.25*W.d22[U] - 0.25*(x*W.d22[P1]+y*W.d22[P2]+z*W.d22[P3]) - 0.5*divP;
  Aij[4] = 0.75*(W.d2[P3]+ W.d3[P2]) - 0.25*W.d23[U] - 0.25*(x*W.d23[P1]+y*W.d23[P2]+z*W.d23[P3]);

  Aij[5] = 1.5*W.d3[P3] - 0.25*W.d33[U] - 0.25*(x*W.d33[P1]+y*W.d33[P2]+z*W.d33[P3]) - 0.5*divP;
  #undef U
  #undef P1
  #undef P2
  #undef P3
}

//---------------------------------------------------------------
//*******           Nonlinear Equations                *********/
//---------------------------------------------------------------
void
NonLinEquations_Aij (int i, int j, int k, int n1, int n2, int n3, 
                 double A, double B, double X, double R,
		 double x, double r, double phi,
		 double y, double z, derivs U, double *values, double time)
{
  DECLARE_CCTK_PARAMETERS;
  double r_plus, r_minus;
  double rhohat, jhat[3];

  r_plus = sqrt ((x - par_b) * (x - par_b) + y * y + z * z);
  r_minus = sqrt ((x + par_b) * (x + par_b) + y * y + z * z);

  getRho ( x, y, z, r_plus, r_minus, &rhohat, jhat);

  // fill in RHS
  values[0] = // P1
    U.d11[0] + U.d22[0] + U.d33[0]; 
  values[1] = // P2
    U.d11[1] + U.d22[1] + U.d33[1]; 
  values[2] = // P3
    U.d11[2] + U.d22[2] + U.d33[2]; 
  values[3] = // U
    U.d11[3] + U.d22[3] + U.d33[3];

  if (rhohat > 0) {
    values[0] = values[0] - 8*Pi * jhat[0];
    values[1] = values[1] - 8*Pi * jhat[1];
    values[2] = values[2] - 8*Pi * jhat[2];
    values[3] = values[3] + 8*Pi * (x*jhat[0] + y*jhat[1] + z*jhat[2]);
  } else {
    // TODO: add something useful here (roland)
  }

}

void
NonLinEquations_Psi (int i, int j, int k, int n1, int n2, int n3, 
                 double A, double B, double X, double R,
		 double x, double r, double phi,
		 double y, double z, derivs U, double *values, double time)
{
  DECLARE_CCTK_PARAMETERS;
  double r_plus, r_minus, psi, psi2, psi4, psi5, psi7;

  r_plus = sqrt ((x - par_b) * (x - par_b) + y * y + z * z);
  r_minus = sqrt ((x + par_b) * (x + par_b) + y * y + z * z);

  psi =
    1. + 0.5 * par_m_plus / r_plus + 0.5 * par_m_minus / r_minus + U.d0[0];
  psi2 = psi * psi;
  psi4 = psi2 * psi2;
  psi5 = psi4 * psi;
  psi7 = psi * psi2 * psi4;

  if(!hydro_field) {

    values[0] =
      U.d11[0] + U.d22[0] + U.d33[0] + 0.125 * BY_KKofxyz (x, y, z, time) / psi7;

  } else {

    double AijAij; 		// square of the total Aij
    double rhohat, jhat[3];
    int indx;

    // calculate some shorthands
    indx = Index (0, i, j, k, 1, n1, n2, n3);
    AijAij = ((double *)NewtonParam)[indx];

    getRho ( x, y, z, r_plus, r_minus, &rhohat, jhat);

    // fill in RHS
    values[0] = // u
      U.d11[0] + U.d22[0] + U.d33[0] + 0.125 * AijAij / psi7;

    if (rhohat > 0) {
    	values[0] = values[0] + 2*Pi * rhohat/(psi2*psi);
    } else {
    	values[0] = values[0] + 2*Pi * rhohat*psi4*psi;
    }

  }

  // Brans Dicke field as in Carroll 4.8 //old: Damour&Farese PRL_70_2220, Sotani&Kokkotas PRD_70_084026
  if(scalar_field) {
    double scalarphi, grad_scalarphi[3], scalarpi, grad_scalarphi_sq, VV;
    // caclInitialScalarFields expects the coordinates of the simulation domain
    calcInitialScalarFields(x - move_origin_x, y, z, &scalarphi, &grad_scalarphi[0], &scalarpi, &VV);
    grad_scalarphi_sq = grad_scalarphi[0]*grad_scalarphi[0]+grad_scalarphi[1]*grad_scalarphi[1]+grad_scalarphi[2]*grad_scalarphi[2];
    values[0] += Pi * ( scalarpi*scalarpi / psi7 + grad_scalarphi_sq * psi + 2.0 * VV * psi5);
  }
}

//---------------------------------------------------------------
//*******               Linear Equations                *********/
//---------------------------------------------------------------
void
LinEquations_Aij (int i, int j, int k, int n1, int n2, int n3, 
              double A, double B, double X, double R,
	      double x, double r, double phi,
	      double y, double z, derivs dU, derivs U, double *values, double time)
{
  DECLARE_CCTK_PARAMETERS;

  // fill in RHS
  values[0] = // P1
    dU.d11[0] + dU.d22[0] + dU.d33[0]; 
  values[1] = // P2
    dU.d11[1] + dU.d22[1] + dU.d33[1]; 
  values[2] = // P3
    dU.d11[2] + dU.d22[2] + dU.d33[2]; 
  values[3] = // U
    dU.d11[3] + dU.d22[3] + dU.d33[3];
} 

void
LinEquations_Psi (int i, int j, int k, int n1, int n2, int n3, 
              double A, double B, double X, double R,
	      double x, double r, double phi,
	      double y, double z, derivs dU, derivs U, double *values, double time)
{
  DECLARE_CCTK_PARAMETERS;
  double r_plus, r_minus, psi, psi2, psi4, psi8;

  r_plus = sqrt ((x - par_b) * (x - par_b) + y * y + z * z);
  r_minus = sqrt ((x + par_b) * (x + par_b) + y * y + z * z);

  psi =
    1. + 0.5 * par_m_plus / r_plus + 0.5 * par_m_minus / r_minus + U.d0[0];
  psi2 = psi * psi;
  psi4 = psi2 * psi2;
  psi8 = psi4 * psi4;

  if (!hydro_field) {

    values[0] = dU.d11[0] + dU.d22[0] + dU.d33[0]
      - 0.875 * BY_KKofxyz (x, y, z, time) / psi8 * dU.d0[0];

  } else {

    double AijAij;	// square of the total Aij
    double rhohat, jhat[3];
    int indx;

    // calculate some shorthands
    indx = Index (0, i, j, k, 1, n1, n2, n3);
    AijAij = ((double *)NewtonParam)[indx];

    getRho ( x, y, z, r_plus, r_minus, &rhohat, jhat);

    // fill in RHS
    values[0] = // u
      dU.d11[0] + dU.d22[0] + dU.d33[0] 
	- 0.875 * AijAij / psi8 * dU.d0[0];

    if (rhohat > 0) {
	values[0] = values[0] - 3 * dU.d0[0] *  2*Pi * rhohat/psi4;
    } else {
	values[0] = values[0] + 5 * dU.d0[0] * 2*Pi * rhohat*psi4;
	// TODO: add something useful here (roland)
    }
    
  } 

  // Brans Dicke field as in Damour&Farese PRL_70_2220, Sotani&Kokkotas PRD_70_084026 
  if(scalar_field) {
    double scalarphi, grad_scalarphi[3], scalarpi, grad_scalarphi_sq, VV;
    // caclInitialScalarFields expects the coordinates of the simulation domain
    calcInitialScalarFields(x - move_origin_x, y, z, &scalarphi, &grad_scalarphi[0], &scalarpi, &VV);
    grad_scalarphi_sq = grad_scalarphi[0]*grad_scalarphi[0]+grad_scalarphi[1]*grad_scalarphi[1]+grad_scalarphi[2]*grad_scalarphi[2];
    values[0] += -7.0*Pi * scalarpi*scalarpi / psi8 * dU.d0[0] + Pi * grad_scalarphi_sq * dU.d0[0] + 10.0 * Pi * VV * psi4 * dU.d0[0];
  }

}

//---------------------------------------------------------------
  
void getRho ( double x, double y, double z, double r_plus, double r_minus, double *rhohat, double jhat[3])
{
  DECLARE_CCTK_PARAMETERS;
  double psi0, psi, psi2, psi4, psi5;
  double chi, chi2, sig2;
  double r, r2, rshell, rshell2;
  double gauschi;

  // Important: Often these are rescaled for better convergence.  If it is,
  // we have to be careful that the sum of the psi factors yields +5.  This 
  // is not important if we're just testing the solver by feeding it a solution
  // for u.
  // Note: If take_deriv=0, du contains junk

  psi0 = 1 + par_m_plus/(2*r_plus) + par_m_minus/(2*r_minus);

  if ( test_u ) {

  	  *rhohat=0;
	  for ( int i = 0; i < cv_ngauss; i++) {

		  // Create a gaussian 
		  r2 = (x-cv_x0[i])*(x-cv_x0[i]) + (y-cv_y0[i])*(y-cv_y0[i]) + (z-cv_z0[i])*(z-cv_z0[i]);
		  r = sqrt(r2); 
		  rshell =  r - cv_r0[i];
		  rshell2 = rshell*rshell;
	          sig2 = cv_sigma[i]*cv_sigma[i];
	          chi2 = rshell2/sig2;
	          chi = sqrt(chi2);
	          gauschi = cv_amp[i]*exp(-0.5*chi2); 

	          psi = psi0 + gauschi;
	          psi2 = psi*psi;
	          psi4 = psi2*psi2;
	          psi5 = psi4*psi;

	          if ( cv_amp[i] > 0 ) {
	               *rhohat += gauschi*(1.0 + 2.0*(rshell/r) - chi2)*psi2*psi / (sig2);
	          } else {
	               *rhohat += gauschi*(1.0 + 2.0*(rshell/r) - chi2) / (sig2 * psi5);
	          }
	  }
	  jhat[0] = jhat[1] = jhat[2] = 0.;
	  

  } else {

	double rho, eps, press, velx, vely, velz;
        double bvecx, bvecy, bvecz;
	double w_lorentz2;

        assert(!scalar_field && "Brans-Dicke scalar fields assume the Einstein frame. Have to transform matter variables (Eq. (5) in Sotani&Kokkots)");

        
	// calcMatterRho returns primitive variables
        // translate coordinates to coordinate system used by evolution
	calcMatterRho( x - move_origin_x, y, z, &psi, &rho, &eps, &press, &velx, &vely, &velz, &bvecx, &bvecy, &bvecz );

	psi = combineInitialGuesses(psi0, psi);

	// convert to quantities measured on the slice (conservative variables)
	if(rho != 0.)
	{
	    w_lorentz2 = 1 / (1 - QUAD(psi) * ( SQR(velx)+SQR(vely)+SQR(velz) ));
	    jhat[0] = rho * w_lorentz2 * (1+eps+press/rho) * velx;
	    jhat[1] = rho * w_lorentz2 * (1+eps+press/rho) * vely;
	    jhat[2] = rho * w_lorentz2 * (1+eps+press/rho) * velz;
	    if(w_lorentz2 < 1)
	        CCTK_VWarn(0,__LINE__,__FILE__,CCTK_THORNSTRING,"Encountered w_lorentz^2 < 1 at (%g,%g,%g): %g.", x,y,z, w_lorentz2);

            // Magnetic energy and pressure contribution shorthands.
            double B2 = QUAD(psi)*(SQR(bvecx) + SQR(bvecy) + SQR(bvecz));
            double vB = QUAD(psi)*(velx*bvecx + vely*bvecy + velz*bvecz);

	    *rhohat = rho * (1+eps+press/rho) * w_lorentz2 - press + (1.-0.5/w_lorentz2)*B2 - 0.5*SQR(vB);

	}
	else
	{
	    jhat[0] = jhat[1] = jhat[2] = 0.;
	    *rhohat = 0.;
	}

	if ( *rhohat < 0 ) {

	    *rhohat /= pow(psi0,5.);

	}  else {

	    // for Gaussian and constant rho profiles we achieve
	    // Rescaling by psi0^8, but dividing by psi0^8 for regularity
	    // by setting psi=1. in MatterLibs
	    *rhohat *= pow(psi,8.);
	    jhat[0] *= pow(psi,10.);
	    jhat[1] *= pow(psi,10.);
	    jhat[2] *= pow(psi,10.);

	}
	
  }

}

//---------------------------------------------------------------
//*******               makeInitialGuess                *********/
//---------------------------------------------------------------
void makeInitialGuess_Aij (double x, double y, double z, double *u, double time)
{
  DECLARE_CCTK_PARAMETERS;

  u[0] = u[1] = u[2] = u[3] = 0;
}
void makeInitialGuess_Psi (double x, double y, double z, double *u, double time)
{
  DECLARE_CCTK_PARAMETERS;

  if(hydro_field)
  {
    double psi;
    // translate coordinates to coordinate system used by evolution
    // Note: we could consider using combineInitialGuesses here, but I doubt it
    // actually makes things converge faster
    psi = calcMatterPsi(x - move_origin_x, y, z);
    u[0] = psi - 1.; // only count the non-flat part
  }
  else
  {
    u[0] = 0.;
  }
}

//---------------------------------------------------------------
//******* Initial guess for puncture psi                *********/
//---------------------------------------------------------------
CCTK_REAL calcPuncturePsi0 (CCTK_REAL x, CCTK_REAL y, CCTK_REAL z)
{
  DECLARE_CCTK_PARAMETERS;

  double r_plus, r_minus;
  double x_unmoved;      // x is in the evolution's coordinate system
  double psi0;

  x_unmoved = x + move_origin_x;
  r_plus = sqrt ((x_unmoved - par_b) * (x_unmoved - par_b) + y * y + z * z);
  r_minus = sqrt ((x_unmoved + par_b) * (x_unmoved + par_b) + y * y + z * z);

  psi0 = 1.;
  if(par_m_plus > 0.)
  {
#ifdef INFINITY
    if(r_plus == 0.)
      psi0 = INFINITY;
    else
#endif
      psi0 += par_m_plus/(2*r_plus);
  }
  if(par_m_minus > 0)
  {
#ifdef INFINITY
    if(r_minus == 0.)
      psi0 = INFINITY;
    else
#endif
      psi0 += par_m_minus/(2*r_minus);
  }

  return psi0;
}


//---------------------------------------------------------------
//******* Combines matter and BH guesses for psi        *********/
//---------------------------------------------------------------
CCTK_REAL combineInitialPsi0Guesses (CCTK_REAL psi_BH, CCTK_REAL psi_matter)
{
  DECLARE_CCTK_PARAMETERS;
  double retval;

  if(CCTK_Equals(combination_method, "matter only"))
  {
    // ensure backwards compatability with versions
    // that did not support combination_method yet
    // of course this fails should I ever actually need the psi contribution
    // from matter alone...
    if(psi0_limit > -1.)
    {
      assert(psi_matter < psi0_limit); // see comment above
      retval = psi_BH + psi_matter - 1.; // == "average psi"
    }
    else
    {
      retval = psi_matter;
    }
  } 
  else if(CCTK_Equals(combination_method, "BH only"))
  {
    retval = psi_BH;
  }
  else if(CCTK_Equals(combination_method, "average psi"))
  {
    retval = psi_BH + psi_matter - 1.;
  }
  else if(CCTK_Equals(combination_method, "average"))
  {
    retval = pow(QUAD(psi_BH) + QUAD(psi_matter) - 1., 0.25);
  }
  else if(CCTK_Equals(combination_method, "maximum"))
  {
    retval = MAX(psi_BH, psi_matter);
  }
  else
  {
    CCTK_VWarn(0,__LINE__,__FILE__,CCTK_THORNSTRING, 
      "Invalid combination method: %s", combination_method);
    retval = 0.; /* NOTREACHED */
  }

  if(psi0_limit > -1. && retval > psi0_limit)
    retval = psi0_limit;

  return retval;
}
