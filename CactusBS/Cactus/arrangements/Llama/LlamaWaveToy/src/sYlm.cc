#include <stdlib.h>
#include <cassert>
#include "cctk.h"



/*
   Calculates the real and imaginary spin-weighted spherical harmonics 
   with given spin s for l,m in spherical coordinates.
   For s=0, we obtain the standard spherical harmonics Ylm.
*/
extern "C" void CCTK_FNAME(sylm_re)(CCTK_INT *sss, CCTK_INT *lll, CCTK_INT *mmm, CCTK_REAL *thetatheta, CCTK_REAL *phiphi, CCTK_REAL *ergebnils);
extern "C" void CCTK_FNAME(sylm_im)(CCTK_INT *sss, CCTK_INT *lll, CCTK_INT *mmm, CCTK_REAL *thetatheta, CCTK_REAL *phiphi, CCTK_REAL *ergebnils);



// pi
const double pi = 4.0*atan(1.0);


inline int fac(int n)
{
   int result = 1;
   
   for (int i=2; i <= n; i++)
      result *= i;
	    
   return result;
}


// coefficient function
inline double Cslm(int s, int l, int m)
{
   assert(l >= 0);

   return sqrt( l*l * (4.0*l*l - 1.0) / ( (l*l - m*m) * (l*l - s*s) ) );
}



// recursion function
double s_lambda_lm(int s, int l, int m, double x)
{
   double Pm;

   assert(m >= 0);
   assert(l >= 0);
   assert(m >= s);

   Pm = pow(-0.5, m);
/*
   if (l == m) // recursion end
      return pow(-2.0, -m) * sqrt( fac(2*m + 1) * 1.0 / ( 4.0*pi * fac(m+s) * fac(m-s) ) * pow(1.0-x, m+s) * pow(1.0+x, m-s) );
   else if (l >= m+2)
      return  (x + s*m * 1.0 / ( l * (l-1.0) ) ) * Cslm(s, l, m) * s_lambda_lm(s, l-1, m, x) 
             - Cslm(s, l, m) * 1.0 / Cslm(s, l-1, m) * s_lambda_lm(s, l-2, m, x);
   else if (l >= m+1)
      return  (x + s*m * 1.0 / ( l * (l-1.0) ) ) * Cslm(s, l, m) * s_lambda_lm(s, l-1, m, x);
   else
      return 0.0;
*/

   if (m !=  s) Pm = Pm * pow(1.0+x, (m-s)*1.0/2);
   if (m != -s) Pm = Pm * pow(1.0-x, (m+s)*1.0/2);
   
   Pm = Pm * sqrt( fac(2*m + 1) * 1.0 / ( 4.0*pi * fac(m+s) * fac(m-s) ) );
   
   if (l == m)
      return Pm;
   
   double Pm1 = (x + s*1.0/(m+1) ) * Cslm(s, m+1, m) * Pm;
   
   if (l == m+1)
      return Pm1;
   else
   {
      double Pn;
      for (int n=m+2; n <= l; n++)
      {
         Pn = (x + s*m * 1.0 / ( n * (n-1.0) ) ) * Cslm(s, n, m) * Pm1 
             - Cslm(s, n, m) * 1.0 / Cslm(s, n-1, m) * Pm;
         Pm = Pm1;
         Pm1 = Pn;
         
      }
      return Pn;
   }
}



void CCTK_FNAME(sylm_re)(CCTK_INT *sss, CCTK_INT *lll, CCTK_INT *mmm, CCTK_REAL *thetatheta, CCTK_REAL *phiphi, CCTK_REAL *ergebnils)
{
   CCTK_REAL result;
   CCTK_INT s, l, m;
   CCTK_REAL Pm = 1.0;

   int ss=*sss;
   int ll=*lll;
   int mm=*mmm;
   double theta=*thetatheta;
   double phi=*phiphi;

   l = ll;
   m = mm;
   s = ss;

   if (l < 0)
      return;
   if (abs(m) > l || l < abs(s))
      return;

   if (abs(mm) < abs(ss))
   {
      s=mm;
      m=ss;
      if ((m+s) % 2)
         Pm  = -Pm;
   }
   
   if (m < 0)
   {
      s=-s;
      m=-m;
      if ((m+s) % 2)
         Pm  = -Pm;
   }

   result = Pm * s_lambda_lm(s, l, m, cos(theta));
   *ergebnils = result * cos(mm*phi);
}


void CCTK_FNAME(sylm_im)(CCTK_INT *sss, CCTK_INT *lll, CCTK_INT *mmm, CCTK_REAL *thetatheta, CCTK_REAL *phiphi, CCTK_REAL *ergebnils)
{
   CCTK_REAL result;
   CCTK_INT s, l, m;
   CCTK_REAL Pm = 1.0;

   int ss=*sss;
   int ll=*lll;
   int mm=*mmm;
   double theta=*thetatheta;
   double phi=*phiphi;


   l = ll;
   m = mm;
   s = ss;

   if (l < 0)
      return;
   if (abs(m) > l || l < abs(s))
      return;

   if (abs(mm) < abs(ss))
   {
      s=mm;
      m=ss;
      if ((m+s) % 2)
         Pm  = -Pm;
   }
   
   if (m < 0)
   {
      s=-s;
      m=-m;
      if ((m+s) % 2)
         Pm  = -Pm;
   }

   result = Pm * s_lambda_lm(s, l, m, cos(theta));
   *ergebnils = result * sin(mm*phi);
}


