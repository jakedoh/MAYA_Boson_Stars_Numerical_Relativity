/* Some basic utilities */
#include "RecoverMHD.h"
#include <cmath>

double spatialdeterminant( double gxx, double gxy, double gxz, double gyy, double gyz, double gzz )
{
  double det = 	- SQR(gxz)*gyy - gxx*SQR(gyz) - SQR(gxy)*gzz  
	  	+ 2*gxy*gxz*gyz + gxx*gyy*gzz;
  return det;
}

void invertspatialdet( double gxx, double gxy, double gxz, double gyy, double gyz, double gzz, 
	double det, double& ugxx, double& ugxy, double& ugxz, double& ugyy, double& ugyz, double& ugzz )
{
  ugxx = (-SQR(gyz) + gyy*gzz)/det;
  ugxy = (gxz*gyz - gxy*gzz)/det;
  ugxz = (-gxz*gyy + gxy*gyz)/det;
  ugyy = (-SQR(gxz) + gxx*gzz)/det;
  ugyz = (gxy*gxz - gxx*gyz)/det;
  ugzz = (-SQR(gxy) + gxx*gyy)/det;
}

bool check_interval( double x_lower, double x_upper, double epsrel, double epsabs)
{ 

  double abs_lower, abs_upper;
  double min_abs, tolerance;
  bool interval_symmetric;

  abs_upper = fabs(x_upper);
  abs_lower = fabs(x_lower);
  if ((x_lower > 0.0 && x_upper > 0.0) || (x_lower < 0.0 && x_upper < 0.0)) {
    min_abs = (abs_lower < abs_upper) ? abs_lower : abs_upper;
  } else {
    min_abs = 0;
  }

  tolerance = epsabs + epsrel * min_abs;

  if ( fabs(x_upper - x_lower) <= tolerance ) {
    interval_symmetric = true;
  } else {
    interval_symmetric = false;
  }

  return interval_symmetric;

}
