 /*@@
   @file      Startup.c
   @date      Wed Apr 19 20:35:04 2000
   @author    Gerd Lanfermann
   @desc 
   Startup.c
   @enddesc 
 @@*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cctk.h"

#include "Ell_DBstructure.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(CactusElliptic_EllBase_Startup_c)

int Ell_RegisterBaseEqTypes(void);


 /*@@
   @routine    Ell_RegisterBaseEqTypes
   @date       Wed Apr 19 20:36:08 2000
   @author     Gerd Lanfermann
   @desc 
       At Startup, EllBase registers the elliptic equation classes for which 
       it provides solvers. Other routines, which may come up with new classes, 
       can registers the classes in their own thorns.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/


int Ell_RegisterBaseEqTypes(void) 
{

  int Ell_RegisterEq(void (*function),const char *);
  void Ell_LinConfMetricRegistry(void (*function),const char *);
  void Ell_LinMetricRegistry(void (*function),const char *);
  void Ell_LinFlatRegistry(void (*function),const char *);
  void Ell_BrBrConfMetricRegistry(void (*function),const char *);
  void Ell_PolyConfMetricRegistry(void (*function),const char *);
  
  int err=0;
  
  err += Ell_RegisterEq(Ell_LinConfMetricRegistry, "Ell_LinConfMetric");
  err += Ell_RegisterEq(Ell_BrBrConfMetricRegistry,"Ell_BrBrConfMetric");
  err += Ell_RegisterEq(Ell_PolyConfMetricRegistry,"Ell_PolyConfMetric");
  err += Ell_RegisterEq(Ell_LinMetricRegistry,     "Ell_LinMetric");
  err += Ell_RegisterEq(Ell_LinFlatRegistry,       "Ell_LinFlat");

  err += Ell_CreateKey(CCTK_VARIABLE_STRING,
                      "EllLinFlat::Bnd");
  err += Ell_CreateKey(CCTK_VARIABLE_STRING,
                      "EllLinConfMetric::Bnd");
  err += Ell_CreateKey(CCTK_VARIABLE_STRING,
                      "EllLinMetric::Bnd");

  err += Ell_CreateKey(CCTK_VARIABLE_STRING,
                      "EllLinFlat::Bnd::Robin");
  err += Ell_CreateKey(CCTK_VARIABLE_STRING,
                      "EllLinConfMetric::Bnd::Robin");
  err += Ell_CreateKey(CCTK_VARIABLE_STRING,
                      "EllLinMetric::Bnd::Robin");
  
  /* Register the variables needed to use these boundaries */
  err += Ell_CreateKey(CCTK_VARIABLE_REAL, 
                      "EllLinConfMetric::Bnd::Robin::inf"); 
  err += Ell_CreateKey(CCTK_VARIABLE_INT,  
                      "EllLinConfMetric::Bnd::Robin::falloff"); 
  err += Ell_CreateKey(CCTK_VARIABLE_REAL, 
                      "EllLinConfMetric::Bnd::Const::V0");

  /* Register the variables needed to use these boundaries */
  err += Ell_CreateKey(CCTK_VARIABLE_REAL, 
                      "EllLinMetric::Bnd::Robin::inf"); 
  err += Ell_CreateKey(CCTK_VARIABLE_INT,  
                      "EllLinMetric::Bnd::Robin::falloff"); 
  err += Ell_CreateKey(CCTK_VARIABLE_REAL, "EllLinMetric::Bnd::Const::V0");

  /* Register the variables needed to use these boundaries */
  err += Ell_CreateKey(CCTK_VARIABLE_REAL, 
                      "EllLinFlat::Bnd::Robin::inf"); 
  err += Ell_CreateKey(CCTK_VARIABLE_INT,  
                      "EllLinFlat::Bnd::Robin::falloff"); 
  err += Ell_CreateKey(CCTK_VARIABLE_REAL, 
                      "EllLinFlat::Bnd::Const::V0");

  if (err<0) CCTK_WARN(1,"Error registering the basic elliptic classes");
  
  return 0;
}
