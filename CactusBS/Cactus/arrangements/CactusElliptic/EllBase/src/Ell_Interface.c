 /*@@
   @header    Ell_Interface.h
   @date      
   @author    Gerd Lanferman
   @desc 
   Elliptic class routines for:

   * Registering the equation class wrapper (the function which is called for
   a specific class of problems by passing all the necessay arguments PLUS
   the name of the desired solver 

   * Equation class wrapper, for each elliptic class 
   Derives the function to call from the passed registration name of the 
   solver "sname".

   @enddesc 
   @version $Header$
 @@*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_FortranString.h"     
#include "StoreNamedData.h"  

#include "EllBase.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(CactusElliptic_EllBase_Ell_Interface_c)

static pNamedData *LinConfMetricSolverDB;
static pNamedData *LinMetricSolverDB;
static pNamedData *LinFlatSolverDB;
static pNamedData *BrBrConfMetricSolverDB;
static pNamedData *PolyConfMetricSolverDB;

int Ell_LinFlatRegistry(void (*function), const char *sname);
int Ell_LinConfMetricRegistry(int (*function), const char *sname);
int Ell_LinMetricRegistry(void (*function), const char *sname);
int Ell_BrBrConfMetricRegistry(void (*function), const char *sname);
int Ell_PolyConfMetricRegistry(void (*function), const char *sname);

int Ell_LinConfMetricSolver(cGH *GH, 
                            int *MetricPsi, 
                            int FieldIndex, 
                            int MIndex, 
                            int NIndex, 
                            CCTK_REAL *AbsTol, 
                            CCTK_REAL *RelTol, 
                            const char *sname);
int Ell_LinMetricSolver(cGH *GH, 
                        int *Metric, 
                        int FieldIndex, 
                        int MIndex, 
                        int NIndex, 
                        CCTK_REAL *AbsTol, 
                        CCTK_REAL *RelTol, 
                        const char *sname);
int Ell_BrBrConfMetricSolver(cGH *GH, 
                             int *MetricPsi, 
                             int FieldIndex, 
                             int MIndex, 
                             int NIndex, 
                             CCTK_REAL *AbsTol, 
                             CCTK_REAL *RelTol, 
                             const char *sname);
int Ell_PolyConfMetricSolver(cGH *GH, 
                             int *MetricPsi, 
                             int FieldIndex, 
                             int *PIndex, 
                             int Pcount, 
                             CCTK_REAL *AbsTol, 
                             CCTK_REAL *RelTol, 
                             const char *sname);
int Ell_LinFlatSolver(cGH *GH, 
                      int FieldIndex, 
                      int MIndex, 
                      int NIndex, 
                      CCTK_REAL *AbsTol, 
                      CCTK_REAL *RelTol, 
                      const char *sname);


void CCTK_FCALL CCTK_FNAME(Ell_LinConfMetricSolver)
     (int *ierr, 
      cGH **GH, 
      int *MetricPsi, 
      int *FieldIndex, 
      int *MIndex, 
      int *NIndex, 
      CCTK_REAL *AbsTol, 
      CCTK_REAL *RelTol, 
      ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME(Ell_LinMetricSolver)
     (int *ierr, 
      cGH **GH, 
      int *Metric, 
      int *FieldIndex, 
      int *MIndex, 
      int *NIndex, 
      CCTK_REAL *AbsTol, 
      CCTK_REAL *RelTol, 
      ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME(Ell_LinFlatSolver)
     (int *ierr, 
      cGH **GH, 
      int *FieldIndex, 
      int *MIndex, 
      int *NIndex, 
      CCTK_REAL *AbsTol, 
      CCTK_REAL *RelTol, 
      ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME(Ell_BrBrConfMetricSolver)
     (int *ierr, 
      cGH **GH, 
      int *MetricPsi, 
      int *FieldIndex, 
      int *MIndex, 
      int *NIndex, 
      CCTK_REAL *AbsTol, 
      CCTK_REAL *RelTol, 
      ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME(Ell_PolyConfMetricSolver)
     (int *ierr, 
      cGH **GH, 
      int *MetricPsi, 
      int *FieldIndex, 
      int *PIndex, 
      int *Pcount, 
      CCTK_REAL *AbsTol, 
      CCTK_REAL *RelTol, 
      ONE_FORTSTRING_ARG);



/*
######################################################
###### Elliptic Equation class: LinEllConfMetric #####
######################################################
*/

int Ell_LinConfMetricRegistry(int (*function), const char *sname) 
{
  int retval;

  if(!GetNamedData(LinConfMetricSolverDB,sname))
  {
    StoreNamedData(&LinConfMetricSolverDB,sname,(int*)function);
    retval = ELL_SUCCESS;
  }
  else
  {
    retval = ELL_SOLVEREXISTS;
  }
  
  return retval;
}

int Ell_LinConfMetricSolver(cGH *GH, 
                            int *MetricPsi, 
                            int FieldIndex, 
                            int MIndex, 
                            int NIndex, 
                            CCTK_REAL *AbsTol, 
                            CCTK_REAL *RelTol, 
                            const char *sname) 
{

  int retval=ELL_SUCCESS;
  int (*fn)(cGH *,int *,int,int,int,CCTK_REAL *,CCTK_REAL *);
 
  fn = (int(*)(cGH*,int*,int,int,int,CCTK_REAL*,CCTK_REAL*))
        (GetNamedData(LinConfMetricSolverDB,sname));

  if (fn)
  {
    retval = fn(GH, MetricPsi, FieldIndex, MIndex, NIndex, AbsTol, RelTol);
  }    
  else  
  { 
    CCTK_WARN(2,"Cannot find solver for LinEllConfMetric");
    retval = ELL_NOSOLVER;
  }

  return retval;

}

void CCTK_FCALL CCTK_FNAME(Ell_LinConfMetricSolver)
     (int *ierr, 
      cGH **GH, 
      int *MetricPsi, 
      int *FieldIndex, 
      int *MIndex, 
      int *NIndex, 
      CCTK_REAL *AbsTol, 
      CCTK_REAL *RelTol, 
      ONE_FORTSTRING_ARG) 
{

  ONE_FORTSTRING_CREATE(sname); 

  *ierr = Ell_LinConfMetricSolver(*GH, 
                                  MetricPsi, 
                                  *FieldIndex, 
                                  *MIndex, 
                                  *NIndex, 
                                  AbsTol, 
                                  RelTol, 
                                  sname);
  free(sname);

}



/*
##################################################
###### Elliptic Equation class: LinEllMetric #####
##################################################
*/

int Ell_LinMetricRegistry(void (*function), const char *sname) 
{
  int retval;

  if(!GetNamedData(LinMetricSolverDB,sname))
  {
    StoreNamedData(&LinMetricSolverDB,sname,(void*)function);
    retval = ELL_SUCCESS;
  }
  else
  {
    retval = ELL_SOLVEREXISTS;
  }
  
  return retval;
}

int Ell_LinMetricSolver(cGH *GH, 
                        int *Metric, 
                        int FieldIndex, 
                        int MIndex, 
                        int NIndex, 
                        CCTK_REAL *AbsTol, 
                        CCTK_REAL *RelTol, 
                        const char *sname) 
{

  int retval=ELL_SUCCESS;
  int (*fn)(cGH *,int *,int,int,int,CCTK_REAL *,CCTK_REAL *);
  
  fn = (int(*)(cGH*,int*,int,int,int,CCTK_REAL*,CCTK_REAL*))
    (GetNamedData(LinMetricSolverDB,sname));

  if (fn)
  {
    retval = fn(GH, Metric, FieldIndex, MIndex, NIndex, AbsTol, RelTol);
  }    
  else  
  { 
    CCTK_WARN(2,"Cannot find solver for LinEllMetric");
    retval = ELL_NOSOLVER;
  }

  return(retval);
}
  
void CCTK_FCALL CCTK_FNAME(Ell_LinMetricSolver)
     (int *ierr, 
      cGH **GH, 
      int *Metric, 
      int *FieldIndex, 
      int *MIndex, 
      int *NIndex, 
      CCTK_REAL *AbsTol, 
      CCTK_REAL *RelTol, 
      ONE_FORTSTRING_ARG) 
{

  ONE_FORTSTRING_CREATE(sname); 

  *ierr = Ell_LinMetricSolver(*GH, 
                              Metric, 
                              *FieldIndex, 
                              *MIndex, 
                              *NIndex, 
                              AbsTol, 
                              RelTol, 
                              sname);
  free(sname);
}
  


/*
################################################
###### Elliptic Equation class: LinEllFlat #####
################################################
*/

int Ell_LinFlatRegistry(void (*function), const char *sname) 
{
  int retval;

  if(!GetNamedData(LinFlatSolverDB,sname))
  {
    StoreNamedData(&LinFlatSolverDB,sname,(void*)function);
    retval = ELL_SUCCESS;
  }
  else
  {
    retval = ELL_SOLVEREXISTS;
  }
  
  return retval;

}

int Ell_LinFlatSolver(cGH *GH, 
                      int FieldIndex, 
                      int MIndex, 
                      int NIndex, 
                      CCTK_REAL *AbsTol, 
                      CCTK_REAL *RelTol, 
                      const char *sname) 
{

  int retval=ELL_SUCCESS;
  int (*fn)(cGH *,int,int,int,CCTK_REAL *,CCTK_REAL *);
  
  fn = (int (*)(cGH *,int,int,int,CCTK_REAL *,CCTK_REAL *)) 
    (GetNamedData(LinFlatSolverDB,sname));

  if (fn)
  {
    retval = fn(GH, FieldIndex, MIndex, NIndex, AbsTol, RelTol);
  }    
  else  
  { 
    CCTK_WARN(2,"Ell_LinFlatSolver: Cannot find solver for LinEllFlat");
    retval = ELL_NOSOLVER;
  }
  
  return retval;
}
  
void CCTK_FCALL CCTK_FNAME(Ell_LinFlatSolver)
     (int *ierr, 
      cGH **GH, 
      int *FieldIndex, 
      int *MIndex, 
      int *NIndex, 
      CCTK_REAL *AbsTol, 
      CCTK_REAL *RelTol, 
      ONE_FORTSTRING_ARG) 
{

  ONE_FORTSTRING_CREATE(sname); 
  *ierr = Ell_LinFlatSolver(*GH, 
                            *FieldIndex, 
                            *MIndex, 
                            *NIndex, 
                            AbsTol, 
                            RelTol, 
                            sname);
  free(sname);

}


/*
####################################################
###### Elliptic Equation class: BrBrConfMetric #####
####################################################
(Brandt-Bruemann Data with conformal metric) 
*/

int Ell_BrBrConfMetricRegistry(void (*function), const char *sname) 
{
  int retval;

  if(!GetNamedData(BrBrConfMetricSolverDB,sname))
  {
    StoreNamedData(&BrBrConfMetricSolverDB,sname, (void*)function);
    retval = ELL_SUCCESS;
  }
  else
  {
    retval = ELL_SOLVEREXISTS;
  }
  
  return retval;

}

int Ell_BrBrConfMetricSolver(cGH *GH, 
                             int *MetricPsi, 
                             int FieldIndex, 
                             int MIndex, 
                             int NIndex, 
                             CCTK_REAL *AbsTol, 
                             CCTK_REAL *RelTol, 
                             const char *sname) 
{

  int retval=ELL_SUCCESS;
  int (*fn)(cGH *,int *,int,int,int,CCTK_REAL *,CCTK_REAL *);
 
  fn = (int(*)(cGH*,int*,int,int,int,CCTK_REAL*,CCTK_REAL*))
        (GetNamedData(BrBrConfMetricSolverDB,sname));

  if (fn)
  {
    retval = fn(GH, MetricPsi, FieldIndex, MIndex, NIndex, AbsTol, RelTol);
  }    
  else  
  { 
    CCTK_WARN(2,"Ell_BrBrConfMetricSolver: Cannot find solver for BrBrConfMetric");
    retval = ELL_NOSOLVER;
  }

  return retval;

}

void CCTK_FCALL CCTK_FNAME(Ell_BrBrConfMetricSolver)
     (int *ierr, 
      cGH **GH, 
      int *MetricPsi, 
      int *FieldIndex, 
      int *MIndex, 
      int *NIndex, 
      CCTK_REAL *AbsTol, 
      CCTK_REAL *RelTol, 
      ONE_FORTSTRING_ARG) 
{

  ONE_FORTSTRING_CREATE(sname); 

  *ierr = Ell_BrBrConfMetricSolver(*GH, 
                                   MetricPsi, 
                                   *FieldIndex, 
                                   *MIndex, 
                                   *NIndex, 
                                   AbsTol, 
                                   RelTol, 
                                   sname);

  free(sname);

}


/*
####################################################
###### Elliptic Equation class: PolyConfMetric #####
####################################################
*/

int Ell_PolyConfMetricRegistry(void (*function), const char *sname) 
{
  int retval;

  if(!GetNamedData(PolyConfMetricSolverDB,sname))
  {
    StoreNamedData(&PolyConfMetricSolverDB,sname, (void*)function);
    retval = ELL_SUCCESS;
  }
  else
  {
    retval = ELL_SOLVEREXISTS;
  }
  
  return retval;
}

int Ell_PolyConfMetricSolver(cGH *GH, 
                             int *MetricPsi, 
                             int FieldIndex, 
                             int *PIndex, 
                             int Pcount, 
                             CCTK_REAL *AbsTol, 
                             CCTK_REAL *RelTol, 
                             const char *sname) 
{

  int retval=ELL_SUCCESS;
  int (*fn)(cGH *,int *,int,int *,int,CCTK_REAL *,CCTK_REAL *);
 
  fn = (int(*)(cGH*,int*,int,int*,int,CCTK_REAL*,CCTK_REAL*))
        (GetNamedData(PolyConfMetricSolverDB,sname));

  if (fn)
  {
    retval = fn(GH, MetricPsi, FieldIndex, PIndex, Pcount, AbsTol, RelTol);
  }    
  else  
  { 
    CCTK_WARN(2,"Cannot find solver for PolyConfMetric");
    retval = ELL_NOSOLVER;
  }
  
  return retval;

}

void CCTK_FCALL CCTK_FNAME(Ell_PolyConfMetricSolver)
     (int *ierr, 
      cGH **GH, 
      int *MetricPsi, 
      int *FieldIndex, 
      int *PIndex, 
      int *Pcount, 
      CCTK_REAL *AbsTol, 
      CCTK_REAL *RelTol, 
      ONE_FORTSTRING_ARG) 
{

  ONE_FORTSTRING_CREATE(sname); 
  *ierr = Ell_PolyConfMetricSolver(*GH, 
                                   MetricPsi, 
                                   *FieldIndex, 
                                   PIndex, 
                                   *Pcount, 
                                   AbsTol, 
                                   RelTol, 
                                   sname);
  free(sname);

}
