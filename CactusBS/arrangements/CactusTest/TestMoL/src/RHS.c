#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

#include <assert.h>
#include <math.h>

/*****************************************************************************
 ************** Local types **************************************************
 *****************************************************************************/
enum RHStype {
  unity, linear, quadratic, cubic, quartic, quintic, sixthorder,
  seventhorder, eighthorder, ninthorder, exponential
};

/*****************************************************************************
 ************** Local data ***************************************************
 *****************************************************************************/
#define COEFF_SLOWGF 0.25
#define COEFF_GF 1.0
#define COEFF_GA 0.5

/*****************************************************************************
 ************** Local routine prototypes *************************************
 *****************************************************************************/
static enum RHStype get_rhs_type(const char *rhstype);
static CCTK_REAL eval_rhs(CCTK_REAL cctk_time,
                          enum RHStype rhs_type,
                          CCTK_REAL coeff);
static CCTK_REAL eval_solution(CCTK_REAL cctk_time,
                               enum RHStype rhs_type,
                               CCTK_REAL coeff);
static int get_mol_slowstep(const cGH *cctkGH);

/*****************************************************************************
 ************** External rtouines  *******************************************
 *****************************************************************************/
void TestMoL_RHS_GF(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTS;

  const enum RHStype rhs_type = get_rhs_type(RHSexpression);

#pragma omp parallel
  CCTK_LOOP3_ALL(TestMoL_RHS, cctkGH, i,j,k)
  {
    const int idx = CCTK_GFINDEX3D(cctkGH, i,j,k);
    rhs_gf[idx] = eval_rhs(cctk_time, rhs_type, COEFF_GF);
  } CCTK_ENDLOOP3_ALL(TestMoL_RHS);

  // TODO: add code to compute the constrained vars
  // TODO: add code that uses the sandr vars to check that MoL correctly saves
  // and restores them
}

void TestMoL_RHSSlow_GF(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTS;

  if(!get_mol_slowstep(cctkGH))
    return;

  const enum RHStype rhsslow_type = get_rhs_type(RHSSlowexpression);

#pragma omp parallel
  CCTK_LOOP3_ALL(TestMoL_RHSSlow, cctkGH, i,j,k)
  {
    const int idx = CCTK_GFINDEX3D(cctkGH, i,j,k);
    rhsslow_gf[idx] = eval_rhs(cctk_time, rhsslow_type, COEFF_SLOWGF);
  } CCTK_ENDLOOP3_ALL(TestMoL_RHSSlow);
}

void TestMoL_RHS_GA(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTS;

  cGroupDynamicData groupdata;

  const int group = CCTK_GroupIndex(CCTK_THORNSTRING "::evolved_ga");
  const int ierr = CCTK_GroupDynamicData(cctkGH, group, &groupdata);
  if(ierr != 0)
  {
    CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                "Failed to obtain information about evolved_ga group: %d",
                ierr);
  }
  assert(groupdata.dim == 1);

  const enum RHStype rhs_type = get_rhs_type(RHSexpression);
  const CCTK_REAL coeff = COEFF_GA;
  for(int i = 0 ; i < groupdata.lsh[0] ; i++)
    rhs_ga[i] = eval_rhs(cctk_time, rhs_type, coeff);

  // TODO: add code to compute the constrained vars
}

void TestMoL_Compare_GF(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTS;

  const enum RHStype rhs_type = get_rhs_type(RHSexpression);
  const enum RHStype rhsslow_type = get_rhs_type(RHSSlowexpression);

#pragma omp parallel
  CCTK_LOOP3_ALL(TestMoL_Compare, cctkGH, i,j,k)
  {
    const int idx = CCTK_GFINDEX3D(cctkGH, i,j,k);
    analytic_gf[idx] = eval_solution(cctk_time, rhs_type, COEFF_GF);
    analyticslow_gf[idx] = eval_solution(cctk_time, rhsslow_type, COEFF_SLOWGF);
    diff_gf[idx] = evolved_gf[idx] - analytic_gf[idx];
    diffslow_gf[idx] = evolvedslow_gf[idx] - analyticslow_gf[idx];
  } CCTK_ENDLOOP3_ALL(TestMoL_Compare);
}

void TestMoL_Compare_GA(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTS;

  cGroupDynamicData groupdata;

  const int group = CCTK_GroupIndex(CCTK_THORNSTRING "::evolved_ga");
  const int ierr = CCTK_GroupDynamicData(cctkGH, group, &groupdata);
  if(ierr != 0)
  {
    CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                "Failed to obtain information about evolved_ga group: %d",
                ierr);
  }
  assert(groupdata.dim == 1);

  const enum RHStype rhs_type = get_rhs_type(RHSexpression);
  const CCTK_REAL coeff = COEFF_GA;
  for(int i = 0 ; i < groupdata.lsh[0] ; i++)
  {
    analytic_ga[i] = eval_solution(cctk_time, rhs_type, coeff);
    diff_ga[i] = evolved_ga[i] - analytic_ga[i];
  }

  // TODO: add code to compute the constrained vars
}

/*****************************************************************************
 ************** Local Routines ***********************************************
 *****************************************************************************/
static enum RHStype get_rhs_type(const char *rhstype)
{
  enum RHStype retval;

  if(CCTK_EQUALS(rhstype, "1"))
    retval = unity;
  else if(CCTK_EQUALS(rhstype, "t"))
    retval = linear;
  else if(CCTK_EQUALS(rhstype, "t**2"))
    retval = quadratic;
  else if(CCTK_EQUALS(rhstype, "t**3"))
    retval = cubic;
  else if(CCTK_EQUALS(rhstype, "t**4"))
    retval = quartic;
  else if(CCTK_EQUALS(rhstype, "t**5"))
    retval = quintic;
  else if(CCTK_EQUALS(rhstype, "t**6"))
    retval = sixthorder;
  else if(CCTK_EQUALS(rhstype, "t**7"))
    retval = seventhorder;
  else if(CCTK_EQUALS(rhstype, "t**8"))
    retval = eighthorder;
  else if(CCTK_EQUALS(rhstype, "t**9"))
    retval = ninthorder;
  else if(CCTK_EQUALS(rhstype, "exp(t)"))
    retval = exponential;
  else
    CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                "Unknown rhs expression '%s'", rhstype);
  return retval;
}


static CCTK_REAL eval_solution(CCTK_REAL cctk_time,
                               enum RHStype rhs_type,
                               CCTK_REAL coeff)
{
  CCTK_REAL retval;

  switch(rhs_type)
  {
    case unity:
      retval = -1. + coeff*cctk_time;
      break;
    case linear:
      retval = -1. + coeff*1./2*pow(cctk_time, 2);
      break;
    case quadratic:
      retval = -1. + coeff*1./3*pow(cctk_time, 3);
      break;
    case cubic:
      retval = -1. + coeff*1./4*pow(cctk_time, 4);
      break;
    case quartic:
      retval = -1. + coeff*1./5*pow(cctk_time, 5);
      break;
    case quintic:
      retval = -1. + coeff*1./6*pow(cctk_time, 6);
      break;
    case sixthorder:
      retval = -1. + coeff*1./7*pow(cctk_time, 7);
      break;
    case seventhorder:
      retval = -1. + coeff*1./8*pow(cctk_time, 8);
      break;
    case eighthorder:
      retval = -1. + coeff*1./9*pow(cctk_time, 9);
      break;
    case ninthorder:
      retval = -1. + coeff*1./10*pow(cctk_time, 10);
      break;
    case exponential:
      retval = -1. + coeff*exp(cctk_time);
      break;
    default:
      CCTK_ERROR("Unexpected rhs type. This should not happen");
      break;
  }

  return retval;
}
static CCTK_REAL eval_rhs(CCTK_REAL cctk_time,
                          enum RHStype rhs_type,
                          CCTK_REAL coeff)
{
  CCTK_REAL retval;

  switch(rhs_type)
  {
    case unity:
      retval = coeff;
      break;
    case linear:
      retval = coeff*cctk_time;
      break;
    case quadratic:
      retval = coeff*pow(cctk_time, 2);
      break;
    case cubic:
      retval = coeff*pow(cctk_time, 3);
      break;
    case quartic:
      retval = coeff*pow(cctk_time, 4);
      break;
    case quintic:
      retval = coeff*pow(cctk_time, 5);
      break;
    case sixthorder:
      retval = coeff*pow(cctk_time, 6);
      break;
    case seventhorder:
      retval = coeff*pow(cctk_time, 7);
      break;
    case eighthorder:
      retval = coeff*pow(cctk_time, 8);
      break;
    case ninthorder:
      retval = coeff*pow(cctk_time, 9);
      break;
    case exponential:
      retval = coeff*exp(cctk_time);
      break;
    default:
      CCTK_ERROR("Unexpected rhs type. This should not happen");
      break;
  }

  return retval;
}

static int get_mol_slowstep(const cGH *cctkGH)
{
  CCTK_INT *slowstep;
  slowstep = CCTK_VarDataPtr(cctkGH, 0, "MoL::MoL_SlowStep");
  assert(slowstep);
  return *slowstep;
}
