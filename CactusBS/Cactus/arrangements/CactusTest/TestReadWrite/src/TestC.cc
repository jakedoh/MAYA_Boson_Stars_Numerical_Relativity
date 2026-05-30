#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Functions.h>
#include <cctk_Parameters.h>
#include <assert.h>
#include <iostream>
#include <string>

#include <type_traits>

extern "C"
void TestReadWrite_SelectBCs (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;

  const char *groupnames[] = {
    "TestReadWriteImp::Var1",
    "TestReadWrite::Var2",
    "TestReadWriteImp::Var3",
    "TestReadWriteImp::TestGroup",
    NULL
  };
  for(int i = 0; groupnames[i]; i++) {
    int ierr = Driver_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                                       groupnames[i], "none");
    assert(!ierr);
  }
}

extern "C"
void TestReadWrite_Reset (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;

  const char *groupnames[] = {
    "TestReadWriteImp::Var1",
    "TestReadWriteImp::rad1",
    "TestReadWriteImp::rad2",
    "TestReadWriteImp::rad3",
    "TestReadWrite::Var2",
    "TestReadWrite::ScalarVar",
    "TestReadWriteImp::Var3",
    "TestReadWriteImp::TestGroup",
    "TestReadWriteImp::UnusedVar",
    "TestReadWriteImp::gridArrayGroup",
    NULL
  };
  for(int i = 0; groupnames[i]; i++) {
    int const firstvarindex = CCTK_FirstVarIndex(groupnames[i]);
    assert(firstvarindex >= 0);
    int const numvars = CCTK_NumVarsInGroup(groupnames[i]);
    assert(numvars >= 0);
    int const tls = CCTK_ActiveTimeLevelsGN(cctkGH, groupnames[i]);
    for(int var = 0; var < numvars; var++) {
      for(int tl = 0; tl < tls; tl++) {
        Driver_SetValidRegion(cctkGH, firstvarindex + var, tl,
                              CCTK_VALID_NOWHERE);
      }
    }
  }
}

extern "C"
void TestReadWrite_TestC_A (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS_CHECKED(TestReadWrite_TestC_A);
  DECLARE_CCTK_PARAMETERS;

  static int first_time = 1;

  CCTK_REAL UnusedVar = -1;
  (void)UnusedVar;

  // test that same function scheduled twice has different accessible variables
  if(first_time) {
    groupVar1[0] = 43;
    groupVar2[0] = 44;

    Var3[0] = 46;
    Var2[0] = x[0];

    groupVar1_p_p[0] = 53;
    groupVar2_p_p[0] = 54;

    gridArray1[0] = 62.;

    rad1[1*radx + 0] = 63.;

    *ScalarVar = 74.;

    int vi;
    *ValidRegionsC = 1;

    vi = CCTK_VarIndex("TestReadWriteImp::groupVar1[0]");
    assert(vi >= 0);
    *ValidRegionsC &=
      (Driver_GetValidRegion(cctkGH, vi, 0) == CCTK_VALID_NOWHERE);

    vi = CCTK_VarIndex("TestReadWriteImp::groupVar1[0]");
    assert(vi >= 0);
    *ValidRegionsC &=
      (Driver_GetValidRegion(cctkGH, vi, 2) == CCTK_VALID_NOWHERE);

    vi = CCTK_VarIndex("TestReadWriteImp::groupVar2[0]");
    assert(vi >= 0);
    *ValidRegionsC &=
      (Driver_GetValidRegion(cctkGH, vi, 2) == CCTK_VALID_NOWHERE);

    vi = CCTK_VarIndex("TestReadWriteImp::Var3[0]");
    assert(vi >= 0);
    *ValidRegionsC &=
      (Driver_GetValidRegion(cctkGH, vi, 0) == CCTK_VALID_NOWHERE);

    vi = CCTK_VarIndex("TestReadWriteImp::GridArray1");
    assert(vi >= 0);
    *ValidRegionsC &=
      (Driver_GetValidRegion(cctkGH, vi, 0) == CCTK_VALID_NOWHERE);

    vi = CCTK_VarIndex("TestReadWrite::ScalarVar");
    assert(vi >= 0);
    *ValidRegionsC &=
      (Driver_GetValidRegion(cctkGH, vi, 0) == CCTK_VALID_NOWHERE);

    first_time = 0;
  } else {
    Var1[0] = 42.;
    Var1_p[0] = 52.;

    int vi;

    vi = CCTK_VarIndex("TestReadWriteImp::Var1");
    assert(vi >= 0);
    *ValidRegionsC &=
      (Driver_GetValidRegion(cctkGH, vi, 0) == CCTK_VALID_NOWHERE);
  }

  CCTK_REAL gridArray1_p;
  (void)gridArray1_p;
  // Fortran test checks size of arrays
  (void)gridArray1_p_p;
}

extern "C"
void TestReadWrite_TestC_B (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS_TestReadWrite_TestC_B;

  Var2[0] = Var3[0] - 1;

  *UnusedVarIsNullC = Var1 == NULL && Var1_p == NULL && gridArray1_p_p == NULL;
  *UnusedGroupIsNullC = groupVar1 == NULL && groupVar2 == NULL;

  int vi;

  vi = CCTK_VarIndex("TestReadWriteImp::groupVar1[0]");
  assert(vi >= 0);
  *ValidRegionsC &=
    (Driver_GetValidRegion(cctkGH, vi, 0) == CCTK_VALID_EXTERIOR);

  vi = CCTK_VarIndex("TestReadWriteImp::groupVar1[0]");
  assert(vi >= 0);
  *ValidRegionsC &=
    (Driver_GetValidRegion(cctkGH, vi, 2) == CCTK_VALID_EVERYWHERE);

  vi = CCTK_VarIndex("TestReadWriteImp::groupVar2[0]");
  assert(vi >= 0);
  *ValidRegionsC &=
    (Driver_GetValidRegion(cctkGH, vi, 2) == CCTK_VALID_EVERYWHERE);

  vi = CCTK_VarIndex("TestReadWriteImp::Var3[0]");
  assert(vi >= 0);
  *ValidRegionsC &=
    (Driver_GetValidRegion(cctkGH, vi, 0) == CCTK_VALID_EVERYWHERE);

  vi = CCTK_VarIndex("TestReadWriteImp::GridArray1");
  assert(vi >= 0);
  *ValidRegionsC &=
    (Driver_GetValidRegion(cctkGH, vi, 0) == CCTK_VALID_INTERIOR);

  vi = CCTK_VarIndex("TestReadWriteImp::Var1");
  assert(vi >= 0);
  *ValidRegionsC &=
    (Driver_GetValidRegion(cctkGH, vi, 0) == CCTK_VALID_EVERYWHERE);

  vi = CCTK_VarIndex("TestReadWrite::ScalarVar");
  assert(vi >= 0);
  *ValidRegionsC &=
    (Driver_GetValidRegion(cctkGH, vi, 0) == CCTK_VALID_EVERYWHERE);

}

extern "C"
void TestReadWrite_TestC_C(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS_TestReadWrite_TestC_C;
  DECLARE_CCTK_PARAMETERS;

  *VarCurrentLevelC = Var1[0] == 42. && Var2[0] == 45. && Var3[0] == 46;
  *GroupCurrentLevelC = groupVar1[0] == 43 && groupVar2[0] == 44;

  *VarPastLevelC = Var1_p[0] == 52.;
  *GroupPastLevelC = groupVar1_p_p[0] == 53 && groupVar2_p_p[0] == 54;

  *ArrayCurrentLevelC = gridArray1[0] == 62. && rad1[1*radx + 0] == 63.;

  *ScalarVarCurrentLevelC = *ScalarVar == 74.;

  static_assert(!std::is_const<std::remove_reference<decltype(*VarCurrentLevelC)>::type>::value,
                "VarCurrentLevelC is const");

  static_assert(std::is_const<std::remove_reference<decltype(*Var1)>::type>::value,
                "Var1) is not const");
  static_assert(std::is_const<std::remove_reference<decltype(*Var1_p)>::type>::value,
                "Var1_p is not const");
  static_assert(std::is_const<std::remove_reference<decltype(*groupVar1)>::type>::value &&
                std::is_const<std::remove_reference<decltype(*groupVar2)>::type>::value,
                "groupVar1 or groupVar2 are not const");
  static_assert(std::is_const<std::remove_reference<decltype(*groupVar1_p_p)>::type>::value &&
                std::is_const<std::remove_reference<decltype(*groupVar2_p_p)>::type>::value,
                "groupVar1_p_p or groupVar2_p_p are not const");
}
