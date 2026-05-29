#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

subroutine TestReadWrite_TestF_A(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS_CHECKED(TestReadWrite_TestF_A)
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  logical, save :: first_time = .true.
  integer :: vi

  ! test that same function scheduled twice has different accessible variables
  if(first_time) then
    Var3(1,1,1,1) = 46.
    Var2(1,1,1) = x(1,1,1)

    groupVar1(1,1,1,1) = 43
    groupVar2(1,1,1,1) = 44

    groupVar1_p_p(1,1,1,1) = 53
    groupVar2_p_p(1,1,1,1) = 54

    gridArray1(1,1) = 62.

    rad1(2,1) = 63.

    ScalarVar = 74.

    if(SIZE(gridArray1_p_p,1) .eq. 12 .and. SIZE(gridArray2_p_p,2) .eq. 13 .and. &
       SIZE(rad1,1) .eq. radx .and. SIZE(rad1,2) .eq. rad_param .and. &
       SIZE(rad2,1) .eq. radx .and. SIZE(rad2,2) .eq. rady .and. &
       SIZE(rad2,3) .eq. rad_param .and. &
       SIZE(rad3,1) .eq. radx .and. SIZE(rad3,2) .eq. rady .and. &
       SIZE(rad3,3) .eq. radz .and. SIZE(rad3,4) .eq. rad_param) then
      ArraySizeF = 1
    else
      ArraySizeF = 0
    end if

    ValidRegionsF = 1

    call CCTK_VarIndex(vi, "TestReadWriteImp::groupVar1[0]")
    if(Driver_GetValidRegion(cctkGH, vi, 0) .ne. CCTK_VALID_NOWHERE) then
      ValidRegionsF = 0
    end if

    call CCTK_VarIndex(vi, "TestReadWriteImp::groupVar1[0]")
    if(Driver_GetValidRegion(cctkGH, vi, 2) .ne. CCTK_VALID_NOWHERE) then
      ValidRegionsF = 0
    end if

    call CCTK_VarIndex(vi, "TestReadWriteImp::groupVar2[0]")
    if(Driver_GetValidRegion(cctkGH, vi, 2) .ne. CCTK_VALID_NOWHERE) then
      ValidRegionsF = 0
    end if

    call CCTK_VarIndex(vi, "TestReadWriteImp::Var3[0]")
    if(Driver_GetValidRegion(cctkGH, vi, 0) .ne. CCTK_VALID_NOWHERE) then
      ValidRegionsF = 0
    end if

    call CCTK_VarIndex(vi, "TestReadWriteImp::GridArray1")
    if(Driver_GetValidRegion(cctkGH, vi, 0) .ne. CCTK_VALID_NOWHERE) then
      ValidRegionsF = 0
    end if

    call CCTK_VarIndex(vi, "TestReadWrite::ScalarVar")
    if(Driver_GetValidRegion(cctkGH, vi, 0) .ne. CCTK_VALID_NOWHERE) then
      ValidRegionsF = 0
    end if

    first_time = .false.
  else
    Var1(1,1,1) = 42.
    Var1_p(1,1,1) = 52.

    call CCTK_VarIndex(vi, "TestReadWriteImp::Var1")
    if(Driver_GetValidRegion(cctkGH, vi, 0) .ne. CCTK_VALID_NOWHERE) then
      ValidRegionsF = 0
    end if
  end if

end subroutine

subroutine TestReadWrite_TestF_B(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS_TESTREADWRITE_TESTF_B
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS

  integer :: vi

  CCTK_POINTER :: NULL

  NULL = CCTK_NullPointer()

  Var2(1,1,1) = Var3(1,1,1,1) - 1.

  if(CCTK_PointerTo(Var1) .eq. NULL .and. CCTK_PointerTo(Var1_p) .eq. NULL .and. &
    CCTK_PointerTo(gridArray1_p_p) .eq. NULL) then
    UnusedVarIsNullF = 1
  else
    UnusedVarIsNullF = 0
  end if
  if(CCTK_PointerTo(groupVar1) .eq. NULL .and. &
     CCTK_PointerTo(groupVar2) .eq. NULL) then
    UnusedGroupIsNullF = 1
  else
    UnusedGroupIsNullF = 0
  end if

  call CCTK_VarIndex(vi, "TestReadWriteImp::groupVar1[0]")
  if(Driver_GetValidRegion(cctkGH, vi, 0) .ne. CCTK_VALID_EXTERIOR) then
    ValidRegionsF = 0
  end if

  call CCTK_VarIndex(vi, "TestReadWriteImp::groupVar1[0]")
  if(Driver_GetValidRegion(cctkGH, vi, 2) .ne. CCTK_VALID_EVERYWHERE) then
    ValidRegionsF = 0
  end if

  call CCTK_VarIndex(vi, "TestReadWriteImp::groupVar2[0]")
  if(Driver_GetValidRegion(cctkGH, vi, 2) .ne. CCTK_VALID_EVERYWHERE) then
    ValidRegionsF = 0
  end if

  call CCTK_VarIndex(vi, "TestReadWriteImp::Var3[0]")
  if(Driver_GetValidRegion(cctkGH, vi, 0) .ne. CCTK_VALID_EVERYWHERE) then
    ValidRegionsF = 0
  end if

  call CCTK_VarIndex(vi, "TestReadWriteImp::GridArray1")
  if(Driver_GetValidRegion(cctkGH, vi, 0) .ne. CCTK_VALID_INTERIOR) then
    ValidRegionsF = 0
  end if

  call CCTK_VarIndex(vi, "TestReadWriteImp::Var1")
  if(Driver_GetValidRegion(cctkGH, vi, 0) .ne. CCTK_VALID_EVERYWHERE) then
    ValidRegionsF = 0
  end if

  call CCTK_VarIndex(vi, "TestReadWrite::ScalarVar")
  if(Driver_GetValidRegion(cctkGH, vi, 0) .ne. CCTK_VALID_EVERYWHERE) then
    ValidRegionsF = 0
  end if
end subroutine

subroutine TestReadWrite_TestF_C(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS_TestReadWrite_TestF_C
  DECLARE_CCTK_PARAMETERS

  if(Var1(1,1,1) .eq. 42. .and. Var2(1,1,1) .eq. 45. .and. Var3(1,1,1,1) .eq. 46.) then
    VarCurrentLevelF = 1
  else
    VarCurrentLevelF = 0
  end if
  if(groupVar1(1,1,1,1) .eq. 43 .and. groupVar2(1,1,1,1) .eq. 44) then
    GroupCurrentLevelF = 1
  else
    GroupCurrentLevelF = 0
  end if

  if(Var1_p(1,1,1) .eq. 52) then
    VarPastLevelF = 1
  else
    VarPastLevelF = 0
  end if
  if(groupVar1_p_p(1,1,1,1) .eq. 53 .and. &
     groupVar2_p_p(1,1,1,1) .eq. 54) then
    GroupPastLevelF = 1
  else
    GroupPastLevelF = 0
  end if

  if(gridArray1(1,1) .eq. 62. .and.  rad1(2,1) .eq. 63.) then
    ArrayCurrentLevelF = 1
  else
    ArrayCurrentLevelF = 0
  end if

  if(ScalarVar .eq. 74.) then
    ScalarVarCurrentLevelF = 1
  else
    ScalarVarCurrentLevelF = 0
  end if
end subroutine

! this routine must never actually execute since it messes with its call signature and passes a CCTK_REAL array to TRIM
subroutine TestReadWrite_TestF_D(CCTK_ARGUMENTS)

  implicit none

  ! a horrible hack to trigger compile errors if READS variables are not
  ! intent(in)
  ! the capitalizationg if "intent" matches the one prodced by rdwr.pl
  ! in the final replacement INTENT must not match any of the variants defined
  ! here
  ! using OPTIONAL for both READS only variables and unused variables means I
  ! cannot distinguish betweent the two, but both are supposed to be
  ! INTENT(IN) so this should be fine
#define cctk_iNteNt_iN INTENT(iN), optional
#define iNteNt(iN) cctk_iNteNt_/**/iN
#define cctk_intent_IN INTENT(IN), optional
#define intent(IN) cctk_intent_/**/IN
  DECLARE_CCTK_ARGUMENTS_TestReadWrite_TestF_D
  DECLARE_CCTK_PARAMETERS

  integer :: idx

  ! these will fail to compile unless the dummy arguments were marked as
  ! optional
  if(PRESENT(Var1)) then
    STOP
  end if

  if(PRESENT(Var1_p)) then
    STOP
  end if

  if(PRESENT(Var2)) then
    STOP
  end if

  if(PRESENT(groupVar1)) then
    STOP
  end if

  if(PRESENT(groupVar2)) then
    STOP
  end if

  if(PRESENT(groupVar1_p_p)) then
    STOP
  end if

  if(PRESENT(groupVar2_p_p)) then
    STOP
  end if

  if(PRESENT(UnusedVar)) then
    STOP
  end if

  idx = INDEX(UnusedVar, "")
end subroutine
