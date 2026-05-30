#include <cctk.h>
#include <cctk_Parameters.h>
#include <util_Table.h>

#include <assert.h>

static void
set_group_tags (int const checkpoint,
                int const persistent,
                int const prolongate,
                char const * restrict const gn);

int
CL_BSSN_SetGroupTags (void)
{
  DECLARE_CCTK_PARAMETERS;
  
  int const checkpoint = timelevels > 1;
  set_group_tags (checkpoint, checkpoint, 1, "ADMBase::metric");
  set_group_tags (checkpoint, checkpoint, 1, "ADMBase::curv");
  set_group_tags (checkpoint, checkpoint, 1, "ADMBase::lapse");
  set_group_tags (checkpoint, checkpoint, 1, "ADMBase::shift");
  set_group_tags (0         , 0         , 0, "ADMBase::dtlapse");
  set_group_tags (checkpoint, checkpoint, 1, "ADMBase::dtshift");
  
  set_group_tags (checkpoint, checkpoint, 0, "CL_BSSN::CL_cons_dlog_confac");
  set_group_tags (checkpoint, checkpoint, 0, "CL_BSSN::CL_cons_dmetric");
  set_group_tags (checkpoint, checkpoint, 0, "CL_BSSN::CL_cons_dlapse");
  set_group_tags (checkpoint, checkpoint, 0, "CL_BSSN::CL_cons_dshift");
  set_group_tags (checkpoint, checkpoint, 0, "CL_BSSN::CL_cons_detg");
  set_group_tags (checkpoint, checkpoint, 0, "CL_BSSN::CL_cons_Gamma");
  set_group_tags (checkpoint, checkpoint, 0, "CL_BSSN::CL_cons_traceA");
  set_group_tags (checkpoint, checkpoint, 0, "CL_BSSN::CL_Ham");
  set_group_tags (checkpoint, checkpoint, 0, "CL_BSSN::CL_mom");
  
  int const rhs_checkpoint = rhs_timelevels > 1;
  set_group_tags (rhs_checkpoint, rhs_checkpoint, 0, "CL_BSSN::CL_log_confacrhs");
  set_group_tags (rhs_checkpoint, rhs_checkpoint, 0, "CL_BSSN::CL_dlog_confacrhs");
  set_group_tags (rhs_checkpoint, rhs_checkpoint, 0, "CL_BSSN::CL_metricrhs");
  set_group_tags (rhs_checkpoint, rhs_checkpoint, 0, "CL_BSSN::CL_dmetricrhs");
  set_group_tags (rhs_checkpoint, rhs_checkpoint, 0, "CL_BSSN::CL_Gammarhs");
  set_group_tags (rhs_checkpoint, rhs_checkpoint, 0, "CL_BSSN::CL_trace_curvrhs");
  set_group_tags (rhs_checkpoint, rhs_checkpoint, 0, "CL_BSSN::CL_curvrhs");
  set_group_tags (rhs_checkpoint, rhs_checkpoint, 0, "CL_BSSN::CL_lapserhs");
  set_group_tags (rhs_checkpoint, rhs_checkpoint, 0, "CL_BSSN::CL_dlapserhs");
  set_group_tags (rhs_checkpoint, rhs_checkpoint, 0, "CL_BSSN::CL_shiftrhs");
  set_group_tags (rhs_checkpoint, rhs_checkpoint, 0, "CL_BSSN::CL_dshiftrhs");
  set_group_tags (rhs_checkpoint, rhs_checkpoint, 0, "CL_BSSN::CL_dtshiftrhs");
  
  return 0;
}

static void
set_group_tags (int const checkpoint,
                int const persistent,
                int const prolongate,
                char const * restrict const gn)
{
  assert (gn);
  
  int const gi = CCTK_GroupIndex (gn);
  assert (gi >= 0);
  
  int const table = CCTK_GroupTagsTableI (gi);
  assert (table >= 0);
  
  if (! checkpoint) {
    int const ierr = Util_TableSetString (table, "no", "Checkpoint");
    assert (! ierr);
  }
  
  if (! persistent) {
    int const ierr = Util_TableSetString (table, "no", "Persistent");
    assert (! ierr);
  }
  
  if (! prolongate) {
    int const iret = Util_TableDeleteKey (table, "ProlongationParameter");
    assert (iret == 0 || iret == UTIL_ERROR_TABLE_NO_SUCH_KEY);
    int const ierr = Util_TableSetString (table, "none", "Prolongation");
    assert (ierr == 0 || ierr == 1); /* 0 means key did not exist before.  1 means key existed before and has now been reset */
  }
}
