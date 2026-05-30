//
// startup routines allocating eg. memory
//

#include "matterlibs.h"
#include "string.h"
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

CCTK_REAL matterlibs_poison;

void MatterLibs_Startup(CCTK_ARGUMENTS)
{
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;

    // generate poison value
    memset(&matterlibs_poison, poison_value, sizeof(matterlibs_poison));
}
