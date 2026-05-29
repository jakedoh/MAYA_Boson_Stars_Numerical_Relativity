#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

void ADMMass_InitLoopCounter(CCTK_ARGUMENTS);

void ADMMass_SetLoopCounter(CCTK_ARGUMENTS);

void ADMMass_Loop(CCTK_ARGUMENTS);


/* Initialise the loop counter */
void ADMMass_InitLoopCounter(CCTK_ARGUMENTS)
{
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

          *ADMMass_LoopCounter = 0;
}

/* Set the loop counter to the value of the parameter ADMMass:ADMMass_number */
void ADMMass_SetLoopCounter(CCTK_ARGUMENTS)
{
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

      *ADMMass_LoopCounter = ADMMass_number;
}

/* Decrements the counter to loop over all radii/distances set */
void ADMMass_Loop(CCTK_ARGUMENTS)
{
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    (*ADMMass_LoopCounter)--;

}
