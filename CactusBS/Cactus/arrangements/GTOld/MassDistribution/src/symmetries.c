#include <string.h>
#include <assert.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "MassDistribution.h"


/*******************************
 ****** External Routines ******
 *******************************/
void query_symmetries(const char *applied_symmetries, int *ncopies, double symsigns[8][3])
{
  char sym_buf[5]; /* space for xyzr\0 */
  const char *syms;

  /* get string representation of symmetries in system */
  if(CCTK_Equals(applied_symmetries, "autoprobe"))
  {
    int nsyms;
    const CCTK_INT *sym_active;
    
    nsyms = 0;
    if((sym_active = CCTK_ParameterGet("reflection_x", "ReflectionSymmetry", NULL)) != NULL && *sym_active)
      sym_buf[nsyms++] = 'x';
    if((sym_active = CCTK_ParameterGet("reflection_y", "ReflectionSymmetry", NULL)) != NULL && *sym_active)
      sym_buf[nsyms++] = 'y';
    if((sym_active = CCTK_ParameterGet("reflection_z", "ReflectionSymmetry", NULL)) != NULL && *sym_active)
      sym_buf[nsyms++] = 'z';
    if(CCTK_IsThornActive("RotatingSymmetry180"))
      sym_buf[nsyms++] = 'r';
    sym_buf[nsyms] = '\0';

    syms = sym_buf;
  }
  else
  {
    syms = applied_symmetries;
  }

  /* create required number symmetry copies */
  symsigns[0][0] = symsigns[0][1] = symsigns[0][2] = 1.0;
  *ncopies = 1;
  for( ; syms[0] != '\0' ; syms++)
  {
    for(int i = 0 ; i < *ncopies ; i++)
    {
      switch(syms[0])
      {
        case 'x':
          symsigns[*ncopies+i][0] = -symsigns[i][0];
          symsigns[*ncopies+i][1] = +symsigns[i][1];
          symsigns[*ncopies+i][2] = +symsigns[i][2];
          break;
        case 'y':
          symsigns[*ncopies+i][0] = +symsigns[i][0];
          symsigns[*ncopies+i][1] = -symsigns[i][1];
          symsigns[*ncopies+i][2] = +symsigns[i][2];
          break;
        case 'z':
          symsigns[*ncopies+i][0] = +symsigns[i][0];
          symsigns[*ncopies+i][1] = +symsigns[i][1];
          symsigns[*ncopies+i][2] = -symsigns[i][2];
          break;
        case 'r':
          symsigns[*ncopies+i][0] = -symsigns[i][0];
          symsigns[*ncopies+i][1] = -symsigns[i][1];
          symsigns[*ncopies+i][2] = +symsigns[i][2];
          break;
        default:
          CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                     "Unkown symmetry condition '%c'.", syms[0]);
          break; /* NOTREACHED */
      } /* switch */
    } /* for i */

    *ncopies *= 2;
    assert(*ncopies <= 8);
  } /* for syms */
}
