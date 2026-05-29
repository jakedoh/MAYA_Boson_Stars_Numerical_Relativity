#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <cctk_File.h>

#include <util_String.h>

#include <stdio.h>
#include <stdlib.h>

void TestPar_SetValues (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS_TestPar_SetValues;
  DECLARE_CCTK_PARAMETERS;
  
  int i;
  FILE *fh = NULL;
  char *buf = NULL;
  
  if(CCTK_MyProc(cctkGH) != 0)
    return;

  if(Util_asprintf(&buf, "%s/ParamEval.asc", out_dir) == 0  || 
     CCTK_CreateDirectory(0777, out_dir) < 0 ||
     (fh = fopen(buf, "w+")) == NULL)
  {
    CCTK_VWarn(CCTK_WARN_ABORT, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Could not open file '%s' for output", buf);
  }

  fprintf(fh,"#  %5s %5s %8s %8s %5s %5s\n","int1","int2","real1","real2","bool1","bool2");
  for(i = 0 ; i < 10 ; i++)
  {
    fprintf(fh, "%2d %5d %5d %8g %8g %5d %5d\n", i, (int)int1[i], (int)int2[i],
            (double)real1[i], (double)real2[i], (int)bool1[i], (int)bool2[i]);
  }

  if(fh)
    fclose(fh);
  if(buf)
    free(buf);
}
