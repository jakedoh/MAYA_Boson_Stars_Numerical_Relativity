
#include "cctk.h"

#include "Slicing.h"

int BaikalVacuum_RegisterSlicing (void)
{
  Einstein_RegisterSlicing ("BaikalVacuum");
  return 0;
}