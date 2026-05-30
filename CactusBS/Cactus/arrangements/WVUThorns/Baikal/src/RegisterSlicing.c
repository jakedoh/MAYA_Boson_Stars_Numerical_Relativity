
#include "cctk.h"

#include "Slicing.h"

int Baikal_RegisterSlicing (void)
{
  Einstein_RegisterSlicing ("Baikal");
  return 0;
}