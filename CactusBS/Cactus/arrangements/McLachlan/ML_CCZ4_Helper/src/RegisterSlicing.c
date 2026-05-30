#include <cctk.h>
#include <Slicing.h>

int ML_CCZ4_RegisterSlicing(void) {
  Einstein_RegisterSlicing("ML_CCZ4");
  return 0;
}
