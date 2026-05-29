/*@@
  @file      EOS_Base.h
  @date      Tue Dec 14 22:18:46 1999
  @author    Tom Goodale
  @desc
  Header file for EOS basic functions
  @enddesc
  @versions $Header$
@@*/

#ifndef _EOS_BASE_H_
#define _EOS_BASE_H_

#ifdef __cplusplus
extern "C" {
#endif

int EOS_RegisterMethod(const char *name);
int EOS_Handle(const char *name);

#define EOS_REGISTER_FUNCTION(x)                                               \
  int EOS_Register##x(int, CCTK_REAL (*func)(CCTK_REAL, CCTK_REAL))
#define EOS_CALL_FUNCTION(x) CCTK_REAL EOS_##x(int, CCTK_REAL, CCTK_REAL)

EOS_REGISTER_FUNCTION(Pressure);
EOS_REGISTER_FUNCTION(SpecificIntEnergy);
EOS_REGISTER_FUNCTION(RestMassDens);
EOS_REGISTER_FUNCTION(DPressByDRho);
EOS_REGISTER_FUNCTION(DPressByDEps);

EOS_CALL_FUNCTION(Pressure);
EOS_CALL_FUNCTION(SpecificIntEnergy);
EOS_CALL_FUNCTION(RestMassDens);
EOS_CALL_FUNCTION(DPressByDRho);
EOS_CALL_FUNCTION(DPressByDEps);

#ifdef __cplusplus
}
#endif

#endif /* EOS_Base.h */
