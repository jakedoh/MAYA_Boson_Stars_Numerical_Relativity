/*@@
  @file      Register.c
  @date      Tue Dec 14 20:55:06 1999
  @author    Tom Goodale
  @desc
  Functions for registering and calling EOS routines
  @enddesc
@@*/

#include <stdlib.h>

#include "cctk.h"
#include "cctk_FortranString.h"

#include "StoreHandledData.h"

#include "EOS_Base.h"

static cHandledData *methods = NULL;

enum call_methods { none, c, f };

typedef struct {
  enum call_methods lang;

  CCTK_REAL (*f)(CCTK_REAL *, CCTK_REAL *);
  CCTK_REAL (*c)(CCTK_REAL, CCTK_REAL);
} func_t;

typedef struct {
  func_t Pressure;
  func_t SpecificIntEnergy;
  func_t RestMassDens;
  func_t DPressByDRho;
  func_t DPressByDEps;
} method_t;

/*@@
  @routine    EOS_RegisterMethod
  @date       Tue Dec 14 22:04:23 1999
  @author     Tom Goodale
  @desc
  Registers an EOS method
  @enddesc
  @calls
  @calledby
  @history

  @endhistory

@@*/
int EOS_RegisterMethod(const char *name) {
  int handle;

  method_t *method;

  method = (method_t *)malloc(sizeof(method_t));

  if (method) {
    method->Pressure.lang = none;
    method->SpecificIntEnergy.lang = none;
    method->RestMassDens.lang = none;
    method->DPressByDRho.lang = none;
    method->DPressByDEps.lang = none;

    handle = Util_NewHandle(&methods, name, (void *)method);
  } else {
    handle = -1;
  }

  return handle;
}

void CCTK_FCALL CCTK_FNAME(EOS_RegisterMethod)(int *handle, ONE_FORTSTRING_ARG);
void CCTK_FCALL
    CCTK_FNAME(EOS_RegisterMethod)(int *handle, ONE_FORTSTRING_ARG) {
  ONE_FORTSTRING_CREATE(name)

  *handle = EOS_RegisterMethod(name);

  free(name);
}

/*@@
  @routine    EOS_Handle
  @date       Tue Dec 14 22:04:47 1999
  @author     Tom Goodale
  @desc
  Gets the handle associated with an EOS method
  @enddesc
  @calls
  @calledby
  @history

  @endhistory

@@*/
int EOS_Handle(const char *name) {
  int handle;
  method_t method;

  handle = Util_GetHandle(methods, name, (void *)&method);

  return handle;
}

void CCTK_FCALL CCTK_FNAME(EOS_Handle)(int *handle, ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME(EOS_Handle)(int *handle, ONE_FORTSTRING_ARG) {
  ONE_FORTSTRING_CREATE(name)

  *handle = EOS_Handle(name);

  free(name);
}

#define REGISTER_FUNCTION(x)                                                   \
  int EOS_Register##x(int handle, CCTK_REAL (*func)(CCTK_REAL, CCTK_REAL)) {   \
    int retval;                                                                \
    method_t *method;                                                          \
                                                                               \
    method = (method_t *)Util_GetHandledData(methods, handle);                 \
                                                                               \
    if (method) {                                                              \
      method->x.lang = c;                                                      \
      method->x.c = func;                                                      \
      retval = 0;                                                              \
    } else {                                                                   \
      retval = 1;                                                              \
    }                                                                          \
                                                                               \
    return retval;                                                             \
  }

/* Functions to register C EOS functions */
REGISTER_FUNCTION(Pressure)
REGISTER_FUNCTION(SpecificIntEnergy)
REGISTER_FUNCTION(RestMassDens)
REGISTER_FUNCTION(DPressByDRho)
REGISTER_FUNCTION(DPressByDEps)

#define REGISTER_FORTRAN_FUNCTION_ARGS                                         \
  (int *retval, int *handle, CCTK_REAL (*func)(CCTK_REAL *, CCTK_REAL *))
#define REGISTER_FORTRAN_FUNCTION(x)                                           \
  REGISTER_FORTRAN_FUNCTION_ARGS {                                             \
    method_t *method;                                                          \
                                                                               \
    method = (method_t *)Util_GetHandledData(methods, *handle);                \
                                                                               \
    if (method) {                                                              \
      method->x.lang = f;                                                      \
      method->x.f = func;                                                      \
      *retval = 0;                                                             \
    } else {                                                                   \
      *retval = 1;                                                             \
    }                                                                          \
  }

/* Functions to register Fortran EOS functions */
void CCTK_FCALL CCTK_FNAME(EOS_RegisterPressure) REGISTER_FORTRAN_FUNCTION_ARGS;
void CCTK_FCALL CCTK_FNAME(EOS_RegisterPressure)
    REGISTER_FORTRAN_FUNCTION(Pressure) void CCTK_FCALL
    CCTK_FNAME(EOS_RegisterSpecificIntEnergy) REGISTER_FORTRAN_FUNCTION_ARGS;
void CCTK_FCALL CCTK_FNAME(EOS_RegisterSpecificIntEnergy)
    REGISTER_FORTRAN_FUNCTION(SpecificIntEnergy) void CCTK_FCALL
    CCTK_FNAME(EOS_RegisterRestMassDens) REGISTER_FORTRAN_FUNCTION_ARGS;
void CCTK_FCALL CCTK_FNAME(EOS_RegisterRestMassDens)
    REGISTER_FORTRAN_FUNCTION(RestMassDens) void CCTK_FCALL
    CCTK_FNAME(EOS_RegisterDPressByDRho) REGISTER_FORTRAN_FUNCTION_ARGS;
void CCTK_FCALL CCTK_FNAME(EOS_RegisterDPressByDRho)
    REGISTER_FORTRAN_FUNCTION(DPressByDRho) void CCTK_FCALL
    CCTK_FNAME(EOS_RegisterDPressByDEps) REGISTER_FORTRAN_FUNCTION_ARGS;
void CCTK_FCALL CCTK_FNAME(EOS_RegisterDPressByDEps)
    REGISTER_FORTRAN_FUNCTION(DPressByDEps)

#define CALL_FUNC(x)                                                           \
  CCTK_REAL EOS_##x(int handle, CCTK_REAL a, CCTK_REAL b) {                    \
    CCTK_REAL retval;                                                          \
    method_t *method;                                                          \
                                                                               \
    method = (method_t *)Util_GetHandledData(methods, handle);                 \
                                                                               \
    if (method && method->x.lang != none) {                                    \
      if (method->x.lang == c) {                                               \
        retval = method->x.c(a, b);                                            \
      } else {                                                                 \
        retval = method->x.f(&a, &b);                                          \
      }                                                                        \
    } else {                                                                   \
      retval = -1;                                                             \
    }                                                                          \
    return retval;                                                             \
  }

    /* Functions to call EOS functions from C */
    CALL_FUNC(Pressure) CALL_FUNC(SpecificIntEnergy) CALL_FUNC(RestMassDens)
        CALL_FUNC(DPressByDRho) CALL_FUNC(DPressByDEps)

#define CALL_FORTRAN_FUNC_ARGS (int *handle, CCTK_REAL *a, CCTK_REAL *b)
#define CALL_FORTRAN_FUNC(x)                                                   \
  CALL_FORTRAN_FUNC_ARGS {                                                     \
    CCTK_REAL retval;                                                          \
    method_t *method;                                                          \
                                                                               \
    method = (method_t *)Util_GetHandledData(methods, *handle);                \
                                                                               \
    if (method && method->x.lang != none) {                                    \
      if (method->x.lang == c) {                                               \
        retval = method->x.c(*a, *b);                                          \
      } else {                                                                 \
        retval = method->x.f(a, b);                                            \
      }                                                                        \
    } else {                                                                   \
      retval = -1;                                                             \
    }                                                                          \
    return retval;                                                             \
  }

    /* Functions to call EOS functions from Fortran*/
    CCTK_REAL CCTK_FCALL CCTK_FNAME(EOS_Pressure) CALL_FORTRAN_FUNC_ARGS;
CCTK_REAL CCTK_FCALL CCTK_FNAME(EOS_Pressure)
    CALL_FORTRAN_FUNC(Pressure) CCTK_REAL CCTK_FCALL
    CCTK_FNAME(EOS_SpecificIntEnergy) CALL_FORTRAN_FUNC_ARGS;
CCTK_REAL CCTK_FCALL CCTK_FNAME(EOS_SpecificIntEnergy)
    CALL_FORTRAN_FUNC(SpecificIntEnergy) CCTK_REAL CCTK_FCALL
    CCTK_FNAME(EOS_RestMassDens) CALL_FORTRAN_FUNC_ARGS;
CCTK_REAL CCTK_FCALL CCTK_FNAME(EOS_RestMassDens)
    CALL_FORTRAN_FUNC(RestMassDens) CCTK_REAL CCTK_FCALL
    CCTK_FNAME(EOS_DPressByDRho) CALL_FORTRAN_FUNC_ARGS;
CCTK_REAL CCTK_FCALL CCTK_FNAME(EOS_DPressByDRho)
    CALL_FORTRAN_FUNC(DPressByDRho) CCTK_REAL CCTK_FCALL
    CCTK_FNAME(EOS_DPressByDEps) CALL_FORTRAN_FUNC_ARGS;
CCTK_REAL CCTK_FCALL CCTK_FNAME(EOS_DPressByDEps)
    CALL_FORTRAN_FUNC(DPressByDEps)
