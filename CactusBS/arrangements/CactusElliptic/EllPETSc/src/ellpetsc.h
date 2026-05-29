/*@@
  @header ellpetsc.h
  @author Gabrielle Allen
  @date 17th September 2001
  @desc
  Include file for thorn EllPETSc
  @enddesc
@@*/

#ifndef _ELLPETSC_H_
#define _ELLPETSC_H_ 1

#define XDM 0
#define XDP 1
#define YDM 2
#define YDP 3
#define ZDM 4
#define ZDP 5

/* PETSc 3.7 requires an "options" arguemnt */
#if ((PETSC_VERSION_MAJOR > 3) || \
     ((PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR >= 7)))
#define PetscOptionsSetValue(iname, value) PetscOptionsSetValue(NULL, iname, value)
#endif

#endif /* _ELLPETSC_H_ */
