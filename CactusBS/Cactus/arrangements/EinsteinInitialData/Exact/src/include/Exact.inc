c Exact.inc -- Misc Fortran 77/90 definitions for the whole thorn
c $Header$

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c
c This thorn is mainly written in a mixture of Fortran 77 and Fortran 90;
c a few routines are written in C.  At present Fortran 90 is used only
c for the "arbitrary slice" option (see the thorn guide for details).
c If this option is not needed, then the Fortran 90 code can all be
c #ifdef-ed out, allowing this thorn to be compiled on a system
c having only Fortran 77 and C compilers (i.e. no Fortran 90 compiler).
c To do this, uncomment the following line:
c #define EXACT_NO_F90

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
