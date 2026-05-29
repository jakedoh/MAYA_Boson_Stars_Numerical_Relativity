
c/*@@
c  @file      EOS_Ideal_Fluid_Startup.F
c  @date      December 1999
c  @author    Mark Miller
c  @desc
c  Startup for EOS_Ideal_Fluid
c  @enddesc
c@@*/

#include "cctk.h"

c/*@@
c  @routine    EOS_Ideal_Fluid_Startup
c  @date       December 1999
c  @author     Mark Miller
c  @desc
c
c  @enddesc
c  @calls      EOS_RegisterMethod EOS_RegisterPressure EOS_RegisterSpecificInternalEnergy EOS_RegisterRestMassDens EOS_RegisterDPressByDRho EOS_RegisterDPressByDEps CCTK_WARN
c  @calledby
c  @history
c
c  @endhistory
c@@*/

      INTEGER FUNCTION EOS_Ideal_Fluid_Startup()
      implicit none
      integer handle,ierr
      external EOS_Ideal_Fluid_Pressure
      external EOS_Ideal_Fluid_SpecificIE
      external EOS_Ideal_Fluid_RestMassDens
      external EOS_Ideal_Fluid_DPressByDRho
      external EOS_Ideal_Fluid_DPressByDEps

      call EOS_RegisterMethod(handle,"Ideal_Fluid")

      IF(handle .ge. 0) THEN

        CALL EOS_RegisterPressure(ierr, handle, EOS_Ideal_Fluid_Pressure);
        CALL EOS_RegisterSpecificIntEnergy(ierr, handle, 
     .                                     EOS_Ideal_Fluid_SpecificIE);
        CALL EOS_RegisterRestMassDens(ierr, handle, 
     .                                EOS_Ideal_Fluid_RestMassDens);
        CALL EOS_RegisterDPressByDRho(ierr, handle, 
     .                                EOS_Ideal_Fluid_DPressByDRho);
        CALL EOS_RegisterDPressByDEps(ierr, handle, 
     .                                EOS_Ideal_Fluid_DPressByDEps);

      ELSE

         call CCTK_WARN(0, "Unable to register an EOS method from fortran 8-(")

      ENDIF

      EOS_Ideal_Fluid_Startup = 0

      end function EOS_Ideal_Fluid_Startup


