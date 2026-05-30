
#include "cctk.h"
#include "cctk_Functions.h"

c     ========================================================================

      SUBROUTINE D3_to_D2(cctkGH,conformal_state,do_ADMmass,do_momentum,do_spin,
     &     Psi_power,origin,myproc,interpolation_operator,interpolation_order,
     &     Dx,Dy,Dz,Psi,
     &     g00,gxx,gxy,gxz,gyy,gyz,gzz,hxx,hxy,hxz,hyy,hyz,hzz,
     &     x,y,z,eta,Nt,Np,nx,ny,nz,theta,phi,
     &     Psis,g00s,gxxs,gxys,gxzs,gyys,gyzs,gzzs,dPsis,dgxxs,dgxys,dgxzs,
     &     dgyys,dgyzs,dgzzs,ADMmass_int1,ADMmass_int2,
     &     momentum_int1,momentum_int2,momentum_int3,
     &     spin_int1,spin_int2,spin_int3,Extract_temp3d)

c     ------------------------------------------------------------------------
c
c     Project the 3-metric and its 1st radial derivatives onto a 2-sphere.
c
c     ------------------------------------------------------------------------

      USE met_rad_der_int
      USE ADMmass_integrand3D_int
      USE momentum_integrand3D_int
      USE spin_integrand3D_int

      IMPLICIT NONE

      DECLARE_CCTK_FUNCTIONS

c     Input variables

      CCTK_POINTER :: cctkGH

      CCTK_INT,INTENT(IN) ::
     &     conformal_state
      INTEGER,INTENT(IN) ::
     &     myproc,Psi_power
      CCTK_INT, INTENT(IN) ::
     &     Nt,Np,nx,ny,nz,do_momentum,do_spin,interpolation_order
      CCTK_REAL,INTENT(IN) ::
     &     origin(3),Dx,Dy,Dz,eta
      CCTK_REAL,INTENT(IN),DIMENSION(Nt) :: theta
      CCTK_REAL,INTENT(IN),DIMENSION(Np) :: phi
      CCTK_REAL,INTENT(IN),DIMENSION(nx) :: x
      CCTK_REAL,INTENT(IN),DIMENSION(ny) :: y
      CCTK_REAL,INTENT(IN),DIMENSION(nz) :: z
      CCTK_REAL,INTENT(IN),DIMENSION(nx,ny,nz) ::
     &     Psi,g00,gxx,gxy,gxz,gyy,gyz,gzz,
     &     hxx,hxy,hxz,hyy,hyz,hzz
      INTEGER,INTENT(IN) :: do_ADMmass(2)
      CCTK_STRING,INTENT(IN) :: 
     &     interpolation_operator
      CCTK_REAL,INTENT(INOUT),DIMENSION(nx,ny,nz) ::
     &     Extract_temp3d

c     Output variables

      CCTK_REAL,INTENT(OUT),DIMENSION(Nt,Np) ::
     &     Psis,g00s,gxxs,gxys,gxzs,gyys,
     &     gyzs,gzzs,dPsis,dgxxs,dgxys,dgxzs,dgyys,dgyzs,dgzzs,
     &     ADMmass_int1,ADMmass_int2,
     &     momentum_int1,momentum_int2,momentum_int3,
     &     spin_int1,spin_int2,spin_int3


c     Local variables, passed on

      LOGICAL ::
     &    err_flag
      INTEGER ::
     &    iorder,npoints
      INTEGER,DIMENSION(Nt,Np) ::
     &    ib,jb,kb
      CCTK_REAL,DIMENSION(Nt,Np) ::
     &    xs,ys,zs,ux,uy,uz,xb,yb,zb

c     Local variables, here only

      INTEGER :: i, j, num_arrays, vindex, ierror
      INTEGER :: interp_handle, param_table_handle, coord_system_handle
      CCTK_POINTER, dimension(3) :: interp_coords
      CCTK_INT,  dimension(8)    :: in_array_indices
      CCTK_POINTER, dimension(8) :: out_arrays
      CCTK_INT, dimension(8)     :: out_array_type_codes

      character(128) options_string
      character(128) operator
      CCTK_INT nchars
      character(len=20) :: error_code

c     ------------------------------------------------------------------------

c     Compute Cartesian coordinates on the surface
c     --------------------------------------------
c      print*,'Entering D3_to_D2'
c      print*,nx,ny,nz,Np,Nt
      DO j = 1, Np
         DO i = 1, Nt

            xs(i,j) = origin(1)+eta*SIN(theta(i))*COS(phi(j))
            ys(i,j) = origin(2)+eta*SIN(theta(i))*SIN(phi(j))
            zs(i,j) = origin(3)+eta*COS(theta(i))

         ENDDO
      ENDDO

c      print*,'Coordinates set'

c     Only do interpolation on one processor
c     --------------------------------------
      SELECT CASE (myproc)

      CASE (0)
         npoints = Np*Nt

      CASE DEFAULT
         npoints = 0

      END SELECT

c      print*,'Number of points set'

c     Get the interpolator, parameter table, and coordinate system handles
c     --------------------------------------------------------------------
      interp_handle = -1
      param_table_handle = -1
      coord_system_handle = -1


      call CCTK_FortranString(nchars,interpolation_operator,operator)
      call CCTK_InterpHandle (interp_handle,operator)
      if (interp_handle .lt. 0) then
        call CCTK_WARN(0,"Cannot get handle for interpolation ! Forgot to activate an implementation providing interpolation operators ??")
      endif

      options_string = "order = " // char(ichar('0') + interpolation_order)
      call Util_TableCreateFromString (param_table_handle, options_string)
      if (param_table_handle .lt. 0) then
        call CCTK_WARN(0,"Cannot create parameter table for interpolator")
      endif

      call CCTK_CoordSystemHandle (coord_system_handle, "cart3d")
      if (coord_system_handle .lt. 0) then
        call CCTK_WARN(0,"Cannot get handle for cart3d coordinate system ! Forgot to activate an implementation providing coordinates ??")
      endif

c      print*,'Handles obtained'

c     fill in the input/output arrays for the interpolator
c     ----------------------------------------------------
      interp_coords(1) = CCTK_PointerTo(xs)
      interp_coords(2) = CCTK_PointerTo(ys)
      interp_coords(3) = CCTK_PointerTo(zs)

      call CCTK_VarIndex(vindex, "extract::g00")
      in_array_indices(1) = vindex
      call CCTK_VarIndex(vindex, "admbase::gxx")
      in_array_indices(2) = vindex
      call CCTK_VarIndex(vindex, "admbase::gxy")
      in_array_indices(3) = vindex
      call CCTK_VarIndex(vindex, "admbase::gxz")
      in_array_indices(4) = vindex
      call CCTK_VarIndex(vindex, "admbase::gyy")
      in_array_indices(5) = vindex
      call CCTK_VarIndex(vindex, "admbase::gyz")
      in_array_indices(6) = vindex
      call CCTK_VarIndex(vindex, "admbase::gzz")
      in_array_indices(7) = vindex
      call CCTK_VarIndex(vindex, "staticconformal::psi")
      in_array_indices(8) = vindex

      out_arrays(1) = CCTK_PointerTo(g00s)
      out_arrays(2) = CCTK_PointerTo(gxxs)
      out_arrays(3) = CCTK_PointerTo(gxys)
      out_arrays(4) = CCTK_PointerTo(gxzs)
      out_arrays(5) = CCTK_PointerTo(gyys)
      out_arrays(6) = CCTK_PointerTo(gyzs)
      out_arrays(7) = CCTK_PointerTo(gzzs)
      out_arrays(8) = CCTK_PointerTo(Psis)

      out_array_type_codes = CCTK_VARIABLE_REAL


c     Project un-physical metric and conformal factor  onto sphere
c     ------------------------------------------------------------
      if (conformal_state > 0) then
        num_arrays = 8
      else
        num_arrays = 7
      end if

c      print*,'1'
c      print*,in_array_indices
c      print*,out_arrays
      call CCTK_InterpGridArrays (ierror, cctkGH, 3, interp_handle,
     $                            param_table_handle, coord_system_handle,
     $                            npoints, CCTK_VARIABLE_REAL, interp_coords,
     $                            num_arrays, in_array_indices,
     $                            num_arrays, out_array_type_codes, out_arrays)
c      print*,'1 done'
      if (ierror < 0) then
        call CCTK_WARN (1, "interpolator call returned an error code")
        write(error_code,'(a7,i13)') 'code = ',ierror
        call CCTK_WARN (1, error_code )
      endif


c     Calculate radial derivatives and project onto sphere
c     ----------------------------------------------------
      call CCTK_VarIndex(vindex, "extract::temp3d")
      in_array_indices(1) = vindex

      if (conformal_state > 0) then
c         print*,'1 1/2'
         CALL met_rad_der(origin,Dx,Dy,Dz,x,y,z,Psi,Extract_temp3d)
c         print*,'1 1/2 done'
         CALL CCTK_SyncGroup(ierror,cctkGH,"extract::temps")

         out_arrays(1) = CCTK_PointerTo(dPsis)
c         print*,'2'
c         print*,in_array_indices(1)
c         print*,out_arrays(1)
         call CCTK_InterpGridArrays (ierror, cctkGH, 3, interp_handle,
     $                               param_table_handle, coord_system_handle,
     $                               npoints, CCTK_VARIABLE_REAL, interp_coords,
     $                               1, in_array_indices(1),
     $                               1, out_array_type_codes, out_arrays(1))
c         print*,'2 done'
         if (ierror < 0) then
           call CCTK_WARN (1, "interpolator call returned an error code")
           write(error_code,'(a7,i13)') 'code = ',ierror
           call CCTK_WARN (1, error_code )
         endif
      end if

c      print*,'3'
      CALL met_rad_der(origin,Dx,Dy,Dz,x,y,z,gxx,Extract_temp3d)
c      print*,'3 done'
c      print*,'4'
      CALL CCTK_SyncGroup(ierror,cctkGH,"extract::temps")
c      print*,'4 done'

      out_arrays(1) = CCTK_PointerTo(dgxxs)
c      print*,'5'
c      print*,in_array_indices(1)
c      print*,out_arrays(1)
      call CCTK_InterpGridArrays (ierror, cctkGH, 3, interp_handle,
     $                            param_table_handle, coord_system_handle,
     $                            npoints, CCTK_VARIABLE_REAL, interp_coords,
     $                            1, in_array_indices(1),
     $                            1, out_array_type_codes, out_arrays(1))
c      print*,'5 done'
      if (ierror < 0) then
        call CCTK_WARN (1, "interpolator call returned an error code")
        write(error_code,'(a7,i13)') 'code = ',ierror
        call CCTK_WARN (1, error_code )
      endif

      CALL met_rad_der(origin,Dx,Dy,Dz,x,y,z,gxy,Extract_temp3d)
      CALL CCTK_SyncGroup(ierror,cctkGH,"extract::temps")

      out_arrays(1) = CCTK_PointerTo(dgxys)
c      print*,'6'
c      print*,in_array_indices(1)
c      print*,out_arrays(1)
      call CCTK_InterpGridArrays (ierror, cctkGH, 3, interp_handle,
     $                            param_table_handle, coord_system_handle,
     $                            npoints, CCTK_VARIABLE_REAL, interp_coords,
     $                            1, in_array_indices(1),
     $                            1, out_array_type_codes, out_arrays(1))
c      print*,'6 done'
      if (ierror < 0) then
        call CCTK_WARN (1, "interpolator call returned an error code")
        write(error_code,'(a7,i13)') 'code = ',ierror
        call CCTK_WARN (1, error_code )
      endif

      CALL met_rad_der(origin,Dx,Dy,Dz,x,y,z,gxz,Extract_temp3d)
      CALL CCTK_SyncGroup(ierror,cctkGH,"extract::temps")

      out_arrays(1) = CCTK_PointerTo(dgxzs)
c      print*,'7'
c      print*,in_array_indices(1)
c      print*,out_arrays(1)
      call CCTK_InterpGridArrays (ierror, cctkGH, 3, interp_handle,
     $                            param_table_handle, coord_system_handle,
     $                            npoints, CCTK_VARIABLE_REAL, interp_coords,
     $                            1, in_array_indices(1),
     $                            1, out_array_type_codes, out_arrays(1))
c      print*,'7 done'
      if (ierror < 0) then
        call CCTK_WARN (1, "interpolator call returned an error code")
        write(error_code,'(a7,i13)') 'code = ',ierror
        call CCTK_WARN (1, error_code )
      endif

      CALL met_rad_der(origin,Dx,Dy,Dz,x,y,z,gyy,Extract_temp3d)
      CALL CCTK_SyncGroup(ierror,cctkGH,"extract::temps")

      out_arrays(1) = CCTK_PointerTo(dgyys)
      call CCTK_InterpGridArrays (ierror, cctkGH, 3, interp_handle,
     $                            param_table_handle, coord_system_handle,
     $                            npoints, CCTK_VARIABLE_REAL, interp_coords,
     $                            1, in_array_indices(1),
     $                            1, out_array_type_codes, out_arrays(1))
      if (ierror < 0) then
        call CCTK_WARN (1, "interpolator call returned an error code")
        write(error_code,'(a7,i13)') 'code = ',ierror
        call CCTK_WARN (1, error_code )
      endif

      CALL met_rad_der(origin,Dx,Dy,Dz,x,y,z,gyz,Extract_temp3d)
      CALL CCTK_SyncGroup(ierror,cctkGH,"extract::temps")

      out_arrays(1) = CCTK_PointerTo(dgyzs)
c      print*,'8'
c      print*,in_array_indices(1)
c      print*,out_arrays(1)
      call CCTK_InterpGridArrays (ierror, cctkGH, 3, interp_handle,
     $                            param_table_handle, coord_system_handle,
     $                            npoints, CCTK_VARIABLE_REAL, interp_coords,
     $                            1, in_array_indices(1),
     $                            1, out_array_type_codes, out_arrays(1))
c      print*,'8 done'
      if (ierror < 0) then
        call CCTK_WARN (1, "interpolator call returned an error code")
        write(error_code,'(a7,i13)') 'code = ',ierror
        call CCTK_WARN (1, error_code )
      endif

      CALL met_rad_der(origin,Dx,Dy,Dz,x,y,z,gzz,Extract_temp3d)
      CALL CCTK_SyncGroup(ierror,cctkGH,"extract::temps")

      out_arrays(1) = CCTK_PointerTo(dgzzs)
c      print*,'9'
c      print*,in_array_indices(1)
c      print*,out_arrays(1)
      call CCTK_InterpGridArrays (ierror, cctkGH, 3, interp_handle,
     $                            param_table_handle, coord_system_handle,
     $                            npoints, CCTK_VARIABLE_REAL, interp_coords,
     $                            1, in_array_indices(1),
     $                            1, out_array_type_codes, out_arrays(1))
c      print*,'9 done'
      if (ierror < 0) then
        call CCTK_WARN (1, "interpolator call returned an error code")
        write(error_code,'(a7,i13)') 'code = ',ierror
        call CCTK_WARN (1, error_code )
      endif


c     Calculate integrands for ADM masses
c     -----------------------------------
c     Standard equation
      IF (do_ADMmass(1) == 1) THEN
c         print*,'10'
c         print*,in_array_indices(1)
c         print*,out_arrays(1)
         CALL ADMmass_integrand3D(origin,Dx,Dy,Dz,x,y,z,gxx,gxy,
     &        gxz,gyy,gyz,gzz,Extract_temp3d,Psi,Psi_power,conformal_state)
         CALL CCTK_SyncGroup(ierror,cctkGH,"extract::temps")
c         print*,'10 done'

         out_arrays(1) = CCTK_PointerTo(ADMmass_int1)
c         print*,'11'
c         print*,in_array_indices(1)
c         print*,out_arrays(1)
         call CCTK_InterpGridArrays (ierror, cctkGH, 3, interp_handle,
     $                               param_table_handle, coord_system_handle,
     $                               npoints, CCTK_VARIABLE_REAL, interp_coords,
     $                               1, in_array_indices(1),
     $                               1, out_array_type_codes, out_arrays(1))
c         print*,'11 done'
         if (ierror < 0) then
           call CCTK_WARN (1, "interpolator call returned an error code")
           write(error_code,'(a7,i13)') 'code = ',ierror
           call CCTK_WARN (1, error_code )
         endif
      END IF

c     Conformal equation
      IF (do_ADMmass(2) == 1) THEN
             ADMmass_int2 = -eta**2/2D0/3.1416D0*dPsis
      ENDIF

c     Calculate integrands for momentum
c     ---------------------------------
      IF (do_momentum==1) THEN

c         print*,'12'
c         print*,in_array_indices(1)
c         print*,out_arrays(1)
         CALL momentum_integrand3D(origin,Dx,Dy,Dz,x,y,z,
     &        gxx,gxy,gxz,gyy,gyz,gzz,hxx,hxy,hxz,hyy,hyz,hzz,
     &        Extract_temp3d,Psi,Psi_power,conformal_state)
         CALL CCTK_SyncGroup(ierror,cctkGH,"extract::temps")
c         print*,'12 done'

         out_arrays(1) = CCTK_PointerTo(momentum_int1)
c         print*,'13'
c         print*,in_array_indices(1)
c         print*,out_arrays(1)
         call CCTK_InterpGridArrays (ierror, cctkGH, 3, interp_handle,
     $                               param_table_handle, coord_system_handle,
     $                               npoints, CCTK_VARIABLE_REAL, interp_coords,
     $                               1, in_array_indices(1),
     $                               1, out_array_type_codes, out_arrays(1))
c         print*,'13 done'
         if (ierror < 0) then
           call CCTK_WARN (1, "interpolator call returned an error code")
           write(error_code,'(a7,i13)') 'code = ',ierror
           call CCTK_WARN (1, error_code )
         endif

c         print*,'14'
c         print*,in_array_indices(1)
c         print*,out_arrays(1)
         CALL momentum_integrand3D(origin,Dx,Dy,Dz,x,y,z,
     &        gxx,gxy,gxz,gyy,gyz,gzz,hxx,hxy,hxz,hyy,hyz,hzz,
     &        Extract_temp3d,Psi,Psi_power,conformal_state)
         CALL CCTK_SyncGroup(ierror,cctkGH,"extract::temps")
c         print*,'14 done'

         out_arrays(1) = CCTK_PointerTo(momentum_int2)
c         print*,'15'
c         print*,in_array_indices(1)
c         print*,out_arrays(1)
         call CCTK_InterpGridArrays (ierror, cctkGH, 3, interp_handle,
     $                               param_table_handle, coord_system_handle,
     $                               npoints, CCTK_VARIABLE_REAL, interp_coords,
     $                               1, in_array_indices(1),
     $                               1, out_array_type_codes, out_arrays(1))
c         print*,'15 done'
         if (ierror < 0) then
           call CCTK_WARN (1, "interpolator call returned an error code")
           write(error_code,'(a7,i13)') 'code = ',ierror
           call CCTK_WARN (1, error_code )
         endif

c         print*,'16'
c         print*,in_array_indices(1)
c         print*,out_arrays(1)
         CALL momentum_integrand3D(origin,Dx,Dy,Dz,x,y,z,
     &        gxx,gxy,gxz,gyy,gyz,gzz,hxx,hxy,hxz,hyy,hyz,hzz,
     &        Extract_temp3d,Psi,Psi_power,conformal_state)
         CALL CCTK_SyncGroup(ierror,cctkGH,"extract::temps")
c         print*,'16 done'

         out_arrays(1) = CCTK_PointerTo(momentum_int3)
c         print*,'17'
c         print*,in_array_indices(1)
c         print*,out_arrays(1)
         call CCTK_InterpGridArrays (ierror, cctkGH, 3, interp_handle,
     $                               param_table_handle, coord_system_handle,
     $                               npoints, CCTK_VARIABLE_REAL, interp_coords,
     $                               1, in_array_indices(1),
     $                               1, out_array_type_codes, out_arrays(1))
c         print*,'17 done'
         if (ierror < 0) then
           call CCTK_WARN (1, "interpolator call returned an error code")
           write(error_code,'(a7,i13)') 'code = ',ierror
           call CCTK_WARN (1, error_code )
         endif
      END IF

c     Calculate integrands for spin
c     -----------------------------
      IF (do_spin==1) THEN

c         print*,'18'
c         print*,in_array_indices(1)
c         print*,out_arrays(1)
         CALL spin_integrand3D(origin,x,y,z,
     &        gxx,gxy,gxz,gyy,gyz,gzz,hxx,hxy,hxz,hyy,hyz,hzz,
     &        Extract_temp3d,Psi,Psi_power,conformal_state)
         CALL CCTK_SyncGroup(ierror,cctkGH,"extract::temps")
c         print*,'18 done'

         out_arrays(1) = CCTK_PointerTo(spin_int1)
c         print*,'19'
c         print*,in_array_indices(1)
c         print*,out_arrays(1)
         call CCTK_InterpGridArrays (ierror, cctkGH, 3, interp_handle,
     $                               param_table_handle, coord_system_handle,
     $                               npoints, CCTK_VARIABLE_REAL, interp_coords,
     $                               1, in_array_indices(1),
     $                               1, out_array_type_codes, out_arrays(1))
c         print*,'19 done'
         if (ierror < 0) then
           call CCTK_WARN (1, "interpolator call returned an error code")
           write(error_code,'(a7,i13)') 'code = ',ierror
           call CCTK_WARN (1, error_code )
         endif

c         print*,'20'
c         print*,in_array_indices(1)
c         print*,out_arrays(1)
         CALL spin_integrand3D(origin,x,y,z,
     &        gxx,gxy,gxz,gyy,gyz,gzz,hxx,hxy,hxz,hyy,hyz,hzz,
     &        Extract_temp3d,Psi,Psi_power,conformal_state)
         CALL CCTK_SyncGroup(ierror,cctkGH,"extract::temps")
c         print*,'20 done'

         out_arrays(1) = CCTK_PointerTo(spin_int2)
c         print*,'21'
c         print*,in_array_indices(1)
c         print*,out_arrays(1)
         call CCTK_InterpGridArrays (ierror, cctkGH, 3, interp_handle,
     $                               param_table_handle, coord_system_handle,
     $                               npoints, CCTK_VARIABLE_REAL, interp_coords,
     $                               1, in_array_indices(1),
     $                               1, out_array_type_codes, out_arrays(1))
c         print*,'21 done'
         if (ierror < 0) then
           call CCTK_WARN (1, "interpolator call returned an error code")
           write(error_code,'(a7,i13)') 'code = ',ierror
           call CCTK_WARN (1, error_code )
         endif

c         print*,'22'
c         print*,in_array_indices(1)
c         print*,out_arrays(1)
         CALL spin_integrand3D(origin,x,y,z,
     &        gxx,gxy,gxz,gyy,gyz,gzz,hxx,hxy,hxz,hyy,hyz,hzz,
     &        Extract_temp3d,Psi,Psi_power,conformal_state)
         CALL CCTK_SyncGroup(ierror,cctkGH,"extract::temps")
c         print*,'22 done'

         out_arrays(1) = CCTK_PointerTo(spin_int3)
c         print*,'23'
c         print*,in_array_indices(1)
c         print*,out_arrays(1)
         call CCTK_InterpGridArrays (ierror, cctkGH, 3, interp_handle,
     $                               param_table_handle, coord_system_handle,
     $                               npoints, CCTK_VARIABLE_REAL, interp_coords,
     $                               1, in_array_indices(1),
     $                               1, out_array_type_codes, out_arrays(1))
c         print*,'23 done'
         if (ierror < 0) then
           call CCTK_WARN (1, "interpolator call returned an error code")
           write(error_code,'(a7,i13)') 'code = ',ierror
           call CCTK_WARN (1, error_code )
         endif

      END IF

      call Util_TableDestroy (ierror, param_table_handle)

      END SUBROUTINE D3_to_D2
