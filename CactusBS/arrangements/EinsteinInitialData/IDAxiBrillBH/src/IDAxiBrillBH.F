c/*@@
c  @file      IDAxiBrillBH.F
c  @date
c  @author
c  @desc
c
c  @enddesc
c@@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

c/*@@
c  @routine    IDAxiBrillBH
c  @date
c  @author
c  @desc
c
c  @enddesc
c  @calls
c  @calledby
c  @history
c
c  @endhistory
c@@*/

      subroutine IDAxiBrillBH(CCTK_ARGUMENTS)

      use cctk

      implicit none

      DECLARE_CCTK_ARGUMENTS_IDAXIBRILLBH
      DECLARE_CCTK_PARAMETERS
      DECLARE_CCTK_FUNCTIONS

      interface
      subroutine mgparm(m,ifd59,id5,id9,idi,idg,imx,jmx)
      implicit CCTK_REAL (a-h,o-z)
      end
      subroutine mg5 (idim,ilower,iupper,jdim,jlower,jupper,
     &     cc,cn,cs,cw,ce,u,rhs,
     &     id5,id9,idi,idg,ifmg,eps,rmax,ier)
      implicit CCTK_REAL (a-h,o-z)
      CCTK_REAL cc(idim,jdim),cn(idim,jdim),cs(idim,jdim),cw(idim,jdim),
     &     ce(idim,jdim),cnw(idim,jdim),cne(idim,jdim),csw(idim,jdim),
     &     cse(idim,jdim),u(idim,jdim),rhs(idim,jdim)
      end
      end interface

      CCTK_REAL, parameter :: one = 1
      CCTK_REAL, parameter :: eps = 1.0e-10
      CCTK_REAL axibheps, rmax, dq, deta
      integer levels,id5,id9,idi,idg,ier
      CCTK_REAL, allocatable :: cc(:,:),ce(:,:),cw(:,:),cn(:,:),cs(:,:),
     $     rhs(:,:),psi2d(:,:),detapsi2d(:,:),dqpsi2d(:,:),
     $     detaetapsi2d(:,:),detaqpsi2d(:,:),dqqpsi2d(:,:)
      CCTK_REAL, allocatable :: etagrd(:),qgrd(:)
      CCTK_REAL, allocatable :: eta(:,:,:),abseta(:,:,:),sign_eta(:,:,:),
     $     q(:,:,:),phi(:,:,:)
      CCTK_REAL, allocatable :: psi2dv(:,:,:),detapsi2dv(:,:,:),
     $     dqpsi2dv(:,:,:),detaetapsi2dv(:,:,:),detaqpsi2dv(:,:,:),
     $     dqqpsi2dv(:,:,:)
      CCTK_REAL error_at_this_grid_point,max_error_in_grid
      CCTK_REAL o1,o2,o3,o4,o5,o6,o7,o8,o9
      CCTK_REAL o10,o11,o12,o13,o14,o15,o16,o17,o18,o19
      CCTK_REAL o20,o21,o22,o23,o24,o25,o26,o27,o28,o29
      CCTK_REAL o30,o31,o32,o33,o34,o35,o36,o37,o38,o39
      CCTK_REAL o40,o41,o42,o43,o44,o45,o46,o47,o48,o49
      CCTK_REAL o50,o51,o52,o53,o54,o55,o56,o57,o58,o59
      CCTK_REAL o60,o61,o62,o63,o64,o65,o66,o67,o68,o69
      CCTK_REAL o70,o71,o72,o73,o74,o75,o76,o77,o78,o79
      CCTK_REAL o80,o81,o82,o83,o84,o85,o86,o87,o88,o89
      CCTK_REAL o90,o91,o92,o93,o94,o95,o96
      integer i22
      CCTK_REAL pi
      CCTK_REAL adm
      CCTK_REAL exp_mhalf_eta, psi3d
      CCTK_INT nx,ny,nz
      integer i,j,k,nquads
      integer npoints,ierror
      integer neb, nqb
      integer posn
      integer io_status

      integer, parameter :: max_string_length = 500

      integer                         param_table_handle, interp_handle
      character(max_string_length) :: message_buffer

      integer                         fstring_length
      character(max_string_length) :: interpolator_name_fstring
      character(max_string_length) :: interpolator_pars_fstring
      character(max_string_length) :: output_psi2D_file_name_fstring

      CCTK_REAL,    dimension(2) :: coord_origin, coord_delta
      CCTK_POINTER, dimension(2) :: interp_coords
      CCTK_POINTER, dimension(6) :: in_arrays, out_arrays
      CCTK_INT,     dimension(2) :: in_array_dims
      CCTK_INT,     dimension(6), parameter :: type_codes = CCTK_VARIABLE_REAL


      pi = 4.0d0*atan(1.0d0)

c     Set up the grid spacings
      nx = cctk_lsh(1)
      ny = cctk_lsh(2)
      nz = cctk_lsh(3)

c
c ***** set up interpolator handle and parameter table *****
c

c
c ... we do this now, rather than waiting till we need them,
c     so any errors (eg user forgot to activate a suitable interpolator thorn)
c     get printed right away, rather than making the user wait for them
c ... n.b. we must first convert our C-style string parameters
c     to Fortran strings
c
      call CCTK_FortranString(fstring_length,
     $                        interpolator_name,
     $                        interpolator_name_fstring)
      call CCTK_InterpHandle(interp_handle,
     $                       interpolator_name_fstring(1:fstring_length))
      if (interp_handle .lt. 0) then
         call CCTK_ERROR("Cannot get interpolator handle!  did you forgot to activate a suitable local-interpolator thorn?")
c     " (fix syntax highlighting)
      endif

c
c build the interpolator parameter table from a suitable string:
c if  interpolator_pars  is a nonempty string, use that, otherwise
c use  "order=n"  where  n  is given by the  interpolation_order
c parameter
c
      call CCTK_FortranString(fstring_length,
     $                        interpolator_pars,
     $                        interpolator_pars_fstring)
      if (fstring_length .eq. 0) then
         interpolator_pars_fstring
     $         = "order=" // char(ichar('0') + interpolation_order)
      endif
      call Util_TableCreateFromString
     $            (param_table_handle,
     $             interpolator_pars_fstring(1:fstring_length))
      if (param_table_handle .lt. 0) then
        write(message_buffer, '(A,I8)')
     $      'failed to create interpolator param table: error code ',
     $      param_table_handle
        call CCTK_WARN(CCTK_WARN_ABORT, message_buffer)
      endif

c
c ***** solve Brill-wave equation on 2D (eta,theta) grid *****
c
c The code uses the abbreviation 'q' for theta.  This is a bit confusing,
c because this is *not* the same quantity as $q(\eta,\theta)$ described
c in the thorn guide (which is is stored in the 3-D array q(nx,ny,nz)).
c
c The (eta,theta) grid spans the range
c	eta   in [0,etamax] + ghost zones   (neb points, spacing deta)
c	theta in [0,pi    ] + ghost zones   (nqb points, spacing dq  )
c
c     2D grid size NE x NQ, plus 2 zones for boundaries
c
c     21/11/00 TR: dont change parameters in place
c                  but keep a copy in local variables
c                  Otherwise the changed parameters cause trouble
c                  after recovery.
c
      neb = ne+2
      nqb = nq+2

      allocate(          cc(neb,nqb))
      allocate(          ce(neb,nqb))
      allocate(          cw(neb,nqb))
      allocate(          cn(neb,nqb))
      allocate(          cs(neb,nqb))
      allocate(         rhs(neb,nqb))
      allocate(       psi2d(neb,nqb))
      allocate(   detapsi2d(neb,nqb))
      allocate(     dqpsi2d(neb,nqb))
      allocate(detaetapsi2d(neb,nqb))
      allocate(  detaqpsi2d(neb,nqb))
      allocate(    dqqpsi2d(neb,nqb))

      allocate(etagrd(neb))
      allocate(  qgrd(nqb))

      allocate(          eta(nx,ny,nz))
      allocate(       abseta(nx,ny,nz))
      allocate(     sign_eta(nx,ny,nz))
      allocate(            q(nx,ny,nz))
      allocate(          phi(nx,ny,nz))
      allocate(       psi2dv(nx,ny,nz))
      allocate(   detapsi2dv(nx,ny,nz))
      allocate(     dqpsi2dv(nx,ny,nz))
      allocate(detaetapsi2dv(nx,ny,nz))
      allocate(  detaqpsi2dv(nx,ny,nz))
      allocate(    dqqpsi2dv(nx,ny,nz))

c Initialize some arrays
      psi2d = 1.0d0
      detapsi2d = 0.0d0

      nquads = 2
      dq = nquads*0.5d0*pi/(nqb-2)
      deta = etamax/(neb-3)

      do j=1,nqb
        qgrd(j) = (j-1.5d0)*dq
        do i=1,neb
          etagrd(i) = (i-2)*deta
#include "EinsteinInitialData/IDAxiBrillBH/src/bhbrill.x"
        enddo
      enddo
c Boundary conditions
      do j=1,nqb
        ce(2,j)=ce(2,j)+cw(2,j)
        cw(2,j)=0.0d0

        cw(neb-1,j)=cw(neb-1,j)+ce(neb-1,j)
        cc(neb-1,j)=cc(neb-1,j)-deta*ce(neb-1,j)
        ce(neb-1,j)=0.0d0

      enddo
      do i=1,neb
        cc(i,2)=cc(i,2)+cs(i,2)
        cs(i,2)=0.0d0
        cc(i,nqb-1)=cc(i,nqb-1)+cn(i,nqb-1)
        cn(i,nqb-1)=0.0d0
      enddo

c Do the solve
      axibheps = error_tolerance
      call CCTK_INFO("Calling axisymmetric solver")

      call mgparm (levels,5,id5,id9,idi,idg,neb,nqb)
      call mg5 (neb,2,neb-1,nqb,2,nqb-1,
     $          cc,cn,cs,cw,ce,psi2d,rhs,
     $          id5,id9,idi,idg,1,axibheps,rmax,ier)

      call CCTK_INFO("Solve complete")

c
c The solution is (hopefully) now available.
c
      if(ier .ne. 0) then
         write(message_buffer, '(A,I8)')
     $      'failed to solve elliptic equation: ier=', ier
         call CCTK_WARN(CCTK_WARN_ABORT, message_buffer)
      end if

      print *,'rmax = ',rmax
      print *,'axibheps = ',axibheps
      print *,'psi2d = ',maxval(psi2d),' ',minval(psi2d)

      max_error_in_grid = 0.0d0
      do j=2,nqb-1
        do i=2,neb-1
          error_at_this_grid_point = rhs(i,j)
     $                               - psi2d(i,j  )*cc(i,j)
     $                               - psi2d(i,j+1)*cn(i,j)
     $                               - psi2d(i,j-1)*cs(i,j)
     $                               - psi2d(i+1,j)*ce(i,j)
     $                               - psi2d(i-1,j)*cw(i,j)
          max_error_in_grid = max(max_error_in_grid,
     $                            abs(error_at_this_grid_point))
        enddo
      enddo
      print *,'Resulting eps =',max_error_in_grid

c Boundary conditions
      do j=1,nqb
         psi2d(1  ,j) = psi2d(3,j)
         psi2d(neb,j) = - deta*psi2d(neb-1,j) + psi2d(neb-2,j)
      enddo

      do i=1,neb
        psi2d(i,1  ) = psi2d(i,2)
        psi2d(i,nqb) = psi2d(i,nqb-1)
      enddo

c derivatives of psi
      do j=2,nqb-1
         do i=2,neb-1
            dqpsi2d  (i,j) = 0.5d0*(psi2d(i,j+1)-psi2d(i,j-1))/dq
            dqqpsi2d (i,j) = (psi2d(i,j+1)+psi2d(i,j-1)-2.0d0*psi2d(i,j))/dq**2
            detapsi2d(i,j) = sinh(0.5d0*etagrd(i))
     $                       + 0.5d0*(psi2d(i+1,j)-psi2d(i-1,j))/deta
            detaetapsi2d(i,j)
     $          = 0.5d0*cosh(0.5d0*etagrd(i))
     $            + (psi2d(i+1,j)+psi2d(i-1,j)-2.0d0*psi2d(i,j))/deta**2
         enddo
      enddo

      do j=1,nqb
        detapsi2d(1,j)=-detapsi2d(3,j)
        detapsi2d(neb,j)=detapsi2d(neb-2,j) ! simplified

        detaetapsi2d(1,j)=detaetapsi2d(3,j)
        detaetapsi2d(neb,j)=detaetapsi2d(neb-2,j) ! simplified...

        dqqpsi2d(1,j)=dqqpsi2d(3,j)
        dqqpsi2d(neb,j)=dqqpsi2d(neb-2,j) ! simplified

        dqpsi2d(1,j)=dqpsi2d(3,j)
        dqpsi2d(neb,j)=-dq*dqpsi2d(neb-1,j)+dqpsi2d(neb-2,j)
      enddo

      do i=1,neb
        detapsi2d(i,1)=detapsi2d(i,2)
        detapsi2d(i,nqb)=detapsi2d(i,nqb-1)

        detaetapsi2d(i,1)=detaetapsi2d(i,2)
        detaetapsi2d(i,nqb)=detaetapsi2d(i,nqb-1)

        dqqpsi2d(i,1)=dqqpsi2d(i,2)
        dqqpsi2d(i,nqb)=dqqpsi2d(i,nqb-1)

        dqpsi2d(i,1)=-dqpsi2d(i,2)
        dqpsi2d(i,nqb)=-dqpsi2d(i,nqb-1)
      enddo

      do j=2,nqb-1
        do i=2,neb-1
          detaqpsi2d(i,j)=0.5d0*(detapsi2d(i,j+1)-detapsi2d(i,j-1))/dq
        enddo
      enddo

      do j=1,nqb
        detaqpsi2d(1,j)=-detaqpsi2d(3,j)
        detaqpsi2d(neb,j)=detaqpsi2d(neb-2,j) ! simplified
      enddo

      do i=1,ne
        detaqpsi2d(i,1)=-detaqpsi2d(i,2)
        detaqpsi2d(i,nqb)=-detaqpsi2d(i,nqb-1)
      enddo

      do j=1,nqb
        psi2d(:,j)=psi2d(:,j)+2.0d0*cosh(0.5d0*etagrd)
      enddo

      if (debug .ge. 6) then
         print *, '### 2-D grid results (= inputs to interpolation) ###'
         print *, 'effective 2-D grid size: neb,nqb =', neb,nqb
         print *, 'debug_{ii,jj} =', debug_ii,debug_jj
         print *, 'at this 2-D grid point...'
         print *, '             eta =', etagrd(debug_ii)
         print *, '       theta [q] =',   qgrd(debug_jj)
         print *, '           psi2d =',        psi2d(debug_ii,debug_jj)
         print *, '       detapsi2d =',    detapsi2d(debug_ii,debug_jj)
         print *, '         dqpsi2d =',      dqpsi2d(debug_ii,debug_jj)
         print *, '    detaetapsi2d =', detaetapsi2d(debug_ii,debug_jj)
         print *, '      detaqpsi2d =',   detaqpsi2d(debug_ii,debug_jj)
         print *, '        dqqpsi2d =',     dqqpsi2d(debug_ii,debug_jj)
      endif

c
c write conformal factor psi on 2D grid to output file if requested
c
      if (output_psi2D .ne. 0) then
         write (message_buffer, '(A,A,A)')
     $         'writing 2D psi to "',
     $            output_psi2D_file_name_fstring(1:fstring_length),
     $         '"'
         call CCTK_INFO(message_buffer)

         call CCTK_FortranString(fstring_length,
     $                           output_psi2D_file_name,
     $                           output_psi2D_file_name_fstring)
         open (9, iostat=io_status, status='replace',
     $         file=output_psi2D_file_name_fstring)
         if (io_status .ne. 0) then 
            write (message_buffer, '(A,A,A,I8)')
     $            'error opening psi2D output file "',
     $            output_psi2D_file_name_fstring(1:fstring_length),
     $            '": io_status=', io_status
            call CCTK_WARN(CCTK_WARN_ABORT, message_buffer)
         endif

         write (9, '(a)')
     $         '# eta          theta (=q)     psi (2D)            psi (3D)'
            do i = 1,neb
            exp_mhalf_eta = exp(-0.5d0 * etagrd(i))
               do j = 1,nqb
               psi3d = psi2d(i,j) * exp_mhalf_eta
               write (9, '(f12.8,a, f12.8,a, g20.10e3,a, g20.10e3)')
     $               etagrd(i), ' ', qgrd(j), ' ',
     $               psi2d(i,j), ' ', psi3d
               end do
            write (9, '(a)') ''
            end do

         close (9, iostat=io_status)
         if (io_status .ne. 0) then 
            write(message_buffer, '(A,I8)')
     $         'error closing psi2D output file: io_status=', io_status
            call CCTK_WARN(CCTK_WARN_ABORT, message_buffer)
         endif
      endif



c
c ***** interpolate psi and its derivatives to the xyz grid positions *****
c ***** and compute the ADM variables there *****
c
c More precisely, at this point we have q, psi, and the psi derivatives
c on the (eta,theta) grid.  We want to interpolate these values to the
c (eta,theta) locations of each of the (x,y,z) grid points.
c
c The (eta,theta) grid only spans the range eta >= 0, so we actually
c interpolate using the (x,y,z) grid points |eta| values, and fix
c things up afterwords.
c

      call CCTK_INFO("interpolating solution to xyz grid points")

      eta = log (r + eps)
      abseta = abs (eta)
      sign_eta = sign (one, eta)
      q = atan2 (sqrt ((x + eps)**2 + y**2), z)
      phi = atan2 (y, x + eps)

      if (debug .ge. 6) then
c posn = (0-origin) 1-D position of this point in interpolator arrays
c        (this is useful because interpolator error messages refer to it)
	 posn = (debug_i-1) + nx*(debug_j-1) + nx*ny*(debug_k-1)
	 print *, '### 3-D interpolation coordinates and inputs ###'
	 print *, '3-D grid size: n[xyz] =', nx,ny,nz
	 print *, 'debug_[ijk] =', debug_i,debug_j,debug_k, '==> posn =', posn
         print *, 'at this 3-D grid point, ...'
	 print *, '        x =',      x(debug_i,debug_j,debug_k)
	 print *, '        y =',      y(debug_i,debug_j,debug_k)
	 print *, '        z =',      z(debug_i,debug_j,debug_k)
	 print *, '      eta =',    eta(debug_i,debug_j,debug_k)
	 print *, '   abseta =', abseta(debug_i,debug_j,debug_k)
	 print *, 'theta [q] =',      q(debug_i,debug_j,debug_k)
	 print *, '      phi =',    phi(debug_i,debug_j,debug_k)
      endif

c set up the interpolator array pointers
      npoints = nx*ny*nz

      coord_origin(1) = etagrd(1)
      coord_origin(2) =   qgrd(1)
      coord_delta(1)  = etagrd(2) - etagrd(1)
      coord_delta(2)  =   qgrd(2) -   qgrd(1)

      if (debug .ge. 6) then
         print *, '### 2-D grid origin:   eta=', coord_origin(1)
         print *, '                     theta=', coord_origin(2)
         print *, '              delta:   eta=', coord_delta(1)
         print *, '                     theta=', coord_delta(2)
      end if

      interp_coords(1) = CCTK_PointerTo(abseta)
      interp_coords(2) = CCTK_PointerTo(q)

      in_array_dims(1) = neb
      in_array_dims(2) = nqb

      in_arrays(1) = CCTK_PointerTo(       psi2d)
      in_arrays(2) = CCTK_PointerTo(   detapsi2d)
      in_arrays(3) = CCTK_PointerTo(     dqpsi2d)
      in_arrays(4) = CCTK_PointerTo(detaetapsi2d)
      in_arrays(5) = CCTK_PointerTo(  detaqpsi2d)
      in_arrays(6) = CCTK_PointerTo(    dqqpsi2d)

      out_arrays(1) = CCTK_PointerTo(       psi2dv)
      out_arrays(2) = CCTK_PointerTo(   detapsi2dv)
      out_arrays(3) = CCTK_PointerTo(     dqpsi2dv)
      out_arrays(4) = CCTK_PointerTo(detaetapsi2dv)
      out_arrays(5) = CCTK_PointerTo(  detaqpsi2dv)
      out_arrays(6) = CCTK_PointerTo(    dqqpsi2dv)

             psi2dv = 1
         detapsi2dv = 0
           dqpsi2dv = 0
      detaetapsi2dv = 0
        detaqpsi2dv = 0
          dqqpsi2dv = 0

      call CCTK_InterpLocalUniform (ierror, 2,
     $                              interp_handle, param_table_handle,
     $                              coord_origin, coord_delta,
     $                              npoints, int(type_codes(1)), interp_coords,
     $                              6, in_array_dims, type_codes, in_arrays,
     $                              6, type_codes, out_arrays)
      if (ierror < 0) then
        write(message_buffer, '(A,I8)')
     $      'error in interpolator: ierror=', ierror
        call CCTK_WARN(CCTK_WARN_ALERT, message_buffer)
      endif

      call Util_TableDestroy (ierror, param_table_handle)

      if (debug .ge. 6) then
         print *, '### interpolation results (at this 3-D grid point) ###'
	 print *, '   psi2dv =', psi2dv(debug_i,debug_j,debug_k)
      end if
c
c ***** compute the ADMBase conformal factor, its derivatives, *****
c ***** metric, and extrinsic curvature from the interpolation output *****
c
      psi = psi2dv * exp (-0.5d0 * eta)
      detapsi2dv = sign_eta * detapsi2dv
      detaqpsi2dv = sign_eta * detaqpsi2dv

      do k=1,nz
         do j=1,ny
            do i=1,nx
c     psix = \partial psi / \partial x / psi
#include "EinsteinInitialData/IDAxiBrillBH/src/psi_1st_deriv.x"

c     psixx = \partial^2\psi / \partial x^2 / psi
#include "EinsteinInitialData/IDAxiBrillBH/src/psi_2nd_deriv.x"
            enddo
         enddo
      enddo

      do k=1,nz
         do j=1,ny
            do i=1,nx
c     Conformal metric
c       gxx = ...
c     Derivatives of the metric (currently commented-out)
c       dxgxx = 1/2 \partial gxx / \partial x
#include "EinsteinInitialData/IDAxiBrillBH/src/gij.x"

               if (r(i,j,k) < 1.0d-10) then
                  
                  psi(i,j,k) = 1
                  
                  psix(i,j,k) = 0
                  psiy(i,j,k) = 0
                  psiz(i,j,k) = 0
                  
                  psixx(i,j,k) = 0
                  psixy(i,j,k) = 0
                  psixz(i,j,k) = 0
                  psiyy(i,j,k) = 0
                  psiyz(i,j,k) = 0
                  psizz(i,j,k) = 0
                  
                  gxx(i,j,k) = 1
                  gxy(i,j,k) = 0
                  gxz(i,j,k) = 0
                  gyy(i,j,k) = 1
                  gyz(i,j,k) = 0
                  gzz(i,j,k) = 1
                  
               end if

            enddo
         enddo
      enddo

c convert to physical metric if StaticConformal is not wanted
      if (generate_StaticConformal_metric .eq. 0) then
         call CCTK_INFO("converting to physical metric")
         call ConfToPhysInPlace(nx, ny, nz,
     $				psi,
     $				gxx, gxy, gxz,
     $				     gyy, gyz,
     $				          gzz)
c        record that we now have a physical metric
         conformal_state = 0
      else
c        record that we computed psi and its 1st and 2nd derivatives
         conformal_state = 3
      end if

c Extrinsic Curvature is identically zero
      kxx = 0.0d0
      kxy = 0.0d0
      kxz = 0.0d0
      kyy = 0.0d0
      kyz = 0.0d0
      kzz = 0.0d0

      if (debug .ge. 6) then
         print *, '### final results (again at this 3-D grid point) ###'
         if (conformal_state .gt. 0) then
            print *, '### ... conformal metric with'
	    print *, '   psi =', psi(debug_i,debug_j,debug_k)
         else
            print *, '### ... physical metric with'
         end if
	 print *, '   gxx =', gxx(debug_i,debug_j,debug_k)
	 print *, '   gxy =', gxy(debug_i,debug_j,debug_k)
	 print *, '   gxz =', gxz(debug_i,debug_j,debug_k)
	 print *, '   gyy =', gyy(debug_i,debug_j,debug_k)
	 print *, '   gyz =', gyz(debug_i,debug_j,debug_k)
	 print *, '   gzz =', gzz(debug_i,debug_j,debug_k)
      endif

c
c ***** diagnostics and final cleanup
c

c ADM mass
      call CCTK_INFO("computing ADM mass")
      i = neb-15
      adm = 0.0d0
      do j=2,nqb-1
        adm = adm
     $        + (psi2d(i,j)-(psi2d(i+1,j)-psi2d(i-1,j))/deta)
     $          *exp(0.5d0*etagrd(i))
      enddo
      adm=adm/(nqb-2)
      print *,'ADM mass: ',adm

      if (CCTK_EQUALS(initial_lapse,"schwarz")) then
         write (*,*)"Initial with schwarzschild-like lapse"
         write (*,*)"using alp = (2r - adm)/(2r+adm)"
         alp = (2.0d0*r - adm)/(2.0d0*r+adm)
      endif

      deallocate(cc,ce,cw,cn,cs,rhs,psi2d,detapsi2d,dqpsi2d,
     $     detaetapsi2d,detaqpsi2d,dqqpsi2d,
     $     etagrd,qgrd,
     $     eta,abseta,sign_eta,q,phi,psi2dv,detapsi2dv,dqpsi2dv,
     $     detaetapsi2dv,detaqpsi2dv,dqqpsi2dv)

      return
      end
