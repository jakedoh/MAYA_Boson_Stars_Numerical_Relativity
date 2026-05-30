! vim: syntax=fortran
#include "cctk.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

module NullNews_Bondi
  implicit none
contains

  subroutine NullNews_News2BondiNews(cctkGH, lsh, dt, qsize, Uo, Un, Uyo, Uyn,&
       zEvolo, zEvoln, zEvolh, News, NewsB, uBondio,&
       uBondi, uBondin, redshiftB, omega, Jh, eth_J,&
       beth_J, eth_U, beth_U, du_J, beta, patch_index,&
       deltao, deltan, nzeta, n_tmp1_cgfn, n_tmp1_cgfs,&
       n_tmp2_cgfn, n_tmp2_cgfs, n_tmp3_cgfn,&
       n_tmp3_cgfs, n_tmp1_rgfn, n_tmp1_rgfs,&
       n_tmp2_rgfn, n_tmp2_rgfs, starting, stereo_mask)              
    implicit none

    CCTK_POINTER,                intent(in) :: cctkGH
    CCTK_INT,     dimension(2),  intent(in) :: lsh
    CCTK_REAL,                   intent(in) :: dt, qsize

    CCTK_INT, dimension(lsh(1),lsh(2)), intent(in) :: stereo_mask

    CCTK_COMPLEX, dimension (lsh(1), lsh(2),2), intent (in)    :: Uo, Un, Jh, eth_J, beth_J
    CCTK_COMPLEX, dimension (lsh(1), lsh(2),2), intent (in)    :: eth_U, beth_u, du_J
    CCTK_COMPLEX, dimension (lsh(1), lsh(2),2), intent (inout) :: Uyo, Uyn
    CCTK_COMPLEX, dimension (lsh(1), lsh(2),2), intent (inout) :: zEvolo, zEvoln, zEvolh, nzeta
    CCTK_COMPLEX, dimension (lsh(1), lsh(2),2), intent (in)    :: News
    CCTK_COMPLEX, dimension (lsh(1), lsh(2),2), intent (out)   :: NewsB
    CCTK_REAL,    dimension (lsh(1), lsh(2),2), intent (in)    :: omega, beta
    CCTK_REAL,    dimension (lsh(1), lsh(2),2), intent (inout) :: uBondi, uBondin, uBondio
    CCTK_REAL,    dimension (lsh(1), lsh(2),2), intent (inout) :: redshiftB
    CCTK_REAL,    dimension (lsh(1), lsh(2),2), intent (inout) :: deltao, deltan
    CCTK_INT,     dimension (lsh(1), lsh(2),2), intent (inout) :: patch_index

    CCTK_COMPLEX, dimension(lsh(1), lsh(2)), intent(inout) ::&
         n_tmp1_cgfn, n_tmp1_cgfs,&
         n_tmp2_cgfn, n_tmp2_cgfs,&
         n_tmp3_cgfn, n_tmp3_cgfs

    CCTK_REAL,    dimension(lsh(1), lsh(2)), intent(inout) ::&
         n_tmp1_rgfn, n_tmp1_rgfs,&
         n_tmp2_rgfn, n_tmp2_rgfs 

    CCTK_INT, intent(inout) :: starting

    CCTK_COMPLEX, dimension (:,:,:), allocatable, save :: Uh, Uyh, Po, Pm
    CCTK_COMPLEX, dimension (:,:,:), allocatable, save :: delta_RHS, delta_RHS_B
    CCTK_REAL,    dimension (:,:,:), allocatable, save :: delta, K, betaB, omegaB

    CCTK_REAL, save :: pii = 3.1415926535897932385d0
    CCTK_INT ::  i, j, ip, ipo, ipm
    CCTK_REAL qh, ph
    CCTK_COMPLEX :: ii  = (0.d0, 1.d0)
    CCTK_INT :: p_this, p_other, l, ierror
    logical, save :: FirstTime = .true.
    ! interpolator declarations
    integer :: N_in_arrays, N_out_arrays, N_interp_points, N_dims
    integer, save :: interp_handle, coord_system_handle, param_table_handle
    CCTK_POINTER, dimension(2) :: interp_coords
    integer                     :: interp_coords_type_code
    CCTK_INT,     dimension(5) :: in_array_indices
    CCTK_POINTER, dimension(5) :: out_arrays
    CCTK_INT,     dimension(5) :: out_array_type_codes
    integer,      dimension(10), save :: gf_indices


    CCTK_COMPLEX, dimension(:), allocatable, save ::&
         Iarr1_this, Iarr1_other,&
         Iarr2_this, Iarr2_other,&
         Iarr4_this, Iarr4_other

    CCTK_REAL, dimension(:), allocatable, save ::&
         Rarr3_this, Rarr3_other,&
         Rarr5_this, Rarr5_other,&
         qcoord_this, qcoord_other,&
         pcoord_this, pcoord_other

    CCTK_INT, dimension(:), allocatable, save ::&
         reindexi_this, reindexi_other,&
         reindexj_this, reindexj_other


    character(len=500) :: message

    DECLARE_CCTK_FUNCTIONS
    DECLARE_CCTK_PARAMETERS

    ! end declarations

    !   write(message, *) "FirstTime is ", FirstTime, " and starting is ", starting
    !   call CCTK_INFO(trim(message))

    if (FirstTime) then
       allocate (Uh(lsh(1), lsh(2),2), Uyh(lsh(1), lsh(2),2),&
            Po(lsh(1), lsh(2),2), Pm(lsh(1), lsh(2),2),&
            delta_RHS(lsh(1), lsh(2),2), delta_RHS_B(lsh(1), lsh(2),2),&
            delta(lsh(1), lsh(2),2), K(lsh(1), lsh(2),2),&
            betaB(lsh(1), lsh(2),2), omegaB(lsh(1), lsh(2),2) )

       allocate(Iarr1_this(lsh(1)*lsh(2)), Iarr1_other(lsh(1)*lsh(2)),&
            Iarr2_this(lsh(1)*lsh(2)), Iarr2_other(lsh(1)*lsh(2)),&
            Iarr4_this(lsh(1)*lsh(2)), Iarr4_other(lsh(1)*lsh(2)),&
            Rarr3_this(lsh(1)*lsh(2)), Rarr3_other(lsh(1)*lsh(2)),&
            Rarr5_this(lsh(1)*lsh(2)), Rarr5_other(lsh(1)*lsh(2)),&
            qcoord_this(lsh(1)*lsh(2)), qcoord_other(lsh(1)*lsh(2)),&
            pcoord_this(lsh(1)*lsh(2)), pcoord_other(lsh(1)*lsh(2)),&
            reindexi_this(lsh(1)*lsh(2)), reindexi_other(lsh(1)*lsh(2)),&
            reindexj_this(lsh(1)*lsh(2)), reindexj_other(lsh(1)*lsh(2))  )

       Uh=0; Uyh=0; Po=0; Pm=0; delta_RHS=0; delta_RHS_B=0
       delta=0; K=0; betaB=0; omegaB=0; Iarr1_this=0; Iarr1_other=0
       Iarr2_this=0; Iarr2_other=0; Iarr4_this=0; Iarr4_other=0
       Rarr3_this=0; Rarr3_other=0
       Rarr5_this=0; Rarr5_other=0; qcoord_this=0; qcoord_other=0
       pcoord_this=0; pcoord_other=0; reindexi_this=0; reindexi_other=0
       reindexj_this=0; reindexj_other=0

       FirstTime = .false.
       if (starting .eq. 1) then
          starting = 0
          Uyo = Uo
          ! correctly mask all point initialy outside equator
          do ip = 1, 2
             do j = 1, lsh(2)
                do i = 1, lsh(1)
                   if (dble(zEvolo(i,j,ip)*conjg(zEvolo(i,j,ip)))>1.0d0) then
                      patch_index(i,j,ip) = 1
                      deltao(i,j,ip) = pii -2.0d0 *atan2(dimag(zEvolo(i,j,ip)), dble(zEvolo(i,j,ip)))
                      zEvolo(i,j,ip) = 1.0d0 / zEvolo(i,j,ip)
                      Uyo(i,j,ip) = Uyo(i,j,ip) * (-zEvolo(i,j,ip) /conjg(zEvolo(i,j,ip)))
                   end if
                end do
             end do
          end do
       end if


       call CCTK_VarIndex(gf_indices(1),  "NullNews::n_tmp1_cgfn")
       call CCTK_VarIndex(gf_indices(2),  "NullNews::n_tmp1_cgfs")
       call CCTK_VarIndex(gf_indices(3),  "NullNews::n_tmp2_cgfn")
       call CCTK_VarIndex(gf_indices(4),  "NullNews::n_tmp2_cgfs")
       call CCTK_VarIndex(gf_indices(5),  "NullNews::n_tmp1_rgfn")
       call CCTK_VarIndex(gf_indices(6),  "NullNews::n_tmp1_rgfs")
       call CCTK_VarIndex(gf_indices(7),  "NullNews::n_tmp3_cgfn")
       call CCTK_VarIndex(gf_indices(8),  "NullNews::n_tmp3_cgfs")
       call CCTK_VarIndex(gf_indices(9),  "NullNews::n_tmp2_rgfn")
       call CCTK_VarIndex(gf_indices(10), "NullNews::n_tmp2_rgfs")

       !      interp_coords_type_code = CCTK_VARIABLE_REAL
       !      out_array_type_codes(1) = CCTK_VARIABLE_COMPLEX
       !      out_array_type_codes(2) = CCTK_VARIABLE_COMPLEX
       !      out_array_type_codes(3) = CCTK_VARIABLE_REAL
       !      out_array_type_codes(4) = CCTK_VARIABLE_COMPLEX
       !      out_array_type_codes(5) = CCTK_VARIABLE_REAL

       interp_handle =-1
       param_table_handle = -1
       coord_system_handle = -1

       if(minval(gf_indices) < 0 ) then
          call CCTK_WARN(0, "Error Obtaining VarIndex")
       end if
       call CCTK_InterpHandle (interp_handle, "generalized polynomial interpolation")
       if (interp_handle < 0) then
          call CCTK_WARN(0,"Interpolation operator not found")
       end if

       call CCTK_CoordSystemHandle (coord_system_handle, "stereo")
       if (coord_system_handle < 0) then
          call CCTK_WARN(0,"Coordinate system 'stereo' not registered")
       end if

       call  Util_TableCreateFromString(param_table_handle, &
             "order = " // char(ichar('0') + news_interp_order))
       if (param_table_handle < 0) then
          call  CCTK_WARN(-1, "Can t create parameter table!")
       end if

    end if ! FirstTime

    ! "old" level of U(y^A) is known, step forward half a dt


    Po = 1.0d0+ zEvolo*conjg(zEvolo)
    zEvoln = zEvolo + 0.25d0 * Po * Uyo * dt
    Uh = 0.5d0 * (Un + Uo)
    Pm = 1.0d0+ zEvoln*conjg(zEvoln)
    ! interpolation of Uh (at mid level) to get Uyh. 
    !   write(message, *) "maxval zbo = ", maxval(abs(zEvolo))
    !   CALL CCTK_INFO(trim(message))
    !   write(message, *) "maxval zbn = ", maxval(abs(zEvoln))
    !   CALL CCTK_INFO(trim(message))
    !   write(message, *) "maxval Uyo = ", maxval(abs(Uyo))
    !   CALL CCTK_INFO(trim(message))

    if (linearized_inertial_frame .eq. 0) then
       Uyh = 0

       do ip = 1, 2

          p_this = 0
          p_other = 0
          do j = 1, lsh(2)
             do i = 1, lsh(1)
                qh = dble(zEvoln(i,j,ip))
                ph = dimag(zEvoln(i,j,ip))
                ipm = mod( ip + patch_index(i,j,ip) - 1 ,2 ) +1
                if (abs(qh) > 1.2 .OR. abs(ph) > 1.2 ) then
                   CALL CCTK_WARN(0,"Abs qh or abs ph > 1.2")
                endif
              if (stereo_mask(i,j).eq.1) then
                if ( ipm .eq. ip ) then
                   p_this = p_this + 1
                   qcoord_this(p_this) = qh
                   pcoord_this(p_this) = ph
                   reindexi_this(p_this) = i
                   reindexj_this(p_this) = j
                else
                   p_other = p_other + 1
                   qcoord_other(p_other) = qh
                   pcoord_other(p_other) = ph
                   reindexi_other(p_other) = i
                   reindexj_other(p_other) = j
                end if
              end if
             end do
          end do
          ipo = mod( ip, 2 ) +1
          n_tmp1_cgfn = Uh(1:lsh(1), 1:lsh(2), ip)  ! not really n and s
          n_tmp1_cgfs = Uh(1:lsh(1), 1:lsh(2), ipo)

          !call CCTK_INFO("About to Interp Uh this")
          interp_coords(1)        = CCTK_PointerTo(qcoord_this)
          interp_coords(2)        = CCTK_PointerTo(pcoord_this)
          in_array_indices(1)     = gf_indices(1)
          out_arrays(1)           = CCTK_PointerTo(Iarr1_this)
          interp_coords_type_code = CCTK_VARIABLE_REAL
          out_array_type_codes(1) = CCTK_VARIABLE_COMPLEX

          N_dims=2
          N_in_arrays=1
          N_out_arrays=1
          N_interp_points = p_this

          call CCTK_InterpGridArrays (ierror, cctkGH, N_dims,&
               interp_handle, param_table_handle, coord_system_handle,&
               N_interp_points, interp_coords_type_code,&
               interp_coords,&
               N_in_arrays, in_array_indices(1:1),&
               N_out_arrays, out_array_type_codes(1:1), out_arrays(1:1))

          !      call CCTK_InterpGV(ierror, cctkGH, interp_handle,&
          !         coord_system_handle, p_this , 1, 1,&
          !        qcoord_this, pcoord_this,&
          !        CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL,&
          !        gf_indices(1), Iarr1_this,&
          !        CCTK_VARIABLE_COMPLEX)

          if (ierror .ne. 0 ) then
             call CCTK_WARN(0, "News Interp Error")
          endif

          !call CCTK_INFO("About to Interp Uh other")
          interp_coords(1)        = CCTK_PointerTo(qcoord_other)
          interp_coords(2)        = CCTK_PointerTo(pcoord_other)
          in_array_indices(2)     = gf_indices(2)
          out_arrays(2)           = CCTK_PointerTo(Iarr1_other)
          interp_coords_type_code = CCTK_VARIABLE_REAL
          out_array_type_codes(2) = CCTK_VARIABLE_COMPLEX

          N_dims=2
          N_in_arrays=1
          N_out_arrays=1
          N_interp_points = p_other

          call CCTK_InterpGridArrays (ierror, cctkGH, N_dims,&
               interp_handle, param_table_handle, coord_system_handle,&
               N_interp_points, interp_coords_type_code,&
               interp_coords,&
               N_in_arrays, in_array_indices(2:2),&
               N_out_arrays, out_array_type_codes(2:2), out_arrays(2:2))

          if (ierror .ne. 0 ) then
             call CCTK_WARN(0, "News Interp Error")
          endif

          do l = 1, p_this
             Uyh(reindexi_this(l), reindexj_this(l), ip) = Iarr1_this(l)
          end do
          do l = 1, p_other
             Uyh(reindexi_other(l), reindexj_other(l), ip) = Iarr1_other(l)
          end do

       end do ! ip

       ! now step forward a full dt, with the values of Uyh (at mid level)
       ! to get a second-order accurate (RK2) value for z(Bondi) at the new level

       zEvoln = zEvolo + .5*Pm * Uyh * dt
       zEvolh = 0.5d0 * (zEvoln + zEvolo)

       ! interpolation of Un (at the new level) to get Uyn. these values 
       ! will be needed the next time around.
       ! shamelessly reuse the scalar variables [qh, ph, qho, pho, ipo]
       ! to mean the coordinates and index at the new level in the loops below.

       ! the uBn, uBo live on the computational grid at scri, and *not* in the
       ! inertial observers grid, so we interpolate to get the change in Bondi
       ! time for those observers.


       Uyn = 0

       do ip = 1, 2

          p_this = 0
          p_other = 0
          do j = 1, lsh(2)
             do i = 1, lsh(1)
                qh = dble(zEvoln(i,j,ip))
                ph = dimag(zEvoln(i,j,ip))
                if (abs(qh) > 1.2 .OR. abs(ph) > 1.2 ) then
                   write(message, *) "CFL_VIOALTION_NEWS: ", &
                        zEvoln(i,j,ip), zEvolo(i,j,ip), Uyh(i,j,ip) * dt
                   CALL CCTK_WARN(0,trim(message))
                endif
                ipm = mod( ip + patch_index(i,j,ip) - 1 ,2 ) +1
              if (stereo_mask(i,j).eq.1) then
                if ( ipm .eq. ip ) then
                   p_this = p_this + 1
                   qcoord_this(p_this) = qh
                   pcoord_this(p_this) = ph
                   reindexi_this(p_this) = i
                   reindexj_this(p_this) = j
                else
                   p_other = p_other + 1
                   qcoord_other(p_other) = qh
                   pcoord_other(p_other) = ph
                   reindexi_other(p_other) = i
                   reindexj_other(p_other) = j
                end if
              end if
             end do
          end do
          ipo = mod( ip, 2 ) +1
          n_tmp1_cgfn = Un(1:lsh(1), 1:lsh(2), ip)  ! not really n and s
          n_tmp1_cgfs = Un(1:lsh(1), 1:lsh(2), ipo)

          !call CCTK_INFO("About to Interp Un this")
          interp_coords(1)        = CCTK_PointerTo(qcoord_this)
          interp_coords(2)        = CCTK_PointerTo(pcoord_this)
          in_array_indices(1)     = gf_indices(1)
          out_arrays(1)           = CCTK_PointerTo(Iarr1_this)
          interp_coords_type_code = CCTK_VARIABLE_REAL
          out_array_type_codes(1) = CCTK_VARIABLE_COMPLEX

          N_dims=2
          N_in_arrays=1
          N_out_arrays=1
          N_interp_points = p_this

          call CCTK_InterpGridArrays (ierror, cctkGH, N_dims,&
               interp_handle, param_table_handle, coord_system_handle,&
               N_interp_points, interp_coords_type_code,&
               interp_coords,&
               N_in_arrays, in_array_indices(1:1),&
               N_out_arrays, out_array_type_codes(1:1), out_arrays(1:1))


          if (ierror .ne. 0 ) then
             call CCTK_WARN(0, "News Interp Error")
          endif

          !call CCTK_INFO("About to Interp Un other")
          interp_coords(1)        = CCTK_PointerTo(qcoord_other)
          interp_coords(2)        = CCTK_PointerTo(pcoord_other)
          in_array_indices(2)     = gf_indices(2)
          out_arrays(2)           = CCTK_PointerTo(Iarr1_other)
          interp_coords_type_code = CCTK_VARIABLE_REAL
          out_array_type_codes(2) = CCTK_VARIABLE_COMPLEX

          N_dims=2
          N_in_arrays=1
          N_out_arrays=1
          N_interp_points = p_other

          call CCTK_InterpGridArrays (ierror, cctkGH, N_dims,&
               interp_handle, param_table_handle, coord_system_handle,&
               N_interp_points, interp_coords_type_code,&
               interp_coords,&
               N_in_arrays, in_array_indices(2:2),&
               N_out_arrays, out_array_type_codes(2:2), out_arrays(2:2))

          if (ierror .ne. 0 ) then
             call CCTK_WARN(0, "News Interp Error")
          endif

          do l = 1, p_this
             Uyn(reindexi_this(l), reindexj_this(l), ip) = Iarr1_this(l)
          end do
          do l = 1, p_other
             Uyn(reindexi_other(l), reindexj_other(l), ip) = Iarr1_other(l)
          end do

       end do ! ip




       K = dsqrt(1.0d0 + dble(Jh*conjg(Jh)))
       delta_RHS =  conjg(du_J)*Jh / (K + 1.0d0) + &
            .5d0 /(K + 1.0d0) * (Jh*(Uh*conjg(eth_J)+&
            conjg(Uh)*conjg(beth_J))) +&
            Jh*conjg(eth_U) + K*beth_U


       ! interpolation of the News, uBondi and redshift (at mid level) 
       ! to the points of the grid of inertial observers at scri.

       Uyh = 0
       delta_RHS_B=0
       betaB = 0
       NewsB = 0
       OmegaB = 1

       do ip = 1, 2
          p_this = 0
          p_other = 0
          do j = 1, lsh(2)
             do i = 1, lsh(1)
                qh = dble(zEvolh(i,j,ip))
                ph = dimag(zEvolh(i,j,ip))
                ipm = mod( ip + patch_index(i,j,ip) - 1 ,2 ) +1
              if (stereo_mask(i,j).eq.1) then
                if ( ipm .eq. ip ) then
                   p_this = p_this + 1
                   qcoord_this(p_this) = qh
                   pcoord_this(p_this) = ph
                   reindexi_this(p_this) = i
                   reindexj_this(p_this) = j
                else
                   p_other = p_other + 1
                   qcoord_other(p_other) = qh
                   pcoord_other(p_other) = ph
                   reindexi_other(p_other) = i
                   reindexj_other(p_other) = j
                end if
              end if
             end do
          end do
          ipo = mod( ip, 2 ) +1

          n_tmp1_cgfn = Uh(1:lsh(1), 1:lsh(2), ip)  ! not really n and s
          n_tmp1_cgfs = Uh(1:lsh(1), 1:lsh(2), ipo)
          n_tmp2_cgfn = delta_RHS(1:lsh(1), 1:lsh(2), ip)
          n_tmp2_cgfs = delta_RHS(1:lsh(1), 1:lsh(2), ipo)
          n_tmp1_rgfn = beta(1:lsh(1), 1:lsh(2), ip)
          n_tmp1_rgfs = beta(1:lsh(1), 1:lsh(2), ipo)
          n_tmp3_cgfn = News(1:lsh(1), 1:lsh(2), ip)
          n_tmp3_cgfs = News(1:lsh(1), 1:lsh(2), ipo)
          n_tmp2_rgfn = omega(1:lsh(1), 1:lsh(2), ip)
          n_tmp2_rgfs = omega(1:lsh(1), 1:lsh(2), ipo)
          !Do We NEED TO CALL CCTK_SyncGroup?

          !call CCTK_INFO("About to Interp massive this")
          interp_coords(1)        = CCTK_PointerTo(qcoord_this)
          interp_coords(2)        = CCTK_PointerTo(pcoord_this)


          in_array_indices(1)     = gf_indices(1)
          in_array_indices(2)     = gf_indices(3)
          in_array_indices(3)     = gf_indices(7)
          out_arrays(1)           = CCTK_PointerTo(Iarr1_this)
          out_arrays(2)           = CCTK_PointerTo(Iarr2_this)
          out_arrays(3)           = CCTK_PointerTo(Iarr4_this)
          interp_coords_type_code = CCTK_VARIABLE_REAL
          out_array_type_codes(1) = CCTK_VARIABLE_COMPLEX
          out_array_type_codes(2) = CCTK_VARIABLE_COMPLEX
          out_array_type_codes(3) = CCTK_VARIABLE_COMPLEX

          N_dims=2
          N_in_arrays=3
          N_out_arrays=3
          N_interp_points = p_this

          call CCTK_InterpGridArrays (ierror, cctkGH, N_dims,&
               interp_handle, param_table_handle, coord_system_handle,&
               N_interp_points, interp_coords_type_code,&
               interp_coords,&
               N_in_arrays, in_array_indices(1:3),&
               N_out_arrays, out_array_type_codes(1:3), out_arrays(1:3))

          if (ierror .ne. 0 ) then
             call CCTK_WARN(0, "News Interp Error")
          endif

          in_array_indices(1)     = gf_indices(5)
          in_array_indices(2)     = gf_indices(9)
          out_arrays(1)           = CCTK_PointerTo(Rarr3_this)
          out_arrays(2)           = CCTK_PointerTo(Rarr5_this)
          interp_coords_type_code = CCTK_VARIABLE_REAL
          out_array_type_codes(1) = CCTK_VARIABLE_REAL
          out_array_type_codes(2) = CCTK_VARIABLE_REAL

          N_dims=2
          N_in_arrays=2
          N_out_arrays=2
          N_interp_points = p_this

          call CCTK_InterpGridArrays (ierror, cctkGH, N_dims,&
               interp_handle, param_table_handle, coord_system_handle,&
               N_interp_points, interp_coords_type_code,&
               interp_coords,&
               N_in_arrays, in_array_indices(1:2),&
               N_out_arrays, out_array_type_codes(1:2), out_arrays(1:2))

          if (ierror .ne. 0 ) then
             call CCTK_WARN(0, "News Interp Error")
          endif

          !call CCTK_INFO("About to Interp massive other")
          interp_coords(1)        = CCTK_PointerTo(qcoord_other)
          interp_coords(2)        = CCTK_PointerTo(pcoord_other)
          in_array_indices(1)     = gf_indices(2)
          in_array_indices(2)     = gf_indices(4)
          in_array_indices(3)     = gf_indices(8)
          out_arrays(1)           = CCTK_PointerTo(Iarr1_other)
          out_arrays(2)           = CCTK_PointerTo(Iarr2_other)
          out_arrays(3)           = CCTK_PointerTo(Iarr4_other)
          interp_coords_type_code = CCTK_VARIABLE_REAL
          out_array_type_codes(1) = CCTK_VARIABLE_COMPLEX
          out_array_type_codes(2) = CCTK_VARIABLE_COMPLEX
          out_array_type_codes(3) = CCTK_VARIABLE_COMPLEX

          N_dims=2
          N_in_arrays=3
          N_out_arrays=3
          N_interp_points = p_other

          call CCTK_InterpGridArrays (ierror, cctkGH, N_dims,&
               interp_handle, param_table_handle, coord_system_handle,&
               N_interp_points, interp_coords_type_code,&
               interp_coords,&
               N_in_arrays, in_array_indices(1:3),&
               N_out_arrays, out_array_type_codes(1:3), out_arrays(1:3))

          if (ierror .ne. 0 ) then
             call CCTK_WARN(0, "News Interp Error")
          endif

          !call CCTK_INFO("About to Interp massive other")
          interp_coords(1)        = CCTK_PointerTo(qcoord_other)
          interp_coords(2)        = CCTK_PointerTo(pcoord_other)
          in_array_indices(1)     = gf_indices(6)
          in_array_indices(2)     = gf_indices(10)
          out_arrays(1)           = CCTK_PointerTo(Rarr3_other)
          out_arrays(2)           = CCTK_PointerTo(Rarr5_other)
          interp_coords_type_code = CCTK_VARIABLE_REAL
          out_array_type_codes(1) = CCTK_VARIABLE_REAL
          out_array_type_codes(2) = CCTK_VARIABLE_REAL

          N_dims=2
          N_in_arrays=2
          N_out_arrays=2
          N_interp_points = p_other

          call CCTK_InterpGridArrays (ierror, cctkGH, N_dims,&
               interp_handle, param_table_handle, coord_system_handle,&
               N_interp_points, interp_coords_type_code,&
               interp_coords,&
               N_in_arrays, in_array_indices(1:2),&
               N_out_arrays, out_array_type_codes(1:2), out_arrays(1:2))

          if (ierror .ne. 0 ) then
             call CCTK_WARN(0, "News Interp Error")
          endif

          do l = 1, p_this
             i = reindexi_this(l)
             j = reindexj_this(l)
             Uyh(i, j, ip)         = Iarr1_this(l)
             delta_RHS_B(i, j, ip) = Iarr2_this(l)
             betaB(i, j, ip)       = Rarr3_this(l)
             NewsB(i, j, ip)       = Iarr4_this(l)
             OmegaB(i, j, ip)      = Rarr5_this(l)

          end do
          do l = 1, p_other
             i = reindexi_other(l)
             j = reindexj_other(l)
             Uyh(i, j, ip)         = Iarr1_other(l)
             delta_RHS_B(i, j, ip) = Iarr2_other(l)
             betaB(i, j, ip)       = Rarr3_other(l)
             NewsB(i, j, ip)       = Iarr4_other(l)
             OmegaB(i, j, ip)      = Rarr5_other(l)
          end do

       end do


       delta_RHS_B = .5d0 * delta_RHS_B  +  uyh*conjg(zEvolh) 
       deltan = deltao + dt * dimag(delta_RHS_B)

    else ! linearized

       zEvoln = nzeta
       Uyh = Uh
       Uyn = Un
       NewsB = News
       deltao  = 0
       betaB = beta
       omegaB = omega
       zEvoln = nzeta
       deltan = 0

    endif ! linearized

    redshiftB = omegaB * exp(2*betaB) 
    uBondin = uBondio + dt * omegaB * exp(2*betaB)
    uBondi = .5d0 * (uBondin + uBondio)


    delta = .5d0 * (deltan + deltao)
    NewsB = NewsB * ( cos(2.0d0*delta) - ii*sin(2.0d0*delta) )

    ! update the time levels of Uy (the "shifted" array of U)

    zEvolo = zEvoln
    Uyo = Uyn
    deltao = deltan
    uBondio = uBondin

    do ip = 1, 2
       do j = 1, lsh(2)
          do i = 1, lsh(1)
             if ( dble(zEvolo(i,j,ip) * conjg(zEvolo(i,j,ip))) > 1.0d0 ) then
                patch_index(i,j,ip) = mod(patch_index(i,j,ip) +1, 2)  ! flip patch (or back, as the case may be)
                ! e^{i delta} does have spinweight with respect to this n/s flipping
                deltao(i,j,ip) = deltao(i,j,ip) -2.0d0 *atan2(dimag(zEvolo(i,j,ip)), dble(zEvolo(i,j,ip)))&
                     + pii - 2 * pii * nint(deltao(i,j,ip)/(2.0d0*pii)) !!!! this phase factor  could be ingored
                zEvolo(i,j,ip) = 1.0d0 / zEvolo(i,j,ip)
                Uyo(i,j,ip) = Uyo(i,j,ip) * (-zEvolo(i,j,ip) /conjg(zEvolo(i,j,ip)))
              
             end if
          end do
       end do
    end do


  end subroutine NullNews_News2BondiNews


end module NullNews_Bondi
