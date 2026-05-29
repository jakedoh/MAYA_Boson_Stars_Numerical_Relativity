! vim: syntax=fortran

  subroutine my_nullify1(imin, imax, gf_1d)

    use NullSHRE_modGFdef 
    implicit none

    CCTK_INT,intent(in) :: imin, imax
    type (gf2d), intent(inout) :: gf_1d(imin:imax)
    integer :: i
    do i = imin, imax
       NULLIFY(gf_1d(i)%d)
    end do

  end subroutine my_nullify1


  subroutine my_nullify1c(imin, imax, gf_1d)

    use NullSHRE_modGFdef 
    implicit none

    CCTK_INT, intent(in) :: imin, imax
    type (gf2dc), intent(inout) :: gf_1d(imin:imax)
    integer :: i
    do i = imin, imax
     NULLIFY(gf_1d(i)%d)
    end do

  end subroutine my_nullify1c


  subroutine my_nullify2(imin, imax, jmin, jmax, gf_2d)

    use NullSHRE_modGFdef
    implicit none

    CCTK_INT, intent(in) :: imin, imax, jmin, jmax
    type (gf2d), intent(inout) :: gf_2d(imin:imax,jmin:jmax)
    integer :: i,j
    do j = jmin, jmax
       do i = imin, imax
          NULLIFY(gf_2d(i,j)%d)
       end do
    end do

  end subroutine my_nullify2


  subroutine my_nullify3(imin, imax, jmin, jmax, kmin, kmax, gf_3d)

    use NullSHRE_modGFdef
    implicit none

    CCTK_INT, intent(in) :: imin, imax, jmin, jmax, kmin, kmax
    type (gf2d), intent(inout) :: gf_3d(imin:imax,jmin:jmax,kmin:kmax)
    integer :: i,j,k
    do k = kmin, kmax
       do j = jmin, jmax
          do i =  imin, imax
             NULLIFY(gf_3d(i,j,k)%d)
          end do
       end do
    end do

  end subroutine my_nullify3


  subroutine nullify_all

    use NullSHRE_modGFdef
    implicit none

    CCTK_INT, parameter :: izero = 0
    integer, parameter :: ik = kind(izero)

    call my_nullify2(1_ik,4_ik,1_ik,4_ik, g_n)
    call my_nullify2(1_ik,4_ik,1_ik,4_ik, g_s)

    call my_nullify2(1_ik,4_ik,1_ik,4_ik, g1_n)
    call my_nullify2(1_ik,4_ik,1_ik,4_ik, g1_s)

    call my_nullify2(1_ik,4_ik,1_ik,4_ik, gup_n)
    call my_nullify2(1_ik,4_ik,1_ik,4_ik, gup_s)

    call my_nullify2(1_ik,4_ik,1_ik,4_ik, j0_n)
    call my_nullify2(1_ik,4_ik,1_ik,4_ik, j0_s)

    call my_nullify2(1_ik,4_ik,1_ik,4_ik, j1_n)
    call my_nullify2(1_ik,4_ik,1_ik,4_ik, j1_s)

    call my_nullify2(1_ik,4_ik,1_ik,4_ik, eta1_n)
    call my_nullify2(1_ik,4_ik,1_ik,4_ik, eta1_s)

    call my_nullify2(1_ik,4_ik,1_ik,4_ik, eta0_n)
    call my_nullify2(1_ik,4_ik,1_ik,4_ik, eta0_s)

    call my_nullify2(1_ik,4_ik,1_ik,4_ik, etaup0_n)
    call my_nullify2(1_ik,4_ik,1_ik,4_ik, etaup0_s)

    call my_nullify2(1_ik,4_ik,1_ik,4_ik, etaup1_n)
    call my_nullify2(1_ik,4_ik,1_ik,4_ik, etaup1_s)

    call my_nullify1c(2_ik,3_ik,qa)

    call my_nullify2(1_ik,3_ik,1_ik,4_ik, dbeta_n)
    call my_nullify2(1_ik,3_ik,1_ik,4_ik, dbeta_s)

    call my_nullify1(1_ik,3_ik,beta_n)
    call my_nullify1(1_ik,3_ik,beta_s)

    call my_nullify1(1_ik,4_ik,dalpha_n)
    call my_nullify1(1_ik,4_ik,dalpha_s)

    call my_nullify1(1_ik,3_ik, sigma_n)
    call my_nullify1(1_ik,3_ik, sigma_s)

    call my_nullify1(1_ik,4_ik,sa_n)
    call my_nullify1(1_ik,4_ik,sa_s)

    call my_nullify1(1_ik,4_ik,na_n)
    call my_nullify1(1_ik,4_ik,na_s)

    call my_nullify1(1_ik,4_ik,ell_n)
    call my_nullify1(1_ik,4_ik,ell_s)

    call my_nullify1(1_ik,4_ik,delld_n)
    call my_nullify1(1_ik,4_ik,delld_s)

    call my_nullify1(1_ik,4_ik,dr0_n)
    call my_nullify1(1_ik,4_ik,dr0_s)

    call my_nullify1(1_ik,2_ik,rl_old_n)
    call my_nullify1(1_ik,2_ik,rl_old_s)

    call my_nullify1(1_ik,4_ik,dr1_n)
    call my_nullify1(1_ik,4_ik,dr1_s)

    call my_nullify2(1_ik,3_ik,2_ik,3_ik,dsigma_n)
    call my_nullify2(1_ik,3_ik,2_ik,3_ik,dsigma_s)

    call my_nullify3(1_ik,4_ik,1_ik,4_ik,1_ik,4_ik,dg_n)
    call my_nullify3(1_ik,4_ik,1_ik,4_ik,1_ik,4_ik,dg_s)

    call my_nullify2(1_ik,4_ik,2_ik,4_ik,dsa_n)
    call my_nullify2(1_ik,4_ik,2_ik,4_ik,dsa_s)

    call my_nullify2(1_ik,4_ik,2_ik,4_ik,dna_n)
    call my_nullify2(1_ik,4_ik,2_ik,4_ik,dna_s)

    call my_nullify3(1_ik,3_ik,2_ik,3_ik,2_ik,3_ik,dj0_n)
    call my_nullify3(1_ik,3_ik,2_ik,3_ik,2_ik,3_ik,dj0_s)

    call my_nullify3(2_ik,3_ik,2_ik,3_ik,2_ik,3_ik,deta0_n)
    call my_nullify3(2_ik,3_ik,2_ik,3_ik,2_ik,3_ik,deta0_s)

    call my_nullify2(1_ik,4_ik,2_ik,3_ik,bracket_n)
    call my_nullify2(1_ik,4_ik,2_ik,3_ik,bracket_s)

    NULLIFY(detg_n%d, detg_s%d,&
         sigma2_n%d, sigma2_s%d,& 
         elld_n%d, elld_s%d,&
         r0_n%d, r0_s%d,&
         alpha_n%d, alpha_s%d,&
         x_wt_n%d, x_wt_s%d,&
         j_wt_n%d, j_wt_s%d,&
         j_l_n%d, j_l_s%d,&
         nu_wt_n%d, nu_wt_s%d,&
         nu_x_n%d, nu_x_s%d,&
         ck_wt_n%d, ck_wt_s%d,&
         ck_x_n%d, ck_x_s%d,&
         beta_wt_n%d, beta_wt_s%d,&
         beta_l_n%d, beta_l_s%d,&
         cb_wt_n%d, cb_wt_s%d,&
         cb_x_n%d, cb_x_s%d,&
         u_wt_n%d, u_wt_s%d,&
         u_l_n%d, u_l_s%d,&
         w_wt_n%d, w_wt_s%d,&
         w_l_n%d, w_l_s%d)

  end subroutine nullify_all
