! vim: syntax=fortran

  subroutine point_v1(lsh, vsize, i_min, i_max, gf_n, gf_s, v1)

    use NullSHRE_modGFdef 
    implicit none
 
    CCTK_INT,            intent(in)    :: lsh(2), vsize, i_min, i_max
    type(gf2d),          intent(inout) :: gf_n(i_min:i_max)
    type(gf2d),          intent(inout) :: gf_s(i_min:i_max)
    CCTK_REAL,   target, intent(in)    :: v1(lsh(1),lsh(2),vsize)

    CCTK_INT :: i, ll

    if (vsize.ne.2*(i_max-i_min+1))&
    call CCTK_WARN(0, "error vsize in pointing library(v1)")

    ll = 1
    do i = i_min, i_max

       gf_n(i)%d => v1(:,:,ll)
       gf_s(i)%d => v1(:,:,ll+1)

       ll = ll+2

    end do

    if (ll.ne.vsize+1)&
    call CCTK_WARN(0, "error ll in pointing library(v1)")

  end subroutine point_v1

  subroutine  point_v1sym(lsh, vsize, i_min, i_max, gf_n, gf_s, v1sym)

    use NullSHRE_modGFdef 
    implicit none

    CCTK_INT,           intent(in)    :: lsh(2), vsize, i_min, i_max
    type(gf2d),         intent(inout) :: gf_n(i_min:i_max, i_min:i_max)
    type(gf2d),         intent(inout) :: gf_s(i_min:i_max, i_min:i_max)
    CCTK_REAL,  target, intent(in)    :: v1sym(lsh(1),lsh(2), vsize)

    CCTK_INT :: i, ll

    if (vsize.ne.2*(i_max-i_min+1))&
    call CCTK_WARN(0, "error vsize in pointing library(v1sym)")
    ll = 1
    do i = i_min, i_max

          gf_n(i,4)%d  => v1sym(:,:,ll)
          gf_s(i,4)%d  => v1sym(:,:,ll+1)

          if (i.ne.4) then
             gf_n(4,i)%d  => v1sym(:,:,ll)
             gf_s(4,i)%d  => v1sym(:,:,ll+1)
          end if

          ll = ll+2

       end do

    if (ll.ne.vsize+1)&
    call CCTK_WARN(0, "error ll in pointing library(v1sym)")

  end subroutine point_v1sym

  subroutine  point_v2(lsh, vsize, i_min, i_max, j_min, j_max, gf_n, gf_s, v2)

    use NullSHRE_modGFdef 
    implicit none

    CCTK_INT,           intent(in)    :: lsh(2), vsize, i_min, i_max, j_min, j_max
    type(gf2d),         intent(inout) :: gf_n(i_min:i_max, j_min:j_max)
    type(gf2d),         intent(inout) :: gf_s(i_min:i_max, j_min:j_max)
    CCTK_REAL,  target, intent(in)    :: v2(lsh(1),lsh(2),vsize)

    CCTK_INT :: i, j, ll
    if (vsize.ne.2*(i_max-i_min+1)*(j_max-j_min+1))&
    call CCTK_WARN(0, "error vsize in pointing library(v2)")
    ll = 1
    do i = i_min, i_max
       do j = j_min, j_max

          gf_n(i,j)%d  => v2(:,:,ll)
          gf_s(i,j)%d  => v2(:,:,ll+1)

          ll = ll+2

       end do
    end do

    if (ll.ne.vsize+1)&
    call CCTK_WARN(0, "error ll in pointing library(v2)")

  end subroutine point_v2

  subroutine  point_v2sym(lsh, vsize, i_min, i_max, gf_n, gf_s, v2sym)

    use NullSHRE_modGFdef 
    implicit none

    CCTK_INT,           intent(in)    :: lsh(2), vsize, i_min, i_max
    type(gf2d),         intent(inout) :: gf_n(i_min:i_max, i_min:i_max)
    type(gf2d),         intent(inout) :: gf_s(i_min:i_max, i_min:i_max)
    CCTK_REAL,  target, intent(in)    :: v2sym(lsh(1),lsh(2), vsize)

    CCTK_INT :: i, j, ll

    if (vsize.ne.(i_max-i_min+1)*(i_max-i_min+2))&
    call CCTK_WARN(0, "error vsize in pointing library(v2sym)")
    ll = 1
    do i = i_min, i_max
       do j = i, i_max

          gf_n(i,j)%d  => v2sym(:,:,ll)
          gf_s(i,j)%d  => v2sym(:,:,ll+1)

          if (i.ne.j) then
             gf_n(j,i)%d  => v2sym(:,:,ll)
             gf_s(j,i)%d  => v2sym(:,:,ll+1)
          end if

          ll = ll+2

       end do
    end do

    if (ll.ne.vsize+1)&
    call CCTK_WARN(0, "error ll in pointing library(v2sym)")

  end subroutine point_v2sym

  subroutine  point_v3_sym12(lsh, vsize, i_min, i_max, k_min, k_max, gf_n, gf_s, v3_sym12)

    use NullSHRE_modGFdef 
    implicit none

    CCTK_INT,           intent(in)    :: lsh(2), vsize, i_min, i_max, k_min, k_max
    type(gf2d),         intent(inout) :: gf_n(i_min:i_max, i_min:i_max, k_min:k_max)
    type(gf2d),         intent(inout) :: gf_s(i_min:i_max, i_min:i_max, k_min:k_max)
    CCTK_REAL,  target, intent(in)    :: v3_sym12(lsh(1),lsh(2), vsize)

    CCTK_INT :: i, j, k, ll

    if (vsize.ne.(i_max-i_min+1)*(i_max-i_min+2)*(k_max-k_min+1))&
    call CCTK_WARN(0, "error vsize in pointing library(v3_sym12)")
    ll = 1
    do i = i_min, i_max
       do j = i, i_max
          do k = k_min, k_max

             gf_n(i,j,k)%d  => v3_sym12(:,:,ll)
             gf_s(i,j,k)%d  => v3_sym12(:,:,ll+1)

             if (i.ne.j) then
                gf_n(j,i,k)%d  => v3_sym12(:,:,ll)
                gf_s(j,i,k)%d  => v3_sym12(:,:,ll+1)
             end if

             ll = ll+2

          end do
       end do
    end do

    if (ll.ne.vsize+1)&
    call CCTK_WARN(0, "error ll in pointing library(v3_sym12)")

  end subroutine point_v3_sym12

  subroutine  point_v3_sym23(lsh, vsize, i_min, i_max, j_min, j_max, gf_n, gf_s, v3_sym23)

    use NullSHRE_modGFdef 
    implicit none

    CCTK_INT,           intent(in)    :: lsh(2), vsize, i_min, i_max, j_min, j_max
    type(gf2d),         intent(inout) :: gf_n(i_min:i_max, j_min:j_max, j_min:j_max)
    type(gf2d),         intent(inout) :: gf_s(i_min:i_max, j_min:j_max, j_min:j_max)
    CCTK_REAL,  target, intent(in)    :: v3_sym23(lsh(1),lsh(2),vsize)

    CCTK_INT :: i, j, k, ll

    if (vsize.ne.(j_max-j_min+1)*(j_max-j_min+2)*(i_max-i_min+1))&
    call CCTK_WARN(0, "error vsize in pointing library(v3_sym23)")
    ll = 1
    do i = i_min, i_max
       do j = j_min, j_max
          do k = j, j_max

             gf_n(i,j,k)%d  => v3_sym23(:,:,ll)
             gf_s(i,j,k)%d  => v3_sym23(:,:,ll+1)

             if (j.ne.k) then
                gf_n(i,k,j)%d  => v3_sym23(:,:,ll)
                gf_s(i,k,j)%d  => v3_sym23(:,:,ll+1)
             end if

             ll = ll+2

          end do
       end do
    end do

    if (ll.ne.vsize+1)&
    call CCTK_WARN(0, "error ll in pointing library(v3_sym23)")

  end subroutine point_v3_sym23
