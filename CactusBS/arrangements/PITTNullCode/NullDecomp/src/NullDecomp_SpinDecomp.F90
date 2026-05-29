/*
  vim: syntax=fortran
*/
#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

module NullDecomp_SpinDecomp
   use  NullDecomp_sYlm
   implicit none
   private

   CCTK_COMPLEX, dimension(:,:,:,:,:), allocatable, save :: stored_Ylm
   CCTK_INT, dimension(:,:,:), allocatable, save :: lm_to_stor
   CCTK_INT, save :: stored = 0
   CCTK_INT, save :: ntotal = 0
   logical, save :: IWasCalled = .false.

   public SpinDecompCoefs
   public SpinDecompRecon
   public SpinDecompFilter

 contains

!--------------------------------------------------------
! PUBLIC ROUTINES:
!--------------------------------------------------------

subroutine SpinDecompCoefs(GH, nq, np, s, z, f, sphcoef)
   use NullDecomp_Vars, only: lmax
   implicit none

   CCTK_POINTER, intent(in) :: GH
   CCTK_INT,                            intent (in)  :: s, nq, np
   CCTK_COMPLEX,   dimension (nq,np)  , intent (in)  :: z
   CCTK_COMPLEX,   dimension (nq,np,2), intent (in)  :: f

   CCTK_COMPLEX,   dimension(:,:),      intent (out) :: sphcoef
 
   CCTK_COMPLEX, dimension(nq, np, 2) :: temp
   
   CCTK_COMPLEX, dimension(s:lmax,-lmax:lmax) :: sphcoef_test

   CCTK_REAL :: norm

   DECLARE_CCTK_FUNCTIONS

   if (CCTK_IsThornActive("NullDecomp") .eq. 0) then
     call CCTK_WARN(0, "You must activate NullDecomp")
   endif

   if (any(shape(sphcoef)  .ne. shape(sphcoef_test) )) then
      call CCTK_WARN(0, "wrong input shape for 'sphcoef' -- must be (s:lmax,-lmax:lmax)")
   end if

   call wt_sphdec(GH, nq, np, s, lmax, z, f, temp, sphcoef, norm)

end subroutine SpinDecompCoefs


subroutine SpinDecompRecon(GH, nq, np, s, z, f, sphcoef)
   use NullDecomp_Vars, only: lmax
   implicit none

   CCTK_POINTER, intent(in) :: GH
   CCTK_INT,                            intent (in)  :: s, nq, np
   CCTK_COMPLEX,   dimension (nq,np)  , intent (in)  :: z
   CCTK_COMPLEX,   dimension (nq,np,2), intent (out) :: f

   CCTK_COMPLEX,   dimension(:,:),      intent (in)  :: sphcoef
 
   CCTK_COMPLEX, dimension(s:lmax,-lmax:lmax) :: sphcoef_test

   DECLARE_CCTK_FUNCTIONS

   if (CCTK_IsThornActive("NullDecomp") .eq. 0) then
     call CCTK_WARN(0, "You must activate NullDecomp")
   endif

   if (any(shape(sphcoef)  .ne. shape(sphcoef_test) )) then
      call CCTK_WARN(0, "wrong input shape for 'sphcoef' -- must be (s:lmax,-lmax:lmax)")
   end if

   call wt_recon(nq, np, s, lmax, z, f, sphcoef)

end subroutine SpinDecompRecon


subroutine SpinDecompFilter(GH, nq, np, s, z, f)
   use NullDecomp_Vars, only: lmax
   implicit none

   CCTK_POINTER, intent(in) :: GH
   CCTK_INT,                               intent (in)  :: s, nq, np
   CCTK_COMPLEX,   dimension (nq,np)  , intent (in)  :: z
   CCTK_COMPLEX,   dimension (nq,np,2), intent (inout)  :: f
 
   CCTK_COMPLEX, dimension(nq, np, 2) :: temp
   CCTK_COMPLEX, dimension(s:lmax,-lmax:lmax) :: sphcoef
   
   CCTK_REAL :: norm

   DECLARE_CCTK_FUNCTIONS

   if (CCTK_IsThornActive("NullDecomp") .eq. 0) then
     call CCTK_WARN(0, "You must activate NullDecomp")
   endif

   call wt_sphdec(GH, nq, np, s, lmax, z, f, temp,&
             sphcoef, norm)

   call wt_recon(nq, np, s, lmax, z, f, sphcoef)



end subroutine SpinDecompFilter

!--------------------------------------------------------
! PRIVATE ROUTINES:
!--------------------------------------------------------

subroutine wt_sphdec(GH, nq, np, s, lmax, z, f, ferr,&
             sphcoef, norm)

   !------------------------------------------------------------------
   ! This routine decomposes a function on the sphere into spherical 
   ! harmonics, returning its coeficients in sphcoef and the remainder
   ! in ferr. Ylm, fint are scratch arrays.
   !------------------------------------------------------------------
   use NullDecomp_Vars, only: lq, uq, lp, up, sum_handle, kern_p, area_p

   implicit none
   CCTK_POINTER, intent(in) :: GH
   CCTK_INT,                               intent (in)  :: s, nq, np, lmax
   CCTK_COMPLEX,   dimension (nq,np)  , intent (in)  :: z
   CCTK_COMPLEX,   dimension (nq,np,2), intent (in)  :: f
   CCTK_COMPLEX,   dimension (nq,np,2), intent (out) :: ferr
   CCTK_COMPLEX,   dimension (s:lmax,-lmax:lmax), intent (out) :: sphcoef
   CCTK_REAL, intent(out) :: norm
 
   CCTK_REAL :: gl_norm
   CCTK_REAL :: gl_sphc_re, gl_sphc_im
   CCTK_REAL :: lcoef_re, lcoef_im
   CCTK_INT :: err
   CCTK_INT :: l, m
   CCTK_COMPLEX,   dimension (nq,np,2) :: Ylm, fint
   CCTK_REAL, parameter :: pi = 3.141592653589793238460
   CCTK_COMPLEX :: t_norm

   DECLARE_CCTK_PARAMETERS

   if (l_max .ne. lmax) then
     call CCTK_WARN(0, "inconsistent call");
   end if

   if (s < 0) then
     call CCTK_WARN(0, "This routine reuqires s>-1 for a stupid reason ... I'll fix it if someone complains")
   endif

   ! we assume (stupidly) that s is positive
   if (s > lmax) then
     call CCTK_WARN(0, "wt_sphdec: invalid lmax")
   endif
   if (abs(s) > 2) then
     call CCTK_WARN(0, "Another stupid limitation ... s<3 required")
   endif

   if (.NOT.IWasCalled .AND. (store_ylms.ne.0) )then
     IWasCalled = .true.
     ntotal = ((lmax+2)*(lmax+1))/2

     allocate(stored_Ylm(0:2,nq,np,2,ntotal))
     allocate(lm_to_stor(0:2,0:lmax, -lmax:lmax))
     lm_to_stor = 0
     stored_Ylm = 1.0d100
   endif

   sphcoef = (0., 0.)

   fint = f * conjg(f) * kern_p


   call wt_intster(uq-lq+1, up-lp+1, area_p(lq:uq,lp:up),& 
         fint(lq:uq,lp:up,1) + fint(lq:uq,lp:up,2), t_norm)
   norm = dble(t_norm)

   call CCTK_ReduceLocScalar(err, GH, -1, sum_handle,&
        norm, gl_norm, CCTK_VARIABLE_REAL)
   if (err .ne. 0) then
     call CCTK_WARN(0,"Error getting sum")
   endif
   norm = gl_norm

   ferr = f

   do l = s, lmax
      do m = -l, l
         if (store_ylms.ne.0) then
           if ( lm_to_stor(s,l,m) .ne. 0) then
             Ylm = stored_Ylm(s,:,:,:,lm_to_stor(s,l,m))
           elseif (stored < ntotal) then
             if (use_rsYlm.ne.0) then
               call rsYlm(s,l, m, nq, np, z, Ylm)
             else
               call sYlm(s,l, m, nq, np, z, Ylm)
             endif
             stored = stored + 1
             lm_to_stor(s,l,m) = stored
             stored_Ylm(s,:,:,:,stored) = Ylm
           else
             if (use_rsYlm.ne.0) then
               call rsYlm(s,l, m, nq, np, z, Ylm)
             else
               call sYlm(s,l, m, nq, np, z, Ylm)
             endif
           endif
         else
           if (use_rsYlm.ne.0) then
             call rsYlm(s,l, m, nq, np, z, Ylm)
           else
             call sYlm(s,l, m, nq, np, z, Ylm)
           endif
         endif
      
         fint = ferr * conjg(Ylm) * kern_p

         call wt_intster(uq-lq+1, up-lp+1, area_p(lq:uq,lp:up),&
             fint(lq:uq,lp:up,1) + fint(lq:uq,lp:up,2), sphcoef(l,m))

         lcoef_re = dble(sphcoef(l,m))
         lcoef_im = dimag(sphcoef(l,m))
         call CCTK_ReduceLocScalar(err, GH, -1, sum_handle,&
            lcoef_re, gl_sphc_re, CCTK_VARIABLE_REAL)
         if (err .ne. 0) then
           call CCTK_WARN(0,"Error getting sum")
         endif
         call CCTK_ReduceLocScalar(err, GH, -1, sum_handle,&
            lcoef_im, gl_sphc_im, CCTK_VARIABLE_REAL)
         if (err .ne. 0) then
           call CCTK_WARN(0,"Error getting sum")
         endif
         sphcoef(l,m) = dcmplx(gl_sphc_re, gl_sphc_im)
         

         ! we could have combined all communcation together if we
         ! did not do the following. After the test phase we should move this
         ! part outside the loop (equivalent at continuum level) and commincate
         ! the enter sphcoef array
         !ferr(lq:uq,lp:up,:) = ferr(lq:uq,lp:up,:) -& 
         !   sphcoef(l,m) * Ylm(lq:uq,lp:up,:) 
         ferr = ferr - sphcoef(l,m) * Ylm

      end do
   end do

end subroutine wt_sphdec


subroutine wt_recon(nq, np, s, lmax, z, f, sphcoef)

   !------------------------------------------------------------------
   ! This routine decomposes a function on the sphere into spherical 
   ! harmonics, returning its coeficients in sphcoef and the remainder
   ! in ferr. Ylm, fint are scratch arrays.
   !------------------------------------------------------------------

   implicit none
   CCTK_INT,                               intent (in)  :: s, nq, np, lmax
   CCTK_COMPLEX,   dimension (nq,np)  , intent (in)  :: z
   CCTK_COMPLEX,   dimension (nq,np,2), intent (out)  :: f
   CCTK_COMPLEX,   dimension (s:lmax,-lmax:lmax), intent (in) :: sphcoef
 
   CCTK_INT :: l,m
   CCTK_COMPLEX,   dimension (nq,np,2) :: Ylm

   DECLARE_CCTK_PARAMETERS

   if (l_max .ne. lmax) then
     call CCTK_WARN(0, "inconsistent call");
   end if


   if (s < 0) then
     call CCTK_WARN(0, "This routine reuqires s>-1 for a stupid reason ... I'll fix it if someone complains")
   endif

   ! we assume (stupidly) that s is positive
   if (s > lmax) then
     call CCTK_WARN(0, "wt_sphdec: invalid lmax")
   endif
   if (abs(s) > 2) then
     call CCTK_WARN(0, "Another stupid limitation ... s<3 required")
   endif

   if (.NOT.IWasCalled .AND. (store_ylms.ne.0) )then
     IWasCalled = .true.
     ntotal = ((lmax+2)*(lmax+1))/2

     allocate(stored_Ylm(0:2,nq,np,2,ntotal))
     allocate(lm_to_stor(0:2,0:lmax, -lmax:lmax))
     lm_to_stor = 0
     stored_Ylm = 1.0d100
   endif

   f = 0

   do l = s, lmax
      do m = -l, l
         if (store_ylms.ne.0) then
           if ( lm_to_stor(s,l,m) .ne. 0) then
             Ylm = stored_Ylm(s,:,:,:,lm_to_stor(s,l,m))
           elseif (stored < ntotal) then
             if (use_rsYlm.ne.0) then
               call rsYlm(s,l, m, nq, np, z, Ylm)
             else
               call sYlm(s,l, m, nq, np, z, Ylm)
             endif
             stored = stored + 1
             lm_to_stor(s,l,m) = stored
             stored_Ylm(s,:,:,:,stored) = Ylm
           else
             if (use_rsYlm.ne.0) then
               call rsYlm(s,l, m, nq, np, z, Ylm)
             else
               call sYlm(s,l, m, nq, np, z, Ylm)
             endif
           endif
         else
           if (use_rsYlm.ne.0) then
             call rsYlm(s,l, m, nq, np, z, Ylm)
           else
             call sYlm(s,l, m, nq, np, z, Ylm)
           endif
         endif
      
         f = f + sphcoef(l,m) * Ylm

      end do
   end do

end subroutine wt_recon



subroutine wt_intster(nq, np, area, field, integral)
  implicit none
  CCTK_INT,          intent (in)  :: nq, np
  CCTK_REAL, intent (in)  :: area(nq,np)
  CCTK_COMPLEX,   intent (in)  :: field(nq,np)
  CCTK_COMPLEX,   intent (out) :: integral 

  CCTK_INT ::  l, k
  CCTK_COMPLEX sum

  sum = 0.
  do l = 1, np-1
     do k = 1, nq-1
        sum = sum + 0.25 * area(k,l) * &
        ( field(k,l) + field(k+1,l) + field(k+1,l+1) + field(k,l+1) )
     end do
  end do
  integral = sum

end subroutine wt_intster

end module NullDecomp_SpinDecomp
