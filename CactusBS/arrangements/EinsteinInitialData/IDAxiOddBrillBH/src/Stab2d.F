      subroutine bicgst2d(cc,ce,cw,cn,cs,u,rhs,
     &             eps,rmax,ier,im,jm)
c
c This routine was lifted from stab.f. Minor modifications have
c been made.
c
      implicit none
c
      integer,intent(in) :: im,jm
      real*8,intent(inout) :: cc(im,jm),cn(im,jm),cs(im,jm),
     $     ce(im,jm),cw(im,jm)
      real*8,intent(out) :: eps
      real*8,intent(out) :: rmax
      real*8 :: u(im,jm),rhs(im,jm)
c     Local variable
      integer ncyc
      integer iscale,i,j,k,ier
c
c******************************************************************************
c
      ncyc = (im-2)*(jm-2)
      
      ier=0

*
*  Determine whether we can diagonally scale the problem to speed
*  convergence. Can only be done if there are no zeros on the main
*  diagonal (ie. central difference coefficient).
*
      iscale=1
      do j = 2,jm-1
         do i = 2,im-1
            if (cc(i,j).eq.0.) iscale = 0
         enddo
      enddo
      
*
*  Do the diagonal scaling if we can.
*
      if (iscale.eq.1) then


         do j = 2,jm-1
            do i = 2,im-1
               rhs(i,j)=rhs(i,j)/cc(i,j)
               cw(i,j)=cw(i,j)/cc(i,j)
               ce(i,j)=ce(i,j)/cc(i,j)
               cs(i,j)=cs(i,j)/cc(i,j)
               cn(i,j)=cn(i,j)/cc(i,j)
               cc(i,j)=1.
            enddo
         enddo
         
      endif
      
*     
*   Now call the bicgstab routine
*     
      
      call bicgstab2d (cc,cn,cs,ce,cw,u,rhs,eps,ncyc,rmax,ier,
     &     im,jm)
      if (rmax.gt.eps) then
         ier = -1
         print*,'Did not converge'
         print*,' maximum residual    = ',rmax
         print*,' tolerance           = ',eps
      endif
      
      return
      end
c     
c******************************************************************
c
      subroutine bicgstab2d (cc,cn,cs,ce,cw,x,r,tol,ncyc,rnorm,
     $     ier,im,jm)
c     
c This routine was lifted from stab.f. Minor modifications have
c been made.
c
      implicit none
c     
      integer,intent(in) :: im,jm
      integer :: ncyc
      real*8,intent(in) :: cc(im*jm),cn(im*jm)
      real*8,intent(in) :: cs(im*jm),ce(im*jm)
      real*8,intent(in) :: cw(im*jm)
      real*8,intent(in) :: tol
      real*8,intent(out):: rnorm
      integer,intent(out) :: ier
      real*8 x(im*jm),r(im*jm)
c     Local variables
      integer :: i,j,k,kk
      real*8 :: p(im*jm),Ap(im*jm),w(im*jm),As(im*jm)
      real*8 :: omega, chi,chi1,chi2, delta, deltap, pp
*
***********************************************************************
*     

      do i = 1,im*jm
         p(i) = 0.
         Ap(i) = 0.
         w(i) = 0.
         As(i) = 0.
      enddo
      
      kk = 0
 10   call usermv2d(cc,cn,cs,ce,cw,x,Ap,im,jm)

      
      do i = 1,im*jm
         r(i) = r(i)-Ap(i)
         p(i) = r(i)
      enddo

c      delta = sum(r)
      delta = 0.

      do i = 1,im*jm
         delta = delta+r(i)
      enddo

      if (delta.eq.0.) then
        ier=-1
        return
      endif

      call usermv2d(cc,cn,cs,ce,cw,p,Ap,im,jm)

c      phi = sum(Ap)
      pp = 0.

      do i = 1,im*jm
         pp = pp+Ap(i)
      enddo
      pp = pp/delta

      if (pp.eq.0.) then
        ier=-1
        return
      endif

c      rnorm = sum(r**2)
      rnorm = 0.

      do i = 1,im*jm
         rnorm = rnorm + r(i)**2
      enddo
      rnorm=sqrt(rnorm)
      
c     Test if initial guess is great (residual less than tolerance)
      if (rnorm .lt. tol)  return
      
 1    continue
      kk = kk + 1

      omega = 1./pp


      do i = 1,im*jm
         w(i) = r(i) - omega*Ap(i)
      enddo

      call usermv2d(cc,cn,cs,ce,cw,w,As,im,jm)

c      chi1 = sum(As*w)
      chi1 = 0.

      do i = 1,im*jm
         chi1 = chi1+As(i)*w(i)
      enddo
c      chi2 = sum(As**2)
      chi2 = 0.

      do i = 1,im*jm
         chi2 = chi2+As(i)**2
      enddo

      chi = chi1/chi2


      do i = 1,im*jm
         r(i) = w(i) - chi*As(i)
         x(i) = x(i) + omega*p(i) + chi*w(i)
      enddo

      deltap = delta
c      delta = sum(r)
      delta = 0.

      do i = 1,im*jm
         delta = delta+r(i)
      enddo

      if (delta.eq.0.) then
        goto 10
      endif


      do i = 1,im*jm
         p(i) = r(i) + (p(i)-chi*Ap(i))*omega*
     &        delta/(deltap*chi)
      enddo
      
      call usermv2d(cc,cn,cs,ce,cw,p,Ap,im,jm)

c      phi = sum(Ap)
      pp = 0.

      do i = 1,im*jm
         pp = pp+Ap(i)
      enddo
      pp=pp/delta

      if (pp.eq.0.) then
        goto 10
      endif

      if (kk.gt.ncyc) then
        print*,' BI-CGStab solver reached maximum nuber of iterations.'
        ier=-1
        return
      endif

c      rnorm = sum(r**2)
      rnorm = 0.

      do i = 1,im*jm
         rnorm = rnorm+r(i)**2
      enddo
      rnorm=sqrt(rnorm)
            
      if (rnorm .gt. tol) goto 1
      return
      end
c
c******************************************************************
c
      subroutine usermv2d(cc,cn,cs,ce,cw,x,y,im,jm)
c
c This routine was lifted from stab.f. Minor modifications have
c been made.
c
c Be careful that the cs are zero on their outer boundary!!
c
      implicit none
c
      integer,intent(in) :: im,jm
      real*8,intent(in) :: cc(im*jm),cn(im*jm)
      real*8,intent(in) :: cs(im*jm),ce(im*jm)
      real*8,intent(in) :: cw(im*jm)
      real*8 :: x(im*jm), y(im*jm)
c     Local variables
      integer :: i, j, k
*
***********************************************************************
*
*
      do i = im+1,im*jm-im-1
         y(i) = cw(i)*x(i-1)
     &        +cc(i)*x(i)
     &        +ce(i)*x(i+1)
     &        +cn(i)*x(i+im)
     &        +cs(i)*x(i-im)
      enddo
c     
      return
      end
      
