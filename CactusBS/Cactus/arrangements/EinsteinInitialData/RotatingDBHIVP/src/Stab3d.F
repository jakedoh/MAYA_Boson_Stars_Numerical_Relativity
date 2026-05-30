      subroutine rbicgst3d(cc,ce,cw,cn,cs,ct,cb,u,rhs,
     &             eps,rmax,ier,im,jm,km)
c
c This routine was lifted from stab.f. Minor modifications have
c been made.
c
      implicit none
c
      integer,intent(in) :: im,jm,km
      real*8,intent(inout) :: cc(im,jm,km),cn(im,jm,km),cs(im,jm,km),
     $     ce(im,jm,km),cw(im,jm,km),ct(im,jm,km),cb(im,jm,km)
      real*8,intent(in) :: eps
      real*8,intent(out) :: rmax
      real*8,intent(inout) :: u(im,jm,km),rhs(im,jm,km)
c     Local variable
      integer ncyc
      integer iscale,i,j,k,ier
c
c******************************************************************************
c

      ncyc = (im-2)*(jm-2)*(km-2)
      
      ier=0

*
*  Determine whether we can diagonally scale the problem to speed
*  convergence. Can only be done if there are no zeros on the main
*  diagonal (ie. central difference coefficient).
*
      iscale=1
      do k = 2,km-1
         do j = 2,jm-1
            do i = 2,im-1
               if (cc(i,j,k).eq.0.) iscale = 0
            enddo
         enddo
      enddo

*
*  Do the diagonal scaling if we can.
*
      if (iscale.eq.1) then


         do k = 2,km-1
            do j = 2,jm-1
               do i = 2,im-1
                  rhs(i,j,k)=rhs(i,j,k)/cc(i,j,k)
                  cb(i,j,k)=cb(i,j,k)/cc(i,j,k)
                  ct(i,j,k)=ct(i,j,k)/cc(i,j,k)
                  cw(i,j,k)=cw(i,j,k)/cc(i,j,k)
                  ce(i,j,k)=ce(i,j,k)/cc(i,j,k)
                  cs(i,j,k)=cs(i,j,k)/cc(i,j,k)
                  cn(i,j,k)=cn(i,j,k)/cc(i,j,k)
                  cc(i,j,k)=1.
               enddo
            enddo
         enddo

      endif

* 
*   Now call the bicgstab routine
*
      
      call rbicgstab (cc,cn,cs,ce,cw,ct,cb,u,rhs,eps,ncyc,rmax,ier,
     &     im,jm,km)
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
      subroutine rbicgstab (cc,cn,cs,ce,cw,ct,cb,x,r,tol,ncyc,rnorm,
     $     ier,im,jm,km)
c
c This routine was lifted from stab.f. Minor modifications have
c been made.
c
      implicit none
c
      integer,intent(in) :: im,jm,km
      integer,intent(in) :: ncyc
      real*8,intent(in)  :: cc(im*jm*km),cn(im*jm*km)
      real*8,intent(in)  :: cs(im*jm*km),ce(im*jm*km)
      real*8,intent(in)  :: cw(im*jm*km),ct(im*jm*km)
      real*8,intent(in)  :: cb(im*jm*km)
      real*8,intent(in)  :: tol
      real*8,intent(out) :: rnorm
      integer,intent(out) :: ier
      real*8,intent(inout) :: x(im*jm*km),r(im*jm*km)
c     Local variables
      integer :: i,kk
      real*8, allocatable :: p(:),Ap(:),w(:),As(:)
      real*8 :: omega, chi,chi1,chi2, delta, deltap, pp
*
***********************************************************************
*
      allocate(p(im*jm*km),Ap(im*jm*km),w(im*jm*km),As(im*jm*km))

      do i = 1,im*jm*km
         p(i) = 0.
         Ap(i) = 0.
         w(i) = 0.
         As(i) = 0.
      enddo

      kk = 0
  10  call rusermv(cc,cn,cs,ce,cw,ct,cb,x,Ap,im,jm,km)

      do i = 1,im*jm*km
         r(i) = r(i)-Ap(i)
         p(i) = r(i)
      enddo

c      delta = sum(r)
      delta = 0.

      do i = 1,im*jm*km
         delta = delta+r(i)
      enddo

      if (delta.eq.0.) then
        ier=-1
        return
      endif

      call rusermv(cc,cn,cs,ce,cw,ct,cb,p,Ap,im,jm,km)

c      phi = sum(Ap)
      pp = 0.

      do i = 1,im*jm*km
         pp = pp+Ap(i)
      enddo
      pp = pp/delta

      if (pp.eq.0.) then
        ier=-1
        return
      endif

c      rnorm = sum(r**2)
      rnorm = 0.

      do i = 1,im*jm*km
         rnorm = rnorm + r(i)**2
      enddo
      rnorm=sqrt(rnorm)

c     Test if initial guess is great (residual less than tolerance)
      if (rnorm .lt. tol)  return
      
 1    continue
      kk = kk + 1

      omega = 1./pp


      do i = 1,im*jm*km
         w(i) = r(i) - omega*Ap(i)
      enddo

      call rusermv(cc,cn,cs,ce,cw,ct,cb,w,As,im,jm,km)

c      chi1 = sum(As*w)
      chi1 = 0.

      do i = 1,im*jm*km
         chi1 = chi1+As(i)*w(i)
      enddo
c      chi2 = sum(As**2)
      chi2 = 0.

      do i = 1,im*jm*km
         chi2 = chi2+As(i)**2
      enddo

      chi = chi1/chi2


      do i = 1,im*jm*km
         r(i) = w(i) - chi*As(i)
         x(i) = x(i) + omega*p(i) + chi*w(i)
      enddo

      deltap = delta
c      delta = sum(r)
      delta = 0.

      do i = 1,im*jm*km
         delta = delta+r(i)
      enddo

      if (delta.eq.0.) then
        goto 10
      endif


      do i = 1,im*jm*km
         p(i) = r(i) + (p(i)-chi*Ap(i))*omega*
     &        delta/(deltap*chi)
      enddo
      
      call rusermv(cc,cn,cs,ce,cw,ct,cb,p,Ap,im,jm,km)

c      phi = sum(Ap)
      pp = 0.

      do i = 1,im*jm*km
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

      do i = 1,im*jm*km
         rnorm = rnorm+r(i)**2
      enddo
      rnorm=sqrt(rnorm)

      if (rnorm .gt. tol) goto 1
c
      deallocate(p,Ap,w,As)
c
      return
      end
c
c******************************************************************
c
      subroutine rusermv(cc,cn,cs,ce,cw,ct,cb,x,y,im,jm,km)
c
c This routine was lifted from stab.f. Minor modifications have
c been made.
c
c Be careful that the cs are zero on their outer boundary!!
c
      implicit none
c
      integer,intent(in) :: im,jm,km
      real*8,intent(in)  :: cc(im*jm*km),cn(im*jm*km)
      real*8,intent(in)  :: cs(im*jm*km),ce(im*jm*km)
      real*8,intent(in)  :: cw(im*jm*km),ct(im*jm*km)
      real*8,intent(in)  :: cb(im*jm*km)
      real*8,intent(inout) :: x(im*jm*km), y(im*jm*km)
c     Local variables
      integer :: i
*
***********************************************************************
*
*
      do i = im*jm+im+1,im*(jm*km-jm-1)
         y(i) = cw(i)*x(i-1)
     &        +cc(i)*x(i)
     &        +ce(i)*x(i+1)
     &        +cn(i)*x(i+im)
     &        +cs(i)*x(i-im)
     &        +ct(i)*x(i+im*jm)
     &        +cb(i)*x(i-im*jm)
      enddo
c     
      return
      end
      
