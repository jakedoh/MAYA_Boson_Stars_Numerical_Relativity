c /*@@
c@routine    interp2d
c@date       Fri Feb 14 08:46:53 1997
c@author     Paul Walker
c@desc 
c      Interpolates from 2D data var with coordinates x and y and
c      sizes nx and ny onto 1D data out with position outx and outy
c      and nout points.
c     <p>
c     This has linear, quadratic and cubic interpolators in it.
c     Or will one day very soon.
c@enddesc 
c@calls     
c@calledby   numerical_axig
c@@*/
      
      subroutine interp2d(var,x,y,nx,ny,out,outx,outy,nout,order)
      implicit none 
      integer nx,ny,nout
      real*8 var(nx,ny), x(nx), y(ny)
      real*8 out(nout),outx(nout),outy(nout)
      integer order
c Interpolation goes from ibelow to ibelow+[1,2,3] depending on order
      integer i,j,ibelow,jbelow,pt
      real*8 xsym,ysym,findx,findy,frac
      real*8 ydir(order+1)
      real*8 ft(10), xt(10)
      real*8 poly2inter, quad_2d, cubic_2d
      real*8 dx, dy, PI
      integer twobhjsad

      PI = 3.14159265

c     Set up the grid spacings
      dx = x(2) - x(1)
      dy = y(2) - y(1)

c Loop over all out points
      do pt=1,nout
         ysym = 1.0D0
         xsym = 1.0D0
c Check bounds
         findx = outx(pt)
         if (findx .lt. x(1)) then 
            write (*,*) "Below inner bound at ",pt,outx(pt)
            STOP
         endif
         if (findx .gt. x(nx)) then
            write (*,*) "Above x bounds at ",pt,outx(pt),x(nx)
            STOP
         endif
         findy = outy(pt)
         if (findy .lt. y(1)) then 
            write (*,*) "Below  y inner bound at ",pt,outy(pt)
            STOP
         endif
         if (findy .gt. y(ny)) then 
            write (*,*) "Below  y inner bound at ",pt,outy(pt)
            STOP
         endif

c Locate ourselves in i,j space
c         do i=1,nx
c            if (x(i) .lt. findx) then
c               ibelow = i
c            endif
c         enddo
c     Assume a regular grid
         ibelow = (findx-x(1))/dx+1

         if (order .eq. 3 .and. ibelow .gt. 1) then
            ibelow = ibelow - 1
         endif
         if (ibelow + order .gt. nx) then
            ibelow = nx - order
         endif

c         do i=1,ny
c            if (y(i) .lt. findy) then
c               jbelow = i
c            endif
c         enddo
c     Assume a regular grid
         jbelow = (findy-y(1))/dy+1
         if (order .eq. 3 .and. jbelow .gt. 1) then
            jbelow = jbelow - 1
         endif
         if (jbelow + order .gt. ny) then
            jbelow = ny - order
         endif

c         write (*,*) "PT :",findx,findy
c         write (*,*) "SYM:",sym
c         write (*,*) "BOUND X ",ibelow,x(ibelow),x(ibelow+1)
c         write (*,*) "BOUND Y ",jbelow,y(jbelow),y(jbelow+1)


c     So do the interpolation
         if (order .eq. 1) then 
c     Interp in the x direction
            frac = (findx-x(ibelow))/(x(ibelow+1)-x(ibelow))
            ydir(1) = frac * var(ibelow+1,jbelow) + 
     $           (1.0 - frac)*var(ibelow,jbelow)
            ydir(2) = frac * var(ibelow+1,jbelow+1) + 
     $           (1.0 - frac)*var(ibelow,jbelow+1)
c     And now in the y direction
            frac = (findy-y(jbelow))/(y(jbelow+1)-y(jbelow))
            out(pt) = xsym * ysym * 
     $           (frac * ydir(2) + (1.0 - frac) * ydir(1))
         else if (order .eq. 2) then
c     Load up for calls to poly2inter
            do j=1,3
               do i=1,3
                  ft(i) = var(ibelow+i-1,jbelow+j-1)
                  xt(i) = x(ibelow+i-1)
               enddo
               ydir(j) = quad_2d(ft,xt(1),dx,findx)
            enddo
            do j=1,3
               xt(j) = y(jbelow+j-1)
            enddo
            out(pt) = xsym * ysym*quad_2d(ydir,xt(1),dy,findy)
         else if (order .eq. 3) then
c     Load up for calls to cubic
            do j=1,4
               do i=1,4
                  ft(i) = var(ibelow+i-1,jbelow+j-1)
                  xt(i) = x(ibelow+i-1)
               enddo
               ydir(j) = cubic_2d(ft,xt(1),dx,findx)
            enddo
            do j=1,4
               xt(j) = y(jbelow+j-1)
            enddo
            out(pt) = xsym * ysym*cubic_2d(ydir,xt(1),dy,findy)
         else
            write (*,*) "ORDER set wrong in interp2d",order
            stop
         endif
      enddo
         
      return
      end
      
      real*8 function linear_2d(f, x0, dx, findx)
      implicit none
      real*8 f(2),x0,dx,findx
      real*8 frac
      
      frac = (findx-x0)/dx
      linear_2d = (frac)*f(2) + (1.0-frac)*f(1)
         
      return
      end
      
      real*8 function quad_2d(f, x0, dx, findx)
      implicit none
      real*8 f(3),x0, dx, findx
      real*8 f0,f1,f2
      real*8 a,b,c, dx2, x02, o2dx2
c     Mathematica tells us
c     -  List(List(Rule(c,(2*dx**2*f0 + 3*dx*f0*x0 - 4*dx*f1*x0 
c     -        + dx*f2*x0 +  f0*x0**2 - 2*f1*x0**2 + f2*x0**2)
c     -        /(2*dx**2)), Rule(b,(-3*dx*f0 + 4*dx*f1 - dx*f2 - 2*f0*x0 
c     -        + 4*f1*x0 - 2*f2*x0)/(2*dx**2)),Rule(a,(f0 - 2*f1 + 
c     -        f2)/(2*dx**2))))
      
      f0 = f(1)
      f1 = f(2)
      f2 = f(3)
      dx2 = dx**2
      x02 = x0**2
      o2dx2 = 1.0D0/(2.0D0*dx2)
      
      c = (2.0D0*dx2*f0 + dx*x0*(3.0D0*f0 - 4.0D0*f1 + f2) + 
     $     x02*(f0 - 2.0D0*f1 + f2))*o2dx2
      b = (dx * (-3.0D0*f0 + 4.0D0*f1 - f2) + x0 * (- 2.0D0*f0 + 
     $     4.0D0*f1 - 2.0D0*f2))*o2dx2
      
      a = (f0 - 2.0D0*f1 + f2)*o2dx2
      
      quad_2d = a*findx**2 + b*findx + c
      
      return
      end
      
      real*8 function cubic_2d(f, x0, dx, findx)
      implicit none
      real*8 a,b,c,d
      real*8 f(4),x0,dx,findx
      
      a = -(f(1)-3.0*f(2)+3.0*f(3)-f(4)) / (6.0*(dx**3))
      
      b = (f(1)-2.0*f(2)+f(3))/(2.0*(dx**2)) +
     $     (f(1)-3.0*f(2)+3.0*f(3)-f(4))*(dx+x0)/(2.0*(dx**3))
      
      c = ((dx**2)*(-11.0*f(1) + 18.0*f(2) -  9.0*f(3) +  2.0*f(4)) +
     $     dx*x0*  (-12.0*f(1) + 30.0*f(2) - 24.0*f(3) +  6.0*f(4)) +
     $     (x0**2)*( -3.0*f(1) +  9.0*f(2) -  9.0*f(3) +  3.0*f(4))) /
     $     (6.0*(dx**3))
      
      d = ((dx**3)*   (  6.0*f(1) ) +
     $     (dx**2)*x0*( 11.0*f(1) - 18.0*f(2) +  9.0*f(3) -  2.0*f(4)) +
     $     (x0**2)*dx*(  6.0*f(1) - 15.0*f(2) + 12.0*f(3) -  3.0*f(4)) +
     $     (x0**3)*   (  1.0*f(1) -  3.0*f(2) +  3.0*f(3) -  1.0*f(4)))/
     $     (6.0*(dx**3))
         
      cubic_2d = ((a*findx + b)*findx + c)*findx + d
            
      return
      end


