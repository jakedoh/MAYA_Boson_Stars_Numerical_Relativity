
      PROGRAM energy

c     -----------------------------------------------------------------
c
c     Calculate energy flux from even parity variable (only real part)
c
c     Now copes with having unequally spaced timesteps, it works out
c     the smallest timestep, and uses this to interpolate to a new
c     array with equally spaced timesteps.
c
c     -----------------------------------------------------------------

      IMPLICIT NONE

      INTEGER,PARAMETER :: nmax=5000
      INTEGER :: i,itotal,it1,it2,Nt
      CHARACTER*50 filename
      REAL*8 t(nmax),Qplus(nmax),dEdt(nmax),E,ddQ(nmax),Qnew(nmax)
      REAL*8 Dtnew,fac,time,t1,t2,dQdt_l,dQdt_r,dtmin,dt
      INTEGER :: ioerror 
      WRITE(6,*)
      WRITE(6,*) "Starting PROGRAM energy ..."
      WRITE(6,*) "---------------------------"

      WRITE(6,131) 
 131  FORMAT("Filename containing Q ? ",$)
      READ(5,*) filename
      WRITE(6,*) "Filename is ... ",filename

      OPEN(UNIT=12,FILE=filename,STATUS="old")

      DO i=1,nmax
         READ(12,*,ERR=101, IOSTAT=ioerror) t(i),Qplus(i)
	 if(ioerror /= 0) goto 101
         itotal = i
      END DO
      WRITE(6,*) "Didn't read all of data file"
      STOP

 101  CLOSE(12) 

c     Give information about the data
      WRITE(6,*) "We have data from t=",t(1)," to ",t(itotal)

c     Find minimum time step
      dtmin = 100000.0
      do i=2,itotal
         dt = t(i)-t(i-1)
         if (dt<dtmin) dtmin = dt
      end do
      WRITE(6,*) "Minimum timestep is ",dtmin

c     Initial time for energy flux
 201  WRITE(6,132) 
 132  FORMAT(" Initial time ? ",$)
      READ(5,*) t1
      IF (t1 < t(1)) GOTO 201
     
c     Calculate initial timeslice
      do i=1,itotal
         if (t(i)<t1) it1=i
      end do
      WRITE(6,*) "Initial timeslice is ...",it1

c     Final time for energy flux
      WRITE(6,133) 
 133  FORMAT(" Final time ? ",$)
 301  READ(5,*) t2
      IF (t2 > t(itotal)) GOTO 301

c     Calculate final timeslice
      do i=1,itotal
         if (t(i)<t2) it2=i
      end do
      WRITE(6,*) "Final timeslice is ...",it2


      Nt = (t2-t1)/dtmin
      WRITE(6,*) "Interpolating to ",Nt," timesteps"
      

c     Do a spline
      Dtnew = Dtmin
      dQdt_l = ( Qplus(2)             - Qplus(1)                ) / 
     &           ( t(2)       - t(1)          )
      dQdt_r = (Qplus(itotal)       - Qplus(itotal-1)       ) /
     &           (t(itotal) - t(itotal-1)) 
      CALL SPLINE(t,Qplus,itotal,dQdt_l,dQdt_r,ddQ)  
      DO i = 1,Nt 
        time = t1+DBLE(i-1)*Dtnew
        CALL SPLINT(t,Qplus,ddQ,itotal,time,Qnew(i))
      ENDDO

c     This is the normalisation used
      fac = 1.0D0/32.0D0/ACOS(-1.0D0)

c     Calculate energy flux
      dEdt(1) = fac*((Qnew(2)-Qnew(1))/Dtnew)**2
      DO i = 2,Nt-1
         dEdt(i) = fac*((Qnew(i+1)-Qnew(i-1))*0.5D0/Dtnew)**2
      END DO
      dEdt(Nt) = fac*((Qnew(Nt)-Qnew(Nt-1))/Dtnew)**2

c     Write results to file
      OPEN(UNIT=21,FILE="dEdt",STATUS="unknown")
      DO i = 1,Nt
         time = t1+DBLE(i-1)*Dtnew
         WRITE(21,*) time,dEdt(i)
      END DO
      CLOSE(21)

c     Use Simpsons rule to calculate integral
      IF ( MOD(Nt,2) /= 0) THEN
         WRITE(6,*) "Using ",Nt," timeslices. Removing one" 
         WRITE(6,*) "to give Simpsons rule an odd number ..."
      END IF

      E = dEdt(1) + 4.0D0*dEdt(Nt-1)+dEdt(Nt)
      DO i=2,Nt-3,2
         E = E + 4.0D0*dEdt(i) + 2.0D0*dEdt(i+1)
      ENDDO
      E = E*Dtnew/3.0D0

c     Write out total energy
      WRITE(6,*)
      WRITE(6,*) "Total energy in this mode is : ",E
      WRITE(6,*) "****************************"

      END PROGRAM energy




      SUBROUTINE spline(x,y,n,yp1,ypn,y2)

      INTEGER, PARAMETER :: nmax=5000

      INTEGER   n,i,k
      REAL*8    yp1,ypn,x(nmax),y(nmax),y2(nmax)
      REAL*8    p,qn,sig,un,u(nmax)

      IF (yp1.GT..99d30) THEN
        y2(1)=0.d0
        u(1)=0.d0
      ELSE
        y2(1)=-0.5d0
        u(1)=(3.d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      ENDIF

      DO i=2,n-1
        sig   = (x(i)-x(i-1))/(x(i+1)-x(i-1))
        p     = sig*y2(i-1)+2.d0
        y2(i) = (sig-1.d0)/p
        u(i)  = (6.d0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))
     &          /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
      ENDDO

      IF (ypn.GT..99d30) then
        qn=0.d0
        un=0.d0
      ELSE
        qn=0.5d0
        un=(3.d0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      ENDIF

      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.d0)

      DO k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
      ENDDO

      RETURN
      END



      SUBROUTINE splint(xa,ya,y2a,n,x,y)
 
      INTEGER n 
      REAL*8  x,y,xa(n),y2a(n),ya(n) 
      INTEGER k,khi,klo 
      REAL*8  a,b,h 

      klo=1 
      khi=n 
1     if (khi-klo.gt.1) then 
        k=(khi+klo)/2 
        if(xa(k).gt.x)then 
          khi=k 
        else 
          klo=k 
        endif 
      goto 1 
      endif 
      h=xa(khi)-xa(klo) 
      if (h.eq.0.d0) pause 'bad xa input in splint' 
      a=(xa(khi)-x)/h 
      b=(x-xa(klo))/h 
      y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h** 
     *2)/6.d0 

      RETURN 
      END 


