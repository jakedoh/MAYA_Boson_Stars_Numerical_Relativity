! vim: syntax=fortran
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!debugging using the Schwarzchild metric solution in IEF coords
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       g(1,1)%d = 1+8*mass*qs**2/pp**2/cr
       g(1,2)%d = 8*mass/cr*qs*(3-2*ip)*ps/pp**2
       g(1,3)%d = -4*mass/cr*qs*(3-2*ip)*(-1+qs**2+ps**2)/pp**2
       g(2,2)%d = 1+8*mass*ps**2/pp**2/cr
       g(2,3)%d = -4*mass*(-1+qs**2+ps**2)*ps/cr/pp**2
       g(3,3)%d = 1+2*mass*(-1+qs**2+ps**2)**2/pp**2/cr

       beta(1)%d = -4*mass*qs/pp/(cr+2*mass)
       beta(2)%d = -4*mass*(3-2*ip)*ps/pp/(cr+2*mass)
       beta(3)%d = 2*mass*(3-2*ip)*(-1+qs**2+ps**2)/pp/(cr+2*mass) 

       alpha%d = 1/sqrt(cr+2*mass)*sqrt(cr) 

       dg(1,1,1)%d = -8*mass*qs**2/cr**2/pp**2 
       dg(1,2,1)%d = -(3-2*ip)*8*mass*qs*ps/cr**2/pp**2
       dg(1,3,1)%d = (3-2*ip)*4*mass*(-1+qs**2+ps**2)*qs/cr**2/pp**2
       dg(2,2,1)%d = -8*mass*ps**2/cr**2/pp**2
       dg(2,3,1)%d = 4*mass*(-1+qs**2+ps**2)*ps/cr**2/pp**2
       dg(3,3,1)%d = -2*mass*(-1+qs**2+ps**2)**2/pp**2/cr**2

       dg(1,1,2)%d = -16*qs*mass*(-1+qs**2-ps**2)/pp**3/cr
       dg(1,2,2)%d = -(3-2*ip)*8*ps*mass*(-1+3*qs**2-ps**2)/pp**3/cr
       dg(1,3,2)%d = (3-2*ip)*4*mass*(1+qs**4-ps**4-6*qs**2)/pp**3/cr
       dg(2,2,2)%d = -32*mass*ps**2/pp**3/cr*qs
       dg(2,3,2)%d = 8*mass*qs*ps*(-3+qs**2+ps**2)/pp**3/cr
       dg(3,3,2)%d = 16*mass*qs*(-1+qs**2+ps**2)/pp**3/cr

       dg(1,1,3)%d = -32*mass/cr*qs**2*ps/pp**3 
       dg(1,2,3)%d = (3-2*ip)*8*qs*mass*(1+qs**2-3*ps**2)/pp**3/cr
       dg(1,3,3)%d = (3-2*ip)*8*mass*qs*ps*(-3+qs**2+ps**2)/pp**3/cr
       dg(2,2,3)%d = 16*ps*mass*(1+qs**2-ps**2)/pp**3/cr
       dg(2,3,3)%d = -4*mass*(6*ps**2-ps**4-1+qs**4)/pp**3/cr
       dg(3,3,3)%d = 16*mass*(-1+qs**2+ps**2)*ps/pp**3/cr

       dg(1,1,4)%d = 0.d0
       dg(1,2,4)%d = 0.d0
       dg(1,3,4)%d = 0.d0
       dg(2,2,4)%d = 0.d0
       dg(2,3,4)%d = 0.d0
       dg(3,3,4)%d = 0.d0 

       dbeta(1,1)%d = 4*mass*qs/(cr+2*mass)**2/pp 
       dbeta(2,1)%d = (3-2*ip)*4*mass*ps/(cr+2*mass)**2/pp
       dbeta(3,1)%d = -(3-2*ip)*2*mass/(cr+2*mass)**2*(-1+qs**2+ps**2)/pp

       dbeta(1,2)%d = 4*mass*(-1+qs**2-ps**2)/(cr+2*mass)/pp**2
       dbeta(2,2)%d = (3-2*ip)*8*mass*qs/(cr+2*mass)/pp**2*ps
       dbeta(3,2)%d = (3-2*ip)*8*mass*qs/(cr+2*mass)/pp**2

       dbeta(1,3)%d = 8*mass*qs/(cr+2*mass)/pp**2*ps
       dbeta(2,3)%d = -(3-2*ip)*4*mass*(1+qs**2-ps**2)/(cr+2*mass)/pp**2
       dbeta(3,3)%d = (3-2*ip)*8*mass*ps/(cr+2*mass)/pp**2

       dbeta(1,4)%d = 0.d0 
       dbeta(2,4)%d = 0.d0 
       dbeta(3,4)%d = 0.d0 

       dalpha(1)%d = mass/sqrt(cr+2*mass)**3/sqrt(cr) 
       dalpha(2)%d = 0.d0
       dalpha(3)%d = 0.d0      
       dalpha(4)%d = 0.d0       
