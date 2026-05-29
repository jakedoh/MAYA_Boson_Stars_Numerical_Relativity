#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"



#define sx(i,j,k) scon(i,j,k,1)
#define sy(i,j,k) scon(i,j,k,2)
#define sz(i,j,k) scon(i,j,k,3)
#define Bconsx(i,j,k) Bcons(i,j,k,1)
#define Bconsy(i,j,k) Bcons(i,j,k,2)
#define Bconsz(i,j,k) Bcons(i,j,k,3)
#define Bvecx(i,j,k) Bvec(i,j,k,1)
#define Bvecy(i,j,k) Bvec(i,j,k,2)
#define Bvecz(i,j,k) Bvec(i,j,k,3)
#define velx(i,j,k) vel(i,j,k,1)
#define vely(i,j,k) vel(i,j,k,2)
#define velz(i,j,k) vel(i,j,k,3)


  subroutine Fluid_Wind(CCTK_ARGUMENTS)

! --->   INITIAL DATA   <---

! This routine calculates the initial data for

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS
  integer i,j,k
  CCTK_REAL zero,third,half,one,two,pi,Mass_my,g_up_xx,r1,det,det_gamma
  real(kind=8) pii
! --->   INITIAL DATA   <---
  CCTK_REAL M0,L0,K1,L1,K2a,L2a,K2b,L2b,K2c,L2c
  CCTK_REAL K2d,L2d,K3a,L3a,K3b,L3b,K3c,L3c,K3d,L3d,K4a
  CCTK_REAL L4a,K4b,L4b,K4c,L4c,K4d,L4d,K4e,L4e,K4f,L4f
  CCTK_REAL K4g,L4g,K4h,L4h,K4i,L4i
  CCTK_REAL M,M1,aa1,qq,a1,a2,m_1,m_2,aa,mm
  CCTK_REAL Z1,Z2,Spar,risco,Lisco,Eisco,aeff,dm,Dpar,eta
  
  
  
one=1.0d0
! --->   NUMBERS   <---

  zero = 0.0d0
  one  = 1.0d0
  half = 0.5d0
  two  = 2.0d0
  pi  = 4.0D0*atan(1.0D0)
 
 
 
 
  M0  =  0.951507d0
  L0  =  0.686710d0
  K1  = -0.051379d0
  L1  =  0.613247d0
  K2a = -0.004804d0
  L2a = -0.145427d0
  K2b = -0.054522d0
  L2b = -0.115689d0
  K2c = -0.000022d0
  L2c = -0.005254d0
  K2d =  1.995246d0
  L2d =  0.801838d0
  K3a =  0.007064d0
  L3a = -0.073839d0
  K3b = -0.017599d0
  L3b =  0.004759d0
  K3c = -0.119175d0
  L3c = -0.078377d0
  K3d =  0.025000d0
  L3d =  1.585809d0
  K4a = -0.068981d0
  L4a = -0.003050d0
  K4b = -0.011383d0
  L4b = -0.002968d0
  K4c = -0.002284d0
  L4c =  0.004364d0
  K4d = -0.165658d0
  L4d = -0.047204d0
  K4e =  0.019403d0
  L4e = -0.053099d0
  K4f =  2.980990d0
  L4f =  0.953458d0
  K4g =  0.020250d0
  L4g = -0.067998d0
  K4h = -0.004091d0
  L4h =  0.001629d0
  K4i =  0.078441d0
  L4i = -0.066693d0
  
  mm=M_total
  M1=Final_Mass
  aa1=Final_Spin
  
  print *, 'are we here?'
  open(1, file = 'brute.dat', status = 'new')
  write(1,*) '# Final Mass= ',M1,'Final_Spin= ',aa1
  write(1,*) '#  mass ratio  ','  a1  ','  a2  ','  1/mass ratio  '


! --->   INITIAL DATA   <---
   do k = 1, cctk_lsh(3)
   do j = 1, cctk_lsh(2)
   do i = 1, cctk_lsh(1)

  !! vel(i,j,k,1) is equivalent to vel^{x}    
  
         a1=y(i,j,k)
         a2=z(i,j,k)
         m_1 = mm * x(i,j,k) / ( 1.0d0 + x(i,j,k) )
         m_2 = mm / ( 1.0d0 + x(i,j,k) )
         dm = ( m_1 - m_2 ) / mm
         eta = m_1 * m_2 / mm**2
         Spar = ( a1 * m_1**2 + a2 * m_2**2 ) / mm**2
         aeff=Spar
         Dpar = (a2*m_2-a1*m_1)/mm
         Z1 = 1.0d0+(1.0d0-aeff**2)**(1.0d0/3.0d0)*((1.0d0+aeff)**(1.0d0/3.0d0)+(1.0d0-aeff)**(1.0d0/3.0d0))
         Z2 = sqrt(3.0d0*aeff**2+Z1**2)
         risco = 3.0d0+Z2-SIGN(one,aeff)*sqrt((3.0d0-Z1)*(3.0d0+Z1+2.0d0*Z2))
         Lisco = 2.0d0/sqrt(3.0d0*risco)*(3.0d0*sqrt(risco)-2.0d0*aeff)
         Eisco = (1.0d0-2.0d0/risco+aeff/risco**(1.5d0))/sqrt(1.0d0-3.0d0/risco+2.0d0*aeff/risco**(1.5d0))

         M=(4.0d0*eta)**2*(M0 + K1*Spar + K2a*Dpar*dm + K2b*Spar**2 + K2c*Dpar**2 + &
         K2d*dm**2 + K3a*Dpar*Spar*dm + K3b*Spar*Dpar**2 + K3c*Spar**3 + K3d*Spar*dm**2 + &
         K4a*Dpar*Spar**2*dm + K4b*Dpar**3*dm + K4c*Dpar**4 + K4d*Spar**4 + &
         K4e*Dpar**2*Spar**2 + K4f*dm**4 + K4g*Dpar*dm**3 + K4h*Dpar**2*dm**2 + K4i*Spar**2*dm**2) + (1.0d0+eta*(Eisco+11.0d0))*dm**6

         aa=(4*eta)**2*(L0 + L1*Spar+L2a*Dpar*dm+L2b*Spar**2+L2c*Dpar**2+L2d*dm**2 + &
         L3a*Dpar*Spar*dm+L3b*Spar*Dpar**2+L3c*Spar**3+L3d*Spar*dm**2+L4a*Dpar*Spar**2*dm+ &
         L4b*Dpar**3*dm+L4c*Dpar**4+L4d*Spar**4+L4e*Dpar**2*Spar**2+L4f*dm**4+& 
         L4g*Dpar*dm**3+L4h*Dpar**2*dm**2+L4i*Spar**2*dm**2) + Spar*(1.0d0+8.0d0*eta)*dm**4+eta*Lisco*dm**6
           IF ( abs((aa-aa1)/aa1) .LE. precision .AND. abs((M-M1)/M1) .LE. precision ) THEN
                               print*,'saving'
                                  write(1,*) 1.0d0/x(i,j,k),y(i,j,k),z(i,j,k),x(i,j,k)
				  print*,'Mass_ratio= ',x(i,j,k),'a1= ',y(i,j,k),'a1= ',z(i,j,k)


           END IF
  
         
        
         rho(i,j,k) = rho_ini_cons

         press(i,j,k) = C_s**2 * ( mygamma-one ) * rho(i,j,k) &
                      / (mygamma * ( mygamma - one ) - C_s**2 * mygamma )

!===============>>>>>> contravariant component of gamma_xx===============  
         g_up_xx = -(gyz(i,j,k)**2 - gyy(i,j,k)*gzz(i,j,k))&
                   /(-(gyy(i,j,k)*gxz(i,j,k)**2 &
                 - 2.0d0*gxy(i,j,k)*gxz(i,j,k)*gyz(i,j,k) &
                 + gxx(i,j,k)*gyz(i,j,k)**2 &
                 + gzz(i,j,k)*gxy(i,j,k)**2 &
                 - gxx(i,j,k)*gyy(i,j,k)*gzz(i,j,k)))
!==============>>>>>>> ===================================================
    if (CCTK_EQUALS(parallel, "no")) then
         vel(i,j,k,1)= v_inf_x
         vel(i,j,k,2)= v_inf_y
         vel(i,j,k,3)= v_inf_z
      
    else if (CCTK_EQUALS(parallel, "yes")) then
         vel(i,j,k,1)= v_inf_x
         vel(i,j,k,2)= v_inf_y
         vel(i,j,k,3)= v_inf_z
    end if

         w_lorentz(i,j,k) = 1.0d0/sqrt(1.0d0-(&
                           + gxx(i,j,k) * vel(i,j,k,1) * vel(i,j,k,1)&
                           + gyy(i,j,k) * vel(i,j,k,2) * vel(i,j,k,2)&
                           + gzz(i,j,k) * vel(i,j,k,3) * vel(i,j,k,3)&
                           + 2.0d0*gxy(i,j,k) * vel(i,j,k,1) * vel(i,j,k,2)&
                           + 2.0d0*gxz(i,j,k) * vel(i,j,k,1) * vel(i,j,k,3)&
                           + 2.0d0*gyz(i,j,k) * vel(i,j,k,2) * vel(i,j,k,3)))  

   end do
   end do
   end do
       close(1)
       eps= press/(( mygamma - one )*rho)

    if CCTK_EQUALS(initial_Bvec,"penner") then

       do k = 1, cctk_lsh(3)
       do j = 1, cctk_lsh(2)
       do i = 1, cctk_lsh(1)
                
            r1 = sqrt(x(i,j,k)**2 + y(i,j,k)**2 + z(i,j,k)**2)

            Bvecx(i,j,k) = x(i,j,k)*((x(i,j,k)*y(i,j,k)*(-(r1**4*(x(i,j,k)**2 &
		+ y(i,j,k)**2)**2) + 2*my_a**2*r1**2*(x(i,j,k)**2 &
		+ y(i,j,k)**2)*((-1 + r1**2)*x(i,j,k)**2 &
		+ (-1 + 2*my_M*r1)*y(i,j,k)**2) + my_a**6*(x(i,j,k)**4 &
		+ 2*x(i,j,k)**2*y(i,j,k)**2) + my_a**4*((-1 + 2*my_M*r1 &
		+ 3*r1**2)*x(i,j,k)**4 + 2*(-1 + 2*my_M*r1 & 
		+ 2*r1**2)*x(i,j,k)**2*y(i,j,k)**2 &
		+ (-1 + 4*my_M*r1)*y(i,j,k)**4)))/(r1*(x(i,j,k)**2 &
		+ y(i,j,k)**2)*(my_a**2*y(i,j,k)**2 + r1**2*(x(i,j,k)**2 &
		+ y(i,j,k)**2))**2) + (-(r1**3*(-1 + 2*r1**2)*(x(i,j,k)**2 & 
		+ y(i,j,k)**2)**2) - my_a**2*r1*(x(i,j,k)**2 &
		+ y(i,j,k)**2)*((-1 + my_M*r1 + 2*r1**2)*x(i,j,k)**2 & 
		+ (-1 + 2*my_M*r1 + 4*r1**2)*y(i,j,k)**2) &
		+ my_a**4*(my_M*y(i,j,k)**2*(x(i,j,k)**2 + 2*y(i,j,k)**2) & 
		- r1*(x(i,j,k)**4 + 2*x(i,j,k)**2*y(i,j,k)**2 & 
		+ 2*y(i,j,k)**4)))/(Sqrt(1 &
		+ x(i,j,k)**2/y(i,j,k)**2)*(my_a**2*y(i,j,k)**2 & 
		+ r1**2*(x(i,j,k)**2 + y(i,j,k)**2))**2) & 
		- (my_a*(1 + (2*my_M*r1*(-my_a**2 &
		+ r1**2))/(r1**2 + (my_a**2*y(i,j,k)**2)/(x(i,j,k)**2 &
		+ y(i,j,k)**2))**2)*z(i,j,k))/((1 &
		+ x(i,j,k)**2/y(i,j,k)**2)**1.5*y(i,j,k)))*B_0





            Bvecy(i,j,k) = x(i,j,k)*((y(i,j,k)**2*(-(r1**4*(x(i,j,k)**2 &
		  + y(i,j,k)**2)**2) + 2*my_a**2*r1**2*(x(i,j,k)**2 &
		  + y(i,j,k)**2)*((-1 + r1**2)*x(i,j,k)**2 &
		  + (-1 + 2*my_M*r1)*y(i,j,k)**2) + my_a**6*(x(i,j,k)**4 &
		  + 2*x(i,j,k)**2*y(i,j,k)**2) +  my_a**4*((-1 + 2*my_M*r1 &
		  + 3*r1**2)*x(i,j,k)**4 + 2*(-1 + 2*my_M*r1 &
		  + 2*r1**2)*x(i,j,k)**2*y(i,j,k)**2 &
		  + (-1 + 4*my_M*r1)*y(i,j,k)**4)))/(r1*(x(i,j,k)**2 &
		  + y(i,j,k)**2)*(my_a**2*y(i,j,k)**2 + r1**2*(x(i,j,k)**2 &
		  + y(i,j,k)**2))**2) + (x(i,j,k)*(r1**3*(-1 &
		  + 2*r1**2)*(x(i,j,k)**2 + y(i,j,k)**2)**2 &
		  + my_a**2*r1*(x(i,j,k)**2 + y(i,j,k)**2)*((-1 + my_M*r1 &
		  + 2*r1**2)*x(i,j,k)**2 + (-1 + 2*my_M*r1 + 4*r1**2)*y(i,j,k)**2) &
		  + my_a**4*(-(my_M*y(i,j,k)**2*(x(i,j,k)**2 + 2*y(i,j,k)**2)) &
		  + r1*(x(i,j,k)**4 + 2*x(i,j,k)**2*y(i,j,k)**2 &
		  + 2*y(i,j,k)**4))))/(Sqrt(1 &
		  + x(i,j,k)**2/y(i,j,k)**2)*y(i,j,k)*(my_a**2*y(i,j,k)**2 &
		  + r1**2*(x(i,j,k)**2 + y(i,j,k)**2))**2) &
		  - (my_a*x(i,j,k)*(1 + (2*my_M*r1*(-my_a**2 &
		  + r1**2))/(r1**2 + (my_a**2*y(i,j,k)**2)/(x(i,j,k)**2 &
		  + y(i,j,k)**2))**2)*z(i,j,k))/((1 &
		  + x(i,j,k)**2/y(i,j,k)**2)**1.5*y(i,j,k)**2))*B_0




           Bvecz(i,j,k) = (x(i,j,k)*((my_a*x(i,j,k)*(1 + (2*my_M*r1*(-my_a**2 &
		+ r1**2))/(r1**2 + (my_a**2*y(i,j,k)**2)/(x(i,j,k)**2 & 
		+ y(i,j,k)**2))**2))/Sqrt(1 + x(i,j,k)**2/y(i,j,k)**2) &
		+ (y(i,j,k)**2*(-(r1**4*(x(i,j,k)**2 + y(i,j,k)**2)**2) &
		+ 2*my_a**2*r1**2*(x(i,j,k)**2 + y(i,j,k)**2)*((-1 &
                + r1**2)*x(i,j,k)**2 & 
		+ (-1 + 2*my_M*r1)*y(i,j,k)**2) + my_a**6*(x(i,j,k)**4 &
		+ 2*x(i,j,k)**2*y(i,j,k)**2) &
		+ my_a**4*((-1 + 2*my_M*r1 + 3*r1**2)*x(i,j,k)**4 &
		+ 2*(-1 + 2*my_M*r1 + 2*r1**2)*x(i,j,k)**2*y(i,j,k)**2 &
		+ (-1 + 4*my_M*r1)*y(i,j,k)**4))*z(i,j,k))/(r1*(x(i,j,k)**2 &
		+ y(i,j,k)**2)*(my_a**2*y(i,j,k)**2 + r1**2*(x(i,j,k)**2 &
		+ y(i,j,k)**2))**2))*B_0)/y(i,j,k)






                  det_gamma = -(gxz(i,j,k)**2*gyy(i,j,k)) + 2.0d0*gxy(i,j,k)*gxz(i,j,k)*gyz(i,j,k) &
                          -gxx(i,j,k)*gyz(i,j,k)**2 &
                          -gxy(i,j,k)**2*gzz(i,j,k) + gxx(i,j,k)*gyy(i,j,k)*gzz(i,j,k)


!              call Prim2ConGenM(GRHydro_eos_handle,gxx(i,j,k),gxy(i,j,k),&
!                   gxz(i,j,k),gyy(i,j,k),gyz(i,j,k),gzz(i,j,k),&
!                   det, dens(i,j,k),sx(i,j,k),sy(i,j,k),sz(i,j,k),&
!                   tau(i,j,k),Bconsx(i,j,k),Bconsy(i,j,k),Bconsz(i,j,k),rho(i,j,k),&
!                   velx(i,j,k),vely(i,j,k),velz(i,j,k),&
!                   eps(i,j,k),press(i,j,k),Bvecx(i,j,k),Bvecy(i,j,k),Bvecz(i,j,k),&
!                   w_lorentz(i,j,k))


        end do
        end do
        end do


    else if CCTK_EQUALS(initial_Bvec,"poloidal_field") then

!   Poloidal Magnetic field implemented as in "General relativistic 
!   simulations of magnetized binary neutron star mergers" - 
!    by Yuk Tung Liu, Stuart L. Shapiro, Zachariah B. Etienne, and 
!   Keisuke Taniguchi - Phys. Rev. D 78, 024012 (2008)  
!        Bvec=0.0d0
       do k = 1, cctk_lsh(3) 
       do j = 1, cctk_lsh(2) 
       do i = 1, cctk_lsh(1) 


               r1 = sqrt(x(i,j,k)**2 + y(i,j,k)**2 + z(i,j,k)**2)
    

              det = sqrt(- gxz(i,j,k)**2*gyy(i,j,k) + 2.0d0*gxy(i,j,k)*gxz(i,j,k)*gyz(i,j,k) &
                      - gxx(i,j,k)*gyz(i,j,k)**2 - gxy(i,j,k)**2*gzz(i,j,k) &
                     + gxx(i,j,k)*gyy(i,j,k)*gzz(i,j,k))
          

                Bvecx(i,j,k) = 3.0d0*A_b*x(i,j,k)*z(i,j,k)/r1**5/det            
                Bvecy(i,j,k) = 3.0d0*A_b*y(i,j,k)*z(i,j,k)/r1**5/det            
                Bvecz(i,j,k) = 3.0d0*A_b*z(i,j,k)/r1**5 - 3.0d0*A_b/r1**3/det        


      end do
      end do
      end do

    else if CCTK_EQUALS(initial_Bvec,"constant_field") then

       do k = 1, cctk_lsh(3) 
       do j = 1, cctk_lsh(2) 
       do i = 1, cctk_lsh(1) 


                Bvecx(i,j,k) = 0.0d0            
                Bvecy(i,j,k) = 0.0d0            
                Bvecz(i,j,k) = A_b        



               

                  
         end do
         end do
         end do          
         print *,  "Yeah....It worked"
                 

    end if 



 
  end subroutine Fluid_Wind
