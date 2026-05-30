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

! --->   NUMBERS   <---

  zero = 0.0d0
  one  = 1.0d0
  half = 0.5d0
  two  = 2.0d0
  pi  = 4.0D0*atan(1.0D0)
 

 
! --->   INITIAL DATA   <---
   do k = 1, cctk_lsh(3)
   do j = 1, cctk_lsh(2)
   do i = 1, cctk_lsh(1)

  !! vel(i,j,k,1) is equivalent to vel^{x}       
        
         rho(i,j,k) = rho_ini_cons

         eps(i,j,k) = A_b * rho(i,j,k)**(mygamma - one) / (mygamma - one)
! A_b = polytropic constant usually 100
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
       press = (mygamma -one) * (rho * eps)
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
