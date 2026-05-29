 /*@@
   @file      Whisky_Con2Prim.F90
   @date      Sat Jan 26 02:49:32 2002
   @author    Luca Baiotti
   @desc 
   A test of the conservative to primitive exchange.
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

 /*@@
   @routine    Whisky_con2primtest
   @date       Sat Jan 26 02:49:58 2002
   @author     Luca Baiotti
   @desc 
   A test of the conservative to primitive variable solver.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

subroutine Whisky_Init_Data_RefinementLevel(CCTK_ARGUMENTS)

  implicit none
  DECLARE_CCTK_ARGUMENTS

  whisky_init_data_reflevel = aint(log10(dble(cctk_levfac(1)))/log10(2.0d0))

end subroutine Whisky_Init_Data_RefinementLevel


subroutine Whisky_con2primtest(CCTK_ARGUMENTS)
  
  implicit  none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer didit,i,j,k,nx,ny,nz,istat,st,ierr
  integer :: strlength=40, maxstrlength=200
  character(len=200) :: fname, dirname
  CCTK_REAL det,uxx,uxy,uxz,uyy,uyz,uzz
  CCTK_REAL dens_send,sx_send,sy_send,sz_send,tau_send
  CCTK_REAL rho_send,velx_send,vely_send,velz_send,eps_send
  CCTK_REAL press_send,w_lorentz_send,x_send,y_send,z_send,r_send
  CCTK_REAL bnx_send, bny_send, bnz_send
  CCTK_REAL bvecx_send, bvecy_send, bvecz_send
  CCTK_REAL weight
  logical epsnegative
  
  call CCTK_WARN(1,"For this test, remember to use a polytropic EoS and to set eos_gamma = 2.0 and eos_k = 100.0")

  nx = cctk_lsh(1)
  ny = cctk_lsh(2)
  nz = cctk_lsh(3)
  
  x_send = 0.0d0
  y_send = 0.0d0
  z_send = 0.0d0
  r_send = 0.0d0
  
  det = 1.0d0
  uxx = 1.0d0
  uyy = 1.0d0 
  uzz = 1.0d0
  uxy = 0.0d0
  uxz = 0.0d0
  uyz = 0.0d0
  weight = 1d0

  dens_send = 1.29047362d0
  sx_send   = 0.166666658d0
  sy_send   = 0.166666658d0
  sz_send   = 0.166666658d0
  tau_send  = 0.484123939d0
  if (whisky_mhd_handle.gt.1) then
     bnx_send  = 0.253289603d0
     bny_send  = 0.253289603d0
     bnz_send  = 0.253289603d0
  else !! Poison
     bnx_send  = -100d0
     bny_send  = -100d0
     bnz_send  = -100d0
  end if
  
  rho_send = 1.0d0
  velx_send = 0.0d0
  vely_send = 0.0d0
  velz_send = 0.0d0
  eps_send = 1.0d-6
  press_send = 6.666666666666667d-7
  w_lorentz_send = 1.0d0
  if ( whisky_mhd_handle.gt.1 ) then
     bvecx_send = 0.0d0
     bvecy_send = 0.0d0
     bvecz_send = 0.0d0
  else !! Poison
     bvecx_send = -100d0
     bvecy_send = -100d0
     bvecz_send = -100d0
  end if

  epsnegative = .false.

  call CCTK_FortranString(strlength,out_dir,dirname)
  if (strlength .gt. maxstrlength) then
     call CCTK_WARN(0,"The output directory string is too long for current Whisky_Init_Data output.")
  end if
  write(fname,*) trim(dirname),'/Con2Prim.info' !! Use trim to remove trailing blanks
  write(*,*) 'Writing Con2Prim output to file: ',fname
  !! BEWARE: CCTK_FortranString also has leading whitespace that needs to be removed via adjustl
  open(unit=10,file=adjustl(fname),status='replace',iostat=st)
  if ( st/=0) then
     write(*,*) 'Unable to open file ',fname
     call CCTK_WARN(0,"Unable to open file for output.")
  end if
  
  write(10,*) 'Con2Prim test: converting to primitive variables.'
  write(10,*) 'Con2Prim test: the conservative variables were'
  write(10,*) '   dens  : ',dens_send
  write(10,*) '   sx    : ',sx_send
  write(10,*) '   sy    : ',sy_send
  write(10,*) '   sz    : ',sz_send
  write(10,*) '   tau   : ',tau_send
  if ( whisky_mhd_handle.gt.1 ) then
     write(10,*) '   Bnx   : ',bnx_send
     write(10,*) '   Bny   : ',bny_send
     write(10,*) '   Bnz   : ',bnz_send
  end if
  write(10,*) '   eps   : ',eps_send, ' (Initial guess)'

  if ( whisky_mhd_handle.gt.1 ) then
     ierr = Con2PrimMHD(whisky_eos_handle,whisky_polytrope_handle,0,& !! 0 signifies not a polytrope
            cctk_iteration,1,4,                                     & !! 1 signifies GST_Newton, 4 hv^2 method
            dens_send,tau_send,sx_send,sy_send,sz_send,bnx_send,bny_send,bnz_send,&
            rho_send,eps_send,velx_send,vely_send,velz_send,bvecx_send,bvecy_send,&
            bvecz_send,w_lorentz_send,press_send,x_send,y_send,z_send,&
            uxx,uxy,uxz,uyy,uyz,uzz,whisky_init_data_reflevel,weight,&
            0,whisky_rho_min,whisky_tau_min)
  else
     call Con2Prim_pt(whisky_eos_handle,dens_send,sx_send,sy_send,sz_send,&
          tau_send,rho_send,velx_send,vely_send,velz_send,&
          eps_send,press_send,w_lorentz_send, &
          uxx,uxy,uxz,uyy,uyz,uzz,det,x_send,y_send,z_send,r_send,&
          epsnegative,whisky_rho_min, whisky_init_data_reflevel, weight)
  end if  

  write(10,*) 'Con2Prim test: the primitive variables are'
  write(10,*) '   primitive variables: '
  write(10,*) '   rho   : ',rho_send
  write(10,*) '   velx  : ',velx_send
  write(10,*) '   vely  : ',vely_send
  write(10,*) '   velz  : ',velz_send
  write(10,*) '   eps   : ',eps_send
  write(10,*) '   press : ',press_send
  write(10,*) '   w_lor : ',w_lorentz_send
  if ( whisky_mhd_handle.gt.1 ) then
     write(10,*) '   B_x   : ',bvecx_send
     write(10,*) '   B_y   : ',bvecy_send
     write(10,*) '   B_z   : ',bvecz_send
  end if
  close(10)

  !!STOP
  write(*,*) 'Con2Prim test: We tested the Conservatives -> Primitives conversion.  Now exiting.  Is this correct?'
  call CCTK_Exit(istat , cctkGH, 0 )
  STOP

  return

end subroutine Whisky_con2primtest
