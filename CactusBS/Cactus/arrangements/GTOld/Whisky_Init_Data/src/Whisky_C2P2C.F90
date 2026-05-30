 /*@@
   @file      Whisky_C2P2C.F90
   @date      Sat Jan 26 02:44:43 2002
   @author    Luca Baiotti
   @desc 
   A test of the conservative <--> primitive variable exchange
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

 /*@@
   @routine    c2p2c
   @date       Sat Jan 26 02:45:19 2002
   @author     Luca Baiotti
   @desc 
   Testing the conservative <--> primitive variable transformations.
   The values before and after should match.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

subroutine c2p2c(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer didit,i,j,k,nx,ny,nz,istat,st,ierr
  integer :: strlength=40, maxstrlength=200
  character(len=200) :: fname, dirname
  CCTK_REAL det,uxx,uxy,uxz,uyy,uyz,uzz
  CCTK_REAL gxx_loc,gxy_loc,gxz_loc,gyy_loc,gyz_loc,gzz_loc
  CCTK_REAL dens_send,sx_send,sy_send,sz_send,tau_send
  CCTK_REAL rho_send,velx_send,vely_send,velz_send,eps_send
  CCTK_REAL press_send,w_lorentz_send,x_send,y_send,z_send,r_send
  CCTK_REAL bnx_send,bny_send,bnz_send,bvecx_send,bvecy_send,bvecz_send
  CCTK_REAL weight
  logical epsnegative

  call CCTK_WARN(1,"This test works only with Ideal_Fluid EoS")

  nx = cctk_lsh(1)
  ny = cctk_lsh(2)
  nz = cctk_lsh(3)
  
  x_send = 0.0d0
  y_send = 0.0d0
  z_send = 0.0d0
  r_send = 0.0d0
  
  gxx_loc = 1.0d0
  gyy_loc = 1.0d0 
  gzz_loc = 1.0d0
  gxy_loc = 0.0d0
  gxz_loc = 0.0d0
  gyz_loc = 0.0d0

  weight = 1.0d0
  
  det = 1.0d0
  
  uxx = 1.0d0
  uyy = 1.0d0 
  uzz = 1.0d0
  uxy = 0.0d0
  uxz = 0.0d0
  uyz = 0.0d0

  if (whisky_mhd_handle.gt.1) then
     bnx_send  = 0.2d0
     bny_send  = 1.3d0
     bnz_send  = 0.1d0
     dens_send = 2.17075d0
     sx_send   = 0.76701d0
     sy_send   = 0.185505d0
     sz_send   = 6.38358d0
     tau_send  = 5.06554d0
     eps_send  = 1.0d-6           
  else !! HD + Poisoned Bfield
     dens_send = 1.29047362d0
     sx_send   = 0.166666658d0
     sy_send   = 0.166666658d0
     sz_send   = 0.166666658d0
     tau_send  = 0.484123939d0
     bnx_send  = -100d0
     bny_send  = -100d0
     bnz_send  = -100d0
     eps_send = 1.0d-6
  end if
  
  press_send = 6.666666666666667d-7
  w_lorentz_send = 1.0d0
  if (whisky_mhd_handle.gt.1) then
     bvecx_send = 0.0d0
     bvecy_send = 0.0d0
     bvecz_send = 0.0d0
  else !! Poison
     bvecx_send = -100.0d0
     bvecy_send = -100.0d0
     bvecz_send = -100.0d0
  end if

  epsnegative = .false.
  
  whisky_rho_min = 1.e-10

  call CCTK_FortranString(strlength,out_dir,dirname)
  if (strlength .gt. maxstrlength) then
     call CCTK_WARN(0,"The output directory string is too long for current Whisky_Init_Data output.")
  end if
  write(fname,*) trim(dirname),'/Con2Prim2Con.info' !! Use trim to remove trailing blanks
  write(*,*) 'Writing Con2Prim2Con output to file: ',fname
  !! BEWARE: CCTK_FortranString also has leading whitespace that needs to be removed via adjustl
  open(unit=10,file=adjustl(fname),status='replace',iostat=st)
  if ( st/=0) then
     write(*,*) 'Unable to open file ',fname
     call CCTK_WARN(0,"Unable to open file for output.")
  end if

  write(10,*) 'C2P2C test: initial values.'
  write(10,*) '   conservative variables: '
  write(10,*) '   dens: ',dens_send
  write(10,*) '   sx  : ',sx_send
  write(10,*) '   sy  : ',sy_send
  write(10,*) '   sz  : ',sz_send
  write(10,*) '   tau : ',tau_send
  if ( whisky_mhd_handle.gt.1 ) then
     write(10,*) '   Bnx   : ',bnx_send
     write(10,*) '   Bny   : ',bny_send
     write(10,*) '   Bnz   : ',bnz_send
  end if
  write(10,*) '   eps : ',eps_send, ' (For initial guess) '
  
  write(10,*) 'C2P2C test: getting the associated primitive variables.'

  if (whisky_mhd_handle.gt.1) then
     ierr = Con2PrimMHD(whisky_eos_handle,whisky_polytrope_handle,0,& !! 0 signifies not a polytrope
          cctk_iteration,1,4,                                       & !! 1 signifies GST_Newton, 4 hv^2 solver
          dens_send,tau_send,sx_send,sy_send,sz_send,bnx_send,bny_send,bnz_send,&
          rho_send,eps_send,velx_send,vely_send,velz_send,bvecx_send,bvecy_send,&
          bvecz_send,w_lorentz_send,press_send,x_send,y_send,z_send,&
          gxx_loc,gxy_loc,gxz_loc,gyy_loc,gyz_loc,gzz_loc,whisky_init_data_reflevel,weight,&
          0,whisky_rho_min,whisky_tau_min)
  else
     call Con2Prim_pt(whisky_eos_handle,dens_send,sx_send,sy_send,sz_send, &
          tau_send,rho_send,velx_send,vely_send,velz_send, &
          eps_send,press_send,w_lorentz_send, &
          uxx,uxy,uxz,uyy,uyz,uzz,det,x_send,y_send,z_send,r_send,&
          epsnegative,whisky_rho_min, whisky_init_data_reflevel,weight)
  end if

  write(10,*) 'C2P2C test: the primitive variables are'
  write(10,*) '   primitive variables: '
  write(10,*) '   rho   : ',rho_send
  write(10,*) '   velx  : ',velx_send
  write(10,*) '   vely  : ',vely_send
  write(10,*) '   velz  : ',velz_send
  write(10,*) '   press : ',press_send
  if ( whisky_mhd_handle.gt.1 ) then
     write(10,*) '   Bx  : ',bvecx_send
     write(10,*) '   By  : ',bvecy_send
     write(10,*) '   Bz  : ',bvecz_send
  end if 
  write(10,*) '   eps   : ',eps_send
  write(10,*) '   W     : ',w_lorentz_send
  
  write(10,*) 'C2P2C test: converting back to conserved variables.'
  if (whisky_mhd_handle.gt.1) then
     call Prim2ConMHD(whisky_eos_handle,gxx_loc, gxy_loc, gxz_loc, gyy_loc, gyz_loc, gzz_loc, det, &
          dens_send, sx_send, sy_send, sz_send, tau_send, rho_send, &
          velx_send, vely_send, velz_send, eps_send, press_send, w_lorentz_send, &
          bnx_send, bny_send, bnz_send, bvecx_send, bvecy_send, bvecz_send)
  else
     call Prim2ConGen(whisky_eos_handle,gxx_loc, gxy_loc, gxz_loc, gyy_loc, gyz_loc, gzz_loc, det, &
          dens_send, sx_send, sy_send, sz_send, tau_send, rho_send, &
          velx_send, vely_send, velz_send, eps_send, press_send, w_lorentz_send) 
  end if
  
  write(10,*) 'C2P2C test: the conserved variables are'
  write(10,*) '   conservative variables: '
  write(10,*) '   dens: ',dens_send
  write(10,*) '   sx  : ',sx_send
  write(10,*) '   sy  : ',sy_send
  write(10,*) '   sz  : ',sz_send
  write(10,*) '   tau : ',tau_send
  if ( whisky_mhd_handle.gt.1 ) then
     write(10,*) '   Bnx   : ',bnx_send
     write(10,*) '   Bny   : ',bny_send
     write(10,*) '   Bnz   : ',bnz_send
  end if
  close(10)

  !!STOP
  write(*,*) 'C2P2C test: We tested the Conservatives -> Primitives -> Conservatives conversion.  Now exiting.  Is this correct?'
  call CCTK_Exit(istat , cctkGH, 0 )

  return

end subroutine c2p2c
