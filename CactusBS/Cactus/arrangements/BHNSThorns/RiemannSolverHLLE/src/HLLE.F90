 /*@@
   @file      Whisky_HLLE.F90
   @date      Sat Jan 26 01:40:14 2002
   @author    Ian Hawke, Pedro Montero, Toni Font
   @desc 
   The HLLE solver. Called from the wrapper function, so works in 
   all directions.
   @enddesc 
 @@*/
   
#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

#include "SpaceMask.h"

 /*@@
   @routine    Whisky_HLLE
   @date       Sat Jan 26 01:41:02 2002
   @author     Ian Hawke, Pedro Montero, Toni Font
   @desc 
   The HLLE solver. Sufficiently simple that its just one big routine.
   @enddesc 
   @calls     
   @calledby   
   @history 
   Altered from Cactus 3 routines originally written by Toni Font.
   @endhistory 

@@*/

subroutine Whisky_HLLE(CCTK_ARGUMENTS)

  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  
  integer :: i, j, k, m
  CCTK_REAL, dimension(5) :: consp,consm_i,fplus,fminus,lamplus
  CCTK_REAL, dimension(5) :: f1,lamminus
  CCTK_REAL, dimension(5) :: qdiff
  CCTK_REAL ::  charmin, charmax, charpm,avg_alp,avg_det
  CCTK_REAL :: gxxh, gxyh, gxzh, gyyh, gyzh, gzzh, uxxh, uxyh, &
       uxzh, uyyh, uyzh, uzzh, avg_beta, usendh, alp_l, alp_r
    
  CCTK_INT :: type_bits, trivial, not_trivial

  if (flux_direction == 1) then
    call SpaceMask_GetTypeBits(type_bits, "Hydro_RiemannProblemX")
    call SpaceMask_GetStateBits(trivial, "Hydro_RiemannProblemX", &
         &"trivial")
  else if (flux_direction == 2) then
    call SpaceMask_GetTypeBits(type_bits, "Hydro_RiemannProblemY")
    call SpaceMask_GetStateBits(trivial, "Hydro_RiemannProblemY", &
         &"trivial")
  else if (flux_direction == 3) then
    call SpaceMask_GetTypeBits(type_bits, "Hydro_RiemannProblemZ")
    call SpaceMask_GetStateBits(trivial, "Hydro_RiemannProblemZ", &
         &"trivial")
  else
    call CCTK_WARN(0, "Flux direction not x,y,z")
  end if

  !$OMP PARALLEL DO &
  !$OMP PRIVATE(i,j,k,m, avg_alp, avg_det, charmin, charmax, charpm,&
  !$OMP         gxxh, gxyh, gxzh, gyyh, gyzh, gzzh, uxxh, uxyh, uxzh,&
  !$OMP         uyyh, uyzh, uzzh, avg_beta, usendh, alp_l, alp_r,&
  !$OMP         consp,consm_i,fplus,fminus,lamplus, f1,lamminus, qdiff) 
  do k = whisky_stencil, cctk_lsh(3) - whisky_stencil
    do j = whisky_stencil, cctk_lsh(2) - whisky_stencil
      do i = whisky_stencil, cctk_lsh(1) - whisky_stencil
        
        f1 = 0.d0
        lamminus = 0.d0
        lamplus = 0.d0
        consp = 0.d0
        consm_i = 0.d0
        fplus = 0.d0
        fminus = 0.d0
        qdiff = 0.d0
        
!!$        Set the left (p for plus) and right (m_i for minus, i+1) states
        
        consp(1)   = densplus(i,j,k) 
        consp(2)   = sxplus(i,j,k)
        consp(3)   = syplus(i,j,k)
        consp(4)   = szplus(i,j,k)
        consp(5)   = tauplus(i,j,k)
        
        consm_i(1) = densminus(i+xoffset,j+yoffset,k+zoffset)
        consm_i(2) = sxminus(i+xoffset,j+yoffset,k+zoffset)
        consm_i(3) = syminus(i+xoffset,j+yoffset,k+zoffset)
        consm_i(4) = szminus(i+xoffset,j+yoffset,k+zoffset)
        consm_i(5) = tauminus(i+xoffset,j+yoffset,k+zoffset) 
          
!!$        Calculate various metric terms here.
!!$        Note also need the average of the lapse at the 
!!$        left and right points.

        if (flux_direction == 1) then
          avg_beta = 0.5d0 * (betax(i+xoffset,j+yoffset,k+zoffset) + &
               betax(i,j,k))
        else if (flux_direction == 2) then
          avg_beta = 0.5d0 * (betay(i+xoffset,j+yoffset,k+zoffset) + &
               betay(i,j,k))
        else if (flux_direction == 3) then
          avg_beta = 0.5d0 * (betaz(i+xoffset,j+yoffset,k+zoffset) + &
               betaz(i,j,k))
        else
          call CCTK_WARN(0, "Flux direction not x,y,z")
        end if

        avg_alp = 0.5 * (alp(i,j,k) + alp(i+xoffset,j+yoffset,k+zoffset))
        
        gxxh = 0.5d0 * (gxx(i+xoffset,j+yoffset,k+zoffset) + gxx(i,j,k))
        gxyh = 0.5d0 * (gxy(i+xoffset,j+yoffset,k+zoffset) + gxy(i,j,k))
        gxzh = 0.5d0 * (gxz(i+xoffset,j+yoffset,k+zoffset) + gxz(i,j,k))
        gyyh = 0.5d0 * (gyy(i+xoffset,j+yoffset,k+zoffset) + gyy(i,j,k))
        gyzh = 0.5d0 * (gyz(i+xoffset,j+yoffset,k+zoffset) + gyz(i,j,k))
        gzzh = 0.5d0 * (gzz(i+xoffset,j+yoffset,k+zoffset) + gzz(i,j,k))

        call SpatialDeterminant(gxxh,gxyh,gxzh,&
             gyyh,gyzh,gzzh,avg_det)

!!$ If the Riemann problem is trivial, just calculate the fluxes from the 
!!$ left state and skip to the next cell
          
        if (SpaceMask_CheckStateBitsF90(space_mask, i, j, k, type_bits, trivial)) then

          if (flux_direction == 1) then
            call num_x_flux(consp(1),consp(2),consp(3),consp(4),consp(5),&
                 f1(1),f1(2),f1(3),&
                 f1(4),f1(5),&
                 velxplus(i,j,k),pressplus(i,j,k),&
                 avg_det,avg_alp,avg_beta)
          else if (flux_direction == 2) then
            call num_x_flux(consp(1),consp(3),consp(4),consp(2),consp(5),&
                 f1(1),f1(3),f1(4),&
                 f1(2),f1(5),&
                 velyplus(i,j,k),pressplus(i,j,k),&
                 avg_det,avg_alp,avg_beta)
          else if (flux_direction == 3) then
            call num_x_flux(consp(1),consp(4),consp(2),consp(3),consp(5),&
                 f1(1),f1(4),f1(2),&
                 f1(3),f1(5),&
                 velzplus(i,j,k),pressplus(i,j,k),&
                 avg_det,avg_alp,avg_beta)
          else
            call CCTK_WARN(0, "Flux direction not x,y,z")
          end if
          
        else !!! The end of this branch is right at the bottom of the routine
            
          call UpperMetric(uxxh, uxyh, uxzh, uyyh, uyzh, uzzh, &
               avg_det,gxxh, gxyh, gxzh, &
               gyyh, gyzh, gzzh)
          
          if (flux_direction == 1) then
            usendh = uxxh
          else if (flux_direction == 2) then
            usendh = uyyh
          else if (flux_direction == 3) then
            usendh = uzzh
          else
            call CCTK_WARN(0, "Flux direction not x,y,z")
          end if
          
!!$        Calculate the jumps in the conserved variables
          
          qdiff(1) = consm_i(1) - consp(1)
          qdiff(2) = consm_i(2) - consp(2)
          qdiff(3) = consm_i(3) - consp(3)
          qdiff(4) = consm_i(4) - consp(4)
          qdiff(5) = consm_i(5) - consp(5)
          
!!$        Eigenvalues and fluxes either side of the cell interface
          
          if (flux_direction == 1) then
            call eigenvalues_hlle(whisky_eos_handle,&
                 rhominus(i+xoffset,j+yoffset,k+zoffset),&
                 velxminus(i+xoffset,j+yoffset,k+zoffset),&
                 velyminus(i+xoffset,j+yoffset,k+zoffset),&
                 velzminus(i+xoffset,j+yoffset,k+zoffset),&
                 epsminus(i+xoffset,j+yoffset,k+zoffset),&
                 w_lorentzminus(i+xoffset,j+yoffset,k+zoffset),&
                 lamminus,gxxh,gxyh,gxzh,gyyh,&
                 gyzh,gzzh,&
                 usendh,avg_det,avg_alp,avg_beta)
            call eigenvalues_hlle(whisky_eos_handle,rhoplus(i,j,k),&
                 velxplus(i,j,k),velyplus(i,j,k),&
                 velzplus(i,j,k),epsplus(i,j,k),w_lorentzplus(i,j,k),&
                 lamplus,gxxh,gxyh,gxzh,gyyh,&
                 gyzh,gzzh,&
                 usendh,avg_det,avg_alp,avg_beta)
            call num_x_flux(consp(1),consp(2),consp(3),consp(4),consp(5),&
                 fplus(1),fplus(2),fplus(3),fplus(4),&
                 fplus(5),velxplus(i,j,k),pressplus(i,j,k),&
                 avg_det,avg_alp,avg_beta)
            call num_x_flux(consm_i(1),consm_i(2),consm_i(3),&
                 consm_i(4),consm_i(5),fminus(1),fminus(2),fminus(3),&
                 fminus(4), fminus(5),&
                 velxminus(i+xoffset,j+yoffset,k+zoffset),&
                 pressminus(i+xoffset,j+yoffset,k+zoffset),&
                 avg_det,avg_alp,avg_beta)
          else if (flux_direction == 2) then
            call eigenvalues_hlle(whisky_eos_handle,&
                 rhominus(i+xoffset,j+yoffset,k+zoffset),&
                 velyminus(i+xoffset,j+yoffset,k+zoffset),&
                 velzminus(i+xoffset,j+yoffset,k+zoffset),&
                 velxminus(i+xoffset,j+yoffset,k+zoffset),&
                 epsminus(i+xoffset,j+yoffset,k+zoffset),&
                 w_lorentzminus(i+xoffset,j+yoffset,k+zoffset),&
                 lamminus,gyyh,gyzh,gxyh,gzzh,&
                 gxzh,gxxh,&
                 usendh,avg_det,avg_alp,avg_beta)
            call eigenvalues_hlle(whisky_eos_handle,rhoplus(i,j,k),&
                 velyplus(i,j,k),velzplus(i,j,k),&
                 velxplus(i,j,k),epsplus(i,j,k),w_lorentzplus(i,j,k),&
                 lamplus,gyyh,gyzh,gxyh,gzzh,&
                 gxzh,gxxh,&
                 usendh,avg_det,avg_alp,avg_beta)
            call num_x_flux(consp(1),consp(3),consp(4),consp(2),consp(5),&
                 fplus(1),fplus(3),fplus(4),fplus(2),&
                 fplus(5),velyplus(i,j,k),pressplus(i,j,k),&
                 avg_det,avg_alp,avg_beta)
            call num_x_flux(consm_i(1),consm_i(3),consm_i(4),&
                 consm_i(2),consm_i(5),fminus(1),fminus(3),fminus(4),&
                 fminus(2), fminus(5),&
                 velyminus(i+xoffset,j+yoffset,k+zoffset),&
                 pressminus(i+xoffset,j+yoffset,k+zoffset),&
                 avg_det,avg_alp,avg_beta)
          else if (flux_direction == 3) then
            call eigenvalues_hlle(whisky_eos_handle,&
                 rhominus(i+xoffset,j+yoffset,k+zoffset),&
                 velzminus(i+xoffset,j+yoffset,k+zoffset),&
                 velxminus(i+xoffset,j+yoffset,k+zoffset),&
                 velyminus(i+xoffset,j+yoffset,k+zoffset),&
                 epsminus(i+xoffset,j+yoffset,k+zoffset),&
                 w_lorentzminus(i+xoffset,j+yoffset,k+zoffset),&
                 lamminus,gzzh,gxzh,gyzh,&
                 gxxh,gxyh,gyyh,&
                 usendh,avg_det,avg_alp,avg_beta)
            call eigenvalues_hlle(whisky_eos_handle,rhoplus(i,j,k),&
                 velzplus(i,j,k),velxplus(i,j,k),&
                 velyplus(i,j,k),epsplus(i,j,k),w_lorentzplus(i,j,k),&
                 lamplus,gzzh,gxzh,gyzh,&
                 gxxh,gxyh,gyyh,&
                 usendh,avg_det,avg_alp,avg_beta)
            call num_x_flux(consp(1),consp(4),consp(2),consp(3),consp(5),&
                 fplus(1),fplus(4),fplus(2),fplus(3),&
                 fplus(5),velzplus(i,j,k),pressplus(i,j,k),avg_det,&
                 avg_alp,avg_beta)
            call num_x_flux(consm_i(1),consm_i(4),consm_i(2),&
                 consm_i(3),consm_i(5),fminus(1),fminus(4),fminus(2),&
                 fminus(3), fminus(5),&
                 velzminus(i+xoffset,j+yoffset,k+zoffset),&
                 pressminus(i+xoffset,j+yoffset,k+zoffset),&
                 avg_det,avg_alp,avg_beta)
          else
            call CCTK_WARN(0, "Flux direction not x,y,z")
          end if
          
!!$        Find minimum and maximum wavespeeds
      
          charmin = min(0.d0, lamplus(1), lamplus(2), lamplus(3), &
               lamplus(4),lamplus(5),  lamminus(1),lamminus(2),lamminus(3),&
               lamminus(4),lamminus(5))  
          
          charmax = max(0.d0, lamplus(1), lamplus(2), lamplus(3), &
               lamplus(4),lamplus(5),  lamminus(1),lamminus(2),lamminus(3),&
               lamminus(4),lamminus(5))
          
          charpm = charmax - charmin
          
!!$        Calculate flux by standard formula
          
          do m = 1,5 
            
            qdiff(m) = consm_i(m) - consp(m)
            
            f1(m) = (charmax * fplus(m) - charmin * fminus(m) + &
                 charmax * charmin * qdiff(m)) / charpm
            
          end do

        end if !!! The end of the SpaceMask check for a trivial RP.
            
        densflux(i, j, k) = f1(1)
        sxflux(i, j, k) = f1(2)
        syflux(i, j, k) = f1(3)
        szflux(i, j, k) = f1(4)
        tauflux(i, j, k) = f1(5)

      end do
    end do
  end do
  !$OMP END PARALLEL DO 
end subroutine Whisky_HLLE
