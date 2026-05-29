 /*@@
   @file      Whisky_Tracer.F90
   @date      Sat Jan 26 01:40:14 2002
   @author    Ian Hawke, Pedro Montero, Toni Font
   @desc 
   The HLLE tracer routine from Whisky.
   @enddesc 
 @@*/
   
#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

#include "SpaceMask.h"

 /*@@
   @routine    Whisky_HLLE_Tracer
   @date       Mon Mar  8 13:47:13 2004
   @author     Ian Hawke
   @desc 
   HLLE just for the tracer.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

subroutine Whisky_HLLE_Tracer(CCTK_ARGUMENTS)

  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  
  integer :: i, j, k,  m
  CCTK_REAL, dimension(number_of_tracers) :: consp,consm_i,fplus,fminus,f1
  CCTK_REAL, dimension(5) :: lamminus,lamplus
  CCTK_REAL, dimension(number_of_tracers) :: qdiff
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
        
        do m=1,number_of_tracers
           consp(m)   = cons_tracerplus(i,j,k,m)         
           consm_i(m) = cons_tracerminus(i+xoffset,j+yoffset,k+zoffset,m)
        enddo
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
          
        qdiff = consm_i - consp
          
!!$        Eigenvalues and fluxes either side of the cell interface
          
        if (flux_direction == 1) then
          call eigenvalues(whisky_eos_handle,&
               rhominus(i+xoffset,j+yoffset,k+zoffset),&
               velxminus(i+xoffset,j+yoffset,k+zoffset),&
               velyminus(i+xoffset,j+yoffset,k+zoffset),&
               velzminus(i+xoffset,j+yoffset,k+zoffset),&
               epsminus(i+xoffset,j+yoffset,k+zoffset),&
               w_lorentzminus(i+xoffset,j+yoffset,k+zoffset),&
               lamminus,gxxh,gxyh,gxzh,gyyh,&
               gyzh,gzzh,&
               usendh,avg_det,avg_alp,avg_beta)
          call eigenvalues(whisky_eos_handle,rhoplus(i,j,k),&
               velxplus(i,j,k),velyplus(i,j,k),&
               velzplus(i,j,k),epsplus(i,j,k),w_lorentzplus(i,j,k),&
               lamplus,gxxh,gxyh,gxzh,gyyh,&
               gyzh,gzzh,&
               usendh,avg_det,avg_alp,avg_beta)
          fplus(:)  = (velxplus(i,j,k)  - avg_beta / avg_alp) * &
               cons_tracerplus(i,j,k,:)
          fminus(:) = (velxminus(i+xoffset,j+yoffset,k+zoffset) - avg_beta / avg_alp) * &
               cons_tracerminus(i+xoffset,j+yoffset,k+zoffset,:)
        else if (flux_direction == 2) then
          call eigenvalues(whisky_eos_handle,&
               rhominus(i+xoffset,j+yoffset,k+zoffset),&
               velyminus(i+xoffset,j+yoffset,k+zoffset),&
               velzminus(i+xoffset,j+yoffset,k+zoffset),&
               velxminus(i+xoffset,j+yoffset,k+zoffset),&
               epsminus(i+xoffset,j+yoffset,k+zoffset),&
               w_lorentzminus(i+xoffset,j+yoffset,k+zoffset),&
               lamminus,gyyh,gyzh,gxyh,gzzh,&
               gxzh,gxxh,&
               usendh,avg_det,avg_alp,avg_beta)
          call eigenvalues(whisky_eos_handle,rhoplus(i,j,k),&
               velyplus(i,j,k),velzplus(i,j,k),&
               velxplus(i,j,k),epsplus(i,j,k),w_lorentzplus(i,j,k),&
               lamplus,gyyh,gyzh,gxyh,gzzh,&
               gxzh,gxxh,&
               usendh,avg_det,avg_alp,avg_beta)
          fplus(:)  = (velyplus(i,j,k)  - avg_beta / avg_alp) * &
               cons_tracerplus(i,j,k,:)
          fminus(:) = (velyminus(i+xoffset,j+yoffset,k+zoffset) - avg_beta / avg_alp) * &
               cons_tracerminus(i+xoffset,j+yoffset,k+zoffset,:)
        else if (flux_direction == 3) then
          call eigenvalues(whisky_eos_handle,&
               rhominus(i+xoffset,j+yoffset,k+zoffset),&
               velzminus(i+xoffset,j+yoffset,k+zoffset),&
               velxminus(i+xoffset,j+yoffset,k+zoffset),&
               velyminus(i+xoffset,j+yoffset,k+zoffset),&
               epsminus(i+xoffset,j+yoffset,k+zoffset),&
               w_lorentzminus(i+xoffset,j+yoffset,k+zoffset),&
               lamminus,gzzh,gxzh,gyzh,&
               gxxh,gxyh,gyyh,&
               usendh,avg_det,avg_alp,avg_beta)
          call eigenvalues(whisky_eos_handle,rhoplus(i,j,k),&
               velzplus(i,j,k),velxplus(i,j,k),&
               velyplus(i,j,k),epsplus(i,j,k),w_lorentzplus(i,j,k),&
               lamplus,gzzh,gxzh,gyzh,&
               gxxh,gxyh,gyyh,&
               usendh,avg_det,avg_alp,avg_beta)
          fplus(:)  = (velzplus(i,j,k) - avg_beta / avg_alp) * &
               cons_tracerplus(i,j,k,:)
          fminus(:) = (velzminus(i+xoffset,j+yoffset,k+zoffset) - avg_beta / avg_alp) * &
               cons_tracerminus(i+xoffset,j+yoffset,k+zoffset,:)
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
          
        do m = 1,number_of_tracers
            
          qdiff(m) = consm_i(m) - consp(m)
          
          f1(m) = (charmax * fplus(m) - charmin * fminus(m) + &
               charmax * charmin * qdiff(m)) / charpm
          
        end do
            
        cons_tracerflux(i, j, k,:) = f1(:)
!!$
!!$        if ( ((flux_direction.eq.3).and.(i.eq.4).and.(j.eq.4)).or.&
!!$             ((flux_direction.eq.2).and.(i.eq.4).and.(k.eq.4)).or.&
!!$             ((flux_direction.eq.1).and.(j.eq.4).and.(k.eq.4))&
!!$             ) then
!!$          write(*,*) flux_direction, i, j, k, f1(1), consm_i(1), consp(1)
!!$        end if
        
      end do
    end do
  end do

      
end subroutine Whisky_HLLE_Tracer
