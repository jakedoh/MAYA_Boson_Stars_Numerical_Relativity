  /*@@
   @file      Whisky_CalcUpdate.F90
   @date      Thu Jan  11 11:03:32 2002
   @author    Ian Hawke
   @desc 
   Calculates the update terms given the fluxes. Moved to here so that 
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"
#include "SpaceMask.h"


 /*@@
   @routine    UpdateCalculation 
   @date       Wed Feb 13 11:03:32 2002
   @author     Ian Hawke
   @desc 
   Calculates the update terms from the fluxes.
   @enddesc 
   @calls     
   @calledby   
   @history 
   Moved out of the Riemann solver routines to make the FishEye /
   weighted flux calculation easier.
   @endhistory 

@@*/


subroutine UpdateCalculation(CCTK_ARGUMENTS)
    
  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT :: i,j,k,itracer
  CCTK_REAL :: idx, alp_l, alp_r

  CCTK_INT :: type_bits, atmosphere, not_atmosphere
  
  call SpaceMask_GetTypeBits(type_bits, "Hydro_Atmosphere")
  call SpaceMask_GetStateBits(atmosphere, "Hydro_Atmosphere",&
                              "in_atmosphere")
  call SpaceMask_GetStateBits(not_atmosphere, "Hydro_Atmosphere",&
                              "not_in_atmosphere")

  idx = 1.d0 / CCTK_DELTA_SPACE(flux_direction)

  if (CCTK_EQUALS(method_type, "RSA FV")) then

    if (use_weighted_fluxes == 0) then
  
      !$OMP PARALLEL DO PRIVATE (i,j,k,itracer, alp_l, alp_r)
      do k = whisky_stencil + 1, cctk_lsh(3) - whisky_stencil
        do j = whisky_stencil + 1, cctk_lsh(2) - whisky_stencil
          do i = whisky_stencil + 1, cctk_lsh(1) - whisky_stencil
            
            alp_l = 0.5d0 * (alp(i,j,k) + &
                 alp(i-xoffset,j-yoffset,k-zoffset))
            alp_r = 0.5d0 * (alp(i,j,k) + &
                 alp(i+xoffset,j+yoffset,k+zoffset))
            
            densrhs(i,j,k) = densrhs(i,j,k) + &
                 (alp_l * densflux(i-xoffset,j-yoffset,k-zoffset) - &
                 alp_r * densflux(i,j,k)) * idx
            sxrhs(i,j,k) = sxrhs(i,j,k) + &
                 (alp_l * sxflux(i-xoffset,j-yoffset,k-zoffset) - &
                 alp_r * sxflux(i,j,k)) * idx
            syrhs(i,j,k) = syrhs(i,j,k) + &
                 (alp_l * syflux(i-xoffset,j-yoffset,k-zoffset) - &
                 alp_r * syflux(i,j,k)) * idx
            szrhs(i,j,k) = szrhs(i,j,k) + &
                 (alp_l * szflux(i-xoffset,j-yoffset,k-zoffset) - &
                 alp_r * szflux(i,j,k)) * idx
            taurhs(i,j,k) = taurhs(i,j,k) + &
                 (alp_l * tauflux(i-xoffset,j-yoffset,k-zoffset) - &
                 alp_r * tauflux(i,j,k)) * idx

            if (whisky_mhd_handle.eq.2) then
               Bnxrhs(i,j,k) = Bnxrhs(i,j,k) + &
                    (alp_l * Bnxflux(i-xoffset,j-yoffset,k-zoffset) - &
                    alp_r * Bnxflux(i,j,k)) * idx
               Bnyrhs(i,j,k) = Bnyrhs(i,j,k) + &
                    (alp_l * Bnyflux(i-xoffset,j-yoffset,k-zoffset) - &
                    alp_r * Bnyflux(i,j,k)) * idx
               Bnzrhs(i,j,k) = Bnzrhs(i,j,k) + &
                    (alp_l * Bnzflux(i-xoffset,j-yoffset,k-zoffset) - &
                    alp_r * Bnzflux(i,j,k)) * idx
               if ( clean_divergence.ne.0 ) then
                  divclean_psirhs(i,j,k) = divclean_psirhs(i,j,k) + &
                     (alp_l * divclean_psiflux(i-xoffset,j-yoffset,k-zoffset) - &
                      alp_r * divclean_psiflux(i,j,k)) * idx * divclean_ch
               end if
            end if

            if (evolve_tracer .ne. 0) then
               do itracer=1,number_of_tracers
                  cons_tracerrhs(i,j,k,itracer) = cons_tracerrhs(i,j,k,itracer) + &
                       (alp_l * cons_tracerflux(i-xoffset,j-yoffset,k-zoffset,itracer) - &
                       alp_r * cons_tracerflux(i,j,k,itracer)) * idx
               enddo
!!$        if ( ((flux_direction.eq.3).and.(i.eq.4).and.(j.eq.4)).or.&
!!$             ((flux_direction.eq.2).and.(i.eq.4).and.(k.eq.4)).or.&
!!$             ((flux_direction.eq.1).and.(j.eq.4).and.(k.eq.4))&
!!$             ) then
!!$          write(*,*) flux_direction, i, j, k, cons_tracerrhs(i,j,k)
!!$        end if
            
            end if
            
            if (wk_atmosphere .eq. 1) then

              if ( (atmosphere_mask(i,j,k) .eq. 1) .or. &
                   (SpaceMask_CheckStateBitsF90(space_mask,i,j,k,type_bits,atmosphere)) ) then

!!$                We are in the atmosphere so the momentum flux must vanish

                sxrhs(i,j,k) = 0.d0
                syrhs(i,j,k) = 0.d0
                szrhs(i,j,k) = 0.d0

!!$                TODO: think about what to do with the B field in atmosphere

                if ( ( (atmosphere_mask(i-1,j  ,k  ) .eq. 1) .and. &
                       (atmosphere_mask(i+1,j  ,k  ) .eq. 1) .and. &
                       (atmosphere_mask(i  ,j-1,k  ) .eq. 1) .and. &
                       (atmosphere_mask(i  ,j+1,k  ) .eq. 1) .and. &
                       (atmosphere_mask(i  ,j  ,k-1) .eq. 1) .and. &
                       (atmosphere_mask(i  ,j  ,k+1) .eq. 1) &
                     ) .or. &
                     ( (SpaceMask_CheckStateBitsF90(space_mask,i-1,j  ,k  ,type_bits,atmosphere)) .and. &
                       (SpaceMask_CheckStateBitsF90(space_mask,i+1,j  ,k  ,type_bits,atmosphere)) .and. &
                       (SpaceMask_CheckStateBitsF90(space_mask,i  ,j-1,k  ,type_bits,atmosphere)) .and. &
                       (SpaceMask_CheckStateBitsF90(space_mask,i  ,j+1,k  ,type_bits,atmosphere)) .and. &
                       (SpaceMask_CheckStateBitsF90(space_mask,i  ,j  ,k-1,type_bits,atmosphere)) .and. &
                       (SpaceMask_CheckStateBitsF90(space_mask,i  ,j  ,k+1,type_bits,atmosphere)) & 
                     ) &
                    ) then

!!$                    All neighbours are also atmosphere so all rhs vanish

                    densrhs(i,j,k) = 0.d0
                    taurhs(i,j,k)  = 0.d0

                end if
              end if

            end if


          enddo
        enddo
      enddo
      !$OMP END PARALLEL DO
      
    else
      
      call CCTK_WARN(0, "Not supported")
      
!!$    do k = whisky_stencil + 1, cctk_lsh(3) - whisky_stencil
!!$      do j = whisky_stencil + 1, cctk_lsh(2) - whisky_stencil
!!$        do i = whisky_stencil + 1, cctk_lsh(1) - whisky_stencil
!!$          
!!$          alp_l = 0.5d0 * (alp(i,j,k) + &
!!$               alp(i-xoffset,j-yoffset,k-zoffset))
!!$          alp_r = 0.5d0 * (alp(i,j,k) + &
!!$               alp(i+xoffset,j+yoffset,k+zoffset))
!!$          
!!$          densrhs(i,j,k) = densrhs(i,j,k) + &
!!$               (alp_l * &
!!$               &cell_surface(i-xoffset,j-yoffset,k-zoffset,flux_direction) * &
!!$               &densflux(i-xoffset,j-yoffset,k-zoffset) - &
!!$               alp_r * &
!!$               &cell_surface(i,j,k,flux_direction) * &
!!$               &densflux(i,j,k)) * idx / cell_volume(i,j,k)
!!$          sxrhs(i,j,k) = sxrhs(i,j,k) + &
!!$               (alp_l * &
!!$               &cell_surface(i-xoffset,j-yoffset,k-zoffset,flux_direction) * &
!!$               &sxflux(i-xoffset,j-yoffset,k-zoffset) - &
!!$               alp_r * &
!!$               &cell_surface(i,j,k,flux_direction) * &
!!$               &sxflux(i,j,k)) * idx / cell_volume(i,j,k)
!!$          syrhs(i,j,k) = syrhs(i,j,k) + &
!!$               (alp_l * &
!!$               &cell_surface(i-xoffset,j-yoffset,k-zoffset,flux_direction) * &
!!$               &syflux(i-xoffset,j-yoffset,k-zoffset) - &
!!$               alp_r * &
!!$               &cell_surface(i,j,k,flux_direction) * &
!!$               &syflux(i,j,k)) * idx / cell_volume(i,j,k)
!!$          szrhs(i,j,k) = szrhs(i,j,k) + &
!!$               (alp_l * &
!!$               &cell_surface(i-xoffset,j-yoffset,k-zoffset,flux_direction) * &
!!$               &szflux(i-xoffset,j-yoffset,k-zoffset) - &
!!$               alp_r * &
!!$               &cell_surface(i,j,k,flux_direction) * &
!!$               &szflux(i,j,k)) * idx / cell_volume(i,j,k)
!!$          taurhs(i,j,k) = taurhs(i,j,k) + &
!!$               (alp_l * &
!!$               &cell_surface(i-xoffset,j-yoffset,k-zoffset,flux_direction) * &
!!$               &tauflux(i-xoffset,j-yoffset,k-zoffset) - &
!!$               alp_r * &
!!$               &cell_surface(i,j,k,flux_direction) * &
!!$               &tauflux(i,j,k)) * idx / cell_volume(i,j,k)
!!$          
!!$        enddo
!!$      enddo
!!$    enddo

    end if
  
  else if (CCTK_EQUALS(method_type, "Flux split FD")) then

    !$OMP PARALLEL DO PRIVATE (i,j,k)
    do k = whisky_stencil + 1, cctk_lsh(3) - whisky_stencil
      do j = whisky_stencil + 1, cctk_lsh(2) - whisky_stencil
        do i = whisky_stencil + 1, cctk_lsh(1) - whisky_stencil
          
          densrhs(i,j,k) = densrhs(i,j,k) + &
               (densflux(i-xoffset,j-yoffset,k-zoffset) - &
                densflux(i,j,k)) * idx
          sxrhs(i,j,k) = sxrhs(i,j,k) + &
               (sxflux(i-xoffset,j-yoffset,k-zoffset) - &
                sxflux(i,j,k)) * idx
          syrhs(i,j,k) = syrhs(i,j,k) + &
               (syflux(i-xoffset,j-yoffset,k-zoffset) - &
                syflux(i,j,k)) * idx
          szrhs(i,j,k) = szrhs(i,j,k) + &
               (szflux(i-xoffset,j-yoffset,k-zoffset) - &
                szflux(i,j,k)) * idx
          taurhs(i,j,k) = taurhs(i,j,k) + &
               (tauflux(i-xoffset,j-yoffset,k-zoffset) - &
                tauflux(i,j,k)) * idx

          if (whisky_mhd_handle.eq.2) then
             Bnxrhs(i,j,k) = Bnxrhs(i,j,k) + &
                  (Bnxflux(i-xoffset,j-yoffset,k-zoffset) - &
                   Bnxflux(i,j,k)) * idx
             Bnyrhs(i,j,k) = Bnyrhs(i,j,k) + &
                  (Bnyflux(i-xoffset,j-yoffset,k-zoffset) - &
                   Bnyflux(i,j,k)) * idx
             Bnzrhs(i,j,k) = Bnzrhs(i,j,k) + &
                  (Bnzflux(i-xoffset,j-yoffset,k-zoffset) - &
                   Bnzflux(i,j,k)) * idx
             if ( clean_divergence.ne.0 ) then
                divclean_psirhs(i,j,k) = divclean_psirhs(i,j,k) + &
                   (divclean_psiflux(i-xoffset,j-yoffset,k-zoffset) - &
                    divclean_psiflux(i,j,k)) * idx * divclean_ch
             end if
          end if
          
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO
      
  end if
  
  return

end subroutine UpdateCalculation
  
