#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine Exact_AddToTmunu(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_LOOP3_ALL_DECLARE(Exact_AddToTmunu)

  integer :: i,j,k

!C temporary variables
  CCTK_REAL  unu,doi,trei,rr,kkkk
  CCTK_REAL  aaaa,bbbb,aaaa1,bbbb1,kkkk1,r3,bass,term1,term2
  CCTK_REAL  unu1, raz, raz2, razsch2, coefsch, pppsch, unusch
  CCTK_REAL  treiori
  CCTK_REAL  star_m, star_r
  CCTK_REAL, dimension(cctk_ash(1),cctk_ash(2),cctk_ash(3)) :: eTxx, eTxy, &
             eTxz, eTyy, eTyz, eTzz, eTtx, eTty, eTtz, eTtt
  CCTK_POINTER pTxx, pTxy, pTxz, pTyy, pTyz, pTzz, pTtx, pTty, pTtz, pTtt
  integer :: vi_eTxx, vi_eTxy, vi_eTxz, vi_eTyy, vi_eTyz, vi_eTzz, vi_eTtx, &
             vi_eTty, vi_eTtz, vi_eTtt

  CCTK_INT, dimension(10) :: variable_list, timelevel_list, where_list
  CCTK_INT :: ierr
  integer :: have_presync

  pointer (pTxx, eTxx)
  pointer (pTxy, eTxy)
  pointer (pTxz, eTxz)
  pointer (pTyy, eTyy)
  pointer (pTyz, eTyz)
  pointer (pTzz, eTzz)
  pointer (pTtx, eTtx)
  pointer (pTty, eTty)
  pointer (pTtz, eTtz)
  pointer (pTtt, eTtt)

  call CCTK_VarIndex(vi_eTxx, "TmunuBase::eTxx")
  call CCTK_VarIndex(vi_eTxy, "TmunuBase::eTxy")
  call CCTK_VarIndex(vi_eTxz, "TmunuBase::eTxz")
  call CCTK_VarIndex(vi_eTyy, "TmunuBase::eTyy")
  call CCTK_VarIndex(vi_eTyz, "TmunuBase::eTyz")
  call CCTK_VarIndex(vi_eTzz, "TmunuBase::eTzz")
  call CCTK_VarIndex(vi_eTtx, "TmunuBase::eTtx")
  call CCTK_VarIndex(vi_eTty, "TmunuBase::eTty")
  call CCTK_VarIndex(vi_eTtz, "TmunuBase::eTtz")
  call CCTK_VarIndex(vi_eTtt, "TmunuBase::eTtt")

  call CCTK_VarDataPtrI(pTxx, cctkGH, 0, vi_eTxx)
  call CCTK_VarDataPtrI(pTxy, cctkGH, 0, vi_eTxy)
  call CCTK_VarDataPtrI(pTxz, cctkGH, 0, vi_eTxz)
  call CCTK_VarDataPtrI(pTyy, cctkGH, 0, vi_eTyy)
  call CCTK_VarDataPtrI(pTyz, cctkGH, 0, vi_eTyz)
  call CCTK_VarDataPtrI(pTzz, cctkGH, 0, vi_eTzz)
  call CCTK_VarDataPtrI(pTtx, cctkGH, 0, vi_eTtx)
  call CCTK_VarDataPtrI(pTty, cctkGH, 0, vi_eTty)
  call CCTK_VarDataPtrI(pTtz, cctkGH, 0, vi_eTtz)
  call CCTK_VarDataPtrI(pTtt, cctkGH, 0, vi_eTtt)

  call CCTK_IsFunctionAliased(have_presync, "Driver_RequireValidData")
  ! tell Cactus that we are modifying Tmunu
  if (have_presync .ne. 0) then
    variable_list = (/ vi_eTxx, vi_eTxy, vi_eTxz, vi_eTyy, vi_eTyz, vi_eTzz, &
                     vi_eTtx, vi_eTty, vi_eTtz, vi_eTtt/)
    timelevel_list = 0
    where_list = CCTK_VALID_EVERYWHERE
    ierr = Driver_RequireValidData(cctkGH, variable_list, timelevel_list, 10, &
                                   where_list)
  end if

  CCTK_LOOP3_ALL(Exact_AddToTmunu, i,j,k)

!C Here we added the matter variables for several of the metrics
!C you can find in "src" directory through the components of the
!C stress-eergy tensor. Being different in specific cases it is
!C necessary to check first which metric is running.  
!C Author : Dumitru Vulcanov (Timisoara, Romania)
!C $Header$

!C Varianta cu un singur param. (raza initiala) pt Rob-Walker

#include "param_defs.inc"

!c
!c FIXME:
!c If we could be certain that this code were always compiled as Fortran 90,
!c the decode here could be done with a case statement, which ought to give
!c better performance than the if-else chain we use now.  But in practice
!c Cactus has enough other overheads that this is not bothering with...
!c

!c
!c ***** KLUDGE *****
!c
!c This code is #include-d into various evolution thorns, and alas does not
!c have direct access to thorn Exact parameters.  Instead, this code must
!c use the restricted-grid-scalar copies of the parameters.  In practice,
!c this means changing "__" to "___" in all parameter names (but not in the
!c #define constants in "param_defs.inc").  See the comments in param.ccl
!c for further information on this.
!c

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccc Minkowski spacetime cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 if     (decoded_exact_model .eq. EXACT__Minkowski) then
!c        no stress-energy tensor in this model
 elseif (decoded_exact_model .eq. EXACT__Minkowski_shift) then
!c        no stress-energy tensor in this model
 elseif (decoded_exact_model .eq. EXACT__Minkowski_funny) then
!c        no stress-energy tensor in this model
 elseif (decoded_exact_model .eq. EXACT__Minkowski_gauge_wave) then
!c        no stress-energy tensor in this model
 elseif (decoded_exact_model .eq. EXACT__Minkowski_shifted_gauge_wave) then
!c        no stress-energy tensor in this model
 elseif (decoded_exact_model .eq. EXACT__Minkowski_conf_wave) then
!c        no stress-energy tensor in this model

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccc black hole spacetimes cccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 elseif (decoded_exact_model .eq. EXACT__Schwarzschild_EF) then
!c        no stress-energy tensor in this model
 elseif (decoded_exact_model .eq. EXACT__Schwarzschild_PG) then
!c        no stress-energy tensor in this model
 elseif (decoded_exact_model .eq. EXACT__Schwarzschild_BL) then
!c        no stress-energy tensor in this model
 elseif (decoded_exact_model .eq. EXACT__Schwarzschild_Novikov) then
!c        no stress-energy tensor in this model
 elseif (decoded_exact_model .eq. EXACT__Kerr_BoyerLindquist) then
!c        no stress-energy tensor in this model
 elseif (decoded_exact_model .eq. EXACT__Kerr_KerrSchild) then
!c        no stress-energy tensor in this model
 elseif (decoded_exact_model .eq. EXACT__Kerr_KerrSchild_spherical) then
!c        no stress-energy tensor in this model

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!c
!c Schwarzschild-Lemaitre spacetime
!c (Schwarzschild black hole with cosmological constant)
!c
 elseif (decoded_exact_model .eq. EXACT__Schwarzschild_Lemaitre) then
   razsch2=x(i,j,k)*x(i,j,k)+y(i,j,k)*y(i,j,k)+z(i,j,k)*z(i,j,k)
   coefsch=-Schwarzschild_Lemaitre___Lambda/(8.0D0*EXACT__pi)
   pppsch=1.0D0-2.0D0*Schwarzschild_Lemaitre___mass/sqrt(razsch2) &
&               -Schwarzschild_Lemaitre___Lambda*razsch2/3.0D0
   unusch=(1.0D0-pppsch)/pppsch/razsch2

   eTtt(i,j,k) =  eTtt(i,j,k)-coefsch*pppsch
   eTtx(i,j,k) =  eTtx(i,j,k)
   eTty(i,j,k) =  eTty(i,j,k)
   eTtz(i,j,k) =  eTtz(i,j,k)
   eTxx(i,j,k)  = eTxx(i,j,k)+coefsch*(1.0D0+x(i,j,k)*x(i,j,k)*unusch)
   eTyy(i,j,k)  = eTyy(i,j,k)+coefsch*(1.0D0+y(i,j,k)*y(i,j,k)*unusch)
   eTzz(i,j,k)  = eTzz(i,j,k)+coefsch*(1.0D0+z(i,j,k)*z(i,j,k)*unusch)
   eTxy(i,j,k)  = eTxy(i,j,k)+coefsch*x(i,j,k)*y(i,j,k)*unusch
   eTxz(i,j,k)  = eTxz(i,j,k)+coefsch*x(i,j,k)*z(i,j,k)*unusch
   eTyz(i,j,k)  = eTyz(i,j,k)+coefsch*y(i,j,k)*z(i,j,k)*unusch

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 elseif (decoded_exact_model .eq. EXACT__multi_BH) then
!c        no stress-energy tensor in this model
 elseif (decoded_exact_model .eq. EXACT__Alvi) then
!c        no stress-energy tensor in this model
 elseif (decoded_exact_model .eq. EXACT__Thorne_fakebinary) then
!c        no stress-energy tensor in this model

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccc cosmological spacetimes cccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!c
!c Lemaitre spacetime
!c
  elseif (decoded_exact_model .eq. EXACT__Lemaitre) then
    unu1 = sqrt(3.0D0*Lemaitre___Lambda) &
&           * CCTK_time * (Lemaitre___kappa+1.0D0) / (2.0D0)
    raz  = Lemaitre___R0*(cosh(unu1) &
&                         + sqrt(1.0D0+8.0D0*EXACT__pi*Lemaitre___epsilon0 &
&                                               /Lemaitre___Lambda) &
&                           *sinh(unu1)) &
&                        **(2.0D0/(3.0D0*Lemaitre___kappa+3.0D0)) 
    raz2 = raz*raz
    treiori = -Lemaitre___Lambda*raz2/8.0D0/EXACT__pi &
&           +Lemaitre___epsilon0*Lemaitre___kappa &
&                              *raz**(-3.0D0*Lemaitre___kappa-1.0D0)

    eTtt(i,j,k) = eTtt(i,j,k) + Lemaitre___Lambda/8.0D0/EXACT__pi &
&              + Lemaitre___epsilon0*raz**(-3.0D0*(Lemaitre___kappa+1.0D0))
    eTxx(i,j,k) = eTxx(i,j,k) + treiori
    eTyy(i,j,k) = eTyy(i,j,k) + treiori
    eTzz(i,j,k) = eTzz(i,j,k) + treiori
    eTxy(i,j,k) = eTxy(i,j,k)
    eTyz(i,j,k) = eTyz(i,j,k)
    eTxz(i,j,k) = eTxz(i,j,k)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccc
!ccc this metric doesnt work and has been moved to ../../archive/
!ccc
!ccc Robertson-Walker spacetime
!ccc
!cc        elseif (decoded_exact_model .eq. EXACT__Robertson_Walker) then
!cc          rr2  = x(i,j,k)*x(i,j,k)+y(i,j,k)*y(i,j,k)+z(i,j,k)*z(i,j,k)
!cc
!cc          if (Robertson_Walker___pressure .gt. 0) then
!cc            aha1  = Robertson_Walker___k * (Robertson_Walker___R0**2)
!cc     &                                  / (8.0D0*EXACT__pi*(raza(i,j,k)**2))
!cc            aha2 = Robertson_Walker___k/(1.0D0 - Robertson_Walker___k*rr2) 
!cc
!cc            eTtt(i,j,k) =  eTtt(i,j,k) + 3.0D0*aha1/(raza(i,j,k)*raza(i,j,k))
!cc            eTxx(i,j,k) =  eTxx(i,j,k) + aha1*(1.0D0 + aha2*x(i,j,k)*x(i,j,k))
!cc            eTyy(i,j,k) =  eTyy(i,j,k) + aha1*(1.0D0 + aha2*y(i,j,k)*y(i,j,k))
!cc            eTzz(i,j,k) =  eTzz(i,j,k) + aha1*(1.0D0 + aha2*z(i,j,k)*z(i,j,k))
!cc            eTxy(i,j,k) =  eTxy(i,j,k) + aha1*aha2*x(i,j,k)*y(i,j,k)
!cc            eTxz(i,j,k) =  eTxz(i,j,k) + aha1*aha2*x(i,j,k)*z(i,j,k)
!cc            eTyz(i,j,k) =  eTyz(i,j,k) + aha1*aha2*y(i,j,k)*y(i,j,k)
!cc          else
!cc            eTtt(i,j,k) = eTtt(i,j,k)+Robertson_Walker___rho * (Robertson_Walker___R0**3)
!cc     &                                      / (raza(i,j,k)**3)
!cc            eTxx(i,j,k) = eTxx(i,j,k)
!cc            eTyy(i,j,k) = eTyy(i,j,k)
!cc            eTzz(i,j,k) = eTzz(i,j,k)
!cc            eTxy(i,j,k) = eTxy(i,j,k)
!cc            eTxz(i,j,k) = eTxz(i,j,k)
!cc            eTyz(i,j,k) = eTyz(i,j,k)
!cc          endif
!cc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!c
!c de Sitter spacetime
!c
  elseif (decoded_exact_model .eq. EXACT__de_Sitter) then
    eTtt(i,j,k) =  eTtt(i,j,k) + 1.0D0/6.0D0/EXACT__pi/(CCTK_time**2)
    eTtx(i,j,k)  = eTtx(i,j,k)
    eTty(i,j,k)  = eTty(i,j,k)
    eTtz(i,j,k)  = eTtz(i,j,k)
    eTxx(i,j,k)  = eTxx(i,j,k)
    eTyy(i,j,k)  = eTyy(i,j,k)
    eTzz(i,j,k)  = eTzz(i,j,k)
    eTxy(i,j,k)  = eTxy(i,j,k)
    eTxz(i,j,k)  = eTxz(i,j,k)
    eTyz(i,j,k)  = eTyz(i,j,k)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!c
!c de Sitter spacetime with cosmological constant
!c
  elseif (decoded_exact_model .eq. EXACT__de_Sitter_Lambda) then
    aaaa = de_Sitter_Lambda___scale/(8.0D0*EXACT__pi)
    bbbb = aaaa*exp(2.0D0*sqrt(de_Sitter_Lambda___scale/3.0D0)*CCTK_time)

    eTtt(i,j,k)  = eTtt(i,j,k) + aaaa
    eTtx(i,j,k)  = eTtx(i,j,k)
    eTty(i,j,k)  = eTty(i,j,k)
    eTtz(i,j,k)  = eTtz(i,j,k)
    eTxx(i,j,k)  = eTxx(i,j,k) - bbbb
    eTyy(i,j,k)  = eTyy(i,j,k) - bbbb
    eTzz(i,j,k)  = eTzz(i,j,k) - bbbb
    eTxy(i,j,k)  = eTxy(i,j,k)
    eTxz(i,j,k)  = eTxz(i,j,k)
    eTyz(i,j,k)  = eTyz(i,j,k)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!c
!c anti-de Sitter spacetime with cosmological constant
!c

  elseif (decoded_exact_model .eq. EXACT__anti_de_Sitter_Lambda) then
    aaaa1 = anti_de_Sitter_Lambda___scale/(8.0D0*EXACT__pi)
    bbbb1 = aaaa1*exp(2.0D0*sqrt(-anti_de_Sitter_Lambda___scale/3.0D0) &
&                           *x(i,j,k))
  
    eTtt(i,j,k) =  eTtt(i,j,k) + bbbb1
    eTtx(i,j,k)  = eTtx(i,j,k)
    eTty(i,j,k)  = eTty(i,j,k)
    eTtz(i,j,k)  = eTtz(i,j,k)
    eTxx(i,j,k)  = eTxx(i,j,k) - aaaa1
    eTyy(i,j,k)  = eTyy(i,j,k) - bbbb1
    eTzz(i,j,k)  = eTzz(i,j,k) - bbbb1
    eTxy(i,j,k)  = eTxy(i,j,k)
    eTxz(i,j,k)  = eTxz(i,j,k)
    eTyz(i,j,k)  = eTyz(i,j,k)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 elseif (decoded_exact_model .eq. EXACT__Bianchi_I) then
!c        no stress-energy tensor in this model
 elseif (decoded_exact_model .eq. EXACT__Goedel) then
!c        no stress-energy tensor in this model

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!c
!c Bertotti spacetime
!c
  elseif (decoded_exact_model .eq. EXACT__Bertotti) then
    bass  = Bertotti___Lambda/(8.0D0*EXACT__pi)
    term1 = bass*exp(2.0D0*sqrt(-Bertotti___Lambda)*x(i,j,k))
    term2 = bass*exp(2.0D0*sqrt(-Bertotti___Lambda)*z(i,j,k))

    eTtt(i,j,k) =  eTtt(i,j,k) + term1
    eTtx(i,j,k)  = eTtx(i,j,k)
    eTty(i,j,k)  = eTty(i,j,k)
    eTtz(i,j,k)  = eTtz(i,j,k)
    eTxx(i,j,k)  = eTxx(i,j,k) - bass
    eTyy(i,j,k)  = eTyy(i,j,k) - term2
    eTzz(i,j,k)  = eTzz(i,j,k) - bass
    eTxy(i,j,k)  = eTxy(i,j,k)
    eTxz(i,j,k)  = eTxz(i,j,k)
    eTyz(i,j,k)  = eTyz(i,j,k)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!c
!c Kasner-like spacetime
!c
  elseif (decoded_exact_model .eq. EXACT__Kasner_like) then

    kkkk=Kasner_like___q*(2.0D0-3.0D0*Kasner_like___q) &
&                        /(8.0D0*EXACT__pi*(CCTK_time**2))

    eTtt(i,j,k)  = eTtt(i,j,k) + kkkk
    eTtx(i,j,k)  = eTtx(i,j,k)
    eTty(i,j,k)  = eTty(i,j,k)
    eTtz(i,j,k)  = eTtz(i,j,k)
    eTxx(i,j,k)  = eTxx(i,j,k) + kkkk*CCTK_time**(2.0D0*Kasner_like___q)
    eTyy(i,j,k)  = eTyy(i,j,k) + kkkk*CCTK_time**(2.0D0*Kasner_like___q)
    eTzz(i,j,k)  = eTzz(i,j,k) + kkkk*CCTK_time**(2.0D0-4.0D0*Kasner_like___q)
    eTxy(i,j,k)  = eTxy(i,j,k)
    eTxz(i,j,k)  = eTxz(i,j,k)
    eTyz(i,j,k)  = eTyz(i,j,k)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 elseif (decoded_exact_model .eq. EXACT__Kasner_axisymmetric) then
!c        no stress-energy tensor in this model

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!c
!c generalized Kasner spacetime
!c
  elseif (decoded_exact_model .eq. EXACT__Kasner_generalized) then

  kkkk1 = (  Kasner_generalized___p1 - Kasner_generalized___p1**2 &
&           + Kasner_generalized___p2 - Kasner_generalized___p2**2 &
&           - Kasner_generalized___p1*Kasner_generalized___p2 ) &
&          / (8.0D0*EXACT__pi*(CCTK_time**2))

  eTtt(i,j,k)  = eTtt(i,j,k) + kkkk1
  eTtx(i,j,k)  = eTtx(i,j,k)
  eTty(i,j,k)  = eTty(i,j,k)
  eTtz(i,j,k)  = eTtz(i,j,k)
  eTxx(i,j,k)  = eTxx(i,j,k)+kkkk1*CCTK_time**(2.0D0*Kasner_generalized___p1)
  eTyy(i,j,k)  = eTyy(i,j,k)+kkkk1*CCTK_time**(2.0D0*Kasner_generalized___p2)
  eTzz(i,j,k)  = eTzz(i,j,k)+kkkk1*CCTK_time**(2.0D0-2.0D0*Kasner_generalized___p1 &
&                                    -2.0D0*Kasner_generalized___p2)
  eTxy(i,j,k)  = eTxy(i,j,k)
  eTxz(i,j,k)  = eTxz(i,j,k)
  eTyz(i,j,k)  = eTyz(i,j,k)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 elseif (decoded_exact_model .eq. EXACT__Gowdy_wave) then
!c        no stress-energy tensor in this model


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 elseif (decoded_exact_model .eq. EXACT__Milne) then
!c        no stress-energy tensor in this model

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccc miscellaneous spacetimes ccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 elseif (decoded_exact_model .eq. EXACT__boost_rotation_symmetric) then
!c        no stress-energy tensor in this model
 elseif (decoded_exact_model .eq. EXACT__bowl) then
!c        no stress-energy tensor in this model

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!c
!c constant density star
!c
  elseif (decoded_exact_model .eq. EXACT__constant_density_star) then
    rr  = sqrt(x(i,j,k)*x(i,j,k)+y(i,j,k)*y(i,j,k)+ &
&                z(i,j,k)*z(i,j,k))
    star_m = constant_density_star___mass
    star_r = constant_density_star___radius

    r3 = star_r**3
    if (rr.le.star_r) then
      unu = 3.0D0*sqrt(1.0D0-2.0D0*star_m/star_r)
      doi = sqrt(1.0D0-2.0D0*star_m*rr*rr/r3)
      trei= star_m*(unu-3.0D0*doi)/(2*EXACT__pi*(unu-doi)*r3)
      eTtt(i,j,k) = eTtt(i,j,k) + 3.0D0*star_m* &
&        (5.0D0-9.0D0*star_m/star_r - unu*doi &
&        -star_m*rr*rr/r3)/(8.0D0*EXACT__pi*r3)
      eTxx(i,j,k) = eTxx(i,j,k) -trei*(1.0D0+2.0D0*star_m*x(i,j,k)*x(i,j,k)/ &
&      (doi*doi*r3))/2.0D0
      eTyy(i,j,k) = eTyy(i,j,k) -trei*(1.0D0+2.0D0*star_m*y(i,j,k)*y(i,j,k)/ &
&      (doi*doi*r3))/2.0D0
      eTzz(i,j,k) = eTzz(i,j,k) -trei*(1.0D0+2.0D0*star_m*z(i,j,k)*z(i,j,k)/ &
&      (doi*doi*r3))/2.0D0
      eTxy(i,j,k) = eTxy(i,j,k) -trei*star_m*x(i,j,k)*y(i,j,k)/(doi*doi*r3)
      eTyz(i,j,k) = eTyz(i,j,k) -trei*star_m*y(i,j,k)*z(i,j,k)/(doi*doi*r3)
      eTxz(i,j,k) = eTxz(i,j,k) -trei*star_m*x(i,j,k)*z(i,j,k)/(doi*doi*r3)
    endif

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  else
    call CCTK_WARN(0,"Unknown value of Exact::decoded_exact_model")
  endif 

  CCTK_ENDLOOP3_ALL(Exact_AddToTmunu)

  ! tell Cactus that we are modifying Tmunu
  if (have_presync .ne. 0) then
    ierr = Driver_NotifyDataModified(cctkGH, variable_list, timelevel_list, &
                                     10, where_list)
  end if

end subroutine
