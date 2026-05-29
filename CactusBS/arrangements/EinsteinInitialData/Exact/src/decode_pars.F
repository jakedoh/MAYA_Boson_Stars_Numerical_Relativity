c/*@@
c @file		decode_pars.F
c @date		Fri Jun  7 19:47:46 CEST 2002
c @author	Jonathan Thornburg <jthorn@aei.mpg.de>
c @desc
c	Decode/copy parameters for this thorn into grid scalars
c	so we can share this with friends
c	so the Calc_Tmunu code in ../include/Scalar_CalcTmunu.inc
c	can use them in computing the stress-energy tensor
c
c	Actually we only have to copy those parameters which are
c	used by the Calc_Tmunu code.  For simplicity we decode
c	exact_model into the integer decoded_exact_model, then
c	copy only the per-model parameters for those models which
c	have stress-energy tensor code in ../include/Scalar_CalcTmunu.inc .
c
c @enddesc
c @version	$Header$
c@@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

#include "param_defs.inc"

c/*@@
c @routine	Exact__decode_pars
c @date		Fri Jun  7 19:47:46 CEST 2002
c @author	Jonathan Thornburg <jthorn@aei.mpg.de>
c @desc
c	Decode/copy parameters for this thorn into grid scalars, so
c	we can share this with friends for the use of the Calc_Tmunu
c	code in ../include/Scalar_CalcTmunu.inc in computing the
c	stress-energy tensor.
c @enddesc
c @version	$Header$
c@@*/
	subroutine Exact__decode_pars(CCTK_ARGUMENTS)
	implicit none
 	DECLARE_CCTK_ARGUMENTS_EXACT__DECODE_PARS
	DECLARE_CCTK_PARAMETERS
	DECLARE_CCTK_FUNCTIONS

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c
c decode  exact_model  into the integer  decoded_exact_model
c

c Minkowski spacetime
	if      (CCTK_Equals(exact_model, "Minkowski") .ne. 0) then
	     decoded_exact_model = EXACT__Minkowski
	elseif (CCTK_Equals(exact_model, "Minkowski/shift") .ne. 0) then
	     decoded_exact_model = EXACT__Minkowski_shift
	elseif (CCTK_Equals(exact_model, "Minkowski/funny") .ne. 0) then
	     decoded_exact_model = EXACT__Minkowski_funny
	elseif (CCTK_Equals(exact_model, "Minkowski/gauge wave") .ne. 0) then
	     decoded_exact_model = EXACT__Minkowski_gauge_wave
	elseif (CCTK_Equals(exact_model, "Minkowski/shifted gauge wave") .ne. 0) then
	     decoded_exact_model = EXACT__Minkowski_shifted_gauge_wave
	elseif (CCTK_Equals(exact_model, "Minkowski/conf wave") .ne. 0) then
	     decoded_exact_model = EXACT__Minkowski_conf_wave

c black hole spacetimes
	elseif (CCTK_Equals(exact_model, "Schwarzschild/EF") .ne. 0) then
	     decoded_exact_model = EXACT__Schwarzschild_EF
	elseif (CCTK_Equals(exact_model, "Schwarzschild/PG") .ne. 0) then
	     decoded_exact_model = EXACT__Schwarzschild_PG
	elseif (CCTK_Equals(exact_model, "Schwarzschild/BL") .ne. 0) then
	     decoded_exact_model = EXACT__Schwarzschild_BL
	elseif (CCTK_Equals(exact_model, "Schwarzschild/Novikov") .ne. 0) then
	     decoded_exact_model = EXACT__Schwarzschild_Novikov
	elseif (CCTK_Equals(exact_model, "Kerr/Boyer-Lindquist") .ne. 0) then
	     decoded_exact_model = EXACT__Kerr_BoyerLindquist
	elseif (CCTK_Equals(exact_model, "Kerr/Kerr-Schild") .ne. 0) then
	     decoded_exact_model = EXACT__Kerr_KerrSchild
	elseif (CCTK_Equals(exact_model, "Kerr/Kerr-Schild/spherical") .ne. 0) then
	     decoded_exact_model = EXACT__Kerr_KerrSchild_spherical
	elseif (CCTK_Equals(exact_model, "Schwarzschild-Lemaitre") .ne. 0) then
	     decoded_exact_model = EXACT__Schwarzschild_Lemaitre
	elseif (CCTK_Equals(exact_model, "multi-BH") .ne. 0) then
	     decoded_exact_model = EXACT__multi_BH
	elseif (CCTK_Equals(exact_model, "Alvi") .ne. 0) then
	     decoded_exact_model = EXACT__Alvi
	elseif (CCTK_Equals(exact_model, "Thorne-fakebinary") .ne. 0) then
	     decoded_exact_model = EXACT__Thorne_fakebinary

c cosmological spacetimes
	elseif (CCTK_Equals(exact_model, "Lemaitre") .ne. 0) then
	     decoded_exact_model = EXACT__Lemaitre
C this metric doesnt work and has been moved to ../archive/
CC	elseif (CCTK_Equals(exact_model, "Robertson-Walker") .ne. 0) then
CC	     decoded_exact_model = EXACT__Robertson_Walker
	elseif (CCTK_Equals(exact_model, "de Sitter") .ne. 0) then
	     decoded_exact_model = EXACT__de_Sitter
	elseif (CCTK_Equals(exact_model, "de Sitter+Lambda") .ne. 0) then
	     decoded_exact_model = EXACT__de_Sitter_Lambda
	elseif (CCTK_Equals(exact_model, "anti-de Sitter+Lambda") .ne. 0) then
	     decoded_exact_model = EXACT__anti_de_Sitter_Lambda
	elseif (CCTK_Equals(exact_model, "Bianchi I") .ne. 0) then
	     decoded_exact_model = EXACT__Bianchi_I
	elseif (CCTK_Equals(exact_model, "Goedel") .ne. 0) then
	     decoded_exact_model = EXACT__Goedel
	elseif (CCTK_Equals(exact_model, "Bertotti") .ne. 0) then
	     decoded_exact_model = EXACT__Bertotti
	elseif (CCTK_Equals(exact_model, "Kasner-like") .ne. 0) then
	     decoded_exact_model = EXACT__Kasner_like
	elseif (CCTK_Equals(exact_model, "Kasner-axisymmetric") .ne. 0) then
	     decoded_exact_model = EXACT__Kasner_axisymmetric
	elseif (CCTK_Equals(exact_model, "Kasner-generalized") .ne. 0) then
	     decoded_exact_model = EXACT__Kasner_generalized
	elseif (CCTK_Equals(exact_model, "Gowdy-wave") .ne. 0) then
	     decoded_exact_model = EXACT__Gowdy_wave
	elseif (CCTK_Equals(exact_model, "Milne") .ne. 0) then
	     decoded_exact_model = EXACT__Milne

c miscellaneous spacetimes
	elseif (CCTK_Equals(exact_model, "boost-rotation symmetric") .ne. 0) then
	     decoded_exact_model = EXACT__boost_rotation_symmetric
	elseif (CCTK_Equals(exact_model, "bowl") .ne. 0) then
	     decoded_exact_model = EXACT__bowl
	elseif (CCTK_Equals(exact_model, "constant density star") .ne. 0) then
	     decoded_exact_model = EXACT__constant_density_star
	else
	     call CCTK_WARN(0, "Unknown exact_model")
	endif

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c
c parameters for Schwarzschild-Lemaitre metric
c (Schwarzschiled black hole with cosmological constant)
c
	Schwarzschild_Lemaitre___Lambda = Schwarzschild_Lemaitre__Lambda
	Schwarzschild_Lemaitre___mass   = Schwarzschild_Lemaitre__mass

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c
c parameters for Lemaitre-type spacetime
c
	Lemaitre___kappa    = Lemaitre__kappa
	Lemaitre___Lambda   = Lemaitre__Lambda
	Lemaitre___epsilon0 = Lemaitre__epsilon0
	Lemaitre___R0       = Lemaitre__R0

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
CCC
CCC this metric doesnt work and has been moved to ../archive/
CCc
CCc parameters for Robertson-Walker spacetime
CCc
CC	Robertson_Walker___R0       = Robertson_Walker__R0
CC	Robertson_Walker___rho      = Robertson_Walker__rho
CC	Robertson_Walker___k        = Robertson_Walker__k
CC	Robertson_Walker___pressure = Robertson_Walker__pressure
CC
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c
c parameters for de Sitter spacetime
c
	de_Sitter___scale = de_Sitter__scale

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c
c parameters for de Sitter spacetime with cosmological constant
c
	de_Sitter_Lambda___scale = de_Sitter_Lambda__scale

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c
c parameters for anti-de Sitter spacetime with cosmological constant
c
	anti_de_Sitter_Lambda___scale = anti_de_Sitter_Lambda__scale

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c
c parameters for Bertotti spacetime
c
	Bertotti___Lambda = Bertotti__Lambda

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c
c parameters for Kasner-like spacetime
c
	Kasner_like___q = Kasner_like__q

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c
c parameters for generalized Kasner spacetime
c
	Kasner_generalized___p1 = Kasner_generalized__p1
	Kasner_generalized___p2 = Kasner_generalized__p2

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c
c parameters for constant density (Schwarzschild) star
c
	constant_density_star___mass   = constant_density_star__mass
	constant_density_star___radius = constant_density_star__radius

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	return
	end
