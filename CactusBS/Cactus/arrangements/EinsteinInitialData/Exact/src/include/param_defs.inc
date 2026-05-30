! param_defs.inc -- integer constants for decoded_exact_model
! $Header$

!
! For reasons explained in our param.ccl file, we decode the  exact_model
! parameter, and all our other string-valued parameters used in computing
! the stress-energy tensor, into integers.  This file contains #define
! definitions for those integers.
!

!
! For each parameter, the value 0 is deliberately *not* a legal value
! for the decoded integer, to help catch bugs where the decoded integer
! is not initialized properly.
!

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!
! ***** definitions for decoded_exact_model *****
!
! These *MUST* be distinct integers, and there must be precisely one
! definition for each value of  exact_model  in  param.ccl .  It is not
! necessary that they be in numerical order, that is just to make things
! look prettier.
!
! N.b. We are counting on the preprocessor being case-sensitive here,
!      since these same names with EXACT changed to Exact, are names
!      of subroutines for the individual metric types!
!

! Minkowski spacetime
#define EXACT__Minkowski			1
#define EXACT__Minkowski_shift			2
#define EXACT__Minkowski_funny			3
#define EXACT__Minkowski_gauge_wave		4
#define EXACT__Minkowski_shifted_gauge_wave	5
#define EXACT__Minkowski_conf_wave		6

! black hole spacetimes
#define EXACT__Schwarzschild_EF			10
#define EXACT__Schwarzschild_PG			11
#define EXACT__Schwarzschild_BL			12
#define EXACT__Schwarzschild_Novikov		13
#define EXACT__Kerr_BoyerLindquist		14
#define EXACT__Kerr_KerrSchild			15
#define EXACT__Kerr_KerrSchild_spherical	16
#define EXACT__Schwarzschild_Lemaitre		17
#define EXACT__multi_BH				18
#define EXACT__Alvi				19
#define EXACT__Thorne_fakebinary		20

! cosmological spacetimes
#define EXACT__Lemaitre			50
!C this metric doesnt work, and has been moved to ./archive/
!C#define EXACT__Robertson_Walker	51
#define EXACT__de_Sitter		52
#define EXACT__de_Sitter_Lambda		53
#define EXACT__anti_de_Sitter_Lambda	54
#define EXACT__Bianchi_I		55
#define EXACT__Goedel			56
#define EXACT__Bertotti			57
#define EXACT__Kasner_like		58
#define EXACT__Kasner_axisymmetric	59
#define EXACT__Kasner_generalized	60
#define EXACT__Gowdy_wave		61
#define EXACT__Milne			62

! miscelaneous spacetimes
#define EXACT__boost_rotation_symmetric	80
#define EXACT__bowl			81
#define EXACT__constant_density_star	82

! math
#define EXACT__pi 3.14159265358979d0
