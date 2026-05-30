#! /usr/bin/perl -i.bak
# $Header$

my $help_msg = <<'EOF';
This perl script converts old parameter files for this thorn to use
the new (June 2002) parameter names.

Usage:
	convert-pars.pl myparfile.par myparfile2.par ...

This will write new parameter files myparfile.par, myparfile2.par, ... ,
renaming the original files to myparfile.par.bak, myparfile2.par.bak, ... .

This script can be fooled by \-line-continues in strange places or
other lexical wierdness, but in practice it should work for almost all
"reasonable" parameter files.
EOF

#
# Assumptions:
# - We assume that the exactmodel/exact_model parameter is specified
#   *before* any of the parameters for the individual model.
# - We assume that all keywords/strings are in "double quotes".
#

################################################################################

if ( (scalar(@ARGV) == 1)
     && (($ARGV[0] eq '-help') || ($ARGV[0] eq '--help')) )
	{
	print $help_msg;
	exit;
	}

########################################

my $exact_model;

	while (<>)
	{
	# we only want to convert these strings if they're Exact:: parameters!
	# n.b. we *do* convert inside comments
	if (/Exact\s*::/i)
		{
		# thorn name
		s/^\s*exact\s*::/Exact::/gi;

		# exact_model
		s/exactmodel/exact_model/gi;
		if (/exact_model/)
			{
			s/minkowski/Minkowski/gi;
			s/flatshift/Minkowski\/shift/gi;
			s/flatfunny/Minkowski\/funny/gi;
			s/gauge_wave/Minkowski\/gauge wave/gi;
			s/finkelstein/Schwarzschild\/EF/gi;
			s/flatschwarz/Schwarzschild\/PG/gi;
			s/novikov/Schwarzschild\/Novikov/gi;
			s/kerr/Kerr\/Boyer-Lindquist/gi;
			s/kerrschild/Kerr\/Kerr-Schild/gi;
			s/cosschwarz/Schwarzschild-Lemaitre/gi;
			s/multiBH/multi-BH/gi;
			s/alvi/Alvi/gi;
			s/fakebinary/Thorne-fakebinary/gi;
			s/lemaitre/Lemaitre/gi;
			s/rob-wal/Robertson-Walker/gi;
			s/desitter/de Sitter/gi;
			s/cosdesitter/de Sitter+Lambda/gi;
			s/cosadesitter/anti-de Sitter+Lambda/gi;
			s/bianchiI/Bianchi I/gi;
			s/godel/Goedel/gi;
			s/cosbertotti/Bertotti/gi;
			s/kasner/Kasner-like/gi;
			s/kasner_axi/Kasner-axisymmetric/gi;
			s/kasner_gen/Kasner-generalized/gi;
			s/milne/Milne/gi;
			s/boostrot/boost-rotation symmetric/gi;
			s/bowl/bowl/gi;
			s/starschwarz/constant density star/gi;
			$exact_model = $_;
			}

		# parameters for Minkowski spacetime
		# with time-dependent shift vector
		s/flatshift_a/Minkowski_shift__amplitude/gi;
		s/flatshift_s/Minkowski_shift__sigma/gi;

		# parameters for Minkowski spacetime
		# with non-trivial spatial coordinates
		s/flatfunny_a/Minkowski_funny__amplitude/gi;
		s/flatfunny_s/Minkowski_funny__sigma/gi;

		# parameters for Minkowski spacetime in gauge-wave coordinates
		s/GW_H/Minkowski_gauge_wave__what_fn/gi;
		if (/Minkowski_gauge_wave__what_fn/i)
			{
			s/"gau"/"Gaussian"/gi;
			}
		s/GW_a/Minkowski_gauge_wave__amplitude/gi;
		s/GW_omega/Minkowski_gauge_wave__omega/gi;
		s/GW_diagonal/Minkowski_gauge_wave__diagonal/gi;
		s/GW_del/Minkowski_gauge_wave__lambda/gi;
		s/GW_phase_shift/Minkowski_gauge_wave__phase/gi;

		# parameters for Schwarzschild metric
		# in Eddington-Finkelstein coordinates
		if ($exact_model =~ /"Schwarzschild\/EF"/)
			{
			s/kerrschild_m/Schwarzschild_EF__mass/gi;
			s/kerrschild_eps/Schwarzschild_EF__epsilon/gi;
			}

		# parameters for Schwarzschild metric
		# in Painleve-Gustrand coordinates
		if ($exact_model =~ /"Schwarzschild\/PG"/)
			{
			s/kerrschild_m/Schwarzschild_PG__mass/gi;
			s/kerrschild_eps/Schwarzschild_PG__epsilon/gi;
			}

		# parameters for Schwarzschild metric in Novikov coordinates
		# there were no old parameters for this model

		# parameters for Schwarzschild-Lemaitre metric
		# (Schwarzschild black hole with cosmological constant)
		s/schwarz_lam/Schwarzschild_Lemaitre__Lambda/gi;
		s/schwarz_mas/Schwarzschild_Lemaitre__mass/gi;

		# parameters for Kerr metric in Boyer-Lindquist coordinates
		s/kerrc_m/Kerr_BoyerLindquist__mass/gi;
		s/kerrc_a/Kerr_BoyerLindquist__spin/gi;

		# parameters for Kerr metric in Kerr-Schild coordinates
		if ($exact_model =~ /"Kerr\/Kerr-Schild"/)
			{
			s/kerrschild_boostv/Kerr_KerrSchild__boost_v/gi;
			s/kerrschild_eps/Kerr_KerrSchild__epsilon/gi;
			s/kerrschild_m/Kerr_KerrSchild__mass/gi;
			s/kerrschild_a/Kerr_KerrSchild__spin/gi;
			}

		# parameters for Majumdar-Papapetrou or Kastor-Traschen
		# maximally-charged multiple black hole solution
		s/kt_nBH/multi_BH__nBH/gi;
		s/kt_hubble/multi_BH__Hubble/gi;
		s/m_bh1/multi_BH__mass1/gi;
		s/co_bh1x/multi_BH__x1/gi;
		s/co_bh1y/multi_BH__y1/gi;
		s/co_bh1z/multi_BH__z1/gi;
		s/m_bh2/multi_BH__mass2/gi;
		s/co_bh2x/multi_BH__x2/gi;
		s/co_bh2y/multi_BH__y2/gi;
		s/co_bh2z/multi_BH__z2/gi;
		s/m_bh3/multi_BH__mass3/gi;
		s/co_bh3x/multi_BH__x3/gi;
		s/co_bh3y/multi_BH__y3/gi;
		s/co_bh3z/multi_BH__z3/gi;
		s/m_bh4/multi_BH__mass4/gi;
		s/co_bh4x/multi_BH__x4/gi;
		s/co_bh4y/multi_BH__y4/gi;
		s/co_bh4z/multi_BH__z4/gi;

		# parameters for Alvi spacetime
		s/alvi_m1/Alvi__mass1/gi;
		s/alvi_m2/Alvi__mass2/gi;
		s/alvi_bb/Alvi__separation/gi;

		# parameters for Thorne's fake binary solution
		s/fakebinary_atype/Thorne_fakebinary__atype/gi;
		s/fakebinary_retarded/Thorne_fakebinary__retarded/gi;
		s/fakebinary_eps/Thorne_fakebinary__epsilon/gi;
		s/fakebinary_a0/Thorne_fakebinary__separation/gi;
		s/fakebinary_omega0/Thorne_fakebinary__Omega0/gi;
		s/fakebinary_m/Thorne_fakebinary__mass/gi;
		s/fakebinary_bround/Thorne_fakebinary__smoothing/gi;

		# parameters for Lemaitre-type spacetime
		s/lemaitre_k/Lemaitre__kappa/gi;
		s/lemaitre_l/Lemaitre__Lambda/gi;
		s/lemaitre_e/Lemaitre__epsilon0/gi;
		s/lemaitre_r/Lemaitre__R0/gi;

		# parameters for Robertson-Walker spacetime
		s/robson_a/Robertson_Walker__R0/gi;
		s/robson_b/Robertson_Walker__rho/gi;
		s/robson_k/Robertson_Walker__k/gi;
		s/robson_model/Robertson_Walker__pressure/gi;
		if (/Robertson_Walker__pressure/)
			{
			s/"pressure"/"true"/gi;
			s/"no_pressure"/"false"/gi;
			}

		# parameters for de Sitter spacetime
		if ($exact_model =~ /"de Sitter"/)
			{
			s/desitt_a/de_Sitter__scale/gi;
			}

		# parameters for de Sitter spacetime with cosmological constant
		if ($exact_model =~ /"de Sitter\+Lambda"/)
			{
			s/desitt_a/de_Sitter_Lambda__scale/gi;
			}

		if ($exact_model =~ /"anti\-de Sitter\+Lambda"/)
			{
			s/desitt_a/anti_de_Sitter_Lambda__scale/gi;
			}

		# parameters for Bianchi type I cosmology
		s/bia/Bianchi_I__scale/gi;

		# parameters for Goedel spacetime
		s/godel_a/Goedel__scale/gi;

		# parameters for Bertotti spacetime
		s/bertotti_lam/Bertotti__Lambda/gi;

		# parameters for Kasner-like spacetime
		s/kasner_q/Kasner_like__q/gi;

		# parameters for axisymmetric Kasner spacetime
		# there are no parameters

		# parameters for generalized Kasner spacetime
		s/kasner_p1/Kasner_generalized__p1/gi;
		s/kasner_p2/Kasner_generalized__p2/gi;

		# parameters for Milne spacetime
		# there are no parameters

		# parameters for boost-rotation symmetricmetric spacetime
		s/boostrotscale/boost_rotation_symmetric__scale/gi;
		s/boostrotstrength/boost_rotation_symmetric__amp/gi;
		s/boostrotsafedistance/boost_rotation_symmetric__min_d/gi;

		# parameters for bowl (bag-of-gold) spacetime (non-Einstein)
		s/bowl_type/bowl__shape/gi;
		if (/bowl__shape/)
			{
			s/"gauss"/"Gaussian"/gi;
			s/"fermi"/"Fermi"/gi;
			}
		s/bowl_evolve/bowl__evolve/gi;
		s/bowl_a/bowl__strength/gi;
		s/bowl_c/bowl__center/gi;
		s/bowl_s/bowl__sigma/gi;
		s/bowl_dx/bowl__x_scale/gi;
		s/bowl_dy/bowl__y_scale/gi;
		s/bowl_dz/bowl__z_scale/gi;
		s/bowl_t0/bowl__t0/gi;
		s/bowl_st/bowl__sigma_t/gi;

		# parameters for constant density (Schwarzschild) star
		s/starschwarz_m/constant_density_star__mass/gi;
		s/starschwarz_r/constant_density_star__radius/gi;
		}

	print;

	# don't want $exact_model to carry over from one par file to another
	if (eof) { $exact_model = undef(); }
	}
