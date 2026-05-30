#!/bin/bash

#Upgrade previously generated Kranc2BSSNChiMatter thorns 
#PATH=$PATH:Kranc/Bin
#kranc Kranc2BSSN.m

# Put altered thorn into different directory (for putting in Tanja arrangements)
thorntoupgrade=Kranc2BSSNChiMatter
thorndir=${thorntoupgrade}CVs

upgradeET=no

rm -rf ${thorndir}
cp -r ${thorntoupgrade} ${thorndir}

#######################
## configuration.ccl ##
#######################

	## Require the thorn ConstraintBdyMatter
	sed '
/GenericFD/ a\\
REQUIRE THORNS\\: ConstraintBdyMatter\\
' ${thorndir}/configuration.ccl > ${thorndir}/configuration.ccl.new
	mv ${thorndir}/configuration.ccl.new ${thorndir}/configuration.ccl

##################
## schedule.ccl ##
##################

	## More storage for constraints (and Einstein Tensor)
	sed '
/STORAGE: scalarconstraints/ s/\[1\]/\[3\]/g
/STORAGE: bssnmom/ s/\[1\]/\[3\]/g
/STORAGE: bssnCcons/ s/\[1\]/\[3\]/g
/STORAGE: EinsteinTensor/ s/\[1\]/\[3\]/g
' ${thorndir}/schedule.ccl > ${thorndir}/schedule.ccl.new
	mv ${thorndir}/schedule.ccl.new ${thorndir}/schedule.ccl

	## Reschedule calc_constraints at POSTINITIAL from ANALYSIS
	sed "
/bssnchimatter_calc_constraints_4th at CCTK_ANALYSIS/ s/at CCTK_ANALYSIS/at CCTK_POSTINITIAL after (MoL_PostStep MoL_PostInitial bssn_evolved_to_adm)/g
/bssnchimatter_calc_constraints_2nd at CCTK_ANALYSIS/ s/at CCTK_ANALYSIS/at CCTK_POSTINITIAL after (MoL_PostStep MoL_PostInitial bssn_evolved_to_adm)/g
" ${thorndir}/schedule.ccl > ${thorndir}/schedule.ccl.new
	mv ${thorndir}/schedule.ccl.new ${thorndir}/schedule.ccl

	## Add constraint calculations at POST_RECOVER_VARIABLES and EVOL 
	sed "
/schedule bssnchimatter_calc_constraints_4th/ a\\
  {\\
    LANG: C\\
    SYNC: bssnCcons_group\\
    SYNC: bssnmom_group\\
    SYNC: scalarconstraints\\
  } \"bssnchimatter_calc_constraints_4th at PostInitial\"\\
\\
  schedule bssnchimatter_calc_constraints_4th at CCTK_POST_RECOVER_VARIABLES after (MoL_PostStep MoL_PostInitial bssn_evolved_to_adm) as bssn_calc_constraints\\
  {\\
    LANG: C\\
    SYNC: bssnCcons_group\\
    SYNC: bssnmom_group\\
    SYNC: scalarconstraints\\
  } \"bssnchimatter_calc_constraints_4th after recovery\"\\
\\
  schedule bssnchimatter_calc_constraints_4th at CCTK_EVOL after (MoL_Evolution bssn_evolved_to_adm) as bssn_calc_constraints
" ${thorndir}/schedule.ccl > ${thorndir}/schedule.ccl.new
	mv ${thorndir}/schedule.ccl.new ${thorndir}/schedule.ccl

	sed "
/schedule bssnchimatter_calc_constraints_2nd/ a\\
  {\\
    LANG: C\\
    SYNC: bssnCcons_group\\
    SYNC: bssnmom_group\\
    SYNC: scalarconstraints\\
  } \"bssnchimatter_calc_constraints_2nd at PostInitial\"\\
\\
  schedule bssnchimatter_calc_constraints_2nd at CCTK_POST_RECOVER_VARIABLES after (MoL_PostStep MoL_PostInitial bssn_evolved_to_adm) as bssn_calc_constraints\\
  {\\
    LANG: C\\
    SYNC: bssnCcons_group\\
    SYNC: bssnmom_group\\
    SYNC: scalarconstraints\\
  } \"bssnchimatter_calc_constraints_2nd after recovery\"\\
\\
  schedule bssnchimatter_calc_constraints_2nd at CCTK_EVOL after (MoL_Evolution bssn_evolved_to_adm) as bssn_calc_constraints
" ${thorndir}/schedule.ccl > ${thorndir}/schedule.ccl.new
	mv ${thorndir}/schedule.ccl.new ${thorndir}/schedule.ccl

	## Add Einstein Tensor calculations in POSTINITIAL and POST_RECOVER_VARIABLES
if [ $upgradeET = "yes" ]; then
	sed "
/schedule bssnchimatter_calc_EinsteinTensor_full4th/ a\\
  {\\
    LANG: C\\
    SYNC: EinsteinTensor\\
  } \"bssnchimatter_calc_EinsteinTensor_full4th during evolution\"\\
\\
  schedule bssnchimatter_calc_EinsteinTensor_full4th at CCTK_POST_RECOVER_VARIABLES after (MoL_PostStep MoL_PostInitial bssn_evolved_to_adm) before bssn_calc_constraints as bssn_EinsteinTensor\\
  {\\
    LANG: C\\
    SYNC: EinsteinTensor\\
  } \"bssnchimatter_calc_EinsteinTensor_full4th after recovery\"\\
\\
  schedule bssnchimatter_calc_EinsteinTensor_full4th at CCTK_POSTINITIAL after (MoL_PostStep MoL_PostInitial bssn_evolved_to_adm) before bssn_calc_constraints as bssn_EinsteinTensor

" ${thorndir}/schedule.ccl > ${thorndir}/schedule.ccl.new
	mv ${thorndir}/schedule.ccl.new ${thorndir}/schedule.ccl

	sed "
/schedule bssnchimatter_calc_EinsteinTensor_Matter/ a\\
  {\\
    LANG: C\\
    SYNC: EinsteinTensor\\
  } \"bssnchimatter_calc_EinsteinTensor_Matter during evolution\"\\
\\
  schedule bssnchimatter_calc_EinsteinTensor_Matter at CCTK_POST_RECOVER_VARIABLES after (MoL_PostStep MoL_PostInitial bssn_evolved_to_adm bssn_EinsteinTensor) before bssn_calc_constraints as bssn_EinsteinTensor_Matter\\
  {\\
    LANG: C\\
    SYNC: EinsteinTensor\\
  } \"bssnchimatter_calc_EinsteinTensor_Matter after recovery\"\\
\\
  schedule bssnchimatter_calc_EinsteinTensor_Matter at CCTK_POSTINITIAL after (MoL_PostStep MoL_PostInitial bssn_evolved_to_adm bssn_EinsteinTensor) before bssn_calc_constraints as bssn_EinsteinTensor_Matter

" ${thorndir}/schedule.ccl > ${thorndir}/schedule.ccl.new
	mv ${thorndir}/schedule.ccl.new ${thorndir}/schedule.ccl
fi

###################
## interface.ccl ##
###################

	## Increase timelevels for constraints to 3
	sed '
s/scalarconstraints type=GF timelevels=1/scalarconstraints type=GF timelevels=3/g
s/bssnmom_group type=GF timelevels=1/bssnmom_group type=GF timelevels=3/g
s/bssnCcons_group type=GF timelevels=1/bssnCcons_group type=GF timelevels=3/g
' ${thorndir}/interface.ccl > ${thorndir}/interface.ccl.new
	mv ${thorndir}/interface.ccl.new ${thorndir}/interface.ccl

if [ $upgradeET = "yes" ]; then
	sed '
s/EinsteinTensor type=GF timelevels=1/EinsteinTensor type=GF timelevels=3/g
' ${thorndir}/interface.ccl > ${thorndir}/interface.ccl.new
	mv ${thorndir}/interface.ccl.new ${thorndir}/interface.ccl
fi


#######################
## src/RegisterMoL.c ##
#######################

	## Register constraints as constrained
	sed "
/return/ i\\
  \\
  ierr += MoLRegisterConstrainedGroup(CCTK_GroupIndex(\"Kranc2BSSNChiMatter\\:\\:scalarconstraints\"));\\
  ierr += MoLRegisterConstrainedGroup(CCTK_GroupIndex(\"Kranc2BSSNChiMatter\\:\\:bssnmom_group\"));\\
  ierr += MoLRegisterConstrainedGroup(CCTK_GroupIndex(\"Kranc2BSSNChiMatter\\:\\:bssnCcons_group\"));\\
" ${thorndir}/src/RegisterMoL.c > ${thorndir}/src/RegisterMoL.c.new
	mv ${thorndir}/src/RegisterMoL.c.new ${thorndir}/src/RegisterMoL.c

######################
## src/Boundaries.c ##
######################

	## Hard-code "flat" BCs for constraints
#	sed "
#/CCTK_INT ierr = 0/ a\\
#\\
#  /* Hard-code no outer boundary for constraints */\\
#  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, -1,\\
#                \"Kranc2BSSNChiMatter\\:\\:scalarconstraints\", \"zero\");\\
#  if (ierr < 0)\\
#     CCTK_WARN(0, \"Failed to register BCs for Kranc2BSSNChiMatter\\:\\:scalarconstraints!\");\\
#  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, -1,\\
#                \"Kranc2BSSNChiMatter\\:\\:bssnmom_group\", \"zero\"); \\
#  if (ierr < 0)\\
#     CCTK_WARN(0, \"Failed to register BCs for Kranc2BSSNChiMatter\\:\\:bssnmom_group!\"); \\
#  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, -1, \\
#                \"Kranc2BSSNChiMatter\\:\\:bssnCcons_group\", \"zero\");\\
#  if (ierr < 0)\\
#     CCTK_WARN(0, \"Failed to register BCs for Kranc2BSSNChiMatter\\:\\:bssnCcons_group!\");\\
#  /***********************************************/\\
#\\
#" ${thorndir}/src/Boundaries.c > ${thorndir}/src/Boundaries.c.new
#	mv ${thorndir}/src/Boundaries.c.new ${thorndir}/src/Boundaries.c

##################################
## src/bssnchimatter_init_gfs.c ##
##################################

	## Add check that ConstraintBdyMatter is active
	sed "
/GenericFD_LoopOverEverything/ i\\
\\
  if ( !CCTK_IsThornActive(\"ConstraintBdyMatter\")) {\\
     CCTK_WARN (CCTK_WARN_ABORT, \"This Kranc thorn requires ConstraintBdyMatter to be active.\");\\
  }\\
" ${thorndir}/src/bssnchimatter_init_gfs.c > ${thorndir}/src/bssnchimatter_init_gfs.c.new
	mv ${thorndir}/src/bssnchimatter_init_gfs.c.new ${thorndir}/src/bssnchimatter_init_gfs.c
