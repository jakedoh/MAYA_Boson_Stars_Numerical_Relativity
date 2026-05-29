#!/bin/bash

rm -rf EinsteinTensor/;

Kranclocal=Kranc/Bin/kranc
$Kranclocal EinsteinTensor.m &> EinsteinTensor.out

echo "----------------- Output ----------------------"
tail EinsteinTensor.out
echo "-----------------------------------------------"

if [ ! -r $thdir ]; then
   echo "Thorn generation failed.  See EinsteinTensor.out"
   exit
fi

echo "Now modifying output."
srcdir=EinsteinTensor/src
thdir=EinsteinTensor

# new versions of Kranc generate .cc files old ones .c files
if [ -r EinsteinTensor/src/RegisterMoL.c ] ; then
  cext=c
elif [ -r EinsteinTensor/src/RegisterMoL.cc ] ; then
  cext=cc
else
  echo >/dev/stderr "Cannot determine filename extension of Kranc generated source code."
  exit 1
fi



################################
## Additional initialization  ##
################################


####################################################
## Remove fake definitions of deltatime, etc.     ##
##   Get from previous timelevels of existing GFs ##
####################################################
srcfiles=( calcTimeDerivs.$cext )
for f in ${srcfiles[*]}
do
	sed '/deltatime =/ s/deltatime = 0/deltatime = CCTK_DELTA_TIME/' $srcdir/$f >$srcdir/$f.new
	mv $srcdir/$f.new $srcdir/$f

	sed '/KTMINUS2 =/ s/KTMINUS2 = 0/KTMINUS2 = K_p_p[index]/' $srcdir/$f >$srcdir/$f.new
	mv $srcdir/$f.new $srcdir/$f
	sed '/KTMINUS =/ s/KTMINUS = 0/KTMINUS = K_p[index]/' $srcdir/$f >$srcdir/$f.new
	mv $srcdir/$f.new $srcdir/$f

	sed '/A11TMINUS2 =/ s/A11TMINUS2 = 0/A11TMINUS2 = A11_p_p[index]/' $srcdir/$f >$srcdir/$f.new
	mv $srcdir/$f.new $srcdir/$f
	sed '/A11TMINUS =/ s/A11TMINUS = 0/A11TMINUS = A11_p[index]/' $srcdir/$f >$srcdir/$f.new
	mv $srcdir/$f.new $srcdir/$f

	sed '/A21TMINUS2 =/ s/A21TMINUS2 = 0/A21TMINUS2 = A21_p_p[index]/' $srcdir/$f >$srcdir/$f.new
	mv $srcdir/$f.new $srcdir/$f
	sed '/A21TMINUS =/ s/A21TMINUS = 0/A21TMINUS = A21_p[index]/' $srcdir/$f >$srcdir/$f.new
	mv $srcdir/$f.new $srcdir/$f

	sed '/A31TMINUS2 =/ s/A31TMINUS2 = 0/A31TMINUS2 = A31_p_p[index]/' $srcdir/$f >$srcdir/$f.new
	mv $srcdir/$f.new $srcdir/$f
	sed '/A31TMINUS =/ s/A31TMINUS = 0/A31TMINUS = A31_p[index]/' $srcdir/$f >$srcdir/$f.new
	mv $srcdir/$f.new $srcdir/$f

	sed '/A22TMINUS2 =/ s/A22TMINUS2 = 0/A22TMINUS2 = A22_p_p[index]/' $srcdir/$f >$srcdir/$f.new
	mv $srcdir/$f.new $srcdir/$f
	sed '/A22TMINUS =/ s/A22TMINUS = 0/A22TMINUS = A22_p[index]/' $srcdir/$f >$srcdir/$f.new
	mv $srcdir/$f.new $srcdir/$f

	sed '/A32TMINUS2 =/ s/A32TMINUS2 = 0/A32TMINUS2 = A32_p_p[index]/' $srcdir/$f >$srcdir/$f.new
	mv $srcdir/$f.new $srcdir/$f
	sed '/A32TMINUS =/ s/A32TMINUS = 0/A32TMINUS = A32_p[index]/' $srcdir/$f >$srcdir/$f.new
	mv $srcdir/$f.new $srcdir/$f

	sed '/A33TMINUS2 =/ s/A33TMINUS2 = 0/A33TMINUS2 = A33_p_p[index]/' $srcdir/$f >$srcdir/$f.new
	mv $srcdir/$f.new $srcdir/$f
	sed '/A33TMINUS =/ s/A33TMINUS = 0/A33TMINUS = A33_p[index]/' $srcdir/$f >$srcdir/$f.new
	mv $srcdir/$f.new $srcdir/$f

	sed '/Gam1TMINUS2 =/ s/Gam1TMINUS2 = 0/Gam1TMINUS2 = Gam1_p_p[index]/' $srcdir/$f >$srcdir/$f.new
	mv $srcdir/$f.new $srcdir/$f
	sed '/Gam1TMINUS =/ s/Gam1TMINUS = 0/Gam1TMINUS = Gam1_p[index]/' $srcdir/$f >$srcdir/$f.new
	mv $srcdir/$f.new $srcdir/$f

	sed '/Gam2TMINUS2 =/ s/Gam2TMINUS2 = 0/Gam2TMINUS2 = Gam2_p_p[index]/' $srcdir/$f >$srcdir/$f.new
	mv $srcdir/$f.new $srcdir/$f
	sed '/Gam2TMINUS =/ s/Gam2TMINUS = 0/Gam2TMINUS = Gam2_p[index]/' $srcdir/$f >$srcdir/$f.new
	mv $srcdir/$f.new $srcdir/$f

	sed '/Gam3TMINUS2 =/ s/Gam3TMINUS2 = 0/Gam3TMINUS2 = Gam3_p_p[index]/' $srcdir/$f >$srcdir/$f.new
	mv $srcdir/$f.new $srcdir/$f
	sed '/Gam3TMINUS =/ s/Gam3TMINUS = 0/Gam3TMINUS = Gam3_p[index]/' $srcdir/$f >$srcdir/$f.new
	mv $srcdir/$f.new $srcdir/$f

done

########################################
## Add assertion for _p_p's existence ##
########################################
f=calcTimeDerivs.$cext
sed "/Loop over the grid points/ i\\
  /* Check whether enough timelevels exist in current scope */\\
  assert( A11_p_p != NULL );\\
  assert( A11_p != NULL );\\
" $srcdir/$f >$srcdir/$f.new
mv $srcdir/$f.new $srcdir/$f

########################################
## Add scheduling of EinsteinTensor in INITIAL ##
########################################
f=schedule.ccl
cat > $thdir/$f.additional << EOF

# Additional schedulings

schedule group createEinsteinTensor at CCTK_POSTINITIAL after (MoL_PostStep MoL_PostInitial EinsteinTensor_time_derivs)
{
} "Create the Einstein Tensor purely from spacetime variables."

schedule group createEinsteinTensor at CCTK_POST_RECOVER_VARIABLES after ( MoL_PostStep EinsteinTensor_time_derivs )
{
} "Create the Einstein Tensor purely from spacetime variables."

schedule group createEinsteinTensor at CCTK_EVOL after ( MoL_Evolution MoL_PostStep EinsteinTensor_time_derivs )
{
} "Create the Einstein Tensor purely from spacetime variables."
EOF
cat $thdir/$f $thdir/$f.additional > $thdir/$f.new
mv $thdir/$f.new $thdir/$f; rm -f $thdir/$f.additional


##########################################################
## Move storage within a matterCheck if in schedule.ccl ##
##########################################################
f=schedule.ccl
sed '/STORAGE: jmatter_group/ d' $thdir/$f >$thdir/$f.new
mv $thdir/$f.new $thdir/$f
sed '/STORAGE: rhoPlusSmatter_group/ d' $thdir/$f >$thdir/$f.new
mv $thdir/$f.new $thdir/$f
sed '/STORAGE: rhoSmatter_group/ d' $thdir/$f >$thdir/$f.new
mv $thdir/$f.new $thdir/$f

sed "/schedule einstein_matter_initGFs/ i\\
\\
  STORAGE: jmatter_group[1]\\
  STORAGE: rhoPlusSmatter_group[1]\\
  STORAGE: rhoSmatter_group[1]\\
\\
" $thdir/$f >$thdir/$f.new
mv $thdir/$f.new $thdir/$f

sed '/./,/^$/!d' $thdir/$f >$thdir/$f.new
mv $thdir/$f.new $thdir/$f

##########################################################
## Put Boundary treatment of  into group ##
##########################################################
f=schedule.ccl
sed '/in MoL_PostStep/ s/in MoL_PostStep/in EinsteinTensor_Boundaries/' $thdir/$f >$thdir/$f.new
mv $thdir/$f.new $thdir/$f

cat > $thdir/$f.additional << EOF

# Schedule as desired

schedule group EinsteinTensor_Boundaries at CCTK_POSTINITIAL after (MoL_PostStep createEinsteinTensor)
{
} "Apply boundaries to "

schedule group EinsteinTensor_Boundaries at CCTK_POSTRESTRICTINITIAL after (MoL_PostStep createEinsteinTensor)
{
} "Apply boundaries to "

schedule group EinsteinTensor_Boundaries at CCTK_POST_RECOVER_VARIABLES after (MoL_PostStep createEinsteinTensor) before SetTmunu
{
} "Apply boundaries to "

schedule group EinsteinTensor_Boundaries at CCTK_POSTRESTRICT after (MoL_PostStep createEinsteinTensor)
{
} "Apply boundaries to "

schedule group EinsteinTensor_Boundaries at CCTK_EVOL after (MoL_Evolution createEinsteinTensor)
{
} "Apply boundaries to "

schedule group EinsteinTensor_Boundaries at CCTK_POSTREGRID after (MoL_PostStep createEinsteinTensor)
{
} "Apply boundaries to "

EOF
cat $thdir/$f $thdir/$f.additional > $thdir/$f.new
mv $thdir/$f.new $thdir/$f; rm -f $thdir/$f.additional

#########################################################
## Make ProjectEinsteinOnto restricted so other thorns ##
##  (e.g. SurfaceECs) can see what is contained in GF  ##
#########################################################
f=param.ccl
sed "/KEYWORD ProjectEinsteinOnto/ i\\
restricted: " $thdir/$f >$thdir/$f.new
mv $thdir/$f.new $thdir/$f

##################################
## Hardcode Outer Boundary type ##
##################################
f=Boundaries.$cext
sed "/ierr = 0/ a\\
\\
   ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 6, -3,\\
        \"EinsteinTensor::et00_group\", \"scalar\");\\
   if (ierr < 0)\\
          CCTK_WARN(-1, \"Failed to register BC for et00_group!\");\\
\\
   ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 6, -3,\\
        \"EinsteinTensor::et0j_group\", \"scalar\");\\
   if (ierr < 0)\\
          CCTK_WARN(-1, \"Failed to register BC for et0j_group!\");\\
\\
   ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 6, -3,\\
        \"EinsteinTensor::etij_group\", \"scalar\");\\
   if (ierr < 0)\\
          CCTK_WARN(-1, \"Failed to register BC for etij_group!\");\\
\\
   //ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 6, -3,\\
   //     \"EinsteinTensor::rhoPlusS_group\", \"scalar\");\\
   //if (ierr < 0)\\
   //       CCTK_WARN(-1, \"Failed to register BC for rhoPlusS_group!\");\\
\\
   //ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 6, -3,\\
   //     \"EinsteinTensor::jBSSN_group\", \"scalar\");\\
   //if (ierr < 0)\\
   //       CCTK_WARN(-1, \"Failed to register BC for jBSSN!\");\\
\\
   //ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 6, -3,\\
   //     \"EinsteinTensor::rhoSAijBSSN_group\", \"scalar\");\\
   //if (ierr < 0)\\
   //       CCTK_WARN(-1, \"Failed to register BC for rhoSAijBSSN!\");\\
\\
   ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 6, -3,\\
        \"EinsteinTensor::dtBSSN\", \"scalar\");\\
   if (ierr < 0)\\
          CCTK_WARN(-1, \"Failed to register BC for dtBSSN!\");\\
\\
   ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 6, -3,\\
        \"EinsteinTensor::dtcA_group\", \"scalar\");\\
   if (ierr < 0)\\
          CCTK_WARN(-1, \"Failed to register BC for dtcA_group!\");\\
\\
   ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 6, -3,\\
        \"EinsteinTensor::dtGam_group\", \"scalar\");\\
   if (ierr < 0)\\
          CCTK_WARN(-1, \"Failed to register BC for dtGam_group!\");\\
\\
  if ( CCTK_EQUALS(matterCheck,\"yes\") ) {\\
   ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 6, -3,\\
        \"EinsteinTensor::jmatter_group\", \"scalar\");\\
   if (ierr < 0)\\
          CCTK_WARN(-1, \"Failed to register BC for jmatter!\");\\
\\
   ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 6, -3,\\
        \"EinsteinTensor::rhoSmatter_group\", \"scalar\");\\
   if (ierr < 0)\\
          CCTK_WARN(-1, \"Failed to register BC for rhoSmatter!\");\\
\\
   ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 6, -3,\\
        \"EinsteinTensor::rhoPlusSmatter_group\", \"scalar\");\\
   if (ierr < 0)\\
          CCTK_WARN(-1, \"Failed to register BC for rhoPlusSmatter_group!\");\\
  } \\
\\
\\
" $srcdir/$f >$srcdir/$f.new
mv $srcdir/$f.new $srcdir/$f

##########################
## Generate README file ##
##########################
f=README
cat > $thdir/$f << EOF
Cactus Code Thorn $thdir 
Thorn Author(s)     : Tanja Bode <tanja.bode@physics.gatech.edu>
Thorn Maintainer(s) : Tanja Bode <tanja.bode@physics.gatech.edu>
--------------------------------------------------------------------------

Purpose of the thorn:

Creates an effective Einstein tensor from the spacetime variables only.
Can create checks for both the stress energy tensor in the form 
used in the spacetime evolution (S^i, rho(n), and Sij^TF) or construct
the Einstein Tensor with respect to t^a for direct comparisons to 
the stress energy tensor.

Some things to note:

  * Since this depends on numerical time derivatives from evolved
    quantities, it is initially junk.
  * It is also very noisy due to the second derivatives taken in 
    calculating.  Give it time to settle before you take anything 
    out of it.
  * Currently inherits from both ADMBase and Kranc2BSSNChiMatter.
    This is historical only.
  * Einstein tensor includes the factor of 8 pi so you can directly 
    compare with Tmunu.

EOF

############################################
## Copy thorn generation scripts to thorn ##
############################################
mkdir $thdir/m
cp EinsteinTensor.m $thdir/m
cp MakeEinsteinTensor.sh $thdir/m

###########################
echo "Done."
###########################
