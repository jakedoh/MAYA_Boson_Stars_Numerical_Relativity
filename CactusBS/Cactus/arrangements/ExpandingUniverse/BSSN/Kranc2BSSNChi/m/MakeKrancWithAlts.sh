#!/bin/bash

# Run it
../../Cactus/Kranc/Bin/kranc Kranc2BSSN.m

# new versions of Kranc generate .cc files old ones .c files
if [ -r Kranc2BSSN/src/RegisterMoL.c ] ; then
  cext=c
elif [ -r Kranc2BSSN/src/RegisterMoL.cc ] ; then
  cext=cc
else
  echo >/dev/stderr "Cannot determine filename extension of Kranc generated source code."
  exit 1
fi

#### Kranc2BSSN*Matter Alterations ####
for Chi in "" Chi ; do
chi=${Chi/C/c}

## Alter interface.ccl ##
gawk "
{print}
/USES FUNCTION MoLRegisterEvolved/ {
found++
print \"\"
print \"CCTK_INT FUNCTION MoLRegisterConstrainedGroup(CCTK_INT IN ConstrainedIndex)\"
print \"USES FUNCTION MoLRegisterConstrainedGroup\"
print \"\"
print \"CCTK_INT FUNCTION MoLRegisterSaveAndRestoreGroup(CCTK_INT IN SandRIndex)\"
print \"USES FUNCTION MoLRegisterSaveAndRestoreGroup\"
}
END {exit(found!=1)}
" Kranc2BSSN${Chi}Matter/interface.ccl > Kranc2BSSN${Chi}Matter/interface.ccl.new
if [ $? -ne 0 ] ; then 
        echo >/dev/stderr "Unable to patch 'Kranc2BSSN${Chi}Matter/interface.ccl'"
        exit 1
fi
mv Kranc2BSSN${Chi}Matter/interface.ccl.new Kranc2BSSN${Chi}Matter/interface.ccl

## Alter param.ccl ##
gawk "
/::MoL_Num_Evolved_Vars/ {
found++
print \"CCTK_INT Kranc2BSSN${Chi}Matter_MaxNumConstrainedVars \\\"Number of constrained variables used by this thorn\\\" ACCUMULATOR-BASE=MethodofLines::MoL_Num_Constrained_Vars\"
print \"{\"
print \"  16:16 :: \\\"Number of constrained variables used by this thorn\\\"\"
print \"} 16\"
print \"\"
print \"CCTK_INT Kranc2BSSN${Chi}Matter_MaxNumSaveAndRestoreVars \\\"Number of save and restore variables used by this thorn\\\" ACCUMULATOR-BASE=MethodofLines::MoL_Num_SaveAndRestore_Vars\"
print \"{\"
print \"  10:10 :: \\\"Number of save and restore variables used by this thorn\\\"\"
print \"} 10\"
print \"\"
print \"restricted:\"
}
{print}
/CCTK_INT MoL_Num_Evolved_Vars/ {
found++
print \"USES CCTK_INT MoL_Num_Constrained_Vars\"
print \"USES CCTK_INT MoL_Num_SaveAndRestore_Vars\"
}
END {exit(found!=2)}
" Kranc2BSSN${Chi}Matter/param.ccl > Kranc2BSSN${Chi}Matter/param.ccl.new
if [ $? -ne 0 ] ; then 
        echo >/dev/stderr "Unable to patch 'Kranc2BSSN${Chi}Matter/param.ccl'"
        exit 1
fi
mv Kranc2BSSN${Chi}Matter/param.ccl.new Kranc2BSSN${Chi}Matter/param.ccl
 
## Alter RegisterMoL.cc? ##
gawk "
/return/ {
found++
print \"\"
print \"  if (CCTK_EQUALS(couple_matter, \\\"yes\\\") || CCTK_EQUALS(recalculate_constraints, \\\"matter\\\")) {\"
print \"    ierr += MoLRegisterSaveAndRestoreGroup(CCTK_GroupIndex(\\\"TmunuBase::stress_energy_scalar\\\"));\"
print \"    ierr += MoLRegisterSaveAndRestoreGroup(CCTK_GroupIndex(\\\"TmunuBase::stress_energy_vector\\\"));\"
print \"    ierr += MoLRegisterSaveAndRestoreGroup(CCTK_GroupIndex(\\\"TmunuBase::stress_energy_tensor\\\"));\"
print \"  }\"
print \"\"
print \"  ierr += MoLRegisterConstrainedGroup(CCTK_GroupIndex(\\\"ADMBase::metric\\\"));\"
print \"  ierr += MoLRegisterConstrainedGroup(CCTK_GroupIndex(\\\"ADMBase::curv\\\"));\"
print \"  ierr += MoLRegisterConstrainedGroup(CCTK_GroupIndex(\\\"ADMBase::lapse\\\"));\"
print \"  ierr += MoLRegisterConstrainedGroup(CCTK_GroupIndex(\\\"ADMBase::shift\\\"));\"
}
{print}
END {exit(found!=1)}
" Kranc2BSSN${Chi}Matter/src/RegisterMoL.$cext > Kranc2BSSN${Chi}Matter/src/RegisterMoL.$cext.new
if [ $? -ne 0 ] ; then 
        echo >/dev/stderr "Unable to patch 'Kranc2BSSN${Chi}Matter/src/RegisterMoL.$cext'"
        exit 1
fi
mv Kranc2BSSN${Chi}Matter/src/RegisterMoL.$cext.new Kranc2BSSN${Chi}Matter/src/RegisterMoL.$cext

done

###################################################
### put header file for externalEtaBeta into src ##
###################################################
for dir in Kranc2BSSN Kranc2BSSNChi Kranc2BSSNMatter Kranc2BSSNChiMatter
do
cat > $dir/src/externaletabeta.h <<EOF 
const CCTK_REAL *externalEtaBeta = (CCTK_REAL *)CCTK_VarDataPtr(cctkGH, 0, "ExternalEtaBeta::etaBeta");
if(externalEtaBeta == NULL)
{
  CCTK_WARN(CCTK_WARN_COMPLAIN,"Could not get data pointer for ExternalEtaBeta::etaBeta. Is ExternalEtaBeta active?");
  /* NOTREACHED */
}
EOF
done


############################################
## Copy thorn generation scripts to thorn ##
############################################
possibly_generated_thorns=( Kranc2BSSN Kranc2BSSNChi Kranc2BSSNMatter Kranc2BSSNChiMatter )
for thdir in ${possibly_generated_thorns[*]}
do
   if [ -d $thdir ]; then
      mkdir -p $thdir/m
      cp -f Kranc2BSSN.m $thdir/m
      cp -f MakeKrancWithAlts.sh $thdir/m
      cp -f $KFile $thdir/m
   fi
done

###########################
echo "Done."
###########################
