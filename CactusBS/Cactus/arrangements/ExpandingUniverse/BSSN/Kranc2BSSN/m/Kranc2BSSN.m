(* ::Package:: *)

Get["KrancThorn`"];

SetEnhancedTimes[False];
SetSourceLanguage["C"];

(**************************************************************************************)
(* variable types *)
(**************************************************************************************)

(* Treatment of dir needs to change for vectorisation! *)
$Assumptions=Element[dir, Integers]

(**************************************************************************************)
(* Derivatives *)
(**************************************************************************************)

derivatives =
{
  PDstandard2nd[i_] -> StandardCenteredDifferenceOperator[1,1,i],
  PDstandard2nd[i_, i_] -> StandardCenteredDifferenceOperator[2,1,i],
  PDstandard2nd[i_, j_] -> StandardCenteredDifferenceOperator[1,1,i] StandardCenteredDifferenceOperator[1,1,j],

  PDstandard4th[i_] -> StandardCenteredDifferenceOperator[1,2,i],
  PDstandard4th[i_, i_] -> StandardCenteredDifferenceOperator[2,2,i],
  PDstandard4th[i_, j_] -> StandardCenteredDifferenceOperator[1,2,i] StandardCenteredDifferenceOperator[1,2,j],

  PDstandard6th[i_] -> StandardCenteredDifferenceOperator[1,3,i],
  PDstandard6th[i_, i_] -> StandardCenteredDifferenceOperator[2,3,i],
  PDstandard6th[i_, j_] -> StandardCenteredDifferenceOperator[1,3,i] StandardCenteredDifferenceOperator[1,3,j],

  PDstandardNth[i_]    -> StandardCenteredDifferenceOperator[1,constraintsfdorder/2,i],
  PDstandardNth[i_,i_] -> StandardCenteredDifferenceOperator[2,constraintsfdorder/2,i],
  PDstandardNth[i_,j_] -> StandardCenteredDifferenceOperator[1,constraintsfdorder/2,i] *
                          StandardCenteredDifferenceOperator[1,constraintsfdorder/2,j],

  PDplus[i_] -> DPlus[i],
  PDminus[i_] -> DMinus[i],
  PDplus[i_,j_] -> DPlus[i] DPlus[j],

  PDonesided2nd[1] -> dir[1] (-shift[1]^(2 dir[1]) + 4 shift[1]^dir[1] - 3 )/(2 spacing[1]),
  PDonesided2nd[2] -> dir[2] (-shift[2]^(2 dir[2]) + 4 shift[2]^dir[2] - 3 )/(2 spacing[2]),
  PDonesided2nd[3] -> dir[3] (-shift[3]^(2 dir[3]) + 4 shift[3]^dir[3] - 3 )/(2 spacing[3]),

  PDlopsided4th[1] -> dir[1] (-3 shift[1]^(-1 dir[1]) - 10 + 18 shift[1]^(1 dir[1]) - 6 shift[1]^(2 dir[1]) + shift[1]^(3 dir[1])) / (12 spacing[1]),
  PDlopsided4th[2] -> dir[2] (-3 shift[2]^(-1 dir[2]) - 10 + 18 shift[2]^(1 dir[2]) - 6 shift[2]^(2 dir[2]) + shift[2]^(3 dir[2])) / (12 spacing[2]),
  PDlopsided4th[3] -> dir[3] (-3 shift[3]^(-1 dir[3]) - 10 + 18 shift[3]^(1 dir[3]) - 6 shift[3]^(2 dir[3]) + shift[3]^(3 dir[3])) / (12 spacing[3]),
(*  ExtrapBound[] -> -shift[dir]^(faceSign 4) + 4 shift[dir]^(faceSign 3) - 6 shift[dir]^(faceSign 2) + 4 shift[dir]^(faceSign 1)*)

  PDlopsided6th[1] -> dir[1] (2 shift[1]^(-2 dir[1]) -24 shift[1]^(-1 dir[1]) - 35 shift[1]^(0) + 80 shift[1]^( dir[1])  -30 shift[1]^( 2 dir[1])+ 8 shift[1]^(3 dir[1])- shift[1]^(4 dir[1])  ) / (60 spacing[1]),
  PDlopsided6th[2] -> dir[2] (2 shift[2]^(-2 dir[2]) -24 shift[2]^(-1 dir[2]) - 35 shift[2]^(0) + 80 shift[2]^( dir[2])  -30 shift[2]^( 2 dir[2])+ 8 shift[2]^(3 dir[2])- shift[2]^(4 dir[2])  ) / (60 spacing[2]),
  PDlopsided6th[3] -> dir[3] (2 shift[3]^(-2 dir[3]) -24 shift[3]^(-1 dir[3]) - 35 shift[3]^(0) + 80 shift[3]^( dir[3])  -30 shift[3]^( 2 dir[3])+ 8 shift[3]^(3 dir[3])- shift[3]^(4 dir[3])  ) / (60 spacing[3])
};

PDN = PDstandardNth;

(**************************************************************************************)
(* Tensors *)
(**************************************************************************************)

(* Register all the tensors that will be used with TensorTools *)
Map[DefineTensor, 
{
  h, hInv, phi, chi, A, K, alpha, Gam, beta, betat, R, Rphi, gamma, bssnmom,  bssnCcons,
  g, k, AInv, DDphi, Dphi, B, DDalpha, divA, gInv,
  gPhys, kPhys, n, dir,
  rho, S, Si, Sij, Tj,
  falpha
}];


(* Set the attributes of all tensor-type quantities for which the
   defaults are not suitable. *)

SetTensorAttribute[phi, TensorWeight, 1/6];
SetTensorAttribute[phi, TensorSpecial, "log"];
SetTensorAttribute[chi, TensorWeight, -2/3];
SetTensorAttribute[A, TensorWeight, -2/3];
SetTensorAttribute[h, TensorWeight, -2/3];
SetTensorAttribute[Gam, TensorWeight, 2/3];
SetTensorAttribute[Gam, TensorSpecial, "gamma"];

(* Register the TensorTools symmetries (this is very simplistic) *)
Map[AssertSymmetricIncreasing,
{
  hInv[ua,ub], AInv[ua,ub]
}];

Map[AssertSymmetricDecreasing, 
{
  h[la,lb], A[la,lb], R[la,lb], Rphi[la,lb],
  g[la,lb], k[la,lb], Sij[la,lb]
}];

AssertSymmetricDecreasing[gamma[ua,lb,lc], lb, lc];

(* Determinants of the metrics in terms of their components
   (Mathematica symbolic expressions) *)
hDet = Det[MatrixOfComponents[h[la,lb]]];
gDet = Det[MatrixOfComponents[g[la,lb]]];

Do[

phiMethod = OddQ[loopIndex];
addMatter = If[loopIndex >= 3, True, False];

fnPrefix = "bssn" <> If[phiMethod,"","chi"] <> If[addMatter,"matter",""];
thornName = "Kranc2BSSN" <> If[phiMethod,"","Chi"] <> If[addMatter,"Matter",""];

(**************************************************************************************)
(* Groups *)
(**************************************************************************************)

(* Cactus group definitions *)

evolvedGroups = Map[CreateGroupFromTensor, {h[li,lj], A[li,lj], If[phiMethod, phi, chi], K, Gam[ui], alpha, beta[ui], betat[ui]}];

admGroups = 
  {{"admbase::metric", {gxx,gxy,gxz,gyy,gyz,gzz}},
   {"admbase::curv", {kxx,kxy,kxz,kyy,kyz,kzz}},
   {"admbase::lapse", {alp}},
   {"admbase::shift", {betax,betay,betaz}}};

evaluatedGroups =
  {{"scalarconstraints", {bssnham, bssnT, bssnD}}, CreateGroupFromTensor[bssnmom[li]], CreateGroupFromTensor[bssnCcons[ui]]};
evaluatedGroups = Map[AddGroupExtra[#, Timelevels -> 3] &, evaluatedGroups];
(*evaluatedGroups = Map[AddGroupExtra[#, Timelevels -> cvTimelevels] &, evaluatedGroups];*)

stressenergyGroups =
  {{"TmunuBase::stress_energy_scalar", {eTtt} },
   {"TmunuBase::stress_energy_vector", {eTtx,eTty,eTtz} },
   {"TmunuBase::stress_energy_tensor", {eTxx,eTxy,eTxz,eTyy,eTyz,eTzz} } };

externalEtaBetaGroups =
   {{"ExternalEtaBeta::etaBeta_group", {externalEtaBeta} } };

declaredGroups = Join[evolvedGroups, evaluatedGroups];
declaredGroupNames = Map[First, declaredGroups];

If[addMatter, groups = Join[declaredGroups, admGroups, stressenergyGroups, externalEtaBetaGroups];, groups = Join[declaredGroups, admGroups, externalEtaBetaGroups];];

(**************************************************************************************)
(* Shorthands *)
(**************************************************************************************)

shorthands = 
{
  deth, invdeth, detg, invdetg, hInv[ua,ub], gamma[ua, lb, lc],
  R[la,lb], Rphi[la,lb], RS, RphiS, em4phi, e4phi, detgmthirdroot,
  twelfth, g[li,lj], k[li,lj], AInv[ua,ub], B[la,lb],
  DDalpha[li,lj], divA[ui], DDphi[li,lj], Dphi[li], gInv[ua,ub], 
  trB, trA, n[ui], dir[ui], faceSign, chiCutoff,
  rho, S, Si[la], Sij[la,lb], Tj[la], pi, muellerTop, muellerBottom, mueller, e2phi,
  expansion, Hubble, falpha , expansionPrime , HubblePrime
};

(**************************************************************************************)
(* Parameters *)
(**************************************************************************************)

realParameters = 
{
  {Name -> m,
   Description -> "Factor of momentum constraint added to Gam evolution. See Beyer and Sarbach 2010. (DEPRECATED -- Hard-coded to 1)", 
   Default -> 1, Steerable -> Always},
  {Name -> lapseAdvection, 
   Description -> "Gauge Parameter, All Lapse: dt_alp ~ lapseAdvection*( beta^i d_i(alp) ) + ...",
   Default -> 0, Steerable -> Always},
  {Name -> harmonicf, 
   Description -> "Gauge Parameter, 1+log lapse: d_t alp = - harmonicf*alp*K + ... ",
   Default -> 2, Steerable -> Always},
  {Name -> gammaDriverLambda, 
   Description -> "Gauge Parameter, NASA/Mueller Shifts: dt_beta^i ~ alp^gammaDriverLambda * betat^i + ...",
   Default -> 1, Steerable -> Always},
  {Name -> betaDotAlphaFactor, 
   Description -> "Gauge Parameter, NASA/Mueller Shifts: dt_beta^i = betaDotAlphaFactor*( alp^gammaDriverLambda betat^i ) + ...",
   Default -> 3/4, Steerable -> Always},
  {Name -> newNASAAdvection, 
   Description -> "Gauge Parameter, NASA/Mueller Shifts: dt_beta^i = newNASAAdvection*( beta^j d_j(beta^i) )",
   Default -> 0, Steerable -> Always},
  {Name -> etaBeta, 
   Description -> "Gauge Parameter, NASA/Mueller Shifts: dt_betat^i ~ etaBeta*betat^i + ...",
   Default -> 0, Steerable -> Always},
  {Name -> chiBeta, 
   Description -> "Gauge Parameter, NASA/Mueller Shifts: dt_betat^i ~ alp^(gammaDriverLapsePower)*chiBeta*d_t(Gam^i) + ...",
   Default -> 0, Steerable -> Always},
  {Name -> gammaDriverLapsePower, 
   Description -> "Gauge Parameter, NASA/Mueller Shifts: dt_betat^i ~ alp^(gammaDriverLapsePower)*chiBeta*d_t(Gam^i) + ...",
   Default -> 0, Steerable -> Always}, 
  {Name -> NASAAdvection, 
   Description -> "Gauge Parameter, NASA/Mueller Shifts: dt_betat^i ~ NASAAdvection*( beta^j*d_j(Gam^i) ) + ...",
   Default -> 0, Steerable -> Always}, 
  {Name -> betatAdvection, 
   Description -> "Gauge Parameter: NASA/Mueller Shifts: dt_betat^i = betatAdvection*( beta^j d_j(betat^i) ) + ...",
   Default -> 0, Steerable -> Always},
  {Name -> addYoTerm, 
   Description -> "Add the Yo term to the evolution of Gam. Required for BBH run stability!",
   Default -> 1,
   Steerable -> Always},
  {Name -> ShiftA,
   Description -> "Gauge Parameter, Mueller shift: muellerBottom ~ (1-chi^(ShiftA/2))^ShiftB", 
   Default -> 1, 
   Steerable -> Always},
  {Name -> Ro,
   Description -> "Gauge Parameter, FLRW Lapse", 
   Default -> 10, 
   Steerable -> Always},
  {Name -> sigma,
   Description -> "Gauge Parameter, FLRW Lapse", 
   Default -> 10, 
   Steerable -> Always},
  {Name -> ShiftB, 
   Description -> "Gauge Parameter, Mueller shift: muellerBottom ~ (1-chi^(ShiftA/2))^ShiftB",
   Default -> 2, Steerable -> Always},
  {Name -> a0,
   Description -> "FLRW parameter, Expansion of the universe at the start of the simulation",
   Default -> 1, Steerable -> Always},
  {Name -> H0,
   Description -> "FLRW parameter, Hubble parameter at the start of the simulation",
   Default -> 1, Steerable -> Always},
  {Name -> t0,
   Description -> "FLRW parameter, time of the start of the simulation",
   Default -> 1, Steerable -> Always}
};

intParameters =
{ 
  {Name -> constraintsfdorder,
   Description -> "Finite-differencing order used in calculating the BSSN constraints -- only used when recalculate_constraints =/= no.",
   Default -> 4, AllowedValues -> {2,4,6},
   Steerable -> Recover},
  {Name -> cvtimelevels,
   Description -> "Number of timelevels for constraint GFs (only valid after ET_2011_10 branches).", 
   Default -> 1, 
   Steerable -> Recover}
};

chiParameters = 
{
  {Name -> chiEps, 
   Description -> "Minimum chi allowed in most denominators.",
   Default -> 0, Steerable -> Always},
  {Name -> chiUniversalCutoff, 
   Description -> "Minimum value for chi enforced after updates",
   Default -> -1, Steerable -> Always}
};

If[Not[phiMethod], realParameters = Join[ realParameters, chiParameters ] ];

matterParameters = 
{
  {Name -> delayGREvolutiontoItn, 
   Description -> "Freeze spacetime until this iteration -- Must be a factor of every_coarse",
   Default -> 0}
};

If[addMatter, realParameters = Join[ realParameters, matterParameters ] ];

(* these are ignored by the thorn but have to be set to something different *)
(* than ADMbase's default of "static" to avoid unneccessary SYNCs from      *)
(* ADMBase's schedule.ccl                                                   *)
evolutionMethodParameters =
{
  {
    Name -> "ADMBase::evolution_method",
    Description -> "Extend ADMBase evolution method",
    AllowedValues -> {thornName}
  },
  {
    Name -> "ADMBase::lapse_evolution_method",
    Description -> "Extend ADMBase lapse evolution method",
    AllowedValues -> {thornName}
  },
  {
    Name -> "ADMBase::shift_evolution_method",
    Description -> "Extend ADMBase shift evolution method",
    AllowedValues -> {thornName}
  },
  {
    Name -> "ADMBase::dtlapse_evolution_method",
    Description -> "d_t(alpha) evolution",
    AllowedValues -> {thornName}
  },
  {
    Name -> "ADMBase::dtshift_evolution_method",
    Description -> "d^2_t(shift) evolution",
    AllowedValues -> {thornName}
  }
};

(**************************************************************************************)
(* Initialize all grid functions to a recognizable value. These should be overwritten *)
(* later. *)
(**************************************************************************************)

INITVALUE = 100;

initGFsCalc = 
{
  Name -> fnPrefix <> "_init_gfs",
  Schedule -> {"at BaseGrid as bssn_initgfs"},
  Equations -> 
  {
    h[li,lj] -> INITVALUE,
    A[li,lj] -> INITVALUE,
    If[phiMethod, phi, chi] -> INITVALUE,
    K -> INITVALUE,
    Gam[ui] -> INITVALUE,
    alpha -> INITVALUE,
    beta[ui] -> INITVALUE,
    betat[ui] -> INITVALUE
  }
};

initConstraintGFs =
{
  Name -> fnPrefix <> "_init_constraints_gfs",
  Schedule -> {"at BaseGrid as bssn_initgfs_constraints"},
  Equations ->
  {
    bssnham -> 0,
    bssnD -> 0,
    bssnT -> 0,
    bssnmom[li] -> 0,
    bssnCcons[ui] -> 0
  } 
};

(**************************************************************************************)
(* ADMBase <--> BSSN variable translation *)
(**************************************************************************************)

admToEvolvedCalc = 
{
  Name -> fnPrefix <> "_adm_to_evolved",
  Schedule -> {"at INITIAL after ADMBase_PostInitial as bssn_adm_to_evolved"},
  Shorthands -> {twelfth, detg, detgmthirdroot, invdetg, g[la,lb],k[la,lb], gInv[ua,ub],em4phi},
  Equations ->
    {
     g11 -> gxx, g21 -> gxy, g22 -> gyy, g31 -> gxz, g32 -> gyz, g33 -> gzz,
     twelfth -> 1/12,
     detg -> gDet,
     detgmthirdroot -> detg^(-1/3),
     invdetg -> 1 / detg,
     gInv[ua,ub] -> invdetg gDet MatrixInverse[g[ua,ub]],
     k11 -> kxx, k21 -> kxy, k22 -> kyy, k31 -> kxz, k32 -> kyz, k33 -> kzz,
     h[la,lb] -> g[la,lb] * detgmthirdroot,
     If[phiMethod, phi -> twelfth Log[detg], chi -> detg^(-1/3)],
     K -> k[la,lb] gInv[ua,ub],
     If[phiMethod, em4phi -> Exp[-4 phi], em4phi -> chi],
     A[la,lb] -> em4phi (k[la,lb] - 1/3 K g[la,lb]),
     alpha -> alp,
     beta1 -> betax, beta2 -> betay, beta3 -> betaz,
     betat1 -> 0, betat2 -> 0, betat3 -> 0,
     Gam[ui] -> 0 (* Really this wants to be populated via extrapolation later. Then
                     we can have two boundary points instead of three. *)

   }
};

(* This is in a separate calculation from the other variables as it
   needs to take finite differences and hence the loop is only in the
   interior. We calculate the Gams based on h rather than g for
   convenience. *)

gammaEvolvedCalc[fdOrder_, PD_, PDadvect_] := 
{
  Name -> fnPrefix <> "_calc_evolved_gammas_" <> fdOrder,
  ConditionalOnKeyword -> {"fd_order", fdOrder},
  Where -> Interior,
  Schedule -> 
  {"at INITIAL after ADMBase_InitialData after bssn_adm_to_evolved as bssn_calc_evolved_gammas"},
  Shorthands -> {deth, invdeth, hInv[ua,ub], gamma[ua,lb,lc], dir[ui]},

  Equations ->
  {
     dir[ui] -> Sign[beta[ui]],
     deth -> hDet, invdeth -> 1 / deth, 
     (* The invdeth hDet is a no-op that is necessary for accurate inverted matrices (according to Ian Hinder). 
	Whether this is necessary in newer (>5.0) versions of Mathematica has not been tested.  *)
     hInv[ua,ub] -> invdeth hDet MatrixInverse[h[ua,ub]],
     gamma[ua,lb,lc] -> 1/2 hInv[ua,ud] (PD[h[lb,ld], lc] + PD[h[lc,ld], lb] - PD[h[lb,lc],ld]),

     Gam[ui] -> hInv[uj,uk] gamma[ui,lj,lk]
  }
};

(* On the boundary we cannot calculate the Gams using a centered
   difference, so we use a one-sided difference pointing into the
   grid.  This is always used with a 2nd order accurate operator as a
   fourth order accurate version requires too large a stencil.*)

gammaEvolvedBoundCalc[PD_] := 
{
  Name -> fnPrefix <> "_calc_evolved_bound_gammas",
  Where -> Boundary,
  Schedule -> 
  {"at INITIAL after ADMBase_InitialData after bssn_adm_to_evolved"},
  Shorthands -> {n1,n2,n3,dir[ui],deth, invdeth, hInv[ua,ub],gamma[ua,lb,lc]},

  Equations ->
  {
    n1 -> -x/r, 
    n2 -> -y/r,
    n3 -> -z/r,
    dir[ui] -> Sign[n[ui]],

    deth -> hDet, invdeth -> 1 / deth, 
    hInv[ua,ub] -> invdeth hDet MatrixInverse[h[ua,ub]],
    gamma[ua,lb,lc] -> 1/2 hInv[ua,ud] (PD[h[lb,ld], lc] + PD[h[lc,ld], lb] - PD[h[lb,lc],ld]),

    Gam[ui] -> hInv[uj,uk] gamma[ui,lj,lk]
  }
};

(*
extrapolateGammasCalc[fdOrder_] := 
{
  Name -> "bssn_extrapolate_gammas_" <> fdOrder,
  ConditionalOnKeyword -> {"fd_order", fdOrder},
  Where -> Boundary,
  Schedule -> 
  {"at INITIAL after bssn_calc_evolved_gammas"},
  Shorthands -> shorthands,

  Equations ->
  {
    faceSign -> 2 face - 1,
    Gam[ui] -> ExtrapBound[Gam[ui]]
  }
}
*)
updateADMParam = 
{
  Name -> "update_admbase_variables",
  Description -> "When to update the ADMBase variables from the BSSN variables",
  Default -> "always",
  AllowedValues -> {"always", "never"},
  Steerable -> Recover
};

evolvedToADMInitialCalc = 
{
  Name -> fnPrefix <> "_evolved_to_adm_initial",
  Schedule -> {"AT Initial as bssn_evolved_to_adm_initial after bssn_adm_to_evolved before HydroBase_Prim2ConInitial"},
  Shorthands -> {e4phi, g[la,lb], k[la,lb]},
  ConditionalOnKeyword -> {"update_admbase_variables", "always"},
  Equations ->
    {
     e4phi -> If[phiMethod, Exp[4 phi], 1/Max[chi,chiEps]],
     g[la,lb] -> h[la,lb] e4phi,
     k[la,lb] -> e4phi A[la,lb] + 1/3 g[la,lb] K,
     gxx -> g11, gxy -> g21, gyy -> g22, gxz -> g31, gyz -> g32, gzz -> g33,
     kxx -> k11, kxy -> k21, kyy -> k22, kxz -> k31, kyz -> k32, kzz -> k33,
     alp -> alpha,
     betax -> beta1,
     betay -> beta2,
     betaz -> beta3
   }
};

evolvedToADMPoststepCalc = 
{
  Name -> fnPrefix <> "_evolved_to_adm",
  Schedule -> {"in MoL_PostStep as bssn_evolved_to_adm after " <> thornName <> "_ApplyBCs before HydroBase_PostStep"},
  Shorthands -> {e4phi, g[la,lb], k[la,lb]},
  ConditionalOnKeyword -> {"update_admbase_variables", "always"},
  Equations ->
    {
     e4phi -> If[phiMethod, Exp[4 phi], 1/Max[chi,chiEps]],
     g[la,lb] -> h[la,lb] e4phi,
     k[la,lb] -> e4phi A[la,lb] + 1/3 g[la,lb] K,
     gxx -> g11, gxy -> g21, gyy -> g22, gxz -> g31, gyz -> g32, gzz -> g33,
     kxx -> k11, kxy -> k21, kyy -> k22, kxz -> k31, kyz -> k32, kzz -> k33,
     alp -> alpha,
     betax -> beta1,
     betay -> beta2,
     betaz -> beta3
   }
};

evolvedToADMCalc = 
{
  Name -> fnPrefix <> "_evolved_to_adm",
  Schedule -> {"at CCTK_EVOL as bssn_evolved_to_adm after MoL_Evolution"},
  Shorthands -> {e4phi, g[la,lb], k[la,lb]},
  ConditionalOnKeyword -> {"update_admbase_variables", "always"},
  Equations ->
    {
     e4phi -> If[phiMethod, Exp[4 phi], 1/Max[chi,chiEps]],
     g[la,lb] -> h[la,lb] e4phi,
     k[la,lb] -> e4phi A[la,lb] + 1/3 g[la,lb] K,
     gxx -> g11, gxy -> g21, gyy -> g22, gxz -> g31, gyz -> g32, gzz -> g33,
     kxx -> k11, kxy -> k21, kyy -> k22, kxz -> k31, kyz -> k32, kzz -> k33,
     alp -> alpha,
     betax -> beta1,
     betay -> beta2,
     betaz -> beta3
   }
};


(**************************************************************************************)
(* BSSN evolution *)
(**************************************************************************************)

(* The evolution equations were originally copied from the
   presentation given in Phys.Rev. D70 (2004) 104004 (gr-qc/0406003).
   A couple of sign errors have been corrected.  CHECK: what about the
   twist / Yo terms? *)

evolveMetricCalc[fdOrder_, PD_, PDadvect_] := 
{
  Name -> fnPrefix <> "_evolve_metric_" <> fdOrder,
  Schedule -> {"in MoL_CalcRHS as bssn_evolve_metric"},
  Where -> Interior,
  ConditionalOnKeyword -> {"fd_order", fdOrder},
  Shorthands -> {dir[ui]},
  CollectList -> {hInv[ui,uj]},
  Equations -> 
  {
    (* Shift advection *)
    dir[ui] -> Sign[beta[ui]],

    If[phiMethod, 
      dot[phi] -> -alpha/6 K  + beta[ui] PDadvect[phi,li]  + 1/6 PD[beta[uj], lj],
      dot[chi] -> -4 chi * (-alpha/6 K   + 1/6 PD[beta[uj], lj]) + beta[ui] PDadvect[chi,li]],

    dot[h[li,lj]] -> -2 alpha A[li,lj] + beta[uk] PDadvect[h[li,lj], lk] + 
                     2 TTSymmetrize[h[lk,li] PD[beta[uk], lj], li,lj] - 2/3 h[li,lj] PD[beta[uk], lk]
  }
};


evolveNonMetricCalc[fdOrder_, PD_, PDadvect_] := 
{
  Name -> fnPrefix <> "_evolve_nonmetric_" <> fdOrder,
  Schedule -> {"in MoL_CalcRHS as bssn_evolve"},
  Where -> Interior,
  ConditionalOnKeyword -> {"fd_order", fdOrder},
  Shorthands -> {deth,invdeth,hInv[ua,ub],gamma[ua,lb,lc],AInv[ua,ub],em4phi,R[li,lj],
                 chiCutoff,DDphi[lj,li],Dphi[li],Rphi[li,lj],DDalpha[lj,li],B[li,lj],trB,
                 dir[ui]},
  CollectList -> {hInv[ui,uj]},
  Equations -> 
  {
    (* Shorthands *)

    deth -> hDet,
    invdeth -> 1 / deth,
    hInv[ua,ub] -> invdeth hDet MatrixInverse[h[ua,ub]],
    gamma[ua, lb, lc] -> 1/2 hInv[ua,ud] (PD[h[lb,ld], lc] + PD[h[lc,ld], lb] - PD[h[lb,lc],ld]),
    AInv[ua,ub] -> hInv[ua,ui] hInv[ub,uj] A[li,lj],

    em4phi -> If[phiMethod, Exp[-4 phi], chi],

    R[li,lj] -> -1/2 hInv[uk,ul] PD[h[li,lj],ll,lk] + TTSymmetrize[h[lk,li] PD[Gam[uk],lj], li,lj] 
                + TTSymmetrize[h[li,lc] gamma[uc,lj,lk] Gam[uk],li,lj]
                + hInv[ul,us](2 TTSymmetrize[gamma[uk,ll,li] h[lj,la] gamma[ua,lk,ls],li,lj]
                              + gamma[uk,li,ls] h[la,lk] gamma[ua,ll,lj]),

    chiCutoff -> If[phiMethod, 0, Max[chi,chiEps]],

    DDphi[lj,li] -> If[phiMethod, PD[phi,lj,li] - gamma[uk,li,lj] PD[phi,lk],
                                  -1/(4 chiCutoff) ( PD[chi,lj,li] - gamma[uk,li,lj] PD[chi,lk] - chiCutoff^-1 PD[chi,li] PD[chi,lj] )],

    Dphi[li] -> If[phiMethod, PD[phi,li], -1/(4 chiCutoff) PD[chi,li]],

    Rphi[li,lj] -> -2 DDphi[lj,li] - 2 h[li,lj] hInv[uk,ul] DDphi[lk,ll] + 4 Dphi[li] Dphi[lj] 
                   - 4 h[li,lj] hInv[uk,ua] Dphi[la] Dphi[lk],

    DDalpha[lj,li] -> PD[alpha,lj,li] - gamma[ul,li,lj] PD[alpha,ll],

    (* This is the bracket that is trace free in the dot[A] equation *)
    B[li,lj] -> alpha R[li,lj] + alpha Rphi[li,lj] - DDalpha[lj,li] 
                + (* This is - in the paper *) 4 TTSymmetrize[Dphi[li] PD[alpha,lj], li,lj],

    trB -> B[li,lj] hInv[ui,uj],


    (* Constraints are calculated here as it is cheaper than recalculating the Ricci 
       in ANALYSIS. Note that the two are not equivalent. *)
    bssnham -> em4phi hInv[ui,uj] (R[li,lj] + Rphi[li,lj]) - A[li,lj] AInv[ui,uj] + 2/3 K^2,

    (* Shift advection *)
    dir[ui] -> Sign[beta[ui]],

    (* Evolution equations *)

    dot[K] -> -em4phi(hInv[ui,uj] DDalpha[li,lj] + (* This is - in the paper *) 
              2 hInv[ui,uj] Dphi[li] PD[alpha,lj])  
              + alpha * (AInv[ui,uj] A[li,lj] + 1/3 K^2 ) + beta[uj] PDadvect[K,lj],

    dot[A[li,lj]] -> em4phi B[li,lj] - em4phi/3 h[li,lj] trB
                     + alpha K A[li,lj] - 2 alpha A[li,lk] A[lm,lj] hInv[um,uk]
                     + beta[uk] PDadvect[A[li,lj], lk] + 2 TTSymmetrize[A[lk,li] PD[beta[uk],lj],li,lj] 
                     - 2/3 A[li,lj] PD[beta[uk],lk],

    dot[Gam[ui]] -> (beta[uj] PDadvect[Gam[ui], lj]
         + hInv[uk,ul] PD[beta[ui],ll,lk] 
         + 1/3 hInv[ui,uj] PD[beta[uk],lj,lk] 
         - Gam[uj] PD[beta[ui],lj] + 2 Gam[ui] PD[beta[uj],lj] / 3 
         - addYoTerm 4/3 (Gam[ui] - hInv[ul,um] gamma[ui,ll,lm]) PD[beta[uj],lj]
         - 2 AInv[ui,uj] PD[alpha,lj]
         + 2 alpha * ( - 2/3 hInv[ui,uj] PD[K,lj]
         + (gamma[ui,lk,ll] AInv[uk,ul] + 6 AInv[ui,uj] Dphi[lj]))  )
  }
};

CoupleMatterParam =
{
  Name -> "couple_matter",
  Description -> "Triggers matter terms in evolving the spacetime",
  Default -> "yes",
  AllowedValues -> {"yes","no"},
  Steerable -> Recover
};

nullifyGREvol =
{
  Name -> fnPrefix <> "_nullify_GREvol",
  Schedule -> {"in MoL_RHSBoundaries after bssn_boundary as bssn_nullify"},
  (*ConditionalOnTextuals -> {"delayGREvolutiontoTime > 0 && cctk_time < delayGREvolutiontoTime"},*)
  ConditionalOnTextuals -> {"delayGREvolutiontoItn > 0 && (cctk_iteration-1) < delayGREvolutiontoItn"},
  Where -> Everywhere,
  Equations ->
  {

  dot[h[li,lj]] -> 0,
  If[phiMethod, 
	dot[phi] -> 0, 
	dot[chi] -> 0],
  dot[K] -> 0,
  dot[A[li,lj]] -> 0,
  dot[Gam[ui]] -> 0,

  dot[alpha] -> 0,
  dot[beta[ui]] -> 0,
  dot[betat[ui]] -> 0 

  }
};

evolveNonMetricCalcMat = 
{
  Name -> fnPrefix <> "_evolve_nonmetric_matter",
  Schedule -> {"in MoL_CalcRHS after bssn_evolve as bssn_evolve_matter"},
  Where -> Interior,
  ConditionalOnKeyword -> {"couple_matter", "yes"},
  Shorthands -> {deth,invdeth,hInv[ua,ub],em4phi, e4phi,pi, rho,Sij[la,lb],Tj[la],Si[la],S},
  Equations -> 
  {

     deth -> hDet,
     invdeth -> 1 / deth,
     hInv[ua,ub] -> invdeth hDet MatrixInverse[h[ua,ub]],
     em4phi -> If[phiMethod, Exp[-4 phi], chi],
     e4phi -> If[phiMethod, Exp[4 phi], 1/Max[chi,chiEps]],

     (* Quantities for matter in the RHS *)
     Sij11 -> eTxx, Sij21 -> eTxy, Sij31 -> eTxz, 
     Sij22 -> eTyy, Sij32 -> eTyz, Sij33 -> eTzz,
     Tj1 -> eTtx, Tj2 -> eTty, Tj3 -> eTtz,
     Si[la] -> (1/alpha) ( - Tj[la] + beta[ub] Sij[la,lb] ),
     S -> em4phi hInv[ua,ub] Sij[la,lb],
     rho -> (1/alpha)^2 ( eTtt - 2 beta[ua] Tj[la] + beta[ua] beta[ub] Sij[la,lb] ),
     pi -> 3.14159265358979323846,

     dot[K] -> dot[K] + 4 pi alpha * ( rho + S ),
     dot[A[li,lj]] -> dot[A[li,lj]] - 8 pi alpha * ( em4phi Sij[li,lj] - 1/3 h[li,lj] S),
     dot[Gam[ui]] -> dot[Gam[ui]] - 16 pi alpha hInv[ui,uj] Si[lj]

  }

};

(**************************************************************************************)
(* Lapse conditions *)
(**************************************************************************************)

lapseParam = 
{
  Name -> "lapse_condition",
  Description -> "Lapse evolution method",
  Default -> "constant",
  AllowedValues -> {"1 + log", "1 + log 2ndcentred", "1 + log full4th", "1 + log 4thcentred", "constant",
    "harmonic", "exact", "1 + log 6thcentred", "1 + log 6th" , "FLRW", "FLRW2", "FLRW3"},
  Steerable -> Recover
};

constantLapseCalc = 
{
  Name -> fnPrefix <> "_constant_lapse",
  Schedule -> {"in MoL_CalcRHS as bssn_lapse"},
  ConditionalOnKeyword -> {"lapse_condition", "constant"},
  Shorthands -> {dir[ui]},
  Where -> Interior,
  Equations -> 
  {
    (* Shift advection *)
    dir[ui] -> Sign[beta[ui]],

    (* WARNING: This will be second order if lapseAdvection != 0 *)
    dot[alpha] -> lapseAdvection beta[ua] PDonesided2nd[alpha,la]
  }
};

exactLapseCalc = 
{
  Name -> fnPrefix <> "_exact_lapse",
  Schedule -> {"in MoL_PostStep after Exact__gauge as bssn_lapse"},
  ConditionalOnKeyword -> {"lapse_condition", "exact"},
  Where -> Everywhere,
  Equations -> 
  {
    alpha -> alp
  }
};

harmonicLapseCalc = 
{
  Name -> fnPrefix <> "_harmonic_lapse",
  Schedule -> {"in MoL_CalcRHS as bssn_lapse"},
  Where -> Interior,
  ConditionalOnKeyword -> {"lapse_condition", "harmonic"},
  Shorthands -> {dir[ui]},
  Equations -> 
  {
    (* Shift advection *)
    dir[ui] -> Sign[beta[ui]],

    (* WARNING: This will be second order if lapseAdvection != 0 *)
    dot[alpha] -> -alpha^2 K + lapseAdvection beta[ua] PDonesided2nd[alpha,la]
  }
};

onePlusLogLapseCalc = 
{
  Name -> fnPrefix <> "_one_plus_log_lapse",
  Schedule -> {"in MoL_CalcRHS as bssn_lapse"},
  Where -> Interior,
  ConditionalOnKeyword -> {"lapse_condition", "1 + log"},
  Shorthands -> {dir[ui]},
  Equations -> 
  {
    (* Shift advection *)
    dir[ui] -> Sign[beta[ui]],

    (* WARNING: This will be second order if lapseAdvection != 0 *)
    dot[alpha] -> - harmonicf alpha K + lapseAdvection beta[ua] PDonesided2nd[alpha,la] 
  }
};

onePlusLogLapseFull4thCalc = 
{
  Name -> fnPrefix <> "_one_plus_log_lapse_full4th",
  Schedule -> {"in MoL_CalcRHS as bssn_lapse"},
  Where -> Interior,
  ConditionalOnKeyword -> {"lapse_condition", "1 + log full4th"},
  Shorthands -> {dir[ui]},
  Equations -> 
  {
    (* Shift advection *)
    dir[ui] -> Sign[beta[ui]],

    dot[alpha] -> - harmonicf alpha K + lapseAdvection beta[ua] PDlopsided4th[alpha,la] 
  }
};

onePlusLogLapse6thCalc = 
{
  Name -> fnPrefix <> "_one_plus_log_lapse_6th",
  Schedule -> {"in MoL_CalcRHS as bssn_lapse"},
  Where -> Interior,
  ConditionalOnKeyword -> {"lapse_condition", "1 + log 6th"},
  Shorthands -> {dir[ui]},
  Equations -> 
  {
    (* Shift advection *)
    dir[ui] -> Sign[beta[ui]],

    dot[alpha] -> - harmonicf alpha K + lapseAdvection beta[ua] PDlopsided6th[alpha,la] 
  }
};

onePlusLogLapse2ndCentredCalc = 
{
  Name -> fnPrefix <> "_one_plus_log_lapse_2ndcentred",
  Schedule -> {"in MoL_CalcRHS as bssn_lapse"},
  Where -> Interior,
  ConditionalOnKeyword -> {"lapse_condition", "1 + log 2ndcentred"},
  Shorthands -> {dir[ui]},
  Equations -> 
  {
    (* Shift advection *)
    dir[ui] -> Sign[beta[ui]],

    dot[alpha] -> - harmonicf alpha K + lapseAdvection beta[ua] PDstandard2nd[alpha,la] 
  }
};

onePlusLogLapse4thCentredCalc = 
{
  Name -> fnPrefix <> "_one_plus_log_lapse_4thcentred",
  Schedule -> {"in MoL_CalcRHS as bssn_lapse"},
  Where -> Interior,
  ConditionalOnKeyword -> {"lapse_condition", "1 + log 4thcentred"},
  Shorthands -> {dir[ui]},
  Equations -> 
  {
    (* Shift advection *)
    dir[ui] -> Sign[beta[ui]],

    dot[alpha] -> - harmonicf alpha K + lapseAdvection beta[ua] PDstandard4th[alpha,la] 
  }
};

onePlusLogLapse6thCentredCalc = 
{
  Name -> fnPrefix <> "_one_plus_log_lapse_6thcentred",
  Schedule -> {"in MoL_CalcRHS as bssn_lapse"},
  Where -> Interior,
  ConditionalOnKeyword -> {"lapse_condition", "1 + log 6thcentred"},
  Shorthands -> {dir[ui]},
  Equations -> 
  {
    (* Shift advection *)
    dir[ui] -> Sign[beta[ui]],

    dot[alpha] -> - harmonicf alpha K + lapseAdvection beta[ua] PDstandard6th[alpha,la] 
  }
};

FLRWLapseCalc = 
{
  Name -> fnPrefix <> "_FLRW_lapse",
  Schedule -> {"in MoL_CalcRHS as bssn_lapse"},
  Where -> Interior,
  ConditionalOnKeyword -> {"lapse_condition", "FLRW"},
  Shorthands -> {dir[ui],falpha},
  Equations -> 
  {
    (* Shift advection *)
    dir[ui] -> Sign[beta[ui]],
    falpha -> 2 (-0.5 Tanh[(r - Ro)/sigma] + 0.5 ) / alpha + (0.5 Tanh[(r - Ro)/sigma] + 0.5) / 3,
    dot[alpha] -> - falpha alpha^2 K + lapseAdvection beta[ua] PDlopsided4th[alpha,la]
  }
};

FLRWLapseCalc2 = 
{
  Name -> fnPrefix <> "_FLRW_lapse2",
  Schedule -> {"in MoL_CalcRHS as bssn_lapse"},
  Where -> Interior,
  ConditionalOnKeyword -> {"lapse_condition", "FLRW2"},
  Shorthands -> {dir[ui],falpha},
  Equations -> 
  {
    (* Shift advection *)
    dir[ui] -> Sign[beta[ui]],
    falpha -> 1 / 3,
    dot[alpha] -> - falpha alpha^2 K + lapseAdvection beta[ua] PDlopsided4th[alpha,la]
  }
};

FLRWLapseCalc3 = 
{
  Name -> fnPrefix <> "_FLRW_lapse3",
  Schedule -> {"in MoL_CalcRHS as bssn_lapse"},
  Where -> Interior,
  ConditionalOnKeyword -> {"lapse_condition", "FLRW3"},
  Shorthands -> {dir[ui],expansion,Hubble},
  Equations -> 
  {
    (* Shift advection *)
    dir[ui] -> Sign[beta[ui]],

    expansion -> a0 (t/t0)^(1),
    Hubble -> H0 expansion^(-2),

    dot[alpha] -> -2 (-0.5 Tanh[(r - Ro)/sigma] + 0.5 ) alpha (K + 3 Hubble) - (0.5 Tanh[(r - Ro)/sigma] + 0.5) alpha^2 K / 3 + lapseAdvection beta[ua] PDlopsided4th[alpha,la]
  }
};


(**************************************************************************************)
(* Shift conditions *)
(**************************************************************************************)

shiftParam = 
{
  Name -> "shift_condition",
  Description -> "Shift evolution method",
  Default -> "constant",
  AllowedValues -> {"constant", "NASA", "NASA4th", "NASAfull4th", "NASA6th", "NASA2ndcentred", "NASA4thcentred", "exact", "Mueller2nd", "Muellerfull4th", "Mueller6th", "Mueller4thcentred", "Mueller2ndcentred", "Mueller6thcentred"},
  Steerable -> Recover
};

constantShiftCalc =
{
  Name -> fnPrefix <> "_constant_shift",
  Schedule -> {"in MoL_CalcRHS as bssn_shift"},
  Where -> Interior,
  ConditionalOnKeyword -> {"shift_condition", "constant"},
  Equations -> 
  {
    dot[beta[ui]] -> 0,
    dot[betat[ui]] -> 0
  }
};

exactShiftCalc =
{
  Name -> fnPrefix <> "_exact_shift",
  Schedule -> {"in MoL_PostStep after Exact__gauge as bssn_shift"},
  Where -> Everywhere,
  ConditionalOnKeyword -> {"shift_condition", "exact"},
  Equations -> 
  {
    beta1 -> betax,
    beta2 -> betay,
    beta3 -> betaz
  }
};

nasaShiftCalc[fdOrder_, PD_, PDadvect_] := 
{
  Name -> fnPrefix <> "_nasa_shift_" <> fdOrder,
  Schedule -> {"in MoL_CalcRHS after bssn_evolve after bssn_evolve_matter as bssn_shift"},
  Where -> Interior,
  Shorthands -> {dir[ui]},
  ConditionalOnKeyword -> {"shift_condition", "NASA" <> fdOrder},
  Equations -> 
  {
    (* Some of these are needed in the computation of GamRHS below *)

    (* Shift advection *)
    dir[ui] -> Sign[beta[ui]],

    dot[beta[ui]] -> betaDotAlphaFactor alpha^gammaDriverLambda betat[ui] + 
      newNASAAdvection beta[ua] PDadvect[beta[ui],la],
    dot[betat[ui]] -> -etaBeta betat[ui] + chiBeta alpha^gammaDriverLapsePower dot[Gam[ui]] - NASAAdvection beta[ua] PDadvect[Gam[ui],la]  + betatAdvection beta[ua] PDadvect[betat[ui],la]
 
  }
};

muellerShiftCalc[fdOrder_, PD_, PDadvect_] :=
{
  Name -> fnPrefix <> "_mueller_shift_" <> fdOrder,
  Schedule -> {"in MoL_CalcRHS after bssn_evolve after bssn_evolve_matter as bssn_shift"},
  Where -> Interior,
  Shorthands -> {dir[ui],deth,invdeth,hInv[ua,ub],e2phi,muellerTop,muellerBottom,mueller},
  ConditionalOnKeyword -> {"shift_condition", "Mueller" <> fdOrder},
  Equations ->
  {
    (* Shift advection *)
    dir[ui] -> Sign[beta[ui]],
    deth -> hDet,
    invdeth -> 1 / deth,
    hInv[ua,ub] -> invdeth hDet MatrixInverse[h[ua,ub]],

    e2phi -> If[phiMethod, Exp[2 phi], 1],

    If[phiMethod,
      muellerTop -> hInv[ua,ub] PDadvect[phi,la] PDadvect[phi,lb],
      muellerTop -> hInv[ua,ub] PDadvect[chi,la] PDadvect[chi,lb]],
    If[phiMethod,
      muellerBottom -> 1/2 e2phi ( 1 - e2phi^(-ShiftA) )^ShiftB,
      muellerBottom -> 2 Sqrt[Max[chi,chiEps]] (1-Max[chi,chiEps]^ShiftA/2)^ShiftB],

    mueller -> etaBeta Sqrt[Abs[muellerTop]] / muellerBottom,

    dot[beta[ui]] -> betaDotAlphaFactor alpha^gammaDriverLambda betat[ui] +
      newNASAAdvection beta[ua] PDadvect[beta[ui],la],
    dot[betat[ui]] -> -mueller betat[ui] + chiBeta alpha^gammaDriverLapsePower dot[Gam[ui]] - NASAAdvection beta[ua] PDadvect[Gam[ui],la]  + betatAdvection beta[ua] PDadvect[betat[ui],la]
  }
};


(* Need to modify Kranc so that we can set multiple conditional
parameters on a calculation. *)

(**************************************************************************************)
(* Boundary conditions *)
(**************************************************************************************)

boundaryParam = 
{
  Name -> "boundary_condition",
  Description -> "My boundary conditions",
  Default -> "radiative",
  AllowedValues -> {"none", "radiative", "FLRWRadiative" , "FLRW", "FLRWConformal", "FLRWConformalRadiative"},
  Steerable -> Recover
};

boundaryCalc =
{
  Name -> fnPrefix <> "_boundary",
  Schedule -> {"in MoL_RHSBoundaries as bssn_boundary"},
  ConditionalOnKeyword -> {"boundary_condition", "radiative"},
  Where -> Boundary,
  Shorthands -> {n1,n2,n3,dir[ui]},
  Equations ->
  {
    n1 -> -x/r, 
    n2 -> -y/r,
    n3 -> -z/r,

    dir[ui] -> Sign[n[ui]],

    (* \partial_t u = - (u - u_0) / r - \partial_r u
       for u_0 some background solution. In this case, Minkowski in Cartesian 
       coordinates. *)
    
    (* WARNING: This is only second order accurate *)

    dot[h[li,lj]] -> -(h[li,lj] - Euc[li,lj]) / r + n[uk] PDonesided2nd[h[li,lj], lk],
    dot[A[li,lj]] -> -(A[li,lj]) / r + n[uk] PDonesided2nd[A[li,lj], lk],
    dot[K] -> -(K) / r + n[uk] PDonesided2nd[K, lk],
    dot[Gam[ui]] -> 0, (* This is conventional *)

    If[phiMethod,
      dot[phi] -> -(phi) / r + n[uk] PDonesided2nd[phi, lk],
      dot[chi] -> -(chi - 1) / r + n[uk] PDonesided2nd[chi, lk]],
      
    dot[alpha] -> -(alpha - 1) / r + n[uk] PDonesided2nd[alpha, lk], 
    dot[beta[ui]] -> -(beta[ui] - 0) / r + n[uk] PDonesided2nd[beta[ui], lk], 
    dot[betat[ui]] -> -(betat[ui] - 0) / r + n[uk] PDonesided2nd[betat[ui], lk],

    (* clear constraints in the boundary region to get rid of poison *)

    bssnham -> 0, bssnD -> 0, bssnT -> 0,
    bssnmom1 -> 0, bssnmom2 -> 0, bssnmom3 -> 0
  }
};

boundaryCalcFLRWRadiative =
{
  Name -> fnPrefix <> "_boundaryFLRWRadiative",
  Schedule -> {"in MoL_RHSBoundaries as bssn_boundary"},
  ConditionalOnKeyword -> {"boundary_condition", "FLRWRadiative"},
  Where -> Boundary,
  Shorthands -> {n1,n2,n3,dir[ui], expansion , Hubble , expansionPrime , HubblePrime},
  Equations ->
  {
    n1 -> -x/r, 
    n2 -> -y/r,
    n3 -> -z/r,

    expansion -> a0 (t/t0)^(1.0/2.0),
    Hubble -> H0 t0/t,
    expansionPrime -> a0 (1.0/(2.0 t0)) (t/t0)^(-1.0/2.0),
    HubblePrime -> -H0 t0/(t^2),

    dir[ui] -> Sign[n[ui]],

    (* \partial_t u = \partial_t u_0 - (u - u_0) / r - \partial_r u
       for u_0 some background solution. In this case, FLRW in Cartesian 
       coordinates. *)
    
    (* WARNING: This is only second order accurate *)

    dot[h[li,lj]] -> -(h[li,lj] - Euc[li,lj]) / r + n[uk] PDonesided2nd[h[li,lj], lk],
    dot[A[li,lj]] -> -(A[li,lj] - 0) / r + n[uk] PDonesided2nd[A[li,lj], lk],
    dot[K] -> -3 HubblePrime -(K + 3 Hubble) / r + n[uk] PDonesided2nd[K, lk],
    dot[Gam[ui]] -> 0, (* This is conventional *)

    If[phiMethod,
      dot[phi] -> -expansionPrime/(2 expansion) -(phi + Log[expansion]/2) / r + n[uk] PDonesided2nd[phi, lk],
      dot[chi] -> -2 expansionPrime/(expansion^3) -(chi - expansion^(-2)) / r + n[uk] PDonesided2nd[chi, lk]],
      
    dot[alpha] -> expansionPrime -(alpha - expansion) / r + n[uk] PDonesided2nd[alpha, lk], 
    dot[beta[ui]] -> -(beta[ui] - 0) / r + n[uk] PDonesided2nd[beta[ui], lk], 
    dot[betat[ui]] -> -(betat[ui] - 0) / r + n[uk] PDonesided2nd[betat[ui], lk],

    (* clear constraints in the boundary region to get rid of poison *)

    bssnham -> 0, bssnD -> 0, bssnT -> 0,
    bssnmom1 -> 0, bssnmom2 -> 0, bssnmom3 -> 0
  }
};


boundaryCalcFLRWConformalRadiative =
{
  Name -> fnPrefix <> "_boundaryFLRWRadiativeConformal",
  Schedule -> {"in MoL_RHSBoundaries as bssn_boundary"},
  ConditionalOnKeyword -> {"boundary_condition", "FLRWConformalRadiative"},
  Where -> Boundary,
  Shorthands -> {n1,n2,n3,dir[ui], expansion , Hubble , expansionPrime , HubblePrime},
  Equations ->
  {
    n1 -> -x/r, 
    n2 -> -y/r,
    n3 -> -z/r,

    expansion -> a0 (t/t0)^(1),
    expansionPrime -> a0/t0,
    Hubble -> H0 expansion^(-2),
    HubblePrime -> -2.0 H0 expansionPrime expansion^(-3.0),
    

    dir[ui] -> Sign[n[ui]],

    (* \partial_t u =  \partial_t u_0 - (u - u_0) / r - \partial_r u
       for u_0 some background solution. In this case, FLRW in Cartesian 
       coordinates. *)
    
    (* WARNING: This is only second order accurate *)

    dot[h[li,lj]] -> -(h[li,lj] - Euc[li,lj]) / r + n[uk] PDonesided2nd[h[li,lj], lk],
    dot[A[li,lj]] -> -(A[li,lj] - 0) / r + n[uk] PDonesided2nd[A[li,lj], lk],
    dot[K] -> -3*HubblePrime -(K + 3*Hubble) / r + n[uk] PDonesided2nd[K, lk],
    dot[Gam[ui]] -> 0, (* This is conventional *)

    If[phiMethod,
      dot[phi] -> -expansionPrime/(expansion*2) -(phi + Log[expansion]/2) / r + n[uk] PDonesided2nd[phi, lk],
      dot[chi] -> -2*expansionPrime/(expansion^3) -(chi - expansion^(-2)) / r + n[uk] PDonesided2nd[chi, lk]],
      
    dot[alpha] -> expansionPrime -(alpha - expansion) / r + n[uk] PDonesided2nd[alpha, lk], 
    dot[beta[ui]] -> -(beta[ui] - 0) / r + n[uk] PDonesided2nd[beta[ui], lk], 
    dot[betat[ui]] -> -(betat[ui] - 0) / r + n[uk] PDonesided2nd[betat[ui], lk],

    (* clear constraints in the boundary region to get rid of poison *)

    bssnham -> 0, bssnD -> 0, bssnT -> 0,
    bssnmom1 -> 0, bssnmom2 -> 0, bssnmom3 -> 0
  }
};



boundaryCalcFLRW =
{
  Name -> fnPrefix <> "_boundaryFLRW",
  Schedule -> {"in MoL_PostStep after " <> thornName <> "_ApplyBCs before bssn_evolved_to_adm as bssn_boundary"},
  ConditionalOnKeyword -> {"boundary_condition", "FLRW"},
  Where -> Boundary,
  Shorthands -> {n1,n2,n3,dir[ui], expansion , Hubble},
  Equations ->
  {
    n1 -> -x/r, 
    n2 -> -y/r,
    n3 -> -z/r,

    expansion -> a0 (t/t0)^(1/2),
    Hubble -> H0 t0/t,

    dir[ui] -> Sign[n[ui]],

    h[li,lj] -> Euc[li,lj],
    A[li,lj] -> 0,
    K -> -3 Hubble,
    dot[Gam[ui]] -> 0,

    If[phiMethod,
      phi -> -Log[expansion]/2,
      chi -> expansion^(-2)],
      
    alpha -> expansion, 
    beta[ui] -> 0, 
    betat[ui] -> 0,

    (* clear constraints in the boundary region to get rid of poison *)

    bssnham -> 0, bssnD -> 0, bssnT -> 0,
    bssnmom1 -> 0, bssnmom2 -> 0, bssnmom3 -> 0
  }
};


boundaryCalcFLRWConformal =
{
  Name -> fnPrefix <> "_boundaryFLRWConformal",
  Schedule -> {"in MoL_PostStep after " <> thornName <> "_ApplyBCs before bssn_evolved_to_adm as bssn_boundary"},
  ConditionalOnKeyword -> {"boundary_condition", "FLRWConformal"},
  Where -> Boundary,
  Shorthands -> {n1,n2,n3,dir[ui], expansion , Hubble},
  Equations ->
  {
    n1 -> -x/r, 
    n2 -> -y/r,
    n3 -> -z/r,

    expansion -> a0 (t/t0)^(1),
    Hubble -> H0 expansion^(-2),

    dir[ui] -> Sign[n[ui]],

    h[li,lj] -> Euc[li,lj],
    A[li,lj] -> 0,
    K -> -3 Hubble,
    dot[Gam[ui]] -> 0,

    If[phiMethod,
      phi -> -Log[expansion]/2,
      chi -> expansion^(-2)],
      
    alpha -> expansion, 
    beta[ui] -> 0, 
    betat[ui] -> 0,

    (* clear constraints in the boundary region to get rid of poison *)

    bssnham -> 0, bssnD -> 0, bssnT -> 0,
    bssnmom1 -> 0, bssnmom2 -> 0, bssnmom3 -> 0
  }
};



(**************************************************************************************)
(* Project the BSSN algebraic constraints *)
(**************************************************************************************)

algconsParam = 
{
  Name -> "algebraic_constraint_projection",
  Description -> "Enforce algebraic constraints: chiUniversalCutoff and tracelessness of A_ij",
  Default -> "psu",
  AllowedValues -> {"none", "psu", "psuboundary"},
  Steerable -> Recover
};

psuProjectConstraintsCalc = 
{
  Name -> fnPrefix <> "_project_constraints_psu",
  Schedule -> {"in MoL_PostStep as bssn_project_constraints"},
  ConditionalOnKeyword -> {"algebraic_constraint_projection", "psu"},
  Shorthands -> {invdeth,hInv[ua,ub],trA},
  Equations ->
  {
    invdeth -> 1 / hDet,
    hInv[ua,ub] -> invdeth hDet MatrixInverse[h[ua,ub]],
    trA ->  A[lm,ln] hInv[um,un],

    If[phiMethod, invdeth -> 1 / hDet , chi -> Max[chi, chiUniversalCutoff]],

    A[la,lb] -> A[la,lb] - 1/3 h[la,lb] trA
  }
};

psuBoundaryProjectConstraintsCalc = 
{
  Name -> fnPrefix <> "_project_constraints_psu_boundary",
  Schedule -> {"in MoL_PostStep as bssn_project_constraints"},
  Where -> Boundary,
  ConditionalOnKeyword -> {"algebraic_constraint_projection", "psuboundary"},
  Shorthands -> {invdeth,hInv[ua,ub],trA},
  Equations ->
  {
    invdeth -> 1 / hDet,
    hInv[ua,ub] -> invdeth hDet MatrixInverse[h[ua,ub]],
    trA ->  A[lm,ln] hInv[um,un],

    A[la,lb] -> A[la,lb] - 1/3 h[la,lb] trA
  }
};

(**************************************************************************************)
(* Constraints calc *)
(**************************************************************************************)

recalcConstraintsParam = 
{
  Name -> "recalculate_constraints",
  Description -> "Recalculate the constraints from BSSN variables. Choose between vacuum or matter. FD order given by parameter constraintsfdorder.",
  Default -> "no",
  AllowedValues -> {"vacuum","matter", "no"},
  Steerable -> Recover
};

vacuumconstraintsCalc = 
{
  Name -> fnPrefix <> "_calc_constraints_vacuum",
  Schedule -> Automatic,
  After -> "ADMBase_SetADMVars",
  Where -> Interior,
  ConditionalOnKeyword -> {"recalculate_constraints", "vacuum"},
  Shorthands -> {deth,invdeth,dir[ui],hInv[ua,ub],gamma[ua,lb,lc],AInv[ua,ub],em4phi,chiCutoff,
                 DDphi[lj,li],Dphi[li],R[li,lj],Rphi[li,lj]},
  Equations ->
  {
        deth -> hDet,
        invdeth -> 1 / deth,
        dir[ui] -> Sign[beta[ui]],
        hInv[ua,ub] -> invdeth hDet MatrixInverse[h[ua,ub]],
        gamma[ua, lb, lc] -> 1/2 hInv[ua,ud] (PDN[h[lb,ld], lc] + PDN[h[lc,ld], lb] - PDN[h[lb,lc],ld]),
        AInv[ua,ub] -> hInv[ua,ui] hInv[ub,uj] A[li,lj],
        em4phi -> If[phiMethod, Exp[-4 phi], chi],
        chiCutoff -> If[phiMethod, 0, Max[chi,chiEps]],


        DDphi[lj,li] -> If[phiMethod, PDN[phi,lj,li] - gamma[uk,li,lj] PDN[phi,lk],
                                      -1/(4 chiCutoff) ( PDN[chi,lj,li] - gamma[uk,li,lj] PDN[chi,lk] 
                                                   - chiCutoff^-1 PDN[chi,li] PDN[chi,lj] )],
        Dphi[li] -> If[phiMethod, PDN[phi,li], -1/(4 chiCutoff) PDN[chi,li]],

        R[li,lj] -> -1/2 hInv[uk,ul] PDN[h[li,lj],ll,lk] + TTSymmetrize[h[lk,li] PDN[Gam[uk],lj], li,lj] 
                    + TTSymmetrize[h[li,lc] gamma[uc,lj,lk] Gam[uk],li,lj]
                    + hInv[ul,us](2 TTSymmetrize[gamma[uk,ll,li] h[lj,la] gamma[ua,lk,ls],li,lj]
                                  + gamma[uk,li,ls] h[la,lk] gamma[ua,ll,lj]),

        Rphi[li,lj] -> -2 DDphi[lj,li] - 2 h[li,lj] hInv[uk,ul] DDphi[lk,ll] + 4 Dphi[li] Dphi[lj] 
                       - 4 h[li,lj] hInv[uk,ua] Dphi[la] Dphi[lk],


        bssnham -> em4phi hInv[ui,uj] (R[li,lj] + Rphi[li,lj]) - A[li,lj] AInv[ui,uj] + 2/3 K^2,
        bssnmom[li] -> 6 hInv[uj,uk] A[lk,li] Dphi[lj]
                      + hInv[um,uk] (PDN[A[lm,li],lk]- gamma[ud,lm,lk] A[ld,li] 
                                                    - gamma[ud,li,lk] A[lm,ld]) 
                      -2/3 PDN[K,li],
	bssnCcons[ui] -> Gam[ui]  - hInv[uj,uk] gamma[ui,lj,lk],
        bssnD -> Log[deth],
        bssnT -> hInv[ui,uj] A[li,lj]
  }
};

matterconstraintsCalc = 
{
  Name -> fnPrefix <> "_calc_constraints_matter",
  Schedule -> Automatic,
  After -> "ADMBase_SetADMVars",
  Where -> Interior, 
  ConditionalOnKeyword -> {"recalculate_constraints", "matter"},
  Shorthands -> {Sij[la,lb],Tj[la],Si[la],rho,pi,deth,invdeth,dir[ui],hInv[ua,ub],
                 gamma[ua,lb,lc],AInv[ua,ub],em4phi,chiCutoff,DDphi[lj,li],Dphi[li],
                 R[li,lj],Rphi[li,lj]},
  Equations ->
  {

     (* Quantities for matter in the RHS *)
     Sij11 -> eTxx, Sij21 -> eTxy, Sij31 -> eTxz, 
     Sij22 -> eTyy, Sij32 -> eTyz, Sij33 -> eTzz,
     Tj1 -> eTtx, Tj2 -> eTty, Tj3 -> eTtz,
     Si[la] -> (1/alpha) ( - Tj[la] + beta[ub] Sij[la,lb] ),
     rho -> (1/alpha)^2 ( eTtt - 2 beta[ua] Tj[la] + beta[ua] beta[ub] Sij[la,lb] ),
     pi -> 3.14159265358979323846,

     deth -> hDet,
     invdeth -> 1 / deth,
     dir[ui] -> Sign[beta[ui]],
     hInv[ua,ub] -> invdeth hDet MatrixInverse[h[ua,ub]],
     gamma[ua, lb, lc] -> 1/2 hInv[ua,ud] (PDN[h[lb,ld], lc] + PDN[h[lc,ld], lb] - PDN[h[lb,lc],ld]),
     AInv[ua,ub] -> hInv[ua,ui] hInv[ub,uj] A[li,lj],
     em4phi -> If[phiMethod, Exp[-4 phi], chi],
     chiCutoff -> If[phiMethod, 0, Max[chi,chiEps]],


     DDphi[lj,li] -> If[phiMethod, PDN[phi,lj,li] - gamma[uk,li,lj] PDN[phi,lk],
                                   -1/(4 chiCutoff) ( PDN[chi,lj,li] - gamma[uk,li,lj] PDN[chi,lk]
                                                - chiCutoff^-1 PDN[chi,li] PDN[chi,lj] )],
     Dphi[li] -> If[phiMethod, PDN[phi,li], -1/(4 chiCutoff) PDN[chi,li]],

     R[li,lj] -> -1/2 hInv[uk,ul] PDN[h[li,lj],ll,lk] + TTSymmetrize[h[lk,li] PDN[Gam[uk],lj], li,lj]
                 + TTSymmetrize[h[li,lc] gamma[uc,lj,lk] Gam[uk],li,lj]
                 + hInv[ul,us](2 TTSymmetrize[gamma[uk,ll,li] h[lj,la] gamma[ua,lk,ls],li,lj]
                               + gamma[uk,li,ls] h[la,lk] gamma[ua,ll,lj]),

     Rphi[li,lj] -> -2 DDphi[lj,li] - 2 h[li,lj] hInv[uk,ul] DDphi[lk,ll] + 4 Dphi[li] Dphi[lj]
                    - 4 h[li,lj] hInv[uk,ua] Dphi[la] Dphi[lk],


     bssnham -> em4phi hInv[ui,uj] (R[li,lj] + Rphi[li,lj]) - A[li,lj] AInv[ui,uj] + 2/3 K^2,
     bssnmom[li] -> 6 hInv[uj,uk] A[lk,li] Dphi[lj]
                   + hInv[um,uk] (PDN[A[lm,li],lk]- gamma[ud,lm,lk] A[ld,li]
                                                 - gamma[ud,li,lk] A[lm,ld])
                   -2/3 PDN[K,li],
     bssnCcons[ui] -> Gam[ui]  - hInv[uj,uk] gamma[ui,lj,lk],
     bssnD -> Log[deth],
     bssnT -> hInv[ui,uj] A[li,lj],

     (* Add matter contributions here *)
     bssnham -> bssnham - 16 pi rho,
     bssnmom[li] -> bssnmom[li] - 8 pi Si[li]

  }
}; 

(**************************************************************************************)
(* Construct the thorn *)
(**************************************************************************************)

fdOrderParam = 
{
  Name -> "fd_order",
  Description -> "Finite differencing scheme order",
  Default -> "2nd",
  AllowedValues -> {"2nd", "2ndcentred", (* "4th", *) "full4th", "4thcentred", "6thcentred", "6th"},
  Steerable -> Recover
};

calculationsVacuum =
{
  initGFsCalc,
  initConstraintGFs,
  admToEvolvedCalc,
  (*evolvedToADMCalc, *)

  gammaEvolvedCalc["2nd", PDstandard2nd, PDonesided2nd],
(*  gammaEvolvedCalc["4th", PDstandard4th, PDonesided2nd],*)
  gammaEvolvedCalc["full4th", PDstandard4th, PDlopsided4th],
  gammaEvolvedCalc["4thcentred", PDstandard4th, PDstandard4th],
  gammaEvolvedCalc["6thcentred", PDstandard6th, PDstandard6th],
  gammaEvolvedCalc["6th", PDstandard6th, PDlopsided6th],
  gammaEvolvedBoundCalc[PDonesided2nd],

  evolveMetricCalc["2nd", PDstandard2nd, PDonesided2nd],
(*  evolveCalc["4th", PDstandard4th, PDonesided2nd],*)
  evolveMetricCalc["full4th", PDstandard4th, PDlopsided4th],
  evolveMetricCalc["2ndcentred", PDstandard2nd, PDstandard2nd],
  evolveMetricCalc["4thcentred", PDstandard4th, PDstandard4th],
  evolveMetricCalc["6thcentred", PDstandard6th, PDstandard6th],
  evolveMetricCalc["6th", PDstandard6th, PDlopsided6th],

  evolveNonMetricCalc["2nd", PDstandard2nd, PDonesided2nd],
(*  evolveCalc["4th", PDstandard4th, PDonesided2nd],*)
  evolveNonMetricCalc["full4th", PDstandard4th, PDlopsided4th],
  evolveNonMetricCalc["2ndcentred", PDstandard2nd, PDstandard2nd],
  evolveNonMetricCalc["4thcentred", PDstandard4th, PDstandard4th],
  evolveNonMetricCalc["6thcentred", PDstandard6th, PDstandard6th],
  evolveNonMetricCalc["6th", PDstandard6th, PDlopsided6th],

  psuProjectConstraintsCalc,
  psuBoundaryProjectConstraintsCalc,
  boundaryCalc,
  boundaryCalcFLRW,
  boundaryCalcFLRWConformal,
  boundaryCalcFLRWRadiative,
  boundaryCalcFLRWConformalRadiative,
  onePlusLogLapseCalc,
  onePlusLogLapseFull4thCalc,
  onePlusLogLapse2ndCentredCalc,
  onePlusLogLapse4thCentredCalc,
  onePlusLogLapse6thCentredCalc, 
  onePlusLogLapse6thCalc, 
  harmonicLapseCalc,
  exactLapseCalc,
  FLRWLapseCalc,  
  FLRWLapseCalc2,  
  FLRWLapseCalc3,  
  exactShiftCalc,
  vacuumconstraintsCalc,
  (*constraintsCalc["4th",PDstandard4th, PDlopsided4th],
  constraintsCalc["2nd",PDstandard2nd, PDonesided2nd],*)

  nasaShiftCalc["2nd", PDstandard2nd, PDonesided2nd],
(*  nasaShiftCalc["4th", PDstandard4th, PDonesided2nd],*)
  nasaShiftCalc["full4th", PDstandard4th, PDlopsided4th],
  nasaShiftCalc["2ndcentred", PDstandard2nd, PDstandard2nd],
  nasaShiftCalc["4thcentred", PDstandard4th, PDstandard4th],
  nasaShiftCalc["6thcentred", PDstandard6th, PDstandard6th],
  nasaShiftCalc["6th", PDstandard6th, PDlopsided6th],
  muellerShiftCalc["2nd", PDstandard2nd, PDonesided2nd]
};

calculationsMatter =
{
  evolvedToADMInitialCalc,
  evolvedToADMPoststepCalc,
  evolveNonMetricCalcMat,
  matterconstraintsCalc,
  nullifyGREvol
};

If[addMatter, calculations = Join[calculationsVacuum,calculationsMatter],calculations = Join[calculationsVacuum,{evolvedToADMCalc}]];
 
keywordParametersVacuum = 
{
  lapseParam, 
  shiftParam,
  fdOrderParam,
  boundaryParam,
  algconsParam,
  recalcConstraintsParam,
  updateADMParam
};

keywordParametersMatter =
{
  CoupleMatterParam
};
 
If[addMatter, keywordParameters = Join[keywordParametersVacuum,keywordParametersMatter], keywordParameters=keywordParametersVacuum];

extendedKeywordParameters = Join[ evolutionMethodParameters ];

CreateKrancThornTT[groups, ".", thornName,
  Calculations -> calculations,
  DeclaredGroups -> declaredGroupNames,
  PartialDerivatives -> derivatives,
  RealParameters -> realParameters,
  IntParameters -> intParameters,
  KeywordParameters -> keywordParameters,
  ExtendedKeywordParameters -> extendedKeywordParameters,
  UseLoopControl -> True,
  (*UseJacobian -> True,*)
  (*UseVectors -> True,*)
  If[addMatter, InheritedImplementations -> {"admbase","TmunuBase"}, InheritedImplementations -> {"admbase"} ]];

, {loopIndex, 1, 4}];
