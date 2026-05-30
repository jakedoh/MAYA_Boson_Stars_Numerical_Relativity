1;2000;0c(* ::Package:: *)

Get["KrancThorn`"];

SetEnhancedTimes[False];
SetSourceLanguage["C"];

(**************************************************************************************)
(* variable types *)
(**************************************************************************************)
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
    BDpi, BDphi, g3
}];


(* Set the attributes of all tensor-type quantities for which the
   defaults are not suitable. *)

SetTensorAttribute[phi, TensorWeight, 1/6];
SetTensorAttribute[chi, TensorWeight, -2/3];
SetTensorAttribute[A, TensorWeight, -2/3];
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
  g[la,lb], k[la,lb], Sij[la,lb], g3[la,lb]
}];

AssertSymmetricDecreasing[gamma[ua,lb,lc], lb, lc];

(* Determinants of the metrics in terms of their components
   (Mathematica symbolic expressions) *)
hDet = Det[MatrixOfComponents[h[la,lb]]];
gDet = Det[MatrixOfComponents[g[la,lb]]];

Do[

phiMethod = Which[loopIndex==1, True, loopIndex==2, False, loopIndex==3, False];
addMatter = If[loopIndex == 3, True, False];

fnPrefix = If[addMatter,"bssnchimatter",If[phiMethod, "bssn", "bssnchi"]];
thornName = If[addMatter,"Kranc2BSSNChiMatter",If[phiMethod, "Kranc2BSSN", "Kranc2BSSNChi"]];

(**************************************************************************************)
(* Groups *)
(**************************************************************************************)

(* Cactus group definitions *)

evolvedGroups = Map[CreateGroupFromTensor, {h[li,lj], A[li,lj], If[phiMethod, phi, chi], K, Gam[ui], alpha, beta[ui], betat[ui], BDphi, BDpi}];

admGroups = 
  {{"admbase::metric", {gxx,gxy,gxz,gyy,gyz,gzz}},
   {"admbase::curv", {kxx,kxy,kxz,kyy,kyz,kzz}},
   {"admbase::lapse", {alp}},
   {"admbase::shift", {betax,betay,betaz}}};

evaluatedGroups =
  {{"scalarconstraints", {bssnham, bssnT, bssnD}}, CreateGroupFromTensor[bssnmom[li]], CreateGroupFromTensor[bssnCcons[li]], CreateGroupFromTensor[g3[li,lj]]};

stressenergyGroups =
  {{"TmunuBase::stress_energy_scalar", {eTtt} },
   {"TmunuBase::stress_energy_vector", {eTtx,eTty,eTtz} },
   {"TmunuBase::stress_energy_tensor", {eTxx,eTxy,eTxz,eTyy,eTyz,eTzz} } };

declaredGroups = Join[evolvedGroups, evaluatedGroups];
declaredGroupNames = Map[First, declaredGroups];

If[addMatter, groups = Join[declaredGroups, admGroups, stressenergyGroups];, groups = Join[declaredGroups, admGroups];];

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
  rho, S, Si[la], Sij[la,lb], Tj[la], pi
};

(**************************************************************************************)
(* Parameters *)
(**************************************************************************************)

realParameters = 
{
  {Name -> m, Default -> 1},
  diss, 
  etaBeta, chiBeta,
  gammaDriverLapsePower, NASAAdvection, lapseAdvection,
  {Name -> betaDotAlphaFactor, Default -> 3/4},
  {Name -> gammaDriverLambda, Default -> 1},
  {Name -> newNASAAdvection, Default -> 0},
  {Name -> betatAdvection, Default -> 0},
  {Name -> chiEps, Default -> 0},
  {Name -> chiUniversalCutoff, Default -> -1},
  {Name -> delayGREvolutiontoTime, Default -> 0}
};

(**************************************************************************************)
(* Initialize all grid functions to a recognizable value. These should be overwritten *)
(* later. *)
(**************************************************************************************)

INITVALUE = "100";

initGFsCalc = 
{
  Name -> fnPrefix <> "_init_gfs",
  Schedule -> {"at BaseGrid"},
  Equations -> 
  {
    h[li,lj] -> INITVALUE,
    A[li,lj] -> INITVALUE,
    If[phiMethod, phi, chi] -> INITVALUE,
    K -> INITVALUE,
    Gam[ui] -> INITVALUE,
    alpha -> INITVALUE,
    beta[ui] -> INITVALUE,
    betat[ui] -> INITVALUE,
    BDphi -> INITVALUE,
    BDpi  -> INITVALUE,
    g3[li,lj] -> INITVALUE
  }
};

(**************************************************************************************)
(* ADMBase <--> BSSN variable translation *)
(**************************************************************************************)

admToEvolvedCalc = 
{
  Name -> fnPrefix <> "_adm_to_evolved",
  Schedule -> {"at INITIAL after ADMBase_PostInitial"},
  Shorthands -> shorthands,
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
     Gam[ui] -> 0, (* Really this wants to be populated via extrapolation later. Then
                     we can have two boundary points instead of three. *)

     
     g3[li,lj] -> g[li,lj],
    BDpi -> 0, 
     (*BDphi  -> 1 + Exp[-(Sqrt[x^2+y^2+z^2]-50)^2/5^2] / 10^6 / Sqrt[x^2+y^2+z^2]*)
    BDphi  -> 1 + Exp[-(Sqrt[x^2+y^2+z^2]-7)^2/25]/100

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
  {"at INITIAL after ADMBase_InitialData after " <> fnPrefix <> "_adm_to_evolved as " <> fnPrefix <> "_calc_evolved_gammas"},
  Shorthands -> shorthands,

  Equations ->
  {
     deth -> hDet, invdeth -> 1 / deth, 
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
  {"at INITIAL after ADMBase_InitialData after " <> fnPrefix <> "_adm_to_evolved"},
  Shorthands -> shorthands,

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
  Default -> "always",
  AllowedValues -> {"always", "never"}
};

evolvedToADMInitialCalc = 
{
  Name -> fnPrefix <> "_evolved_to_adm_initial",
  Schedule -> {"AT Initial as evolved_to_adm_initial after " <> fnPrefix <> "_adm_to_evolved before HydroBase_Prim2ConInitial"},
  Shorthands -> shorthands,
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
     betaz -> beta3,

     g3[la,lb] -> g[la,lb]
   }
};

evolvedToADMPoststepCalc = 
{
  Name -> fnPrefix <> "_evolved_to_adm",
  Schedule -> {"in MoL_PostStep as evolved_to_adm after " <> thornName <> "_ApplyBCs before HydroBase_PostStep"},
  Shorthands -> shorthands,
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
     betaz -> beta3,


    g3[li,lj] -> g[li,lj]
   }
};

evolvedToADMCalc = 
{
  Name -> fnPrefix <> "_evolved_to_adm",
  Schedule -> {"at CCTK_EVOL as evolved_to_adm after MoL_Evolution"},
  Shorthands -> shorthands,
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
     betaz -> beta3,

     g3[li,lj] -> g[li,lj]
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
  Shorthands -> shorthands,
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
  Shorthands -> shorthands,
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
         - 4/3 (Gam[ui] - hInv[ul,um] gamma[ui,ll,lm]) PD[beta[uj],lj]
         - 2 AInv[ui,uj] PD[alpha,lj]
         + 2 alpha * ( - 2/3 hInv[ui,uj] PD[K,lj]
		       + (gamma[ui,lk,ll] AInv[uk,ul] + 6 AInv[ui,uj] Dphi[lj]))  ),



     (* Brans-Dickie scalar field *)
    e4phi -> If[phiMethod, Exp[4 phi], 1/Max[chi,chiEps]],
    g[li,lj]    -> h[li,lj] e4phi,
    detg -> gDet,
    invdetg -> 1 / detg,
    gInv[ua,ub] -> invdetg gDet MatrixInverse[g[ua,ub]],

   

    dot[BDphi] -> alpha BDpi + beta[ui] PD[BDphi,li],
    dot[BDpi]  -> beta[ui] PD[BDpi,li] + alpha K BDpi
    + alpha gInv[ui,uj] PD[BDphi,li,lj]
    - 0.5 alpha gInv[ui,uj] gInv[ul,um] PD[g3[ll,lm],lj] PD[BDphi,li]
     + gInv[ui,uj] PD[alpha,lj] PD[BDphi,li] 

    (*  dot[BDphi] -> BDpi,
     dot[BDpi]  -> Euc[ua,ub] PD[BDphi,la,lb]*)
    
  }
};

CoupleMatterParam =
{
  Name -> "couple_matter",
  Default -> "yes",
  AllowedValues -> {"yes","no"}
};

nullifyGREvol =
{
  Name -> fnPrefix <> "_nullify_GREvol",
  Schedule -> {"in MoL_PostRHS after dissipation_add as bssn_nullify "},
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
  Shorthands -> shorthands,
  Equations -> 
  {

     deth -> hDet,
     invdeth -> 1 / deth,
     hInv[ua,ub] -> invdeth hDet MatrixInverse[h[ua,ub]],
     em4phi -> If[phiMethod, Exp[-4 phi], chi],
     e4phi -> If[phiMethod, Exp[4 phi], 1/Max[chi,chiEps]],

     (* auxiliary variables *)
     betaSqr -> e4phi h[li,lj] beta[ui] beta[uj],
     dBDphiSqr -> dot[BDphi]^2/alpha^2 - 2 BDpi dot[BDphi]/alpha
     + em4phi hInv[ui,uj] PD[BDphi,li] PD[BDphi,lj],
     dotdotBDphi -> alpha dot[BDpi] + dot[alpha] BDpi + dot[beta[ui]] PD[BDphi,li]
     + beta[ui] ( alpha PD[BDpi,li] + PD[alpha,li] BDpi
		  + PD[beta[uj],li] PD[BDphi,lj] 
		  + beta[uj] PD[BDphi,li,lj] ),

     (* Christoffel symbols in terms of 3+1 quantities *)
     (* First index is upper the others at lower \Gamma^a_{bc} *)
     (* Expressions from Appendix B, M. Alcubierre book *)
     k[la,lb] -> e4phi A[la,lb] + 1/3 g[la,lb] K,
     Gam000 -> ( dot[alpha] + beta[ui] PD[alpha,li] 
		 - beta[ui] beta[uj] k[li,lj] ) / alpha,
     Gami00[ui] -> alpha PD[alpha,lj] em4phi hInv[uj,ui]
     - 2 alpha beta[uj] k[lj,lm] em4phi hInv[um,ui]
     - beta[ui] ( dot[alpha] + beta[um] PD[alpha,lm] 
		  - beta[um] beta[un] k[lm,ln] ) / alpha + dot[beta[ui]] 
     + beta[uj] ( PD[beta[ui],lj] + beta[uk] *
		  ( gamma[ui,lj,lk] +
		    If[phiMethod,
		       2 ( Euc[ui,lk] PD[phi,lj] +
			   Euc[ui,lj] PD[phi,lk] -
			   hInv[ui,ul] PD[phi,ll] h[lj,lk] ),
		    - 1/(2 chi) ( Euc[ui,lk] PD[chi,lj] +
				  Euc[ui,lj] PD[chi,lk] -
				  hInv[ui,ul] PD[chi,ll] h[lj,lk] )] )),
     Gam00i[li] -> ( PD[alpha,li] - beta[um] k[li,lm] ) / alpha,
     Gamj0i[uj,li] -> -beta[uj] ( PD[alpha,li] - beta[un] k[li,ln] ) / alpha
     - alpha em4phi hInv[uj,un] k[ln,li] + PD[beta[ui],lj] + beta[uk] *
		  ( gamma[uj,li,lk] +
		    If[phiMethod,
		       2 ( Euc[uj,lk] PD[phi,li] +
			   Euc[uj,li] PD[phi,lk] -
			   hInv[uj,ul] PD[phi,ll] h[li,lk] ),
		    - 1/(2 chi) ( Euc[uj,lk] PD[chi,li] +
				  Euc[uj,li] PD[chi,lk] -
				  hInv[uj,ul] PD[chi,ll] h[li,lk] )] ),
		  
		    
					    
			   



     (* Components due to the scalar Brans-Dicke field *)
     BDTtt -> BDomega/BDphi^2 ( dot[BDphi]^2 - 1/2 (-alpha^2 + betaSqr) dBDphiSqr)
     + 1/BDphi (dotdotBDphi - Gam000 dot[BDphi] - Gami00[ui] PD[BDphi,li] ),

     BDTt[li] -> BDomega/BDphi^2 (dot[BDphi] PD[BDphi,li] - 
				  1/2 e4phi h[li,lj] beta[uj] dBDphiSqr ) +
		  1/BDphi ( PD[alpha,li] BDpi + alpha PD[BDpi,li] +
			    PD[beta[uj],li] PD[BDphi,lj] +
			    beta[uj] PD[BDphi,lj,li] -
			    Gam00i[li] dot[BDphi] -
			    Gamj0i[uj,li] PD[BDphi,lj] ),


     BDTxx -> ,
     BDTxy -> ,
     BDTxz -> ,
     BDTyy -> ,
     BDTyz -> ,
     BDTzz -> ,
     

     (* Quantities for matter in the RHS *)
     (* This comes from whatever non gravity matter fields are being evolved *)
     Sij11 -> eTxx, Sij21 -> eTxy, Sij31 -> eTxz, 
     Sij22 -> eTyy, Sij32 -> eTyz, Sij33 -> eTzz,
     Tj1 -> eTtx, Tj2 -> eTty, Tj3 -> eTtz,

     (* rescaling of T_{ab} by Brans-Dicke field *)
     Sij[la,lb] -> Sij[la,lb] / BDphi,
     Tj[la] -> Tj[la] / BDphi,

    
     (* ij components of scalar field T_{ab} *)
     BDTij[li,lj] -> BDomega/BDphi^2 ( PD[BDphi,li] PD[BDphi,lj] -1/2 e4phi h[li,lj] PDBDphi2 )
     + 1/BDphi ( 

     Si[la] -> (1/alpha) ( - Tj[la] + beta[ub] Sij[la,lb] ),
     S -> em4phi hInv[ua,ub] Sij[la,lb],
     rho -> (1/alpha)^2 ( eTtt/BDphi - 2 beta[ua] Tj[la] + beta[ua] beta[ub] Sij[la,lb] ),
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
  Default -> "constant",
  AllowedValues -> {"1 + log", "1 + log full4th", "1 + log 4thcentred", "constant",
    "harmonic", "exact", "1 + log 6thcentred", "1 + log 6th"}
};

constantLapseCalc = 
{
  Name -> fnPrefix <> "_constant_lapse",
  Schedule -> {"in MoL_CalcRHS"},
  ConditionalOnKeyword -> {"lapse_condition", "constant"},
  Shorthands -> shorthands,
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
  Schedule -> {"in MoL_PostStep after Exact__gauge"},
  ConditionalOnKeyword -> {"lapse_condition", "exact"},
  Shorthands -> shorthands,
  Where -> Everywhere,
  Equations -> 
  {
    alpha -> alp
  }
};

harmonicLapseCalc = 
{
  Name -> fnPrefix <> "_harmonic_lapse",
  Schedule -> {"in MoL_CalcRHS"},
  Where -> Interior,
  ConditionalOnKeyword -> {"lapse_condition", "harmonic"},
  Shorthands -> shorthands,
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
  Schedule -> {"in MoL_CalcRHS"},
  Where -> Interior,
  ConditionalOnKeyword -> {"lapse_condition", "1 + log"},
  Shorthands -> shorthands,
  Equations -> 
  {
    (* Shift advection *)
    dir[ui] -> Sign[beta[ui]],

    (* WARNING: This will be second order if lapseAdvection != 0 *)
    dot[alpha] -> -2 alpha K + lapseAdvection beta[ua] PDonesided2nd[alpha,la] 
  }
};

onePlusLogLapseFull4thCalc = 
{
  Name -> fnPrefix <> "_one_plus_log_lapse_full4th",
  Schedule -> {"in MoL_CalcRHS"},
  Where -> Interior,
  ConditionalOnKeyword -> {"lapse_condition", "1 + log full4th"},
  Shorthands -> shorthands,
  Equations -> 
  {
    (* Shift advection *)
    dir[ui] -> Sign[beta[ui]],

    dot[alpha] -> -2 alpha K + lapseAdvection beta[ua] PDlopsided4th[alpha,la] 
  }
};

onePlusLogLapse6thCalc = 
{
  Name -> fnPrefix <> "_one_plus_log_lapse_6th",
  Schedule -> {"in MoL_CalcRHS"},
  Where -> Interior,
  ConditionalOnKeyword -> {"lapse_condition", "1 + log 6th"},
  Shorthands -> shorthands,
  Equations -> 
  {
    (* Shift advection *)
    dir[ui] -> Sign[beta[ui]],

    dot[alpha] -> -2 alpha K + lapseAdvection beta[ua] PDlopsided6th[alpha,la] 
  }
};

onePlusLogLapse4thCentredCalc = 
{
  Name -> fnPrefix <> "_one_plus_log_lapse_4thcentred",
  Schedule -> {"in MoL_CalcRHS"},
  Where -> Interior,
  ConditionalOnKeyword -> {"lapse_condition", "1 + log 4thcentred"},
  Shorthands -> shorthands,
  Equations -> 
  {
    (* Shift advection *)
    dir[ui] -> Sign[beta[ui]],

    dot[alpha] -> -2 alpha K + lapseAdvection beta[ua] PDstandard4th[alpha,la] 
  }
};

onePlusLogLapse6thCentredCalc = 
{
  Name -> fnPrefix <> "_one_plus_log_lapse_6thcentred",
  Schedule -> {"in MoL_CalcRHS"},
  Where -> Interior,
  ConditionalOnKeyword -> {"lapse_condition", "1 + log 6thcentred"},
  Shorthands -> shorthands,
  Equations -> 
  {
    (* Shift advection *)
    dir[ui] -> Sign[beta[ui]],

    dot[alpha] -> -2 alpha K + lapseAdvection beta[ua] PDstandard6th[alpha,la] 
  }
};

(**************************************************************************************)
(* Shift conditions *)
(**************************************************************************************)

shiftParam = 
{
  Name -> "shift_condition",
  Default -> "constant",
  AllowedValues -> {"constant", "NASA", "NASA4th", "NASAfull4th", "NASA6th", "NASA4thcentred", "exact"}
};

constantShiftCalc =
{
  Name -> fnPrefix <> "_constant_shift",
  Schedule -> {"in MoL_CalcRHS"},
  Where -> Interior,
  Shorthands -> shorthands,
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
  Schedule -> {"in MoL_PostStep after Exact__gauge"},
  Where -> Everywhere,
  Shorthands -> shorthands,
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
  Schedule -> {"in MoL_CalcRHS after bssn_evolve after bssn_evolve_matter"},
  Where -> Interior,
  Shorthands -> shorthands,
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

(* Need to modify Kranc so that we can set multiple conditional
parameters on a calculation. *)

(**************************************************************************************)
(* Boundary conditions *)
(**************************************************************************************)

boundaryParam = 
{
  Name -> "boundary_condition",
  Default -> "radiative",
  AllowedValues -> {"none", "radiative"}
};

boundaryCalc =
{
  Name -> fnPrefix <> "_boundary",
  Schedule -> {"in MoL_RHSBoundaries"},
  ConditionalOnKeyword -> {"boundary_condition", "radiative"},
  Where -> Boundary,
  Shorthands -> shorthands,
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
    bssnmom1 -> 0, bssnmom2 -> 0, bssnmom3 -> 0,

    (* Brans-Dicke boundary conditions *)
    BDpi  -> 0,
    BDphi -> 1


  }
};

(**************************************************************************************)
(* Project the BSSN algebraic constraints *)
(**************************************************************************************)

algconsParam = 
{
  Name -> "algebraic_constraint_projection",
  Default -> "psu",
  AllowedValues -> {"none", "psu", "psuboundary"}
};

psuProjectConstraintsCalc = 
{
  Name -> fnPrefix <> "_project_constraints_psu",
  Schedule -> {"in MoL_PostStep"},
  ConditionalOnKeyword -> {"algebraic_constraint_projection", "psu"},
  Shorthands -> shorthands,
  Equations ->
  {
    invdeth -> 1 / hDet,
    hInv[ua,ub] -> invdeth hDet MatrixInverse[h[ua,ub]],
    trA ->  A[lm,ln] hInv[um,un],

    If[phiMethod, invdeth -> 1 / hDet , chi -> Max[chi, chiUniversalCutoff]],

    A[la,lb] -> A[la,lb] - 1/3 h[la,lb] trA,


    (* needed for Brans-Dicke *)
     e4phi -> If[phiMethod, Exp[4 phi], 1/Max[chi,chiEps]],
     g[la,lb] -> h[la,lb] e4phi,
     g3[li,lj] -> g[li,lj]

  }
};

psuBoundaryProjectConstraintsCalc = 
{
  Name -> fnPrefix <> "_project_constraints_psu_boundary",
  Schedule -> {"in MoL_PostStep"},
  Where -> Boundary,
  ConditionalOnKeyword -> {"algebraic_constraint_projection", "psuboundary"},
  Shorthands -> shorthands,
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
  Default -> "no",
  AllowedValues -> {"4th", "no"}
};

constraintsCalc = 
{
  Name -> fnPrefix <> "_calc_constraints",
  Schedule -> {"at CCTK_ANALYSIS as bssn_calc_constraints"},
  Where -> Interior,
  ConditionalOnKeyword -> {"recalculate_constraints", "4th"},
  Shorthands -> shorthands,
  Equations ->
  {
        deth -> hDet,
        invdeth -> 1 / deth,
        hInv[ua,ub] -> invdeth hDet MatrixInverse[h[ua,ub]],
        gamma[ua, lb, lc] -> 1/2 hInv[ua,ud] (PDstandard4th[h[lb,ld], lc] + PDstandard4th[h[lc,ld], lb] - PDstandard4th[h[lb,lc],ld]),
        AInv[ua,ub] -> hInv[ua,ui] hInv[ub,uj] A[li,lj],
        em4phi -> If[phiMethod, Exp[-4 phi], chi],
        chiCutoff -> If[phiMethod, 0, Max[chi,chiEps]],


        DDphi[lj,li] -> If[phiMethod, PDstandard4th[phi,lj,li] - gamma[uk,li,lj] PDstandard4th[phi,lk],
                                      -1/(4 chiCutoff) ( PDstandard4th[chi,lj,li] - gamma[uk,li,lj] PDstandard4th[chi,lk] 
                                                   - chiCutoff^-1 PDstandard4th[chi,li] PDstandard4th[chi,lj] )],
        Dphi[li] -> If[phiMethod, PDstandard4th[phi,li], -1/(4 chiCutoff) PDstandard4th[chi,li]],

        R[li,lj] -> -1/2 hInv[uk,ul] PDstandard4th[h[li,lj],ll,lk] + TTSymmetrize[h[lk,li] PDstandard4th[Gam[uk],lj], li,lj] 
                    + TTSymmetrize[h[li,lc] gamma[uc,lj,lk] Gam[uk],li,lj]
                    + hInv[ul,us](2 TTSymmetrize[gamma[uk,ll,li] h[lj,la] gamma[ua,lk,ls],li,lj]
                                  + gamma[uk,li,ls] h[la,lk] gamma[ua,ll,lj]),

        Rphi[li,lj] -> -2 DDphi[lj,li] - 2 h[li,lj] hInv[uk,ul] DDphi[lk,ll] + 4 Dphi[li] Dphi[lj] 
                       - 4 h[li,lj] hInv[uk,ua] Dphi[la] Dphi[lk],


        bssnham -> em4phi hInv[ui,uj] (R[li,lj] + Rphi[li,lj]) - A[li,lj] AInv[ui,uj] + 2/3 K^2,
        bssnmom[li] -> 6 hInv[uj,uk] A[lk,li] Dphi[lj]
                      + hInv[um,uk] (PDstandard4th[A[lm,li],lk]- gamma[ud,lm,lk] A[ld,li] 
                                                    - gamma[ud,li,lk] A[lm,ld]) 
                      -2/3 PDstandard4th[K,li],
	bssnCcons[ui] -> Gam[ui]  - hInv[uj,uk] gamma[ui,lj,lk],
        bssnD -> Log[deth],
        bssnT -> hInv[ui,uj] A[li,lj]
  }
};

MatterConstraintsParam =
{
  Name -> "matter_constraints",
  Default -> "no",
  AllowedValues -> {"yes","no"}
};

matconstraintsCalc = 
{
  Name -> fnPrefix <> "_calc_matter_constraints",
  Schedule -> {"at CCTK_ANALYSIS after bssn_calc_constraints"},
  Where -> Interior,
  ConditionalOnKeyword -> {"matter_constraints", "yes"},
  Shorthands -> shorthands,
  Equations ->
  {

     (* Quantities for matter in the RHS *)
     Sij11 -> eTxx, Sij21 -> eTxy, Sij31 -> eTxz, 
     Sij22 -> eTyy, Sij32 -> eTyz, Sij33 -> eTzz,
     Tj1 -> eTtx, Tj2 -> eTty, Tj3 -> eTtz,
     Si[la] -> (1/alpha) ( - Tj[la] + beta[ub] Sij[la,lb] ),
     rho -> (1/alpha)^2 ( eTtt - 2 beta[ua] Tj[la] + beta[ua] beta[ub] Sij[la,lb] ),
     pi -> 3.14159265358979323846,

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
  Default -> "2nd",
  AllowedValues -> {"2nd", (* "4th", *) "full4th", "4thcentred", "6thcentred", "6th"}
};

calculationsVacuum =
{
  initGFsCalc,
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
  evolveMetricCalc["4thcentred", PDstandard4th, PDstandard4th],
  evolveMetricCalc["6thcentred", PDstandard6th, PDstandard6th],
  evolveMetricCalc["6th", PDstandard6th, PDlopsided6th],

  evolveNonMetricCalc["2nd", PDstandard2nd, PDonesided2nd],
(*  evolveCalc["4th", PDstandard4th, PDonesided2nd],*)
  evolveNonMetricCalc["full4th", PDstandard4th, PDlopsided4th],
  evolveNonMetricCalc["4thcentred", PDstandard4th, PDstandard4th],
  evolveNonMetricCalc["6thcentred", PDstandard6th, PDstandard6th],
  evolveNonMetricCalc["6th", PDstandard6th, PDlopsided6th],

  psuProjectConstraintsCalc,
  psuBoundaryProjectConstraintsCalc,
  boundaryCalc,
  onePlusLogLapseCalc,
  onePlusLogLapseFull4thCalc,
  onePlusLogLapse4thCentredCalc,
  onePlusLogLapse6thCentredCalc, 
  onePlusLogLapse6thCalc, 
  harmonicLapseCalc,
  exactLapseCalc,
  exactShiftCalc,
  constraintsCalc,
  nullifyGREvol,

  nasaShiftCalc["2nd", PDstandard2nd, PDonesided2nd],
(*  nasaShiftCalc["4th", PDstandard4th, PDonesided2nd],*)
  nasaShiftCalc["full4th", PDstandard4th, PDlopsided4th],
  nasaShiftCalc["4thcentred", PDstandard4th, PDstandard4th],
  nasaShiftCalc["6thcentred", PDstandard6th, PDstandard6th],
  nasaShiftCalc["6th", PDstandard6th, PDlopsided6th]
};

If[addMatter, calculations = Join[calculationsVacuum,{evolvedToADMInitialCalc,evolvedToADMPoststepCalc,evolveNonMetricCalcMat,matconstraintsCalc}],calculations = Join[calculationsVacuum,{evolvedToADMCalc}]];

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

If[addMatter, keywordParameters = Join[keywordParametersVacuum,{CoupleMatterParam,MatterConstraintsParam}], keywordParameters=keywordParametersVacuum];

CreateKrancThornTT[groups, ".", thornName,
  Calculations -> calculations,
  DeclaredGroups -> declaredGroupNames,
  PartialDerivatives -> derivatives,
  RealParameters -> realParameters,
  KeywordParameters -> keywordParameters,
  If[addMatter, InheritedImplementations -> {"admbase","TmunuBase"}, InheritedImplementations -> {"admbase"} ]];

, {loopIndex, 1, 3}];
