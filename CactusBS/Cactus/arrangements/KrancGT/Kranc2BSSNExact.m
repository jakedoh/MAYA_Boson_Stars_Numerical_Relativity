
Get["KrancThorn`"];

SetEnhancedTimes[False];
SetSourceLanguage["C"];

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

  PDplus[i_] -> DPlus[i],
  PDminus[i_] -> DMinus[i],
  PDplus[i_,j_] -> DPlus[i] DPlus[j],

  PDonesided2nd[1] -> dir[1] (-shift[1]^(2 dir[1]) + 4 shift[1]^dir[1] - 3 )/(2 spacing[1]),
  PDonesided2nd[2] -> dir[2] (-shift[2]^(2 dir[2]) + 4 shift[2]^dir[2] - 3 )/(2 spacing[2]),
  PDonesided2nd[3] -> dir[3] (-shift[3]^(2 dir[3]) + 4 shift[3]^dir[3] - 3 )/(2 spacing[3]),

  PDlopsided4th[1] -> dir[1] (-3 shift[1]^(-1 dir[1]) - 10 + 18 shift[1]^(1 dir[1]) - 6 shift[1]^(2 dir[1]) + shift[1]^(3 dir[1])) / (12 spacing[1]),
  PDlopsided4th[2] -> dir[2] (-3 shift[2]^(-1 dir[2]) - 10 + 18 shift[2]^(1 dir[2]) - 6 shift[2]^(2 dir[2]) + shift[2]^(3 dir[2])) / (12 spacing[2]),
  PDlopsided4th[3] -> dir[3] (-3 shift[3]^(-1 dir[3]) - 10 + 18 shift[3]^(1 dir[3]) - 6 shift[3]^(2 dir[3]) + shift[3]^(3 dir[3])) / (12 spacing[3])
};

(**************************************************************************************)
(* Tensors *)
(**************************************************************************************)

(* Register all the tensors that will be used with TensorTools *)
Map[DefineTensor, 
{
  h, hInv, phi, chi, A, K, alpha, Gam, beta, betat, gamma, 
  g, k, AInv, gInv, 
  n, dir, 
  hExact, phiExact, chiExact, AExact, KExact, alphaExact, GamExact, betaExact, betatExact,
  hError, phiError, chiError, AError, KError, alphaError, GamError, betaError, betatError
}];


(* Set the attributes of all tensor-type quantities for which the
   defaults are not suitable. *)

SetTensorAttribute[phiExact, TensorWeight, 1/6];
SetTensorAttribute[chiExact, TensorWeight, -2/3];
SetTensorAttribute[AExact, TensorWeight, -2/3];
SetTensorAttribute[GamExact, TensorWeight, 2/3];
SetTensorAttribute[GamExact, TensorSpecial, "gamma"];

(* Register the TensorTools symmetries (this is very simplistic) *)
Map[AssertSymmetricIncreasing,
{
  hInv[ua,ub], AInv[ua,ub]
}];

Map[AssertSymmetricDecreasing, 
{
  h[la,lb], hInv[ua,ub], A[la,lb], g[la,lb], k[la,lb], hExact[la,lb], AExact[la,lb]
}];

AssertSymmetricDecreasing[gamma[ua,lb,lc], lb, lc];

(* Determinants of the metrics in terms of their components
   (Mathematica symbolic expressions) *)
hDet = Det[MatrixOfComponents[h[la,lb]]];
hExactDet = Det[MatrixOfComponents[hExact[la,lb]]];
gDet = Det[MatrixOfComponents[g[la,lb]]];

Do[

phiMethod = If[loopIndex == 1, True, False];

fnPrefix = If[phiMethod, "bssnexact", "bssnchiexact"];

(**************************************************************************************)
(* Groups *)
(**************************************************************************************)

(* Cactus group definitions *)

evolvedGroups = Map[CreateGroupFromTensor, {h[li,lj], A[li,lj], If[phiMethod, phi, chi], K, Gam[ui], alpha, beta[ui], betat[ui]}];

evolvedGroups = qualifyGroups[evolvedGroups, Map[groupName, evolvedGroups], "Kranc2BSSN", "Kranc2BSSNExact"];

exactGroups = Map[CreateGroupFromTensor, {hExact[li,lj], AExact[li,lj], If[phiMethod, phiExact, chiExact], KExact, GamExact[ui], alphaExact, betaExact[ui], betatExact[ui]}];

admGroups = 
  {{"admbase::metric", {gxx,gxy,gxz,gyy,gyz,gzz}},
   {"admbase::curv", {kxx,kxy,kxz,kyy,kyz,kzz}},
   {"admbase::lapse", {alp}},
   {"admbase::shift", {betax,betay,betaz}}};

declaredGroups = Join[exactGroups];
declaredGroupNames = Map[First, declaredGroups];

groups = Join[declaredGroups, admGroups, evolvedGroups];

(**************************************************************************************)
(* Shorthands *)
(**************************************************************************************)

shorthands = 
{
  deth, invdeth, detg, invdetg, hInv[ua,ub], gamma[ua, lb, lc],
  em4phi, e4phi, detgmthirdroot,
  third, twelfth, g[li,lj], k[li,lj], AInv[ua,ub], gInv[ua,ub], half, 
  n[ui], dir[ui]
};

(**************************************************************************************)
(* Parameters *)
(**************************************************************************************)

realParameters = 
{
};

(**************************************************************************************)
(* Initialize all grid functions to a recognizable value. These should be overwritten *)
(* later. *)
(**************************************************************************************)

INITVALUE = "100";

initGFsCalc = 
{
  Name -> fnPrefix <> "_init_gfs",
  Schedule -> {"in ADMBase_InitialData"},
  Equations -> 
  {
    hExact[li,lj] -> INITVALUE,
    AExact[li,lj] -> INITVALUE,
    If[phiMethod,  phiExact -> INITVALUE,
                   chiExact -> INITVALUE],
    KExact -> INITVALUE,
    GamExact[ui] -> INITVALUE,
    alphaExact -> INITVALUE,
    betaExact[ui] -> INITVALUE,
    betatExact[ui] -> INITVALUE
  }
};

(**************************************************************************************)
(* ADMBase <--> BSSN variable translation *)
(**************************************************************************************)

admToExactCalc = 
{
  Name -> fnPrefix <> "_adm_to_exact",
  Schedule -> {"at ANALYSIS"},
  Where -> Everywhere,
  Shorthands -> shorthands,
  Equations ->
    {
     g11 -> gxx, g21 -> gxy, g22 -> gyy, g31 -> gxz, g32 -> gyz, g33 -> gzz,
     third -> 1/3, twelfth -> 1/12,
     detg -> gDet,
     detgmthirdroot -> detg^(-third),
     invdetg -> 1 / detg,
     gInv[ua,ub] -> invdetg gDet MatrixInverse[g[ua,ub]],
     k11 -> kxx, k21 -> kxy, k22 -> kyy, k31 -> kxz, k32 -> kyz, k33 -> kzz,
     hExact[la,lb] -> g[la,lb] * detgmthirdroot,
     If[phiMethod, phiExact -> twelfth Log[detg], chiExact -> detg^(-1/3)],
     KExact -> k[la,lb] gInv[ua,ub],
     If[phiMethod, em4phi -> Exp[-4 phiExact], em4phi -> chiExact],
     AExact[la,lb] -> em4phi (k[la,lb] - third K g[la,lb]),
     alphaExact -> alp,
     betaExact1 -> betax, betaExact2 -> betay, betaExact3 -> betaz,
     betatExact1 -> 0, betatExact2 -> 0, betatExact3 -> 0,
     GamExact[ui] -> 0 
   }
};

(* This is in a separate calculation from the other variables as it
   needs to take finite differences and hence the loop is only in the
   interior. We calculate the Gams based on h rather than g for
   convenience. *)

(*gammaExactCalc[fdOrder_, PD_, PDadvect_] := 
{
  Name -> "bssn_calc_exact_gammas_" <> fdOrder,
  ConditionalOnKeyword -> {"fd_order", fdOrder},
  Where -> Interior,
  Schedule -> 
  {"at ANALYSIS after bssn_adm_to_exact as bssn_calc_exact_gammas"},
  Shorthands -> shorthands,

  Equations ->
  {
     deth -> hExactDet, invdeth -> 1 / deth, half -> 1/2,
     hInv[ua,ub] -> invdeth hExactDet MatrixInverse[hExact[ua,ub]],
     gamma[ua,lb,lc] -> half hInv[ua,ud] (PD[hExact[lb,ld], lc] + PD[hExact[lc,ld], lb] - PD[hExact[lb,lc],ld]),

     GamExact[ui] -> hInv[uj,uk] gamma[ui,lj,lk]
  }
}

(* On the boundary we cannot calculate the Gams using a centered
   difference, so we use a one-sided difference pointing into the
   grid.  This is always used with a 2nd order accurate operator as a
   fourth order accurate version requires too large a stencil.*)

gammaExactBoundCalc[PD_] := 
{
  Name -> "bssn_calc_exact_bound_gammas",
  Where -> Boundary,
  Schedule -> 
  {"at ANALYSIS after bssn_calc_exact_gammas"},
  Shorthands -> shorthands,

  Equations ->
  {
    n1 -> -x/r, 
    n2 -> -y/r,
    n3 -> -z/r,
    dir[ui] -> Sign[n[ui]],

    deth -> hExactDet, invdeth -> 1 / deth, half -> 1/2,
    hInv[ua,ub] -> invdeth hDet MatrixInverse[hExact[ua,ub]],
    gamma[ua,lb,lc] -> half hInv[ua,ud] (PD[hExact[lb,ld], lc] + PD[hExact[lc,ld], lb] - PD[hExact[lb,lc],ld]),

    GamExact[ui] -> hInv[uj,uk] gamma[ui,lj,lk]
  }
}
*)

(**************************************************************************************)
(* Construct the thorn *)
(**************************************************************************************)

fdOrderParam = 
{
  Name -> "fd_order",
  Default -> "2nd",
  AllowedValues -> {"2nd", "4th", "full4th"}
};

calculations = 
{
  initGFsCalc,
  admToExactCalc (*,

  gammaExactCalc["2nd", PDstandard2nd, PDonesided2nd],
  gammaExactCalc["4th", PDstandard4th, PDonesided2nd],
  gammaExactCalc["full4th", PDstandard4th, PDlopsided4th],
  gammaExactBoundCalc[PDonesided2nd] *)
};

(*keywordParameters =
{
  fdOrderParam
};*)

CreateKrancThornTT[groups, ".", If[phiMethod, "Kranc2BSSNExact", "Kranc2BSSNChiExact"],
  Calculations -> calculations,
  DeclaredGroups -> declaredGroupNames,
  PartialDerivatives -> derivatives,
  RealParameters -> realParameters,
(*  KeywordParameters -> keywordParameters,*)
  InheritedImplementations -> {"admbase", If[phiMethod, "Kranc2BSSN", "Kranc2BSSNChi"]}];

, {loopIndex, 1, 2}];
