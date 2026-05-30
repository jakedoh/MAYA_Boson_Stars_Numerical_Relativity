
Get["KrancThorn`"];

(*SetDebugLevel[InfoFull];*)
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

  PDonesided2nd[1] -> dir[1] (-shift[1]^(2 dir[1]) + 4 shift[1]^dir[1] - 3 )/(2 spacing[1]),
  PDonesided2nd[2] -> dir[2] (-shift[2]^(2 dir[2]) + 4 shift[2]^dir[2] - 3 )/(2 spacing[2]),
  PDonesided2nd[3] -> dir[3] (-shift[3]^(2 dir[3]) + 4 shift[3]^dir[3] - 3 )/(2 spacing[3]),

  PDonesided4th[1] -> dir[1] (- 3 shift[1]^(4 dir[1]) + 16 shift[1]^(3 dir[1]) - 36 shift[1]^(2 dir[1]) + 48 shift[1]^dir[1] - 25 ) / (12 spacing[1]),
  PDonesided4th[2] -> dir[2] (- 3 shift[2]^(4 dir[2]) + 16 shift[2]^(3 dir[2]) - 36 shift[2]^(2 dir[2]) + 48 shift[2]^dir[2] - 25 ) / (12 spacing[2]),
  PDonesided4th[3] -> dir[3] (- 3 shift[3]^(4 dir[3]) + 16 shift[3]^(3 dir[3]) - 36 shift[3]^(2 dir[3]) + 48 shift[3]^dir[3] - 25 ) / (12 spacing[3]),

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
  beta, n, dir, swSigma, swEta, DDsigma, energy,
  Sij, Si, g, h, swPhi, chi, k, K, em4swPhi,
  gInv, hInv, deth, detg, invdetg, alpha, Ddetg
}];

(* Register the TensorTools symmetries (this is very simplistic) *)

Map[AssertSymmetricDecreasing, 
{
  Sij[la,lb], g[la,lb], k[la,lb], gInv[ua,ub], h[la,lb], hInv[ua,ub]
}];


(* Determinants of the metrics in terms of their components
   (Mathematica symbolic expressions) *)
gDet = Det[MatrixOfComponents[g[la,lb]]];
hDet = Det[MatrixOfComponents[h[la,lb]]];

(**************************************************************************************)
(* Groups *)
(**************************************************************************************)

(* Cactus group definitions *)

wave = Map[CreateGroupFromTensor, {swSigma, swEta, energy, detg, g[li,lj], beta[ui], alpha, K, gInv[ui,uj], h[li,lj], deth, hInv[li,lj], swPhi}];

fnPrefix = "ScalarField";

admGroups =  
  {{"admbase::metric", {gxx,gxy,gxz,gyy,gyz,gzz}},
   {"admbase::curv", {kxx,kxy,kxz,kyy,kyz,kzz}},
   {"admbase::lapse", {alp}},
   {"admbase::shift", {betax,betay,betaz}}};

matterGroups =
  {{"MatterBase::S_group", {S}},
   {"MatterBase::Si_group", {Si1,Si2,Si3}},
   {"MatterBase::Sij_group", {Sij11,Sij21,Sij31,Sij22,Sij32,Sij33}},
   {"MatterBase::rho_group", {rho}}};

declaredGroups = Join[wave];
declaredGroupNames = Map[First, declaredGroups];

groups = Join[declaredGroups, admGroups, matterGroups];

(**************************************************************************************)
(* Shorthands *)
(**************************************************************************************)

shorthands = { 
  invdetg, dir[ui], n[ui], DDsigma[li,lj], half, ADVsigma, ADVeta, rSqr, radi, k[la,lb], Ddetg[li], em4swPhi
};

(**************************************************************************************)
(* Parameters *)
(**************************************************************************************)

realParameters = 
{
  Amp, lambda
};

(**************************************************************************************)
(* Initialize all grid functions to a recognizable value. These should be overwritten *)
(* later. *)
(**************************************************************************************)

INITVALUE = "666";

initGFsCalc = {
  Name -> fnPrefix <> "_init_gfs",
  Schedule -> {"at INITIAL"},
  Equations -> {
    swSigma -> INITVALUE,
    swEta -> INITVALUE
  } 
};

(**************************************************************************************)
(* ADMBase <--> BSSN variable translation *)
(**************************************************************************************)

admToADMCalc =
{
  Name -> fnPrefix <> "_adm_to_adm",
  Schedule -> {"at INITIAL", "in MoL_PostStep as wave_ADMupdate"},
  Shorthands -> shorthands,
  Equations ->
  {
    g11 -> gxx, g21 -> gxy, g22 -> gyy, g31 -> gxz, g32 -> gyz, g33 -> gzz,
    k11 -> kxx, k21 -> kxy, k22 -> kyy, k31 -> kxz, k32 -> kyz, k33 -> kzz,
    beta1 -> betax, beta2 -> betay, beta3 -> betaz,
    alpha -> alp,
    K -> k[la,lb] gInv[ua,ub]
  }
};

admToBSSNCalc =
{
  Name -> fnPrefix <> "_adm_to_bssn",
  Schedule -> {"at INITIAL", "in MoL_PostStep after wave_ADMupdate as wave_BSSNupdate"},
  Shorthands -> shorthands,
  Equations ->
  {
    detg -> gDet,
    swPhi -> Log[detg]/12,
    h[la,lb] -> g[la,lb] * detg^(-1/3)
  }
};

MatterCoupleParam =
{
  Name -> "couple_matter",
  Default -> "never",
  AllowedValues -> {"always", "never"}
};

setMatterTermsADM[PD_, PDadvect_] :=
{
  Name -> fnPrefix <> "_set_matter_terms_adm",
  Schedule -> {"in MoL_PostStep after(wave_ADMupdate MatterBase_Init) as wave_updateMatter"},
  ConditionalOnKeyword -> {"decomp_type", "ADM"},
  Shorthands-> shorthands,
  Where -> Interior,
  Equations ->
  { 
    detg -> gDet,
    invdetg -> 1 / detg,
    gInv[ua,ub] -> invdetg gDet MatrixInverse[g[ua,ub]],
    ADVsigma -> beta[ui] PD[swSigma,li],
    
    rho -> rho + ( swEta^2 + gInv[ui,uj] PD[swSigma,li] PD[swSigma,lj] )/2,
    Si[li] -> Si[li] + ( -swEta ) PD[swSigma,li],
    Sij[li,lj] -> Sij[li,lj] + PD[swSigma,li] PD[swSigma,lj] - g[li,lj] ( gInv[uk,um] PD[swSigma,lk] PD[swSigma,lm] - swEta^2 )/2,
    S -> gInv[ui,uj] Sij[li,lj] 
  }
 };


setMatterTermsBSSN[PD_, PDadvect_] :=
{
  Name -> fnPrefix <> "_set_matter_terms_bssn",
  Schedule -> {"in MoL_PostStep after(wave_BSSNupdate MatterBase_Init) as wave_updateMatter"},
  ConditionalOnKeyword -> {"decomp_type", "BSSN"},
  Shorthands-> shorthands,
  Where -> Interior,
  Equations ->
  { 
    detg -> gDet,
    invdetg -> 1 / detg,
    gInv[ua,ub] -> invdetg gDet MatrixInverse[g[ua,ub]],
    ADVsigma -> beta[ui] PD[swSigma,li],
    em4swPhi -> Exp[-4 swPhi],

    (* the rho that I'm adding to should be 0 or set by another matter thorn
       I assume that any other matter thorn will appropriately set rho *)
    rho -> rho +  ( swEta^2 + gInv[ui,uj] PD[swSigma,li] PD[swSigma,lj] )/2,
    Si[li] -> Si[li] + ( -swEta ) PD[swSigma,li],
    Sij[li,lj] -> Sij[li,lj] + PD[swSigma,li] PD[swSigma,lj] - g[li,lj] ( gInv[uk,um] PD[swSigma,lk] PD[swSigma,lm] - swEta^2 )/2,
    S -> gInv[ui,uj] Sij[li,lj]
  }
 };

killMatterTerms =
{
  Name -> fnPrefix <> "kill_matter_terms",
  Schedule -> {"in MoL_PostStep after wave_updateMatter"},
  ConditionalOnKeyword -> {"couple_matter", "never"},
  Equations ->
  {
    rho -> 0,
    Si[li] -> 0,
    Sij[li,lj] -> 0,
    S -> 0
  }
};

(**************************************************************************************)
(* Scalar wave equation *)
(**************************************************************************************)
waveCalcADM[PD_, PDadvect_] =
{
  Name -> fnPrefix <> "_wave_adm",
  Schedule -> {"in MoL_CalcRHS as wave"},
  ConditionalOnKeyword -> {"decomp_type", "ADM"},
  Shorthands -> shorthands,
  Where -> Interior,
  Equations -> 
  {
    
     detg -> gDet,
     invdetg -> 1 / detg,
     gInv[ua,ub] -> invdetg gDet MatrixInverse[g[ua,ub]],

     (* Shorthands *)

       
     half -> 0.5,
     DDsigma[li,lj] -> PD[swSigma,li,lj], 
     Ddetg[li] -> detg gInv[uj,uk] PD[g[lj,lk],li],
     ADVsigma -> beta[uk] PD[swSigma,lk],
     
     (* Shift advection *)
     dir[ui] -> Sign[beta[ui]],

     dot[swSigma] -> alpha swEta + ADVsigma,
     dot[swEta] -> alpha * gInv[ui,uj] DDsigma[li,lj] + alpha *( PD[gInv[ui, uj],li] + gInv[ui,uj] Ddetg[li]/ (2 detg) ) PD[swSigma,lj] + alpha K swEta + gInv[ui,uj] PD[alpha,li] PD[swSigma,lj] 
  }
};

waveCalcBSSN[PD_, PDadvect_] =
{
  Name -> fnPrefix <> "_wave_bssn",
  Schedule -> {"in MoL_CalcRHS as wave"},
  ConditionalOnKeyword -> {"decomp_type", "BSSN"},
  Shorthands -> shorthands,
  Where -> Interior,
  Equations -> 
  { 
    deth -> hDet,
    hInv[ua,ub] -> hDet/deth MatrixInverse[h[ua,ub]],
    ADVsigma ->  beta[uk] PD[swSigma,lk],
    half -> 0.5,
    ADVeta ->  beta[uk] PD[swEta,lk],
    em4swPhi -> Exp[-4 swPhi],

    dir[ui] -> Sign[beta[ui]],

    dot[swSigma] -> alpha swEta + ADVsigma,
    dot[swEta] -> alpha em4swPhi * ( PD[hInv[ui,uj],li] PD[swSigma,lj] + hInv[ui,uj] PD[swSigma,li,lj] ) + 2 alpha em4swPhi hInv[ui,uj] PD[swPhi,li] PD[swSigma,lj] + alpha K swEta + em4swPhi hInv[ui,uj] PD[alpha,li] PD[swSigma,lj]
  }
};
(* In these initial data bits I multiply by em4swPhi to get the physical field *)
setGaussianID =
{
  Name -> fnPrefix <> "_set_gauss_id",
  Schedule -> {"at INITIAL"},
  ConditionalOnKeyword -> {"matter_initial_data", "gaussian"},
  AllowedSymbols -> {gAmp, gLambda},
  Shorthands -> shorthands,
  Equations ->
  {
    swSigma -> 0,
    swEta -> -4 gAmp gLambda Exp[-(x^2 + y^2 +z^2) gLambda]
  }
};

setPhysGaussianID =
{
  Name -> fnPrefix <> "_set_phys_gauss_id",
  Schedule -> {"at INITIAL"},
  ConditionalOnKeyword -> {"matter_initial_data", "phys_gaussian"},
  AllowedSymbols -> {gAmp, gLambda},
  Shorthands -> shorthands,
  Equations ->
  {
    
    swSigma -> 2 gAmp Exp[-(x^2+y^2+z^2) gLambda],
    swEta -> 0
  }
};


(**************************************************************************************)
(* Boundary conditions *)
(**************************************************************************************)

BoundaryParam = 
{
  Name -> "apply_radiative_bcs",
  Default -> "never",
  AllowedValues -> {"always", "never"}
};

DecompType =
{
  Name -> "decomp_type",
  Default -> "BSSN",
  AllowedValues -> {"BSSN","ADM"}
};

boundaryCalc =
{
  Name -> fnPrefix <> "boundary",
  Schedule -> {"in MoL_RHSBoundaries"},
  ConditionalOnKeyword -> {"apply_radiative_bcs", "always"},
  Where -> Boundary,
  Shorthands -> shorthands,
  Equations ->
  {
    (* n is the inward pointing normal, needed for dir in the 
       definition of the one-sided stencil. The boundary conditions
       then require a change of sign. *)
    n1 -> -x/r, 
    n2 -> -y/r,
    n3 -> -z/r,

    dir[ui] -> Sign[n[ui]],

    (* \partial_t u = - (u - u_0) / r - \partial_r u
       for u_0 some background solution. In this case, Minkowski in Cartesian 
       coordinates.  See: gr-qc/0203102 *)
    
    dot[swSigma] -> -(swSigma - 0) / r + n[uk] PDonesided2nd[swSigma, lk],
    dot[swEta] -> -(swEta - 0) / r + n[uk] PDonesided2nd[swEta, lk] 
  }
};

(**************************************************************************************)
(* Construct the thorn *)
(**************************************************************************************)

fdOrderParam = 
{
  Name -> "fd_order",
  Default -> "4th",
  AllowedValues -> {"2nd", "4th", "full4th"}
};

calculations = 
{
  initGFsCalc,
  admToADMCalc,
  admToBSSNCalc,
  killMatterTerms,
  setGaussianID,
  setPhysGaussianID,

  waveCalcADM[PDstandard4th, PDlopsided4th],
  waveCalcBSSN[PDstandard4th, PDlopsided4th],
  setMatterTermsADM[PDstandard4th, PDlopsided4th],
  setMatterTermsBSSN[PDstandard4th, PDlopsided4th],
  
  boundaryCalc
  (*energyCalc*)
};

keywordParameters = 
{
  fdOrderParam,
  BoundaryParam,
  DecompType,
  MatterCoupleParam
};

inheritedRealParameters = 
{
  "MatterBase::gAmp",
  "MatterBase::gLambda"
};

inheritedKeywordParameters =
{
  "MatterBase::matter_initial_data"
};

CreateKrancThornTT[groups,".", "ScalarField",
  Calculations -> calculations,
  DeclaredGroups -> declaredGroupNames,
  PartialDerivatives -> derivatives,
  RealParameters -> realParameters,
  KeywordParameters -> keywordParameters,
  InheritedImplementations -> {"admbase","MatterBase"},
  InheritedRealParameters -> inheritedRealParameters,
  InheritedKeywordParameters -> inheritedKeywordParameters
  ];
