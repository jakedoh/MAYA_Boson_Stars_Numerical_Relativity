$Path = Join[$Path, {"../../../repos/Kranc/Tools/CodeGen",
                     "../../../repos/Kranc/Tools/MathematicaMisc"}];

Get["KrancThorn`"];

SetEnhancedTimes[False];
SetSourceLanguage["C"];



(******************************************************************************)
(* Derivatives *)
(******************************************************************************)

KD = KroneckerDelta;

partialDerivatives =
{
  PDstandardNth[i_]    -> StandardCenteredDifferenceOperator[1,derivOrder/2,i],
  PDstandardNth[i_,i_] -> StandardCenteredDifferenceOperator[2,derivOrder/2,i],
  PDstandardNth[i_,j_] -> StandardCenteredDifferenceOperator[1,derivOrder/2,i] *
                          StandardCenteredDifferenceOperator[1,derivOrder/2,j],
  PDdissipationNth[i_] ->
    spacing[i]^(derivOrder+1) / 2^(derivOrder+2) *
    StandardCenteredDifferenceOperator[derivOrder+2,derivOrder/2+1,i],
};

ResetJacobians;

DefineJacobian[PD, FD, J, dJ];



(******************************************************************************)
(* Tensors *)
(******************************************************************************)

DefineTensor[normal];
DefineTensor[tangentA];
DefineTensor[tangentB];
DefineTensor[dir];

AssertSymmetricIncreasing[admg[lalb]];
AssertSymmetricIncreasing[admK[lalb]];
AssertSymmetricIncreasing [dJ[ua,lb,lc], lb, lc];
AssertSymmetricIncreasing [G[ua,lb,lc], lb, lc];
AssertSymmetricIncreasing [Gtl[la,lb,lc], lb, lc];
AssertSymmetricIncreasing [Gt[ua,lb,lc], lb, lc];
AssertSymmetricIncreasing [gK[la,lb,lc], la, lb];

DefineConnection [CD, PD, G];
DefineConnection [CDt, PD, Gt];



(******************************************************************************)
(* Expressions *)
(******************************************************************************)

pi = N[Pi,40];

detgExpr  = Det [MatrixOfComponents [g [la,lb]]];
ddetgExpr[la_] =
  Sum [D[Det[MatrixOfComponents[g[la, lb]]], X] PD[X, la],
       {X, Union[Flatten[MatrixOfComponents[g[la, lb]]]]}];

detgtExpr = Det [MatrixOfComponents [gt[la,lb]]];
ddetgtExpr[la_] =
  Sum [D[Det[MatrixOfComponents[gt[la, lb]]], X] PD[X, la],
       {X, Union[Flatten[MatrixOfComponents[gt[la, lb]]]]}];



(******************************************************************************)
(* Groups *)
(******************************************************************************)

declaredGroups = {
  SetGroupName [CreateGroupFromTensor [phi      ], "log_confac"],
  SetGroupName [CreateGroupFromTensor [gt[la,lb]], "metric"    ],
  SetGroupName [CreateGroupFromTensor [Xt[ua]   ], "Gamma"     ]
};

declaredGroupNames = Map[First, declaredGroups];

extraGroups =
  {{"Grid::coordinates", {x, y, z, r}},
   {"ADMBase::metric",  {gxx, gxy, gxz, gyy, gyz, gzz}},
   {"ADMBase::curv",    {kxx, kxy, kxz, kyy, kyz, kzz}}}
};

allGroups = Join[declaredGroups, extraGroups];



(******************************************************************************)
(* Initial data *)
(******************************************************************************)

initialCalc =
{
  Name -> BSSN <> "_Minkowski",
  Schedule -> {"IN ADMBase_InitialData"},
  ConditionalOnKeyword -> {"my_initial_data", "Minkowski"},
  Equations -> 
  {
    phi       -> IfThen [conformalMethod, 1, 0],
    gt[la,lb] -> KD[la,lb],
    trK       -> 0,
    At[la,lb] -> 0,
    Xt[ua]    -> 0,
    alpha     -> 1,
    A         -> 0,
    beta[ua]  -> 0,
    B[ua]     -> 0
  }
};



(******************************************************************************)
(* Convert from ADMBase *)
(******************************************************************************)

convertFromADMBaseCalc =
{
  Name -> BSSN <> "_convertFromADMBase",
  Schedule -> {"AT initial AFTER ADMBase_PostInitial"},
  ConditionalOnKeyword -> {"my_initial_data", "ADMBase"},
  Shorthands -> {g[la,lb], detg, gu[ua,ub], em4phi},
  Equations -> 
  {
    g[la,lb]  -> admg[la,lb],
    detg      -> detgExpr,
    gu[ua,ub] -> 1/detg detgExpr MatrixInverse [g[ua,ub]],
    
    phi       -> IfThen [conformalMethod, detg^(-1/6), Log[detg]/12],
    em4phi    -> IfThen [conformalMethod, phi^2, Exp[-4 phi]],
    gt[la,lb] -> em4phi g[la,lb],
    
    trK       -> gu[ua,ub] admK[la,lb],
    At[la,lb] -> em4phi (admK[la,lb] - (1/3) g[la,lb] trK),
    
    alpha     -> admalpha,
    
    beta[ua]  -> admbeta[ua]
  }
};

convertFromADMBaseGammaCalc =
{
  Name -> BSSN <> "_convertFromADMBaseGamma",
  Schedule -> {"AT initial AFTER " <> BSSN <> "_convertFromADMBase"},
  ConditionalOnKeyword -> {"my_initial_data", "ADMBase"},
  (*
  Where -> InteriorNoSync,
  *)
  (* Do not synchronise right after this routine; instead, synchronise
     after extrapolating *)
  Where -> Interior,
  (* Synchronise after this routine, so that the refinement boundaries
     are set correctly before extrapolating.  (We will need to
     synchronise again after extrapolating because extrapolation does
     not fill ghost zones, but this is irrelevant here.)  *)
  Shorthands -> {dir[ua],
                 detgt, gtu[ua,ub], Gt[ua,lb,lc], theta},
  Equations -> 
  {
    dir[ua] -> Sign[beta[ua]],
    
    detgt        -> 1 (* detgtExpr *),
    gtu[ua,ub]   -> 1/detgt detgtExpr MatrixInverse [gt[ua,ub]],
    Gt[ua,lb,lc] -> 1/2 gtu[ua,ud]
                    (PD[gt[lb,ld],lc] + PD[gt[lc,ld],lb] - PD[gt[lb,lc],ld]),
    Xt[ua] -> gtu[ub,uc] Gt[ua,lb,lc],
    
(*
    A -> - admdtalpha / (harmonicF alpha^harmonicN) (LapseAdvectionCoeff - 1),
*)
    (* If LapseACoeff=0, then A is not evolved, in the sense that it
       does not influence the time evolution of other variables.  *)
    A -> IfThen [LapseACoeff != 0,
                 1 / (- harmonicF alpha^harmonicN)
                 (+ admdtalpha
                  - LapseAdvectionCoeff beta[ua] PDua[alpha,la]
                  - LapseAdvectionCoeff Abs[beta[ua]] PDus[alpha,la]),
                 0],
    
    theta -> thetaExpr,
    
    (* If ShiftBCoeff=0 or theta ShiftGammaCoeff=0, then B^i is not
       evolved, in the sense that it does not influence the time
       evolution of other variables.  *)
    B[ua] -> IfThen [ShiftGammaCoeff ShiftBCoeff != 0,
                     1 / (theta ShiftGammaCoeff)
                     (+ admdtbeta[ua]
                      - ShiftAdvectionCoeff beta[ub] PDua[beta[ua],lb]
                      - ShiftAdvectionCoeff Abs[beta[ub]] PDus[beta[ua],lb]),
                     0]
  }
};



(******************************************************************************)
(* Implementations *)
(******************************************************************************)

inheritedImplementations = {"ADMBase", "TmunuBase"};



(******************************************************************************)
(* Parameters *)
(******************************************************************************)

keywordParameters =
{
  {
    Name -> "my_initial_data",
    (* Visibility -> "restricted", *)
    (* Description -> "ddd", *)
    AllowedValues -> {"ADMBase", "Minkowski"},
    Default -> "ADMBase"
  },
  {
    Name -> "my_initial_boundary_condition",
    Visibility -> "restricted",
    (* Description -> "ddd", *)
    AllowedValues -> {"none"},
    Default -> "none"
  }
};

realParameters =
{
  {
    Name -> LapseACoeff,
    Description -> "Whether to evolve A in time",
    Default -> 0
  },
  {
    Name -> harmonicF,
    Description -> "d/dt alpha = - f alpha^n K   (harmonic=1, 1+log=2)",
    Default -> 1
  }
};



(******************************************************************************)
(* Construct the thorns *)
(******************************************************************************)

calculations =
{
  initialCalc,
  convertFromADMBaseCalc,
  convertFromADMBaseGammaCalc
};

CreateKrancThornTT [allGroups, ".", BSSN,
  Calculations -> calculations,
  DeclaredGroups -> declaredGroupNames,
  PartialDerivatives -> derivatives,
  EvolutionTimelevels -> evolutionTimelevels,
  DefaultEvolutionTimelevels -> 3,
  UseLoopControl -> True,
  InheritedImplementations -> inheritedImplementations,
  InheritedKeywordParameters -> inheritedKeywordParameters,
  ExtendedKeywordParameters -> extendedKeywordParameters,
  KeywordParameters -> keywordParameters,
  IntParameters -> intParameters,
  RealParameters -> realParameters
];
