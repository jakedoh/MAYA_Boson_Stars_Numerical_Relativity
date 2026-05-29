
SetEnhancedTimes[False];

(******************************************************************************)
(* Options *)
(******************************************************************************)

(* useJacobian: True or False *)
useJacobian = True;

(* timelevels: 2 or 3
   (keep this at 3; this is better chosen with a run-time parameter) *)
evolutionTimelevels = 3;

(* matter: 0 or 1 *)
addMatter = 1;



prefix = "ML_";
suffix =
  (* If [useJacobian, "_MP", ""] <> *)
  (* If [derivOrder!=4, "_O" <> ToString[derivOrder], ""] <> *)
  If [evolutionTimelevels!=3, "_TL" <> ToString[evolutionTimelevels], ""] <>
  (* If [addMatter!=0, "_M", ""] <> *)
  "";

ADM = prefix <> "ADM" <> suffix;

(******************************************************************************)
(* Derivatives *)
(******************************************************************************)

KD = KroneckerDelta;

derivatives =
{
  PDstandardNth[i_]     -> StandardCenteredDifferenceOperator[1,fdOrder/2,i],
  PDstandardNth[i_, i_] -> StandardCenteredDifferenceOperator[2,fdOrder/2,i],
  PDstandardNth[i_, j_] -> StandardCenteredDifferenceOperator[1,fdOrder/2,i]
                           StandardCenteredDifferenceOperator[1,fdOrder/2,j]
};

PD = PDstandardNth;



(******************************************************************************)
(* Tensors *)
(******************************************************************************)

(* Register the tensor quantities with the TensorTools package *)
Map [DefineTensor,
     {g, K, alpha, beta, H, M, detg, gu, G, R, trR, Km, trK,
      T00, T0, T, rho, S}];

Map [AssertSymmetricIncreasing,
     {g[la,lb], K[la,lb], R[la,lb],
      T[la,lb]}];
AssertSymmetricIncreasing [G[ua,lb,lc], lb, lc];

DefineConnection [CD, PD, G];

Map [DefineTensor,
     {gxx, gxy, gxz, gyy, gyz, gzz,
      kxx, kxy, kxz, kyy, kyz, kzz,
      alp,
      dtalp,
      betax, betay, betaz,
      dtbetax, dtbetay, dtbetaz,
      eTtt,
      eTtx, eTty, eTtz,
      eTxx, eTxy, eTxz, eTyy, eTyz, eTzz}];

(******************************************************************************)
(* Expressions *)
(******************************************************************************)

detgExpr  = Det [MatrixOfComponents [g [la,lb]]];
ddetgExpr[la_] =
  Sum [D[Det[MatrixOfComponents[g[la, lb]]], X] PD[X, la],
       {X, Union[Flatten[MatrixOfComponents[g[la, lb]]]]}];

(******************************************************************************)
(* Groups *)
(******************************************************************************)

SetGroupTimelevels[g_,tl_] = Join[g, {Timelevels -> tl}];

evolvedGroups =
  {SetGroupName [CreateGroupFromTensor [g[la,lb]], prefix <> "metric"],
   SetGroupName [CreateGroupFromTensor [K[la,lb]], prefix <> "curv"  ],
   SetGroupName [CreateGroupFromTensor [alpha   ], prefix <> "lapse" ],
   SetGroupName [CreateGroupFromTensor [beta[ua]], prefix <> "shift" ]};
evaluatedGroups =
  {SetGroupName [CreateGroupFromTensor [H    ], prefix <> "Ham"],
   SetGroupName [CreateGroupFromTensor [M[la]], prefix <> "mom"]};

declaredGroups = Join [evolvedGroups, evaluatedGroups];
declaredGroupNames = Map [First, declaredGroups];



extraGroups =
  {{"ADMBase::metric",   {gxx, gxy, gxz, gyy, gyz, gzz}},
   {"ADMBase::curv",     {kxx, kxy, kxz, kyy, kyz, kzz}},
   {"ADMBase::lapse",    {alp}},
   {"ADMBase::dtlapse",  {dtalp}},
   {"ADMBase::shift",    {betax, betay, betaz}},
   {"ADMBase::dtshift",  {dtbetax, dtbetay, dtbetaz}},
   {"TmunuBase::stress_energy_scalar", {eTtt}},
   {"TmunuBase::stress_energy_vector", {eTtx, eTty, eTtz}},
   {"TmunuBase::stress_energy_tensor", {eTxx, eTxy, eTxz, eTyy, eTyz, eTzz}}
};



groups = Join [declaredGroups, extraGroups];

(******************************************************************************)
(* Initial data *)
(******************************************************************************)

initialCalc =
{
  Name -> ADM <> "_Minkowski",
  Schedule -> {"IN ADMBase_InitialData"},
  ConditionalOnKeyword -> {"my_initial_data", "Minkowski"},
  Equations -> 
  {
    g[la,lb] -> KD[la,lb],
    K[la,lb] -> 0,
    alpha    -> 1,
    beta[ua] -> 0
  }
};

(******************************************************************************)
(* Convert from ADMBase *)
(******************************************************************************)

convertFromADMBaseCalc =
{
  Name -> ADM <> "_convertFromADMBase",
  Schedule -> {"AT initial AFTER ADMBase_PostInitial"},
  ConditionalOnKeyword -> {"my_initial_data", "ADMBase"},
  Equations -> 
  {
    g11   -> gxx,
    g12   -> gxy,
    g13   -> gxz,
    g22   -> gyy,
    g23   -> gyz,
    g33   -> gzz,
    K11   -> kxx,
    K12   -> kxy,
    K13   -> kxz,
    K22   -> kyy,
    K23   -> kyz,
    K33   -> kzz,
    (* TODO: this is incomplete; it ignores dtalp and dtbeta^i *)
    alpha -> alp,
    beta1 -> betax,
    beta2 -> betay,
    beta3 -> betaz
  }
};

(******************************************************************************)
(* Convert to ADMBase *)
(******************************************************************************)

convertToADMBaseCalc =
{
  Name -> ADM <> "_convertToADMBase",
  Schedule -> {"IN MoL_PostStep AFTER " <> ADM <> "_ApplyBCs"},
  Equations -> 
  {
    gxx     -> g11,
    gxy     -> g12,
    gxz     -> g13,
    gyy     -> g22,
    gyz     -> g23,
    gzz     -> g33,
    kxx     -> K11,
    kxy     -> K12,
    kxz     -> K13,
    kyy     -> K22,
    kyz     -> K23,
    kzz     -> K33,
    (* TODO: this is wrong; it sets dtalp and dtbeta^i incorrectly *)
    alp     -> alpha,
    dtalp   -> 0,
    betax   -> beta1,
    betay   -> beta2,
    betaz   -> beta3,
    dtbetax -> 0,
    dtbetay -> 0,
    dtbetaz -> 0
  }
};

(******************************************************************************)
(* Evolution equations *)
(******************************************************************************)

evolCalc =
{
  Name -> ADM <> "_RHS",
  Schedule -> {"IN MoL_CalcRHS", "AT analysis"},
  Where -> Interior,
  Shorthands -> {detg, gu[ua,ub], G[ua,lb,lc], R[la,lb], Km[ua,lb], trK},
  Equations -> 
  {
    detg -> detgExpr,
    gu[ua,ub] -> 1/detg detgExpr MatrixInverse [g[ua,ub]],
    G[ua,lb,lc] -> 1/2 gu[ua,ud]
                   (PD[g[lb,ld],lc] + PD[g[lc,ld],lb] - PD[g[lb,lc],ld]),
    R[la,lb] -> G[uc,ld,la] G[ud,lc,lb] - G[uc,la,lb] G[ud,lc,ld]
                + 1/2 gu[uc,ud] (- PD[g[lc,ld],la,lb] + PD[g[lc,la],ld,lb]
                                 - PD[g[la,lb],lc,ld] + PD[g[ld,lb],lc,la]),
    Km[ua,lb] -> gu[ua,uc] K[lc,lb],
    trK -> Km[ua,la],

    dot[g[la,lb]] -> -2 alpha K[la,lb]
                     + Lie[g[la,lb], beta],
    dot[K[la,lb]] -> - CD[alpha,la,lb]
                     + alpha * (R[la,lb] + K[la,lb] trK - 2 K[la,lc] Km[uc,lb])
                     + Lie[K[la,lb], beta],
    dot[alpha]    -> 0,
    dot[beta[ua]] -> 0
  }
};

(******************************************************************************)
(* Boundary conditions *)
(******************************************************************************)

boundaryCalc =
{
  Name -> ADM <> "_boundary",
  Schedule -> {"IN MoL_PostStep"},
  ConditionalOnKeyword -> {"my_boundary_condition", "Minkowski"},
  Where -> BoundaryWithGhosts,
  Equations -> 
  {
    g[la,lb] -> KD[la,lb],
    K[la,lb] -> 0,
    alpha    -> 1,
    beta[ua] -> 0
  }
};

(******************************************************************************)
(* Constraint equations *)
(******************************************************************************)

constraintsCalc =
{
  Name -> ADM <> "_constraints",
  Schedule -> {"AT analysis"},
  Where -> Interior,
  Shorthands -> {detg, gu[ua,ub], G[ua,lb,lc], R[la,lb], trR, Km[ua,lb], trK},
  Equations -> 
  {
    detg -> detgExpr,
    gu[ua,ub] -> 1/detg detgExpr MatrixInverse[g[ua,ub]],
    G[ua,lb,lc] -> 1/2 gu[ua,ud]
                   (PD[g[lb,ld],lc] + PD[g[lc,ld],lb] - PD[g[lb,lc],ld]),
    R[la,lb] -> gu[us,ur](G[um,la,lr] G[uk,ls,lb] g[lk,lm] - G[um,la,lb] G[uk,ls,lr] g[lk,lm])
                + 1/2 gu[uc,ud] (- PD[g[lc,ld],la,lb] + PD[g[lc,la],ld,lb]
                                 - PD[g[la,lb],lc,ld] + PD[g[ld,lb],lc,la]),
    trR -> R[la,lb] gu[ua,ub],
    Km[ua,lb] -> gu[ua,uc] K[lc,lb],
    trK -> Km[ua,la],

    H -> trR - Km[ua,lb] Km[ub,la] + trK^2,
    M[la] -> gu[ub,uc] (CD[K[lc,la], lb] - CD[K[lc,lb], la])
  }
};

constraintsBoundaryCalc =
{
  Name -> ADM <> "_constraints_boundary",
  Schedule -> {"AT analysis AFTER " <> ADM <> "_constraints"},
  Where -> BoundaryWithGhosts,
  Equations -> 
  {
    H     -> 0,
    M[la] -> 0
  }
};

(******************************************************************************)
(* Implementations *)
(******************************************************************************)

inheritedImplementations =
  Join[{"ADMBase"},
       If [addMatter!=0, {"TmunuBase"}, {}]];

(******************************************************************************)
(* Parameters *)
(******************************************************************************)

extendedKeywordParameters =
{
  {
    Name -> "ADMBase::evolution_method",
    AllowedValues -> {BSSN}
  },
  {
    Name -> "ADMBase::lapse_evolution_method",
    AllowedValues -> {BSSN}
  },
  {
    Name -> "ADMBase::shift_evolution_method",
    AllowedValues -> {BSSN}
  }
};

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
    Name -> "my_boundary_condition",
    (* Visibility -> "restricted", *)
    (* Description -> "ddd", *)
    AllowedValues -> {"none", "Minkowski"},
    Default -> "none"
  }
};

intParameters =
{
  {
    Name -> fdOrder,
    Default -> 4,
    AllowedValues -> {2,4,6,8}
  }
};

(******************************************************************************)
(* Construct the thorns *)
(******************************************************************************)

calculations =
{
  initialCalc,
  convertFromADMBaseCalc,
  evolCalc,
  boundaryCalc,
  convertToADMBaseCalc,
  constraintsCalc,
  constraintsBoundaryCalc
};

CreateKrancThornTT [groups, ".", ADM,
  Calculations -> calculations,
  DeclaredGroups -> declaredGroupNames,
  PartialDerivatives -> derivatives,
  EvolutionTimelevels -> evolutionTimelevels,
  UseJacobian -> useJacobian,
  UseLoopControl -> True,
  UseVectors -> True,
  InheritedImplementations -> inheritedImplementations,
  KeywordParameters -> keywordParameters,
  IntParameters -> intParameters
];
