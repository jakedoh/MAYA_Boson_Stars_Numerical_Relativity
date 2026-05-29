
SetEnhancedTimes[False];

(******************************************************************************)
(* Derivatives *)
(******************************************************************************)

derivOrder = 4;

derivatives =
  {PDstandardNth[i_]     ->
   StandardCenteredDifferenceOperator[1,derivOrder/2,i],
   PDstandardNth[i_, i_] ->
   StandardCenteredDifferenceOperator[2,derivOrder/2,i],
   PDstandardNth[i_, j_] ->
   StandardCenteredDifferenceOperator[1,derivOrder/2,i]
   StandardCenteredDifferenceOperator[1,derivOrder/2,j]};

PD = PDstandardNth;

evolutionTimelevels = 4;

KD = KroneckerDelta;

(******************************************************************************)
(* Tensors *)
(******************************************************************************)

(* Register the tensor quantities with the TensorTools package *)
Map[DefineTensor, {u, rho, eps}];

(******************************************************************************)
(* Groups *)
(******************************************************************************)

evolvedGroups =
  {SetGroupName[CreateGroupFromTensor[u  ], "WT_u"  ],
   SetGroupName[CreateGroupFromTensor[rho], "WT_rho"]};
evaluatedGroups =
  {AddGroupTag[
      SetGroupName[CreateGroupFromTensor[eps], "WT_eps"],
      "Prolongation" -> "None"]};

declaredGroups = Join[evolvedGroups, evaluatedGroups];
declaredGroupNames = Map[First, declaredGroups];

extraGroups =
  {{"grid::coordinates", {x, y, z, r}}};

groups = declaredGroups;

(******************************************************************************)
(* Initial data *)
(******************************************************************************)

initialCalcGaussian =
  {Name                 -> "WT_Gaussian",
   Schedule             -> {"AT initial"},
   ConditionalOnKeyword -> {"initial_data", "Gaussian"},
   Equations            -> 
   {u   -> amplitude Exp[-(1/2) (r/width)^2],
    rho -> 0}};

initialCalcStanding =
  {Name                 -> "WT_Standing",
   Schedule             -> {"AT initial"},
   ConditionalOnKeyword -> {"initial_data", "Standing"},
   Shorthands           -> {kvec, omega},
   Equations            -> 
   {kvec  -> Pi / width,
    omega -> Sqrt[3 kvec^2],
    u     -> amplitude Cos[kvec x] Cos[kvec y] Cos[kvec z] Cos[omega t],
    rho   -> amplitude Cos[kvec x] Cos[kvec y] Cos[kvec z] Sin[omega t] (-omega)}};

(******************************************************************************)
(* Evolution equations *)
(******************************************************************************)

evolCalc =
  {Name      -> "WT_RHS",
   Schedule  -> {"IN MoL_CalcRHS" (*"AT analysis"*)},
   Where     -> Interior,
   Equations -> 
   {dot[u]   -> rho,
    dot[rho] -> KD[ua,ub] PD[u,la,lb]}};

(******************************************************************************)
(* Boundary conditions *)
(******************************************************************************)

boundaryCalc =
  {Name      -> "WT_Dirichlet",
   Schedule  -> {"IN MoL_CalcRHS", "AT analysis"},
   Where     -> Boundary,
   Equations -> 
   {dot[u]   -> 0,
    dot[rho] -> 0}};

(******************************************************************************)
(* Energy *)
(******************************************************************************)

energyCalc =
  {Name      -> "WT_Energy",
   Schedule  -> {"AT analysis"},
   Where     -> Interior,
   Equations -> 
   {eps -> 1/2 (rho^2 + KD[ua,ub] PD[u,la] PD[u,lb])}};

energyBoundaryCalc =
  {Name      -> "WT_EnergyBoundary",
   Schedule  -> {"AT analysis"},
   Where     -> Boundary,
   Equations -> 
   {eps -> 0}};

(******************************************************************************)
(* Parameters *)
(******************************************************************************)

keywordParameters =
  {{Name          -> "initial_data",
    AllowedValues -> {"Gaussian", "Standing"},
    Default       -> "Gaussian"}};

realParameters =
  {{Name        -> amplitude,
    Description -> "Amplitude of initial Gaussian",
    Default     -> 1},
   {Name        -> width,
    Description -> "Width of initial Gaussian",
    Default     -> 1}};

(******************************************************************************)
(* Construct the thorns *)
(******************************************************************************)

calculations = 
  {initialCalcGaussian, initialCalcStanding,
   evolCalc, boundaryCalc,
   energyCalc, energyBoundaryCalc};

CreateKrancThornTT[
  groups, ".", "ML_WaveToy",
  Calculations        -> calculations,
  DeclaredGroups      -> declaredGroupNames,
  PartialDerivatives  -> derivatives,
  UseLoopControl      -> True,
  UseVectors          -> True,
  EvolutionTimelevels -> evolutionTimelevels,
  KeywordParameters   -> keywordParameters,
  RealParameters      -> realParameters];

renameCalc[calc_, oldprefix_, newprefix_] := Module[
  {name, newname, newcalc},
  name = lookup[calc, Name];
  newname = StringReplace[name,
                          RegularExpression["^" <> oldprefix] -> newprefix];
  newcalc = mapReplace[calc, Name, newname];
  newcalc];

CreateKrancThornTT[
  groups, ".", "ML_WaveToy_CL",
  Calculations        -> Map[renameCalc[#, "WT_", "WT_CL_"]&, calculations],
  DeclaredGroups      -> declaredGroupNames,
  PartialDerivatives  -> derivatives,
  UseLoopControl      -> True,
  UseOpenCL           -> True,
  UseVectors          -> True,
  EvolutionTimelevels -> evolutionTimelevels,
  KeywordParameters   -> keywordParameters,
  RealParameters      -> realParameters];
