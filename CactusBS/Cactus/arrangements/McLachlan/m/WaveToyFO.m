
SetEnhancedTimes[False];

(******************************************************************************)
(* Derivatives *)
(******************************************************************************)

derivOrder = 4;

derivatives =
{
  PDstandardNth[i_]     -> StandardCenteredDifferenceOperator[1,derivOrder/2,i],
  PDstandardNth[i_, i_] -> StandardCenteredDifferenceOperator[2,derivOrder/2,i],
  PDstandardNth[i_, j_] -> StandardCenteredDifferenceOperator[1,derivOrder/2,i]
                           StandardCenteredDifferenceOperator[1,derivOrder/2,j]
(*
  PDstandardNth[i_, i_, i_] ->
    StandardCenteredDifferenceOperator[3,derivOrder/2,i],
  PDstandardNth[i_, i_, j_] ->
    StandardCenteredDifferenceOperator[2,derivOrder/2,i]
    StandardCenteredDifferenceOperator[1,derivOrder/2,j],
  PDstandardNth[i_, j_, i_] ->
    StandardCenteredDifferenceOperator[2,derivOrder/2,i]
    StandardCenteredDifferenceOperator[1,derivOrder/2,j],
  PDstandardNth[j_, i_, i_] ->
    StandardCenteredDifferenceOperator[2,derivOrder/2,i]
    StandardCenteredDifferenceOperator[1,derivOrder/2,j],
  PDstandardNth[i_, j_, k_] ->
    StandardCenteredDifferenceOperator[1,derivOrder/2,i]
    StandardCenteredDifferenceOperator[1,derivOrder/2,j]
    StandardCenteredDifferenceOperator[1,derivOrder/2,k]
*)
};

PD = PDstandardNth;

(* timelevels *)
evolutionTimelevels = 2;

KD = KroneckerDelta;

(******************************************************************************)
(* Tensors *)
(******************************************************************************)

(* Register the tensor quantities with the TensorTools package *)
Map [DefineTensor, {u, rho, v, w}];

(******************************************************************************)
(* Groups *)
(******************************************************************************)

evolvedGroups =
  {SetGroupName [CreateGroupFromTensor [u    ], "WT_u"  ],
   SetGroupName [CreateGroupFromTensor [v[la]], "WT_v"  ],
   SetGroupName [CreateGroupFromTensor [rho  ], "WT_rho"]};
evaluatedGroups =
  {AddGroupTag[
      SetGroupName [CreateGroupFromTensor [w[ua]], "WT_w"],
      "Prolongation"->"None"]};

declaredGroups = Join [evolvedGroups, evaluatedGroups];
declaredGroupNames = Map [First, declaredGroups];

groups = declaredGroups;

(******************************************************************************)
(* Initial data *)
(******************************************************************************)

initialCalc =
{
  Name -> "WTFO_Gaussian",
  Schedule -> {"AT initial"},
  (* Where -> Boundary, *)
  (* Where -> Interior, *)
  Equations -> 
  {
    u -> 0,
    v[la] -> 0,
    rho -> 0
  }
};

(******************************************************************************)
(* Evolution equations *)
(******************************************************************************)

evolCalc =
{
  Name -> "WTFO_RHS",
  Schedule -> {"IN MoL_CalcRHS", "AT analysis"},
  Where -> Interior,
  Equations -> 
  {
    dot[u] -> rho,
    dot[rho] -> KD[ua,ub] PD[v[la],lb],
    dot[v[la]] -> PD[rho,la]
  }
};

(******************************************************************************)
(* Constraint equations *)
(******************************************************************************)

constraintsCalc =
{
  Name -> "WTFO_constraints",
  Schedule -> {"AT analysis"},
  Where -> Interior,
  Equations -> 
  {
    w[ua] -> Eps[ua,ub,uc] PD[v[lb],lc]
  }
};

(******************************************************************************)
(* Construct the thorns *)
(******************************************************************************)

calculations = 
{
  initialCalc,
  evolCalc,
  constraintsCalc
};

CreateKrancThornTT [groups, ".", "ML_WaveToyFO",
  Calculations -> calculations,
  DeclaredGroups -> declaredGroupNames,
  PartialDerivatives -> derivatives,
  UseLoopControl -> True,
  EvolutionTimelevels -> evolutionTimelevels
];
