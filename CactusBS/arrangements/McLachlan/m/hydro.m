$Path = Join[$Path, {"../../../repos/Kranc/Tools/CodeGen",
                     "../../../repos/Kranc/Tools/MathematicaMisc"}];

Get["KrancThorn`"];

(*SetDebugLevel[InfoFull];*)

SetEnhancedTimes[False];

(****************************************************************************
 Derivatives
****************************************************************************)

derivOrder = 2;

KD = KroneckerDelta;

derivatives =
{
  PDstandardNth[i_]     -> StandardCenteredDifferenceOperator[1,derivOrder/2,i],
  PDstandardNth[i_, i_] -> StandardCenteredDifferenceOperator[2,derivOrder/2,i],
  PDstandardNth[i_, j_] -> StandardCenteredDifferenceOperator[1,derivOrder/2,i]
                           StandardCenteredDifferenceOperator[1,derivOrder/2,j]
};

FD = PDstandardNth;
ResetJacobians;
DefineJacobian[PD, FD, KD, Zero3];

(****************************************************************************
 Tensors 
****************************************************************************)

(* Register all the tensors that will be used with TensorTools *)
Map[DefineTensor, 
{
  rho, vel, eps, press, vol,
  mass, mom, ene,
  massflux, momflux, eneflux
}];

(* Determinants of the metrics in terms of their components
  (Mathematica symbolic expressions) *)
gDet = Det[MatrixOfComponents[g[la,lb]]];

(****************************************************************************
 Groups
****************************************************************************)

SetGroupTimelevels[g_,tl_] = Join[g, {Timelevels -> tl}];

evolvedGroups = {
  CreateGroupFromTensor [mass   ],             (* mass *)
  CreateGroupFromTensor [mom[la]],             (* momentum *)
  CreateGroupFromTensor [ene    ]              (* total energy *)
};
evaluatedGroups = {
  CreateGroupFromTensor [rho    ],             (* mass density *)
  CreateGroupFromTensor [vel[ua]],             (* velocity *)
  CreateGroupFromTensor [eps    ],             (* specific internal energy *)
  CreateGroupFromTensor [press  ],             (* pressure *)
  CreateGroupFromTensor [vol    ],             (* volume *)
  CreateGroupFromTensor [massflux[ua]],        (* mass flux *)
  CreateGroupFromTensor [momflux[la,ub]],      (* momentum flux *)
  CreateGroupFromTensor [eneflux[ua]]          (* energy flux *)
};

declaredGroups = Join [evolvedGroups, evaluatedGroups];
declaredGroupNames = Map [First, declaredGroups];

groups = Join[declaredGroups];

(******************************************************************************)
(* Initial data *)
(******************************************************************************)

vacuumCalc =
{
  Name -> "hydro_vacuum",
  Schedule -> {"IN ADMBase_InitialData"},
  ConditionalOnKeyword -> {"initial_data", "vacuum"},
  Equations -> 
  {
    rho     -> 0,
    vel[ua] -> 0,
    eps     -> 0
  }
};

soundWaveCalc =
{
  Name -> "hydro_soundWave",
  Schedule -> {"IN ADMBase_InitialData"},
  ConditionalOnKeyword -> {"initial_data", "sound wave"},
  Equations -> 
  {
    rho     -> 1.0,
    vel[ua] -> A Sin[2 Pi x / L],
    eps     -> 1.0
  }
};

(******************************************************************************)
(* Convert from primitive to conserved variables *)
(******************************************************************************)

prim2conCalc =
{
  Name -> "hydro_prim2con",
  Schedule -> {"AT initial AFTER ADMBase_PostInitial"},
  Equations -> 
  {
    vol     -> h^3,
    mass    -> vol rho,
    mom[la] -> mass Euc[la,lb] vel[ub],
    ene     -> (1/2) mass Euc[la,lb] vel[ub] vel[ua] + mass eps
  }
};

(******************************************************************************)
(* Convert from conserved to primitive variables *)
(******************************************************************************)

con2primCalc =
{
  Name -> "hydro_con2prim",
  Schedule -> {"IN hydro_con2primGroup"},
  Equations -> 
  {
    rho     -> mass / vol,
    vel[ua] -> Euc[ua,ub] mom[lb] / mass,
    eps     -> ene / mass - (1/2) Euc[la,lc] vel[ua] vel[uc],
    
    press   -> Gamma rho eps
  }
};

(******************************************************************************)
(* Evolution equations *)
(******************************************************************************)

fluxCalc =
{
  Name -> "hydro_fluxes",
  Schedule -> {"IN hydro_evolCalcGroup"},
  Equations -> 
  {
    massflux[ua]    -> mass vel[ua],
    momflux [la,ub] -> mom[la] vel[ub] + KD[la,ub] press / vol,
    eneflux [ua]    -> (ene + press / vol) vel[ua]
  }
};

evolCalc =
{
  Name -> "hydro_RHS",
  Schedule -> {"IN hydro_evolCalcGroup AFTER hydro_fluxes"},
  Where -> Interior,
  Equations -> 
  {
    (* dt rho + div rho v = 0 *)
    (* dt m + div m v = 0 *)
    dot[mass] -> - PD[massflux[ua],la],
    
    (* dt pi + div (pi v + P) = 0 *)
    (* dt p + div (p v + P/V) = 0 *)
    dot[mom[la]] -> - PD[momflux[la,ub],lb],
    
    (* dt tau + div (tau v + P) = 0 *)
    (* dt t + div (t v + P/V) = 0 *)
    dot[ene] -> - PD[eneflux[ua],la]
  }
};

(******************************************************************************)
(* Parameters *)
(******************************************************************************)

keywordParameters =
{
  {
    Name -> "initial_data",
    (* Visibility -> "restricted", *)
    (* Description -> "ddd", *)
    AllowedValues -> {"vacuum", "sound wave"},
    Default -> "vacuum"
  }
};

realParameters =
{
  {
    Name -> h,
    Description -> "grid spacing",
    Default -> 0.01
  },
  {
    Name -> A,
    Description -> "sound wave amplitude",
    Default -> 0.001
  },
  {
    Name -> L,
    Description -> "sound wave wavelength",
    Default -> 1
  },
  {
    Name -> alpha,
    Description -> "artificial viscosity coefficient",
    Default -> 0
  },
  {
    Name -> Gamma,
    Description -> "polytropic exponent",
    Default -> 4/3
  }
};

(******************************************************************************)
(* Construct the thorns *)
(******************************************************************************)

calculations =
{
  vacuumCalc,
  soundWaveCalc,
  prim2conCalc,
  con2primCalc,
  evolCalc
};

CreateKrancThornTT [groups, ".", "ML_hydro",
  Calculations -> calculations,
  DeclaredGroups -> declaredGroupNames,
  PartialDerivatives -> derivatives,
  UseLoopControl -> True,
  KeywordParameters -> keywordParameters,
  RealParameters -> realParameters
];
