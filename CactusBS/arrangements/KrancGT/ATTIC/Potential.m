
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
  n, dir, Phi, Psi, DDPhi, energy, V
}];


(**************************************************************************************)
(* Groups *)
(**************************************************************************************)

(* Cactus group definitions *)

wave = Map[CreateGroupFromTensor, {Phi, Psi, energy}];
potential = Map[CreateGroupFromTensor, {V}];

declaredGroups = Join[wave,potential];
declaredGroupNames = Map[First, declaredGroups];

groups = Join[declaredGroups];

(**************************************************************************************)
(* Shorthands *)
(**************************************************************************************)

shorthands = 
{
  dir[ui], n[ui], DDPhi, half, zred, xmin, Vmin 
};

(**************************************************************************************)
(* Parameters *)
(**************************************************************************************)

realParameters = 
{
  mass, eps, x0, sigma, ell
};

(**************************************************************************************)
(* Initialize all grid functions to a recognizable value. These should be overwritten *)
(* later. *)
(**************************************************************************************)

INITVALUE = "100";

initScalarGFsCalc = 
{
  Name -> "init_scalargfs",
  Schedule -> {"at INITIAL"},
  Equations -> 
  {
    Phi -> INITVALUE,
    Psi -> INITVALUE,
    V -> INITVALUE
  }
};

(**************************************************************************************)
(* Initial data and potential *)
(**************************************************************************************)

IDParam =
{
  Name -> "initial_data",
  Default -> "pulse",
  AllowedValues -> {"pulse"}
};

setIDPulse = 
{
  Name -> "set_ID_Pulse",
  Schedule -> {"at INITIAL"},
  ConditionalOnKeyword -> {"initial_data", "pulse"},
  Shorthands -> shorthands,
  Equations ->
    {
     Phi -> Exp[(-(x-x0)^2)/sigma^2],
     Psi -> - 2 (x-x0) / sigma^2 Phi
   }
};

setPotential = 
{
  Name -> "set_Potential",
  Schedule -> {"at INITIAL"},
  Shorthands -> shorthands,
  Equations ->
    {
      xmin -> 1.5 * (ell*(ell+1)-1)/(ell*(ell+1)) + Sqrt[(1.5 * (ell*(ell+1)-1)/(ell*(ell+1)))^2 + 8/(ell*(ell+1))], (* Location of the Schwarzschild potential peak *)
      Vmin -> (1-2*mass/xmin) * (ell*(ell+1)/xmin^2 + 2/xmin^3), (* Value of the Schwarzschild potential at the peak *)
      V -> (Cosh[(x-xmin)/(2*2.75*mass)])^2, (* Eckart potential; 2.75 constant from discussion following equation (6) in Phys.Lett.A 21, 251-254 (1996)*)
      V -> Vmin/V
    }
};

(**************************************************************************************)
(* Wave equation *)
(**************************************************************************************)


waveCalc[fdOrder_, PD_, PDadvect_] := 
{
  Name -> "wave_" <> fdOrder,
  Schedule -> {"in MoL_CalcRHS as wave"},
  ConditionalOnKeyword -> {"fd_order", fdOrder},
  Shorthands -> shorthands,
  Where -> Interior,
  Equations -> 
  {
     (* Shorthands *)

     DDPhi -> PD[Phi,1,1], (* second derivative of Phi wrt x *)

     dot[Phi] -> Psi,
     dot[Psi] -> DDPhi - V Phi
  }
};


(**************************************************************************************)
(* Boundary conditions *)
(**************************************************************************************)

BoundaryParam = 
{
  Name -> "apply_radiative_bcs",
  Default -> "no",
  AllowedValues -> {"yes", "no"}
};

boundaryCalc =
{
  Name -> "boundary",
  Schedule -> {"in MoL_RHSBoundaries"},
  ConditionalOnKeyword -> {"apply_radiative_bcs", "yes"},
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
    
    dot[Phi] -> -(Phi - 0) / r + n[uk] PDonesided2nd[Phi, lk],
    dot[Psi] -> -(Psi - 0) / r + n[uk] PDonesided2nd[Psi, lk] 
  }
};

(**************************************************************************************)
(* Energy *)
(**************************************************************************************)

energyCalc = 
{
  Name -> "calc_energy",
  Schedule -> {"at CCTK_ANALYSIS"},
  Shorthands -> shorthands,
  Where -> Interior,
  Equations ->
  {
     DDPhi -> KroneckerDelta[ui,uj] PDstandard4th[Phi,li,lj], 

     energy -> DDPhi + Psi^2
  }
};

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
  initScalarGFsCalc,
  setIDPulse,
  setPotential,

  waveCalc["2nd", PDstandard2nd, PDonesided2nd],
  waveCalc["4th", PDstandard4th, PDonesided2nd],
  waveCalc["full4th", PDstandard4th, PDlopsided4th],

  boundaryCalc,
  energyCalc
};

keywordParameters = 
{
  fdOrderParam,
  IDParam,
  BoundaryParam
};

CreateKrancThornTT[groups, ".", "Potential", 
  Calculations -> calculations,
  DeclaredGroups -> declaredGroupNames,
  PartialDerivatives -> derivatives,
  RealParameters -> realParameters,
  KeywordParameters -> keywordParameters];
