
Get["KrancThorn`"];

SetEnhancedTimes[False];
SetSourceLanguage["C"];


(**************************************************************************************)
(* Tensors *)
(**************************************************************************************)

(* Register all the tensors that will be used with TensorTools *)
Map[DefineTensor, 
{
  h, hInv, phi, A, K, alpha, Gam, beta, betat, R, Rphi, gamma, bssnmom,  bssnCcons,
  g, k, AInv, DDphi, B, DDalpha, divA, gInv, divh,
  gPhys, kPhys, n, dir, dampG, dampDD, dampDDD, dampDG, ell
}];



Map[AssertSymmetricDecreasing, 
{
  h[la,lb], hInv[ua,ub], A[la,lb], R[la,lb], Rphi[la,lb],
  g[la,lb], k[la,lb]
}];


(**************************************************************************************)
(* Groups *)
(**************************************************************************************)

(* Cactus group definitions *)


admGroups = 
  {{"admbase::metric", {gxx,gxy,gxz,gyy,gyz,gzz}},
   {"admbase::curv", {kxx,kxy,kxz,kyy,kyz,kzz}},
   {"admbase::lapse", {alp}},
   {"admbase::shift", {betax,betay,betaz}}};

declaredGroups = Join[];
declaredGroupNames = Map[First, declaredGroups];

groups = Join[declaredGroups, admGroups];

(**************************************************************************************)
(* Shorthands *)
(**************************************************************************************)

shorthands = 
{
  ell[li], ss, rr
};

(**************************************************************************************)
(* Parameters *)
(**************************************************************************************)

realParameters = 
{
  {Name -> eps, Default -> 1e-3},
  x0,y0,z0, mass,{Name -> pp, Default -> 4}
};

intParameters = 
{

}

k11=kxx; k21=kxy; k22=kyy; k31=kxz; k32=kyz; k33=kzz;
g11=gxx; g21=gxy; g22=gyy; g31=gxz; g32=gyz; g33=gzz;

kerrSchildCalc = 
{
  Name -> "kerr_schild",
  Schedule -> {"in ADMBase_InitialData"},
  Shorthands -> shorthands,
  Equations -> 
  {
    rr -> (x^2 + y^2 + z^2)^(1/2),
    ss -> (rr^pp + eps^pp)^(1/pp),
    ell1 -> x/ss,
    ell2 -> y/ss,
    ell3 -> z/ss,
    g[li,lj] -> Euc[li,lj] + 2 mass / ss ell[li] ell[lj],
    alp -> (1 + 2 mass / ss)^(-1/2),
    betax -> -2 mass / ss ell1 (1+2 mass / ss) ^(-1),
    betay -> -2 mass / ss ell2 (1+2 mass / ss) ^(-1),
    betaz -> -2 mass / ss ell3 (1+2 mass / ss) ^(-1),
    k[li,lj] -> 2 mass / ss^4 (1+2 mass / ss)^(-1/2) (ss^2 Euc[li,lj] - (2+mass/ss) ss^2 ell[li] ell[lj])
  }
};

calculations = 
{
  kerrSchildCalc
};

keywordParameters = 
{
};

CreateKrancThornTT[groups, ".", "KerrSchild", 
  Calculations -> calculations,
  DeclaredGroups -> declaredGroupNames,
(*  PartialDerivatives -> derivatives,*)
  IntParameters -> intParameters,
  RealParameters -> realParameters,
  KeywordParameters -> keywordParameters,
  InheritedImplementations -> {"admbase"}];
