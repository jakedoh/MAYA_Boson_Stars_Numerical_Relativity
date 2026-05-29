(* ::Package:: *)

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
  beta, n, dir, exetaBeta, gInv, g
}];


(* Set the attributes of all tensor-type quantities for which the
   defaults are not suitable. *)

(* Register the TensorTools symmetries (this is very simplistic) *)
Map[AssertSymmetricIncreasing,
{
  gInv[ua,ub]
}];

Map[AssertSymmetricDecreasing, 
{
  g[la,lb]
}];

(* Determinants of the metrics in terms of their components
   (Mathematica symbolic expressions) *)
gDet = Det[MatrixOfComponents[g[la,lb]]];

fnPrefix = "eta_beta"; 
thornName = "ExternalEtaBeta";

(**************************************************************************************)
(* Groups *)
(**************************************************************************************)

(* Cactus group definitions *)
admGroups = 
  {{"admbase::metric", {gxx,gxy,gxz,gyy,gyz,gzz}},
   {"admbase::shift", {betax,betay,betaz}}};

declaredGroups = Map[CreateGroupFromTensor, {exetaBeta}];
declaredGroupNames = Map[First, declaredGroups];

groups = Join[declaredGroups, admGroups];

(**************************************************************************************)
(* Substitutions *)
(**************************************************************************************)

(* Use the AdmBase variable names *)
g11=gxx;  g21=gxy;  g22=gyy;  g31=gxz;  g32=gyz;  g33=gzz;
beta1=betax; beta2=betay; beta3=betaz;

(**************************************************************************************)
(* Shorthands *)
(**************************************************************************************)

shorthands = 
{
  detg, invdetg, gInv[ua,ub],
  n[ui], dir[ui], S,
  Mtot, A, C1, C2,
  xx1, xx2, yy1, yy2, zz1, zz2,
  r1hat, r2hat, top, bottom, rr1, rr2
};

(**************************************************************************************)
(* Parameters *)
(**************************************************************************************)

intParameters =
{
  {Name -> surfindex1, Default -> 0},
  {Name -> surfindex2, Default -> 1}
};

realParameters =
{
  {Name -> w1, 
   Description -> "weight in front of f(r) in denominator of BH+ contribution, w1*f(r1)^n",
   Default -> 0, Steerable -> Always},
  {Name -> w2, 
   Description -> "weight in front of f(r) in denominator of BH- contribution, w2*f(r2)^n",
   Default -> 0, Steerable -> Always},
  {Name -> n, 
   Description -> "function f(r) raised to the power n, wi*f(ri)^n",
   Default -> 0, Steerable -> Always},
  {Name -> Mplus, 
   Description -> "horizon mass of the BH at positive x position",
   Default -> 0.5, Steerable -> Always},
  {Name -> Mminus, 
   Description -> "horizon mass of the BH at negative x position",
   Default -> 0.5, Steerable -> Always},
  {Name -> sigma, 
   Description -> "width of the tanh wall.  Used when rescaling of exetaBeta is done via tanh",
   Default -> 1, Steerable -> Always},
  {Name -> Rcutoff, 
   Description -> "Distance from origin that the rescaled exetaBeta= 1/2 exetaBeta_BG.  Only used when rescaling",
   Default -> 1, Steerable -> Always}
};

(* Keyword Params *)

typeParam =
{
  Name -> "etaBeta_form", 
  Description -> "Form of the function in the denominator of exetaBeta",
  Default -> "polynomial", AllowedValues -> {"polynomial", "exponential"} (*, "gridfunction"}*)
};

rescaleParam =
{
  Name -> "rescale_etaBeta", 
  Description -> "Resale etaBeta to zero as r->infinity using either a polynomial form or hyperbolic tangent",
  Default -> "no", AllowedValues -> {"tanh", "yes", "no"}
};

(*fdOrderParam =
{
  Name -> "fd_order",
  Description -> "Finite differencing scheme order",
  Default -> "2nd",
  AllowedValues -> {"2nd", "2ndcentred", "full4th", "4thcentred", "6thcentred", "6th"}
};*)

(*boundaryParam =
{
  Name -> "boundary_condition",
  Default -> "radiative",
  AllowedValues -> {"none", "radiative"}
};*)

(*evolveParam =
{
  Name -> "evolve_etaBeta",
  Default -> "no",
  AllowedValues -> {"no", "yes"}
};*)

(*initEtaBetaRHS =
{
  Name -> "init_etaBeta_rhs",
  Default -> "no",
  AllowedValues -> {"no", "yes"}
};*)

keywordParameters =      
{
  rescaleParam,
  typeParam
  (*fdOrderParam,*)
  (*
  boundaryParam,
  evolveParam,
  initEtaBetaRHS*)
};

(**************************************************************************************)
(* Initialize all grid functions to a recognizable value. These should be overwritten *)
(* later. *)
(**************************************************************************************)

INITVALUE = "33";

initGFsCalc = 
{
  Name -> fnPrefix <> "_init_etaBeta",
  Schedule -> {"at BaseGrid as etaBeta_initgfs"},
  Where -> Everywhere,
  Equations -> 
  {
    exetaBeta -> INITVALUE
  }
};

(*initGFsRHSCalc =
{
  Name -> fnPrefix <> "_init_etaBeta_rhs",
  Schedule -> {"in MoL_CalcRHS as etaBeta_initgfs before etaBeta_calc_etaBeta"},
  Where -> Everywhere,
  ConditionalOnKeywords -> {{"evolve_etaBeta","no"},{"init_etaBeta_rhs","yes"}},
  Equations ->
  {
    exetaBeta -> INITVALUE
  }
};*)

(**************************)
(* Evolution Calculations *)
(**************************)

(* Initial Calculation for evolved... uses exetaBeta = 1.0 *)
(*InitCalc[fdOrder_, PD_, PDadvect_] :=
{
  Name -> fnPrefix <> "_init_calc_" <> fdOrder,
  Schedule -> {"at INITIAL after ADMBase_PostInitial as etaBeta_init_calc"},
  ConditionalOnKeywords -> {{"fd_order", fdOrder},{"evolve_etaBeta", "yes"}},
  Where -> Interior,
  Shorthands -> shorthands,
  Equations ->
    {
     (* placeholder for xx1 thru zz2 *)
     xx1 -> 0, yy1 -> 0, zz1 -> 0,
     xx2 -> 0, yy2 -> 0, zz2 -> 0,

     rr1 -> Sqrt[ (xx1 - x)^2 + (yy1 - y)^2 + (zz1 - z)^2 ],
     rr2 -> Sqrt[ (xx2 - x)^2 + (yy2 - y)^2 + (zz2 - z)^2 ],

     detg -> gDet,
     invdetg -> 1.0 / detg,
     gInv[ua,ub] -> invdetg gDet MatrixInverse[g[ua,ub]],

     top -> gInv[ui,uj] gInv[ua,ub] gInv[uc,ud] PD[g[la,lb],li] PD[g[lc,ld],lj],
     bottom -> 6 (1 - detg^(-1.0/6.0))^2,

     S -> ( 1.31241 Sqrt[Abs[top]] / bottom ) ( Rcutoff^2 / ( r^2 + Rcutoff^2 ) ),
     exetaBeta -> 2.0
    }
};

(* etaBeta evolution *)
evolveEtaBetaCalc[fdOrder_, PD_, PDadvect_] := 
{
  Name -> fnPrefix <> "_evolve_etaBeta_" <> fdOrder,
  Schedule -> {"in MoL_CalcRHS as etaBeta_evolve_etaBeta"},
  Where -> Interior,
  ConditionalOnKeywords -> {{"fd_order", fdOrder},{"evolve_etaBeta", "yes"}},
  Shorthands -> shorthands,
  Equations -> 
  {
    (* placeholder for xx1 thru zz2 *)
     xx1 -> 0, yy1 -> 0, zz1 -> 0,
     xx2 -> 0, yy2 -> 0, zz2 -> 0,

     rr1 -> Sqrt[ (xx1 - x)^2 + (yy1 - y)^2 + (zz1 - z)^2 ],
     rr2 -> Sqrt[ (xx2 - x)^2 + (yy2 - y)^2 + (zz2 - z)^2 ],

     detg -> gDet,
     invdetg -> 1.0 / detg,
     gInv[ua,ub] -> invdetg gDet MatrixInverse[g[ua,ub]],

     top -> gInv[ui,uj] gInv[ua,ub] gInv[uc,ud] PD[g[la,lb],li] PD[g[lc,ld],lj],
     bottom -> 6 (1 - detg^(-1.0/6.0))^2,

     S -> ( 1.31241 Sqrt[Abs[top]]/bottom ) ( Rcutoff^2 / ( r^2 + Rcutoff^2 ) ),

    (* Shift advection *)
    dir[ui] -> Sign[beta[ui]],

    (* Mtot = 1 *)
    Mtot -> 1.0,

    dot[exetaBeta] -> beta[ui] PDadvect[exetaBeta,li] + 1/Mtot ( -exetaBeta + S )
                   
  }
};*)

(**********************************************)
(* Grid function form of etaBeta Calculations *)
(**********************************************)

(* Initial calculation of grid function Mueller term from arxiv: *)
(*EtaBetaGFInitCalc[fdOrder_, PD_, PDadvect_] :=
{
  Name -> fnPrefix <> "_calc_etaBeta_gf_init_" <> fdOrder,
  Schedule -> {"at INITIAL after ADMBase_PostInitial as etaBeta_calc_etaBeta"},
  ConditionalOnKeywords -> {{"fd_order", fdOrder},{"etaBeta_form", "gridfunction"}},
  Where -> Interior,
  Shorthands -> shorthands,
  Equations ->
    {
     detg -> gDet,
     invdetg -> 1.0 / detg,
     gInv[ua,ub] -> invdetg gDet MatrixInverse[g[ua,ub]],

     top -> gInv[ui,uj] gInv[ua,ub] gInv[uc,ud] PD[g[la,lb],li] PD[g[lc,ld],lj],
     bottom -> 6 (1 - detg^(-a/6.0))^b,

     exetaBeta -> 1.31241 Sqrt[Abs[top]]/bottom
   }
};

(* Same calculation scheduled in MoL_CalcRHS before the shift evolution *)
EtaBetaGFCalc[fdOrder_, PD_, PDadvect_] :=
{
  Name -> fnPrefix <> "_calc_wetaBeta_gf_" <> fdOrder,
  Schedule -> {"in MoL_CalcRHS before bssn_shift as etaBeta_calc_etaBeta"},
  ConditionalOnKeywords -> {{"fd_order", fdOrder},{"etaBeta_form", "gridfunction"}},
  Where -> Interior,
  Shorthands -> shorthands,
  Equations ->
   {
     detg -> gDet,
     invdetg -> 1.0 / detg,
     gInv[ua,ub] -> invdetg gDet MatrixInverse[g[ua,ub]],

     top -> gInv[ui,uj] gInv[ua,ub] gInv[uc,ud] PD[g[la,lb],li] PD[g[lc,ld],lj],
     bottom -> 6 (1 - detg^(-a/6.0))^b,

     exetaBeta -> 1.31241 Sqrt[Abs[top]]/bottom
   }
};*)

(***********************************)
(* Calculate etaBeta algebraically *)
(***********************************)

(* polynomial form of etaBeta from arxiv: 1003.4681*)
EtaBetaPolyCalc :=
{
  Name -> fnPrefix <> "_calc_etaBeta_poly",
  Schedule -> {"at INITIAL after ShiftTracker_Set_SphericalSurface as etaBeta_calc_etaBeta","in MoL_CalcRHS before bssn_shift as etaBeta_calc_etaBeta"},
  Where -> Everywhere,
  ConditionalOnKeyword -> {"etaBeta_form", "polynomial"},
  Shorthands -> shorthands,
  Equations ->
  {
     (* calculate consts *)
     Mtot -> Mplus + Mminus,
     A  -> 2.0/Mtot,
     C1 -> 1.0/Mplus - A,
     C2 -> 1.0/Mminus - A,

     (* placeholder for xx1 thru zz2 *)
     xx1 -> 0, yy1 -> 0, zz1 -> 0,
     xx2 -> 0, yy2 -> 0, zz2 -> 0,

     (* calculate rhati = |rveci - rvec|/|rvec1 - rvec2| *)
     r1hat -> Sqrt[ (xx1-x)^2+(yy1-y)^2+(zz1-z)^2 ] / Sqrt[ (xx1-xx2)^2+(yy1-yy2)^2+(zz1-zz2)^2 ],
     r2hat -> Sqrt[ (xx2-x)^2+(yy2-y)^2+(zz2-z)^2 ] / Sqrt[ (xx1-xx2)^2+(yy1-yy2)^2+(zz1-zz2)^2 ],

     (* calculate etaBeta *)
     exetaBeta -> A + C1/( 1 + w1 (r1hat^2)^n ) + C2/( 1 + w2 (r2hat^2)^n )
  }
};

(* initial polynomial calculation *)
(*EtaBetaPolyInitCalc :=
{
  Name -> fnPrefix <> "_calc_etaBeta_poly_init",
  Schedule -> {"at INITIAL after ShiftTracker_Set_SphericalSurface as etaBeta_calc_etaBeta"},
  Where -> Everywhere,
  ConditionalOnKeyword -> {"etaBeta_form", "polynomial"},
  Shorthands -> shorthands,
  Equations ->
  {
     (* calculate consts *)
     Mtot -> Mplus + Mminus,
     A  -> 2.0/Mtot,
     C1 -> 1.0/Mplus - A,
     C2 -> 1.0/Mminus - A,

     (* placeholder for xx1 thru zz2 *)
     xx1 -> 0, yy1 -> 0, zz1 -> 0,
     xx2 -> 0, yy2 -> 0, zz2 -> 0,

     (* calculate rhati = |rveci - rvec|/|rvec1 - rvec2| *)
     r1hat -> Sqrt[ (xx1-x)^2+(yy1-y)^2+(zz1-z)^2 ] / Sqrt[ (xx1-xx2)^2+(yy1-yy2)^2+(zz1-zz2)^2 ],
     r2hat -> Sqrt[ (xx2-x)^2+(yy2-y)^2+(zz2-z)^2 ] / Sqrt[ (xx1-xx2)^2+(yy1-yy2)^2+(zz1-zz2)^2 ],

     (* calculate etaBeta *)
     exetaBeta -> A + C1/( 1 + w1 (r1hat^2)^n ) + C2/( 1 + w2 (r2hat^2)^n )
  }
};*)

(* exponential form of etaBeta *)
EtaBetaExpCalc :=
{
  Name -> fnPrefix <> "_calc_etaBeta_exp",
  Schedule -> {"in MoL_CalcRHS before bssn_shift as etaBeta_calc_etaBeta","at INITIAL after ShiftTracker_Set_SphericalSurface as etaBeta_calc_etaBeta"},
  Where -> Everywhere,
  ConditionalOnKeyword -> {"etaBeta_form", "exponential"},
  Shorthands -> shorthands,
  Equations ->
  {
     (* calculate consts *)
     Mtot -> Mplus + Mminus,
     A 	-> 2.0/Mtot,
     C1 -> 1.0/Mplus - A,
     C2 -> 1.0/Mminus - A,

     (* placeholder for xx1 thru zz2 *)
     xx1 -> 0, yy1 -> 0, zz1 -> 0,
     xx2 -> 0, yy2 -> 0, zz2 -> 0,

     (* calculate rhati = |rveci - rvec|/|rvec1 - rvec2| *)
     r1hat -> Sqrt[ (xx1-x)^2+(yy1-y)^2+(zz1-z)^2 ] / Sqrt[ (xx1-xx2)^2+(yy1-yy2)^2+(zz1-zz2)^2 ],
     r2hat -> Sqrt[ (xx2-x)^2+(yy2-y)^2+(zz2-z)^2 ] / Sqrt[ (xx1-xx2)^2+(yy1-yy2)^2+(zz1-zz2)^2 ],

     (* calculate etaBeta *)
     exetaBeta -> A + C1 Exp[-w1 (r1hat^2)^n] + C2 Exp[-w2 (r2hat^2)^n]
  }
};

(* initial calculation of etaBeta as an exponential *)
(*EtaBetaExpInitCalc :=
{
  Name -> fnPrefix <> "_calc_etaBeta_exp_init",
  Schedule -> {"at INITIAL after ShiftTracker_Set_SphericalSurface as etaBeta_calc_etaBeta"},
  Where -> Everywhere,
  ConditionalOnKeyword -> {"etaBeta_form", "exponential"},
  Shorthands -> shorthands,
  Equations ->
  {
     (* calculate consts *)
     Mtot -> Mplus + Mminus,
     A  -> 2.0/Mtot,
     C1 -> 1.0/Mplus - A,
     C2 -> 1.0/Mminus - A,

     (* placeholder for xx1 thru zz2 *)
     xx1 -> 0, yy1 -> 0, zz1 -> 0,
     xx2 -> 0, yy2 -> 0, zz2 -> 0,

     (* calculate rhati = |rveci - rvec|/|rvec1 - rvec2| *)
     r1hat -> Sqrt[ (xx1-x)^2+(yy1-y)^2+(zz1-z)^2 ] / Sqrt[ (xx1-xx2)^2+(yy1-yy2)^2+(zz1-zz2)^2 ],
     r2hat -> Sqrt[ (xx2-x)^2+(yy2-y)^2+(zz2-z)^2 ] / Sqrt[ (xx1-xx2)^2+(yy1-yy2)^2+(zz1-zz2)^2 ],

     (* calculate etaBeta *)
     exetaBeta -> A + C1 Exp[-w1 (r1hat^2)^n] + C2 Exp[-w2 (r2hat^2)^n]
  }
};*)

(**************************************************)
(* rescale etaBeta to fall off to 0 at boundaries *)
(**************************************************)

(*RescaleEtaBetaInitCalc :=
{
  Name -> fnPrefix <> "_rescale_etaBeta_init_calc",
  Schedule -> {"at INITIAL after etaBeta_calc_etaBeta as rescale_etaBeta"},
  Where -> Everywhere,
  ConditionalOnKeyword -> {"rescale_etaBeta", "yes"},
  Shorthands -> shorthands,
  Equations ->
  {
     (* calculate etaBeta *)
     exetaBeta -> exetaBeta Rcutoff^2 / ( r^2 + Rcutoff^2 ) 
  }
};*)

RescaleEtaBetaCalc :=
{
  Name -> fnPrefix <> "_rescale_etaBeta_calc",
  Schedule -> {"at INITIAL after etaBeta_calc_etaBeta as rescale_etaBeta","in MoL_CalcRHS after etaBeta_calc_etaBeta as rescale_etaBeta"},
  Where -> Everywhere,
  ConditionalOnKeyword -> {"rescale_etaBeta", "yes"},
  Shorthands -> shorthands,
  Equations ->
  {
     (* calculate etaBeta *)
     exetaBeta -> exetaBeta Rcutoff^2 / ( r^2 + Rcutoff^2 ) 
  }
};

RescaleEtaBetaTanhCalc :=
{
  Name -> fnPrefix <> "_rescale_etaBeta_tanh_calc",
  Schedule -> {"at INITIAL after etaBeta_calc_etaBeta as rescale_tanh_etaBeta","in MoL_CalcRHS after etaBeta_calc_etaBeta as rescale_tanh_etaBeta"},
  Where -> Everywhere,
  ConditionalOnKeyword -> {"rescale_etaBeta", "tanh"},
  Shorthands -> shorthands,
  Equations ->
  {
     (* calculate etaBeta *)
     exetaBeta -> exetaBeta * ( -0.5 Tanh[ ( r - Rcutoff ) / sigma ] + 0.5 )
  }
};


(***********************)
(* Boundary conditions *)
(***********************)

(*boundaryCalc =
{
  Name -> fnPrefix <> "_boundary",
  Schedule -> {"in MoL_RHSBoundaries as etaBeta_boundary"},
  ConditionalOnKeywords -> {{"boundary_condition", "radiative"},{"evolve_etaBeta","yes"}},
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

    (* exetaBeta -> 0 as r -> inf *)
    dot[exetaBeta] -> -exetaBeta / r + n[uk] PDonesided2nd[exetaBeta, lk]

  }
};*)

(**************************************************************************************)
(* Construct the thorn *)
(**************************************************************************************)

calculations =
{
  (*initGFsCalc,
  initGFsRHSCalc,*)
  (*InitCalc["2nd", PDstandard2nd, PDonesided2nd],
  InitCalc["2ndcentred", PDstandard2nd, PDstandard2nd],
  InitCalc["4thcentred", PDstandard4th, PDstandard4th],
  InitCalc["6thcentred", PDstandard6th, PDstandard6th],
  EtaBetaGFCalc["2nd", PDstandard2nd, PDonesided2nd],
  EtaBetaGFCalc["2ndcentred", PDstandard2nd, PDstandard2nd],
  EtaBetaGFCalc["4thcentred", PDstandard4th, PDstandard4th],
  EtaBetaGFCalc["6thcentred", PDstandard6th, PDstandard6th],
  EtaBetaGFInitCalc["2nd", PDstandard2nd, PDonesided2nd],
  EtaBetaGFInitCalc["2ndcentred", PDstandard2nd, PDstandard2nd],
  EtaBetaGFInitCalc["4thcentred", PDstandard4th, PDstandard4th],
  EtaBetaGFInitCalc["6thcentred", PDstandard6th, PDstandard6th],
  evolveEtaBetaCalc["2nd", PDstandard2nd, PDonesided2nd],
  evolveEtaBetaCalc["2ndcentred", PDstandard2nd, PDstandard2nd],
  evolveEtaBetaCalc["4thcentred", PDstandard4th, PDstandard4th],
  evolveEtaBetaCalc["6thcentred", PDstandard6th, PDstandard6th],*)
  (*InitCalc["full4th", PDstandard4th, PDlopsided4th],
  InitCalc["6th", PDstandard6th, PDlopsided6th],*)
  (*EtaBetaGFCalc["full4th", PDstandard4th, PDlopsided4th],
  EtaBetaGFCalc["6th", PDstandard6th, PDlopsided6th],
  EtaBetaGFInitCalc["full4th", PDstandard4th, PDlopsided4th],
  EtaBetaGFInitCalc["6th", PDstandard6th, PDlopsided6th],*)
  (*evolveEtaBetaCalc["full4th", PDstandard4th, PDlopsided4th],
  evolveEtaBetaCalc["6th", PDstandard6th, PDlopsided6th],*)
  (*EtaBetaPolyInitCalc,*)
  (*EtaBetaExpInitCalc,*)
  (*boundaryCalc*)
  RescaleEtaBetaCalc,
  RescaleEtaBetaTanhCalc,
  EtaBetaPolyCalc,
  EtaBetaExpCalc
};

 
CreateKrancThornTT[groups, ".", thornName,
  Calculations -> calculations,
  DeclaredGroups -> declaredGroupNames,
  PartialDerivatives -> derivatives,
  IntParameters -> intParameters,
  RealParameters -> realParameters,
  KeywordParameters -> keywordParameters,
  UseLoopControl -> True,
  InheritedImplementations -> {"admbase","sphericalsurface"}];

