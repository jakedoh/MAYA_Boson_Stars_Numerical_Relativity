(* ::Package:: *)

Get["KrancThorn`"];

SetEnhancedTimes[False];
SetSourceLanguage["C"];

(**************************************************************************************)
(* variable types *)
(*this notebook is used to generate the evolution of multi scalar field based on prd97,023526(2018)*)
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
  n, dir, beta, 
  g, gInv, gLocal, gLocalInv, 
  kcurv, kLocal, KtraceLocal, ktrace,
  Gam3, Gam0ij, Gamkij, (*define the variables of metric and Connections*)
  Tmn, rhoE, phiamp, TmnLocal, ttrace,
  pia, phia, pib, phib, SFSi, dtphia, dtphib, (*vpotential,*)(*define the variables of scalar fields*)
  (*ginvphiaa, ginvphiab, ginvphibb, gvphiaa, gvphiab, gvphibb, Gamaaa, Gamaab, Gamabb, (*define the curverd target space of the scalar fields*)*)
  SFTtt, SFTtj, SFTij
}];

Map[AssertSymmetricIncreasing,
{
  gInv[ua,ub], gLocalInv[ua,ub]
}];

Map[AssertSymmetricDecreasing, 
{
  g[la,lb], gLocal[la,lb], kcurv[la,lb], kLocal[la,lb], SFTij[la,lb], Tmn[la,lb], TmnLocal[la,lb], Gam0ij[la,lb]
}];

AssertSymmetricDecreasing[Gam3[ua,lb,lc], lb, lc];
AssertSymmetricDecreasing[Gamkij[ua,lb,lc], lb, lc];

(* Determinants of the metrics in terms of their components
   (Mathematica symbolic expressions) *)
gDet = Det[MatrixOfComponents[gLocal[la,lb]]];


fnPrefix = "EvolveScalarFields";
thornName = "EvolveScalarFields";

(**************************************************************************************)
(* Groups *)
(**************************************************************************************)

(* Cactus group definitions *)

evolvedGroups = Map[CreateGroupFromTensor, {phia, pia, phib, pib(*, vpotential*)}];

admGroups = 
  {{"admbase::metric", {gxx,gxy,gxz,gyy,gyz,gzz}},
   {"admbase::curv",  {kxx,kxy,kxz,kyy,kyz,kzz}},
   {"admbase::lapse", {alp}},
   {"admbase::shift", {betax,betay,betaz}}};

evaluatedGroups = Map[CreateGroupFromTensor,{gLocal[li,lj], ttrace, SFTtt, rhoE, phiamp}];
evaluatedGroups = Map[AddGroupExtra[#, Timelevels -> 3] &, evaluatedGroups];

stressenergyGroups =
  {{"TmunuBase::stress_energy_scalar", {eTtt} },
   {"TmunuBase::stress_energy_vector", {eTtx,eTty,eTtz} },
   {"TmunuBase::stress_energy_tensor", {eTxx,eTxy,eTxz,eTyy,eTyz,eTzz} } };

declaredGroups = Join[evolvedGroups, evaluatedGroups];
declaredGroupNames = Map[First, declaredGroups];

groups = Join[declaredGroups, admGroups, stressenergyGroups];

(**************************************************************************************)
(* Shorthands *) (*Temporary variables which will be used in this calculation*)
(**************************************************************************************)

shorthands = 
{
  ktrace, 
  n[ui], dir[ui], detg, invdetg, 
  gLocalInv[ui,uj], gInv[ui,uj],
  (*TttLocal, TtxLocal, TtyLocal, TtzLocal, TmnLocal[li,lj],*)
  betaSqr, dphiaSqr, dphibSqr, dotpia, dotpib,
  (*factornew, ginvphiaa, ginvphiab, ginvphibb, gvphiaa, gvphiab, gvphibb,*) (*define the metric fo curved field space*)
  (*Gamaaa, Gamaab, Gamabb, Gambaa, Gambab, Gambbb,*) (*define the Christoff symbols fo curved field space*)
   SFTtj[li], SFTij[li,lj],  SFSi[li], 
  Gam3[ua,lb,lc], Gam0ij[li,lj], Gamkij[uk,li,lj], 
  kLocal[li,lj], Ktrace, dtphia, derivVa, dtphib, derivVb, STpotential, phiaRS, piaRS, phibRS, pibRS
};

(* Use the names of shift vector from ADMBase *)
alpha=alp;
beta1=betax; beta2=betay; beta3=betaz;
kcurv11=kxx; kcurv21=kxy; kcurv22=kyy; kcurv31=kxz; kcurv32=kyz; kcurv33=kzz; 
g11=gxx; g21=gxy; g22=gyy; g31=gxz; g32=gyz; g33=gzz;

(* Use the names from TmunuBase*)
Tmn11=eTxx; Tmn21=eTxy; Tmn22=eTyy; Tmn31=eTxz; Tmn32=eTyz; Tmn33=eTzz; 

(**************************************************************************************)
(* Parameters *)
(**************************************************************************************)

realParameters = 
{
  {Name -> phiamass, Default -> 0 },
  {Name -> phibmass, Default -> 0 },
  {Name -> vv0, Default -> 0 },
  {Name -> vv1, Default -> 0 },
  {Name -> phivsign, Default -> 1 },
  {Name -> r0apotential, Default -> 0 },
  {Name -> controlmasspositiveornegative, Default -> 1 },
  {Name -> phiaBG, Default -> 0},
  {Name -> piaBG, Default -> 0},
  {Name -> phibBG, Default -> 0},
  {Name -> pibBG, Default -> 0},
  {Name -> gcouple, Default -> 0},
  {Name -> lambdaphia, Default -> 0},
  {Name -> lambdaphib, Default -> 0},
  {Name -> lambdaphian, Default -> 0},
  {Name -> lambdaphibn, Default -> 0},
  {Name -> phiav, Default -> 0},
  {Name -> phibv, Default -> 0},
  {Name -> Vi0, Default -> 1},
  {Name -> phi00, Default -> 1}
};

(* If this is set to 'no', then a function that puts zeros in  TmunuTraceMatterLocal is scheduled. *)
computeScalarFieldTmunuParam =
{
  Name -> "couple_scalar_field_to_gravity",
  Default -> "yes",
  AllowedValues -> {"yes", "no"}
};

fdOrderParam =
{
  Name -> "fd_order",
  Default -> "full4th",
  AllowedValues -> {"2nd", "full4th", "4thcentred", "6thcentred", "6th"}
};

(**************************************************************************************)
(* Initialize all grid functions to a recognizable value. These should be overwritten *)
(* later. *)
(**************************************************************************************)

INITVALUE = "100";

initGFsCalc = 
{
  Name -> fnPrefix <> "_init_gfs",
  Schedule -> {"at BaseGrid"},
  Equations -> 
  {
    phia      -> 0,
    pia       -> 0,
    phib      -> 0,
    pib       -> 0,
    rhoE       -> 0,
    phiamp     -> 0
    (*vpotential-> 0*)
  }
};

(**************************************************************************************)
(* scalar function eq. (19) of prd97,023526(2018)*)
(**************************************************************************************)
phivxyz=vv0*(vv1+phivsign*Tanh[r-r0apotential]);
(**************************************************************************************)
(* potential functions *)
(**************************************************************************************)
Vnone[phia_, phib_]:= 0;
Vaxion[phia_, phib_]:= lambdaphian/8 (phia^2-phivxyz^2)^2 + lambdaphibn/8 (phib^2-phivxyz^2)^2;
Vonlymass[phia_, phib_]:= 1/2*controlmasspositiveornegative*phiamass^2*phia^2 + 1/2*controlmasspositiveornegative*phibmass^2*phib^2;
Vdoublewell[phia_, phib_]:= 1/2 phiamass^2 phia^2 + 1/2 phibmass^2 phib^2 +  gcouple/2 phia^2 phib^2  + lambdaphian/4 (phia^2-phiav^2)^2 + lambdaphibn/4 (phib^2-phibv^2)^2; 
VphiFourthall[phia_, phib_]:= 1/2 phiamass^2 phia^2 + 1/2 phibmass^2 phib^2+ lambdaphia/4 phia^4 + lambdaphib/4 phib^4 + gcouple/2 phia^2 phib^2 + lambdaphian/4 (phia^2-phiav^2)^2  + lambdaphibn/4 (phib^2-phibv^2)^2; 
Vphiinflation[phia_, phib_]:= Vi0 + lambdaphia/4 (phia^2-phi00^2)^2  + lambdaphibn/4 (phib^2-phibv^2)^2; 
Vphilinear[phia_, phib_] := Vi0 + lambdaphia (phia) + lambdaphibn (phib); 

Vprimea[V_,x_,y_]:= D[V[x, y], x];
Vprimeb[V_,x_,y_]:= D[V[x, y], y];
potentialParam =
{
  Name -> "potential_type",
  Default -> "doublewell",
  AllowedValues -> {"none","onlymass","axion","doublewell","phiFourth", "inflation", "linear"}
};



(**************************************************************************************)
(* Evolution of the scalar field*)
(**************************************************************************************)

evolveScalarField[fdOrder_, PD_, PDadvect_, potential_, V_] := 
{
  Name -> fnPrefix<>"_EvolveScalarField_"<>fdOrder<>"_potential_"<>potential,
  Schedule -> {"in MoL_CalcRHS as EvolveScalarField"},
  Where -> Interior,
  ConditionalOnKeywords -> {{"fd_order", fdOrder},{"potential_type", potential}},
  Shorthands -> shorthands,
  (*  CollectList -> {gLocalInv[ui,uj]}, *)
  Equations -> 
  {
    dir[ui] -> Sign[beta[ui]],
    (*pi -> 3.14159265358979323846,*)
 
    gInv[ui,uj] -> MatrixInverse[g[ui,uj]],
    ktrace -> kcurv[li,lj] gInv[ui, uj],

    Gam3[ua, lb, lc] -> 1/2 gInv[ua,ud] (PD[g[lb,ld], lc] + PD[g[lc,ld], lb] - PD[g[lb,lc],ld]),
      
    derivVa -> Vprimea[V,phia,phib],
    derivVb -> Vprimeb[V,phia,phib],
    (*vpotential\[Rule]  V[phia,phib],*)
    (* scalar fielda evolution eq. (8) of prd97,023526(2018)*)    
    dot[phia] -> - alpha pia + beta[ui] PDadvect[phia,li], (* partial time-derivative of STphi; STpi is the Lie derivative of STphi along n*)
    dot[pia]  -> beta[ui] PDadvect[pia,li] + alpha ktrace pia
      - alpha gInv[ui,uj] PD[phia,li,lj] + alpha gInv[ui,uj] Gam3[uk,li,lj] PD[phia,lk] - gInv[ui,uj] PD[alpha,lj] PD[phia,li]
      + alpha derivVa,     
    
    (* scalar fielda evolution eq. (8) of prd97,023526(2018)*)    
    dot[phib] -> - alpha pib + beta[ui] PDadvect[phib,li], (* partial time-derivative of STphi; STpi is the Lie derivative of STphi along n*)
    dot[pib]  -> beta[ui] PDadvect[pib,li] + alpha ktrace pib
      - alpha gInv[ui,uj] PD[phib,li,lj] + alpha gInv[ui,uj] Gam3[uk,li,lj] PD[phib,lk] - gInv[ui,uj] PD[alpha,lj] PD[phib,li]
      + alpha derivVb 

  }
};

(*AddPotentialToRHS[potential_,V_] :=
{
  Name -> fnPrefix<>"_AddPotentialToRHS_"<>potential,
  Schedule -> {"in MoL_CalcRHS after EvolveScalarField as AddPotentialToScalarFieldRHS"},
  Where -> Interior,
  ConditionalOnKeyword -> {"potential_type", potential},
  Shorthands -> shorthands,
  Equations  -> 
  {
    factornew -> 2 STA[phia,phib] + 6 STcouplea ^2 phia^2 + 6 STcoupleb ^2 phib^2,
    ginvphiaa -> (2 STA[phia,phib]) / STplankmass^2  (2 STA[phia,phib]+6 STcoupleb ^2 phib^2)/factornew,
    ginvphiab -> -(2 STA[phia,phib]) /STplankmass^2 (6 STcouplea STcoupleb phia phib)/(factornew),
    ginvphibb -> (2 STA[phia,phib]) / STplankmass^2  (2 STA[phia,phib]+6 STcouplea ^2 phia^2)/factornew,
    derivVa -> Vprimea[V,phia,phib],
    derivVb -> Vprimeb[V,phia,phib],
    dot[pia]  -> dot[pia] - alpha ginvphiaa derivVa - alpha ginvphiab derivVb,
    dot[pib]  -> dot[pib] - alpha ginvphiab derivVa - alpha ginvphibb derivVb
  }
};*)

ScalarFieldTmunuCalc[fdOrder_,PD_, PDadvect_, potential_, V_] := 
{
  Name -> fnPrefix <> "_ScalarFieldTmunu_" <> potential <> "_" <> fdOrder,
  (*Schedule -> Automatic,*)
  (*After -> "ADMBase_SetADMVars",*)
  Schedule -> {"in AddToTmunu as ScalarFieldTmunu after Whisky_SetTmunu"},
  Where -> Interior,
  ConditionalOnKeywords -> {{"fd_order", fdOrder},{"couple_scalar_field_to_gravity","yes"},{"potential_type",potential}},
  Shorthands -> shorthands,
  Equations -> 
  {
     dir[ui] -> Sign[beta[ui]],
     (*pi -> 3.14159265358979323846,*)
     
     detg -> gDet,
     invdetg -> 1 / detg,
     gInv[ui,uj] -> MatrixInverse[g[ui,uj]],
     
     (* time derivative from eq of motion *)
     dtphia -> -alpha pia + beta[ui] PDadvect[phia,li],
     dtphib -> -alpha pib + beta[ui] PDadvect[phib,li],
     (* auxiliary variables *)
     betaSqr -> g[li,lj] beta[ui] beta[uj],
     dphiaSqr ->  - pia^2 + gInv[ui,uj] PD[phia,li] PD[phia,lj], 
     dphibSqr ->  - pib^2 + gInv[ui,uj] PD[phib,li] PD[phib,lj], 
    (* eq. (A3) of prd97,023526(2018)*)                          

     (* Components due to the scalar field *)
     (* eq. (7) of prd97,023526(2018)*) 
     SFTtt -> dtphia^2 + dtphib^2  - (-alpha^2 + betaSqr)( 1/2 ( dphiaSqr + dphibSqr) +  V[phia,phib] ),
     SFTtj[li] -> dtphia PD[phia,li] + dtphib PD[phib,li] - g[li,lj] beta[uj]( 1/2 ( dphiaSqr + dphibSqr) +  V[phia,phib] ),
     SFTij[li,lj] -> PD[phia,li] PD[phia,lj] + PD[phib,li] PD[phib,lj] - g[li,lj] (1/2 ( dphiaSqr + dphibSqr)+  V[phia,phib] ),
     (*======================================================================================*)
     (* Trace of the stress energy tensor of matter. At this point assuming this is the last *)
     (* function that adds something to Tmunu, all the matter contributions are in place.    *)
     (* Therefore the following trace is only from the matter fields.                        *)
     (* We compute it here and store it to use it later in the evolution equations           *)
     (*======================================================================================*)

     
     eTtt -> eTtt + SFTtt  ,
     eTtx -> eTtx + SFTtj1 ,
     eTty -> eTty + SFTtj2 ,
     eTtz -> eTtz + SFTtj3 ,
     eTxx -> eTxx + SFTij11,
     eTxy -> eTxy + SFTij21,
     eTxz -> eTxz + SFTij31,
     eTyy -> eTyy + SFTij22,
     eTyz -> eTyz + SFTij32,
     eTzz -> eTzz + SFTij33
     
  }
};

ScalarFieldrhoCalc[fdOrder_,PD_, PDadvect_, potential_, V_] := 
{
  Name -> fnPrefix <> "_ScalarFieldrho_" <> potential <> "_" <> fdOrder,
  Schedule -> Automatic,
  After -> "ADMBase_SetADMVars",
  (*Schedule -> {"in AddToTmunu as ScalarFieldTmunu after Whisky_SetTmunu"},*)
  Where -> Interior,
  ConditionalOnKeywords -> {{"fd_order", fdOrder},{"couple_scalar_field_to_gravity","yes"},{"potential_type",potential}},
  Shorthands -> shorthands,
  Equations -> 
  {
     dir[ui] -> Sign[beta[ui]],
     (*pi -> 3.14159265358979323846,*)
     detg -> gDet,
     invdetg -> 1 / detg,
     gInv[ui,uj] -> MatrixInverse[g[ui,uj]],
     (* time derivative from eq of motion *)
     dtphia -> -alpha pia + beta[ui] PDadvect[phia,li],
     dtphib -> -alpha pib + beta[ui] PDadvect[phib,li],
     (* auxiliary variables *)
     betaSqr -> g[li,lj] beta[ui] beta[uj],
     dphiaSqr ->  - pia^2 + gInv[ui,uj] PD[phia,li] PD[phia,lj], 
     dphibSqr ->  - pib^2 + gInv[ui,uj] PD[phib,li] PD[phib,lj], 

     SFTtt -> dtphia^2 + dtphib^2  - (-alpha^2 + betaSqr)( 1/2 ( dphiaSqr + dphibSqr) +  V[phia,phib] ),
     SFTtj[li] -> dtphia PD[phia,li] + dtphib PD[phib,li] - g[li,lj] beta[uj]( 1/2 ( dphiaSqr + dphibSqr) +  V[phia,phib] ),
     SFTij[li,lj] -> PD[phia,li] PD[phia,lj] + PD[phib,li] PD[phib,lj] - g[li,lj] (1/2 ( dphiaSqr + dphibSqr)+  V[phia,phib] ),

     rhoE -> (1/alpha)^2 ( SFTtt - 2 beta[ui] SFTtj[li] + beta[ui] beta[uj] SFTij[li,lj] )
     
  }
};

ScalarFieldAmplitudeCalc[fdOrder_,PD_, PDadvect_, potential_, V_] :=
{
  Name -> fnPrefix <> "_ScalarFieldAmplitude_" <> potential <> "_" <> fdOrder,
  Schedule -> Automatic,
  After -> "ADMBase_SetADMVars",
  (*Schedule -> {"in AddToTmunu as ScalarFieldTmunu after Whisky_SetTmunu"},*)
  Where -> Interior,
  ConditionalOnKeywords -> {{"fd_order", fdOrder},{"couple_scalar_field_to_gravity","yes"},{"potential_type",potential}},
  Shorthands -> shorthands,
  Equations ->
  {
     
      phiamp -> (phia^2 + phib^2)^(1/2)

  }
};

(**************************************************************************************)
(* Boundary conditions *)
(**************************************************************************************)

boundaryParam = 
{
  Name -> "boundary_condition",
  Default -> "radiative",
  AllowedValues -> {"none", "radiative", "phi_of_t"}
};

boundaryCalc =
{
  Name -> fnPrefix <> "_boundaryRadiative",
  Schedule -> {"in MoL_RHSBoundaries as ScalarFieldBoundaries"},
  ConditionalOnKeyword -> {"boundary_condition", "radiative"},
  Where -> Boundary,
  Shorthands -> shorthands,
  Equations  -> 
  {
    n1 -> -x/r, 
    n2 -> -y/r,
    n3 -> -z/r,

    dir[ui] -> Sign[n[ui]],

    (* \partial_t u = - (u - u_0) / r - \partial_r u,
       for u_0 some background solution. In this case, Minkowski in Cartesian 
       coordinates. *)
    
    (* WARNING: This is only second order accurate *)

    (* Scalar Field boundary conditions *)
    dot[pia]  -> -(pia - piaBG)/r + n[uk] PDonesided2nd[pia,lk],
    dot[phia] -> -(phia-phiaBG)/r + n[uk] PDonesided2nd[phia,lk],
    
    dot[pib]  -> -(pib - pibBG)/r + n[uk] PDonesided2nd[pib,lk],
    dot[phib] -> -(phib-phibBG)/r + n[uk] PDonesided2nd[phib,lk],

    (* clear rho and scalar field amplitude in boundaries to get rid of poison *)
    rhoE -> 0,
    phiamp -> 0

  }
};

boundaryPhiOfTCalc =
{
  Name -> fnPrefix <> "_boundaryPhiOfT",
  Schedule -> {"in MoL_RHSBoundaries as ScalarFieldBoundaries"},
  ConditionalOnKeyword -> {"boundary_condition", "phi_of_t"},
  Where -> Boundary,
  Shorthands -> shorthands,
  Equations  -> 
  {
    n1 -> -x/r,
    n2 -> -y/r,
    n3 -> -z/r,

    dir[ui] -> Sign[n[ui]],

    (* \partial_t u = - (u - u_0) / r - \partial_r u,
       for u_0 some background solution. In this case, Minkowski in Cartesian 
       coordinates. *)

    (* WARNING: This is only second order accurate *)

    (* Scalar Field boundary conditions *)
    dot[pia]  -> -(pia - piaBG  )/r + n[uk] PDonesided2nd[pia,lk],
    dot[phia] -> -(phia- piaBG*t)/r + n[uk] PDonesided2nd[phia,lk],
    
    dot[pib]  -> -(pib - pibBG  )/r + n[uk] PDonesided2nd[pib,lk],
    dot[phib] -> -(phib- pibBG*t)/r + n[uk] PDonesided2nd[phib,lk],

    (* clear rho and scalar field amplitude in boundaries to get rid of poison *)
    rhoE -> 0,
    phiamp -> 0

  }
};

(**************************************************************************************)
(* Construct the thorn *)
(**************************************************************************************)



calculations =
{
  initGFsCalc,
  evolveScalarField["full4th", PDstandard4th, PDlopsided4th,"none",Vnone],
  evolveScalarField["6th", PDstandard6th, PDlopsided6th,"none",Vnone],
  evolveScalarField["full4th", PDstandard4th, PDlopsided4th,"onlymass",Vonlymass],
  evolveScalarField["6th", PDstandard6th, PDlopsided6th,"onlymass",Vonlymass],
  evolveScalarField["full4th", PDstandard4th, PDlopsided4th,"axion",Vaxion],
  evolveScalarField["6th", PDstandard6th, PDlopsided6th,"axion",Vaxion],
  evolveScalarField["full4th", PDstandard4th, PDlopsided4th,"doublewell",Vdoublewell],
  evolveScalarField["6th", PDstandard6th, PDlopsided6th,"doublewell",Vdoublewell],
  evolveScalarField["full4th", PDstandard4th, PDlopsided4th,"phiFourth",VphiFourthall],
  evolveScalarField["6th", PDstandard6th, PDlopsided6th,"phiFourth",VphiFourthall],
  evolveScalarField["full4th", PDstandard4th, PDlopsided4th,"inflation",Vphiinflation],
  evolveScalarField["6th", PDstandard6th, PDlopsided6th,"inflation",Vphiinflation],
  evolveScalarField["full4th", PDstandard4th, PDlopsided4th,"linear",Vphilinear],
  evolveScalarField["6th", PDstandard6th, PDlopsided6th,"linear",Vphilinear],
  
(*  AddPotentialToRHS["none",Vnone],
  AddPotentialToRHS["onlymass",Vonlymass],
  AddPotentialToRHS["phiFourth",VphiFourthall],
  AddPotentialToRHS["doublewell",Vdoublewell],
  AddPotentialToRHS["inflation",Vphiinflation],
  AddPotentialToRHS["linear",Vphilinear],*)

  ScalarFieldTmunuCalc["full4th", PDstandard4th,PDlopsided4th,"none",Vnone],
  ScalarFieldTmunuCalc["6th", PDstandard6th,PDlopsided6th,"none",Vnone],
  ScalarFieldTmunuCalc["full4th", PDstandard4th,PDlopsided4th,"onlymass",Vonlymass],
  ScalarFieldTmunuCalc["6th", PDstandard6th,PDlopsided6th,"onlymass",Vonlymass],
  ScalarFieldTmunuCalc["full4th", PDstandard4th,PDlopsided4th,"axion",Vaxion],
  ScalarFieldTmunuCalc["6th", PDstandard6th,PDlopsided6th,"axion",Vaxion],
  ScalarFieldTmunuCalc["full4th", PDstandard4th,PDlopsided4th,"doublewell",Vdoublewell],
  ScalarFieldTmunuCalc["6th", PDstandard6th,PDlopsided6th,"doublewell",Vdoublewell],
  ScalarFieldTmunuCalc["full4th", PDstandard4th,PDlopsided4th,"phiFourth",VphiFourthall],
  ScalarFieldTmunuCalc["6th", PDstandard6th,PDlopsided6th,"phiFourth",VphiFourthall],
  ScalarFieldTmunuCalc["full4th", PDstandard4th,PDlopsided4th,"inflation",Vphiinflation],
  ScalarFieldTmunuCalc["6th", PDstandard6th,PDlopsided6th,"inflation",Vphiinflation],
  ScalarFieldTmunuCalc["full4th", PDstandard4th,PDlopsided4th,"linear",Vphilinear],
  ScalarFieldTmunuCalc["6th", PDstandard6th,PDlopsided6th,"linear",Vphilinear],


  ScalarFieldrhoCalc["full4th", PDstandard4th,PDlopsided4th,"none",Vnone],
  ScalarFieldrhoCalc["6th", PDstandard6th,PDlopsided6th,"none",Vnone],
  ScalarFieldrhoCalc["full4th", PDstandard4th,PDlopsided4th,"onlymass",Vonlymass],
  ScalarFieldrhoCalc["6th", PDstandard6th,PDlopsided6th,"onlymass",Vonlymass],
  ScalarFieldrhoCalc["full4th", PDstandard4th,PDlopsided4th,"axion",Vaxion],
  ScalarFieldrhoCalc["6th", PDstandard6th,PDlopsided6th,"axion",Vaxion],
  ScalarFieldrhoCalc["full4th", PDstandard4th,PDlopsided4th,"doublewell",Vdoublewell],
  ScalarFieldrhoCalc["6th", PDstandard6th,PDlopsided6th,"doublewell",Vdoublewell],
  ScalarFieldrhoCalc["full4th", PDstandard4th,PDlopsided4th,"phiFourth",VphiFourthall],
  ScalarFieldrhoCalc["6th", PDstandard6th,PDlopsided6th,"phiFourth",VphiFourthall],
  ScalarFieldrhoCalc["full4th", PDstandard4th,PDlopsided4th,"inflation",Vphiinflation],
  ScalarFieldrhoCalc["6th", PDstandard6th,PDlopsided6th,"inflation",Vphiinflation],
  ScalarFieldrhoCalc["full4th", PDstandard4th,PDlopsided4th,"linear",Vphilinear],
  ScalarFieldrhoCalc["6th", PDstandard6th,PDlopsided6th,"linear",Vphilinear],  

  ScalarFieldAmplitudeCalc["full4th", PDstandard4th,PDlopsided4th,"none",Vnone],
  ScalarFieldAmplitudeCalc["6th", PDstandard6th,PDlopsided6th,"none",Vnone],
  ScalarFieldAmplitudeCalc["full4th", PDstandard4th,PDlopsided4th,"onlymass",Vonlymass],
  ScalarFieldAmplitudeCalc["6th", PDstandard6th,PDlopsided6th,"onlymass",Vonlymass],
  ScalarFieldAmplitudeCalc["full4th", PDstandard4th,PDlopsided4th,"axion",Vaxion],
  ScalarFieldAmplitudeCalc["6th", PDstandard6th,PDlopsided6th,"axion",Vaxion],
  ScalarFieldAmplitudeCalc["full4th", PDstandard4th,PDlopsided4th,"doublewell",Vdoublewell],
  ScalarFieldAmplitudeCalc["6th", PDstandard6th,PDlopsided6th,"doublewell",Vdoublewell],
  ScalarFieldAmplitudeCalc["full4th", PDstandard4th,PDlopsided4th,"phiFourth",VphiFourthall],
  ScalarFieldAmplitudeCalc["6th", PDstandard6th,PDlopsided6th,"phiFourth",VphiFourthall],
  ScalarFieldAmplitudeCalc["full4th", PDstandard4th,PDlopsided4th,"inflation",Vphiinflation],
  ScalarFieldAmplitudeCalc["6th", PDstandard6th,PDlopsided6th,"inflation",Vphiinflation],
  ScalarFieldAmplitudeCalc["full4th", PDstandard4th,PDlopsided4th,"linear",Vphilinear],
  ScalarFieldAmplitudeCalc["6th", PDstandard6th,PDlopsided6th,"linear",Vphilinear],

  boundaryCalc,
  boundaryPhiOfTCalc
};



keywordParameters = 
{
  fdOrderParam,
  boundaryParam,
  potentialParam,
  computeScalarFieldTmunuParam
};

(*extendedKeywordParameters =
{
  ADMBaseIDParam,
  ADMBaseLapseIDParam,
  ADMBaseShiftIDParam
};*)

CreateKrancThornTT[groups, ".", thornName,
  Calculations -> calculations,
  DeclaredGroups -> declaredGroupNames,
  PartialDerivatives -> derivatives,
  RealParameters -> realParameters,
  KeywordParameters -> keywordParameters,
  InheritedImplementations -> {"admbase","TmunuBase"}];


