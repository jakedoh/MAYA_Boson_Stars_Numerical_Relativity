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
  n, dir, beta, kcurv, Tmn, betaD, rhoE, rhoJ,
  g, gLocal, alpLocal, betaLocal, Gam3, Gam0ij, Gamkij,
  STpi, STphi, STSi, dtSTphi, 
  STTtj, STTij, (*STphid,*) (*STpid,*)
  exPhidd, exPhid, exPid, alpha,
  TmunuTraceMatterLocal, gLocalInv, Tmn, TmnLocal, kLocal, gInv, KtraceLocal, ktrace, ttrace
}];

Map[AssertSymmetricIncreasing,
{
  gLocalInv[ua,ub], gInv[ua,ub]
}];

Map[AssertSymmetricDecreasing, 
{
  gLocal[la,lb], kcurv[la,lb], kLocal[la,lb], STTij[la,lb], exPhidd[la,lb], Tmn[la,lb], Gam0ij[la,lb],
  TmnLocal[la,lb], g[la,lb]
}];

AssertSymmetricDecreasing[Gam3[ua,lb,lc], lb, lc];
AssertSymmetricDecreasing[Gamkij[ua,lb,lc], lb, lc];

(* Determinants of the metrics in terms of their components
   (Mathematica symbolic expressions) *)
gDet = Det[MatrixOfComponents[gLocal[la,lb]]];


fnPrefix = "ScalarTensor";
thornName = "ScalarTensor";

(**************************************************************************************)
(* Groups *)
(**************************************************************************************)

(* Cactus group definitions *)

evolvedGroups = Map[CreateGroupFromTensor, {STphi, STpi(*, dtSTphi, STphid[li]*)}];

admGroups = 
  {{"admbase::metric", {gxx,gxy,gxz,gyy,gyz,gzz}},
   {"admbase::curv",  {kxx,kxy,kxz,kyy,kyz,kzz}},
   {"admbase::lapse", {alp}},
   {"admbase::shift", {betax,betay,betaz}}};

evaluatedGroups = Map[CreateGroupFromTensor,{gLocal[li,lj], alpLocal, KtraceLocal,(*STpid[li],*) (*aofphi,*) TmunuTraceMatterLocal,rhoE,rhoJ}];
evaluatedGroups = Map[AddGroupExtra[#, Timelevels -> 3] &, evaluatedGroups];

stressenergyGroups =
  {{"TmunuBase::stress_energy_scalar", {eTtt} },
   {"TmunuBase::stress_energy_vector", {eTtx,eTty,eTtz} },
   {"TmunuBase::stress_energy_tensor", {eTxx,eTxy,eTxz,eTyy,eTyz,eTzz} } };

declaredGroups = Join[evolvedGroups, evaluatedGroups];
declaredGroupNames = Map[First, declaredGroups];

groups = Join[declaredGroups, admGroups, stressenergyGroups];

(**************************************************************************************)
(* Shorthands *)
(**************************************************************************************)

shorthands = 
{
  ktrace, ttrace,
  n[ui], dir[ui],  detg, invdetg, 
  gLocalInv[ui,uj], betaD[li], gInv[ui,uj],
  TttLocal, TtxLocal, TtyLocal, TtzLocal, TmnLocal[li,lj],
  betaSqr, dSTphiSqr, dotSTpi, pi,
  STTtt, STTtj[li], STTij[li,lj], STrho, STSi[li],
  exPhi, exPi, exPhid[li], exPhidd[li,lj], exPid[li], sourcePhi, sourcePi, exPhit, exPit,
  Gam3[ua,lb,lc], Gam0ij[li,lj], Gamkij[uk,li,lj], betaLocal[ui],
  ConfFactorA, ConfFactorA2, ConfFactorAm2, kLocal[li,lj], Ktrace, aofphi, dtSTphi, derivV, phiRS, piRS
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
  {Name -> STomega, Default -> 10},
  {Name -> STBeta, Default -> 0 },
  {Name -> STalpha0, Default -> 0 },
  {Name -> STamplitude, Default -> 0},
  {Name -> r0, Default -> 0},
  {Name -> STsigma, Default -> 1},
  {Name -> kx, Default -> 1},
  {Name -> ky, Default -> 1},
  {Name -> kz, Default -> 1},
  {Name -> x00, Default -> 0},
  {Name -> y00, Default -> 0},
  {Name -> z00, Default -> 0},
  {Name -> fieldMass, Default -> 0 },
  {Name -> phiBG, Default -> 0},
  {Name -> piBG, Default -> 0},
  {Name -> Vlambda, Default -> 0},
  {Name -> Veps, Default -> 0},
  {Name -> Veta, Default -> 0},
  {Name -> V0, Default -> 0}
};

(* Brans Dicke scalar field or normal scalar field *)
(*systemTypeParam =
{
  Name -> "system_type",
  Default -> "BransDicke",
  AllowedValues -> {"BransDicke", "ScalarField"}
};*)

calcRhoParam =
{
  Name -> "calculate_rho",
  Default -> "no",
  AllowedValues -> {"yes", "no"}
};

(*(* form of A(phi) *)
AofphiParam =
{
  Name -> "A_of_phi",
  Default -> "exponential phi squared",
  AllowedValues -> {"exponential phi squared", "ST"}
};*)

(* If this is set to 'no', then a function that puts zeros in  TmunuTraceMatterLocal is scheduled. *)
computeScalarTmunuParam =
{
  Name -> "couple_scalar_field_to_matter",
  Default -> "yes",
  AllowedValues -> {"yes", "no"}
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
  AllowedValues -> {"2nd", (* "4th", *) "full4th", "4thcentred", "6thcentred", "6th"}
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
    STphi      -> INITVALUE,
    STpi       -> INITVALUE,
    rhoE       -> 0,
    rhoJ       -> 0
    (*dtSTphi    -> INITVALUE*)
    (*STphid[li] -> INITVALUE*)
  }
};

(**************************************************************************************)
(* conformal factors *)
(**************************************************************************************)

STA[phi_]:= Exp[ STalpha0 (phi - phiBG) + 0.5 STBeta (phi^2-phiBG^2) ];
aofphi[A_,x_]:= D[Log[A[phi]],phi]/.phi->x;

(**************************************************************************************)
(* potential functions *)
(**************************************************************************************)
(*VcosmoSF[phi_]:= Vlambda/8*(phi^2-Veta^2)^2;*)
(*VcosmoST[phi_]:= Vlambda/8*( (2 Sqrt[Pi] phi)^2 - ( 2 Sqrt[Pi] Veta)^2 )^2;*)
(*VcosmoST[phi_]:= Vlambda/8*( (2 Sqrt[Pi] phi)^2 - ( 2 Sqrt[Pi] Veta)^2 )^2 * ( (2 Sqrt[Pi] phi)^2 - ( 2 Sqrt[Pi] phiBG)^2 )^2;*)
(*VcosmoST[phi_]:= Vlambda/8*( phi^2 - Veta^2 )^2 * ( phi^2 - phiBG^2 )^2;*)
Vfourth[phi_]:= Vlambda/8*( phi^2 - phiBG^2 )^2;
Vsixth[phi_] := Vlambda/8*( phi^2 - phiBG^2 )^2 * phi^2;
Veighth[phi_] := Vlambda/8*( phi^2 - phiBG^2 )^2 * ( phi^2 - Veta^2 )^2; 
Vbubble[phi_] := Vlambda/8*( ( phi^2 - phiBG^2 )^2 - 4*Veps*phiBG^3*(phi - phiBG) );
Vexpand[phi_] := Vlambda*phi^2 * ( 1/6 phi^4 - 1/4 ( 1 + Veps^2 ) Veta^2 phi^2 + 1/2 Veps^2 Veta^4 ) + V0;
Vnone[phi_]  := 0;
Vprime[V_,x_]:= D[V[phi],phi]/.phi->x;

potentialParam =
{
  Name -> "potential_type",
  Default -> "none",
  AllowedValues -> {"phiFourth", "phiSixth", "phiEighth", "phiBubble", "phiExpand", "none"}
};

(**************************************************************************************)
(* Scalar Field Initial Data *)
(**************************************************************************************)
(*
initialDataParam = 
{
  Name -> "initialData_type",
  Default -> "Gaussian",
  AllowedValues -> {"Gaussian", "none", "TravelingPulse", "static"}
};

ep=0.0000000001;
IDtemp  = 1 + STamplitude Exp[-(Sqrt[(x-x00)^2+(y-y00)^2+(z-z00)^2+ep]-r0)^2/STsigma^2];
IDtempx = D[IDtemp,x];
IDtempy = D[IDtemp,y];
IDtempz = D[IDtemp,z];

initialDataCalc1 = 
{
  Name -> fnPrefix <> "_InitialDataGaussian",
  Schedule -> {"in ADMBase_PostInitial as STInitialData"},
  ConditionalOnKeyword -> {"initialData_type", "Gaussian"},
  Shorthands -> shorthands,
  Equations ->
  {
    STphi  -> STamplitude Exp[-(Sqrt[(x-x00)^2+(y-y00)^2+(z-z00)^2+ep]-r0)^2/STsigma^2],

    (* This is necesary if we want dot[STphi] to be zero at t=0 *)
    STpi    -> -(betax IDtempx + betay IDtempy + betaz IDtempz )/alp
   }
};

(* Initial data from the static single BH solution *)
initialDataStatSolCalc = 
{
  Name -> fnPrefix <> "_InitiDataStaticBransDicke",
  Schedule -> {"in ADMBase_PostInitial as STInitialData"},
  ConditionalOnKeyword -> {"initialData_type", "static"},
  Shorthands -> shorthands,
  Equations ->
  {
    STphi  -> STamplitude/Sqrt[(x-1.168642873)^2+y^2+z^2+0.1] + STamplitude/Sqrt[(x+1.168642873)^2+y^2+z^2+0.1],

    STpi    -> 0
   }
};


(* Exact solution to test the code *)
(* these definitions are harmless, I will leave them here uncommented *)
lx=ly=lz=2*10.0; (* this actually depends on your parameter file *)

exxPhi  = Cos[t] Cos[2 Pi x/lx] Cos[2 Pi y/ly] Cos[2 Pi z/lz];
exxPhit = D[exxPhi, t];
exxPhix = D[exxPhi, x];
exxPhiy = D[exxPhi, y];
exxPhiz = D[exxPhi, z];

exxPi   = ( exxPhit - betax exxPhix - betay exxPhiy - betaz exxPhiz ) / alp;
exxPit  = D[exxPi, t];
exxPix  = D[exxPi, x];
exxPiy  = D[exxPi, y];
exxPiz  = D[exxPi, z];

exxPhixx = D[exxPhi, x, x];
exxPhiyx = D[exxPhi, y, x];
exxPhizx = D[exxPhi, z, x];
exxPhiyy = D[exxPhi, y, y];
exxPhizy = D[exxPhi, z, y];
exxPhizz = D[exxPhi, z, z];

IDtemp2  = STamplitude Exp[-Sqrt[(x-x00)^2+(y-y00)^2+(z-z00)^2+ep]^2/STsigma^2];
IDtempx2 = D[IDtemp2,x];
IDtempy2 = D[IDtemp2,y];
IDtempz2 = D[IDtemp2,z];

initialDataCalc2 = 
{
  Name -> fnPrefix <> "_InitialDataTravelingPulse",
  (*  Schedule -> {"at INITIAL after ADMBase_PostInitial after bssn_evolved_to_adm_initial"},*)
  Schedule -> {"in ADMBase_PostInitial"},
  ConditionalOnKeyword -> {"initialData_type", "TravelingPulse"},
  Shorthands -> shorthands,
  Equations ->
  {
    (* The next line enables initial data to test the code
     whit an ad hoc exact solution *)
    (* STphi -> exxPhi, STpi -> exxPi *)



    STphi  -> STamplitude Exp[-Sqrt[(x-x00)^2+(y-y00)^2+(z-z00)^2+ep]^2/STsigma^2],
    (* This is necesary if we want dot[STphi] to be zero at t=0 *)
    STpi -> -(betax IDtempx + betay IDtempy + betaz IDtempz )/alpLocal 
		   - 2 STamplitude Exp[-Sqrt[(x-x00)^2+(y-y00)^2+(z-z00)^2+ep]^2/STsigma^2] (x kx + y ky + z kz)/STsigma^2/alpLocal

   }
};


ADMBaseIDParam =
{
  Name -> "ADMBase::initial_data",
  Default -> {};
  AllowedValues -> {"StaticBransDickeBH"}
};

initialDataStaticST =
{
  Name -> fnPrefix <> "_StaticBransDickeBH",
  Schedule -> {"in ADMBase_InitialData"},
  ConditionalOnKeyword -> {"initial_data", "StaticBransDickeBH"},
  Shorthands -> shorthands,
  Equations ->
  {
    gxx -> 1
  }

};

ADMBaseLapseIDParam =
{
  Name -> "ADMBase::initial_lapse",
  Default -> {};
  AllowedValues -> {"StaticBransDickeBH"}
};

initialDataStaticST =
{
  Name -> fnPrefix <> "_StaticBransDickeBH",
  Schedule -> {"in ADMBase_InitialGauge"},
  ConditionalOnKeyword -> {"initial_lapse", "StaticBransDickeBH"},
  Shorthands -> shorthands,
  Equations ->
  {
    alpLocal -> Sqrt[1 - 2/r (1 + 1)]
  }

};

   
ADMBaseShiftIDParam =
{
  Name -> "ADMBase::initial_shift",
  Default -> {};
  AllowedValues -> {"StaticBransDickeBH"}
};

initialDataStaticST =
{
  Name -> fnPrefix <> "_StaticBransDickeBH",
  Schedule -> {"in ADMBase_InitialGauge"},
  ConditionalOnKeyword -> {"initial_shift", "StaticBransDickeBH"},
  Shorthands -> shorthands,
  Equations ->
  {
    beta[ui] -> 0
  }

};


(*computeDtSTphiCalc[fdOrder_, PD_, PDadvect_] := 
{
  Name -> fnPrefix<>"_DtSTphi_"<>fdOrder,
  Schedule -> {"AT Initial after ADMBase_PostInitial before EinsteinToJordanFrame as ScalarFieldDerivatives"},
	       (*"in MoL_PostStep before EinsteinToJordanFrame after bssn_evolved_to_adm as ScalarFieldDerivatives"}, *)
  (* "at CCTK_EVOL after MoL_Evolution before EinsteinToJordanFrame as ScalarFieldDerivatives"},*)
  ConditionalOnKeyword -> {"fd_order", fdOrder},
  Shorthands -> shorthands,
  Equations -> 
  {
    (* This function equation was inside evolution but it only goes through the Interior *)
    (* points. Here we go over Everywhere so there are no nans in the boundaries *)

      (*ConfFactorA   -> Exp[ 0.5 STBeta STphi^2 ],

      STphid[li] -> PD[STphi,li],*)
      dir[ui] -> Sign[beta[ui]],
      dtSTphi   -> alpLocal STpi + beta[ui] PDadvect[STphi,li]
  }
};*)

(*
"AT Initial after STInitialData as ST_ADMToLocal", 
	       "in MoL_PostStep before STEvolveScalarField as ST_ADMToLocal", 
	       "at CCTK_EVOL after MoL_Evolution before STEvolveScalarField as ST_ADMToLocal"
 *)
*)

(**************************************************************************************)
(* Handling of stored grid functions and conformal transformation*)
(**************************************************************************************)

ADMToLocalCalc = 
{
  Name -> fnPrefix <> "_ADMToLocal",
  Schedule -> {"AT Initial after ADMBase_PostInitial before EinsteinToJordanFrame as ScalarFieldADMToLocal_initial", 
	       "in MoL_PostStep after bssn_evolved_to_adm before EinsteinToJordanFrame as ScalarFieldADMToLocal", 
	       "at CCTK_EVOL after MoL_Evolution before EinsteinToJordanFrame as ScalarFieldADMToLocal"},
  (*ConditionalOnKeyword -> {"system_type", "BransDicke"},*)
  Shorthands -> shorthands,
  Equations ->
   {
     alpLocal  -> alpha,
     gLocal[li, lj] -> g[li,lj],
     kLocal[li,lj] -> kcurv[li,lj],

     detg -> gDet,
     invdetg -> 1/detg,
     gLocalInv[ui,uj] -> invdetg gDet MatrixInverse[gLocal[ui,uj]],

     KtraceLocal -> kLocal[li,lj] gLocalInv[ui, uj]
   }
};

EinsteinToJordanCalc=
{
  Name -> fnPrefix <> "_EinsteinToJordan",
  Schedule ->  {"AT Initial as EinsteinToJordanFrame after bssn_evolved_to_adm_initial before HydroBase_Prim2ConInitial",
                "in MoL_PostStep as EinsteinToJordanFrame after bssn_evolved_to_adm before HydroBase_PostStep",
                "at CCTK_EVOL as EinsteinToJordanFrame after bssn_evolved_to_adm" },
  Shorthands -> shorthands,
  (*ConditionalOnKeyword -> {"system_type", "BransDicke"},*)
  Equations ->
  {
    (*ConfFactorA -> Exp[ 0.5 STBeta STphi^2 ],
    aofphi -> STBeta STphi,*)
    ConfFactorA -> STA[STphi],
    aofphi -> aofphi[STA,STphi],

    alp -> ConfFactorA alp,
    (* factor of 2 sqrt(pi) comes with STpi due to rescaling, g_ij is still in the einstein frame *)
    kcurv[li,lj] -> ConfFactorA ( kcurv[li,lj] - g[li,lj] STpi aofphi ),
    g[li, lj] -> ConfFactorA^2 g[li,lj]
  }
};

(**************************************************************************************)
(* Evolution of the scalar field*)
(**************************************************************************************)

evolveScalarField[fdOrder_, PD_, PDadvect_] := 
{
  Name -> fnPrefix<>"_EvolveScalarField_"<>fdOrder,
  Schedule -> {"in MoL_CalcRHS as EvolveScalarField"},
  Where -> Interior,
  ConditionalOnKeyword -> {"fd_order", fdOrder},
  Shorthands -> shorthands,
  (*  CollectList -> {gLocalInv[ui,uj]}, *)
  Equations -> 
  {
    dir[ui] -> Sign[beta[ui]],
    pi -> 3.14159265358979323846,
 
(*
    detg -> gDet,
    invdetg -> 1 / detg,
*)
    gInv[ui,uj] -> MatrixInverse[g[ui,uj]],
    ktrace -> kcurv[li,lj] gInv[ui, uj],
    (*kLocal[li,lj] -> kcurv[li,lj] / ConfFactorA + STBeta STphi / (ConfFactorA alp)*
    ( gLocal[li,lj] ConfFactorA^2 dtSTphi - betaD[lj] PD[STphi,li] - betaD[li] PD[STphi,lj] ), *)
    (*kLocal[li,lj] -> ( kcurv[li,lj]  + gLocal[li,lj] STpi aofphi ) / ConfFactorA,*)

    Gam3[ua, lb, lc] -> 1/2 gInv[ua,ud] (PD[g[lb,ld], lc] + PD[g[lc,ld], lb] - PD[g[lb,lc],ld]),
     
	(*ConfFactorA2  -> Exp[ STBeta STphi^2 ],*)
     ConfFactorA2 -> STA[STphi]^2,
     TttLocal -> eTtt ConfFactorA2,
     TtxLocal -> eTtx ConfFactorA2,
     TtyLocal -> eTty ConfFactorA2,
     TtzLocal -> eTtz ConfFactorA2,
     TmnLocal[li,lj] -> Tmn[li,lj] ConfFactorA2,
     ttrace -> -TttLocal/alpha^2 + 2/alpha^2 (TtxLocal beta1 + TtyLocal beta2 + TtzLocal beta3 )
                           -beta[ui] beta[uj] TmnLocal[li,lj]/alpha^2 + gInv[ui,uj] TmnLocal[li,lj] ,
    
    (* scalar field evolution *)    
    dot[STphi] -> alpha STpi + beta[ui] PDadvect[STphi,li], (* partial time-derivative of STphi; STpi is the Lie derivative of STphi along n*)
    dot[STpi]  -> beta[ui] PDadvect[STpi,li] + alpha ktrace STpi
      + alpha gInv[ui,uj] PD[STphi,li,lj] (*comment this line for Michael orig*) 
      - alpha gInv[ui,uj] Gam3[uk,li,lj] PD[STphi,lk] + gInv[ui,uj] PD[alpha,lj] PD[STphi,li]
      + alpha aofphi[STA,STphi] ttrace
(*
     - 0.5 alpha gInv[ui,uj] gInv[ul,um] PD[g[ll,lm],lj] PD[STphi,li]
     uncomment this line for Michael orig; indices are already corrected
*)     
    
(*
    + alpLocal gLocalInv[ui,uj] PD[STphi,li,lj] 
    (*- 0.5 alpLocal gLocalInv[ui,uj] gLocalInv[ul,um] PD[gLocal[ll,lm],lj] PD[STphi,li] *)
    - alpLocal gLocalInv[ui,uj] Gam3[uk,li,lj] PD[STphi,lk]
*)
	(*0 + gLocalInv[ui,uj] PD[alpLocal,lj] PD[STphi,li]  *)
                 (*- alpLocal fieldMass^2 STphi*)

    (*,dot[STphi] -> STpi,
     dot[STpi]  -> Euc[ua,ub] PD[STphi,la,lb] *)

  }
};

AddPotentialToRHS[potential_,V_] :=
{
  Name -> fnPrefix<>"_AddPotentialToRHS_"<>potential,
  Schedule -> {"in MoL_CalcRHS after EvolveScalarField as AddPotentialToScalarFieldRHS"},
  Where -> Interior,
  ConditionalOnKeyword -> {"potential_type", potential},
  Shorthands -> shorthands,
  Equations ->
  {
    derivV -> Vprime[V,STphi],
    dot[STpi]  -> dot[STpi] - derivV
  }
};

SetTmunuTraceToZeroCalc =
{
  Name -> fnPrefix <> "_ZeroTmunuTrace",
  Schedule -> {"in AddToTmunu as ScalarFieldTmunu after Whisky_SetTmunu"},
  Where -> Everywhere,
  ConditionalOnKeyword -> {"couple_scalar_field_to_matter", "no"},
  Shorthands -> shorthands,
  Equations -> {
    TmunuTraceMatterLocal -> 0
  }
};

ScalarFieldTmunuCalc[fdOrder_,PD_, PDadvect_, potential_, V_] := 
{
  Name -> fnPrefix <> "_ScalarFieldTmunu_" <> potential <> "_" <> fdOrder,
  Schedule -> {"in AddToTmunu as ScalarFieldTmunu after Whisky_SetTmunu"},
  Where -> Interior,
  (*ConditionalOnKeywords -> {{"fd_order", fdOrder},{"couple_scalar_field","yes"},{"system_type","BransDicke"},{"potential_type",potential}},*)
  ConditionalOnKeywords -> {{"fd_order", fdOrder},{"couple_scalar_field_to_gravity","yes"},{"potential_type",potential}},
  Shorthands -> shorthands,
  Equations -> 
  {
     dir[ui] -> Sign[beta[ui]],
     pi -> 3.14159265358979323846,
     
     detg -> gDet,
     invdetg -> 1 / detg,
     gLocalInv[ua,ub] -> invdetg gDet MatrixInverse[gLocal[ua,ub]],
     betaLocal[ui] -> beta[ui],

     (* time derivative from eq of motion *)
     dtSTphi -> alpLocal STpi + betaLocal[ui] PDadvect[STphi,li],

     (* auxiliary variables *)
     betaSqr -> gLocal[li,lj] betaLocal[ui] betaLocal[uj],
     (*dSTphiSqr -> dtSTphi^2/alpLocal^2 - 2 STpi dtSTphi/alpLocal + gLocalInv[ui,uj] PD[STphi,li] PD[STphi,lj],*)
     dSTphiSqr -> gLocalInv[ui,uj] PD[STphi,li] PD[STphi,lj] - STpi^2, 
     
     (* Components due to the scalar field *)
     STTtt -> dtSTphi^2 - 1/2 (-alpLocal^2 + betaSqr) ( dSTphiSqr + fieldMass^2 STphi^2 + 2 V[STphi] ),
     STTtj[li] -> dtSTphi PD[STphi,li] - 1/2 gLocal[li,lj] betaLocal[uj] (dSTphiSqr + fieldMass^2 STphi^2 + 2 V[STphi] ),
     STTij[li,lj] -> PD[STphi,li] PD[STphi,lj] - 1/2 gLocal[li,lj] (dSTphiSqr + fieldMass^2 STphi^2 + 2 V[STphi] ),

     (****************************************************)
     (* Compute Tmunu in Einstein Frame *)
     (*ConfFactorA2  -> Exp[ STBeta STphi^2 ],*)
     ConfFactorA2 -> STA[STphi]^2,
     TttLocal -> eTtt ConfFactorA2,
     TtxLocal -> eTtx ConfFactorA2,
     TtyLocal -> eTty ConfFactorA2,
     TtzLocal -> eTtz ConfFactorA2,
     TmnLocal[li,lj] -> Tmn[li,lj] ConfFactorA2,
     (****************************************************)
     
     (*======================================================================================*)
     (* Trace of the stress energy tensor of matter. At this point assuming this is the last *)
     (* function that adds something to Tmunu, all the matter contributions are in place.    *)
     (* Therefore the following trace is only from the matter fields.                        *)
     (* We compute it here and store it to use it later in the evolution equations           *)
     (*======================================================================================*)
     TmunuTraceMatterLocal -> -TttLocal/alpLocal^2 + 2/alpLocal^2 (TtxLocal betaLocal1 + TtyLocal betaLocal2 + TtzLocal betaLocal3 )
                           -betaLocal[ui] betaLocal[uj] TmnLocal[li,lj]/alpLocal^2 + gLocalInv[ui,uj] TmnLocal[li,lj] ,

     (* Add to TmunuBase variables *)
     (* STTtt, STTtj, STTij are in the Einstein Frame. Here we make sure that the scalar *)
     (* field is the last contribution to the total Tmunu and transform Tmunu to the     *)
     (* Einstein frame. In this way we minimize modification the the BSSN thorn          *)
         
     eTtt -> TttLocal   + STTtt  ,
     eTtx -> TtxLocal   + STTtj1 ,
     eTty -> TtyLocal   + STTtj2 ,
     eTtz -> TtzLocal   + STTtj3 ,
     eTxx -> TmnLocal11 + STTij11,
     eTxy -> TmnLocal21 + STTij21,
     eTxz -> TmnLocal31 + STTij31,
     eTyy -> TmnLocal22 + STTij22,
     eTyz -> TmnLocal32 + STTij32,
     eTzz -> TmnLocal33 + STTij33
     
  }
};

(*ScalarFieldTmunuCalc[fdOrder_,PD_, PDadvect_, potential_,V_] :=
{
  Name -> fnPrefix <> "_ScalarFieldTmunu_" <> potential <> "_" <> fdOrder,
  Schedule -> {"in AddToTmunu as ScalarFieldTmunu after Whisky_SetTmunu"},
  Where -> Interior,
  ConditionalOnKeywords -> {{"fd_order", fdOrder},{"couple_scalar_field","yes"},{"system_type","ScalarField"},{"potential_type",potential}},
  Shorthands -> shorthands,
  Equations ->
  {
     dir[ui] -> Sign[beta[ui]],
     pi -> 3.14159265358979323846,

     detg -> gDet,
     invdetg -> 1 / detg,
     gLocalInv[ua,ub] -> invdetg gDet MatrixInverse[gLocal[ua,ub]],
     betaLocal[ui] -> beta[ui],

     (* time derivative from eq of motion *)
     dtSTphi -> alpLocal STpi + betaLocal[ui] PDadvect[STphi,li],


     (* auxiliary variables *)
     betaSqr -> gLocal[li,lj] betaLocal[ui] betaLocal[uj],
     (*dSTphiSqr -> dtSTphi^2/alpLocal^2 - 2 STpi dtSTphi/alpLocal + gLocalInv[ui,uj] PD[STphi,li] PD[STphi,lj],*)
     dSTphiSqr -> gLocalInv[ui,uj] PD[STphi,li] PD[STphi,lj] - STpi^2,
     
     (* Components due to the scalar field *)
     STTtt -> dtSTphi^2 - 1/2 (-alpLocal^2 + betaSqr) ( dSTphiSqr + fieldMass^2 STphi^2 + 2 V[STphi] ),
     STTtj[li] -> dtSTphi PD[STphi,li] - 1/2 gLocal[li,lj] betaLocal[uj] (dSTphiSqr + fieldMass^2 STphi^2 + 2 V[STphi]),
     STTij[li,lj] -> PD[STphi,li] PD[STphi,lj] - 1/2 gLocal[li,lj] (dSTphiSqr + fieldMass^2 STphi^2 + 2 V[STphi]),

     (****************************************************)
     TttLocal -> eTtt,
     TtxLocal -> eTtx,
     TtyLocal -> eTty,
     TtzLocal -> eTtz,
     TmnLocal[li,lj] -> Tmn[li,lj],
     (****************************************************)
     
     (*======================================================================================*)
     (* Trace of the stress energy tensor of matter. At this point assuming this is the last *)
     (* function that adds something to Tmunu, all the matter contributions are in place.    *)
     (* Therefore the following trace is only from the matter fields.                        *)
     (* We compute it here and store it to use it later in the evolution equations           *)
     (*======================================================================================*)
     (*TmunuTraceMatterLocal -> -TttLocal/alpLocal^2 + 2/alpLocal^2 (TtxLocal betaLocal1 + TtyLocal betaLocal2 + TtzLocal betaLocal3 )
                           -betaLocal[ui] betaLocal[uj] TmnLocal[li,lj]/alpLocal^2 + gLocalInv[ui,uj] TmnLocal[li,lj],*)
     TmunuTraceMatterLocal -> 0,

     (* Add to TmunuBase variables *)
     (* STTtt, STTtj, STTij are in the Einstein Frame. Here we make sure that the scalar *)
     (* field is the last contribution to the total Tmunu and transform Tmunu to the     *)
     (* Einstein frame. In this way we minimize modification the the BSSN thorn          *)

     eTtt -> TttLocal   + STTtt  ,
     eTtx -> TtxLocal   + STTtj1 ,
     eTty -> TtyLocal   + STTtj2 ,
     eTtz -> TtzLocal   + STTtj3 ,
     eTxx -> TmnLocal11 + STTij11,
     eTxy -> TmnLocal21 + STTij21,
     eTxz -> TmnLocal31 + STTij31,
     eTyy -> TmnLocal22 + STTij22,
     eTyz -> TmnLocal32 + STTij32,
     eTzz -> TmnLocal33 + STTij33

  }
};*)

rhoCalc[fdOrder_, PD_, PDadvect_, potential_, V_] :=
{
  Name -> fnPrefix <> "_calc_rho_" <> potential <> "_" <> fdOrder,
  Schedule -> {"at CCTK_ANALYSIS as ScalarFieldRhoCalc"},
  Where -> Interior,
  ConditionalOnKeywords -> {{"fd_order", fdOrder},{"calculate_rho","yes"},{"potential_type",potential}},
  Shorthands -> shorthands,
  Equations ->
  {
     dir[ui] -> Sign[beta[ui]],
     pi -> 3.14159265358979323846,

     detg -> gDet,
     invdetg -> 1 / detg,
     gLocalInv[ua,ub] -> invdetg gDet MatrixInverse[gLocal[ua,ub]],
     betaLocal[ui] -> beta[ui],

     (* rho in the Einstein frame *)
     rhoE -> 0.5 STpi^2 + 0.5 gLocalInv[ui,uj] PD[STphi,li] PD[STphi,lj] + V[STphi],

     (* rho in the Jordan frame *)
     ConfFactorA -> STA[STphi],
     rhoJ -> ConfFactorA^(-4) rhoE

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
  Equations ->
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
    dot[STpi]  -> -(STpi - piBG)/r + n[uk] PDonesided2nd[STpi,lk],
    dot[STphi] -> -(STphi-phiBG)/r + n[uk] PDonesided2nd[STphi,lk],

    (*dot[dtSTphi] -> 0,
    dot[STphid[li]] -> 0*)

    (* clear rho in boundaries to get rid of poison *)
    rhoE -> 0, rhoJ -> 0

  }
};

boundaryPhiOfTCalc =
{
  Name -> fnPrefix <> "_boundaryPhiOfT",
  Schedule -> {"in MoL_RHSBoundaries as ScalarFieldBoundaries"},
  ConditionalOnKeyword -> {"boundary_condition", "phi_of_t"},
  Where -> Boundary,
  Shorthands -> shorthands,
  Equations ->
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
    dot[STpi]  -> -(STpi - piBG  )/r + n[uk] PDonesided2nd[STpi,lk],
    dot[STphi] -> -(STphi- piBG*t)/r + n[uk] PDonesided2nd[STphi,lk],

    (*dot[dtSTphi] -> 0,
    dot[STphid[li]] -> 0*)

    (* clear rho in boundaries to get rid of poison *)
    rhoE -> 0, rhoJ -> 0

  }
};

(**************************************************************************************)
(* Construct the thorn *)
(**************************************************************************************)



calculations =
{
  initGFsCalc,
  (*initialDataCalc1,
  initialDataCalc2,
  initialDataStaticST,
  initialDataStatSolCalc,*)
  ADMToLocalCalc,

  evolveScalarField["full4th", PDstandard4th, PDlopsided4th],
  evolveScalarField["6th", PDstandard6th, PDlopsided6th],

  AddPotentialToRHS["phiFourth",Vfourth],
  AddPotentialToRHS["phiSixth" ,Vsixth],
  AddPotentialToRHS["phiEighth",Veighth],
  AddPotentialToRHS["phiBubble",Vbubble],
  AddPotentialToRHS["phiExpand",Vexpand],

  SetTmunuTraceToZeroCalc,

  ScalarFieldTmunuCalc["full4th", PDstandard4th,PDlopsided4th,"none",Vnone],
  ScalarFieldTmunuCalc["6th", PDstandard6th,PDlopsided6th,"none",Vnone],
  ScalarFieldTmunuCalc["full4th", PDstandard4th,PDlopsided4th,"phiFourth",Vfourth],
  ScalarFieldTmunuCalc["6th", PDstandard6th,PDlopsided6th,"phiFourth",Vfourth],
  ScalarFieldTmunuCalc["full4th", PDstandard4th,PDlopsided4th,"phiSixth",Vsixth],
  ScalarFieldTmunuCalc["6th", PDstandard6th,PDlopsided6th,"phiSixth",Vsixth],
  ScalarFieldTmunuCalc["full4th", PDstandard4th,PDlopsided4th,"phiEighth",Veighth],
  ScalarFieldTmunuCalc["6th", PDstandard6th,PDlopsided6th,"phiEighth",Veighth],
  ScalarFieldTmunuCalc["full4th", PDstandard4th,PDlopsided4th,"phiBubble",Vbubble],
  ScalarFieldTmunuCalc["6th", PDstandard6th,PDlopsided6th,"phiBubble",Vbubble],
  ScalarFieldTmunuCalc["full4th", PDstandard4th,PDlopsided4th,"phiExpand",Vexpand],
  ScalarFieldTmunuCalc["6th", PDstandard6th,PDlopsided6th,"phiExpand",Vexpand],

  (*
  BransDickeTmunuCalc["full4th", PDstandard4th,PDlopsided4th,"none",Vnone],
  BransDickeTmunuCalc["6th", PDstandard6th,PDlopsided6th,"none",Vnone],
  BransDickeTmunuCalc["full4th", PDstandard4th,PDlopsided4th,"phiFourth",Vfourth],
  BransDickeTmunuCalc["6th", PDstandard6th,PDlopsided6th,"phiFourth",Vfourth],
  BransDickeTmunuCalc["full4th", PDstandard4th,PDlopsided4th,"phiSixth",Vsixth],
  BransDickeTmunuCalc["6th", PDstandard6th,PDlopsided6th,"phiSixth",Vsixth],
  BransDickeTmunuCalc["full4th", PDstandard4th,PDlopsided4th,"phiEighth",Veighth],
  BransDickeTmunuCalc["6th", PDstandard6th,PDlopsided6th,"phiEighth",Veighth],
  *)

  rhoCalc["full4th", PDstandard4th,PDlopsided4th,"none",Vnone],
  rhoCalc["6th", PDstandard6th,PDlopsided6th,"none",Vnone],
  rhoCalc["full4th", PDstandard4th,PDlopsided4th,"phiFourth",Vfourth],
  rhoCalc["6th", PDstandard6th,PDlopsided6th,"phiFourth",Vfourth],
  rhoCalc["full4th", PDstandard4th,PDlopsided4th,"phiSixth",Vsixth],
  rhoCalc["6th", PDstandard6th,PDlopsided6th,"phiSixth",Vsixth],
  rhoCalc["full4th", PDstandard4th,PDlopsided4th,"phiEighth",Veighth],
  rhoCalc["6th", PDstandard6th,PDlopsided6th,"phiEighth",Veighth],
  rhoCalc["full4th", PDstandard4th,PDlopsided4th,"phiBubble",Vbubble],
  rhoCalc["6th", PDstandard6th,PDlopsided6th,"phiBubble",Vbubble],
  rhoCalc["full4th", PDstandard4th,PDlopsided4th,"phiExpand",Vexpand],
  rhoCalc["6th", PDstandard6th,PDlopsided6th,"phiExpand",Vexpand],

  boundaryCalc,
  boundaryPhiOfTCalc,
  EinsteinToJordanCalc
};



keywordParameters = 
{
  fdOrderParam,
  boundaryParam,
  (*initialDataParam,*)
  potentialParam,
  computeScalarTmunuParam,
  computeScalarFieldTmunuParam,
  (*systemTypeParam,*)
  calcRhoParam
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
  (*ExtendedKeywordParameters -> extendedKeywordParameters,*)
  InheritedImplementations -> {"admbase","TmunuBase"}];


