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
  n, dir, beta, kcurv, Tmn, betaD,
  g3EF, alpEF, betaEF, Gam3, Gam0ij, Gamkij,
  BDpi, BDphi, BDSi, dtBDphi,
  BDTtj, BDTij, BDphid, BDpid,
  exPhidd, exPhid, exPid,
  TmunuTraceMatterEF, g3InvEF, TmnEF, kcurvEF, g3Inv
}];

Map[AssertSymmetricIncreasing,
{
  g3InvEF[ua,ub]
}];

Map[AssertSymmetricDecreasing, 
{
  g3EF[la,lb], kcurv[la,lb], kcurvEF[la,lb], BDTij[la,lb], exPhidd[la,lb], Tmn[la,lb], Gam0ij[la,lb],
  TmnEF[la,lb]
}];


AssertSymmetricDecreasing[Gam3[ua,lb,lc], lb, lc];
AssertSymmetricDecreasing[Gamkij[ua,lb,lc], lb, lc];

(* Determinants of the metrics in terms of their components
   (Mathematica symbolic expressions) *)
gDet = Det[MatrixOfComponents[g3EF[la,lb]]];


fnPrefix = "BransDicke";
thornName = "BransDicke";

(**************************************************************************************)
(* Groups *)
(**************************************************************************************)

(* Cactus group definitions *)

evolvedGroups = Map[CreateGroupFromTensor, {BDphi, BDpi, dtBDphi, BDphid[li]}];


admGroups = 
  {{"admbase::metric", {gxx,gxy,gxz,gyy,gyz,gzz}},
   {"admbase::curv",  {kxx,kxy,kxz,kyy,kyz,kzz}},
   {"admbase::lapse", {alp}},
   {"admbase::shift", {betax,betay,betaz}}};

evaluatedGroups = Map[CreateGroupFromTensor,{g3EF[li,lj], alpEF, BDpid[li], TmunuTraceMatterEF}];

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
  n[ui], dir[ui],  detg, invdetg, 
  g3InvEF[ui,uj], betaEF[ui], betaD[li],
  TttEF, TtxEF, TtyEF, TtzEF, TmnEF[li,lj],
  betaSqr, dBDphiSqr, dotBDpi, pi,
  BDTtt, BDTtj[li], BDTij[li,lj], BDrho, BDSi[li],
  exPhi, exPi, exPhid[li], exPhidd[li,lj], exPid[li], sourcePhi, sourcePi, exPhit, exPit,
  Gam3[ua,lb,lc], Gam0ij[li,lj], Gamkij[uk,li,lj],
  ConfFactorA, ConfFactorA2, ConfFactorAm2, kcurvEF[li,lj], KtraceEF
};


(* Use the names of shift vector from ADMBase *)
beta1=betax; beta2=betay; beta3=betaz;
kcurv11=kxx; kcurv21=kxy; kcurv22=kyy; kcurv31=kxz; kcurv32=kyz; kcurv33=kzz; 

(* Use the names from TmunuBase*)
Tmn11=eTxx; Tmn21=eTxy; Tmn22=eTyy; Tmn31=eTxz; Tmn32=eTyz; Tmn33=eTzz; 




(**************************************************************************************)
(* Parameters *)
(**************************************************************************************)

realParameters = 
{
  {Name -> BDomega, Default -> 10},
  {Name -> BDBeta, Default -> 0 },
  {Name -> BDamplitude, Default -> 0},
  {Name -> r0, Default -> 0},
  {Name -> BDsigma, Default -> 1},
  {Name -> kx, Default -> 1},
  {Name -> ky, Default -> 1},
  {Name -> kz, Default -> 1},
  {Name -> x00, Default -> 0},
  {Name -> y00, Default -> 0},
  {Name -> z00, Default -> 0},
  {Name -> fieldMass, Default -> 0 }
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
    BDphi      -> INITVALUE,
    BDpi       -> INITVALUE,
    dtBDphi    -> INITVALUE,
    BDphid[li] -> INITVALUE
  }
};

(**************************************************************************************)
(* Scalar Field Initial Data *)
(**************************************************************************************)

initialDataParam = 
{
  Name -> "initialData_type",
  Default -> "Gaussian",
  AllowedValues -> {"Gaussian", "none", "TravelingPulse", "static"}
};




ep=0.0000000001;
IDtemp  = 1 + BDamplitude Exp[-(Sqrt[(x-x00)^2+(y-y00)^2+(z-z00)^2+ep]-r0)^2/BDsigma^2];
IDtempx = D[IDtemp,x];
IDtempy = D[IDtemp,y];
IDtempz = D[IDtemp,z];

initialDataCalc1 = 
{
  Name -> fnPrefix <> "_InitialDataGaussian",
  Schedule -> {"in ADMBase_PostInitial as BDInitialData"},
  ConditionalOnKeyword -> {"initialData_type", "Gaussian"},
  Shorthands -> shorthands,
  Equations ->
  {
    BDphi  -> BDamplitude Exp[-(Sqrt[(x-x00)^2+(y-y00)^2+(z-z00)^2+ep]-r0)^2/BDsigma^2],

    (* This is necesary if we want dot[BDphi] to be zero at t=0 *)
    BDpi    -> -(betax IDtempx + betay IDtempy + betaz IDtempz )/alp
   }
};

(* Initial data from the static single BH solution *)
initialDataStatSolCalc = 
{
  Name -> fnPrefix <> "_InitiDataStaticBransDicke",
  Schedule -> {"in ADMBase_PostInitial as BDInitialData"},
  ConditionalOnKeyword -> {"initialData_type", "static"},
  Shorthands -> shorthands,
  Equations ->
  {
    BDphi  -> BDamplitude/Sqrt[(x-1.168642873)^2+y^2+z^2+0.1] + BDamplitude/Sqrt[(x+1.168642873)^2+y^2+z^2+0.1],

    BDpi    -> 0
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

IDtemp2  = BDamplitude Exp[-Sqrt[(x-x00)^2+(y-y00)^2+(z-z00)^2+ep]^2/BDsigma^2];
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
    (* BDphi -> exxPhi, BDpi -> exxPi *)



    BDphi  -> BDamplitude Exp[-Sqrt[(x-x00)^2+(y-y00)^2+(z-z00)^2+ep]^2/BDsigma^2],
    (* This is necesary if we want dot[BDphi] to be zero at t=0 *)
    BDpi -> -(betax IDtempx + betay IDtempy + betaz IDtempz )/alpEF 
		   - 2 BDamplitude Exp[-Sqrt[(x-x00)^2+(y-y00)^2+(z-z00)^2+ep]^2/BDsigma^2] (x kx + y ky + z kz)/BDsigma^2/alpEF

   }
};


ADMBaseIDParam =
{
  Name -> "ADMBase::initial_data",
  Default -> {};
  AllowedValues -> {"StaticBransDickeBH"}
};

initialDataStaticBD =
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

initialDataStaticBD =
{
  Name -> fnPrefix <> "_StaticBransDickeBH",
  Schedule -> {"in ADMBase_InitialGauge"},
  ConditionalOnKeyword -> {"initial_lapse", "StaticBransDickeBH"},
  Shorthands -> shorthands,
  Equations ->
  {
    alpEF -> Sqrt[1 - 2/r (1 + 1)]
  }

};

   
ADMBaseShiftIDParam =
{
  Name -> "ADMBase::initial_shift",
  Default -> {};
  AllowedValues -> {"StaticBransDickeBH"}
};

initialDataStaticBD =
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


computeDtBDphiCalc[fdOrder_, PD_] := 
{
  Name -> fnPrefix<>"_DtBDphi_"<>fdOrder,
  Schedule -> {"AT Initial after ADMBase_PostInitial before EinsteinToJordanFrame as BDFieldDerivatives"},
	       (*"in MoL_PostStep before EinsteinToJordanFrame after bssn_evolved_to_adm as BDFieldDerivatives"}, *)
  (* "at CCTK_EVOL after MoL_Evolution before EinsteinToJordanFrame as BDFieldDerivatives"},*)
  ConditionalOnKeyword -> {"fd_order", fdOrder},
  Shorthands -> shorthands,
  Equations -> 
  {
    (* This function equation was inside evolution but it only goes through the Interior *)
    (* points. Here we go over Everywhere so there are no nans in the boundaries *)

      ConfFactorA   -> Exp[ 0.5 BDBeta BDphi^2 ],

      BDphid[li] -> PD[BDphi,li],
      dtBDphi   -> (alp/ConfFactorA) BDpi + beta[ui] BDphid[li]
  }
};

(*
"AT Initial after BDInitialData as BD_ADMToLocal", 
	       "in MoL_PostStep before BDEvolveScalarField as BD_ADMToLocal", 
	       "at CCTK_EVOL after MoL_Evolution before BDEvolveScalarField as BD_ADMToLocal"
 *)

ADMToLocalCalc = 
{
  Name -> fnPrefix <> "_ADMToLocal",
  Schedule -> {"AT Initial after ADMBase_PostInitial as BD_ADMToLocal", 
	       "in MoL_PostStep after EinsteinToJordanFrame as BD_ADMToLocal", 
	       "at CCTK_EVOL after MoL_Evolution as BD_ADMToLocal"},

  Shorthands -> shorthands,
  Equations ->
   {
     ConfFactorAm2 -> Exp[ -BDBeta BDphi^2 ],
     ConfFactorA   -> Exp[ 0.5 BDBeta BDphi^2 ],


     alpEF  -> alp / ConfFactorA,
     g3EF11 -> gxx ConfFactorAm2,
     g3EF21 -> gxy ConfFactorAm2,
     g3EF31 -> gxz ConfFactorAm2,
     g3EF22 -> gyy ConfFactorAm2,
     g3EF32 -> gyz ConfFactorAm2,
     g3EF33 -> gzz ConfFactorAm2
   }
};






evolveScalarField[fdOrder_, PD_] := 
{
  Name -> fnPrefix<>"_EvolveScalarField_"<>fdOrder,
  Schedule -> {"in MoL_CalcRHS as BDEvolveScalarField"},
  Where -> Interior,
  ConditionalOnKeyword -> {"fd_order", fdOrder},
  Shorthands -> shorthands,
  CollectList -> {g3Inv[ui,uj]},
  Equations -> 
  {
    (* Shorthands *)

    pi -> 3.14159265358979323846,
 
    (* Transform ADMBase variables to Einstein Frame metric *)
    ConfFactorA  -> Exp[ 0.5 BDBeta BDphi^2 ],
    betaD[li]  -> g3EF[li,lj] beta[uj],

    BDphid[li] -> PD[BDphi,li],
    dtBDphi   -> alpEF BDpi + beta[ui] BDphid[li],
    
    detg -> gDet,
    invdetg -> 1 / detg,
    g3InvEF[ui,uj] -> invdetg gDet MatrixInverse[g3EF[ui,uj]],
    kcurvEF[li,lj] -> kcurv[li,lj] / ConfFactorA + BDBeta BDphi / (ConfFactorA alp)*
    ( g3EF[li,lj] ConfFactorA^2 dtBDphi - betaD[lj] PD[BDphi,li] - betaD[li] PD[BDphi,lj] ), 


    
    KtraceEF -> kcurvEF[la,lb] g3InvEF[ua,ub],

 

    (* Brans-Dickie scalar field *)    

    dot[BDphi] -> alpEF BDpi + beta[ui] PD[BDphi,li],
    dot[BDpi]  -> 
    beta[ui] PD[BDpi,li]
    + alpEF KtraceEF BDpi
    + alpEF g3InvEF[ui,uj] PD[BDphi,li,lj] 
    - 0.5 alpEF g3InvEF[ui,uj] g3InvEF[ul,um] PD[g3EF[ll,lm],lj] PD[BDphi,li] 
    + g3InvEF[ui,uj] PD[alpEF,lj] PD[BDphi,li] 
     + 4 alpEF pi BDBeta BDphi  TmunuTraceMatterEF 

        
                 (*- alpEF fieldMass^2 BDphi*)

    
    (*,dot[BDphi] -> BDpi,
     dot[BDpi]  -> Euc[ua,ub] PD[BDphi,la,lb] *)
    
  }
};






(* If this is set to 'no', then a function that puts zeros in  TmunuTraceMatterEF is scheduled. *)
computeScalarTmunuParam = 
{
  Name -> "BDTmunuEinsteinFrame",
  Default -> "yes",
  AllowedValues -> {"yes", "no"}
};

fdOrderParam = 
{
  Name -> "fd_order",
  Default -> "full4th",
  AllowedValues -> {"2nd", (* "4th", *) "full4th", "4thcentred", "6thcentred", "6th"}
};


SetTmunuTraceToZeroCalc =
{
  Name -> fnPrefix <> "_ZeroTmunuTrace",
  Schedule -> {"in AddToTmunu as ScalarFieldTmunu after Whisky_SetTmunu"},
  Where -> Everywhere,
  ConditionalOnKeyword -> {"BDTmunuEinsteinFrame", "no"},
  Shorthands -> shorthands,
  Equations -> {
    TmunuTraceMatterEF -> 0
  }
};

scalarFieldTmunuCalc[fdOrder_,PD_] := 
{
  Name -> fnPrefix <> "_ScalarFieldTmunu_" <> fdOrder,
  Schedule -> {"in AddToTmunu as ScalarFieldTmunu after Whisky_SetTmunu"},
  Where -> Interior,
  ConditionalOnKeywords -> {{"BDTmunuEinsteinFrame", "yes"},{"fd_order", fdOrder}},
  Shorthands -> shorthands,
  Equations -> 
  {

     pi -> 3.14159265358979323846,
     ConfFactorA2  -> Exp[ BDBeta BDphi^2 ],
     
     detg -> gDet,
     invdetg -> 1 / detg,
     g3InvEF[ua,ub] -> invdetg gDet MatrixInverse[g3EF[ua,ub]],
     betaEF[ui] -> beta[ui],

     (* time derivative from eq of motion *)
     dtBDphi -> alpEF BDpi + betaEF[ui] PD[BDphi,li],


     (* auxiliary variables *)
     betaSqr -> g3EF[li,lj] betaEF[ui] betaEF[uj],
     dBDphiSqr -> dtBDphi^2/alpEF^2 - 2 BDpi dtBDphi/alpEF + g3InvEF[ui,uj] PD[BDphi,li] PD[BDphi,lj],
     

  
		    
     (* Components due to the scalar field *)
     BDTtt -> dtBDphi^2 - 1/2 (-alpEF^2 + betaSqr) ( dBDphiSqr + fieldMass^2 BDphi^2),
     BDTtj[li] -> dtBDphi PD[BDphi,li] - 1/2 g3EF[li,lj] betaEF[uj] (dBDphiSqr + fieldMass^2 BDphi^2),
     BDTij[li,lj] -> PD[BDphi,li] PD[BDphi,lj] - 1/2 g3EF[li,lj] (dBDphiSqr + fieldMass^2 BDphi^2),


     (* Compute Tmunu in Einstein Frame *)
     TttEF -> eTtt ConfFactorA2,
     TtxEF -> eTtx ConfFactorA2,
     TtyEF -> eTty ConfFactorA2,
     TtzEF -> eTtz ConfFactorA2,
     TmnEF[li,lj] -> Tmn[li,lj] ConfFactorA2,
     

     (*======================================================================================*)
     (* Trace of the stress energy tensor of matter. At this point assuming this is the last *)
     (* function that adds something to Tmunu, all the matter contributions are in place.    *)
     (* Therefore the following trace is only from the matter fields.                        *)
     (* We compute it here and store it to use it later in the evolution equations           *)
     (*======================================================================================*)
     TmunuTraceMatterEF -> -TttEF/alpEF^2 + 2/alpEF^2 (TtxEF betaEF1 + TtyEF betaEF2 + TtzEF betaEF3 )
                           -betaEF[ui] betaEF[uj] TmnEF[li,lj]/alpEF^2 + g3InvEF[ui,uj] TmnEF[li,lj],

		   

     (* Add to TmunuBase variables *)
     (* BDTtt, BDTtj, BDTij are in the Einstein Frame. Here we make sure that the scalar *)
     (* field is the last contribution to the total Tmunu and transform Tmunu to the     *)
     (* Einstein frame. In this way we minimize modification the the BSSN thorn          *)

         
     eTtt -> TttEF   + BDTtt  /(4 pi),
     eTtx -> TtxEF   + BDTtj1 /(4 pi),
     eTty -> TtyEF   + BDTtj2 /(4 pi),
     eTtz -> TtzEF   + BDTtj3 /(4 pi),
     eTxx -> TmnEF11 + BDTij11/(4 pi),
     eTxy -> TmnEF21 + BDTij21/(4 pi),
     eTxz -> TmnEF31 + BDTij31/(4 pi),
     eTyy -> TmnEF22 + BDTij22/(4 pi),
     eTyz -> TmnEF32 + BDTij32/(4 pi),
     eTzz -> TmnEF33 + BDTij33/(4 pi)
     
  }

};



computeBDScalarTmunuParam = 
{
  Name -> "BDTmunuJordanFrame",
  Default -> "yes",
  AllowedValues -> {"yes", "no"}
};



(*=================================================================================*)
(* This function is obsolete now that calculations are done in the Einstein frame  *)
(*=================================================================================*)
(*
 bransDickeTmunuCal[fdOrder_,PD_] := 
{
  Name -> fnPrefix <> "_BDScalarFieldTmunu_" <> fdOrder,
  Schedule -> {"in AddToTmunu as ScalarFieldTmunu"},
  Where -> Interior, 
  ConditionalOnKeywords -> {{"BDTmunuJordanFrame", "yes"},{"fd_order", fdOrder}},
  Shorthands -> shorthands,
  Equations -> `
  {

     pi -> 3.14159265358979323846,
     detg -> gDet,
     invdetg -> 1 / detg,
     g3Inv[ua,ub] -> invdetg gDet MatrixInverse[g3[ua,ub]],
     KtraceEF -> kcurv[la,lb] g3Inv[ua,ub],
 
     Gam3[ua,lb,lc] -> 1/2 g3Inv[ua,ud] (PD[g3[lb,ld], lc] + PD[g3[lc,ld], lb] - PD[g3[lb,lc],ld]),


     (*======================================================================================*)
     (* Trace of the stress energy tensor of matter. At this point assuming this is the last *)
     (* function that adds something to Tmunu, all the matter contributions are in place.    *)
     (* Therefore the following trace is only from the matter fields.                        *)
     (*======================================================================================*)
     TmunuTraceMatter -> -eTtt/alp^2 + 2/alp^2 (eTtx beta1 + eTty beta2 + eTtz beta3 )
		         -beta[ui] beta[uj] Tmn[li,lj]/alp^2 + g3Inv[ui,uj] Tmn[li,lj],


     (* time derivative from eq of motion *)
     dotBDphi -> alp BDpi + beta[ui] PD[BDphi,li],
     dotBDpi  -> beta[ui] PD[BDpi,li] + alp KtraceEF BDpi
                 + alp g3Inv[ui,uj] PD[BDphi,li,lj]
                 - 0.5 alp g3Inv[ui,uj] g3Inv[ul,um] PD[g3[ll,lm],lj] PD[BDphi,li]
                 + g3Inv[ui,uj] PD[alp,lj] PD[BDphi,li] 
		 - alp fieldMass^2 BDphi
                 - 8 alp pi / (3 + 2 BDomega) TmunuTraceMatter,



     (* auxiliary variables *)
     betaSqr -> g3[li,lj] beta[ui] beta[uj],
     dBDphiSqr -> dotBDphi^2/alp^2 - 2 BDpi dotBDphi/alp
		      + g3Inv[ui,uj] PD[BDphi,li] PD[BDphi,lj],
    


  
     (* Christoffel symbols in terms of 3+1 quantities, 4-geometry physical metric *)
     (* First index is upper the others at lower \Gamma^a_{bc} *)
     (* Expressions from Appendix B, M. Alcubierre book *)
     Gam0ij[li,lj]    -> -kcurv[li,lj]/alp,
     Gamkij[uk,li,lj] -> Gam3[uk,li,lj] + beta[uk] kcurv[li,lj]/alp,
    
		  
		    
     (* Components due to the scalar Brans-Dicke field *)
     BDrho ->  BDomega/BDphi^2 ( BDpi^2 + 1/2 dBDphiSqr )
               + 1/(BDphi alp) * ( dotBDpi - beta[ui] PD[BDpi,li] - g3Inv[ui,uj] PD[BDphi,li] PD[alp,lj] )
	       - 1/BDphi * 8 pi TmunuTraceMatter/(3 + 2 BDomega),

     BDSi[li] -> -BDomega/BDphi^2 BDpi PD[BDphi,li]
                 -1/BDphi * ( PD[BDpi,li] + g3Inv[uj,uk] PD[BDphi,lj] kcurv[lk,li]),
 

     BDTij[li,lj] -> BDomega/BDphi^2 (PD[BDphi,li] PD[BDphi,lj] - 1/2 g3[li,lj] dBDphiSqr)
                     + 1/BDphi * ( PD[BDphi,li,lj] - Gam0ij[li,lj] dotBDphi - Gamkij[uk,li,lj] PD[BDphi,lk] )
		     - 1/BDphi * g3[li,lj] * 8 pi TmunuTraceMatter/(3 + 2 BDomega),

     (* Now we construct the coordinate components of Tmunu *)
     BDTtj[li] -> beta[uj] BDTij[li,lj] - alp BDSi[li],     
     BDTtt -> alp^2 BDrho + 2 beta[ui] BDTtj[li] - beta[ui] beta[uj] BDTij[li,lj],
		   

     (* Add the the TmunuBase variables *)
     eTtt -> eTtt/BDphi + BDTtt  /(8 pi),
     eTtx -> eTtx/BDphi + BDTtj1 /(8 pi),
     eTty -> eTty/BDphi + BDTtj2 /(8 pi),
     eTtz -> eTtz/BDphi + BDTtj3 /(8 pi),
     eTxx -> eTxx/BDphi + BDTij11/(8 pi),
     eTxy -> eTxy/BDphi + BDTij21/(8 pi),
     eTxz -> eTxz/BDphi + BDTij31/(8 pi),
     eTyy -> eTyy/BDphi + BDTij22/(8 pi),
     eTyz -> eTyz/BDphi + BDTij32/(8 pi),
     eTzz -> eTzz/BDphi + BDTij33/(8 pi)

  }

};
 *)


(**************************************************************************************)
(* Boundary conditions *)
(**************************************************************************************)

boundaryParam = 
{
  Name -> "boundary_condition",
  Default -> "radiative",
  AllowedValues -> {"none", "radiative"}
};

boundaryCalc =
{
  Name -> fnPrefix <> "_boundary",
  Schedule -> {"in MoL_RHSBoundaries as BDScalarBoundaries"},
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


    (* Brans-Dicke boundary conditions *)
    dot[BDpi]  -> -BDpi/r + n[uk] PDonesided2nd[BDpi,lk],
    dot[BDphi] -> -(BDphi-0)/r + n[uk] PDonesided2nd[BDphi,lk],


    dot[dtBDphi] -> 0,
    dot[BDphid[li]] -> 0

  }
};




(**************************************************************************************)
(* Construct the thorn *)
(**************************************************************************************)



calculations =
{
  initGFsCalc,
  initialDataCalc1,
  initialDataCalc2,
  (*initialDataDtBDphiCalc,*)
  initialDataStaticBD,
  initialDataStatSolCalc,
  ADMToLocalCalc,
  evolveScalarField["full4th", PDstandard4th],
  evolveScalarField["6th", PDstandard6th],

  computeDtBDphiCalc["full4th", PDstandard4th],
  computeDtBDphiCalc["6th", PDstandard6th],
  SetTmunuTraceToZeroCalc,
  scalarFieldTmunuCalc["full4th", PDstandard4th],
  scalarFieldTmunuCalc["6th", PDstandard6th],
  boundaryCalc
};



keywordParameters = 
{
  fdOrderParam,
  boundaryParam,
  initialDataParam,
  computeScalarTmunuParam,
  computeBDScalarTmunuParam
};

extendedKeywordParameters =
{
  ADMBaseIDParam,
  ADMBaseLapseIDParam,
  ADMBaseShiftIDParam
};

CreateKrancThornTT[groups, ".", thornName,
  Calculations -> calculations,
  DeclaredGroups -> declaredGroupNames,
  PartialDerivatives -> derivatives,
  RealParameters -> realParameters,
  KeywordParameters -> keywordParameters,
  ExtendedKeywordParameters -> extendedKeywordParameters,
  InheritedImplementations -> {"admbase","TmunuBase"}];


