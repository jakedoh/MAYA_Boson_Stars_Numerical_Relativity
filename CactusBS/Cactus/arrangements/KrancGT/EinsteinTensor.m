$RecursionLimit = 1000;

Get["KrancThorn`"];

(*SetDebugLevel[InfoFull];*)

SetEnhancedTimes[False];
SetSourceLanguage["C"];

(****************************************************************************
  Derivatives
 ****************************************************************************)

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

(****************************************************************************
  Tensors 
 ****************************************************************************)

(* Register all the tensors that will be used with TensorTools *)
Map[DefineTensor, 
{
  et00, et0j, etij, (* Einstein Tensor *)
  g00, g0, g, gInv, k, hgamma, beta, dtbeta, betal, dtk, dt2beta, dt2alp, dtchi, dtTrK, dtcA, dtGam, dir,
  dtgij, dtg0j, dkg00, dkg0j, dkgij, (* 1st derivs of 4-metric *)
  dtgInvij, dtgInv0j, dtgInv00, dkgInv00, dkgInv0j, dkgInvij, dkhgamma, (* 1st deriv of gInv *)
  dldkg00, dldkg0j, dldkgij, (* 2nd deriv of 4-metric *)
  Gam, rhoSAijBSSN, h, hInv, Dphi, DDphi, chgamma, Rh, RGam, Rphi, ThreeRicciBSSN, A, (* BSSN *)
  DDalpha,
  rhoPlusS, ThreeRicci, Sij, Si, S, rho, Tj, jBSSN, MomConstraint, tnn, tni, tij,
  rhoSmatter, rhoPlusSmatter, jmatter, scrG, MscrG, TrMscrG, SijTF, rhat, SrrTF
}];

(* Register the TensorTools symmetries (this is very simplistic) *)
Map[AssertSymmetricIncreasing,
{
  gInv[ua,ub], dtgInvij[ua,ub], hInv[ua,ub]
}];

Map[AssertSymmetricDecreasing, 
{
  k[la,lb], g[la,lb], etij[la,lb], dtgij[la,lb],
  dtk[la,lb], dtcA[la,lb], dtGam[ua],
  ThreeRicci[la,lb], Sij[la,lb], A[la,lb],
  h[la,lb], Rh[la,lb], RGam[la,lb], Rphi[la,lb], ThreeRicciBSSN[la,lb],
  rhoSAijBSSN[la,lb], rhoSmatter[la,lb], SijTF[la,lb], MscrG[la,lb], tij[la,lb]
}];

AssertSymmetricDecreasing[hgamma[ua,lb,lc], lb, lc];
(* Symmetries: 1st derivatives *)
AssertSymmetricDecreasing[dkgij[li,lj,lk], li, lj];
AssertSymmetricIncreasing[dkgInvij[ui,uj,lk], ui, uj];
AssertSymmetricDecreasing[dkhgamma[um,li,lj,lk], li, lj];
(* Symmetries: 2nd derivatives *)
AssertSymmetricIncreasing[dldkgij[li,lj,ll,lk], li, lj];
AssertSymmetricIncreasing[dldkgij[li,lj,ll,lk], lk, ll];
AssertSymmetricIncreasing[dldkg0j[lj,ll,lk], lk, ll];
AssertSymmetricIncreasing[dldkg00[ll,lk], lk, ll];
(* BSSN *)
AssertSymmetricDecreasing[chgamma[ui,lj,lk], lj, lk];

(* Determinants of the metrics in terms of their components
   (Mathematica symbolic expressions) *)
gDet = Det[MatrixOfComponents[g[la,lb]]];
hDet = Det[MatrixOfComponents[h[la,lb]]];
INITVALUE = "42";

(****************************************************************************
  Groups
 ****************************************************************************)

(* Cactus group definitions *)

einsteinTensorGroups = Map[CreateGroupFromTensor,{et00,et0j[la],etij[la,lb]}];
einsteinTensorGroups = Map[AddGroupExtra[#, Timelevels -> 3] &, einsteinTensorGroups];

dtBSSNGroups = { {"dtBSSN",{dtTrK}}, CreateGroupFromTensor[dtcA[la,lb]], CreateGroupFromTensor[dtGam[ua]] };

(* with mattercheck *)
matterGroups = Map[CreateGroupFromTensor,{rhoPlusSmatter,rhoSmatter[la,lb],jmatter[ua]}];

(* For other analysis *)
extraGroups = Map[ CreateGroupFromTensor,{SrrTF} ];

admGroups = 
  {{"admbase::metric", {gxx,gxy,gxz,gyy,gyz,gzz}},
   {"admbase::curv", {kxx,kxy,kxz,kyy,kyz,kzz}},
   {"admbase::lapse", {alp}},
   {"admbase::shift", {betax,betay,betaz}}};

bssnGroups =
  {{"Kranc2BSSNChiMatter::Gam_group", {Gam1, Gam2, Gam3}},
   {"Kranc2BSSNChiMatter::chi_group", {chi}},
   {"Kranc2BSSNChiMatter::K_group", {K}},
   {"Kranc2BSSNChiMatter::h_group", {h11,h21,h31,h22,h32,h33}},
   {"Kranc2BSSNChiMatter::A_group", {A11,A21,A31,A22,A32,A33}},
   {"Kranc2BSSNChiMatter::bssnCcons_group", {bssnCcons1,bssnCcons2,bssnCcons3}}};

stressenergyGroups =
  {{"TmunuBase::stress_energy_scalar", {eTtt} },
   {"TmunuBase::stress_energy_vector", {eTtx,eTty,eTtz} },
   {"TmunuBase::stress_energy_tensor", {eTxx,eTxy,eTxz,eTyy,eTyz,eTzz} } };

(* declaredGroups = Join[einsteinTensorGroups,dtBSSNGroups,rhoSGroups,matterGroups]; *)
declaredGroups = Join[einsteinTensorGroups,dtBSSNGroups,matterGroups,extraGroups];
declaredGroupNames = Map[First, declaredGroups];

groups = Join[declaredGroups, admGroups, bssnGroups, stressenergyGroups ];

realParameters = {
	{Name -> chiEps, Default -> 0.0001, Steerable -> Always }
};

(****************************************************************************
  Shorthands
 ****************************************************************************)

shorthands = 
{
  detg, invdetg, g00, g0[la], gInv[ua,ub], betal[la], hgamma[ua,lb,lc], dir[ui], 
  chiCutoff, deth, invdeth, chgamma[ua,lb,lc], hInv[ua,ub], DDphi[la,lb], Dphi[la], Rh[la,lb], RGam[la,lb], Rphi[la,lb], ThreeRicciBSSN[la,lb], DDalpha[la,lb], (* for BSSN work *)
  dtgij[li,lj], dtg0j[lj], dtg00, dkg00[lk], dkg0j[lj,lk], dkgij[li,lj,lk], (* 1st deriv of 4-metric *)
  dtgInvij[ui,uj], dtgInv0j[uj], dtgInv00, dkgInv00[lk], dkgInv0j[uj,lk], dkgInvij[ui,uj,lk], (* 1st deriv of gInv *)
  (*dthgamma[um,li,lj],*) dkhgamma[um,li,lj,lk],
  dldkg00[lk,ll], dldkg0j[lj,lk,ll], dldkgij[li,lj,lk,ll], (* 2nd deriv of 4-metric *)
  deltatime, Gam1TMINUS, Gam1TMINUS2, Gam2TMINUS, Gam2TMINUS2, Gam3TMINUS, Gam3TMINUS2, (* TMINUS -> Previous timelevels -> _p. Replace in post-Kranc processing. *)
  A11TMINUS, A11TMINUS2, A21TMINUS, A21TMINUS2, A31TMINUS, A31TMINUS2, A22TMINUS, A22TMINUS2, A32TMINUS, A32TMINUS2, A33TMINUS, A33TMINUS2,
  KTMINUS, KTMINUS2, chiTMINUS, chiTMINUS2,
  ThreeRicci[la,lb], Si[la], S, rho, tSijTF, MomConstraint[la], SijTF[la,lb], MscrG[la,lb], TrMscrG, rhat[ui],
  rhoPlusS, rhoSAijBSSN[li,lj], jBSSN[ui], tnn, tni[li], tij[li,lj]
};

k11=kxx; k21=kxy; k22=kyy; k31=kxz; k32=kyz; k33=kzz;
g11=gxx; g21=gxy; g22=gyy; g31=gxz; g32=gyz; g33=gzz;
beta1=betax; beta2=betay; beta3=betaz;
scrG1=bssnCcons1; scrG2=bssnCcons2; scrG3=bssnCcons3;

(* For comparisons with matter-calculated values *)
Sij11=eTxx; Sij21=eTxy; Sij31=eTxz; Sij22=eTyy; Sij32=eTyz; Sij33=eTzz;
Tj1=eTtx; Tj2=eTty; Tj3=eTtz;

einsteinInit = 
{
  Name -> "einstein_initGFs",
  Schedule -> {"at CCTK_BASEGRID"},
  Where -> Everywhere,
  Equations ->
  {

	et00 -> 0,
	et0j[li] -> 0,
	etij[li,lj] -> 0,

	dtcA[li,lj] -> 0,
	dtGam[ui] -> 0,
	dtTrK -> 0,

	SrrTF -> 0
  }
};

calcTimeDerivs =
{
  Name -> "calcTimeDerivs",
  Schedule -> {"at CCTK_EVOL after (MoL_Evolution bssn_evolved_to_adm) as EinsteinTensor_time_derivs"},
  Where -> Interior,
  Shorthands -> shorthands,
  Equations ->
  {

	(* After generating the thorn, replace <var>TMINUS by <var>_p and
	   <var>TMINUS2 by <var>_p_p since this usage is not compatible with
	   Kranc-generated thorns. *)

	(* Need dt of K, Aij, Gam *)
	deltatime -> 0,

	A11TMINUS -> 0,
	A11TMINUS2 -> 0,
	A21TMINUS -> 0,
	A21TMINUS2 -> 0,
	A31TMINUS -> 0,
	A31TMINUS2 -> 0,
	A22TMINUS -> 0,
	A22TMINUS2 -> 0,
	A32TMINUS -> 0,
	A32TMINUS2 -> 0,
	A33TMINUS -> 0,
	A33TMINUS2 -> 0,
	dtcA11 -> -( 4*A11TMINUS - 3*A11 - A11TMINUS2 ) / (2 deltatime), 
	dtcA21 -> -( 4*A21TMINUS - 3*A21 - A21TMINUS2 ) / (2 deltatime), 
	dtcA31 -> -( 4*A31TMINUS - 3*A31 - A31TMINUS2 ) / (2 deltatime), 
	dtcA22 -> -( 4*A22TMINUS - 3*A22 - A22TMINUS2 ) / (2 deltatime), 
	dtcA32 -> -( 4*A32TMINUS - 3*A32 - A32TMINUS2 ) / (2 deltatime), 
	dtcA33 -> -( 4*A33TMINUS - 3*A33 - A33TMINUS2 ) / (2 deltatime),

	Gam1TMINUS -> 0, 
	Gam1TMINUS2 -> 0, 
	Gam2TMINUS -> 0, 
	Gam2TMINUS2 -> 0, 
	Gam3TMINUS -> 0, 
	Gam3TMINUS2 -> 0, 
	dtGam1 -> -( 4*Gam1TMINUS - 3*Gam1 - Gam1TMINUS2 ) / (2 deltatime), 
	dtGam2 -> -( 4*Gam2TMINUS - 3*Gam2 - Gam2TMINUS2 ) / (2 deltatime), 
	dtGam3 -> -( 4*Gam3TMINUS - 3*Gam3 - Gam3TMINUS2 ) / (2 deltatime), 

	KTMINUS -> 0,
	KTMINUS2 -> 0,
	dtTrK -> - ( 4*KTMINUS - 3*K - KTMINUS2 ) / (2 deltatime)

  }
};

fdOrderParam = 
{
  Name -> "fd_order",
  Default -> "6th",
  AllowedValues -> {"2nd", "4th", "6th"}
};

einsteinTensorCalc[fdOrder_, PD_, PDadvect_] := 
{
  Name -> "einsteinTensor_calc_" <> fdOrder,
  Schedule -> {"in createEinsteinTensor as ET_effective_Calc"},
  Where -> Interior,
  ConditionalOnKeyword -> {"fd_order", fdOrder},
  Shorthands -> shorthands,
  Equations ->
  {

	dir[ui] -> Sign[beta[ui]],

	(* Getting the current 4-metric in order *)
        detg -> gDet,
        invdetg -> 1 / detg,
        gInv[ua,ub] -> invdetg gDet MatrixInverse[g[ua,ub]],
        hgamma[ua, lb, lc] -> 1/2 gInv[ua,ud] (PD[g[lb,ld], lc] + 
          PD[g[lc,ld], lb] - PD[g[lb,lc],ld]),

	(* Required for BSSN derivations *)
	deth -> hDet,
	invdeth -> 1 / deth,
	hInv[ua,ub] -> invdeth hDet MatrixInverse[h[ua,ub]],
	chiCutoff -> Max[chi,chiEps],
	chgamma[ua, lb, lc] -> 1/2 hInv[ua,ud] (PD[h[lb,ld], lc] + PD[h[lc,ld], lb] - PD[h[lb,lc],ld]),

	(* Evaluate the tensor density which introduces scrG^i into the evolution equations *)
	MscrG[li,lj] -> - 1/(8 Pi) ( TTSymmetrize[ h[lk,li] PD[ scrG[uk],lj ] + scrG[uk] h[li,lm] chgamma[um,lj,lk], li,lj] ) , 
	TrMscrG -> gInv[ui,uj] MscrG[li,lj],

	(* Calculate the 3-Ricci using the BSSN method which includes the evolved Gam^i   
	     In the MwM world, R(Gam^i) -> R(h) + 8 Pi MscrG				  
	                    so R(h) -> R(Gam^i) - 8 Pi MscrG				*)
	DDphi[lj,li] -> -1/(4 chiCutoff) ( PD[chi,lj,li] - chgamma[uk,li,lj] PD[chi,lk] - chiCutoff^-1 PD[chi,li] PD[chi,lj] ),
	Dphi[li] -> -1/(4 chiCutoff) PD[chi,li],
	RGam[li,lj] -> -1/2 hInv[uk,ul] PD[h[li,lj],ll,lk] + TTSymmetrize[h[lk,li] PD[Gam[uk],lj], li,lj]
		+ TTSymmetrize[h[li,lc] chgamma[uc,lj,lk] Gam[uk],li,lj]
		+ hInv[ul,us](2 TTSymmetrize[chgamma[uk,ll,li] h[lj,la] chgamma[ua,lk,ls],li,lj]
                              + chgamma[uk,li,ls] h[la,lk] chgamma[ua,ll,lj]),
	Rphi[li,lj] -> -2 DDphi[lj,li] - 2 h[li,lj] hInv[uk,ul] DDphi[lk,ll] + 4 Dphi[li] Dphi[lj]
                   - 4 h[li,lj] hInv[uk,ua] Dphi[la] Dphi[lk],
	ThreeRicciBSSN[li,lj] -> ( RGam[li,lj] - 8 Pi MscrG[li,lj] ) + Rphi[li,lj],
	DDalpha[li,lj] -> PD[alp,li,lj] - hgamma[uk,li,lj] PD[alp,lk],

	(*********************************************************************************)
	(* Calculation of Matter or Matter+scrG contributions to the evolution equations *)
	(* from the spacetime variables 						 *)
	(*********************************************************************************)

	(* Calculate (S_ij - M_ij)^TF from dt(Aij), imposing trace-free   
	      -- Should vanish for vacuum and MwM evolutions     	*)
	rhoSAijBSSN[li,lj] -> - 1/(8 Pi alp chiCutoff) ( dtcA[li,lj] (* L_t cAij *)
           - ( beta[um] PDadvect[ A[li,lj],lm] + 2 TTSymmetrize[ A[li,lm] PD[beta[um],lj], li,lj] - 2/3 A[li,lj] PD[beta[um],lm] ) ) (* L_B cAij *)
	   + 1/(8 Pi) ( ThreeRicciBSSN[li,lj] - (1/3) h[li,lj] hInv[um,uk] ThreeRicciBSSN[lm,lk] ) (* R^TF *)
	   - 1/(8 Pi alp) ( DDalpha[li,lj] - (1/3) h[li,lj] hInv[um,uk] DDalpha[lm,lk] ) (* DDalp^TF *)
	   + 1/(8 Pi chiCutoff) ( K A[li,lj] - 2 A[li,lk] hInv[uk,um] A[lm,lj] ),
	tSijTF -> chi hInv[ui,uj] rhoSAijBSSN[li,lj], (* Store for post-run check *)
	rhoSAijBSSN[li,lj] -> rhoSAijBSSN[li,lj] - (1/3) g[li,lj] tSijTF,

	(* Calculate the rho + S from dt(K)   
	      -- Should vanish for vacuum and MwM evolutions  *)
	rhoPlusS -> 1/(4 Pi alp) ( dtTrK - beta[um] PDadvect[K,lm] )
	   + chi/(4 Pi alp) ( hInv[ui,uj] PD[alp,li,lj] + ( hgamma[uk,lk,lm] hInv[um,ui] - hInv[ui,uk] gInv[um,uj] PD[g[lk,lm],lj] ) PD[alp,li] )
	   - chi^2/(4 Pi) hInv[ui,uk] hInv[uj,um] k[li,lj] k[lk,lm],

	(* Calculate j^i+(chi/16 Pi alp)*dt(scrG) from dtGam assuming standard BSSN    
	      -- Should vanish for vacuum and MwM evolutions  *)
	jBSSN[ui] -> - chi/( 16 Pi alp ) ( dtGam[ui] (* - L_t Gam *)
	   + beta[uj] PDadvect[Gam[ui],lj] - Gam[uj] PD[beta[ui],lj] + 2/3 Gam[ui] PD[beta[uj],lj] (* - L_B Gam *)
	   - 4/3 scrG[ui] PD[beta[uj],lj] 	(* Yo term *)
           - hInv[uj,uk] PD[ beta[ui], lk,lj] - 1/3 hInv[ui,uj] PD[beta[uk],lk,lj]
           + 2 hInv[ui,uk] hInv[uj,um] A[lk,lm] PD[alp,lj] )
           + chi/( 8 Pi ) ( chgamma[ui,lj,lk] hInv[uj,um] hInv[uk,ua] A[lm,la] 
               + 6 hInv[ui,uk] hInv[uj,um] A[lk,lm] Dphi[lj] - 2/3 hInv[ui,uj] PD[K,lj] ),

	(****************************************************************************************)
	(* Generation of quantities to create the Einstein tensor from the constraint equations *)
	(****************************************************************************************)

	(* Generate rho=T_nn from vacuum Hamiltonian constraint since we only calculate (rho + S) above *)
	rho -> 1/(16 Pi) ( chi hInv[ui,uj] ThreeRicciBSSN[li,lj] - A[li,lj] hInv[uj,uk] hInv[ui,um] A[lk,lm] + 2/3 K^2 ), 

	(* Generate j_i=-T_ni from the vacuum Momentum constraint *)
	MomConstraint[li] -> chi/(8 Pi) (6 hInv[uj,uk] A[lk,li] Dphi[lj]
                      + hInv[um,uk] (PD[A[lm,li],lk]- chgamma[ud,lm,lk] A[ld,li]
                                                    - chgamma[ud,li,lk] A[lm,ld])
                      -2/3 PD[K,li]),

	(* Generate Sij^TF from (S_ij-M_ij)^TF = rhoSAijBSSN *)
	SijTF[li,lj] -> (MscrG[li,lj]-1/3 g[li,lj] TrMscrG),
	(*SijTF[li,lj] -> (MscrG[li,lj]-1/3 g[li,lj] TrMscrG) + rhoSAijBSSN[li,lj],*)

	(***********************************************)
	(* Fill in Einstein Tensor with respect to n^a *)
	(***********************************************)

	(* This is G_nn = 8 Pi T_nn, direct from the Hamiltonian constraint *)
	et00 -> 8 Pi rho,

	(* This is G_nj = 8 Pi T_nj, direct from the Momentum constraint *)
	et0j[li] -> 8 Pi ( - MomConstraint[li] ),

	(* Reconstruct S_ij from rho+S, rho, and SijTF derived from dt(Aij) *)
	etij[li,lj] -> 8 Pi ( SijTF[li,lj] - 1/3 g[li,lj] rho ),
	(* etij[li,lj] -> 8 Pi ( SijTF[li,lj] - 1/3 g[li,lj] ( rhoPlusS - rho ) *)

	(**************************************)
	(* Fill in extra monitoring variables *)
	(**************************************)

	rhat1 -> x/Sqrt[x^2 + y^2 + z^2], 
	rhat2 -> y/Sqrt[x^2 + y^2 + z^2], 
	rhat3 -> z/Sqrt[x^2 + y^2 + z^2], 
	SrrTF -> SijTF[li,lj] rhat[ui] rhat[uj]

  }
};

einsteinProjectionParam = 
{
  Name -> "ProjectEinsteinOnto",
  Default -> "normal observer",
  AllowedValues -> {"normal observer", "time observer"}
};

projectEinsteinOnTCalc = 
{
  Name -> "project_Einstein_on_t",
  Schedule -> {"in createEinsteinTensor after ET_effective_Calc"},
  Where -> Interior,
  ConditionalOnKeyword -> {"ProjectEinsteinOnto","time observer"},
  Shorthands -> shorthands,
  Equations ->
  {

	(* Getting the current 4-metric in order *)
        detg -> gDet,
        invdetg -> 1 / detg,
        gInv[ua,ub] -> invdetg gDet MatrixInverse[g[ua,ub]],

	(* Required for BSSN derivations *)
	deth -> hDet,
	invdeth -> 1 / deth,
	hInv[ua,ub] -> invdeth hDet MatrixInverse[h[ua,ub]],

        (* The EinsteinTensor contains the normal observer. Store as shorthands and convert to time observer. *)
        tnn -> et00,
        tni[li] -> et0j[li],
        tij[li,lj] -> etij[li,lj],
 
        (* Refill effective stress energy tensor, in terms of t^a = alp n^a + beta^a.
           This only effects time-time and time-space components. *)
        et0j[li] -> - alp tni[li] + beta[uj] tij[li,lj], 
        et00 -> alp^2 tnn - 2 alp tni[li] beta[ui] + beta[ui] beta[uj] tij[li,lj]

  }
};

matterCheckParam = 
{
  Name -> "matterCheck",
  Default -> "no",
  AllowedValues -> {"yes", "no"}
};

matterCheckCalc =
{
  Name -> "fill_matter_from_tmunu",
  Schedule -> {"in createEinsteinTensor after project_Einstein_on_t"},
  Where -> Everywhere,
  ConditionalOnKeyword -> {"matterCheck","yes"},
  Shorthands -> shorthands,
  Equations ->
  {

	(* Getting the current 4-metric in order *)
        detg -> gDet,
        invdetg -> 1 / detg,
        gInv[ua,ub] -> invdetg gDet MatrixInverse[g[ua,ub]],

	(* Generate what it should be from current matter: from t^a to n^a basis, for comparison *)
        Si[la] -> ( - Tj[la] + beta[ub] Sij[la,lb] ) / alp, (* - g_ik n_j T^kj *)
        S -> gInv[ua,ub] Sij[la,lb],
        rho -> (1/alp^2) ( eTtt - 2 beta[ua] Tj[la] + beta[ua] beta[ub] Sij[la,lb] ),

	rhoSmatter[li,lj] -> ( Sij[li,lj] - (1/3) g[li,lj] S ),
	jmatter[ui] -> ( beta[uk] gInv[ui,uj] Sij[lj,lk] - gInv[ui,uj] Tj[lj] - beta[ui] rho ) / alp,
	rhoPlusSmatter -> rho + S,

	(* Rescale effective Einstein Tensor by 8 Pi, so it contains Tmunu instead of Gmunu *)
	etij[li,lj] -> etij[li,lj]/(8 Pi),
	et0j[li] -> et0j[li]/(8 Pi),
	et00 -> et00/(8 Pi)

  }
};

einsteinMatterInit = 
{
  Name -> "einstein_matter_initGFs",
  Schedule -> {"at CCTK_BASEGRID"},
  Where -> Everywhere,
  ConditionalOnKeyword -> { "matterCheck" ,"yes"},
  Equations ->
  {

	rhoPlusSmatter -> 0,
	rhoSmatter[li,lj] -> 0,
	jmatter[ui] -> 0
  }
};


(****************************************************************************
  Construct the thorn
 ****************************************************************************)

keywordParameters = 
{
  fdOrderParam, matterCheckParam, einsteinProjectionParam 
};

calculations = 
{
  einsteinInit,
  einsteinMatterInit,
  calcTimeDerivs,
  einsteinTensorCalc["2nd", PDstandard2nd, PDonesided2nd],
  einsteinTensorCalc["4th", PDstandard4th, PDlopsided4th],
  einsteinTensorCalc["6th", PDstandard6th, PDlopsided6th],
  projectEinsteinOnTCalc,
  matterCheckCalc 
};

CreateKrancThornTT[groups, ".", "EinsteinTensor", 
  Calculations -> calculations,
  DeclaredGroups -> declaredGroupNames,
  PartialDerivatives -> derivatives,
  KeywordParameters -> keywordParameters,
  RealParameters -> realParameters,
  UseLoopControl -> True,
  InheritedImplementations -> {"TmunuBase","ADMBase","Kranc2BSSNChiMatter"}];
