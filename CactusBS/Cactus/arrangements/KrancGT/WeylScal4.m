$Path = Join[$Path, Map[ToString, ReadList[StringToStream[Environment["MATHPATH"]], Word, WordSeparators->{":"}]]];

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
  PDstandard2nd[i_, j_] -> StandardCenteredDifferenceOperator[1,1,i] *
    StandardCenteredDifferenceOperator[1,1,j],

  PDstandard4th[i_] -> StandardCenteredDifferenceOperator[1,2,i],
  PDstandard4th[i_, i_] -> StandardCenteredDifferenceOperator[2,2,i],
  PDstandard4th[i_, j_] -> StandardCenteredDifferenceOperator[1,2,i] *
    StandardCenteredDifferenceOperator[1,2,j],

  PDplus[i_] -> DPlus[i],
  PDminus[i_] -> DMinus[i],
  PDplus[i_,j_] -> DPlus[i] DPlus[j],

  PDonesided2nd[1] -> dir[1] (-shift[1]^(2 dir[1]) + 4 shift[1]^dir[1] - 3 )/
    (2 spacing[1]),
  PDonesided2nd[2] -> dir[2] (-shift[2]^(2 dir[2]) + 4 shift[2]^dir[2] - 3 )/
    (2 spacing[2]),
  PDonesided2nd[3] -> dir[3] (-shift[3]^(2 dir[3]) + 4 shift[3]^dir[3] - 3 )/
    (2 spacing[3])
};

(*PD = PDstandard2nd;*)

(****************************************************************************
  Tensors 
 ****************************************************************************)

(* Register all the tensors that will be used with TensorTools *)
Map[DefineTensor, 
{
  R, gamma,g, gInv, k, ltet, n, rm, im, rmbar, imbar, tn, va, vb, vc,
  wa, wb, wc, ea, eb, ec,
  Ro, Rojo, R4p, Psi0r, Psi0i, Psi1r, Psi1i, Psi2r, Psi2i, Psi3r,
  Psi3i, Psi4r,Psi4i
}];

(* Psi0,2,4 behave as (pseudo)scalars *)
SetTensorAttribute[Psi0r, TensorParity, +1];
SetTensorAttribute[Psi2r, TensorParity, +1];
SetTensorAttribute[Psi4r, TensorParity, +1];
SetTensorAttribute[Psi0i, TensorParity, -1];
SetTensorAttribute[Psi2i, TensorParity, -1];
SetTensorAttribute[Psi4i, TensorParity, -1];
(* these are 'special' and require a patch to the symmetry thorns *)
SetTensorAttribute[Psi1r, TensorManualCartesianParities, {1,1,-1}];
SetTensorAttribute[Psi3r, TensorManualCartesianParities, {-1,-1,1}];
SetTensorAttribute[Psi1i, TensorManualCartesianParities, {-1,-1,1}];
SetTensorAttribute[Psi3i, TensorManualCartesianParities, {-1,-1,1}];

(* Register the TensorTools symmetries (this is very simplistic) *)
Map[AssertSymmetricDecreasing, 
{
  k[la,lb], g[la,lb]
}];

AssertSymmetricDecreasing[gamma[ua,lb,lc], lb, lc];

(* Determinants of the metrics in terms of their components
   (Mathematica symbolic expressions) *)
gDet = Det[MatrixOfComponents[g[la,lb]]];

(****************************************************************************
  Groups
 ****************************************************************************)

(* Cactus group definitions *)

scalars = { (*Psi0r, Psi0i, Psi1r, Psi1i, Psi2r, Psi2i, Psi3r, Psi3i, *) Psi4r,Psi4i};

scalarGroups = Map[CreateGroupFromTensor, scalars];

(* We need extra timelevels so that interpolation onto the extraction sphere
   works properly with mesh refinement *)
scalarGroups = Map[AddGroupExtra[#, Timelevels -> 3] &, scalarGroups];

admGroups = 
  {{"admbase::metric", {gxx,gxy,gxz,gyy,gyz,gzz}},
   {"admbase::curv", {kxx,kxy,kxz,kyy,kyz,kzz}},
   {"admbase::lapse", {alp}},
   {"admbase::shift", {betax,betay,betaz}}};

declaredGroups = scalarGroups;
declaredGroupNames = Map[First, declaredGroups];

groups = Join[declaredGroups, admGroups];

(****************************************************************************
  Shorthands
 ****************************************************************************)

shorthands = 
{
  gamma[ua,lb,lc], R[la,lb,lc,ld], invdetg, detg, third, detgmthirdroot,
  gInv[ua,ub], Ro[la,lb,lc], Rojo[la,lb], R4p[li,lj,lk,ll],
  omega11, omega22, omega33, omega12, omega13, omega23, va[ua], vb[ua], vc[ua],
  wa[ua], wb[ua], wc[ua], ea[ua], eb[ua], ec[ua],
  tn[ua], nn, ltet[ua],n[ua],rm[ua],im[ua],rmbar[ua],imbar[ua], isqrt2, xmoved,
  ymoved, zmoved
};

k11=kxx; k21=kxy; k22=kyy; k31=kxz; k32=kyz; k33=kzz;
g11=gxx; g21=gxy; g22=gyy; g31=gxz; g32=gyz; g33=gzz;

realParameters = {{Name -> offset, Default -> 10^(-15)},xorig,yorig,zorig};

PsisCalc[fdOrder_, PD_] := 
{
  Name -> "psis_calc_" <> fdOrder,
  Where -> Interior,
  ConditionalOnKeyword -> {"fd_order", fdOrder},
  Shorthands -> shorthands,
  Equations ->
  {
        detg -> gDet,
        invdetg -> 1 / detg,
        gInv[ua,ub] -> invdetg gDet MatrixInverse[g[ua,ub]],
        gamma[ua, lb, lc] -> 1/2 gInv[ua,ud] (PD[g[lb,ld], lc] + 
          PD[g[lc,ld], lb] - PD[g[lb,lc],ld]),

(****************************************************************************
  Offset the origin
 ****************************************************************************)

	xmoved -> x - xorig,
	ymoved -> y - yorig,
	zmoved -> z - zorig,

(****************************************************************************
  Compute the local tetrad
 ****************************************************************************)

	va1 -> -ymoved, va2 -> xmoved+offset, va3 -> 0,
	vb1 -> xmoved+offset, vb2 -> ymoved, vb3 -> zmoved,
	vc[ua] -> Sqrt[detg] gInv[ua,ud] Eps[ld,lb,lc] va[ub] vb[uc],

	(* Orthonormalize using Gram-Schmidt*)

        wa[ua] -> va[ua],
        omega11 -> wa[ua] wa[ub] g[la,lb],
        ea[ua] -> wa[ua] / Sqrt[omega11],

        omega12 -> ea[ua] vb[ub] g[la,lb],
        wb[ua] -> vb[ua] - omega12 ea[ua],
        omega22 -> wb[ua] wb[ub] g[la,lb],
        eb[ua] -> wb[ua]/Sqrt[omega22],

        omega13 -> ea[ua] vc[ub] g[la,lb],
        omega23 -> eb[ua] vc[ub] g[la,lb],
        wc[ua] -> vc[ua] - omega13 ea[ua] - omega23 eb[ua],
        omega33 -> wc[ua] wc[ub] g[la,lb],
        ec[ua] -> wc[ua]/Sqrt[omega33],

	(* Create Spatial Portion of Null Tetrad *)
	isqrt2  ->  0.7071067811865475244,
	ltet[ua] -> isqrt2 eb[ua],
	n[ua] -> - isqrt2 eb[ua],
	rm[ua] -> isqrt2 ec[ua],
	im[ua] -> isqrt2 ea[ua],
	rmbar[ua] -> isqrt2 ec[ua],
	imbar[ua] -> -isqrt2 ea[ua],

	(* nn here is the projection of both l^a and n^a with u^a (the time-like unit
	   vector normal to the hypersurface).  We do NOT save the t component of the 
	   tetrads in this code to avoid unnecessary factors of lapse and shift. *)
	nn -> isqrt2,


(****************************************************************************
  Compute the NP pseudoscalars
 ****************************************************************************)

	(* Calculate the relevant Riemann Quantities *)
		
        (* The 3-Riemann *)                
	R[la,lb,lc,ld] -> 1/2 ( PD[g[la,ld],lc,lb] + PD[g[lb,lc],ld,la] )
			- 1/2 ( PD[g[la,lc],lb,ld] + PD[g[lb,ld],la,lc] )
			+ g[lj,le] gamma[uj,lb,lc] gamma[ue,la,ld] 	
			- g[lj,le] gamma[uj,lb,ld] gamma[ue,la,lc],

        (* The 4-Riemann projected into the slice on all its indices. 
           The Gauss equation. *)
        R4p[li,lj,lk,ll] -> R[li,lj,lk,ll] + 
          2 AntiSymmetrize[k[li,lk] k[ll,lj], lk, ll],

        (* The 4-Riemann projected in the unit normal direction on one 
           index, then into the slice on the remaining indices. The Codazzi
           equation. *)
	Ro[lj,lk,ll] ->  - 2 AntiSymmetrize[ PD[k[lj,lk],ll], lk,ll]
			 - 2 AntiSymmetrize[ gamma[up,lj,lk] k[ll,lp], lk,ll],

        (* The 4-Riemann projected in the unit normal direction on two
           indices, and into the slice on the remaining two. *)
	Rojo[lj,ll] ->  gInv[uc,ud] (R[lj,lc,ll,ld] )   
                        - k[lj,lp] gInv[up,ud] k[ld,ll]
			+ k[lc,ld] gInv[uc,ud] k[lj,ll],

	(* Calculate End Quantities
		NOTE: In writing this, I assume m[0]=0!! to save lots of work *)

	Psi4r -> R4p[li,lj,lk,ll] n[ui] n[uk] 
                   ( rmbar[uj] rmbar[ul] - imbar[uj] imbar[ul] )
		 + 2 Ro[lj,lk,ll] n[uk] nn 
                   ( rmbar[uj] rmbar[ul] - imbar[uj] imbar[ul] )
		 + Rojo[lj,ll] nn nn ( rmbar[uj] rmbar[ul] - imbar[uj] imbar[ul] 
              (* + terms in mbar^0 == 0*) ),

	Psi4i -> R4p[la,lb,lc,ld] n[ua] n[uc] ( - rm[ub] im[ud] - im[ub] rm[ud] )
		 + 2 Ro[la,lb,lc] n[ub] nn ( - rm[ua] im[uc] - im[ua] rm[uc] )
		 + Rojo[la,lb] nn nn ( - rm[ua] im[ub] - im[ua] rm[ub] )  (*,


	Psi3r -> R4p[la,lb,lc,ld] ltet[ua] n[ub] rm[uc] n[ud]
		 + Ro[la,lb,lc] ( nn (n[ua]-ltet[ua]) rm[ub] n[uc]
                                      - nn rm[ua] ltet[ub] n[uc] )
		 - Rojo[la,lb] nn (n[ua]-ltet[ua]) nn rm[ub],

	Psi3i -> - R4p[la,lb,lc,ld] ltet[ua] n[ub] im[uc] n[ud]
		 - Ro[la,lb,lc] ( nn (n[ua]-ltet[ua]) im[ub] n[uc] 
                                  - nn im[ua] ltet[ub] n[uc] )
		 + Rojo[la,lb] nn (n[ua]-ltet[ua]) nn im[ub],

	Psi2r -> R4p[la,lb,lc,ld] ltet[ua] n[ud] (rm[ub] rm[uc] + im[ub] im[uc])
		 + Ro[la,lb,lc] nn ( n[uc] (rm[ua] rm[ub] + im[ua] im[ub]) 
                                     - ltet[ub] (rm[ua] rm[uc] + im[ua] im[uc]) )
		 - Rojo[la,lb] nn nn (rm[ua] rm[ub] + im[ua] im[ub]),

	Psi2i -> R4p[la,lb,lc,ld] ltet[ua] n[ud] (im[ub] rm[uc] - rm[ub] im[uc])
		 + Ro[la,lb,lc] nn ( n[uc] (im[ua] rm[ub] - rm[ua] im[ub]) 
                                     - ltet[ub] (rm[ua] im[uc] - im[ua] rm[uc]) )
		 - Rojo[la,lb] nn nn (im[ua] rm[ub] - rm[ua] im[ub]),

	Psi1r -> R4p[la,lb,lc,ld] n[ua] ltet[ub] rm[uc] ltet[ud]
		 + Ro[la,lb,lc] ( nn ltet[ua] rm[ub] ltet[uc] 
                                  - nn rm[ua] n[ub] ltet[uc] 
                                  - nn n[ua] rm[ub] ltet[uc] )
		 + Rojo[la,lb] nn nn ( n[ua] rm[ub] - ltet[ua] rm[ub] ),

	Psi1i -> R4p[la,lb,lc,ld] n[ua] ltet[ub] im[uc] ltet[ud]
		 + Ro[la,lb,lc] ( nn ltet[ua] im[ub] ltet[uc] 
                                  - nn im[ua] n[ub] ltet[uc] 
                                  - nn n[ua] im[ub] ltet[uc] )
		 + Rojo[la,lb] nn nn ( n[ua] im[ub] - ltet[ua] im[ub] ),

	Psi0r -> R4p[la,lb,lc,ld] ltet[ua] ltet[uc] (rm[ub] rm[ud] - im[ub] im[ud])
		 + 2 Ro[la,lb,lc] nn ltet[ub] (rm[ua] rm[uc] - im[ua] im[uc])
		 + Rojo[la,lb] nn nn (rm[ua] rm[ub] - im[ua] im[ub]),

	Psi0i -> R4p[la,lb,lc,ld] ltet[ua] ltet[uc] (rm[ub] im[ud] + im[ub] rm[ud])
		 + 2 Ro[la,lb,lc] nn ltet[ub] (rm[ua] im[uc] + im[ua] rm[uc])
		 + Rojo[la,lb] nn nn (rm[ua] im[ub] + im[ua] rm[ub]) *)
  }
};


(****************************************************************************
  Construct the thorn
 ****************************************************************************)

fdOrderParam = 
{
  Name -> "fd_order",
  Default -> "2nd",
  AllowedValues -> {"2nd", "4th"}
};

keywordParameters = 
{
  fdOrderParam
};

calculations = 
{
  PsisCalc["2nd", PDstandard2nd],
  PsisCalc["4th", PDstandard4th]
};

CreateKrancThornTT[groups, ".", "WeylScal4", 
  Calculations -> calculations,
  DeclaredGroups -> declaredGroupNames,
  PartialDerivatives -> derivatives,
  KeywordParameters -> keywordParameters,
  RealParameters -> realParameters,
  UseLoopControl -> True,
  InheritedImplementations -> {"admbase", "methodoflines"}];
