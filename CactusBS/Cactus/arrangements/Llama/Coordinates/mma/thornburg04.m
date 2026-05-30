(* ::Package:: *)

(* thornburg04.m

Calculate the jacobians for the thornburg04 coordinate system.

This Mathematica script generates header files which contain the
components of the jacobians required for transforming local
thornburg04 coordinates to a global coordinate basis. The generated
header files need to be copied to the
  Coordinates/src/jacobian
directory.

This script makes use of the Format.m package, downloaded from
  http://library.wolfram.com/infocenter/MathSource/60/

*)

<< "Format.m"

WriteHeaderFile[FName_, CurvCoords_] := Module[
  {Coord, Jac, DJac, CJac, CDJac, FP},
  Coord = {xp, yp, zp};
  Jac = Simplify[Outer[D, CurvCoords, Coord]];
  DJac = Simplify[Outer[D, Jac, Coord]];
  CJac = CAssign[J, Jac, AssignOptimize->False];
  CDJac = CAssign[dJ, DJac, AssignOptimize->False];
  FP = OpenWrite[FName<>"_J.hh"];
  WriteString[FP, "// This file has been generated automatically.  Do not change it manually.\n"];
  WriteString[FP, CJac];
  WriteString[FP, "\n"];
  Close[FP];
  FP = OpenWrite[FName<>"_dJ.hh"];
  WriteString[FP, "// This file has been generated automatically.  Do not change it manually.\n"];
  WriteString[FP, CDJac];
  WriteString[FP, "\n"];
  Close[FP]
]

WriteRFile[FName_, Rtilde_] := Module[{},
  FP = OpenWrite[FName<>".hh"];
  WriteString[FP, "// This file has been generated automatically.  Do not change it manually.\n"];
  WriteString[FP, CAssign[R, Rtilde, AssignOptimize->False]];
  WriteString[FP, "\n"];
  Close[FP]
]

rp[x_, y_, z_] = Sqrt[x^2 + y^2 + z^2]


Rg1[Rl_] = ((h1+h0)/(2 h0)) (Rl-(Rmax+Rmin)/2)+((1-(h0+h1)/(2 h0))Sqrt[1+(16 (-Rmax-Rmin)^2)/(Rmax-Rmin)^2](Rmax-Rmin)^2)/(32 (-Rmax-Rmin))Sqrt[(1+((Rl-(Rmax+Rmin)/2)/((Rmax-Rmin)/8))^2)]

Rg[Rl_] = Rg1[Rl-Rstart] - Rg1[0]+Rstart

Rl[Rg_] = Simplify[Solve[Rg[Rl] == Rg, Rl]][[2]][[1]][[2]]


WriteRFile["thornburg04_local_to_global", Rg[Rl]]
WriteRFile["thornburg04_global_to_local", Rl[Rg]]

PatchXPlus = {ArcTan[zp/xp], ArcTan[yp/xp], Rl[rp[xp,yp,zp]]}
WriteHeaderFile["thornburg04_x_plus", PatchXPlus]

PatchYPlus = {ArcTan[zp/yp], ArcTan[xp/yp], Rl[rp[xp,yp,zp]]}
WriteHeaderFile["thornburg04_y_plus", PatchYPlus]

PatchZPlus = {ArcTan[yp/zp], ArcTan[xp/zp], Rl[rp[xp,yp,zp]]}
WriteHeaderFile["thornburg04_z_plus", PatchZPlus]

Print["Done."]



