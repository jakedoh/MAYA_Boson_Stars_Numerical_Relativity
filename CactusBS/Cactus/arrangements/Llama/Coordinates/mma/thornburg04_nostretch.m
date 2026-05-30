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

PatchXPlus = {ArcTan[zp/xp], ArcTan[yp/xp], rp[xp,yp,zp]}
WriteHeaderFile["thornburg04_nostretch_x_plus", PatchXPlus]

PatchYPlus = {ArcTan[zp/yp], ArcTan[xp/yp], rp[xp,yp,zp]}
WriteHeaderFile["thornburg04_nostretch_y_plus", PatchYPlus]

PatchZPlus = {ArcTan[yp/zp], ArcTan[xp/zp], rp[xp,yp,zp]}
WriteHeaderFile["thornburg04_nostretch_z_plus", PatchZPlus]

Print["Done."]
