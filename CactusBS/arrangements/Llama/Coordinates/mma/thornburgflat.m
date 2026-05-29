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
  DJac = Outer[D, Jac, Coord];
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

R[x_,y_,z_] = Sqrt[x^2+y^2+z^2]

tgaxp[x_,y_,z_] = z/x
tgayp[x_,y_,z_] = z/y
tgazp[x_,y_,z_] = y/z

tgbxp[x_,y_,z_] = y/x
tgbyp[x_,y_,z_] = x/y
tgbzp[x_,y_,z_] = x/z

cxp[x_,y_,z_] = (ri-rm)/(R[x,y,z]*ri/x-rm) (R[x,y,z]-rm)+rm
cyp[x_,y_,z_] = (ri-rm)/(R[x,y,z]*ri/y-rm) (R[x,y,z]-rm)+rm
czp[x_,y_,z_] = (ri-rm)/(R[x,y,z]*ri/z-rm) (R[x,y,z]-rm)+rm

cxm[x_,y_,z_] = (ri-rm)/(R[x,y,z]*ri/(-x)-rm) (R[x,y,z]-rm)+rm
cym[x_,y_,z_] = (ri-rm)/(R[x,y,z]*ri/(-y)-rm) (R[x,y,z]-rm)+rm
czm[x_,y_,z_] = (ri-rm)/(R[x,y,z]*ri/(-z)-rm) (R[x,y,z]-rm)+rm


sol1xp = Solve[{tgaxp[x,y,z]==tga, tgbxp[x,y,z]==tgb, cxp[x,y,z]==c}, {x,y,z}][[1]] // Simplify
sol2xp = Solve[{tgaxp[x,y,z]==tga, tgbxp[x,y,z]==tgb, cxp[x,y,z]==c}, {x,y,z}][[2]] // Simplify
sol1xm = Solve[{tgaxp[x,y,z]==tga, tgbxp[x,y,z]==tgb, cxm[x,y,z]==c}, {x,y,z}][[1]] // Simplify
sol2xm = Solve[{tgaxp[x,y,z]==tga, tgbxp[x,y,z]==tgb, cxm[x,y,z]==c}, {x,y,z}][[2]] // Simplify
                                                                                   
sol1yp = Solve[{tgayp[x,y,z]==tga, tgbyp[x,y,z]==tgb, cyp[x,y,z]==c}, {x,y,z}][[1]] // Simplify
sol2yp = Solve[{tgayp[x,y,z]==tga, tgbyp[x,y,z]==tgb, cyp[x,y,z]==c}, {x,y,z}][[2]] // Simplify
sol1ym = Solve[{tgayp[x,y,z]==tga, tgbyp[x,y,z]==tgb, cym[x,y,z]==c}, {x,y,z}][[1]] // Simplify
sol2ym = Solve[{tgayp[x,y,z]==tga, tgbyp[x,y,z]==tgb, cym[x,y,z]==c}, {x,y,z}][[2]] // Simplify
                                                                                  
sol1zp = Solve[{tgazp[x,y,z]==tga, tgbzp[x,y,z]==tgb, czp[x,y,z]==c}, {x,y,z}][[1]] // Simplify
sol2zp = Solve[{tgazp[x,y,z]==tga, tgbzp[x,y,z]==tgb, czp[x,y,z]==c}, {x,y,z}][[2]] // Simplify
sol1zm = Solve[{tgazp[x,y,z]==tga, tgbzp[x,y,z]==tgb, czm[x,y,z]==c}, {x,y,z}][[1]] // Simplify
sol2zm = Solve[{tgazp[x,y,z]==tga, tgbzp[x,y,z]==tgb, czm[x,y,z]==c}, {x,y,z}][[2]] // Simplify

Rglobalxp[tga_,tgb_,c_] = Refine[R[x,y,z] /. sol2xp, c-ri>0] //Simplify
Rglobalxm[tga_,tgb_,c_] = Refine[R[x,y,z] /. sol1xm, c-ri>0] //Simplify
                                                           
Rglobalyp[tga_,tgb_,c_] = Refine[R[x,y,z] /. sol2yp, c-ri>0] //Simplify
Rglobalym[tga_,tgb_,c_] = Refine[R[x,y,z] /. sol1ym, c-ri>0] //Simplify
                                                          
Rglobalzp[tga_,tgb_,c_] = Refine[R[x,y,z] /. sol2zp, c-ri>0] //Simplify
Rglobalzm[tga_,tgb_,c_] = Refine[R[x,y,z] /. sol1zm, c-ri>0] //Simplify


WriteRFile["thornburg04_flat_x_plus_global_to_local", cxp[x,y,z]]
WriteRFile["thornburg04_flat_y_plus_global_to_local", cyp[x,y,z]]
WriteRFile["thornburg04_flat_z_plus_global_to_local", czp[x,y,z]]

WriteRFile["thornburg04_flat_x_minus_global_to_local", cxm[x,y,z]]
WriteRFile["thornburg04_flat_y_minus_global_to_local", cym[x,y,z]]
WriteRFile["thornburg04_flat_z_minus_global_to_local", czm[x,y,z]]


WriteRFile["thornburg04_flat_x_plus_local_to_global", Rglobalxp[xp,yp,zp]]
WriteRFile["thornburg04_flat_y_plus_local_to_global", Rglobalyp[xp,yp,zp]]
WriteRFile["thornburg04_flat_z_plus_local_to_global", Rglobalzp[xp,yp,zp]]

WriteRFile["thornburg04_flat_x_minus_local_to_global", Rglobalxm[xp,yp,zp]]
WriteRFile["thornburg04_flat_y_minus_local_to_global", Rglobalym[xp,yp,zp]]
WriteRFile["thornburg04_flat_z_minus_local_to_global", Rglobalzm[xp,yp,zp]]



PatchXPlus = {ArcTan[zp/xp], ArcTan[yp/xp], cxp[xp,yp,zp]}
WriteHeaderFile["thornburg04_flat_x_plus", PatchXPlus]

PatchXMinus = {ArcTan[zp/xp], ArcTan[yp/xp], cxm[xp,yp,zp]}
WriteHeaderFile["thornburg04_flat_x_minus", PatchXMinus]

PatchYPlus = {ArcTan[zp/yp], ArcTan[xp/yp], cyp[xp,yp,zp]}
WriteHeaderFile["thornburg04_flat_y_plus", PatchYPlus]

PatchYMinus = {ArcTan[zp/yp], ArcTan[xp/yp], cym[xp,yp,zp]}
WriteHeaderFile["thornburg04_flat_y_minus", PatchYMinus]

PatchZPlus = {ArcTan[yp/zp], ArcTan[xp/zp], czp[xp,yp,zp]}
WriteHeaderFile["thornburg04_flat_z_plus", PatchZPlus]

PatchZMinus = {ArcTan[yp/zp], ArcTan[xp/zp], czm[xp,yp,zp]}
WriteHeaderFile["thornburg04_flat_z_minus", PatchZMinus]

Print["Done."]
