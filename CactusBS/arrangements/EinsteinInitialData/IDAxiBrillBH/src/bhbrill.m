$Path = Union[$Path,{"~/SetTensor"}];
Needs["SetTensor`"];

Dimension = 3;
x[1] = eta; x[2] = q; x[3] = phi;
qf[eta_,q_] := amp (Exp[-(eta-eta0)^2/sigma^2]+Exp[-(eta+eta0)^2/sigma^2]) Sin[q]^n

md = {
{Exp[2 qf[eta,q]],0,0},
{0,Exp[2 qf[eta,q]],0},
{0,0,Sin[q]^2}} psi2d[eta,q]^4;
InitializeMetric[md];

Clear[exc];
DefineTensor[exc];
SetTensor[exc[la,lb],{{0,0,0},{0,0,0},{0,0,0}}];

tmp = RicciR[la,lb] Metricg[ua,ub]+exc[la,lb] Metricg[ua,ub] exc[lc,ld] Metricg[uc,ud]-
  exc[la,lb] exc[lc,ld] Metricg[ua,uc] Metricg[ub,ud];
tmp = RicciToAffine[tmp];
tmp = EvalMT[tmp];
tmp = ExpandAll[-Exp[2 qf[eta,q]] psi2d[eta,q]^5/8 tmp]
sav=tmp
tmp = SubFun[sav,psi2d[eta,q],2 Cosh[eta/2]+psi2d[eta,q]]

(* Make the stencil... *)

stencil = ExpandAll[tmp /. {
 D[psi2d[eta,q],eta]->(psi2d[i+1,j]-psi2d[i-1,j])/(2 deta),
 D[psi2d[eta,q],eta,eta]->(psi2d[i+1,j]+psi2d[i-1,j]-2 psi2d[i,j])/(deta deta),
 D[psi2d[eta,q],q]->(psi2d[i,j+1]-psi2d[i,j-1])/(2 dq),
 D[psi2d[eta,q],q,q]->(psi2d[i,j+1]+psi2d[i,j-1]-2 psi2d[i,j])/(dq dq),
 psi2d[eta,q]->psi2d[i,j]
 }];

cn = Coefficient[stencil,psi2d[i,j+1]]
cs = Coefficient[stencil,psi2d[i,j-1]]
ce = Coefficient[stencil,psi2d[i+1,j]]
cw = Coefficient[stencil,psi2d[i-1,j]]
cc = Coefficient[stencil,psi2d[i,j]]
rhs = -SubFun[tmp,psi2d[eta,q],0]

FortranOutputOfDepList = "(i,j)";
$FortranReplace = Union[{
      "UND"->"_",
      "(eta,q)"->"(i,j)"
}];
fd = FortranOpen["bhbrill.x"];
FortranWrite[fd,Cn[i,j],cn ];
FortranWrite[fd,Cs[i,j],cs ];
FortranWrite[fd,Cw[i,j],cw ];
FortranWrite[fd,Cc[i,j],cc ];
FortranWrite[fd,Ce[i,j],ce ];
FortranWrite[fd,Rhs[i,j],rhs ];
FortranClose[fd];

(* Next part, write out conformal g's and d's *)


xv = Exp[eta] Sin[q] Cos[phi];
yv = Exp[eta] Sin[q] Sin[phi];
zv = Exp[eta] Cos[q];

mc = Table[ D[ {xv,yv,zv}[[i]], {eta,q,phi}[[j]] ],{i,1,3},{j,1,3}];
mci = Simplify[Inverse[mc]];

Clear[mct];
DefineTensor[mct,{{1,2},1}];
Iter[mct[ua,lb],
  mct[ua,lb]=mc[[ua,-lb]];
];

Clear[mcti];
DefineTensor[mcti,{{1,2},1}];
Iter[mcti[ua,lb],
  mcti[ua,lb]=mci[[ua,-lb]];
];

gijtmp = Exp[2 eta]/psi2d[eta,q]^4 Metricg[lc,ld] mcti[uc,la] mcti[ud,lb]

Clear[i2];
DefineTensor[i2,{{2,1},1}];

fd = FortranOpen["gij.x"];
Iter[i2[ua,ub],
  v1 = {x,y,z}[[ua]];
  v2 = {x,y,z}[[ub]];
  metv = ToExpression["g"<>ToString[v1]<>ToString[v2]<>"[i,j,k]"];
  gg[v1,v2]=Simplify[EvalMT[gijtmp,{la->-ua,lb->-ub}]];
  FortranWrite[fd,metv,gg[v1,v2]];
  For[ii=1,ii<=3,ii++,
    v3 = {x,y,z}[[ii]];
    dmetv = ToExpression["d"<>ToString[v3]<>ToString[metv]];
    res = OD[gg[v1,v2],lc] mcti[uc,ld]/2;
    res = EvalMT[res,ld-> -ii];
    res = Simplify[res];
    FortranWrite[fd,dmetv,res];
  ];
];
FortranClose[fd];

$FortranReplace = {
  "UND"->"_",
  "(eta,q)"->""
};

fd = FortranOpen["psi_1st_deriv.x"];
For[ii=1,ii<=3,ii++,
  v1 = {x,y,z}[[ii]];
  psv =ToExpression["psi"<>ToString[v1]<>"[i,j,k]"];
  rhs = CD[Exp[-eta/2] psi2dv[eta,q],lc] mcti[uc,la];
  rhs = EvalMT[rhs,{la->-ii}]/(Exp[-eta/2] psi2dv[eta,q]);
  FortranWrite[fd,psv,rhs];
];
FortranClose[fd];

fd = FortranOpen["psi_2nd_deriv.x"];
Iter[i2[ua,ub],
  v1 = {x,y,z}[[ua]];
  v2 = {x,y,z}[[ub]];
  psv = ToExpression["psi"<>ToString[v1]<>ToString[v2]<>"[i,j,k]"];
  rhs = OD[OD[Exp[-eta/2] psi2dv[eta,q],lc] mcti[uc,la],ld] mcti[ud,lb];
  rhs = EvalMT[rhs,{la->-ua,lb->-ub}]/(Exp[-eta/2] psi2dv[eta,q]);
  FortranWrite[fd,psv,rhs];
];
FortranClose[fd];
