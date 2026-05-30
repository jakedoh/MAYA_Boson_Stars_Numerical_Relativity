$Path = Union[$Path,{"~/SetTensor"}];
Needs["SetTensor`"];

Dimension = 3;
x[1] = eta; x[2] = q; x[3] = phi
qf[eta_,q_,phi_] := amp (Exp[-(eta-eta0)^2/sigma^2]+Exp[-(eta+eta0)^2/sigma^2]) Sin[q]^n (1+c Cos[phi]^2)

md = {
{Exp[2 qf[eta,q,phi]],0,0},
{0,Exp[2 qf[eta,q,phi]],0},
{0,0,Sin[q]^2}} psisph[eta,q,phi]^4;
InitializeMetric[md];

Clear[exc];
DefineTensor[exc];
SetTensor[exc[la,lb],{{0,0,0},{0,0,0},{0,0,0}}];

tmp = RicciR[la,lb] Metricg[ua,ub]+exc[la,lb] Metricg[ua,ub] exc[lc,ld] Metricg[uc,ud]-
  exc[la,lb] exc[lc,ld] Metricg[ua,uc] Metricg[ub,ud];
tmp = RicciToAffine[tmp];
tmp = EvalMT[tmp];
tmp = ExpandAll[-Exp[2 qf[eta,q,phi]] psisph[eta,q,phi]^5/8 tmp]
sav=tmp
tmp = SubFun[sav,psisph[eta,q,phi],2 Cosh[eta/2]+psisph[eta,q,phi]]

(* Make the stencil... *)

stencil = ExpandAll[tmp /. {
 D[psisph[eta,q,phi],eta]->(psisph[i+1,j,k]-psisph[i-1,j,k])/(2 deta),
 D[psisph[eta,q,phi],eta,eta]->(psisph[i+1,j,k]+psisph[i-1,j,k]-2 psisph[i,j,k])/(deta deta),
 D[psisph[eta,q,phi],q]->(psisph[i,j+1,k]-psisph[i,j-1,k])/(2 dq),
 D[psisph[eta,q,phi],q,q]->(psisph[i,j+1,k]+psisph[i,j-1,k]-2 psisph[i,j,k])/(dq dq),
 D[psisph[eta,q,phi],phi]->(psisph[i,j,k+1]-psisph[i,j,k-1])/(2 dphi),
 D[psisph[eta,q,phi],phi,phi]->(psisph[i,j,k+1]+psisph[i,j,k-1]-2 psisph[i,j,k])/(dphi dphi),
 psisph[eta,q,phi]->psisph[i,j,k]
 }];

ac = Coefficient[stencil,psisph[i,j,k]]
an = Coefficient[stencil,psisph[i+1,j,k]]
as = Coefficient[stencil,psisph[i-1,j,k]]
ae = Coefficient[stencil,psisph[i,j,k+1]]
aw = Coefficient[stencil,psisph[i,j,k-1]]
aq = Coefficient[stencil,psisph[i,j+1,k]]
ab = Coefficient[stencil,psisph[i,j-1,k]]
rhs = -SubFun[tmp,psisph[eta,q,phi],0]

FortranOutputOfDepList = "(i,j,k)";
$FortranReplace = Union[{
      "UND"->"_",
      "(eta,q,phi)"->"(i,j,k)"
}];
fd = FortranOpen["bhbrill3d.x"];
FortranWrite[fd,An[i,j,k],an ];
FortranWrite[fd,As[i,j,k],as ];
FortranWrite[fd,Ae[i,j,k],ae ];
FortranWrite[fd,Aw[i,j,k],aw ];
FortranWrite[fd,Aq[i,j,k],aq ];
FortranWrite[fd,Ab[i,j,k],ab ];
FortranWrite[fd,Ac[i,j,k],ac ];
FortranWrite[fd,Rhs[i,j,k],rhs ];
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

gijtmp = Exp[2 eta]/psisph[eta,q,phi]^4 Metricg[lc,ld] mcti[uc,la] mcti[ud,lb]

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
  "(eta,q,phi)"->""
};

fd = FortranOpen["psi_1st_deriv.x"];
For[ii=1,ii<=3,ii++,
  v1 = {x,y,z}[[ii]];
  psv =ToExpression["psi"<>ToString[v1]<>"[i,j,k]"];
  rhs = CD[Exp[-eta/2] psi3d[eta,q,phi],lc] mcti[uc,la];
  rhs = EvalMT[rhs,{la->-ii}]/(Exp[-eta/2] psi3d[eta,q,phi]);
  FortranWrite[fd,psv,rhs];
];
FortranClose[fd];

fd = FortranOpen["psi_2nd_deriv.x"];
Iter[i2[ua,ub],
  v1 = {x,y,z}[[ua]];
  v2 = {x,y,z}[[ub]];
  psv = ToExpression["psi"<>ToString[v1]<>ToString[v2]<>"[i,j,k]"];
  rhs = OD[OD[Exp[-eta/2] psi3d[eta,q,phi],lc] mcti[uc,la],ld] mcti[ud,lb];
  rhs = EvalMT[rhs,{la->-ua,lb->-ub}]/(Exp[-eta/2] psi3d[eta,q,phi]);
  FortranWrite[fd,psv,rhs];
];
FortranClose[fd];
