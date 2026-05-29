
/* Copyright 2013 Peter Diener, Nils Dorband, Roland Haas, Ian Hinder,
Christian Ott, Denis Pollney, Thomas Radke, Christian Reisswig, Erik
Schnetter, Barry Wardell and Burkhard Zink

This file is part of Llama.

Llama is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free
Software Foundation, either version 2 of the License, or (at your
option) any later version.

Llama is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with Llama.  If not, see <http://www.gnu.org/licenses/>. */

// This file has been generated automatically.  Do not change it manually.
dJ[0][0][0]=zp*(-(xp*xp*(yp*yp))+yp*yp*(zp*zp)-2.*pow(xp,4.)+pow(yp,4.))*pow\
(xp*xp+yp*yp,-1.5)*pow(xp*xp+yp*yp+zp*zp,-2.);
dJ[0][0][1]=-(xp*yp*zp*(3.*(xp*xp)+3.*(yp*yp)+zp*zp)*pow(xp*xp+yp*yp,-1.5)*p\
ow(xp*xp+yp*yp+zp*zp,-2.));
dJ[0][0][2]=(1.*xp*(xp*xp+yp*yp-zp*zp)*pow(xp*xp+yp*yp+zp*zp,-2.))/sqrt(xp*x\
p+yp*yp);
dJ[0][1][0]=-(xp*yp*zp*(3.*(xp*xp)+3.*(yp*yp)+zp*zp)*pow(xp*xp+yp*yp,-1.5)*p\
ow(xp*xp+yp*yp+zp*zp,-2.));
dJ[0][1][1]=zp*(xp*xp*(-(yp*yp)+zp*zp)+pow(xp,4.)-2.*pow(yp,4.))*pow(xp*xp+y\
p*yp,-1.5)*pow(xp*xp+yp*yp+zp*zp,-2.);
dJ[0][1][2]=(1.*yp*(xp*xp+yp*yp-zp*zp)*pow(xp*xp+yp*yp+zp*zp,-2.))/sqrt(xp*x\
p+yp*yp);
dJ[0][2][0]=(1.*xp*(xp*xp+yp*yp-zp*zp)*pow(xp*xp+yp*yp+zp*zp,-2.))/sqrt(xp*x\
p+yp*yp);
dJ[0][2][1]=(1.*yp*(xp*xp+yp*yp-zp*zp)*pow(xp*xp+yp*yp+zp*zp,-2.))/sqrt(xp*x\
p+yp*yp);
dJ[0][2][2]=2.*zp*pow(xp*xp+yp*yp+zp*zp,-2.)*sqrt(xp*xp+yp*yp);
dJ[1][0][0]=2.*xp*yp*pow(xp*xp+yp*yp,-2.);
dJ[1][0][1]=(-(xp*xp)+yp*yp)*pow(xp*xp+yp*yp,-2.);
dJ[1][0][2]=0.;
dJ[1][1][0]=(-(xp*xp)+yp*yp)*pow(xp*xp+yp*yp,-2.);
dJ[1][1][1]=-2.*xp*yp*pow(xp*xp+yp*yp,-2.);
dJ[1][1][2]=0.;
dJ[1][2][0]=0.;
dJ[1][2][1]=0.;
dJ[1][2][2]=0.;
dJ[2][0][0]=(yp*yp+zp*zp)*pow(xp*xp+yp*yp+zp*zp,-1.5);
dJ[2][0][1]=-(xp*yp*pow(xp*xp+yp*yp+zp*zp,-1.5));
dJ[2][0][2]=-(xp*zp*pow(xp*xp+yp*yp+zp*zp,-1.5));
dJ[2][1][0]=-(xp*yp*pow(xp*xp+yp*yp+zp*zp,-1.5));
dJ[2][1][1]=(xp*xp+zp*zp)*pow(xp*xp+yp*yp+zp*zp,-1.5);
dJ[2][1][2]=-(yp*zp*pow(xp*xp+yp*yp+zp*zp,-1.5));
dJ[2][2][0]=-(xp*zp*pow(xp*xp+yp*yp+zp*zp,-1.5));
dJ[2][2][1]=-(yp*zp*pow(xp*xp+yp*yp+zp*zp,-1.5));
dJ[2][2][2]=(xp*xp+yp*yp)*pow(xp*xp+yp*yp+zp*zp,-1.5);
