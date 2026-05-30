
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
dJ[0][0][0]=0.;
dJ[0][0][1]=0.;
dJ[0][0][2]=0.;
dJ[0][1][0]=0.;
dJ[0][1][1]=2.*yp*zp*pow(yp*yp+zp*zp,-2.);
dJ[0][1][2]=-(1/(yp*yp+zp*zp))+2.*(zp*zp)*pow(yp*yp+zp*zp,-2.);
dJ[0][2][0]=0.;
dJ[0][2][1]=1/(yp*yp+zp*zp)-2.*(yp*yp)*pow(yp*yp+zp*zp,-2.);
dJ[0][2][2]=-2.*yp*zp*pow(yp*yp+zp*zp,-2.);
dJ[1][0][0]=-2.*xp*yp*pow(xp*xp+yp*yp,-2.);
dJ[1][0][1]=1/(xp*xp+yp*yp)-2.*(yp*yp)*pow(xp*xp+yp*yp,-2.);
dJ[1][0][2]=0.;
dJ[1][1][0]=-(1/(xp*xp+yp*yp))+2.*(xp*xp)*pow(xp*xp+yp*yp,-2.);
dJ[1][1][1]=2.*xp*yp*pow(xp*xp+yp*yp,-2.);
dJ[1][1][2]=0.;
dJ[1][2][0]=0.;
dJ[1][2][1]=0.;
dJ[1][2][2]=0.;
dJ[2][0][0]=(2.*ri*(ri-rm)*rm*yp*(ri+yp)*(xp*xp)*pow(rm*yp+ri*sqrt(xp*xp+yp*\
yp+zp*zp),-3.))/(xp*xp+yp*yp+zp*zp)+(ri-rm)*rm*yp*(ri+yp)*(xp*xp)*pow(xp*xp\
+yp*yp+zp*zp,-1.5)*pow(rm*yp+ri*sqrt(xp*xp+yp*yp+zp*zp),-2.)-(1.*(ri-rm)*rm\
*yp*(ri+yp)*pow(rm*yp+ri*sqrt(xp*xp+yp*yp+zp*zp),-2.))/sqrt(xp*xp+yp*yp+zp*\
zp);
dJ[2][0][1]=(ri-rm)*rm*xp*(ri+yp)*(yp*yp)*pow(xp*xp+yp*yp+zp*zp,-1.5)*pow(rm\
*yp+ri*sqrt(xp*xp+yp*yp+zp*zp),-2.)-(1.*(ri-rm)*rm*xp*yp*pow(rm*yp+ri*sqrt(\
xp*xp+yp*yp+zp*zp),-2.))/sqrt(xp*xp+yp*yp+zp*zp)-(1.*(ri-rm)*rm*xp*(ri+yp)*\
pow(rm*yp+ri*sqrt(xp*xp+yp*yp+zp*zp),-2.))/sqrt(xp*xp+yp*yp+zp*zp)+(2.*(ri-\
rm)*rm*xp*yp*(ri+yp)*pow(rm*yp+ri*sqrt(xp*xp+yp*yp+zp*zp),-3.)*(rm+(1.*ri*y\
p)/sqrt(xp*xp+yp*yp+zp*zp)))/sqrt(xp*xp+yp*yp+zp*zp);
dJ[2][0][2]=(2.*ri*(ri-rm)*rm*xp*yp*(ri+yp)*zp*pow(rm*yp+ri*sqrt(xp*xp+yp*yp\
+zp*zp),-3.))/(xp*xp+yp*yp+zp*zp)+(ri-rm)*rm*xp*yp*(ri+yp)*zp*pow(xp*xp+yp*\
yp+zp*zp,-1.5)*pow(rm*yp+ri*sqrt(xp*xp+yp*yp+zp*zp),-2.);
dJ[2][1][0]=(2.*ri*(ri-rm)*xp*(rm*pow(yp,3.)+ri*(-(rm*(xp*xp+zp*zp))+pow(xp*\
xp+yp*yp+zp*zp,1.5)))*pow(rm*yp+ri*sqrt(xp*xp+yp*yp+zp*zp),-3.))/(xp*xp+yp*\
yp+zp*zp)+(ri-rm)*xp*pow(xp*xp+yp*yp+zp*zp,-1.5)*(rm*pow(yp,3.)+ri*(-(rm*(x\
p*xp+zp*zp))+pow(xp*xp+yp*yp+zp*zp,1.5)))*pow(rm*yp+ri*sqrt(xp*xp+yp*yp+zp*\
zp),-2.)-(1.*ri*(ri-rm)*pow(rm*yp+ri*sqrt(xp*xp+yp*yp+zp*zp),-2.)*(-2.*rm*x\
p+3.*xp*sqrt(xp*xp+yp*yp+zp*zp)))/sqrt(xp*xp+yp*yp+zp*zp);
dJ[2][1][1]=(ri-rm)*yp*pow(xp*xp+yp*yp+zp*zp,-1.5)*(rm*pow(yp,3.)+ri*(-(rm*(\
xp*xp+zp*zp))+pow(xp*xp+yp*yp+zp*zp,1.5)))*pow(rm*yp+ri*sqrt(xp*xp+yp*yp+zp\
*zp),-2.)+(2.*(ri-rm)*(rm*pow(yp,3.)+ri*(-(rm*(xp*xp+zp*zp))+pow(xp*xp+yp*y\
p+zp*zp,1.5)))*pow(rm*yp+ri*sqrt(xp*xp+yp*yp+zp*zp),-3.)*(rm+(1.*ri*yp)/sqr\
t(xp*xp+yp*yp+zp*zp)))/sqrt(xp*xp+yp*yp+zp*zp)-(1.*(ri-rm)*pow(rm*yp+ri*sqr\
t(xp*xp+yp*yp+zp*zp),-2.)*(3.*rm*(yp*yp)+3.*ri*yp*sqrt(xp*xp+yp*yp+zp*zp)))\
/sqrt(xp*xp+yp*yp+zp*zp);
dJ[2][1][2]=(2.*ri*(ri-rm)*zp*(rm*pow(yp,3.)+ri*(-(rm*(xp*xp+zp*zp))+pow(xp*\
xp+yp*yp+zp*zp,1.5)))*pow(rm*yp+ri*sqrt(xp*xp+yp*yp+zp*zp),-3.))/(xp*xp+yp*\
yp+zp*zp)+(ri-rm)*zp*pow(xp*xp+yp*yp+zp*zp,-1.5)*(rm*pow(yp,3.)+ri*(-(rm*(x\
p*xp+zp*zp))+pow(xp*xp+yp*yp+zp*zp,1.5)))*pow(rm*yp+ri*sqrt(xp*xp+yp*yp+zp*\
zp),-2.)-(1.*ri*(ri-rm)*pow(rm*yp+ri*sqrt(xp*xp+yp*yp+zp*zp),-2.)*(-2.*rm*z\
p+3.*zp*sqrt(xp*xp+yp*yp+zp*zp)))/sqrt(xp*xp+yp*yp+zp*zp);
dJ[2][2][0]=(2.*ri*(ri-rm)*rm*xp*yp*(ri+yp)*zp*pow(rm*yp+ri*sqrt(xp*xp+yp*yp\
+zp*zp),-3.))/(xp*xp+yp*yp+zp*zp)+(ri-rm)*rm*xp*yp*(ri+yp)*zp*pow(xp*xp+yp*\
yp+zp*zp,-1.5)*pow(rm*yp+ri*sqrt(xp*xp+yp*yp+zp*zp),-2.);
dJ[2][2][1]=(ri-rm)*rm*(ri+yp)*zp*(yp*yp)*pow(xp*xp+yp*yp+zp*zp,-1.5)*pow(rm\
*yp+ri*sqrt(xp*xp+yp*yp+zp*zp),-2.)-(1.*(ri-rm)*rm*yp*zp*pow(rm*yp+ri*sqrt(\
xp*xp+yp*yp+zp*zp),-2.))/sqrt(xp*xp+yp*yp+zp*zp)-(1.*(ri-rm)*rm*(ri+yp)*zp*\
pow(rm*yp+ri*sqrt(xp*xp+yp*yp+zp*zp),-2.))/sqrt(xp*xp+yp*yp+zp*zp)+(2.*(ri-\
rm)*rm*yp*(ri+yp)*zp*pow(rm*yp+ri*sqrt(xp*xp+yp*yp+zp*zp),-3.)*(rm+(1.*ri*y\
p)/sqrt(xp*xp+yp*yp+zp*zp)))/sqrt(xp*xp+yp*yp+zp*zp);
dJ[2][2][2]=(2.*ri*(ri-rm)*rm*yp*(ri+yp)*(zp*zp)*pow(rm*yp+ri*sqrt(xp*xp+yp*\
yp+zp*zp),-3.))/(xp*xp+yp*yp+zp*zp)+(ri-rm)*rm*yp*(ri+yp)*(zp*zp)*pow(xp*xp\
+yp*yp+zp*zp,-1.5)*pow(rm*yp+ri*sqrt(xp*xp+yp*yp+zp*zp),-2.)-(1.*(ri-rm)*rm\
*yp*(ri+yp)*pow(rm*yp+ri*sqrt(xp*xp+yp*yp+zp*zp),-2.))/sqrt(xp*xp+yp*yp+zp*\
zp);
