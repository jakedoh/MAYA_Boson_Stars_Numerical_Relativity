
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
dJ[0][0][0]=2.*xp*zp*pow(xp*xp+zp*zp,-2.);
dJ[0][0][1]=0.;
dJ[0][0][2]=-(1/(xp*xp+zp*zp))+2.*(zp*zp)*pow(xp*xp+zp*zp,-2.);
dJ[0][1][0]=0.;
dJ[0][1][1]=0.;
dJ[0][1][2]=0.;
dJ[0][2][0]=1/(xp*xp+zp*zp)-2.*(xp*xp)*pow(xp*xp+zp*zp,-2.);
dJ[0][2][1]=0.;
dJ[0][2][2]=-2.*xp*zp*pow(xp*xp+zp*zp,-2.);
dJ[1][0][0]=2.*xp*yp*pow(xp*xp+yp*yp,-2.);
dJ[1][0][1]=-(1/(xp*xp+yp*yp))+2.*(yp*yp)*pow(xp*xp+yp*yp,-2.);
dJ[1][0][2]=0.;
dJ[1][1][0]=1/(xp*xp+yp*yp)-2.*(xp*xp)*pow(xp*xp+yp*yp,-2.);
dJ[1][1][1]=-2.*xp*yp*pow(xp*xp+yp*yp,-2.);
dJ[1][1][2]=0.;
dJ[1][2][0]=0.;
dJ[1][2][1]=0.;
dJ[1][2][2]=0.;
dJ[2][0][0]=-((ri-rm)*xp*pow(xp*xp+yp*yp+zp*zp,-1.5)*(-(rm*(ri*(yp*yp+zp*zp)\
+pow(xp,3.)))+ri*pow(xp*xp+yp*yp+zp*zp,1.5))*pow(rm*xp-ri*sqrt(xp*xp+yp*yp+\
zp*zp),-2.))-(2.*(ri-rm)*(-(rm*(ri*(yp*yp+zp*zp)+pow(xp,3.)))+ri*pow(xp*xp+\
yp*yp+zp*zp,1.5))*pow(rm*xp-ri*sqrt(xp*xp+yp*yp+zp*zp),-3.)*(rm-(1.*ri*xp)/\
sqrt(xp*xp+yp*yp+zp*zp)))/sqrt(xp*xp+yp*yp+zp*zp)+(1.*(ri-rm)*pow(rm*xp-ri*\
sqrt(xp*xp+yp*yp+zp*zp),-2.)*(-3.*rm*(xp*xp)+3.*ri*xp*sqrt(xp*xp+yp*yp+zp*z\
p)))/sqrt(xp*xp+yp*yp+zp*zp);
dJ[2][0][1]=(2.*ri*(ri-rm)*yp*(-(rm*(ri*(yp*yp+zp*zp)+pow(xp,3.)))+ri*pow(xp\
*xp+yp*yp+zp*zp,1.5))*pow(rm*xp-ri*sqrt(xp*xp+yp*yp+zp*zp),-3.))/(xp*xp+yp*\
yp+zp*zp)-(ri-rm)*yp*pow(xp*xp+yp*yp+zp*zp,-1.5)*(-(rm*(ri*(yp*yp+zp*zp)+po\
w(xp,3.)))+ri*pow(xp*xp+yp*yp+zp*zp,1.5))*pow(rm*xp-ri*sqrt(xp*xp+yp*yp+zp*\
zp),-2.)+(1.*(ri-rm)*pow(rm*xp-ri*sqrt(xp*xp+yp*yp+zp*zp),-2.)*(-2.*ri*rm*y\
p+3.*ri*yp*sqrt(xp*xp+yp*yp+zp*zp)))/sqrt(xp*xp+yp*yp+zp*zp);
dJ[2][0][2]=(2.*ri*(ri-rm)*zp*(-(rm*(ri*(yp*yp+zp*zp)+pow(xp,3.)))+ri*pow(xp\
*xp+yp*yp+zp*zp,1.5))*pow(rm*xp-ri*sqrt(xp*xp+yp*yp+zp*zp),-3.))/(xp*xp+yp*\
yp+zp*zp)-(ri-rm)*zp*pow(xp*xp+yp*yp+zp*zp,-1.5)*(-(rm*(ri*(yp*yp+zp*zp)+po\
w(xp,3.)))+ri*pow(xp*xp+yp*yp+zp*zp,1.5))*pow(rm*xp-ri*sqrt(xp*xp+yp*yp+zp*\
zp),-2.)+(1.*(ri-rm)*pow(rm*xp-ri*sqrt(xp*xp+yp*yp+zp*zp),-2.)*(-2.*ri*rm*z\
p+3.*ri*zp*sqrt(xp*xp+yp*yp+zp*zp)))/sqrt(xp*xp+yp*yp+zp*zp);
dJ[2][1][0]=-((ri-rm)*rm*(ri-xp)*yp*(xp*xp)*pow(xp*xp+yp*yp+zp*zp,-1.5)*pow(\
rm*xp-ri*sqrt(xp*xp+yp*yp+zp*zp),-2.))+(1.*(ri-rm)*rm*(ri-xp)*yp*pow(rm*xp-\
ri*sqrt(xp*xp+yp*yp+zp*zp),-2.))/sqrt(xp*xp+yp*yp+zp*zp)-(1.*(ri-rm)*rm*xp*\
yp*pow(rm*xp-ri*sqrt(xp*xp+yp*yp+zp*zp),-2.))/sqrt(xp*xp+yp*yp+zp*zp)-(2.*(\
ri-rm)*rm*(ri-xp)*xp*yp*pow(rm*xp-ri*sqrt(xp*xp+yp*yp+zp*zp),-3.)*(rm-(1.*r\
i*xp)/sqrt(xp*xp+yp*yp+zp*zp)))/sqrt(xp*xp+yp*yp+zp*zp);
dJ[2][1][1]=(2.*ri*(ri-rm)*rm*(ri-xp)*xp*(yp*yp)*pow(rm*xp-ri*sqrt(xp*xp+yp*\
yp+zp*zp),-3.))/(xp*xp+yp*yp+zp*zp)-(ri-rm)*rm*(ri-xp)*xp*(yp*yp)*pow(xp*xp\
+yp*yp+zp*zp,-1.5)*pow(rm*xp-ri*sqrt(xp*xp+yp*yp+zp*zp),-2.)+(1.*(ri-rm)*rm\
*(ri-xp)*xp*pow(rm*xp-ri*sqrt(xp*xp+yp*yp+zp*zp),-2.))/sqrt(xp*xp+yp*yp+zp*\
zp);
dJ[2][1][2]=(2.*ri*(ri-rm)*rm*(ri-xp)*xp*yp*zp*pow(rm*xp-ri*sqrt(xp*xp+yp*yp\
+zp*zp),-3.))/(xp*xp+yp*yp+zp*zp)-(ri-rm)*rm*(ri-xp)*xp*yp*zp*pow(xp*xp+yp*\
yp+zp*zp,-1.5)*pow(rm*xp-ri*sqrt(xp*xp+yp*yp+zp*zp),-2.);
dJ[2][2][0]=-((ri-rm)*rm*(ri-xp)*zp*(xp*xp)*pow(xp*xp+yp*yp+zp*zp,-1.5)*pow(\
rm*xp-ri*sqrt(xp*xp+yp*yp+zp*zp),-2.))+(1.*(ri-rm)*rm*(ri-xp)*zp*pow(rm*xp-\
ri*sqrt(xp*xp+yp*yp+zp*zp),-2.))/sqrt(xp*xp+yp*yp+zp*zp)-(1.*(ri-rm)*rm*xp*\
zp*pow(rm*xp-ri*sqrt(xp*xp+yp*yp+zp*zp),-2.))/sqrt(xp*xp+yp*yp+zp*zp)-(2.*(\
ri-rm)*rm*(ri-xp)*xp*zp*pow(rm*xp-ri*sqrt(xp*xp+yp*yp+zp*zp),-3.)*(rm-(1.*r\
i*xp)/sqrt(xp*xp+yp*yp+zp*zp)))/sqrt(xp*xp+yp*yp+zp*zp);
dJ[2][2][1]=(2.*ri*(ri-rm)*rm*(ri-xp)*xp*yp*zp*pow(rm*xp-ri*sqrt(xp*xp+yp*yp\
+zp*zp),-3.))/(xp*xp+yp*yp+zp*zp)-(ri-rm)*rm*(ri-xp)*xp*yp*zp*pow(xp*xp+yp*\
yp+zp*zp,-1.5)*pow(rm*xp-ri*sqrt(xp*xp+yp*yp+zp*zp),-2.);
dJ[2][2][2]=(2.*ri*(ri-rm)*rm*(ri-xp)*xp*(zp*zp)*pow(rm*xp-ri*sqrt(xp*xp+yp*\
yp+zp*zp),-3.))/(xp*xp+yp*yp+zp*zp)-(ri-rm)*rm*(ri-xp)*xp*(zp*zp)*pow(xp*xp\
+yp*yp+zp*zp,-1.5)*pow(rm*xp-ri*sqrt(xp*xp+yp*yp+zp*zp),-2.)+(1.*(ri-rm)*rm\
*(ri-xp)*xp*pow(rm*xp-ri*sqrt(xp*xp+yp*yp+zp*zp),-2.))/sqrt(xp*xp+yp*yp+zp*\
zp);
