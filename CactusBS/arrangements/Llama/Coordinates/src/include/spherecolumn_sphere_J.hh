
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
J[0][0]=(1.*xp*zp)/((xp*xp+yp*yp+zp*zp)*sqrt(xp*xp+yp*yp));
J[0][1]=(1.*yp*zp)/((xp*xp+yp*yp+zp*zp)*sqrt(xp*xp+yp*yp));
J[0][2]=-(sqrt(xp*xp+yp*yp)/(xp*xp+yp*yp+zp*zp));
J[1][0]=-(yp/(xp*xp+yp*yp));
J[1][1]=xp/(xp*xp+yp*yp);
J[1][2]=0.;
J[2][0]=(-0.5*h0*((128.*h0*Rmax*Rmin*xp)/sqrt(xp*xp+yp*yp+zp*zp)+(128.*h\
1*Rmax*Rmin*xp)/sqrt(xp*xp+yp*yp+zp*zp)+(64.*h0*xp*(Rmax*Rmax))/sqrt(xp\
*xp+yp*yp+zp*zp)+(64.*h1*xp*(Rmax*Rmax))/sqrt(xp*xp+yp*yp+zp*zp)+(64.*h\
0*xp*(Rmin*Rmin))/sqrt(xp*xp+yp*yp+zp*zp)+(64.*h1*xp*(Rmin*Rmin))/sqrt(\
xp*xp+yp*yp+zp*zp)-(1.*pow(h0-h1,2.)*(47.*Rmin*(Rmax*Rmax)+47.*Rmax*(Rm\
in*Rmin)+17.*pow(Rmax,3.)+17.*pow(Rmin,3.))*(h0*(128.*(Rmax+Rmin)*xp-(2\
.*xp*(62.*Rmax*Rmin+33.*(Rmax*Rmax)+33.*(Rmin*Rmin)))/sqrt(xp*xp+yp*yp+\
zp*zp))+(2.*h1*xp*pow(Rmax-Rmin,2.))/sqrt(xp*xp+yp*yp+zp*zp)))/sqrt(h0*\
pow(h0-h1,2.)*(47.*Rmin*(Rmax*Rmax)+47.*Rmax*(Rmin*Rmin)+17.*pow(Rmax,3\
.)+17.*pow(Rmin,3.))*(2.*h1*pow(Rmax-Rmin,2.)*sqrt(xp*xp+yp*yp+zp*zp)+h\
0*(47.*Rmin*(Rmax*Rmax)+47.*Rmax*(Rmin*Rmin)+64.*(Rmax+Rmin)*(xp*xp+yp*\
yp+zp*zp)+17.*pow(Rmax,3.)+17.*pow(Rmin,3.)-2.*(62.*Rmax*Rmin+33.*(Rmax\
*Rmax)+33.*(Rmin*Rmin))*sqrt(xp*xp+yp*yp+zp*zp))))))/(-2.*h0*h1*(62.*Rm\
ax*Rmin+33.*(Rmax*Rmax)+33.*(Rmin*Rmin))+h0*h0*pow(Rmax-Rmin,2.)+h1*h1*\
pow(Rmax-Rmin,2.));
J[2][1]=(-0.5*h0*((128.*h0*Rmax*Rmin*yp)/sqrt(xp*xp+yp*yp+zp*zp)+(128.*h\
1*Rmax*Rmin*yp)/sqrt(xp*xp+yp*yp+zp*zp)+(64.*h0*yp*(Rmax*Rmax))/sqrt(xp\
*xp+yp*yp+zp*zp)+(64.*h1*yp*(Rmax*Rmax))/sqrt(xp*xp+yp*yp+zp*zp)+(64.*h\
0*yp*(Rmin*Rmin))/sqrt(xp*xp+yp*yp+zp*zp)+(64.*h1*yp*(Rmin*Rmin))/sqrt(\
xp*xp+yp*yp+zp*zp)-(1.*pow(h0-h1,2.)*(47.*Rmin*(Rmax*Rmax)+47.*Rmax*(Rm\
in*Rmin)+17.*pow(Rmax,3.)+17.*pow(Rmin,3.))*(h0*(128.*(Rmax+Rmin)*yp-(2\
.*yp*(62.*Rmax*Rmin+33.*(Rmax*Rmax)+33.*(Rmin*Rmin)))/sqrt(xp*xp+yp*yp+\
zp*zp))+(2.*h1*yp*pow(Rmax-Rmin,2.))/sqrt(xp*xp+yp*yp+zp*zp)))/sqrt(h0*\
pow(h0-h1,2.)*(47.*Rmin*(Rmax*Rmax)+47.*Rmax*(Rmin*Rmin)+17.*pow(Rmax,3\
.)+17.*pow(Rmin,3.))*(2.*h1*pow(Rmax-Rmin,2.)*sqrt(xp*xp+yp*yp+zp*zp)+h\
0*(47.*Rmin*(Rmax*Rmax)+47.*Rmax*(Rmin*Rmin)+64.*(Rmax+Rmin)*(xp*xp+yp*\
yp+zp*zp)+17.*pow(Rmax,3.)+17.*pow(Rmin,3.)-2.*(62.*Rmax*Rmin+33.*(Rmax\
*Rmax)+33.*(Rmin*Rmin))*sqrt(xp*xp+yp*yp+zp*zp))))))/(-2.*h0*h1*(62.*Rm\
ax*Rmin+33.*(Rmax*Rmax)+33.*(Rmin*Rmin))+h0*h0*pow(Rmax-Rmin,2.)+h1*h1*\
pow(Rmax-Rmin,2.));
J[2][2]=(-0.5*h0*((128.*h0*Rmax*Rmin*zp)/sqrt(xp*xp+yp*yp+zp*zp)+(128.*h\
1*Rmax*Rmin*zp)/sqrt(xp*xp+yp*yp+zp*zp)+(64.*h0*zp*(Rmax*Rmax))/sqrt(xp\
*xp+yp*yp+zp*zp)+(64.*h1*zp*(Rmax*Rmax))/sqrt(xp*xp+yp*yp+zp*zp)+(64.*h\
0*zp*(Rmin*Rmin))/sqrt(xp*xp+yp*yp+zp*zp)+(64.*h1*zp*(Rmin*Rmin))/sqrt(\
xp*xp+yp*yp+zp*zp)-(1.*pow(h0-h1,2.)*(47.*Rmin*(Rmax*Rmax)+47.*Rmax*(Rm\
in*Rmin)+17.*pow(Rmax,3.)+17.*pow(Rmin,3.))*(h0*(128.*(Rmax+Rmin)*zp-(2\
.*zp*(62.*Rmax*Rmin+33.*(Rmax*Rmax)+33.*(Rmin*Rmin)))/sqrt(xp*xp+yp*yp+\
zp*zp))+(2.*h1*zp*pow(Rmax-Rmin,2.))/sqrt(xp*xp+yp*yp+zp*zp)))/sqrt(h0*\
pow(h0-h1,2.)*(47.*Rmin*(Rmax*Rmax)+47.*Rmax*(Rmin*Rmin)+17.*pow(Rmax,3\
.)+17.*pow(Rmin,3.))*(2.*h1*pow(Rmax-Rmin,2.)*sqrt(xp*xp+yp*yp+zp*zp)+h\
0*(47.*Rmin*(Rmax*Rmax)+47.*Rmax*(Rmin*Rmin)+64.*(Rmax+Rmin)*(xp*xp+yp*\
yp+zp*zp)+17.*pow(Rmax,3.)+17.*pow(Rmin,3.)-2.*(62.*Rmax*Rmin+33.*(Rmax\
*Rmax)+33.*(Rmin*Rmin))*sqrt(xp*xp+yp*yp+zp*zp))))))/(-2.*h0*h1*(62.*Rm\
ax*Rmin+33.*(Rmax*Rmax)+33.*(Rmin*Rmin))+h0*h0*pow(Rmax-Rmin,2.)+h1*h1*\
pow(Rmax-Rmin,2.));
