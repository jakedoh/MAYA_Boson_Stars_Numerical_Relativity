
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
dJ[0][1][1]=-2.*yp*zp*pow(yp*yp+zp*zp,-2.);
dJ[0][1][2]=(yp*yp-zp*zp)*pow(yp*yp+zp*zp,-2.);
dJ[0][2][0]=0.;
dJ[0][2][1]=(yp*yp-zp*zp)*pow(yp*yp+zp*zp,-2.);
dJ[0][2][2]=2.*yp*zp*pow(yp*yp+zp*zp,-2.);
dJ[1][0][0]=-2.*xp*zp*pow(xp*xp+zp*zp,-2.);
dJ[1][0][1]=0.;
dJ[1][0][2]=(xp*xp-zp*zp)*pow(xp*xp+zp*zp,-2.);
dJ[1][1][0]=0.;
dJ[1][1][1]=0.;
dJ[1][1][2]=0.;
dJ[1][2][0]=(xp*xp-zp*zp)*pow(xp*xp+zp*zp,-2.);
dJ[1][2][1]=0.;
dJ[1][2][2]=2.*xp*zp*pow(xp*xp+zp*zp,-2.);
dJ[2][0][0]=(yp*yp+zp*zp)*pow(xp*xp+yp*yp+zp*zp,-1.5);
dJ[2][0][1]=-(xp*yp*pow(xp*xp+yp*yp+zp*zp,-1.5));
dJ[2][0][2]=-(xp*zp*pow(xp*xp+yp*yp+zp*zp,-1.5));
dJ[2][1][0]=-(xp*yp*pow(xp*xp+yp*yp+zp*zp,-1.5));
dJ[2][1][1]=(xp*xp+zp*zp)*pow(xp*xp+yp*yp+zp*zp,-1.5);
dJ[2][1][2]=-(yp*zp*pow(xp*xp+yp*yp+zp*zp,-1.5));
dJ[2][2][0]=-(xp*zp*pow(xp*xp+yp*yp+zp*zp,-1.5));
dJ[2][2][1]=-(yp*zp*pow(xp*xp+yp*yp+zp*zp,-1.5));
dJ[2][2][2]=(xp*xp+yp*yp)*pow(xp*xp+yp*yp+zp*zp,-1.5);
