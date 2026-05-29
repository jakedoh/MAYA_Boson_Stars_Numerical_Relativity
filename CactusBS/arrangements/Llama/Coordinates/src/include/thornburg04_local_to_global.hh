
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
R=(-0.25*(h0+h1)*(-Rmax-Rmin))/h0+(0.5*(h0+h1)*(Rl+0.5*(-Rmax-R\
min)-Rstart))/h0+Rstart-(0.03125*(1.-(0.5*(h0+h1))/h0)*(1.+16.\
*pow(-Rmax-Rmin,2.)*pow(Rmax-Rmin,-2.))*pow(Rmax-Rmin,2.))/(-R\
max-Rmin)+(0.03125*(1.-(0.5*(h0+h1))/h0)*pow(Rmax-Rmin,2.)*sqr\
t(1.+16.*pow(-Rmax-Rmin,2.)*pow(Rmax-Rmin,-2.))*sqrt(1.+64.*po\
w(Rmax-Rmin,-2.)*pow(Rl+0.5*(-Rmax-Rmin)-Rstart,2.)))/(-Rmax-R\
min);
