
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

if (*general_coordinates)
  {
    dadx = J11[ijk];
    dady = J12[ijk];
    dadz = J13[ijk];
    dbdx = J21[ijk];
    dbdy = J22[ijk];
    dbdz = J23[ijk];
    dcdx = J31[ijk];
    dcdy = J32[ijk];
    dcdz = J33[ijk];
    
  }


g_diff(cctkGH, gxx, &dg[0][0][0], &dg[0][0][1], &dg[0][0][2],
       *general_coordinates,
       dadx, dbdx, dcdx, dady, dbdy, dcdy, dadz, dbdz, dcdz,
       i, j, k, ni, nj, nk,
       imin[0], imax[0], imin[1], imax[1], imin[2], imax[2], q[0], q[1], q[2],
       ihx, ihy, ihz);

g_diff(cctkGH, gxy, &dg[0][1][0], &dg[0][1][1], &dg[0][1][2],
       *general_coordinates,
       dadx, dbdx, dcdx, dady, dbdy, dcdy, dadz, dbdz, dcdz,
       i, j, k, ni, nj, nk,
       imin[0], imax[0], imin[1], imax[1], imin[2], imax[2], q[0], q[1], q[2],
       ihx, ihy, ihz);

g_diff(cctkGH, gxz, &dg[0][2][0], &dg[0][2][1], &dg[0][2][2],
       *general_coordinates,
       dadx, dbdx, dcdx, dady, dbdy, dcdy, dadz, dbdz, dcdz,
       i, j, k, ni, nj, nk,
       imin[0], imax[0], imin[1], imax[1], imin[2], imax[2], q[0], q[1], q[2],
       ihx, ihy, ihz);

g_diff(cctkGH, gyy, &dg[1][1][0], &dg[1][1][1], &dg[1][1][2],
       *general_coordinates,
       dadx, dbdx, dcdx, dady, dbdy, dcdy, dadz, dbdz, dcdz,
       i, j, k, ni, nj, nk,
       imin[0], imax[0], imin[1], imax[1], imin[2], imax[2], q[0], q[1], q[2],
       ihx, ihy, ihz);

g_diff(cctkGH, gyz, &dg[1][2][0], &dg[1][2][1], &dg[1][2][2],
       *general_coordinates,
       dadx, dbdx, dcdx, dady, dbdy, dcdy, dadz, dbdz, dcdz,
       i, j, k, ni, nj, nk,
       imin[0], imax[0], imin[1], imax[1], imin[2], imax[2], q[0], q[1], q[2],
       ihx, ihy, ihz);

g_diff(cctkGH, gzz, &dg[2][2][0], &dg[2][2][1], &dg[2][2][2],
       *general_coordinates,
       dadx, dbdx, dcdx, dady, dbdy, dcdy, dadz, dbdz, dcdz,
       i, j, k, ni, nj, nk,
       imin[0], imax[0], imin[1], imax[1], imin[2], imax[2], q[0], q[1], q[2],
       ihx, ihy, ihz);

g_diff(cctkGH, betax, &dxBetax, &dyBetax, &dzBetax,
       *general_coordinates,
       dadx, dbdx, dcdx, dady, dbdy, dcdy, dadz, dbdz, dcdz,
       i, j, k, ni, nj, nk,
       imin[0], imax[0], imin[1], imax[1], imin[2], imax[2], q[0], q[1], q[2],
       ihx, ihy, ihz);

g_diff(cctkGH, betay, &dxBetay, &dyBetay, &dzBetay,
       *general_coordinates,
       dadx, dbdx, dcdx, dady, dbdy, dcdy, dadz, dbdz, dcdz,
       i, j, k, ni, nj, nk,
       imin[0], imax[0], imin[1], imax[1], imin[2], imax[2], q[0], q[1], q[2],
       ihx, ihy, ihz);

g_diff(cctkGH, betaz, &dxBetaz, &dyBetaz, &dzBetaz,
       *general_coordinates,
       dadx, dbdx, dcdx, dady, dbdy, dcdy, dadz, dbdz, dcdz,
       i, j, k, ni, nj, nk,
       imin[0], imax[0], imin[1], imax[1], imin[2], imax[2], q[0], q[1], q[2],
       ihx, ihy, ihz);

g_diff(cctkGH, alp, &dxalp, &dyalp, &dzalp,
       *general_coordinates,
       dadx, dbdx, dcdx, dady, dbdy, dcdy, dadz, dbdz, dcdz,
       i, j, k, ni, nj, nk,
       imin[0], imax[0], imin[1], imax[1], imin[2], imax[2], q[0], q[1], q[2],
       ihx, ihy, ihz);

