
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

/*
** Jacobian.h
** 
*/

#ifndef   	JACOBIAN_H_
# define   	JACOBIAN_H_


/*
 * apply_jacobian()
 *
 * Applies a Jacobian to a single index tensor or first derivative.
 *
 * The Jacobian components should be specified in the form:
 *
 * J[i][a] = dx^i/dy^a.
 *
 */
static inline __attribute__((always_inline)) void
apply_jacobian(const cGH* restrict const cctkGH,
               const CCTK_INT nvars,
               CCTK_REAL (*dvars)[3],
               CCTK_REAL (*xdvars)[3],
               CCTK_REAL J[3][3])
{
  for (int v=0; v<nvars; ++v)
  {
     xdvars[v][0] = J[0][0] * dvars[v][0]
                  + J[1][0] * dvars[v][1]
                  + J[2][0] * dvars[v][2];
     xdvars[v][1] = J[0][1] * dvars[v][0]
                  + J[1][1] * dvars[v][1]
                  + J[2][1] * dvars[v][2];
     xdvars[v][2] = J[0][2] * dvars[v][0]
                  + J[1][2] * dvars[v][1]
                  + J[2][2] * dvars[v][2];
  }
/*      for (int a=0; a<3; ++a)
        {
          xdvars[v][a] = 0.0;
          for (int i=0; i<3; ++i)
            xdvars[v][a] += J[i][a] * dvars[v][i];
        }
  } */
  return;
}

/*
 * apply_jacobian2()
 *
 * Applies a Jacobian to second derivatives.
 *
 * As input, this routine requires the standard Jacobians,
 *
 *   J[i][a] = dx^i/dy^a,
 *
 * as well as their derivatives,
 *
 *   K[i][a][b] = d^2 x^i / dy^a dy^b.
 *
 * Note that we take into account the index symmetry in [a,b],
 * and instead use a multi-index which ranges from 1...6:
 *   K[i][a][b] -> K[i][c]
 * where
 *   c = [1,2,3,4,5,6] = [[1,1], [1,2], [1,3], [2,2], [2,3], [3,3]]
 *
 */
static inline __attribute__((always_inline)) void
apply_jacobian2(const cGH* restrict const cctkGH,
                const CCTK_INT nvars,
                CCTK_REAL (*dvars)[3],
                CCTK_REAL (*ddvars)[6],
                CCTK_REAL (*xddvars)[6],
                CCTK_REAL J[3][3],
                CCTK_REAL K[3][6])
{
  for (int v=0; v<nvars; ++v)
  {
    xddvars[v][0] =       J[0][0] * J[0][0] * ddvars[v][0]
                  + 2.0 * J[0][0] * J[1][0] * ddvars[v][1]
                  + 2.0 * J[0][0] * J[2][0] * ddvars[v][2]
                  +       J[1][0] * J[1][0] * ddvars[v][3]
                  + 2.0 * J[1][0] * J[2][0] * ddvars[v][4]
                  +       J[2][0] * J[2][0] * ddvars[v][5]
                  + K[0][0] * dvars[v][0]
                  + K[1][0] * dvars[v][1]
                  + K[2][0] * dvars[v][2];
    xddvars[v][1] =   J[0][0] * J[0][1] * ddvars[v][0]
                  + ( J[0][0] * J[1][1] + J[1][0] * J[0][1] ) * ddvars[v][1]
                  + ( J[0][0] * J[2][1] + J[2][0] * J[0][1] ) * ddvars[v][2]
                  +   J[1][0] * J[1][1] * ddvars[v][3]
                  + ( J[1][0] * J[2][1] + J[2][0] * J[1][1] ) * ddvars[v][4]
                  +   J[2][0] * J[2][1] * ddvars[v][5]
                  + K[0][1] * dvars[v][0]
                  + K[1][1] * dvars[v][1]
                  + K[2][1] * dvars[v][2];
    xddvars[v][2] =   J[0][0] * J[0][2] * ddvars[v][0]
                  + ( J[0][0] * J[1][2] + J[1][0] * J[0][2] ) * ddvars[v][1]
                  + ( J[0][0] * J[2][2] + J[2][0] * J[0][2] ) * ddvars[v][2]
                  +   J[1][0] * J[1][2] * ddvars[v][3]
                  + ( J[1][0] * J[2][2] + J[2][0] * J[1][2] ) * ddvars[v][4]
                  +   J[2][0] * J[2][2] * ddvars[v][5]
                  + K[0][2] * dvars[v][0]
                  + K[1][2] * dvars[v][1]
                  + K[2][2] * dvars[v][2];
    xddvars[v][3] =       J[0][1] * J[0][1] * ddvars[v][0]
                  + 2.0 * J[0][1] * J[1][1] * ddvars[v][1]
                  + 2.0 * J[0][1] * J[2][1] * ddvars[v][2]
                  +       J[1][1] * J[1][1] * ddvars[v][3]
                  + 2.0 * J[1][1] * J[2][1] * ddvars[v][4]
                  +       J[2][1] * J[2][1] * ddvars[v][5]
                  + K[0][3] * dvars[v][0]
                  + K[1][3] * dvars[v][1]
                  + K[2][3] * dvars[v][2];
    xddvars[v][4] =   J[0][1] * J[0][2] * ddvars[v][0]
                  + ( J[0][1] * J[1][2] + J[1][1] * J[0][2] ) * ddvars[v][1]
                  + ( J[0][1] * J[2][2] + J[2][1] * J[0][2] ) * ddvars[v][2]
                  +   J[1][1] * J[1][2] * ddvars[v][3]
                  + ( J[1][1] * J[2][2] + J[2][1] * J[1][2] ) * ddvars[v][4]
                  +   J[2][1] * J[2][2] * ddvars[v][5]
                  + K[0][4] * dvars[v][0]
                  + K[1][4] * dvars[v][1]
                  + K[2][4] * dvars[v][2];
    xddvars[v][5] =       J[0][2] * J[0][2] * ddvars[v][0]
                  + 2.0 * J[0][2] * J[1][2] * ddvars[v][1]
                  + 2.0 * J[0][2] * J[2][2] * ddvars[v][2]
                  +       J[1][2] * J[1][2] * ddvars[v][3]
                  + 2.0 * J[1][2] * J[2][2] * ddvars[v][4]
                  +       J[2][2] * J[2][2] * ddvars[v][5]
                  + K[0][5] * dvars[v][0]
                  + K[1][5] * dvars[v][1]
                  + K[2][5] * dvars[v][2];
  }
//    {
//      int gcmpt = 0;
//      for (int a=0; a<3; ++a)
//        {
//          for (int b=a; b<3; ++b)
//            {
//              /*
//               * First derivative terms.
//               */
//              xddvars[v][gcmpt] = 0;
//
//              /*
//              for (int i=0; i<3; ++i)
//                xddvars[v][gcmpt] += K[i][gcmpt] * dvars[v][i];
//              */
//              
//              /*
//               * Second derivative terms, upper diagonal.
//               */
//              int lcmpt = 0;
//              for (int i=0; i<3; ++i)
//                {
//                  for (int j=i; j<3; ++j)
//                    {
//                      xddvars[v][gcmpt] += J[i][a]*J[j][b] * ddvars[v][lcmpt];
//                      ++lcmpt;
//                    }
//                }
// 
//              /*
//               * Second derivative terms, lower diagonal.
//               */
//              lcmpt = 1;
//              for (int i=0; i<2; ++i)
//                {
//                  for (int j=i+1; j<3; ++j)
//                    {
//                      xddvars[v][gcmpt] += 2.0 * J[i][a]*J[j][b] 
//			* ddvars[v][lcmpt];
//                      lcmpt*=2;
//                    }
//                }
//              ++gcmpt;
//            }
//        }
//    }

  return;
}

#endif 	    /* !JACOBIAN_H_ */
