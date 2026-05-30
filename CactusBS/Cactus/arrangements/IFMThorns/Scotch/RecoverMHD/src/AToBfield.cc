#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"
#include "RecoverMHD.h"

// Second order f.d. for volume-centered quantity from volume-centered quantity
#define DIFF_2(q, iplus, iminus, inv_delta) \
              ( 0.5 * ( q[iplus] - q[iminus] ) * inv_delta )

// Fourth order f.d. for volume-centered quantity from volume-centered quantity
#define DIFF_4(q, ipp, ip, im, imm, inv_delta ) \
              (((-q[ipp] + 8.*q[ip] - 8.*q[im] + q[imm]) / 12.) *inv_delta )
 
/* Populate the B-field everywhere given a populated A-field
   Don't forget to sync afterwards. */

extern "C" void RecoverMHD_Avec2Bvec(CCTK_ARGUMENTS)
{

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  CCTK_REAL inv_dx, inv_dy, inv_dz;
  //CCTK_INT  nx, ny, nz;

  //nx = cctk_lsh[0];
  //ny = cctk_lsh[1];
  //nz = cctk_lsh[2];
  inv_dx = 1./CCTK_DELTA_SPACE(0);
  inv_dy = 1./CCTK_DELTA_SPACE(1);
  inv_dz = 1./CCTK_DELTA_SPACE(2);

  CCTK_LOOP3_INT( calcB_from_A_RecoverMHD, cctkGH, i,j,k )
//  for ( int k = whisky_stencil-1; k < nz-whisky_stencil; k++ )
//    for ( int j = whisky_stencil-1; j < ny-whisky_stencil; j++ )
//      for ( int i = whisky_stencil-1; i < nx-whisky_stencil; i++ )
      {

         
          CCTK_REAL Az_y, Ay_z, Ax_z, Az_x, Ay_x, Ax_y;
          CCTK_REAL det, inv_sqdet;

          CCTK_INT local_spatial_order = spatial_order;
          //if ( spatial_order > 2 ) { /* Check we're not too close to some boundary */
          //   for ( int kk = k-2; kk <= k+2; kk++ )
          //       for ( int jj = j-2; jj <= j+2; jj++ )
          //           for ( int ii = i-2; ii <= i+2; ii++ ) {
          //               CCTK_INT stencil_idx = CCTK_GFINDEX3D( cctkGH, ii, jj, kk );
          //               if ( emask[stencil_idx] < 0.75 )
          //                  local_spatial_order = 2;
          //   }
          //}

          /* Identify indicies */
          CCTK_INT idx   = CCTK_GFINDEX3D( cctkGH, i, j, k );
          CCTK_INT x_idx = CCTK_VECTGFINDEX3D( cctkGH, i, j, k, 0 );
          CCTK_INT y_idx = CCTK_VECTGFINDEX3D( cctkGH, i, j, k, 1 );
          CCTK_INT z_idx = CCTK_VECTGFINDEX3D( cctkGH, i, j, k, 2 );

          // x +- 1
          CCTK_INT y_xplus1 = CCTK_VECTGFINDEX3D( cctkGH, i+1, j, k, 1 );
          CCTK_INT y_xmin1  = CCTK_VECTGFINDEX3D( cctkGH, i-1, j, k, 1 );
          CCTK_INT z_xplus1 = CCTK_VECTGFINDEX3D( cctkGH, i+1, j, k, 2 );
          CCTK_INT z_xmin1  = CCTK_VECTGFINDEX3D( cctkGH, i-1, j, k, 2 );

          // y +- 1
          CCTK_INT x_yplus1 = CCTK_VECTGFINDEX3D( cctkGH, i, j+1, k, 0 );
          CCTK_INT x_ymin1  = CCTK_VECTGFINDEX3D( cctkGH, i, j-1, k, 0 );
          CCTK_INT z_yplus1 = CCTK_VECTGFINDEX3D( cctkGH, i, j+1, k, 2 );
          CCTK_INT z_ymin1  = CCTK_VECTGFINDEX3D( cctkGH, i, j-1, k, 2 );

          // z +- 1
          CCTK_INT x_zplus1 = CCTK_VECTGFINDEX3D( cctkGH, i, j, k+1, 0 );
          CCTK_INT x_zmin1  = CCTK_VECTGFINDEX3D( cctkGH, i, j, k-1, 0 );
          CCTK_INT y_zplus1 = CCTK_VECTGFINDEX3D( cctkGH, i, j, k+1, 1 );
          CCTK_INT y_zmin1  = CCTK_VECTGFINDEX3D( cctkGH, i, j, k-1, 1 );

          if (local_spatial_order == 2) {

             Az_y = DIFF_2(Avec, z_yplus1, z_ymin1, inv_dy);
             Ay_z = DIFF_2(Avec, y_zplus1, y_zmin1, inv_dz);
             Ax_z = DIFF_2(Avec, x_zplus1, x_zmin1, inv_dz);
             Az_x = DIFF_2(Avec, z_xplus1, z_xmin1, inv_dx);
             Ay_x = DIFF_2(Avec, y_xplus1, y_xmin1, inv_dx);
             Ax_y = DIFF_2(Avec, x_yplus1, x_ymin1, inv_dy);

          } else {

             // x +- 2
             CCTK_INT y_xplus2 = CCTK_VECTGFINDEX3D( cctkGH, i+2, j, k, 1 );
             CCTK_INT y_xmin2  = CCTK_VECTGFINDEX3D( cctkGH, i-2, j, k, 1 );
             CCTK_INT z_xplus2 = CCTK_VECTGFINDEX3D( cctkGH, i+2, j, k, 2 );
             CCTK_INT z_xmin2  = CCTK_VECTGFINDEX3D( cctkGH, i-2, j, k, 2 );
             // y +- 2
             CCTK_INT x_yplus2 = CCTK_VECTGFINDEX3D( cctkGH, i, j+2, k, 0 );
             CCTK_INT x_ymin2  = CCTK_VECTGFINDEX3D( cctkGH, i, j-2, k, 0 );
             CCTK_INT z_yplus2 = CCTK_VECTGFINDEX3D( cctkGH, i, j+2, k, 2 );
             CCTK_INT z_ymin2  = CCTK_VECTGFINDEX3D( cctkGH, i, j-2, k, 2 );
             // z +- 2
             CCTK_INT x_zplus2 = CCTK_VECTGFINDEX3D( cctkGH, i, j, k+2, 0 );
             CCTK_INT x_zmin2  = CCTK_VECTGFINDEX3D( cctkGH, i, j, k-2, 0 );
             CCTK_INT y_zplus2 = CCTK_VECTGFINDEX3D( cctkGH, i, j, k+2, 1 );
             CCTK_INT y_zmin2  = CCTK_VECTGFINDEX3D( cctkGH, i, j, k-2, 1 );

             Az_y = DIFF_4(Avec, z_yplus2, z_yplus1, z_ymin1, z_ymin2, inv_dy);
             Ay_z = DIFF_4(Avec, y_zplus2, y_zplus1, y_zmin1, y_zmin2, inv_dz);
             Ax_z = DIFF_4(Avec, x_zplus2, x_zplus1, x_zmin1, x_zmin2, inv_dz);
             Az_x = DIFF_4(Avec, z_xplus2, z_xplus1, z_xmin1, z_xmin2, inv_dx);
             Ay_x = DIFF_4(Avec, y_xplus2, y_xplus1, y_xmin1, y_xmin2, inv_dx);
             Ax_y = DIFF_4(Avec, x_yplus2, x_yplus1, x_ymin1, x_ymin2, inv_dy);
          }

          det = spatialdeterminant(gxx[idx],gxy[idx],gxz[idx],gyy[idx],gyz[idx],gzz[idx]);
          inv_sqdet = 1./sqrt(det);

          /* Set both B^i and Bn[xyz]! */
          CCTK_REAL BnxL = Az_y - Ay_z;
          CCTK_REAL BnyL = Ax_z - Az_x;
          CCTK_REAL BnzL = Ay_x - Ax_y;
          Bvec[x_idx] = inv_sqdet*BnxL;
          Bvec[y_idx] = inv_sqdet*BnyL;
          Bvec[z_idx] = inv_sqdet*BnzL;
          Bnx[idx] = BnxL; 
          Bny[idx] = BnyL; 
          Bnz[idx] = BnzL; 

      }
      CCTK_ENDLOOP3_INT(calcB_from_A_RecoverMHD);
}

