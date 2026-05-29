C slice_normal.inc -- common code for computing n^A and dx^A/dt,
C                     used in ../slice_data.F and ../slice_evolve.F
C $Header$

C     Fill in remaining components of metric.
      do m=1,4
         do n=1,m-1
            gu(m,n) = gu(n,m)
            gd(m,n) = gd(n,m)
         end do
      end do

C     Calculate normal vector with respect to the slice, with indices down,
C     and not yet normalized.
C     n_A = epsilon_ABCD X^A_,1 X^B_,2 X^C_,3.
C     Easiest to hardwired this with an editor.
               nd(1) = s1d(2,1)*s1d(3,2)*s1d(4,3)
     $              + s1d(3,1)*s1d(4,2)*s1d(2,3)
     $              + s1d(4,1)*s1d(2,2)*s1d(3,3)
     $              - s1d(2,1)*s1d(4,2)*s1d(3,3)
     $              - s1d(4,1)*s1d(3,2)*s1d(2,3)
     $              - s1d(3,1)*s1d(2,2)*s1d(4,3)

               nd(2) = - s1d(3,1)*s1d(4,2)*s1d(1,3)
     $              - s1d(4,1)*s1d(1,2)*s1d(3,3)
     $              -s1d(1,1)*s1d(3,2)*s1d(4,3)
     $              + s1d(3,1)*s1d(1,2)*s1d(4,3)
     $              + s1d(1,1)*s1d(4,2)*s1d(3,3)
     $              + s1d(4,1)*s1d(3,2)*s1d(1,3)

               nd(3) = s1d(4,1)*s1d(1,2)*s1d(2,3)
     $              + s1d(1,1)*s1d(2,2)*s1d(4,3)
     $              + s1d(2,1)*s1d(4,2)*s1d(1,3)
     $              - s1d(4,1)*s1d(2,2)*s1d(1,3)
     $              - s1d(2,1)*s1d(1,2)*s1d(4,3)
     $              - s1d(1,1)*s1d(4,2)*s1d(2,3)

               nd(4) = s1d(1,1)*s1d(2,2)*s1d(3,3)
     $              + s1d(2,1)*s1d(3,2)*s1d(1,3)
     $              + s1d(3,1)*s1d(1,2)*s1d(2,3)
     $              - s1d(1,1)*s1d(3,2)*s1d(2,3)
     $              - s1d(3,1)*s1d(2,2)*s1d(1,3)
     $              - s1d(2,1)*s1d(1,2)*s1d(3,3)

C     Normalize normal vector and raise its index, both with the exact
C     solution metric.
C     - N^2 = g^AB n_A n_B).
               norm = 0.d0
               do m=1,4
                  do n=1,4
                     norm = norm + nd(m) * nd(n) * gu(m,n)
                  end do
               end do
C     Debugging test.
               if (norm .ge. 0.d0) then
                  write(6,*) 'slice normal no longer timelike'
                  write(6,*) 'grid point', i, j, k
                  write(6,*) 'x^i', x(i,j,k), y(i,j,k), z(i,j,k)
                  write(6,*) 'x^A', slicex(i,j,k), slicey(i,j,k), 
     $                 slicez(i,j,k), slicet(i,j,k)
                  write(6,*) 'n_A', nd(1), nd(2), nd(3), nd(4)
                  write(6,*) 'n_A n^A', norm
                  write(6,*)  'dx^A/dx^i'
                  do m=1,4
                     do n=1,3
                        write(6,*) m, n, s1d(m,n)
                     end do
                  end do
                  write(6,*)  'g^AB'
                  do m=1,4
                     do n=1,4
                        write(6,*) m, n, gu(m,n)
                     end do
                  end do
                  call CCTK_WARN (0, "aborting")
               end if
C     And now n^A = - 1/N g^AB n_B. (Sign so that it is future-pointing.)
               do m=1,4
                  nu(m) = 0.d0
                  do n=1,4
                     nu(m) = nu(m) + gu(m,n) * nd(n)
                  end do
                  nu(m) = - nu(m) / sqrt(-norm)
               end do

C     dX^A/dt = alpha n^A + beta^i X^A_,i. Store as a GF.
               slicetmp2x(i,j,k) = alp(i,j,k) * nu(1)
     $              + betax(i,j,k) * s1d(1,1)
     $              + betay(i,j,k) * s1d(1,2)
     $              + betaz(i,j,k) * s1d(1,3) 

               slicetmp2y(i,j,k) = alp(i,j,k) * nu(2)
     $              + betax(i,j,k) * s1d(2,1)
     $              + betay(i,j,k) * s1d(2,2)
     $              + betaz(i,j,k) * s1d(2,3) 

               slicetmp2z(i,j,k) = alp(i,j,k) * nu(3)
     $              + betax(i,j,k) * s1d(3,1)
     $              + betay(i,j,k) * s1d(3,2)
     $              + betaz(i,j,k) * s1d(3,3) 

               slicetmp2t(i,j,k) = alp(i,j,k) * nu(4)
     $              + betax(i,j,k) * s1d(4,1)
     $              + betay(i,j,k) * s1d(4,2)
     $              + betaz(i,j,k) * s1d(4,3) 


