! Finding isosurface(s).
! $Header$

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

  module areac
  contains
  subroutine area_computation(x_coord, y_coord, z_coord, t_index, &
                              g11, g12, g13, &
                              g22, g23, g33, &
                              npoints, l, area )

  use EHFinder_mod
  implicit none
  CCTK_INT,intent(in) :: npoints, l
  CCTK_REAL,dimension(:),pointer :: x_coord, y_coord, z_coord
  CCTK_INT, dimension(:),pointer :: t_index
  CCTK_REAL, dimension(npoints), intent(in) :: g11,g12,g13,g22,g23,g33
  CCTK_REAL::ds1,ds2,ds3,ds,dx1,dx2,dx3,dx_1,dx_2,dx_3,average_g11,average_g12,average_g13,average_g22,average_g33,average_g23,inverse_g11,inverse_g12,inverse_g13,inverse_g22,inverse_g33,inverse_g23,g
  CCTK_INT::i,point1,point2,point3
  CCTK_REAL,intent(out)::area

  print*, npoints, l
  print*, t_index(1:l)
  print*, 'x = ',x_coord
  print*, 'g11 = ',g11
  do i=1,l,3
!computing the determinant of the locale averaged metric (the metric is equal with the average of the metric in the corners of the triangle
!we also use point1, point2, point3 as references from t_index that gives the corners of the current triangle  
    point1=t_index(i)
    point2=t_index(i+1)
    point3=t_index(i+2)
    average_g11=(g11(point1)+g11(point2)+g11(point3))/3
    average_g12=(g12(point1)+g12(point2)+g12(point3))/3
    average_g13=(g13(point1)+g13(point2)+g13(point3))/3
    average_g22=(g22(point1)+g22(point2)+g22(point3))/3
    average_g23=(g23(point1)+g23(point2)+g23(point3))/3
    average_g33=(g33(point1)+g33(point2)+g33(point3))/3
    g=average_g11*(average_g22*average_g33-average_g23*average_g23)-average_g12*(average_g12*average_g33-average_g23*average_g13)+average_g13*(average_g12*average_g23-average_g22*average_g13)
!computing inverse of the g-matrix
    inverse_g11=g*(average_g22*average_g33-average_g23*average_g23)
    inverse_g22=g*(average_g11*average_g33-average_g13*average_g13)
    inverse_g33=g*(average_g11*average_g22-average_g12*average_g12)
    inverse_g12=g*(average_g13*average_g23-average_g12*average_g33)
    inverse_g13=g*(average_g12*average_g23-average_g13*average_g22)
    inverse_g23=g*(average_g13*average_g12-average_g11-average_g23)
!computing length elements and area of an elementar triangle
    dx1=x_coord(point2)-x_coord(point1)
    dx2=y_coord(point2)-y_coord(point1)
    dx3=z_coord(point2)-z_coord(point1)
    dx_1=x_coord(point3)-x_coord(point1)
    dx_2=y_coord(point3)-y_coord(point1)
    dx_3=z_coord(point3)-z_coord(point1)
    ds1=0.25*(g**0.5)*(dx2*dx_3-dx3*dx_2)
    ds2=0.25*(g**0.5)*(dx3*dx_1-dx1*dx_3)
    ds3=0.25*(g**0.5)*(dx1*dx_2-dx2*dx_1)
    ds=0.5*(inverse_g11*ds1*ds1+inverse_g22*ds2*ds2+inverse_g33*ds3*ds3+2*inverse_g12*ds1*ds2+2*inverse_g13*ds1*ds3+2*inverse_g23*ds2*ds3)
    area=area+ds
  end do
  end subroutine area_computation
  end module


 module triangular
 contains
subroutine triangularization(n_coords,n_triangles,Nx,Ny,Nz,x_coord,y_coord,z_coord,t_index,grid_value,grid_x,grid_y,grid_z)
!subroutine triangularization(l,contor,nx,ny,nz,x_coord,y_coord,z_coord,t_index,f(:,:,:,1),x(:,1,1),y(1,:,1),z(1,1,:))

  use IsoSurface_mod
  implicit none 
!declaring variables that will be used in the program
  CCTK_INT,intent(in)::Nx,Ny,Nz
  CCTK_REAL:: u,ABx,ABy,ABz,BCx,BCy,BCz
  CCTK_REAL,dimension(10)::temp_x,temp_y,temp_z
  CCTK_REAL,dimension(Nx,Ny,Nz),intent(in)::grid_value
  CCTK_REAL,dimension(Nx),intent(in)::grid_x
  CCTK_REAL,dimension(Ny),intent(in)::grid_y
  CCTK_REAL,dimension(Nz),intent(in)::grid_z
  CCTK_INT::jj,kk,sum,nsize,boolean_finder1,boolean_finder2,boolean_finder3,t,equality,index2,j,i,k,l,contor
  CCTK_INT,intent(out)::n_coords,n_triangles
  CCTK_REAL, dimension(:),pointer:: x_coord,y_coord,z_coord
  CCTK_INT, dimension (:),pointer::t_index
  CCTK_INT,dimension(:),allocatable::resize_t_index

!we store table with triangles in a 2 dimensional constant array using reshape function; Fortran stores column-like and we need row-like
  CCTK_INT, dimension (256,16):: a

  integer, dimension (16,30):: b1(16,30)=reshape&
((/-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
0, 1, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
1, 8, 3, 9, 8, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
0, 8, 3, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
9, 2, 10, 0, 2, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
2, 8, 3, 2, 10, 8, 10, 9, 8, -1, -1, -1, -1, -1, -1, -1,&
3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
0, 11, 2, 8, 11, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
1, 9, 0, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
1, 11, 2, 1, 9, 11, 9, 8, 11, -1, -1, -1, -1, -1, -1, -1,&
3, 10, 1, 11, 10, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
0, 10, 1, 0, 8, 10, 8, 11, 10, -1, -1, -1, -1, -1, -1, -1,&
3, 9, 0, 3, 11, 9, 11, 10, 9, -1, -1, -1, -1, -1, -1, -1,&
9, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
4, 3, 0, 7, 3, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
0, 1, 9, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
4, 1, 9, 4, 7, 1, 7, 3, 1, -1, -1, -1, -1, -1, -1, -1,&
1, 2, 10, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
3, 4, 7, 3, 0, 4, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1,&
9, 2, 10, 9, 0, 2, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1,&
2, 10, 9, 2, 9, 7, 2, 7, 3, 7, 9, 4, -1, -1, -1, -1,&
8, 4, 7, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
11, 4, 7, 11, 2, 4, 2, 0, 4, -1, -1, -1, -1, -1, -1, -1,&
9, 0, 1, 8, 4, 7, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1,&
4, 7, 11, 9, 4, 11, 9, 11, 2, 9, 2, 1, -1, -1, -1, -1,&
3, 10, 1, 3, 11, 10, 7, 8, 4, -1, -1, -1, -1, -1, -1, -1,&
1, 11, 10, 1, 4, 11, 1, 0, 4, 7, 11, 4, -1, -1, -1, -1/),(/16,30/))





integer, dimension (16,30):: b2(16,30)=reshape&
((/4, 7, 8, 9, 0, 11, 9, 11, 10, 11, 0, 3, -1, -1, -1, -1,&
4, 7, 11, 4, 11, 9, 9, 11, 10, -1, -1, -1, -1, -1, -1, -1,&
9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
9, 5, 4, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
0, 5, 4, 1, 5, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
8, 5, 4, 8, 3, 5, 3, 1, 5, -1, -1, -1, -1, -1, -1, -1,&
1, 2, 10, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
3, 0, 8, 1, 2, 10, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1,&
5, 2, 10, 5, 4, 2, 4, 0, 2, -1, -1, -1, -1, -1, -1, -1,&
2, 10, 5, 3, 2, 5, 3, 5, 4, 3, 4, 8, -1, -1, -1, -1,&
9, 5, 4, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
0, 11, 2, 0, 8, 11, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1,&
0, 5, 4, 0, 1, 5, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1,&
2, 1, 5, 2, 5, 8, 2, 8, 11, 4, 8, 5, -1, -1, -1, -1,&
10, 3, 11, 10, 1, 3, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1,&
4, 9, 5, 0, 8, 1, 8, 10, 1, 8, 11, 10, -1, -1, -1, -1,&
5, 4, 0, 5, 0, 11, 5, 11, 10, 11, 0, 3, -1, -1, -1, -1,&
5, 4, 8, 5, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1,&
9, 7, 8, 5, 7, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
9, 3, 0, 9, 5, 3, 5, 7, 3, -1, -1, -1, -1, -1, -1, -1,&
0, 7, 8, 0, 1, 7, 1, 5, 7, -1, -1, -1, -1, -1, -1, -1,&
1, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
9, 7, 8, 9, 5, 7, 10, 1, 2, -1, -1, -1, -1, -1, -1, -1,&
10, 1, 2, 9, 5, 0, 5, 3, 0, 5, 7, 3, -1, -1, -1, -1,&
8, 0, 2, 8, 2, 5, 8, 5, 7, 10, 5, 2, -1, -1, -1, -1,&
2, 10, 5, 2, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1,&
7, 9, 5, 7, 8, 9, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1,&
9, 5, 7, 9, 7, 2, 9, 2, 0, 2, 7, 11, -1, -1, -1, -1,&
2, 3, 11, 0, 1, 8, 1, 7, 8, 1, 5, 7, -1, -1, -1, -1,&
11, 2, 1, 11, 1, 7, 7, 1, 5, -1, -1, -1, -1, -1, -1, -1/),(/16,30/))



integer, dimension (16,30):: b3(16,30)=reshape&
((/9, 5, 8, 8, 5, 7, 10, 1, 3, 10, 3, 11, -1, -1, -1, -1,&
5, 7, 0, 5, 0, 9, 7, 11, 0, 1, 0, 10, 11, 10, 0, -1,&
11, 10, 0, 11, 0, 3, 10, 5, 0, 8, 0, 7, 5, 7, 0, -1,&
11, 10, 5, 7, 11, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
0, 8, 3, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
9, 0, 1, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
1, 8, 3, 1, 9, 8, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1,&
1, 6, 5, 2, 6, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
1, 6, 5, 1, 2, 6, 3, 0, 8, -1, -1, -1, -1, -1, -1, -1,&
9, 6, 5, 9, 0, 6, 0, 2, 6, -1, -1, -1, -1, -1, -1, -1,&
5, 9, 8, 5, 8, 2, 5, 2, 6, 3, 2, 8, -1, -1, -1, -1,&
2, 3, 11, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
11, 0, 8, 11, 2, 0, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1,&
0, 1, 9, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1,&
5, 10, 6, 1, 9, 2, 9, 11, 2, 9, 8, 11, -1, -1, -1, -1,&
6, 3, 11, 6, 5, 3, 5, 1, 3, -1, -1, -1, -1, -1, -1, -1,&
0, 8, 11, 0, 11, 5, 0, 5, 1, 5, 11, 6, -1, -1, -1, -1,&
3, 11, 6, 0, 3, 6, 0, 6, 5, 0, 5, 9, -1, -1, -1, -1,&
6, 5, 9, 6, 9, 11, 11, 9, 8, -1, -1, -1, -1, -1, -1, -1,&
5, 10, 6, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
4, 3, 0, 4, 7, 3, 6, 5, 10, -1, -1, -1, -1, -1, -1, -1,&
1, 9, 0, 5, 10, 6, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1,&
10, 6, 5, 1, 9, 7, 1, 7, 3, 7, 9, 4, -1, -1, -1, -1,&
6, 1, 2, 6, 5, 1, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1,&
1, 2, 5, 5, 2, 6, 3, 0, 4, 3, 4, 7, -1, -1, -1, -1,&
8, 4, 7, 9, 0, 5, 0, 6, 5, 0, 2, 6, -1, -1, -1, -1,&
7, 3, 9, 7, 9, 4, 3, 2, 9, 5, 9, 6, 2, 6, 9, -1,&
3, 11, 2, 7, 8, 4, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1,&
5, 10, 6, 4, 7, 2, 4, 2, 0, 2, 7, 11, -1, -1, -1, -1/),(/16,30/))



integer, dimension (16,30):: b4(16,30)=reshape&
((/0, 1, 9, 4, 7, 8, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1,&
9, 2, 1, 9, 11, 2, 9, 4, 11, 7, 11, 4, 5, 10, 6, -1,&
8, 4, 7, 3, 11, 5, 3, 5, 1, 5, 11, 6, -1, -1, -1, -1,&
5, 1, 11, 5, 11, 6, 1, 0, 11, 7, 11, 4, 0, 4, 11, -1,&
0, 5, 9, 0, 6, 5, 0, 3, 6, 11, 6, 3, 8, 4, 7, -1,&
6, 5, 9, 6, 9, 11, 4, 7, 9, 7, 11, 9, -1, -1, -1, -1,&
10, 4, 9, 6, 4, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
4, 10, 6, 4, 9, 10, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1,&
10, 0, 1, 10, 6, 0, 6, 4, 0, -1, -1, -1, -1, -1, -1, -1,&
8, 3, 1, 8, 1, 6, 8, 6, 4, 6, 1, 10, -1, -1, -1, -1,&
1, 4, 9, 1, 2, 4, 2, 6, 4, -1, -1, -1, -1, -1, -1, -1,&
3, 0, 8, 1, 2, 9, 2, 4, 9, 2, 6, 4, -1, -1, -1, -1,&
0, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
8, 3, 2, 8, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1,&
10, 4, 9, 10, 6, 4, 11, 2, 3, -1, -1, -1, -1, -1, -1, -1,&
0, 8, 2, 2, 8, 11, 4, 9, 10, 4, 10, 6, -1, -1, -1, -1,&
3, 11, 2, 0, 1, 6, 0, 6, 4, 6, 1, 10, -1, -1, -1, -1,&
6, 4, 1, 6, 1, 10, 4, 8, 1, 2, 1, 11, 8, 11, 1, -1,&
9, 6, 4, 9, 3, 6, 9, 1, 3, 11, 6, 3, -1, -1, -1, -1,&
8, 11, 1, 8, 1, 0, 11, 6, 1, 9, 1, 4, 6, 4, 1, -1,&
3, 11, 6, 3, 6, 0, 0, 6, 4, -1, -1, -1, -1, -1, -1, -1,&
6, 4, 8, 11, 6, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
7, 10, 6, 7, 8, 10, 8, 9, 10, -1, -1, -1, -1, -1, -1, -1,&
0, 7, 3, 0, 10, 7, 0, 9, 10, 6, 7, 10, -1, -1, -1, -1,&
10, 6, 7, 1, 10, 7, 1, 7, 8, 1, 8, 0, -1, -1, -1, -1,&
10, 6, 7, 10, 7, 1, 1, 7, 3, -1, -1, -1, -1, -1, -1, -1,&
1, 2, 6, 1, 6, 8, 1, 8, 9, 8, 6, 7, -1, -1, -1, -1,&
2, 6, 9, 2, 9, 1, 6, 7, 9, 0, 9, 3, 7, 3, 9, -1,&
7, 8, 0, 7, 0, 6, 6, 0, 2, -1, -1, -1, -1, -1, -1, -1,&
7, 3, 2, 6, 7, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1/),(/16,30/))




integer, dimension (16,30):: b5(16,30)=reshape&
((/2, 3, 11, 10, 6, 8, 10, 8, 9, 8, 6, 7, -1, -1, -1, -1,&
2, 0, 7, 2, 7, 11, 0, 9, 7, 6, 7, 10, 9, 10, 7, -1,&
1, 8, 0, 1, 7, 8, 1, 10, 7, 6, 7, 10, 2, 3, 11, -1,&
11, 2, 1, 11, 1, 7, 10, 6, 1, 6, 7, 1, -1, -1, -1, -1,&
8, 9, 6, 8, 6, 7, 9, 1, 6, 11, 6, 3, 1, 3, 6, -1,&
0, 9, 1, 11, 6, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
7, 8, 0, 7, 0, 6, 3, 11, 0, 11, 6, 0, -1, -1, -1, -1,&
7, 11, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
3, 0, 8, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
0, 1, 9, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
8, 1, 9, 8, 3, 1, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1,&
10, 1, 2, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
1, 2, 10, 3, 0, 8, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1,&
2, 9, 0, 2, 10, 9, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1,&
6, 11, 7, 2, 10, 3, 10, 8, 3, 10, 9, 8, -1, -1, -1, -1,&
7, 2, 3, 6, 2, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
7, 0, 8, 7, 6, 0, 6, 2, 0, -1, -1, -1, -1, -1, -1, -1,&
2, 7, 6, 2, 3, 7, 0, 1, 9, -1, -1, -1, -1, -1, -1, -1,&
1, 6, 2, 1, 8, 6, 1, 9, 8, 8, 7, 6, -1, -1, -1, -1,&
10, 7, 6, 10, 1, 7, 1, 3, 7, -1, -1, -1, -1, -1, -1, -1,&
10, 7, 6, 1, 7, 10, 1, 8, 7, 1, 0, 8, -1, -1, -1, -1,&
0, 3, 7, 0, 7, 10, 0, 10, 9, 6, 10, 7, -1, -1, -1, -1,&
7, 6, 10, 7, 10, 8, 8, 10, 9, -1, -1, -1, -1, -1, -1, -1,&
6, 8, 4, 11, 8, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
3, 6, 11, 3, 0, 6, 0, 4, 6, -1, -1, -1, -1, -1, -1, -1,&
8, 6, 11, 8, 4, 6, 9, 0, 1, -1, -1, -1, -1, -1, -1, -1,&
9, 4, 6, 9, 6, 3, 9, 3, 1, 11, 3, 6, -1, -1, -1, -1,&
6, 8, 4, 6, 11, 8, 2, 10, 1, -1, -1, -1, -1, -1, -1, -1,&
1, 2, 10, 3, 0, 11, 0, 6, 11, 0, 4, 6, -1, -1, -1, -1/),(/16,30/))



integer, dimension (16,30):: b6(16,30)=reshape&
((/4, 11, 8, 4, 6, 11, 0, 2, 9, 2, 10, 9, -1, -1, -1, -1,&
10, 9, 3, 10, 3, 2, 9, 4, 3, 11, 3, 6, 4, 6, 3, -1,&
8, 2, 3, 8, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1,&
0, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
1, 9, 0, 2, 3, 4, 2, 4, 6, 4, 3, 8, -1, -1, -1, -1,&
1, 9, 4, 1, 4, 2, 2, 4, 6, -1, -1, -1, -1, -1, -1, -1,&
8, 1, 3, 8, 6, 1, 8, 4, 6, 6, 10, 1, -1, -1, -1, -1,&
10, 1, 0, 10, 0, 6, 6, 0, 4, -1, -1, -1, -1, -1, -1, -1,&
4, 6, 3, 4, 3, 8, 6, 10, 3, 0, 3, 9, 10, 9, 3, -1,&
10, 9, 4, 6, 10, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
4, 9, 5, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
0, 8, 3, 4, 9, 5, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1,&
5, 0, 1, 5, 4, 0, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1,&
11, 7, 6, 8, 3, 4, 3, 5, 4, 3, 1, 5, -1, -1, -1, -1,&
9, 5, 4, 10, 1, 2, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1,&
6, 11, 7, 1, 2, 10, 0, 8, 3, 4, 9, 5, -1, -1, -1, -1,&
7, 6, 11, 5, 4, 10, 4, 2, 10, 4, 0, 2, -1, -1, -1, -1,&
3, 4, 8, 3, 5, 4, 3, 2, 5, 10, 5, 2, 11, 7, 6, -1,&
7, 2, 3, 7, 6, 2, 5, 4, 9, -1, -1, -1, -1, -1, -1, -1,&
9, 5, 4, 0, 8, 6, 0, 6, 2, 6, 8, 7, -1, -1, -1, -1,&
3, 6, 2, 3, 7, 6, 1, 5, 0, 5, 4, 0, -1, -1, -1, -1,&
6, 2, 8, 6, 8, 7, 2, 1, 8, 4, 8, 5, 1, 5, 8, -1,&
9, 5, 4, 10, 1, 6, 1, 7, 6, 1, 3, 7, -1, -1, -1, -1,&
1, 6, 10, 1, 7, 6, 1, 0, 7, 8, 7, 0, 9, 5, 4, -1,&
4, 0, 10, 4, 10, 5, 0, 3, 10, 6, 10, 7, 3, 7, 10, -1,&
7, 6, 10, 7, 10, 8, 5, 4, 10, 4, 8, 10, -1, -1, -1, -1,&
6, 9, 5, 6, 11, 9, 11, 8, 9, -1, -1, -1, -1, -1, -1, -1,&
3, 6, 11, 0, 6, 3, 0, 5, 6, 0, 9, 5, -1, -1, -1, -1,&
0, 11, 8, 0, 5, 11, 0, 1, 5, 5, 6, 11, -1, -1, -1, -1,&
6, 11, 3, 6, 3, 5, 5, 3, 1, -1, -1, -1, -1, -1, -1, -1/),(/16,30/))


integer, dimension (16,30):: b7(16,30)=reshape&
((/1, 2, 10, 9, 5, 11, 9, 11, 8, 11, 5, 6, -1, -1, -1, -1,&
0, 11, 3, 0, 6, 11, 0, 9, 6, 5, 6, 9, 1, 2, 10, -1,&
11, 8, 5, 11, 5, 6, 8, 0, 5, 10, 5, 2, 0, 2, 5, -1,&
6, 11, 3, 6, 3, 5, 2, 10, 3, 10, 5, 3, -1, -1, -1, -1,&
5, 8, 9, 5, 2, 8, 5, 6, 2, 3, 8, 2, -1, -1, -1, -1,&
9, 5, 6, 9, 6, 0, 0, 6, 2, -1, -1, -1, -1, -1, -1, -1,&
1, 5, 8, 1, 8, 0, 5, 6, 8, 3, 8, 2, 6, 2, 8, -1,&
1, 5, 6, 2, 1, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
1, 3, 6, 1, 6, 10, 3, 8, 6, 5, 6, 9, 8, 9, 6, -1,&
10, 1, 0, 10, 0, 6, 9, 5, 0, 5, 6, 0, -1, -1, -1, -1,&
0, 3, 8, 5, 6, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
10, 5, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
11, 5, 10, 7, 5, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
11, 5, 10, 11, 7, 5, 8, 3, 0, -1, -1, -1, -1, -1, -1, -1,&
5, 11, 7, 5, 10, 11, 1, 9, 0, -1, -1, -1, -1, -1, -1, -1,&
10, 7, 5, 10, 11, 7, 9, 8, 1, 8, 3, 1, -1, -1, -1, -1,&
11, 1, 2, 11, 7, 1, 7, 5, 1, -1, -1, -1, -1, -1, -1, -1,&
0, 8, 3, 1, 2, 7, 1, 7, 5, 7, 2, 11, -1, -1, -1, -1,&
9, 7, 5, 9, 2, 7, 9, 0, 2, 2, 11, 7, -1, -1, -1, -1,&
7, 5, 2, 7, 2, 11, 5, 9, 2, 3, 2, 8, 9, 8, 2, -1,&
2, 5, 10, 2, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1,&
8, 2, 0, 8, 5, 2, 8, 7, 5, 10, 2, 5, -1, -1, -1, -1,&
9, 0, 1, 5, 10, 3, 5, 3, 7, 3, 10, 2, -1, -1, -1, -1,&
9, 8, 2, 9, 2, 1, 8, 7, 2, 10, 2, 5, 7, 5, 2, -1,&
1, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
0, 8, 7, 0, 7, 1, 1, 7, 5, -1, -1, -1, -1, -1, -1, -1,&
9, 0, 3, 9, 3, 5, 5, 3, 7, -1, -1, -1, -1, -1, -1, -1,&
9, 8, 7, 5, 9, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
5, 8, 4, 5, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1,&
5, 0, 4, 5, 11, 0, 5, 10, 11, 11, 3, 0, -1, -1, -1, -1/),(/16,30/))


integer, dimension (16,30):: b8(16,30)=reshape&
((/0, 1, 9, 8, 4, 10, 8, 10, 11, 10, 4, 5, -1, -1, -1, -1,&
10, 11, 4, 10, 4, 5, 11, 3, 4, 9, 4, 1, 3, 1, 4, -1,&
2, 5, 1, 2, 8, 5, 2, 11, 8, 4, 5, 8, -1, -1, -1, -1,&
0, 4, 11, 0, 11, 3, 4, 5, 11, 2, 11, 1, 5, 1, 11, -1,&
0, 2, 5, 0, 5, 9, 2, 11, 5, 4, 5, 8, 11, 8, 5, -1,&
9, 4, 5, 2, 11, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
2, 5, 10, 3, 5, 2, 3, 4, 5, 3, 8, 4, -1, -1, -1, -1,&
5, 10, 2, 5, 2, 4, 4, 2, 0, -1, -1, -1, -1, -1, -1, -1,&
3, 10, 2, 3, 5, 10, 3, 8, 5, 4, 5, 8, 0, 1, 9, -1,&
5, 10, 2, 5, 2, 4, 1, 9, 2, 9, 4, 2, -1, -1, -1, -1,&
8, 4, 5, 8, 5, 3, 3, 5, 1, -1, -1, -1, -1, -1, -1, -1,&
0, 4, 5, 1, 0, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
8, 4, 5, 8, 5, 3, 9, 0, 5, 0, 3, 5, -1, -1, -1, -1,&
9, 4, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
4, 11, 7, 4, 9, 11, 9, 10, 11, -1, -1, -1, -1, -1, -1, -1,&
0, 8, 3, 4, 9, 7, 9, 11, 7, 9, 10, 11, -1, -1, -1, -1,&
1, 10, 11, 1, 11, 4, 1, 4, 0, 7, 4, 11, -1, -1, -1, -1,&
3, 1, 4, 3, 4, 8, 1, 10, 4, 7, 4, 11, 10, 11, 4, -1,&
4, 11, 7, 9, 11, 4, 9, 2, 11, 9, 1, 2, -1, -1, -1, -1,&
9, 7, 4, 9, 11, 7, 9, 1, 11, 2, 11, 1, 0, 8, 3, -1,&
11, 7, 4, 11, 4, 2, 2, 4, 0, -1, -1, -1, -1, -1, -1, -1,&
11, 7, 4, 11, 4, 2, 8, 3, 4, 3, 2, 4, -1, -1, -1, -1,&
2, 9, 10, 2, 7, 9, 2, 3, 7, 7, 4, 9, -1, -1, -1, -1,&
9, 10, 7, 9, 7, 4, 10, 2, 7, 8, 7, 0, 2, 0, 7, -1,&
3, 7, 10, 3, 10, 2, 7, 4, 10, 1, 10, 0, 4, 0, 10, -1,&
1, 10, 2, 8, 7, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
4, 9, 1, 4, 1, 7, 7, 1, 3, -1, -1, -1, -1, -1, -1, -1,&
4, 9, 1, 4, 1, 7, 0, 8, 1, 8, 7, 1, -1, -1, -1, -1,&
4, 0, 3, 7, 4, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
4, 8, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1/),(/16,30/))



integer, dimension (16,16):: b9(16,16)=reshape&
((/9, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
3, 0, 9, 3, 9, 11, 11, 9, 10, -1, -1, -1, -1, -1, -1, -1,&
0, 1, 10, 0, 10, 8, 8, 10, 11, -1, -1, -1, -1, -1, -1, -1,&
3, 1, 10, 11, 3, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
1, 2, 11, 1, 11, 9, 9, 11, 8, -1, -1, -1, -1, -1, -1, -1,&
3, 0, 9, 3, 9, 11, 1, 2, 9, 2, 11, 9, -1, -1, -1, -1,&
0, 2, 11, 8, 0, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
3, 2, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
2, 3, 8, 2, 8, 10, 10, 8, 9, -1, -1, -1, -1, -1, -1, -1,&
9, 10, 2, 0, 9, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
2, 3, 8, 2, 8, 10, 0, 1, 8, 1, 10, 8, -1, -1, -1, -1,&
1, 10, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
1, 3, 8, 9, 1, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
0, 9, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
0, 3, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1/),(/16,16/))

  do j=1,30
    do i=1,16
      a(j,i)=b1(i,j)
    end do 
  end do

  jj=0
  do j=31,60
    jj=jj+1
      do i=1,16
        a(j,i)=b2(i,jj)
      end do 
  end do

  jj=0
  do j=61,90
    jj=jj+1
      do i=1,16
        a(j,i)=b3(i,jj)
      end do 
  end do

  jj=0
  do j=91,120
    jj=jj+1
      do i=1,16
        a(j,i)=b4(i,jj)
      end do 
  end do

  jj=0
  do j=121,150
    jj=jj+1
      do i=1,16
        a(j,i)=b5(i,jj)
      end do 
  end do

  jj=0
  do j=151,180
    jj=jj+1
      do i=1,16
        a(j,i)=b6(i,jj)
      end do 
  end do

  jj=0
  do j=181,210
    jj=jj+1
      do i=1,16
        a(j,i)=b7(i,jj)
      end do 
  end do

  jj=0
  do j=211,240
    jj=jj+1
      do i=1,16
        a(j,i)=b8(i,jj)
      end do 
  end do

  jj=0
  do j=241,257
    jj=jj+1
      do i=1,16
        a(j,i)=b9(i,jj)
      end do 
  end do

  
  !allocating initial space for coordinate and triangle index arrays.These arrays will be extended if needed during the execution of the program using allocate and deallocate functions.
  nsize=1000
  allocate(x_coord(nsize),y_coord(nsize),z_coord(nsize),t_index(5*nsize))

  !marching through the grid and computing "sum" for each individual cube. "l" is the contor for the vectors that keeps the coordinates of intersection points. we will also produce a vector called t_index that tells which points build which triangles.
  l=0
  contor=0
  do k=1,Nz-1
    do j=1,Ny-1
      do i=1,Nx-1                                 
        sum=1
          if (grid_value(i,j,k)<0) then 
            sum=sum+2**0
          end if
          if (grid_value(i+1,j,k)<0) then 
            sum=sum+2**3
          end if
          if (grid_value(i+1,j+1,k)<0) then 
            sum=sum+2**2
          end if
          if (grid_value(i,j+1,k)<0) then 
            sum=sum+2**1
          end if
          if (grid_value(i+1,j+1,k+1)<0) then 
            sum=sum+2**6
          end if
          if (grid_value(i+1,j,k+1)<0) then 
            sum=sum+2**7
          end if
          if (grid_value(i,j,k+1)<0) then 
            sum=sum+2**4
          end if
          if (grid_value(i,j+1,k+1)<0) then 
            sum=sum+2**5
          end if

          if (a(sum,1)>-1) then
            do jj=1,15,3
              if (a(sum,jj)>-1) then
                index2=0
                  do kk=jj,jj+2
                    index2=index2+1
                    if (a(sum,kk)==0) then 
                      u=abs(grid_value(i,j,k))*(grid_y(j+1)-grid_y(j))/(abs(grid_value(i,j,k))+abs(grid_value(i,j+1,k)))
                      temp_x(index2)=grid_x(i)
                      temp_y(index2)=grid_y(j)+u
                      temp_z(index2)=grid_z(k)
                    end if
 
                    if (a(sum,kk)==1) then 
                      u=abs(grid_value(i,j+1,k))*(grid_x(i+1)-grid_x(i))/(abs(grid_value(i,j+1,k))+abs(grid_value(i+1,j+1,k)))
                      temp_x(index2)=grid_x(i)+u
                      temp_y(index2)=grid_y(j+1)
                      temp_z(index2)=grid_z(k)            
                    end if

                    if (a(sum,kk)==2) then 
                      u=abs(grid_value(i+1,j,k))*(grid_y(j+1)-grid_y(j))/(abs(grid_value(i+1,j+1,k))+abs(grid_value(i+1,j,k)))
                      temp_x(index2)=grid_x(i+1)
                      temp_y(index2)=grid_y(j)+u
                      temp_z(index2)=grid_z(k)
                    end if

                    if (a(sum,kk)==3) then
                      u=abs(grid_value(i,j,k))*(grid_x(i+1)-grid_x(i))/(abs(grid_value(i+1,j,k))+abs(grid_value(i,j,k)))
                      temp_x(index2)=grid_x(i)+u
                      temp_y(index2)=grid_y(j)
                      temp_z(index2)=grid_z(k)
                    end if

                    if (a(sum,kk)==4) then 
                      u=abs(grid_value(i,j,k+1))*(grid_y(j+1)-grid_y(j))/(abs(grid_value(i,j,k+1))+abs(grid_value(i,j+1,k+1)))
                      temp_x(index2)=grid_x(i)
                      temp_y(index2)=grid_y(j)+u
                      temp_z(index2)=grid_z(k+1)
                    end if
 
                    if (a(sum,kk)==5) then 
                      u=abs(grid_value(i,j+1,k+1))*(grid_x(i+1)-grid_x(i))/(abs(grid_value(i,j+1,k+1))+abs(grid_value(i+1,j+1,k+1)))
                      temp_x(index2)=grid_x(i)+u
                      temp_y(index2)=grid_y(j+1)
                      temp_z(index2)=grid_z(k+1)        
                    end if

                    if (a(sum,kk)==6) then 
                      u=abs(grid_value(i+1,j,k+1))*(grid_y(j+1)-grid_y(j))/(abs(grid_value(i+1,j+1,k+1))+abs(grid_value(i+1,j,k+1)))
                      temp_x(index2)=grid_x(i+1)
                      temp_y(index2)=grid_y(j)+u
                      temp_z(index2)=grid_z(k+1)
                    end if
  
                    if (a(sum,kk)==7) then
                      u=abs(grid_value(i,j,k+1))*(grid_x(i+1)-grid_x(i))/(abs(grid_value(i,j,k+1))+abs(grid_value(i+1,j,k+1)))
                      temp_x(index2)=grid_x(i)+u
                      temp_y(index2)=grid_y(j)
                      temp_z(index2)=grid_z(k+1)           
                    end if

                    if (a(sum,kk)==8) then 
                      u=abs(grid_value(i,j,k))*(grid_z(k+1)-grid_z(k))/(abs(grid_value(i,j,k))+abs(grid_value(i,j,k+1)))
                      temp_x(index2)=grid_x(i)
                      temp_y(index2)=grid_y(j)
                      temp_z(index2)=grid_z(k)+u
                    end if

                    if (a(sum,kk)==9) then 
                      u=abs(grid_value(i,j+1,k))*(grid_z(k+1)-grid_z(k))/(abs(grid_value(i,j+1,k))+abs(grid_value(i,j+1,k+1)))
                      temp_x(index2)=grid_x(i)
                      temp_y(index2)=grid_y(j+1)
                      temp_z(index2)=grid_z(k)+u
                    end if

                    if (a(sum,kk)==10) then 
                      u=abs(grid_value(i+1,j+1,k))*(grid_z(k+1)-grid_z(k))/(abs(grid_value(i+1,j+1,k))+abs(grid_value(i+1,j+1,k+1)))
                      temp_x(index2)=grid_x(i+1)
                      temp_y(index2)=grid_y(j+1)
                      temp_z(index2)=grid_z(k)+u
                    end if

                    if (a(sum,kk)==11) then 
                      u=abs(grid_value(i+1,j,k))*(grid_z(k+1)-grid_z(k))/(abs(grid_value(i+1,j,k))+abs(grid_value(i+1,j,k+1)))
                      temp_x(index2)=grid_x(i+1)
                      temp_y(index2)=grid_y(j)
                      temp_z(index2)=grid_z(k)+u
                    end if
                  end do

                ABx=temp_x(3)-temp_x(2)
                BCx=temp_x(2)-temp_x(1)
                ABy=temp_y(3)-temp_y(2)
                BCy=temp_y(2)-temp_y(1)
                ABz=temp_z(3)-temp_z(2)
                BCz=temp_z(2)-temp_z(1)
                equality=1
                if ((abs((ABx*BCy)-(ABy*BCx))<1e-10).and.(abs((ABx*BCz)-(ABz*BCx))<1e-10).and.(abs((ABy*BCz)-(ABz*BCy))<1e-10)) then 
                  equality=0
                end if
                if (equality==1) then
                  boolean_finder1=0
                  boolean_finder2=0
                  boolean_finder3=0
                  if (l>1) then
                    do t=1,l
                      if ((abs(temp_x(1)-x_coord(t))<1e-10).and.(abs(temp_y(1)-y_coord(t))<1e-10).and.(abs(temp_z(1)-z_coord(t))<1e-10)) then
                        boolean_finder1=1
                        contor=contor+1
                        t_index(contor)=t-1
                      end if

                      if ((abs(temp_x(2)-x_coord(t))<1e-10).and.(abs(temp_y(2)-y_coord(t))<1e-10).and.(abs(temp_z(2)-z_coord(t))<1e-10)) then
                        boolean_finder2=1
                        contor=contor+1
                        t_index(contor)=t-1
                      end if
                      if ((abs(temp_x(3)-x_coord(t))<1e-10).and.(abs(temp_y(3)-y_coord(t))<1e-10).and.(abs(temp_z(3)-z_coord(t))<1e-10)) then
                        boolean_finder3=1
                        contor=contor+1
                        t_index(contor)=t-1
                      end if
                    end do
                  end if
                  if (boolean_finder1==0) then
                    l=l+1
                    x_coord(l)=temp_x(1)
                    y_coord(l)=temp_y(1)
                    z_coord(l)=temp_z(1)
                    contor=contor+1
                    t_index(contor)=l-1
                  end if
                  if (boolean_finder2==0) then
                    l=l+1
                    x_coord(l)=temp_x(2)
                    y_coord(l)=temp_y(2)
                    z_coord(l)=temp_z(2)
                    contor=contor+1
                    t_index(contor)=l-1
                  end if
                  if (boolean_finder3==0) then
                    l=l+1
                    x_coord(l)=temp_x(3)
                    y_coord(l)=temp_y(3)
                    z_coord(l)=temp_z(3)
                    contor=contor+1
                    t_index(contor)=l-1
                  end if
                end if
               

              if ((nsize-l)<20) then
                allocate(resize_coord_x(nsize),resize_coord_y(nsize), &
                         resize_coord_z(nsize),resize_t_index(nsize*5))
                resize_coord_x(1:l)=x_coord(1:l)
                resize_coord_y(1:l)=y_coord(1:l)
                resize_coord_z(1:l)=z_coord(1:l)
                resize_t_index(1:contor)=t_index(1:contor)
                deallocate(x_coord,y_coord,z_coord,t_index)
                nsize=nsize*2
                allocate(x_coord(nsize),y_coord(nsize),z_coord(nsize),t_index(nsize*5))
                x_coord(1:l)=resize_coord_x(1:l)
                y_coord(1:l)=resize_coord_y(1:l)
                z_coord(1:l)=resize_coord_z(1:l)
                t_index(1:contor)=resize_t_index(1:contor)
                deallocate(resize_coord_x,resize_coord_y,resize_coord_z,resize_t_index)
              end if
         
              end if
            end do
          end if
        end do
      end do
    end do 
!writing data to the files-temporary output method
  open (unit=9,file='output_x.txt',status='replace',action='write')
  open (unit=10,file='output_y.txt',status='replace',action='write')
  open (unit=11,file='output_z.txt',status='replace',action='write')
  open (unit=12,file='triangle_test.txt',status='replace',action='write')  
  n_triangles=contor
  n_coords=l
  do contor=1,n_triangles
    write(12,*)t_index(contor) 
  end do
  do l=1,n_coords
    write(9,*)x_coord(l)
    write(10,*)y_coord(l)
    write(11,*)z_coord(l)
  end do
  close(9)
  close(10)
  close(11)
  close(12)
  print*,'t_index = ',t_index(1:n_triangles)
  end subroutine triangularization
  end module

module interpa
contains
subroutine interpolation_area(cctkGH,npoints,x_coord,y_coord,z_coord,l,t_index,area)
  
  use areac
  use EHFinder_mod
  implicit none
  DECLARE_CCTK_PARAMETERS
!  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS 

  CCTK_POINTER, intent(IN) :: cctkGH
  CCTK_INT, intent(IN) :: npoints, l
  CCTK_REAL, dimension(:), pointer :: x_coord,y_coord,z_coord
  CCTK_INT, dimension(:), pointer :: t_index
  CCTK_REAL, intent(OUT) :: area
  character(len=200) :: area_interp
  character(len=128) :: warn_message
  CCTK_INT :: area_interp_len
  character(len=7) :: area_order
  CCTK_INT :: interp_handle, table_handle, coord_system_handle
  CCTK_POINTER, dimension(3) :: interp_coords
  CCTK_POINTER, dimension(6) :: out_array
  CCTK_INT, dimension(6) :: in_array
  CCTK_REAL, dimension(:), allocatable :: metric_xx, metric_yy, metric_zz, metric_xy, metric_yz, metric_xz
  
  CCTK_INT, dimension(6), parameter :: out_types = (/ CCTK_VARIABLE_REAL, &
                                                      CCTK_VARIABLE_REAL, &
                                                      CCTK_VARIABLE_REAL, &
                                                      CCTK_VARIABLE_REAL, &
                                                      CCTK_VARIABLE_REAL, &
                                                      CCTK_VARIABLE_REAL /)
  ! Convert the interpolator parameter into a fortran string.
  call CCTK_FortranString ( area_interp_len, area_interpolator, &
                                             area_interp )

!allocation of metric arrays
  allocate(metric_xx(npoints), metric_yy(npoints), metric_zz(npoints), metric_xy(npoints), metric_xz(npoints), metric_yz(npoints))

! Get an interpolation handle.
  call CCTK_InterpHandle ( interp_handle, area_interp(1:area_interp_len) )

  if ( interp_handle .lt. 0 ) then
    warn_message = 'Cannot get handle for interpolation. '
    warn_message = trim(warn_message)//'Forgot to activate an implementation '
    warn_message = trim(warn_message)//'providing interpolation operators?'
    call CCTK_WARN( 0, trim(warn_message) )
  end if

! Write the interpolation order parameter into a string
  write(area_order,'(a6,i1)') 'order=',area_interpolation_order

! Create a table from the interpolation order string.
  call Util_TableCreateFromString ( table_handle, area_order )
  if ( table_handle .lt. 0 ) then
    call CCTK_WARN( 0, 'Cannot create parameter table for interpolator' )
  end if

! Get a handle for the coordinate system.
  call CCTK_CoordSystemHandle ( coord_system_handle, 'cart3d' )
  if ( coord_system_handle .lt. 0) then
    warn_message = 'Cannot get handle for cart3d coordinate system. '
    warn_message = trim(warn_message)//'Forgot to activate an implementation '
    warn_message = trim(warn_message)//'providing coordinates?'
    call CCTK_WARN( 0, trim(warn_message) )
  endif

! Get the pointers to the interpolation points.
  interp_coords(1) = CCTK_PointerTo(x_coord)
  interp_coords(2) = CCTK_PointerTo(y_coord)
  interp_coords(3) = CCTK_PointerTo(z_coord)

! Get the pointers to the interpolation return arrays.
  out_array(1) = CCTK_PointerTo(metric_xx)
  out_array(2) = CCTK_PointerTo(metric_xy)
  out_array(3) = CCTK_PointerTo(metric_xz)
  out_array(4) = CCTK_PointerTo(metric_yy)
  out_array(5) = CCTK_PointerTo(metric_yz)
  out_array(6) = CCTK_PointerTo(metric_zz)
!  out_array(7) = CCTK_PointerTo(psii)
! Get the indices for the grid functions to be interpolated.
  call CCTK_VarIndex ( in_array(1), 'admbase::gxx' )
  call CCTK_VarIndex ( in_array(2), 'admbase::gxy' )
  call CCTK_VarIndex ( in_array(3), 'admbase::gxz' )
  call CCTK_VarIndex ( in_array(4), 'admbase::gyy' )
  call CCTK_VarIndex ( in_array(5), 'admbase::gyz' )
  call CCTK_VarIndex ( in_array(6), 'admbase::gzz' )

  call CCTK_InterpGridArrays ( ierr, cctkGH, 3, interp_handle, &
                                 table_handle, coord_system_handle, &
                                 npoints, CCTK_VARIABLE_REAL, &
                                 interp_coords, 6, in_array(1:6), &
                                 6, out_types(1:6), out_array(1:6) )
  
  call area_computation(x_coord, y_coord, z_coord, t_index, &
                        metric_xx, metric_xy, metric_xz, &
                        metric_yy, metric_yz, metric_zz, &
                        npoints, l, area )

  deallocate(metric_xx, metric_yy, metric_zz, metric_xz, metric_yz, metric_xy)

  end subroutine interpolation_area
  end module

subroutine EHFinder_IsoSurface(CCTK_ARGUMENTS)

  use interpa
  use EHFinder_mod
  use IsoSurface_mod
  use triangular

  implicit none

  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS

  CCTK_REAL, dimension(:), pointer :: x_coord, y_coord, z_coord
  CCTK_INT, dimension(:), pointer :: t_index
  CCTK_INT :: l,contor,p,nsize
  
! Find index ranges for points excluding ghost zones and symmetry points for
! the 3D grid functions.
#include "include/physical_part.h"

  call CCTK_INFO ('Entered IsoSurface')
  
  nx = cctk_lsh(1); ny = cctk_lsh(2); nz = cctk_lsh(3)

  nsize=1000
  allocate(x_coord(nsize),y_coord(nsize),z_coord(nsize),t_index(5*nsize))

  do p = 1, eh_number_level_sets
    call triangularization(l,contor,nx,ny,nz,x_coord,y_coord,z_coord, &
                           t_index,f(:,:,:,p),x(:,1,1),y(1,:,1),z(1,1,:))

    print*, 'contor = ', t_index(1:contor)
    print*,'The surface consists of', contor/3,'triangles and',l,'coordinates'

    call interpolation_area(cctkGH, l, x_coord,y_coord,z_coord, contor, &
                            t_index, eh_area(1,p) )

    print*,'The area is ', eh_area(1,p)   
  end do

  

end subroutine EHFinder_IsoSurface

  
