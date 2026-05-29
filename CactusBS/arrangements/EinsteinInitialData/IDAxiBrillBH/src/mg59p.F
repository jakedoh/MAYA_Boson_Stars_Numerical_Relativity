#include "cctk.h"
c----------------------------------------------------------------------
      subroutine mgparm(m,ifd59,id5,id9,idi,idg,imx,jmx)
      implicit CCTK_REAL (a-h,o-z)
c----------------------------------------------------------------------
c  Given imx, jmx and ifd59 (See comments in mgsu2), mgparm calculates
c  the number of grids that will be needed, and the dimensions
*  needed for the coefficient, right hand side and solution arrays
*  to store values on all grid levels.
c .....................................................................
*
****  Parameters
*
*  INTEGERS:
*     m :
*         This is the number of grid levels that the multigrid
*         routine will use.
*
*     ifmg :
*          ifmg = 0 - The full multigrid algorithm is not used to
*                     obtain a good initial guess on the fine grid.
*                     (use this if you can provide a good initial guess)
*          ifmg = 1 - The full multigrid algorithm is used to obtain a
*                     good initial guess on the fine grid.
*
*     id5 :
*          Dimension of the arrays ac,aw,as,ae,an,q and f. id5 is the
*          total number of grid points on the finest grid and all
*          coarser grids.
*
*     id9 :
*          Dimension of the arrays asw,ase,ane,anw. If ifd59=5 then
*          id9=idi.  If ifd59=9 then id9=id5.
*          (NOTE: This routine specifically written for 9-point)
*
*     idi :
*          Dimension of the work arrays pu and pd. idi is the total
*          number of grid points on all of the coarser grids.
*
*     idg :
*          Dimension of the work array gam. It is set to the value im,
*          the number of grid points in the i-direction on the finest
*          grid.
*
*     imx,jmx :
*          These are the number of points in the i and j directions
*          including the two ficticious points.
*
      parameter(n5=1,n9=2,ni=3,jm=4,i9=5,j9=6,ifd=7,jred=8)
      dimension np2(20,8)
      iq5=1
      iq9=1
      iqi=1
      m=1
      np2(m,1)=jmx
      np2(m,2)=3
   10 if(np2(m,1).le.3) go to 20
      m=m+1
      np2(m,1)=np2(m-1,1)/2+1
      if(np2(m-1,2).eq.2.and.mod(np2(m-1,1),2).eq.1)
     + np2(m,1)=np2(m,1)+1
      np2(m,2)=2
      go to 10
   20 do 30 k=1,m
      np2(m-k+1,jm)=np2(k,1)
   30 np2(m-k+1,jred)=np2(k,2)
      do 40 k=m,1,-1
      ktot=imx*np2(k,jm)
      np2(k,n5)=iq5
      iq5=iq5+ktot
      np2(k,n9)=iq9
      if(k.lt.m.or.ifd59.eq.9) iq9=iq9+ktot
      np2(k,ni)=iqi
   40 if(k.lt.m) iqi=iqi+ktot
      do 50 k=1,m
      np2(k,i9)=imx
      np2(k,j9)=np2(k,jm)
   50 np2(k,ifd)=9
         if(ifd59.eq.5) then
      np2(m,i9)=1
      np2(m,j9)=1
      np2(m,ifd)=5
         endif
      id5=iq5-1
      id9=iq9-1
      idi=iqi-1
      idg=imx
      return
      end
*
*
*
************************************************************************
*
      subroutine mg9 (idim,ilower,iupper,jdim,jlower,jupper,
     &             cc,cn,cs,cw,ce,cnw,cne,csw,cse,u,rhs,
     &             id5,id9,idi,idg,ifmg,eps,rmax,ier)
      implicit CCTK_REAL (a-h,o-z)
*
************************************************************************
*
*
*     This routine is a wrapper for the multigrid solver.  The input
*  arrays are the finite difference coefficients on the 2D grid with
*  the boundary conditions absorbed into them.
*
****  Parameters
*
*  INTEGERS:
*     idim,jdim :
*          These define sizes of the coefficient arrays passed to this
*          routine.
*
*     ilower,iupper,
*     jlower,jupper :
*          These are the indices of the computational grid that
*          correspond to upper and lower index limits for points on
*          which computations are actually done.
*
*     id5 :
*          Dimension of the arrays ac,aw,as,ae,an,q and f. id5 is the
*          total number of grid points on the finest grid and all
*          coarser grids.
*
*     id9 :
*          Dimension of the arrays asw,ase,ane,anw. If ifd59=5 then
*          id9=idi.  If ifd59=9 then id9=id5.
*          (NOTE: This routine specifically written for 9-point)
*
*     idi :
*          Dimension of the work arrays pu and pd. idi is the total
*          number of grid points on all of the coarser grids.
*
*     idg :
*          Dimension of the work array gam. It is set to the value im,
*          the number of grid points in the i-direction on the finest
*          grid.
*
*     ifmg :
*          ifmg = 0 - The full multigrid algorithm is not used to
*                     obtain a good initial guess on the fine grid.
*                     (use this if you can provide a good initial guess)
*          ifmg = 1 - The full multigrid algorithm is used to obtain a
*                     good initial guess on the fine grid.
*
*     ier :
*          This is an error flag with the following possible return
*          values:
*              ier =  0 => solver converged without error
*              ier = -1 => solver did not converge
*
*  REALS:
*     cc(idim,jdim),cn(idim,jdim),
*     cs(idim,jdim),cw(idim,jdim),ce(idim,jdim),
*     cnw(idim,jdim),cne(idim,jdim),
*     csw(idim,jdim),cse(idim,jdim) :
*          These are the finite difference coefficient arrays for a
*          nine-point stencil on the two dimensional grid as follows:
*
*
*               cnw    cn    cne
*            ^
*            |  cw     cc    ce
* increasing |
*  j values  |  csw    cs    cse
*   (theta)  |
*             --------> increasing i values (eta)
*
*  Ie. the coresspondence is : nw : i-1,j+1
*                              n  : i,j+1
*                              ne : i+1,j+1
*                              w  : i-1,j
*                              c  : i,j
*                              e  : i+1,j
*                              sw : i-1,j-1
*                              s  : i,j-1
*                              se : i+1,j-1
*
*     u :
*          Input: this contains the initial guess to the solution of the
*                 equation
*          Output: This contains the final approximation to the solution
*                 determined by the multigrid solver.
*     
*     rhs :
*          This array contains the values of the right hand side of the
*          equation at every point on the rwo dimensional grid.
*
*     eps :
*          eps > 0.  The maximum norm of the residual is calculated at
*                    the end of each multigrid cycle. The algorithm is
*                    terminated when this maximum becomes less than eps
*                    or when the maximum number of iterations (see ncyc)
*                    is exceeded.  It is up to the user to provide a
*                    meaningfull tolerance criteria for the particular
*                    problem being solved.
*
*     rmax:
*          This is the final value of the residual calcu;ated.
*
************************************************************************
*
*

      integer idim,ilower,iupper,jdim,jlower,jupper,ifmg
      integer id5,id9,idi,idg,ier
      CCTK_REAL cc(idim,jdim),cn(idim,jdim),cs(idim,jdim),cw(idim,jdim),
     &     ce(idim,jdim),cnw(idim,jdim),cne(idim,jdim),csw(idim,jdim),
     &     cse(idim,jdim),u(idim,jdim),rhs(idim,jdim)
      CCTK_REAL eps,rmax

*
************************************************************************
*
*  Variable definitions:
*
*  Integers:
*     ifd59 :
*          ifd59 = 5 - means a 5-point finite difference stencil
*                      (ac,an,as,aw,ae) is defined on the finest grid.
*          ifd59 = 9 - means a 9-point finite difference stencil
*                      (ac,an,as,aw,ae,anw,ane,asw,ase) is defined on
*                      the finest grid by the user.
*                      (NOTE: This routine specifically written
*                             for 9-point)
*
*     ncyc :
*          The maximum number of multigrid "v"-cycles to be used. If
*          the maximum norm of the residual is not less than tol at the
*          end of ncyc cycles, the algorithm is terminated.
*          (NOTE: ncyc <= 40 )
*          
*     np(20,8) :
*          Input: When the iskip=1,-1 or -2 option is used, np2 is
*                 assumed to contain the grid information for umgs2.
*          Output: When the iskip=0 option is used, the grid
*                 information for umgs2 is returned in np2.
*          (NOTE: This is only useful for multiple instance problems)
*
*     iskip :
*          iskip = 0 - The coarse grid information, coarse grid
*                      operators and interpolation coefficients are
*                      calculated by umgs2.  This information is stored
*                      in the arrays ac, aw, as, asw, ase, pu, pd, np2
*                      and the variable.
*          iskip = 1 - The calculation of the coarse grid information,
*                      coarse grid operators and interpolation
*                      coefficients is skipped.  This option would be
*                      used when umgs2 has been called with iskip=0 and
*                      is being called again to solve a system of
*                      equations with the same matrix. This would be
*                      the case in, say, parabolic problems with time
*                      independent coefficients.
*          iskip =-1 - The set up of pointers (ugrdfn) is skipped.
*                      Coarse grid operators and interpolation
*                      coefficients are calculated and the given matrix
*                      equation is solved.  This option would be used
*                      when umgs2 has been called with iskip=0 and is
*                      being called again to solve a system of equations
*                      with a different matrix of the same dimensions.
*                      This would be the case for, say, parabolic
*                      problems with time dependent coefficients.
*          iskip =-2 - The set up of pointers (ugrdfn) is skipped.
*                      Coarse grid operators and interpolation
*                      coefficients are calculated and returned.
*                      No matrix solve.
*
*     ipc :
*          ipc = 0 or 1.
*          ipc is a multigrid parameter which determines the type of
*          interpolation to be used.  Usually ipc=1 is best.  However,
*          if the boundary condition equations have been absorbed into
*          the interior equations then ipc=0 can be used which results
*          in a slightly more efficient algorithm.
*
*     nman :
*          nman = 0 usually.
*          nman = 1 signals that the fine grid equations are singulari
*          for the case when homogeneous Neumann boundary conditions are
*          applied along the entire boundary.  In this case, the
*          difference equations are singular and the condition that the
*          integral of q over the domain be zero is added to the set of
*          difference equations.  This condition is satisfied by adding
*          the appropriate constant vector to q on the fine grid.  It is
*          assumed, in this case, that a well-defined problem has been
*          given to mgss2, i.e. the integral of f over the domain is
*          zero.
*
*     im :
*          The number of grid points in the x-direction (including two
*          ficticious points)
*     jm :
*          The number of grid points in the y-direction (including two
*          ficticious points)
*
*     linp :
*          This is a dummy argument left over from the authors
*          development of the code
*             Use:  common /io/ linp,lout
*
*     lout :
*          lout = unit number of output file into which the maximum norm
*          of the residual after each multigrid v-cycle" is printed.
*             Use:  common /io/ linp,lout
*
*     iscale :
*          Flag to indicate whether problem can be diagonally scaled to
*          speed convergence of the multigrid solver.
*
*     REALS:
*
*     ac(id5),an(id5),as(id5),aw(id5),ae(id5),
*     anw(id9),ane(id9),asw(id9),ase(id9) :
*          Input: ac, an, as, aw, ae, anw, ane, asw and ase contain the
*                 stencil coefficients for the difference operator on
*                 the finest grid. When the iskip=1 option is used,
*                 these arrays also are assumed to contain the coarse
*                 grid difference stencil coeficients.
*          Output: when the iskip=0 option is used, the coarse grid
*                 stencil coeficients are returned in ac, an, as, aw,
*                 ae, anw, ane, asw and ase.
*
*     ru(idi),rd(idi),rc(idi) :
*          Real work arrays.
*
*     pu(idi),pd(idi),pc(idi) :
*          Real work arrays.
*          Input: when the iskip=1 option is used, these arrays are
*                 assumed to contain the interpolation coefficients used
*                 in the semi-coarsening multigrid algorithm.
*          Output: when the iskip=0 option is used, the interpolation
*                 coeficients are returned in pu and pd.
*
*     f(id5) :
*          f contains the right hand side vector of the matrix
*          equation to be solved by umgs2.
*
*     q(id5) :
*          If ifmg=0, q contains the initial guess on the fine grid.
*          If ifmg=1, the initial guess on the fine grid is determined
*                     by the full multigrid process and the value of
*                     q on input to umgs2 not used.
*
*     tol :
*          tol > 0.  The maximum norm of the residual is calculated at
*                    the end of each multigrid cycle. The algorithm is
*                    terminated when this maximum becomes less than tol
*                    or when the maximum number of iterations (see ncyc)
*                    is exceeded.  It is up to the user to provide a
*                    meaningfull tolerance criteria for the particular
*                    problem being solved.
*          tol = 0.  Perform ncyc multigrid cycles.  Calculate and print
*                    the maximum norm of the residual after each cycle.
*          tol =-1.  Perform ncyc multigrid cycles.  The maximum norm of
*                    the final residual is calculated and returned in
*                    the variable rmax in the calling list of umgs2.
*          tol =-2.  Perform ncyc multigrid cycles.  The maximum norm of
*                    the residual is not calculated.
*
*     rmax :
*          If tol.ge.-1., the final residual norm is returned in rmax.
*
************************************************************************
*
*

      integer ncyc,ifd59
      parameter (ncyc=40,ifd59=9)

      integer np2(20,8)
      integer iskip,ipc,nman
      parameter (iskip=0,ipc=1,nman=0)
      integer irc,irurd,im,jm
      integer linp,lout
      common /io/ linp,lout

      CCTK_REAL ac(id5),an(id5),as(id5),aw(id5),ae(id5),
     &     anw(id9),ane(id9),asw(id9),ase(id9),
     &     q(id5),f(id5),gam(idg)

      CCTK_REAL ru(idi),rd(idi),rc(idi),pu(idi),pd(idi),pc(idi)

      CCTK_REAL tol

      integer iscale,itry


      ier=0

*
* Set some parameters for multigrid solver
*
      lout=6
c      rewind(unit=lout)
      irc=0
      irurd=0
      im=iupper-ilower+3
      jm=jupper-jlower+3
      tol=eps

*
*  Set up coefficients into vectors with correct indexing
*
      do 110 j=jlower,jupper
      do 100 i=ilower,iupper
      n=(j+1-jlower)*im + i+1-(ilower-1)
      ac(n)=cc(i,j)
      an(n)=cn(i,j)
      as(n)=cs(i,j)
      aw(n)=cw(i,j)
      ae(n)=ce(i,j)
      anw(n)=cnw(i,j)
      ane(n)=cne(i,j)
      asw(n)=csw(i,j)
      ase(n)=cse(i,j)
      q(n)=u(i,j)
      f(n)=rhs(i,j)
  100 continue
  110 continue


*
*  Determine whether we can diagonal scale the problem to speed
*  convergence. Can only be done if there are no zeros on the main
*  diagonal (ie. central difference coefficient).
*
      iscale=1
      do 200 j=jlower,jupper
      do 205 i=ilower,iupper
      n=(j+1-jlower)*im + i+1-(ilower-1)
      if (ac(n) .eq. 0.) then
        iscale=0
      endif
  205 continue
  200 continue

*
*  Do the diagonal scaling if we can.
*
      if (iscale.eq.1) then

      do 210 j=jlower,jupper
      do 215 i=ilower,iupper
        n=(j+1-jlower)*im + i+1-(ilower-1)
        f(n)=f(n)/ac(n)
        ase(n)=ase(n)/ac(n)
        asw(n)=asw(n)/ac(n)
        ane(n)=ane(n)/ac(n)
        anw(n)=anw(n)/ac(n)
        ae(n)=ae(n)/ac(n)
        aw(n)=aw(n)/ac(n)
        as(n)=as(n)/ac(n)
        an(n)=an(n)/ac(n)
        ac(n)=1.
  215 continue
  210 continue

      endif

c 
c   Now call the multigrid routine

      itry=0
      
 1122 call umgs2(
     . ac,aw,as,ae,an,asw,ase,ane,anw,q,f,pu,pd,pc,ru,rd,rc,gam,np2,
     . ifd59,ifmg,ncyc,tol,nman,im,jm,id5,id9,idi,m,iskip,rmax,
     . ipc,irc,irurd)

      if ((rmax.gt.tol).and.(itry.le.5)) then
        itry=itry+1
        print*,"Retry #",itry
        goto 1122
      endif

      if (rmax.gt.tol) then
        ier = -1
        print*,"Did not converge"
        print*," maximum residual    = ",rmax
        print*," tolerance           = ",tol
      endif

*
*  Convert the solution back to the 2D array form
*
      do 510 j=jlower,jupper
      do 500 i=ilower,iupper
      n=(j+1-jlower)*im + i+1-(ilower-1)
      u(i,j)=q(n)
  500 continue
  510 continue

      return
      end



*
*
*
************************************************************************
*
      subroutine mg5 (idim,ilower,iupper,jdim,jlower,jupper,
     &             cc,cn,cs,cw,ce,u,rhs,
     &             id5,id9,idi,idg,ifmg,eps,rmax,ier)
      implicit CCTK_REAL (a-h,o-z)
*
************************************************************************
*
*
*     This routine is a wrapper for the multigrid solver.  The input
*  arrays are the finite difference coefficients on the 2D grid with
*  the boundary conditions absorbed into them.
*
****  Parameters
*
*  INTEGERS:
*     idim,jdim :
*          These define sizes of the coefficient arrays passed to this
*          routine.
*
*     ilower,iupper,
*     jlower,jupper :
*          These are the indices of the computational grid that
*          correspond to upper and lower index limits for points on
*          which computations are actually done.
*
*     id5 :
*          Dimension of the arrays ac,aw,as,ae,an,q and f. id5 is the
*          total number of grid points on the finest grid and all
*          coarser grids.
*
*     id9 :
*          Dimension of the arrays asw,ase,ane,anw. If ifd59=5 then
*          id9=idi.  If ifd59=9 then id9=id5.
*          (NOTE: This routine specifically written for 5-point)
*
*     idi :
*          Dimension of the work arrays pu and pd. idi is the total
*          number of grid points on all of the coarser grids.
*
*     idg :
*          Dimension of the work array gam. It is set to the value im,
*          the number of grid points in the i-direction on the finest
*          grid.
*
*     ifmg :
*          ifmg = 0 - The full multigrid algorithm is not used to
*                     obtain a good initial guess on the fine grid.
*                     (use this if you can provide a good initial guess)
*          ifmg = 1 - The full multigrid algorithm is used to obtain a
*                     good initial guess on the fine grid.
*
*     ier :
*          This is an error flag with the following possible return
*          values:
*              ier =  0 => solver converged without error
*              ier = -1 => solver did not converge
*
*  REALS:
*     cc(idim,jdim),cn(idim,jdim),
*     cs(idim,jdim),cw(idim,jdim),ce(idim,jdim):
*          These are the finite difference coefficient arrays for a
*          five-point stencil on the two dimensional grid as follows:
*
*
*                      cn       
*            ^
*            |  cw     cc    ce
* increasing |
*  j values  |         cs       
*   (theta)  |
*             --------> increasing i values (eta)
*
*  Ie. the coresspondence is : n  : i,j+1
*                              ne : i+1,j+1
*                              c  : i,j
*                              e  : i+1,j
*                              s  : i,j-1
*
*     u :
*          Input: this contains the initial guess to the solution of the
*                 equation
*          Output: This contains the final approximation to the solution
*                 determined by the multigrid solver.
*     
*     rhs :
*          This array contains the values of the right hand side of the
*          equation at every point on the rwo dimensional grid.
*
*     eps :
*          eps > 0.  The maximum norm of the residual is calculated at
*                    the end of each multigrid cycle. The algorithm is
*                    terminated when this maximum becomes less than tol
*                    or when the maximum number of iterations (see ncyc)
*                    is exceeded.  It is up to the user to provide a
*                    meaningfull tolerance criteria for the particular
*                    problem being solved.
*     rmax:
*          This is the final value of the residual calcu;ated.
*
************************************************************************
*
*

      integer idim,ilower,iupper,jdim,jlower,jupper,ifmg
      integer id5,id9,idi,idg,ier
      CCTK_REAL cc(idim,jdim),cn(idim,jdim),cs(idim,jdim),cw(idim,jdim),
     &     ce(idim,jdim),u(idim,jdim),rhs(idim,jdim)
      CCTK_REAL eps,rmax

*
************************************************************************
*
*  Variable definitions:
*
*  Integers:
*     ifd59 :
*          ifd59 = 5 - means a 5-point finite difference stencil
*                      (ac,an,as,aw,ae) is defined on the finest grid.
*          ifd59 = 9 - means a 9-point finite difference stencil
*                      (ac,an,as,aw,ae,anw,ane,asw,ase) is defined on
*                      the finest grid by the user.
*                      (NOTE: This routine specifically written
*                             for 5-point)
*
*     ncyc :
*          The maximum number of multigrid "v"-cycles to be used. If
*          the maximum norm of the residual is not less than tol at the
*          end of ncyc cycles, the algorithm is terminated.
*          (NOTE: ncyc <= 40 )
*          
*     np(20,8) :
*          Input: When the iskip=1,-1 or -2 option is used, np2 is
*                 assumed to contain the grid information for umgs2.
*          Output: When the iskip=0 option is used, the grid
*                 information for umgs2 is returned in np2.
*          (NOTE: This is only useful for multiple instance problems)
*
*     iskip :
*          iskip = 0 - The coarse grid information, coarse grid
*                      operators and interpolation coefficients are
*                      calculated by umgs2.  This information is stored
*                      in the arrays ac, aw, as, asw, ase, pu, pd, np2
*                      and the variable.
*          iskip = 1 - The calculation of the coarse grid information,
*                      coarse grid operators and interpolation
*                      coefficients is skipped.  This option would be
*                      used when umgs2 has been called with iskip=0 and
*                      is being called again to solve a system of
*                      equations with the same matrix. This would be
*                      the case in, say, parabolic problems with time
*                      independent coefficients.
*          iskip =-1 - The set up of pointers (ugrdfn) is skipped.
*                      Coarse grid operators and interpolation
*                      coefficients are calculated and the given matrix
*                      equation is solved.  This option would be used
*                      when umgs2 has been called with iskip=0 and is
*                      being called again to solve a system of equations
*                      with a different matrix of the same dimensions.
*                      This would be the case for, say, parabolic
*                      problems with time dependent coefficients.
*          iskip =-2 - The set up of pointers (ugrdfn) is skipped.
*                      Coarse grid operators and interpolation
*                      coefficients are calculated and returned.
*                      No matrix solve.
*
*     ipc :
*          ipc = 0 or 1.
*          ipc is a multigrid parameter which determines the type of
*          interpolation to be used.  Usually ipc=1 is best.  However,
*          if the boundary condition equations have been absorbed into
*          the interior equations then ipc=0 can be used which results
*          in a slightly more efficient algorithm.
*
*     nman :
*          nman = 0 usually.
*          nman =1 signals that the fine grid equations are singular for
*          the case when homogeneous Neumann boundary conditions are
*          applied along the entire boundary.  In this case, the
*          difference equations are singular and the condition that the
*          integral of q over the domain be zero is added to the set of
*          difference equations.  This condition is satisfied by adding
*          the appropriate constant vector to q on the fine grid.  It is
*          assumed, in this case, that a well-defined problem has been
*          given to mgss2, i.e. the integral of f over the domain is
*          zero.
*
*     im :
*          The number of grid points in the x-direction (including two
*          ficticious points)
*     jm :
*          The number of grid points in the y-direction (including two
*          ficticious points)
*
*     linp :
*          This is a dummy argument left over from the authors
*          development of the code
*             Use:  common /io/ linp,lout
*
*     lout :
*          lout = unit number of output file into which the maximum norm
*          of the residual after each multigrid v-cycle" is printed.
*             Use:  common /io/ linp,lout
*
*     iscale :
*          Flag to indicate whether problem can be diagonally scaled to
*          speed convergence of the multigrid solver.
*
*     REALS:
*
*     ac(id5),an(id5),as(id5),aw(id5),ae(id5),
*     anw(id9),ane(id9),asw(id9),ase(id9) :
*          Input: ac, an, as, aw, ae, anw, ane, asw and ase contain the
*                 stencil coefficients for the difference operator on
*                 the finest grid. When the iskip=1 option is used,
*                 these arrays also are assumed to contain the coarse
*                 grid difference stencil coeficients.
*          Output: when the iskip=0 option is used, the coarse grid
*                 stencil coeficients are returned in ac, an, as, aw,
*                 ae, anw, ane, asw and ase.
*
*     ru(idi),rd(idi),rc(idi) :
*          Real work arrays.
*
*     pu(idi),pd(idi),pc(idi) :
*          Real work arrays.
*          Input: when the iskip=1 option is used, these arrays are
*                 assumed to contain the interpolation coefficients used
*                 in the semi-coarsening multigrid algorithm.
*          Output: when the iskip=0 option is used, the interpolation
*                 coeficients are returned in pu and pd.
*
*     f(id5) :
*          f contains the right hand side vector of the matrix
*          equation to be solved by umgs2.
*
*     q(id5) :
*          If ifmg=0, q contains the initial guess on the fine grid.
*          If ifmg=1, the initial guess on the fine grid is determined
*                     by the full multigrid process and the value of
*                     q on input to umgs2 not used.
*
*     tol :
*          tol > 0.  The maximum norm of the residual is calculated at
*                    the end of each multigrid cycle. The algorithm is
*                    terminated when this maximum becomes less than tol
*                    or when the maximum number of iterations (see ncyc)
*                    is exceeded.  It is up to the user to provide a
*                    meaningfull tolerance criteria for the particular
*                    problem being solved.
*          tol = 0.  Perform ncyc multigrid cycles.  Calculate and print
*                    the maximum norm of the residual after each cycle.
*          tol =-1.  Perform ncyc multigrid cycles.  The maximum norm of
*                    the final residual is calculated and returned in
*                    the variable rmax in the calling list of umgs2.
*          tol =-2.  Perform ncyc multigrid cycles.  The maximum norm of
*                    the residual is not calculated.
*
************************************************************************
*
*

      integer ncyc,ifd59
      parameter (ncyc=40,ifd59=5)

*     integer id5,id9,idi,idg
*  This is for a 103x28 grid
*     parameter(id5=6695,id9=3811,idi=3811,idg=103)
*  This is for a 203x56 grid
*     parameter(id5=24969,id9=13601,idi=13601,idg=203)
*  This is for a 403x118 grid
*     parameter(id5=100750,id9=52793,idi=52793,idg=403)

      integer np2(20,8)
      integer iskip,ipc,nman
      parameter (iskip=0,ipc=1,nman=0)
      integer irc,irurd,im,jm
      integer linp,lout
      common /io/ linp,lout

      CCTK_REAL ac(id5),an(id5),as(id5),aw(id5),ae(id5),
     &     anw(id9),ane(id9),asw(id9),ase(id9),
     &     q(id5),f(id5),gam(idg)

      CCTK_REAL ru(idi),rd(idi),rc(idi),pu(idi),pd(idi),pc(idi)

      CCTK_REAL tol

      integer iscale,itry

c clear all arrays
      ac = 0.
      an = 0.
      as = 0.
      aw = 0.
      ae = 0.
      q = 0.
      f = 0.
      gam = 0.

      ier=0
*
* Set some parameters for multigrid solver
*
      lout=6
c      rewind(unit=lout)
      irc=0
      irurd=0
      im=iupper-ilower+3
      jm=jupper-jlower+3
      tol=eps


*
*  Set up coefficients into vectors with correct indexing
*
      do 110 j=jlower,jupper
      do 100 i=ilower,iupper
      n=(j+1-jlower)*im + i+1-(ilower-1)
      ac(n)=cc(i,j)
      an(n)=cn(i,j)
      as(n)=cs(i,j)
      aw(n)=cw(i,j)
      ae(n)=ce(i,j)
      anw(n)=0.
      ane(n)=0.
      asw(n)=0.
      ase(n)=0.
      q(n)=u(i,j)
      f(n)=rhs(i,j)
  100 continue
  110 continue


*
*  Determine whether we can diagonal scale the problem to speed
*  convergence. Can only be done if there are no zeros on the main
*  diagonal (ie. central difference coefficient).
*
      iscale=1
      do 200 j=jlower,jupper
      do 205 i=ilower,iupper
      n=(j+1-jlower)*im + i+1-(ilower-1)
      if (ac(n) .eq. 0.) then
        iscale=0
      endif
  205 continue
  200 continue

*
*  Do the diagonal scaling if we can.
*
      if (iscale.eq.1) then

      do 210 j=jlower,jupper
      do 215 i=ilower,iupper
        n=(j+1-jlower)*im + i+1-(ilower-1)
        f(n)=f(n)/ac(n)
        ae(n)=ae(n)/ac(n)
        aw(n)=aw(n)/ac(n)
        as(n)=as(n)/ac(n)
        an(n)=an(n)/ac(n)
        ac(n)=1.
  215 continue
  210 continue

      endif

c 
c   Now call the multigrid routine

      itry=0
      
 1122 call umgs2(
     + ac,aw,as,ae,an,asw,ase,ane,anw,q,f,pu,pd,pc,ru,rd,rc,gam,np2,
     + ifd59,ifmg,ncyc,tol,nman,im,jm,id5,id9,idi,m,iskip,rmax,
     + ipc,irc,irurd)


      if ((rmax.gt.tol).and.(itry.le.5)) then
        itry=itry+1
        print*,"Retry #",itry
        goto 1122
      endif

      if (rmax.gt.tol) then
        ier = -1
        print*,"Did not converge"
        print*," maximum residual    = ",rmax
        print*," tolerance           = ",tol
      endif

*
*  Convert the solution back to the 2D array form
*
      do 510 j=jlower,jupper
      do 500 i=ilower,iupper
      n=(j+1-jlower)*im + i+1-(ilower-1)
      u(i,j)=q(n)
  500 continue
  510 continue

      return
      end
