#include "cctk.h"
c----------------------------------------------------------------------
      subroutine umgs2(
     + ac,aw,as,ae,an,asw,ase,ane,anw,q,f,pu,pd,pc,ru,rd,rc,gam,np2,
     + ifd59,ifmg,ncyc,tol,nman,im,jm,id5,id9,idi,m,iskip,rmax,
     + ipc,irc,irurd)
      implicit CCTK_REAL(a-h,o-z)
c----------------------------------------------------------------------
cdir$ noinline
c** SUBROUTINE UMGS2
c
c** COPYRIGHT: Ecodynamics Research Associates, Inc.
c
c** Date written: June, 1990
c** Author:  Steve Schaffer
c            Mathematics Department
c            New Mexico Tech
c            Socorro, NM  87801
c            505-835-5811
c
c** DESCRIPTION:
c     umgs2 is a black box symmetric matrix solver.  It is written
c     in unsymmetric storage mode and can be used to solve mildly
c     nonsymmetric problems.  The user provides a matrix and right hand
c     side vector corresponding to a 5 or 9 point finite difference/
c     finite volume discretization of a symmetric second order PDE.
c     umgs2 will construct a sequence of coarse grids and coarse
c     grid operators and then solve the matrix equation using a
c     y-direction semi-coarsening multigrid algorithm and return
c     the solution vector.  If a sequence of matrix problems are
c     to be solved using the same matrix, computational time can
c     be saved by skipping the construction of the coarse grid
c     information in subsequent calls to umgs2. 
c
c     The matrix on the finest grid is stored in the arrays ac,aw,as,
c     an,asw,ase,ane and anw.  The difference stencil at the point
c     (i,j) given by
c
c         nw  n  ne        anw(i,j)   an(i,j)   ane(i,j)
c          w  c  e    =    aw(i,j)    ac(i,j)   aw(i,j)
c         sw  s  se        asw(i,j)   as(i,j)   ase(i,j)
c
c     If the difference stencil on the fine grid is a 5 point stencil
c     then the arrays asw,ase,ane,anw are not used and the
c     stencil is given by
c
c             n                     an(i,j)
c          w  c  e    =    aw(i,j)  ac(i,j)   ae(i,j)
c             s                     as(i,j)
c
c     However, asw,ase,ane,anw still need to be dimensioned (by id9)
c     in the calling program as they are used in the coarse grid
c     calculations.
c
c** STORAGE
c     It is assumed that a set of ficticious points have been defined
c     along the entire boundary.  These points have nothing to do with
c     the solution and are used for programming convenience and
c     vectorization purposes.  Storage is allocated for the stencil
c     elements ac,aw,as,asw,ase, the solution vector, q, and the
c     right hand side vector, f, at these ficticious points.  The
c     stencils at these ficticious points and all stencil connections
c     to them are set to zero in the subroutine useta which is called
c     by umgs2.  The computational grid is depicted by
c
c         x   x   x   x              x   x   x   x
c
c         x   *   *   *              *   *   *   x
c
c         x   *   *   *              *   *   *   x
c    .
c    .
c    .
c         x   *   *   *              *   *   *   x
c
c         x   *   *   *              *   *   *   x
c
c         x   x   x   x              x   x   x   x
c
c     where x depicts the ficticious points and * depicts the interior
c     points.  The total storage requirements for the fine grid problem
c     is then 5*im*jm for 5 point stencils and 7*im*jm for 9 point
c     stencils.  The total storage requirements for the multigrid
c     solution is approximately 2 to 3 times that of the storage
c     requirements of the fine grid problem. (See DIMENSION PARAMETERS).
c     Note:  The first im*jm elements of the arrays ac,aw,as,[asw,ase],
c     q and f correspond to the finest grid.
c
c** DIMENSION PARAMETERS
c     The arrays ac,aw,ae,asw,ase,q,f,pu,pd and pc are dimensioned as one
c     dimensional arrays in the calling program.  They are dimensioned
c     as two dimensional arrays in the working subroutines.  The one
c     dimensional storage of the arrays, say q, follows: n=(j-1)*jm+i,
c     where n is the element location in the one dimensional storage of
c     q corresponding to the (i,j)th element of the two dimensional
c     storage of q and jm is the number of grid points in the j
c     direction (including the two ficticious points).
c
c     The dimension parameters are id5, id9, idi and idg.  They can be
c     determined by running the companion program MSS2DIM.F.
c      id5 - Integer variable.
c            Dimension of the arrays ac,aw,as,ae,an,q and f in the
c            calling program.  id5 is the total number of grid points
c            on the finest grid and all coarser grids.
c      id9 - Integer variable.
c            Dimension of the arrays asw,ase,ane,anw in the calling
c            program.  If ifd59=5 then id9=idi.  If ifd59=9 then
c            id9=id5.
c      idi - Integer variable.
c            Dimension of the work arrays pu and pd in the calling
c            program.  idi is the total number of grid points on all
c            of the coarser grids.
c      idg - Integer variable.
c            Dimension of the work array gam in the calling program.
c            It is set to the value im, the number of grid points
c            in the i-direction on the finest grid.
c
c** INPUT
c            (Note: all variable types are set implicitly)
c ac,aw,as
c    ae,an - Real arrays.  Dimensioned (id5) in calling program.
c            See comments in DESCRIPTION and DIMENSION PARAMETERS.
c  asw,ase
c  ane,anw - Real arrays.  Dimensioned (id9) in calling program.
c            See comments in DESCRIPTION and DIMENSION PARAMETERS.
c        f - Real array.  Dimensioned (id5) in calling program.
c            f contains the right hand side vector of the matrix
c            equation to be solved by umgs2.
c        q - Real array.  Dimensioned (id5) in calling program.
c            If ifmg=0, q contains the initial guess on the fine
c            grid.  If ifmg=1, the initial guess on the fine grid
c            is determined by the full multigrid process and the
c            value of q on input to umgs2 not used.
c   ifd59 -  Integer variable.
c            =5 - means a 5-point finite difference stencil (ac,aw and
c                 as) is defined on the finest grid by the user.
c            =9 - means a 9-point finite difference stencil (ac,aw,as,
c                 asw, ase) is defined on the finest grid by the user.
c    ifmg - Integer variable.
c           =0 - The full multigrid algorithm is not used to obtain a
c                good initial guess on the fine grid.
c           =1 - The full multigrid algorithm is used to obtain a good
c                initial guess on the fine grid.
c    ncyc - Integer variable.
c           The maximum number of multigrid v-cycles to be used.
c           If the maximum norm of the residual is not less than tol
c           at the end of ncyc cycles, the algorithm is terminated.
c     tol - Real variable.
c           >0 - The maximum norm of the residual is calculated at the
c           end of each multigrid cycle.  The algorithm is terminated
c           when this maximum becomes less than tol or when the maximum
c           number of iterations (see ncyc) is exceeded.  It is up to
c           the user to provide a meaningfull tolerance criteria for
c           the particular problem being solved.
c           =0 - Perform ncyc multigrid cycles.  Calculate and print
c           the maximum norm of the residual after each cycle.
c           =-1. - Perform ncyc multigrid cycles.  The maximum norm of
c           the final residual is calculated and returned in the
c           variable rmax in the calling list of umgs2.
c           =-2. - Perform ncyc multigrid cycles.  The maximum norm of
c           the residual is not calculated.
c   iskip - Integer variable.
c           =0 - The coarse grid information, coarse grid operators
c                and interpolation coefficients are calculated by
c                umgs2.  This information is stored in the arrays
c                ac, aw, as, asw, ase, pu, pd, np2 and the variable m
c                and returned to the calling program.
c           =1 - The calculation of the coarse grid information, coarse
c                grid operators and interpolation coefficients is
c                skipped.  This option would be used when umgs2 has
c                been called with iskip=0 and is being called again
c                to solve a system of equations with the same matrix.
c                This would be the case in, say, parabolic problems
c                with time independent coefficients.
c           =-1 -The set up of pointers (ugrdfn) is skipped.  Coarse grid
c                operators and interpolation coefficients are calculated
c                and the given matrix equation is solved.  This option
c                would be used when umgs2 has been called with iskip=0
c                and is being called again to solve a system of
c                equations with a different matrix of the same
c                dimensions.  This would be the case for, say,
c                parabolic problems with time dependent coefficients.
c           =-2 -The set up of pointers (ugrdfn) is skipped.  Coarse grid
c                operators and interpolation coefficients are calculated
c                and returned to the calling program.  No matrix solve.
c     ipc - Integer variable.
c           =0 or 1.
c           ipc is a multigrid parameter which determines the type of
c           interpolation to be used.  Usually ipc=1 is best.  However, if
c           the boundary contition equations have been absorbed into the
c           interior equations then ipc=0 can be used which results in a
c           slightly more efficient algorithm.
c    nman - Integer variable.
c           =0 usually.
c           =1 signals that the fine grid equations are singular for
c           the case when homogeneous Neumann boundary conditions are
c           applied along the entire boundary.  In this case, the
c           difference equations are singular and the condition that
c           the integral of q over the domain be zero is added to the
c           set of difference equations.  This condition is satisfied
c           by adding the appropriate constant vector to q on the fine
c           grid.  It is assumed, in this case, that a well-defined
c           problem has been given to mgss2, i.e. the integral of f
c           over the domain is zero.
c      im - Integer variable.
c           The number of grid points in the x-direction (including two
c           ficticious points)
c      jm - Integer variable.
c           The number of grid points in the y-direction (including two
c           ficticious points)
c    lout - Integer variable.
c           = unit number of output file into which the maximum norm
c             of the residual after each multigrid v-cycle is printed.
c             Use:  common /iout/ lout
c
c** INPUT/OUTPUT
c       q - Real array.  Dimensioned (id5)
c           On input, if ifmg=0, q contains the initial guess on the
c           finest grid for umgs2.  On output, q contains the final
c           solution on the finest grid.
c  ac-anw - Real arrays.  See DIMENSION.
c           On input, ac, aw, as, [asw and ase] contain the stencil
c           coefficients for the difference operator on the finest
c           grid.  When the iskip=1 option is used, these arrays
c           also are assumed to contain the coarse grid difference
c           stencil coeficients.
c           On output, when the iskip=0 option is used, the coarse
c           grid stencil coeficients are returned in ac - ase.
c
c ru,rd,rc - Real work arrays.  Dimensioned (idi)
c
c pu,pd,pc - Real work arrays.  Dimensioned (idi).
c           On input, when the iskip=1 option is used, these arrays
c           are assumed to contain the interpolation coefficients
c           used in the semi-coarsening multigrid algorithm.
c           On output, when the iskip=0 option is used, the
c           interpolation coeficients are returned in pu and pd.
c     np2 - Integer work array.  Dimensioned np2(20,8).
c           On input, when the iskip=1,-1 or -2 option is used, np2 is
c           assumed to contain the grid information for umgs2.
c           On output, when the iskip=0 option is used, the grid
c           information for umgs2 is returned in np2.
c** OUTPUT
c    rmax - If tol.ge.-1., the final residual norm is returned in rmax.
c
c** SUBROUTINES CALLED BY UMGS2
c
c    - ugrdfn, ukey, uintad, urelax, urscal, ursrhs, useta
c
c** END OF DESCRIPTION OF UMGS2
c .....................................................................
      CCTK_REAL ac(id5),aw(id5),as(id5),ae(id5),an(id5),asw(id9),
     + ase(id9),ane(id9),anw(id9),q(id5),f(id5)
      CCTK_REAL pu(idi),pd(idi),pc(idi),gam(im)
      integer np2(20,8)
      CCTK_REAL ru(idi),rd(idi),rc(idi)
      CCTK_REAL resid(0:40),confac(0:40)
      common /io/ linp,lout
      CCTK_REAL one
      parameter (one=1)
c clear arrays
      resid = 0.
      confac = 0.
      np2 = -1
c
c-time           tsu0=second()
         if(iskip.eq.0) then
      call ugrdfn(m,ifd59,is5,is9,isi,np2,im,jm)
         iquit=0
           if(m.gt.20) then
      iquit=1
      write(lout,*) ' m=',m,' > 20 - np2 is dimensioned np2(m=20,8)'
           endif
           if(is5.gt.id5) then
      iquit=1
      write(lout,*) ' id5=',id5,' too small.  Should be set to',is5
           endif
           if(is9.gt.id9) then
      iquit=1
      write(lout,*) ' id9=',id9,' too small.  Should be set to',is9
           endif
           if(isi.gt.idi) then
      iquit=1
      write(lout,*) ' idi=',idi,' too small.  Should be set to',isi
            endif
            if(is5.lt.2*im*jm) then
      iquit=1
      write(lout,*) ' id5.lt.2*im*jm can cause problems in useta'
      write(lout,*) ' this can be remedied by setting id5 larger'
            endif
      if(iquit.eq.1) return
         endif
         if(iskip.le.0) then
c     ----------  interpolation and coarse grid operators -----------
      do 5 k=m-1,1,-1
cdir$ inline 
      call ukey(k+1,np2,n5,n9,ni,jm,i9,j9,ifd,jr)
      call ukey(k,np2,n5c,n9c,nic,jmc,i9c,j9c,ifdc,jrc)
cdir$ noinline 
      if(k.eq.m-1) n5cqf=n5c
    5 call useta(
     + ac(n5),aw(n5),as(n5),ae(n5),an(n5),asw(n9),ase(n9),
     + ane(n9),anw(n9),ac(n5c),aw(n5c),as(n5c),ae(n5c),an(n5c),asw(n9c),
     + ase(n9c),ane(n9c),anw(n9c),pu(nic),pd(nic),pc(nic),ru(nic),
     + rd(nic),rc(nic),q(n5cqf),f(n5cqf),gam,
     + im,jm,jmc,ifd,i9,j9,nman,k+1,m,jr,ipc,irc,irurd)
         endif
      if(iskip.eq.-2) return
c
         if(ifmg.ge.1) then
      do 6 k=m-1,1,-1
cdir$ inline 
      call ukey(k+1,np2,n5,n9,ni,jm,i9,j9,ifd,jr)
      call ukey(k,np2,n5c,n9c,nic,jmc,i9c,j9c,ifdc,jrc)
cdir$ noinline 
c     TODO 2009-01-06 Erik Schnetter <schnetter@cct.lsu.edu>:
c     I notice that this call does not correspond to the subroutine
c     definition below.  The call passes three arguments too many.
c     Maybe the three arguments pu(nic),pd(nic),pc(nic) should be
c     omitted?
c      call CCTK_WARN (CCTK_WARN_ABORT, "Trying to find out whether this line is reached")
c    6 call ursrhs(f(n5),f(n5c),pu(nic),pd(nic),pc(nic),ru(nic),
c     + rd(nic),rc(nic),im,jm,jmc,m,k+1,jr,irc)
    6 call ursrhs(f(n5),f(n5c),ru(nic),rd(nic),rc(nic),
     + im,jm,jmc,m,k+1,jr,irc)
         endif
c-time            tsu1=second()
c-time            write(lout,*) ' time for setup =',tsu1-tsu0
      l=1
      if(ifmg.eq.0) l=m
      k=l
      mcyc=0
      rmaxo=1.
c     ----------   begin multigrid cycling  ----------------------------
c
      if(l.eq.1) go to 20
cdir$ inline 
   10 call ukey(k,np2,n5,n9,ni,jm,i9,j9,ifd,jr)
cdir$ noinline 
      call urelax(
     + ac(n5),aw(n5),as(n5),ae(n5),an(n5),asw(n9),ase(n9),
     + ane(n9),anw(n9),f(n5),q(n5),gam,
     + im,jm,i9,j9,ifd,nman,k,m,jr,0,0)
cdir$ inline 
      call ukey(k-1,np2,n5c,n9c,nic,jmc,i9c,j9c,ifdc,jrc)
cdir$ noinline 
      call urscal(
     + ac(n5),aw(n5),as(n5),ae(n5),an(n5),asw(n9),ase(n9),
     + ane(n9),anw(n9),q(n5),f(n5),f(n5c),q(n5c),rc(nic),
     + im,jm,jmc,ifd,i9,j9,k,m,jr,tol,rmax,ipc,irc)
      if(k.eq.m.and.rmax.lt.tol) go to 60
         if(k.eq.m.and.tol.ge.-.5) then
      if(rmaxo.ne.0.) then
        rate=rmax/rmaxo
      else
        rate=1.
      endif
      rmaxo=rmax
      if(mcyc.eq.0) rmax0=rmax
      resid(mcyc)=rmax
      confac(mcyc)=rate
         endif
      if(tol.eq.-.5) write(lout,*) ' down ',k,rmax
      k=k-1
      if(k.gt.1) go to 10
c     ---------  solve coarsest grid  ----------------------------------
c
cdir$ inline 
   20 call ukey(1,np2,n5,n9,ni,jm,i9,j9,ifd,jr)
cdir$ noinline 
      call urelax(
     + ac(n5),aw(n5),as(n5),ae(n5),an(n5),asw(n9),ase(n9),
     + ane(n9),anw(n9),f(n5),q(n5),gam,
     + im,jm,i9,j9,ifd,nman,k,m,jr,0,0)
         if(l.eq.1) go to 40
c     ----------  interpolate correction to next finer grid  -----------
c
   30 k=k+1
cdir$ inline 
      call ukey(k,np2,n5,n9,ni,jm,i9,j9,ifd,jr)
      call ukey(k-1,np2,n5c,n9c,nic,jmc,i9c,j9c,ifdc,jrc)
cdir$ noinline 
      call uintad(
     + q(n5),q(n5c),pu(nic),pd(nic),im,jm,jmc,1,jr,ipc)
      call urelax(
     + ac(n5),aw(n5),as(n5),ae(n5),an(n5),asw(n9),ase(n9),
     + ane(n9),anw(n9),f(n5),q(n5),gam,
     + im,jm,i9,j9,ifd,nman,k,m,jr,0,0)
         if(tol.eq.-.5) then
      call urscal(
     + ac(n5),aw(n5),as(n5),ae(n5),an(n5),asw(n9),ase(n9),
     + ane(n9),anw(n9),q(n5),f(n5),f(n5c),q(n5c),rc(nic),
     + im,jm,jmc,ifd,i9,j9,k,m,jr,tol,rmax,ipc,irc)
      write(lout,*) ' up   ',k,rmax
         endif
      if(k.lt.l) go to 30
      if(l.eq.m) go to 50
c     ----------  interpolate solution to new finest grid l+1 in fmg  ----
c
   40 l=l+1
      k=l
cdir$ inline 
      call ukey(l,np2,n5,n9,ni,jm,i9,j9,ifd,jr)
      call ukey(l-1,np2,n5c,n9c,nic,jmc,i9c,j9c,ifdc,jrc)
cdir$ noinline 
      call uintad(
     + q(n5),q(n5c),pu(nic),pd(nic),im,jm,jmc,0,jr,0)
      go to 10
c
   50 if(nman.eq.1) call uneuman(q(n5),im,jm)
      mcyc=mcyc+1
c     ----------  Cycle ncyc times on grid m  ----------------------------
      if(mcyc.lt.ncyc) go to 10
c-time            tmg1=second()
c-time            write(lout,*) ' time in ',ncyc,' cycles =',tmg1-tsu1
c
c     ----------  print out final residual and work units  ---------------
        if(tol.ge.-1.) then
cdir$ inline 
      call ukey(m,np2,n5,n9,ni,jm,i9,j9,ifd,jr)
      call ukey(m-1,np2,n5c,n9c,nic,jmc,i9c,j9c,ifdc,jrc)
cdir$ noinline 
      call urscal(
     + ac(n5),aw(n5),as(n5),ae(n5),an(n5),asw(n9),ase(n9),
     + ane(n9),anw(n9),q(n5),f(n5),f(n5c),q(n5c),rc(nic),
     + im,jm,jmc,ifd,i9,j9,k,m,jr,one,rmax,ipc,irc)
      resid(mcyc)=rmax
      confac(mcyc)=rmax/rmaxo
      nb=0
      ne=min(6,mcyc)
 2029 write(lout,2033) (mc,mc=nb,ne)
      write(lout,2032) (resid(mc),mc=nb,ne)
      write(lout,2031) (confac(mc),mc=nb,ne)
      nb=ne+1
      ne=ne+min(6,mcyc-ne)
      if(nb.le.ne) go to 2029
      fconfac=(rmax/rmax0)**(1./dble(mcyc))
      write(lout,2034) fconfac
 2034 format(30x,6(1h*)/,' average convergence factor =',f7.3,/,
     + 30x,6(1h*))
 2033 format(7(4x,i2,4x))
 2031 format(7(1x,f9.3))
 2032 format(7(1x,e9.3))
         endif
      return
   60 write(lout,1003) mcyc,resid(mcyc),tol
      return
 1003 format(' cyc=',i2,' max(res)=',1pe8.2/
     + ' tolerance condition tol=',1pe8.2,' satisfied')
      end
c----------------------------------------------------------------------
      subroutine ugrdfn(m,ifd59,is5,is9,isi,np2,imx,jmx)
      implicit CCTK_REAL(a-h,o-z)
c----------------------------------------------------------------------
cdir$ noinline
c  Given imx, jmx and ifd59 (See comments in mgss2), ugrdfn calculates
c  the number of grids that will be needed.  Pointers into the arrays
c  ac, aw, as, asw, ase, q, f, pu, pd, pc, ru, rd and rc and the size
c  of each grid is calculated and stored in the array np2.  The 
c  subroutine ukey is called to retrieve the grid information.
c .....................................................................
      parameter(n5=1,n9=2,ni=3,jm=4,i9=5,j9=6,ifd=7,jred=8)
      integer np2(20,8)
      common /cs/ icorstr,iprint
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
      is5=iq5-1
      is9=iq9-1
      isi=iqi-1
      return
      end
c----------------------------------------------------------------------
      subroutine ukey(k,np2,nn5,nn9,nni,jjm,ii9,jj9,iifd,jjred)
      implicit CCTK_REAL(a-h,o-z)
c----------------------------------------------------------------------
c  Returns the grid pointers and dimension variables for grid k.  The
c  information is stored in the array np2.
c......................................................................
      parameter(n5=1,n9=2,ni=3,jm=4,i9=5,j9=6,ifd=7,jred=8)
      integer np2(20,8)
      nn5=np2(k,n5)
      nn9=np2(k,n9)
      nni=np2(k,ni)
      jjm=np2(k,jm)
      ii9=np2(k,i9)
      jj9=np2(k,j9)
      iifd=np2(k,ifd)
      jjred=np2(k,jred)
      return
      end
c----------------------------------------------------------------------
      subroutine uintad(q,qc,pu,pd,im,jm,jmc,iadd,jred,ipc)
      implicit CCTK_REAL (a-h,o-z)
c----------------------------------------------------------------------
c  iadd=1:
c  Interpolates and adds the coarse grid (kf-1) correction, qc, to the
c  fine grid (kf) approximation, q, at the black y-lines.
c  iadd=0:
c  In the full multigrid algorithm, the solution to the coarse grid
c  (kf-1) difference equation is interpolated to the fine grid (kf)
c  to be used as the initial guess vector for kf=2,3,...,m.
c  Interpolation is at black y-lines only.
c .....................................................................
      CCTK_REAL q(im,jm),qc(im,jmc),pu(im,jmc),pd(im,jmc)
      im1=im-1
      jm1=jm-1
      jblack=5-jred
c                                    add correction to next finer grid
         if(iadd.eq.1) then
      jc=3-jred
      do 10 j=jblack,jm1,2
      jc=jc+1
      do 10 i=2,im1
   10 q(i,j)=q(i,j)+pd(i,jc)*qc(i,jc)+pu(i,jc)*qc(i,jc+1)
c
c                       interpolate solution to next finer grid in fmg
         else
      jc=3-jred
      do 40 j=jblack,jm1,2
      jc=jc+1
      do 40 i=2,im1
   40 q(i,j)=pd(i,jc)*qc(i,jc)+pu(i,jc)*qc(i,jc+1)
         endif
      return
      end
c----------------------------------------------------------------------
      subroutine uneuman(q,im,jm)
      implicit CCTK_REAL (a-h,o-z)
c----------------------------------------------------------------------
c  For problems with homogeneous Neumann boundary contitions, the
c  condition that the integral of q over the domain be zero is added
c  to the set of difference equations in order to obtain a unique
c  solution.
c......................................................................
      CCTK_REAL q(im,jm)
      im1=im-1
      jm1=jm-1
      con=0.
      do 10 j=2,jm1
      do 10 i=2,im1
   10 con=con+q(i,j)
      con=con/((im-2)*(jm-2))
      do 20 j=2,jm1
      do 20 i=2,im1
   20 q(i,j)=q(i,j)-con
      return
      end
c----------------------------------------------------------------------
      subroutine urelax(ac,aw,as,ae,an,asw,ase,ane,anw,f,q,gam,
     + im,jm,i9,j9,ifd,nman,k,m,jred,ipc,iprcud)
      implicit CCTK_REAL (a-h,o-z)
c----------------------------------------------------------------------
c  Performs red/black x-line relaxation.  The Thomas algorithm is used
c  to solve the tridiagonal matrices.
c** INPUT -
c       ac-anw=  finite difference operator coeficients
c            q=  initial approximation
c            f=  right hand side vector
c        im,jm=  the number of grid points in the x,y-directions
c        i9,j9=  the i,j-dimensions of the arrays asw,ase
c          ifd=  5 or 9 - the size of the stencil
c         nman-  =0 usually.
c                =1 signals that the fine grid equations are singular
c                for the case when Neumann boundary conditions are
c                applied along the entire boundary.  In this case, the
c                equations on the coarsest grid (consisting of a single
c                line of unknowns) is a singular tridiagonal system
c                and the Thomas algorithm is modified on this grid to
c                obtain a solution with an arbitrary constant vector
c                component.  This constant vector is removed on the
c                finest grid by the call to subroutine uneuman.
c** OUTPUT -
c            q=  final approximation after a red/black relaxation sweep
c .....................................................................
      CCTK_REAL ac(im,jm),aw(im,jm),as(im,jm),ae(im,jm),an(im,jm),
     + asw(i9,j9),ase(i9,j9),ane(i9,j9),anw(i9,j9)
      CCTK_REAL f(im,jm),q(im,jm),gam(im)
      jm1=jm-1
      im1=im-1
      im2=im-2
      jblack=5-jred
c                                              usual red/black relaxatio
      nrel=2
      jrb=jred
c                                             ipc ..brbr relaxation swee
         if(iprcud.eq.1) then
      nrel=ipc
      if(mod(ipc,2).eq.0) jrb=jblack
c                                     1 black relax for calc g pu,pd,ru,
         elseif(iprcud.eq.2) then
      nrel=1
      jrb=jblack
         endif
c
c
      do 109 nrr=1,nrel
         if(jrb.eq.jblack) then
c                                                             black rela
         if(jblack.le.jm1) then
         if(iprcud.ne.2) then
c
      do 110 j=jblack,jm1,2
      do 110 i=2,im1
  110 q(i,j)=f(i,j)-as(i,j)*q(i,j-1)-an(i,j)*q(i,j+1)
         if(ifd.eq.9) then
      do 120 j=jblack,jm1,2
      do 120 i=2,im1
  120 q(i,j)=q(i,j)-asw(i,j)*q(i-1,j-1)-ase(i,j)*q(i+1,j-1)-
     + anw(i,j)*q(i-1,j+1)-ane(i,j)*q(i+1,j+1)
         endif
         endif
c                                                black tridiagonal solve
c**
c**  Moved calculation of loop 129 from loop 130 for vectorization
c**  on vector machines (ie. Cray)
c**  By: John Towns 2/6/92
c**
      do 129 j=jblack,jm1,2        
  129 q(2,j)=q(2,j)/ac(2,j)

c**
c**  Changed bet=(quantity) to bet=1./(quantity) to trade two divisions
c**  for one division and two multiplies (more efficient on all
c**  machines)
c**  By: John Towns 2/6/92
c**

c**
c**  Cray compiler directives to parallelize tridiagonal solve.
c**  By: John Towns 4/13/92
c**

cmic$ parallel  private(bet,gam,i)
cmic$1shared(ac,ae,aw,q,jblack,jm1,im1,im2)
cmic$ do parallel
      do 130 j=jblack,jm1,2
      bet=1./ac(2,j)
      do 140 i=3,im1
      gam(i)=ae(i-1,j)*bet
      bet=1./(ac(i,j)-aw(i,j)*gam(i))
  140 q(i,j)=(q(i,j)-aw(i,j)*q(i-1,j))*bet
      do 150 i=im2,2,-1
  150 q(i,j)=q(i,j)-gam(i+1)*q(i+1,j)
  130 continue
cmic$ end do
cmic$ end parallel
         endif
c                                                             red relax
         else
c
      do 210 j=jred,jm1,2
      do 210 i=2,im1
  210 q(i,j)=f(i,j)-as(i,j)*q(i,j-1)-an(i,j)*q(i,j+1)
         if(ifd.eq.9) then
      do 220 j=jred,jm1,2
      do 220 i=2,im1
  220 q(i,j)=q(i,j)-asw(i,j)*q(i-1,j-1)-ase(i,j)*q(i+1,j-1)-
     + anw(i,j)*q(i-1,j+1)-ane(i,j)*q(i+1,j+1)
         endif
c                                                      tridiagonal solve
c                          nman=1 ==> avoid singularity on coarsest grid
      imm=im1
          if(nman.eq.1.and.k.eq.1) then
      imm=im-2
      q(im1,2)=0.
      gam(im1)=0.
          endif
c
c**
c**  Moved calculation of loop 229 from loop 230 for vectorization
c**  on vector machines (ie. Cray)
c**  By: John Towns 2/6/92
c**
      do 229 j=jred,jm1,2
  229 q(2,j)=q(2,j)/ac(2,j)

c**
c**  Changed bet=(quantity) to bet=1./(quantity) to trade two divisions
c**  for one division and two multiplies (more efficient on all
c**  machines)
c**  By: John Towns 2/6/92
c**

c**
c**  Cray compiler directives to parallelize tridiagonal solve.
c**  By: John Towns 4/13/92
c**

cmic$ parallel  private(bet,gam,i) 
cmic$1shared(ac,ae,aw,q,jred,jm1,im2,imm)
cmic$ do parallel
      do 230 j=jred,jm1,2
      bet=1./ac(2,j)
      do 240 i=3,imm
      gam(i)=ae(i-1,j)*bet
      bet=1./(ac(i,j)-aw(i,j)*gam(i))
  240 q(i,j)=(q(i,j)-aw(i,j)*q(i-1,j))*bet
      do 250 i=im2,2,-1
  250 q(i,j)=q(i,j)-gam(i+1)*q(i+1,j)
  230 continue
cmic$ end do
cmic$ end parallel
         endif
      jrb=5-jrb
  109 continue
      return
      end
c----------------------------------------------------------------------
      subroutine urscal(
     + ac,aw,as,ae,an,asw,ase,ane,anw,q,f,fc,qc,rc,
     + im,jm,jmc,ifd,i9,j9,kf,m,jred,tol,rmax,ipc,irc)
      implicit CCTK_REAL (a-h,o-z)
c----------------------------------------------------------------------
c  Defines the grid kf-1 right hand side, fc, as the restriction of the
c  grid kf residual.  The restriction operator is the transpose of the
c  interpolation operator.  Note:  The grid kf residual is zero at the
c  black lines (j-direction) as a result of red/black relaxation.
c  Thus, the restriction is simple injection.  The initial guess, qc,
c  for the coarse grid correction equation is set to zero.  The
c  maximum norm of the residual is calculated and returned in rmax.
c......................................................................
      CCTK_REAL ac(im,jm),aw(im,jm),as(im,jm),ae(im,jm),an(im,jm),
     + asw(i9,j9),ase(i9,j9),ane(i9,j9),anw(i9,j9)
      CCTK_REAL f(im,jm),q(im,jm),fc(im,jmc),qc(im,jmc)
      CCTK_REAL rc(im,jmc)
      rmax=0.
      im1=im-1
      jm1=jm-1
      jmc1=jmc-1
      jc=1
      do 10 j=jred,jm1,2
      jc=jc+1
      do 10 i=2,im1
   10 fc(i,jc)=f(i,j)-as(i,j)*q(i,j-1)-an(i,j)*q(i,j+1)-
     + aw(i,j)*q(i-1,j)-ae(i,j)*q(i+1,j)-ac(i,j)*q(i,j)
         if(ifd.eq.9) then
      jc=1
      do 20 j=jred,jm1,2
      jc=jc+1
      do 20 i=2,im1
   20 fc(i,jc)=fc(i,jc)-asw(i,j)*q(i-1,j-1)-ane(i,j)*q(i+1,j+1)-
     + ase(i,j)*q(i+1,j-1)-anw(i,j)*q(i-1,j+1)
         endif
c                                           zero out qc as initial guess
      do 25 jc=1,jmc
      do 25 i=1,im
   25 qc(i,jc)=0.
c                                        if kf=m calculate residual norm
         if((kf.eq.m.and.tol.ge.0.).or.tol.eq.-.5) then
      do 30 jc=2,jmc1
      do 30 i=2,im1
      resmax=abs(fc(i,jc))
   30 if(resmax.gt.rmax) rmax=resmax
         endif
c                                                 weight rhs if irc.ge.1
         if(irc.eq.1.and.ipc.ge.1) then
      do 40 jc=2,jmc1
      do 40 i=2,im1
   40 fc(i,jc)=rc(i,jc)*fc(i,jc)
         endif
c
      return
      end
c----------------------------------------------------------------------
      subroutine ursrhs(f,fc,ru,rd,rc,im,jm,jmc,m,kf,jred,irc)
      implicit CCTK_REAL (a-h,o-z)
c----------------------------------------------------------------------
c  Restricts the right hand side vector on grid kf onto grid kf-1 when
c  the full multigrid (ifmg>0) option is used.  The restriction operator
c  is NOT necessarily the transpose of the interpolation operator.
c......................................................................
      CCTK_REAL f(im,jm),fc(im,jmc),ru(im,jmc),rd(im,jmc),rc(im,jmc)
      jm1=jm-1
      im1=im-1
      jc=1
         if(irc.eq.0) then
      do 10 j=jred,jm1,2
      jc=jc+1
      do 10 i=2,im1
   10 fc(i,jc)=ru(i,jc-1)*f(i,j-1)+rd(i,jc)*f(i,j+1)+f(i,j)
         else
      do 20 j=jred,jm1,2
      jc=jc+1
      do 20 i=2,im1
   20 fc(i,jc)=ru(i,jc-1)*f(i,j-1)+rd(i,jc)*f(i,j+1)+
     + rc(i,jc)*f(i,j)
         endif
      return
      end
c----------------------------------------------------------------------
      subroutine useta(
     + ac,aw,as,ae,an,asw,ase,ane,anw,acc,awc,asc,aec,
     + anc,aswc,asec,anec,anwc,pu,pd,pc,ru,rd,rc,qw,fw,gam,
     + im,jm,jmc,ifd,i9,j9,nman,kf,m,jred,ipc,irc,irurd)
      implicit CCTK_REAL (a-h,o-z)
c----------------------------------------------------------------------
cdir$ noinline
c     Calculates the interpolation coefficients from grid kf-1 to
c     grid kf and the coarse grid operator on grid kf-1.
c** INPUT -
c    ac - anw = fine grid (kf) array stencil coeficients
c            m=  total number of grids
c           kf=  grid number of the fine grid
c          ifd=  the size of the fine grid stencil (= 5 or 9)
c        i9,j9=  the i,j-dimensions of the arrays asw,ase
c        qw,fw=  coarse grid portions of q and f used for work arrays here
c   (See comments in MGSS2 for details)
c** OUTPUT -
c  acc - anwc = coarse grid (kf-1) array stencil coeficients
c        pu,pd=  arrays of interpolation coefficients from grid kf-1
c                to grid kf
c .....................................................................
      CCTK_REAL ac(im,jm),aw(im,jm),as(im,jm),ae(im,jm),an(im,jm),
     + asw(i9,j9),ase(i9,j9),ane(i9,j9),anw(i9,j9),
     + ru(im,jmc),rd(im,jmc),rc(im,jmc),
     + pu(im,jmc),pd(im,jmc),pc(im,jmc),gam(im)
      CCTK_REAL acc(im,jmc),awc(im,jmc),asc(im,jmc),aec(im,jmc),
     + anc(im,jmc),aswc(im,jmc),asec(im,jmc),anec(im,jmc),anwc(im,jmc)
      CCTK_REAL qw(im,jm),fw(im,jm)
      common /io/ linp,lout
      common /prsol/ iprsol
c
      pcscale=.001
c
      im1=im-1
      jm1=jm-1
      jmc1=jmc-1
      jblack=5-jred
c                           zeroing out connections to fictitious points
      do 1 j=1,jm
      do 2 i=1,im,im1
      ac(i,j)=0.
      aw(i,j)=0.
      as(i,j)=0.
      ae(i,j)=0.
    2 an(i,j)=0.
      aw(2,j)=0.
    1 ae(im1,j)=0.
      do 3 i=1,im
      do 4 j=1,jm,jm1
      ac(i,j)=0.
      aw(i,j)=0.
      as(i,j)=0.
      ae(i,j)=0.
    4 an(i,j)=0.
      as(i,2)=0.
    3 an(i,jm1)=0.
         if(ifd.eq.9) then
      do 5 j=1,jm
      do 6 i=1,im,im1
      asw(i,j)=0.
      ase(i,j)=0.
      ane(i,j)=0.
    6 anw(i,j)=0.
      asw(2,j)=0.
      anw(2,j)=0.
      ase(im1,j)=0.
    5 ane(im1,j)=0.
      do 7 i=1,im
      do 8 j=1,jm,jm1
      asw(i,j)=0.
      ase(i,j)=0.
      ane(i,j)=0.
    8 anw(i,j)=0.
      ase(i,2)=0.
      asw(i,2)=0.
      ane(i,jm1)=0.
    7 anw(i,jm1)=0.
      endif
c
      do 9 jc=1,jmc
      do 9 i=1,im
      pc(i,jc)=0.
      pu(i,jc)=0.
      pd(i,jc)=0.
      rc(i,jc)=0.
      ru(i,jc)=0.
    9 rd(i,jc)=0.
c
c                               calculation of interpolation coeficients
c

c                                                              define pc
         if(ipc.ge.1) then
      do 20 j=1,jm
      do 20 i=1,im
      fw(i,j)=0.
   20 qw(i,j)=1.
c
      call urelax(ac,aw,as,ae,an,asw,ase,ane,anw,fw,qw,gam,im,jm,
     + i9,j9,ifd,nman,kf,m,jred,ipc,1)
c                                                              scale pc
      pcmax=0.
      jc=1
      do 40 j=jred,jm1,2
      jc=jc+1
      do 40 i=2,im1
      pc(i,jc)=qw(i,j)
   40 pcmax=max(pcmax,abs(qw(i,j)))
      do 50 jc=2,jmc1
      do 50 i=2,im1
      if(pc(i,jc).eq.0.) pc(i,jc)=pcscale
   50 pc(i,jc)=pc(i,jc)/pcmax
c
         else
      do 55 jc=2,jmc1
      do 55 i=2,im1
   55 pc(i,jc)=1.
         endif
c
c                                                              define pu
      jc=3-jred
      do 60 j=jblack,jm1,2
      jc=jc+1
         if(ipc.eq.0) then
      do 70 i=2,im1
   70 qw(i,j)=-an(i,j)
         if(ifd.eq.9) then
      do 80 i=2,im1
   80 qw(i,j)=qw(i,j)-ane(i,j)-anw(i,j)
         endif
         else
      do 90 i=2,im1
   90 qw(i,j)=-an(i,j)*pc(i,jc+1)
         if(ifd.eq.9) then
      do 100 i=2,im1
  100 qw(i,j)=qw(i,j)-ane(i,j)*pc(i+1,jc+1)-anw(i,j)*pc(i-1,jc+1)
         endif
         endif
   60 continue
c                                                          solve for pu
      call urelax(ac,aw,as,ae,an,asw,ase,ane,anw,fw,qw,gam,im,jm,
     + i9,j9,ifd,nman,kf,m,jred,ipc,2)
c

      jc=3-jred
      do 102 j=jblack,jm1,2
      jc=jc+1
         if(j.lt.jm1) then
      do 103 i=2,im1
  103 pu(i,jc)=qw(i,j)
         endif
  102 continue
c
c                                                              define pd
      jc=3-jred
      do 106 j=jblack,jm1,2
      jc=jc+1
         if(ipc.eq.0) then
      do 130 i=2,im1
  130 qw(i,j)=-as(i,j)
         if(ifd.eq.9) then
      do 140 i=2,im1
  140 qw(i,j)=qw(i,j)-ase(i,j)-asw(i,j)
         endif
c
         else
c
      do 150 i=2,im1
  150 qw(i,j)=-as(i,j)*pc(i,jc)
         if(ifd.eq.9) then
      do 160 i=2,im1
  160 qw(i,j)=qw(i,j)-ase(i,j)*pc(i+1,jc)-asw(i,j)*pc(i-1,jc)
         endif
         endif
  106 continue
c                                                          solve for pd
      call urelax(ac,aw,as,ae,an,asw,ase,ane,anw,fw,qw,gam,im,jm,
     + i9,j9,ifd,nman,kf,m,jred,ipc,2)
c
      jc=3-jred
      do 105 j=jblack,jm1,2
      jc=jc+1
         if(j.gt.2) then
      do 104 i=2,im1
  104 pd(i,jc)=qw(i,j)
         endif
  105 continue
c
c                                            define restriction operator
c
c                                                              define rc
         if(irc.eq.1) then
      do 500 jc=2,jmc1
      do 500 i=2,im1
  500 rc(i,jc)=pc(i,jc)
         else
      do 502 jc=2,jmc1
      do 502 i=2,im1
  502 rc(i,jc)=1.
         endif
c   
c                                           compute qw = -Cb(inv) * eb*
         if(irurd.ge.1) then
      jc=3-jred
         if(irurd.eq.1) then
      do 560 j=jblack,jm1,2
      jc=jc+1
      do 560 i=2,im1
  560 qw(i,j)=1.
         elseif(irurd.eq.2) then
      do 561 j=jblack,jm1,2
      jc=jc+1
      do 561 i=2,im1
  561 qw(i,j)=(pd(i,jc)*pc(i,jc)+pu(i,jc)*pc(i,jc+1))
         endif
c
      call urelax(ac,aw,as,ae,an,asw,ase,ane,anw,fw,qw,gam,im,jm,
     + i9,j9,ifd,nman,kf,m,jred,ipc,2)
c
      jc=3-jred
      do 566 j=jblack,jm1,2
      jc=jc+1
c                                           compute ru = -b(j+1) * qw
         if(j.lt.jm1) then
      do 570 i=2,im1
  570 ru(i,jc)=-as(i,j+1)*qw(i,j)
         if(ifd.eq.9) then
      do 580 i=2,im1
  580 ru(i,jc)=ru(i,jc)-ase(i,j+1)*qw(i+1,j)-asw(i,j+1)*qw(i-1,j)
         endif
         endif
c                             compute rd = -a(j-1) * c(j)(inv) * qw
         if(j.gt.2) then
      do 650 i=2,im1
  650 rd(i,jc)=-an(i,j-1)*qw(i,j)
         if(ifd.eq.9) then
      do 660 i=2,im1
  660 rd(i,jc)=rd(i,jc)-ane(i,j-1)*qw(i+1,j)-anw(i,j-1)*qw(i-1,j)
         endif
         endif
  566 continue
c
         else
c                                              else set ru=pu and rd=pd
      jc=3-jred
      do 670 j=jblack,jm1,2
      jc=jc+1
      do 670 i=2,im1
      ru(i,jc)=pu(i,jc)
  670 rd(i,jc)=pd(i,jc)
         endif
c
c                                   calculating the coarse grid operator
c
         if(ipc+irc+irurd.eq.0) then
      j=jred-2
      do 200 jc=2,jmc1
      j=j+2
      do 200 i=2,im1
      acc(i,jc)=ac(i,j)+an(i,j-1)*pu(i,jc-1)+as(i,j+1)*pd(i,jc)+
     + pu(i,jc-1)*(as(i,j)+ac(i,j-1)*pu(i,jc-1))+
     + pd(i,jc)*(an(i,j)+ac(i,j+1)*pd(i,jc))
      awc(i,jc)=aw(i,j)+pd(i-1,jc)*aw(i,j+1)*pd(i,jc)+pu(i-1,jc-1)*
     + aw(i,j-1)*pu(i,jc-1)
      asc(i,jc)=as(i,j)*pd(i,jc-1)+pu(i,jc-1)*(as(i,j-1)+
     + ac(i,j-1)*pd(i,jc-1))
      aec(i,jc)=ae(i,j)+pd(i+1,jc)*ae(i,j+1)*pd(i,jc)+pu(i+1,jc-1)*
     + ae(i,j-1)*pu(i,jc-1)
      anc(i,jc)=an(i,j)*pu(i,jc)+pd(i,jc)*(an(i,j+1)+
     + ac(i,j+1)*pu(i,jc))
      aswc(i,jc)=pd(i-1,jc-1)*aw(i,j-1)*pu(i,jc-1)
      asec(i,jc)=pd(i+1,jc-1)*ae(i,j-1)*pu(i,jc-1)
      anec(i,jc)=pu(i+1,jc)*ae(i,j+1)*pd(i,jc)
  200 anwc(i,jc)=pu(i-1,jc)*aw(i,j+1)*pd(i,jc)
         if(ifd.eq.9) then
      j=jred-2
      do 210 jc=2,jmc1
      j=j+2
      do 210 i=2,im1
      awc(i,jc)=awc(i,jc)+asw(i,j+1)*pd(i,jc)+anw(i,j-1)*pu(i,jc-1)+
     + pd(i-1,jc)*anw(i,j)+pu(i-1,jc-1)*asw(i,j)
      aec(i,jc)=aec(i,jc)+ase(i,j+1)*pd(i,jc)+ane(i,j-1)*pu(i,jc-1)+
     + pd(i+1,jc)*ane(i,j)+pu(i+1,jc-1)*ase(i,j)
      aswc(i,jc)=aswc(i,jc)+asw(i,j-1)*pu(i,jc-1)+pd(i-1,jc-1)*asw(i,j)
      asec(i,jc)=asec(i,jc)+ase(i,j-1)*pu(i,jc-1)+pd(i+1,jc-1)*ase(i,j)
      anec(i,jc)=anec(i,jc)+ane(i,j+1)*pd(i,jc)+pu(i+1,jc)*ane(i,j)
  210 anwc(i,jc)=anwc(i,jc)+anw(i,j+1)*pd(i,jc)+pu(i-1,jc)*anw(i,j)
         endif
c
         else
c
      j=jred-2
      do 300 jc=2,jmc1
      j=j+2
      do 300 i=2,im1
      acc(i,jc)=rc(i,jc)*(ac(i,j)*pc(i,jc)+
     +                    as(i,j)*pu(i,jc-1)+
     +                    an(i,j)*pd(i,jc))+
     +        ru(i,jc-1)*(ac(i,j-1)*pu(i,jc-1)+
     +                    an(i,j-1)*pc(i,jc))+
     +          rd(i,jc)*(ac(i,j+1)*pd(i,jc)+
     +                    as(i,j+1)*pc(i,jc))
      awc(i,jc)=rc(i,jc)*aw(i,j)*pc(i-1,jc)+
     +          rd(i,jc)*aw(i,j+1)*pd(i-1,jc)+
     +        ru(i,jc-1)*aw(i,j-1)*pu(i-1,jc-1)
      asc(i,jc)=rc(i,jc)*as(i,j)*pd(i,jc-1)+
     +       ru(i,jc-1)*(ac(i,j-1)*pd(i,jc-1)+
     +                   as(i,j-1)*pc(i,jc-1))
      aec(i,jc)=rc(i,jc)*ae(i,j)*pc(i+1,jc)+
     +          rd(i,jc)*ae(i,j+1)*pd(i+1,jc)+
     +        ru(i,jc-1)*ae(i,j-1)*pu(i+1,jc-1)
      anc(i,jc)=rc(i,jc)*an(i,j)*pu(i,jc)+
     +         rd(i,jc)*(ac(i,j+1)*pu(i,jc)+
     +                   an(i,j+1)*pc(i,jc+1))
      aswc(i,jc)=ru(i,jc-1)*aw(i,j-1)*pd(i-1,jc-1)
      asec(i,jc)=ru(i,jc-1)*ae(i,j-1)*pd(i+1,jc-1)
      anec(i,jc)=rd(i,jc)*ae(i,j+1)*pu(i+1,jc)
  300 anwc(i,jc)=rd(i,jc)*aw(i,j+1)*pu(i-1,jc)
         if(ifd.eq.9) then
      j=jred-2
      do 310 jc=2,jmc1
      j=j+2
      do 310 i=2,im1
      awc(i,jc)=awc(i,jc)+(rd(i,jc)*asw(i,j+1)+
     +                   ru(i,jc-1)*anw(i,j-1))*pc(i-1,jc)+
     +                    rc(i,jc)*(anw(i,j)*pd(i-1,jc)+
     +                              asw(i,j)*pu(i-1,jc-1))
      aec(i,jc)=aec(i,jc)+(rd(i,jc)*ase(i,j+1)+
     +                   ru(i,jc-1)*ane(i,j-1))*pc(i+1,jc)+
     +                    rc(i,jc)*(ane(i,j)*pd(i+1,jc)+
     +                              ase(i,j)*pu(i+1,jc-1))
      aswc(i,jc)=aswc(i,jc)+ru(i,jc-1)*asw(i,j-1)*pc(i-1,jc-1)+
     +                        rc(i,jc)*asw(i,j)*pd(i-1,jc-1)
      asec(i,jc)=asec(i,jc)+ru(i,jc-1)*ase(i,j-1)*pc(i+1,jc-1)+
     +                        rc(i,jc)*ase(i,j)*pd(i+1,jc-1)
      anec(i,jc)=anec(i,jc)+rd(i,jc)*ane(i,j+1)*pc(i+1,jc+1)+
     +                      rc(i,jc)*ane(i,j)*pu(i+1,jc)
  310 anwc(i,jc)=anwc(i,jc)+rd(i,jc)*anw(i,j+1)*pc(i-1,jc+1)+
     +                      rc(i,jc)*anw(i,j)*pu(i-1,jc)
         endif
         endif
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
         if(iprsol.eq.2.and.kf.eq.m.and.ifd.eq.5) then
      do 111 j=2,jm1
      do 111 i=2,im1
               write(lout,1010) im,jm,an(i,j)
               write(lout,1012) i,j,aw(i,j),ac(i,j),ae(i,j)
  111          write(lout,1011) as(i,j)
 1010 format(2(1x,i2),14x,f12.5)
 1011 format(20x,f12.5)
 1012 format(2(1x,i2),3(1x,f12.5))
          endif
         if(iprsol.eq.2.and.kf.eq.m.and.ifd.eq.9) then
      do 115 j=2,jm1
      do 115 i=2,im1
               write(lout,1017) im,jm,anw(i,j),an(i,j),ane(i,j)
               write(lout,1017) i,j,aw(i,j),ac(i,j),ae(i,j)
  115          write(lout,1016) asw(i,j),as(i,j),ase(i,j)
 1016  format(6x,3(1x,f12.5))
 1017  format(2(1x,i2),3(1x,f12.5))
          endif
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      if(iprsol.eq.2.or.iprsol.eq.3) then
      write(lout,*) ' coarse grid kf-1 ==',kf-1
      do 211 jc=2,jmc1
      do 211 i=2,im1
      write(lout,2015) im,jmc,anwc(i,jc),anc(i,jc),anec(i,jc),pd(i,jc),
     + rd(i,jc)
      write(lout,2013) i,jc,awc(i,jc),acc(i,jc),aec(i,jc),pc(i,jc),
     + rc(i,jc)
      write(lout,2014) aswc(i,jc),asc(i,jc),asec(i,jc),pu(i,jc-1),
     + ru(i,jc-1)
  211 write(lout,*) '                          m, kf-1=',m,kf-1
 2014 format(9x,3(1x,f11.5),3x,2(f11.5,1x))
 2013 format(3x,2(1x,i2),3(1x,f11.5),3x,2(f11.5,1x))
 2015 format(3x,2(1x,i2),3(1x,f11.5),3x,2(f11.5,1x))
          endif
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      return
      end
c**********************************************************************
      subroutine uoutpt(q,im,jm)
      implicit CCTK_REAL (a-h,o-z)
c**********************************************************************
c  Sample output subroutine.  Prints out the values of q at the
c  interior points of the finest grid.
c**********************************************************************
      common /io/ linp,lout
      CCTK_REAL q(im,jm)
      im1=im-1
      jm1=jm-1
      ie=1
   20 ib=ie+1
      ie=ib+min(5,im1-ib)
      do 10 j=jm1,2,-1
   10 write(lout,100) j,(q(i,j),i=ib,ie)
      if(ie.lt.im1) go to 20
  100 format(1x,i2,1x,6(1x,f10.4))
      return
      end
c**********************************************************************
      subroutine uputf(ac,aw,as,ae,an,f,nx,ny,
     + lo,nxd,nyd,i32su)
      implicit CCTK_REAL (a-h,o-z)
c**********************************************************************
      CCTK_REAL ac(lo:nxd,lo:nyd),aw(lo:nxd,lo:nyd),
     + as(lo:nxd,lo:nyd),ae(lo:nxd,lo:nyd),
     + an(lo:nxd,lo:nyd),f(lo:nxd,lo:nyd)
      CCTK_REAL a(20),b(20),ab(20)
      integer il(20),ir(20),jb(20),jt(20)
      integer ibc(4)
      common /io/ linp,lout
c
c dcell is the value assigned to the diagonal element of a dead
c cell, i.e. a cell that has 0 conections to all its neighbors.
c icendif determines the differencing scheme for the first order terms
c icendif=0 - central differencing, =1 - forward differencing.
c
      dcell = 1.
      do 321 j=lo,nyd
      do 321 i=lo,nxd
      ac(i,j)=0.
      aw(i,j)=0.
      as(i,j)=0.
      ae(i,j)=0.
  321 an(i,j)=0.
c
      nx1=nx+1
      ny1=ny+1
      read(linp,*) icendif
      read(linp,*) ibc(1),ibc(2),ibc(3),ibc(4)
      read(linp,*) nreg
      hx=1./nx
      hy=1./ny
      hx2=hx*hx
      hy2=hy*hy
      hxy2=hx2*hy2
      dcell=hxy2*dcell
      write(lout,1011) ibc(1),ibc(2),ibc(3),ibc(4),icendif,dcell,hx,hy
 1011 format(' ibc_l ibc_b ibc_r ibc_t   icendiff   dcell   hx      hy'/
     + 4x,i1,5x,i1,5x,i1,5x,i1,7x,i1,6x,e8.2,2x,f6.4,2x,f6.4/)
c
      do 10 irg=1,nreg
      read(linp,*) il(irg),ir(irg),jb(irg),jt(irg)
      read(linp,*) xk,yk,sreg,freg
      read(linp,*) a(irg),b(irg),ab(irg)
      write(lout,1000) il(irg),ir(irg),jb(irg),jt(irg),xk,yk,sreg,freg
 1000 format(1x,i3,',',i3,' X ',i3,',',i3,2x,1pe8.2,1x,1pe8.2,2x,
     + 1pe8.1,1x,1pe8.1)
      write(lout,1001) a(irg),b(irg),ab(irg)
 1001 format(17x,' a=',1pe10.3,'  b=',1pe10.3,' ab=',1pe10.3)
      xk=xk*hy2
      yk=yk*hx2
      sreg=sreg*hxy2
      freg=freg*hxy2
      a(irg)=a(irg)*hx2*hy
      b(irg)=b(irg)*hx*hy2
      ab(irg)=ab(irg)*hx*hy
         if(icendif.eq.0) then
      a(irg)=a(irg)/2.
      b(irg)=b(irg)/2.
      ab(irg)=ab(irg)/4.
         endif
      if(il(irg).eq.1) il(irg)=0
      if(ir(irg).eq.nx) ir(irg)=nx1
      if(jb(irg).eq.1) jb(irg)=0
      if(jt(irg).eq.ny) jt(irg)=ny1
      do 20 i=il(irg),ir(irg)
      do 20 j=jb(irg),jt(irg)
      aw(i,j)=xk
      as(i,j)=yk
      ac(i,j)=sreg
   20 f(i,j)=freg
   10 continue
      write(lout,*) ' - - - - - - - - - - - - - - - - - - - - - - - - -'
c                             defining coeficients by harmonic averaging
      do 30 i=1,nx
      asio=as(i,0)
      do 30 j=1,ny1
      aa=as(i,j)*asio
         if(aa.gt.0.) then
      t=2.*aa/(as(i,j)+asio)
      asio=as(i,j)
      as(i,j)=t
         else
      asio=as(i,j)
      as(i,j)=0.
         endif
   30 continue
      do 40 j=1,ny
      awoj=aw(0,j)
      do 40 i=1,nx1
      aa=aw(i,j)*awoj
         if(aa.gt.0.) then
      t=2.*aa/(aw(i,j)+awoj)
      awoj=aw(i,j)
      aw(i,j)=t
         else
      awoj=aw(i,j)
      aw(i,j)=0.
         endif
   40 continue
      do 45 i=0,nx
      do 45 j=0,ny
      ae(i,j)=aw(i+1,j)
   45 an(i,j)=as(i,j+1)
      do 50 i=1,nx
      do 50 j=1,ny
      ac(i,j)=ac(i,j)-aw(i,j)-as(i,j)-ae(i,j)-
     + an(i,j)
   50 if(ac(i,j).eq.0.) ac(i,j)=dcell
c                                        adding on the unsymmetric terms
      do 51 irg=1,nreg
c                                           icendif=0 ==> central diff g
         if(icendif.eq.0) then
      do 52 i=il(irg),ir(irg)
      do 52 j=jb(irg),jt(irg)
      aw(i,j)=aw(i,j)-a(irg)
      ae(i,j)=ae(i,j)+a(irg)
      an(i,j)=an(i,j)+b(irg)
   52 as(i,j)=as(i,j)-b(irg)
c                                        icendif=1 ==> upstream diff s
         elseif(icendif.eq.1) then
      do 54 i=il(irg),ir(irg)
      do 54 j=jb(irg),jt(irg)
      ac(i,j)=ac(i,j)-a(irg)
      ae(i,j)=ae(i,j)+a(irg)
      an(i,j)=an(i,j)+b(irg)
   54 ac(i,j)=ac(i,j)-b(irg)
         endif
   51 continue
c                           set boundary conditions for 5 point operator
      do 60 j=1,ny
c                                                          left boundary
      ae(0,j)=aw(1,j)
         if(ibc(1).eq.1) then
      ac(0,j)=aw(1,j)
      f(0,j)=2.*aw(1,j)*0.0
         elseif(ibc(1).eq.2) then
      ac(0,j)=-aw(1,j)
      f(0,j)=hx*aw(1,j)*0.0
         endif
         if(ac(0,j).eq.0.) then
      ae(0,j)=0.
      ac(0,j)=dcell
      f(0,j)=0.
         endif
c                                                         right boundary
      aw(nx1,j)=ae(nx,j)
         if(ibc(3).eq.1) then
      ac(nx1,j)=ae(nx,j)      
      f(nx1,j)=2.*ae(nx,j)*0.
         elseif(ibc(3).eq.2) then
      ac(nx1,j)=-ae(nx,j)
      f(nx1,j)=hx*ae(nx,j)*0.
         endif
         if(ac(nx1,j).eq.0.) then
      aw(nx1,j)=0.
      ac(nx1,j)=dcell
      f(nx1,j)=0.
         endif
   60 continue
c
      do 80 i=1,nx
c                                                         lower boundary
      an(i,0)=as(i,1)
         if(ibc(2).eq.1) then
      ac(i,0)=as(i,1)
      f(i,0)=2.*as(i,1)*0.
         elseif(ibc(2).eq.2) then
      ac(i,0)=-as(i,1)
      f(i,0)=hy*as(i,1)*0.
         endif
         if(ac(i,0).eq.0.) then
      an(i,0)=0.
      ac(i,0)=dcell
      f(i,0)=0.
         endif
c                                                         upper boundary
      as(i,ny1)=an(i,ny)
         if(ibc(4).eq.1) then
      ac(i,ny1)=an(i,ny)
      f(i,ny1)=2.*an(i,ny)*0.
         elseif(ibc(4).eq.2) then
      ac(i,ny1)=-an(i,ny)
      f(i,ny1)=2.*an(i,ny)*0.
         endif
         if(ac(i,ny1).eq.0.) then
      as(i,ny1)=0.
      ac(i,ny1)=dcell
      f(i,ny1)=0.
         endif
   80 continue
c                     connections between ghost boundary points zeroed
      do 83 j=1,ny1
      as(0,j)=0.
      an(0,j-1)=0.
      as(nx1,j)=0.
   83 an(nx1,j-1)=0.
      do 86 i=1,nx1
      aw(i,0)=0.
      ae(i-1,0)=0.
      aw(i,ny1)=0.
   86 ae(i-1,ny1)=0.
c                                        corner stencils and rhs defined
         if(i32su.eq.32) then
      do 90 j=0,ny1,ny1
      do 90 i=0,nx1,nx1
      ac(i,j)=dcell
      aw(i,j)=0.
      ae(i,j)=0.
      as(i,j)=0.
      an(i,j)=0.
   90 f(i,j)=0.
         endif
c                                i32su=22 - boundary conditions absorbed
         if(i32su.eq.22) then
      do 100 j=1,ny
      awac=aw(1,j)/ac(0,j)
      ac(1,j)=ac(1,j)-awac*ae(0,j)
      aw(1,j)=0.
      f(1,j)=f(1,j)-awac*f(0,j)
      ac(0,j)=0.
      ae(0,j)=0.
      f(0,j)=0
      awac=aw(nx1,j)/ac(nx1,j)
      ac(nx,j)=ac(nx,j)-awac*ae(nx1,j)
      ae(nx,j)=0.
      f(nx,j)=f(nx,j)-awac*f(nx1,j)
      ac(nx1,j)=0.
      aw(nx1,j)=0.
  100 f(nx1,j)=0.
c
      do 110 i=1,nx
      asac=as(i,1)/ac(i,0)
      ac(i,1)=ac(i,1)-asac*an(i,0)
      as(i,1)=0.
      f(i,1)=f(i,1)-asac*f(i,0)
      ac(i,0)=0.
      an(i,0)=0.
      f(i,0)=0.
      anac=an(i,ny)/ac(i,ny1)
      ac(i,ny)=ac(i,ny)-anac*as(i,ny1)
      an(i,ny)=0.
      f(i,ny)=f(i,ny)-anac*f(i,ny1)
      ac(i,ny1)=0.
      as(i,ny1)=0.
  110 f(i,ny1)=0.
         endif
      return
      end
