c/*@@
c  @file      Extract.F
c  @date      28th October 1997
c  @author    Gab Allen
c  @desc      Waveform extraction
c
c  @enddesc
c@@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

c/*@@
c  @routine    Extract
c  @date       28th October 1997
c  @author     Gab Allen
c  @desc       Entry point for waveform extraction routines
c
c  @enddesc
c  @calls      D3_extract
c  @calledby   No idea
c  @history
c
c  @endhistory
c@@*/


      SUBROUTINE Extract(CCTK_ARGUMENTS)

      USE D3_extract_int

      IMPLICIT NONE

      DECLARE_CCTK_ARGUMENTS
      DECLARE_CCTK_PARAMETERS
      DECLARE_CCTK_FUNCTIONS

c     Non-Cactus input variables for D3_extract

      INTEGER ::
     &     igrid,lmode,mmode,Psi_power=1
      INTEGER ::
     &     do_ADMmass(2)
      CCTK_REAL ::
     &     orig(3),radius
      CCTK_REAL ::
     &     x_1d(cctk_lsh(1)),y_1d(cctk_lsh(2)),z_1d(cctk_lsh(3))


c     Output variables from D3_extract

      CCTK_REAL ::
     &     dtaudt,ADMmass(2),mass,rsch,momentum(3),spin(3)
      CCTK_REAL,ALLOCATABLE ::
     &     Qodd(:,:,:), Qeven(:,:,:)


c     Local variables

      INTEGER ::
     &     ix,iy,iz,fn1,fn2,lmin,lmax,mmin,mmax,lstep,mstep,
     &     il,im,it,ndet,idet,ioutput
      INTEGER,SAVE :: openfile = 1
      INTEGER,PARAMETER ::
     &     max_detectors = 9,number_timecoords = 2
      CCTK_REAL,PARAMETER ::
     &     two = 2.0D0
      CCTK_REAL ::
     &     time,r1,r2,dr,r_max,xmin,xmax,ymin,ymax,zmin,zmax
      CCTK_REAL ::
     &     detector(max_detectors)
      CHARACTER*100 ::
     &     out1,out2,massfile,rschfile,ADMmassfile1,ADMmassfile2,
     &     momentumfile1,momentumfile2,momentumfile3,spinfile1,spinfile2,spinfile3
      CHARACTER*3 ::
     &     timestring

      LOGICAL :: use_proper,use_coordinate,do_output
      CCTK_REAL,DIMENSION(10),SAVE :: proper_time
      CCTK_REAL,SAVE :: last_time
      INTEGER, save :: open_file_level(10) = 0

      INTEGER myproc,ierr
      CCTK_REAL dx,dy,dz
      CCTK_INT :: nx, ny, nz

      character*80 infoline

c     ------------------------------------------------------------------

      myproc = CCTK_MyProc(cctkGH)
      dx = cctk_delta_space(1)
      dy = cctk_delta_space(2)
      dz = cctk_delta_space(3)

      nx = cctk_lsh(1)
      ny = cctk_lsh(2)
      nz = cctk_lsh(3)

c     ------------------------------------------------------------------
c
c     0. Check to see if Extract should have been called
c
c     ------------------------------------------------------------------

      IF ((itout) .LE. 0) RETURN
      IF(MOD(int(cctk_iteration),int(itout)) .NE. 0) RETURN

      if (verbose == 1) then
         write(infoline,'(A24,G14.7)')
     &        'Calling Extract at time ',cctk_time
         call CCTK_INFO(infoline)
      end if

c     ------------------------------------------------------------------
c
c     1. Initial Stuff
c
c     ------------------------------------------------------------------


c     Get the number of polar divisions

      IF (MOD(int(Nt),2) == 1 .OR. MOD(int(Np),2) == 1 .OR.
     &     Nt < 0 .OR. Np < 0) THEN
        call CCTK_WARN(0,"Error in Nt or Np in Extract")
      END IF


c     Get the origin of spherical symmetry

      orig(1) = origin_x
      orig(2) = origin_y
      orig(3) = origin_z


c     Set the value of igrid

      IF (CCTK_EQUALS(domain,"octant")) THEN
         igrid = 1
      ELSEIF (CCTK_EQUALS(domain,"bitant")) THEN
         igrid = 3
      ELSE
         igrid = 0
      END IF

#ifdef THORN_CARTOON_2D
      IF (contains("cartoon_active","yes") .NE. 0) THEN
         igrid = 2
         Np = 1
      END IF
#endif

c     Create 1D coordinate arrays

      do ix = 1, nx
        x_1d(ix) = x(ix,1,1)
      enddo
      do iy = 1, ny
        y_1d(iy) = y(1,iy,1)
      enddo
      do iz = 1, nz
        z_1d(iz) = z(1,1,iz)
      enddo

c     See if the ADM mass should be calculated and output
      IF (doADMmass == 1) THEN
         do_ADMmass(1) = 1
         IF (conformal_state == 1) THEN
            do_ADMmass(2) = 1
         ELSE
            do_ADMmass(2) = 0
         ENDIF
      ELSE
         do_ADMmass(:) = 0
      END IF

c     Set array containing gtt ... for the moment I am lazy and
c     just make it the lapse and do not include the shift!
      g00 = alp

c     What kind of timecoordinate do I want to use
      if (CCTK_EQUALS(timecoord,"both")) then
         use_proper = .true.
         use_coordinate = .true.
         if (cctk_iteration == 0) then
            proper_time = 0
            last_time = cctk_time
         end if
      elseif (CCTK_EQUALS(timecoord,"proper")) then
         use_proper = .true.
         use_coordinate = .false.
         if (cctk_iteration == 0) then
            proper_time = 0
            last_time = cctk_time
         end if
      elseif (CCTK_EQUALS(timecoord,"coord")) then
         use_proper = .false.
         use_coordinate = .true.
      else
         use_proper = .false.
         use_coordinate = .true.
      endif

c     ------------------------------------------------------------------
c
c     2. Sort out which modes to extract
c
c     ------------------------------------------------------------------

      lmode  = l_mode

      IF (all_modes == 1) THEN
         lmin = 2
         lmax = lmode
c     Get m_mode in this case too, or the T3E has a junk value (PW)
         mmode = 0
         ALLOCATE(Qodd(2,2:lmode,0:lmode),Qeven(2,2:lmode,0:lmode))
      ELSE
         lmin = lmode
         lmax = lmode
         mmode = m_mode
         ALLOCATE(Qodd(2,2:lmode,0:mmode),Qeven(2,2:lmode,0:mmode))
      ENDIF
c      print*,'Extract Qodd'
c      print*,lbound(Qodd,1),lbound(Qodd,2),lbound(Qodd,3)
c      print*,size(Qodd,1),size(Qodd,2),size(Qodd,3)


c     Check do not have bad values for the modes

      IF (lmode < 2 .OR. mmode > lmode .OR. mmode < 0) THEN
         call CCTK_WARN(0,"Error in lmode,mmode in Extract")
      END IF
      IF (lmode > 9) THEN
         call CCTK_WARN(4,"Filenames will probably be crazy in Extract")
      END IF

c     If using an octant, do not do odd modes

      IF (igrid == 1 .OR. igrid == 2) THEN
         lstep = 2 ; mstep = 2
      ELSE
         lstep = 1 ; mstep = 1
      END IF

c     If we have cartoon, then jump past all the m-modes in loop
      IF (igrid == 2) mstep = lmode+1


c     ------------------------------------------------------------------
c
c     3. Find maximum radius of extraction on grid
c
c     ------------------------------------------------------------------

      call CCTK_CoordRange(ierr,cctkGH,xmin,xmax,-1,"x","cart3d")
      call CCTK_CoordRange(ierr,cctkGH,ymin,ymax,-1,"y","cart3d")
      call CCTK_CoordRange(ierr,cctkGH,zmin,zmax,-1,"z","cart3d")

      IF (igrid == 2) THEN
         r_max = MIN(xmax-two*Dx-orig(1),
     &        zmax-two*Dz-orig(3))
      ELSE IF (igrid == 1) THEN
         r_max = MIN(zmax-two*Dx-orig(1),
     &        ymax-two*Dy-orig(2),
     &        zmax-two*Dz-orig(3))
      ELSE IF (igrid == 3) THEN
         r_max = MIN(xmax-two*Dx-orig(1),
     &        ymax-two*Dy-orig(2),
     &        zmax-two*Dz-orig(3),
     &        ABS(xmin+two*Dx-orig(1)),
     &        ABS(ymin+two*Dy-orig(2)))
      ELSE
         r_max = MIN(xmax-two*Dx-orig(1),
     &        ymax-two*Dy-orig(2),
     &        zmax-two*Dz-orig(3),
     &        ABS(xmin+two*Dx-orig(1)),
     &        ABS(ymin+two*Dy-orig(2)),
     &        ABS(zmin+two*Dz-orig(3)))
      END IF




c     ------------------------------------------------------------------
c
c     4. Extract Cauchy initial data for a linear wave equation
c
c     ------------------------------------------------------------------

      test_Cauchy: IF (Cauchy == 1) THEN

      it = Cauchy_timestep


c     check_timestep: IF (cctk_iteration == it) THEN

      check_timestep: IF (open_file_level(1) .eq. 0) THEN

c           Open output files
            write (*,*) "Open output file on level ", 1
            open_file_level(1)=1;

            test_myproc1: IF (CCTK_MyProc(cctkGH) == 0 ) THEN

               out1 = "rsch_ini.rl"
               call Extract_open(cctkGH,0,out1,77)

               out1 = "mass_ini.rl"
               call Extract_open(cctkGH,0,out1,78)

               IF (do_ADMmass(1) == 1) THEN
                  out1 = "ADMmass_ini.rl"
                  call Extract_open(cctkGH,0,out1,79)
               ENDIF
               IF (do_ADMmass(2) == 1) THEN
                  out1 = "ADMmassc_ini.rl"
                  call Extract_open(cctkGH,0,out1,80)
               ENDIF

               IF (do_momentum == 1) THEN
                  out1 = "momentum_x_ini.rl"
                  call Extract_open(cctkGH,0,out1,781)
                  out1 = "momentum_y_ini.rl"
                  call Extract_open(cctkGH,0,out1,782)
                  out1 = "momentum_z_ini.rl"
                  call Extract_open(cctkGH,0,out1,783)
               ENDIF
               IF (do_spin == 1) THEN
                  out1 = "spin_x_ini.rl"
                  call Extract_open(cctkGH,0,out1,784)
                  out1 = "spin_y_ini.rl"
                  call Extract_open(cctkGH,0,out1,785)
                  out1 = "spin_z_ini.rl"
                  call Extract_open(cctkGH,0,out1,786)
               ENDIF

               loop_l1: DO il = lmin,lmax,lstep

                  IF (all_modes == 0) THEN
                     mmin = mmode ; mmax = mmode
                  ELSE
                     mmin = 0 ; mmax = il
                  END IF

                  loop_m1: DO im = mmin,mmax,lstep

                     fn1 = il*10+im
                     out1 = "Qeven_ini_"//
     &                    CHAR(il+48)//CHAR(im+48)//".rl"
                     call Extract_open(cctkGH,0,out1,fn1)

c                    Only print odd-parity if full grid
                     IF (igrid == 0) THEN
                        fn2 = 100+il*10+im
                        out2 = "Qodd_ini_"//
     &                       CHAR(il+48)//CHAR(im+48)//".rl"
                        call Extract_open(cctkGH,0,out2,fn2)
                     END IF

                  END DO loop_m1

               END DO loop_l1

            END IF test_myproc1


c           Find range of extraction radii

            r1 = Cauchy_r1
            r2 = r_max
            dr = Cauchy_dr


c           Do extraction at each radius

            IF (verbose == 1) THEN
               WRITE(*,*) "Extracting Cauchy initial data"
               WRITE(*,*) "  r = ",r1," to ",r2," step ",dr
            END IF

            radius = r1

            extract_at_each_radius: DO WHILE (radius < r2)

c               print*,'Calling D3_extract'
c               print*,nx,ny,nz,Nt,Np
               CALL D3_extract(cctkGH,conformal_state,do_ADMmass,do_momentum,do_spin,
     &           igrid,orig,myproc,interpolation_operator,interpolation_order,
     &           Nt,Np,nx,ny,nz,all_modes,lmode,
     &           mmode,x_1d,y_1d,z_1d,Dx,Dy,Dz,Psi_power,Psi,g00,
     &           gxx,gxy,gxz,gyy,gyz,gzz,kxx,kxy,kxz,kyy,kyz,kzz,
     &           radius,ADMmass,momentum,spin,mass,rsch,
     &           Qodd,Qeven,temp3d,dtaudt)
c               print*,Qodd,Qeven
c               print*,'Exited D3_extract'

               IF (verbose == 1) THEN
                  WRITE(*,*) "Extracted at r =",radius
                  WRITE(*,*) " Sch radius/mass =",rsch,mass
               ENDIF

c              Write to file at each radius

               test_myproc2: IF (CCTK_MyProc(cctkGH) == 0) THEN

                  WRITE( 77,*) radius,rsch
                  WRITE( 78,*) radius,mass
                  IF (do_ADMmass(1) == 1) WRITE( 79,*) radius,ADMmass(1)
                  IF (do_ADMmass(2) == 1) WRITE( 80,*) radius,ADMmass(2)
                  IF (do_momentum == 1) then
                     WRITE(781,*) radius,momentum(1)
                     WRITE(782,*) radius,momentum(2)
                     WRITE(783,*) radius,momentum(3)
                  end if

                  IF (do_spin == 1) THEN
                     WRITE(784,*) radius,spin(1)
                     WRITE(785,*) radius,spin(2)
                     WRITE(786,*) radius,spin(3)
                  END IF

                  loop_l2: DO il = lmin,lmax,lstep

                     IF (all_modes == 0) THEN
                        mmin = mmode ; mmax = mmode
                     ELSE
                        mmin = 0 ; mmax = il
                     END IF

                     loop_m2: DO im = mmin,mmax,mstep

                        fn1 = il*10+im
                        WRITE(fn1,*) rsch,Qeven(:,il,im)

c                       Only print odd-parity if full grid
                        IF (igrid == 0) THEN
                           fn2 = 100+il*10+im
                           WRITE(fn2,*) rsch,Qeven(:,il,im)
                        END IF

                     END DO loop_m2

                  END DO loop_l2

               END IF test_myproc2

               radius = radius+dr

            END DO extract_at_each_radius


c           Now close all the files

            test_myproc3: IF (myproc == 0) THEN

               CLOSE( 77)   ! Schwarzschild radius
               CLOSE( 78)   ! Schwarzschild mass
               IF (do_ADMmass(1) == 1) CLOSE( 79)   ! ADM mass
               IF (do_ADMmass(2) == 1) CLOSE( 80)   ! ADM mass

               IF (do_momentum == 1) THEN
                  CLOSE(781)    ! momentum
                  CLOSE(782)    ! momentum
                  CLOSE(783)    ! momentum
               end if

               IF (do_spin == 1) THEN
                  CLOSE(784)    ! spin
                  CLOSE(785)    ! spin
                  CLOSE(786)    ! spin
               END IF

               loop_l3: DO il = lmin,lmax,lstep

                  IF (all_modes == 0) THEN
                     mmin = mmode ; mmax = mmode
                  ELSE
                     mmin = 0 ; mmax = il
                  END IF

                  loop_m3: DO im = mmin,mmax,mstep

                     fn1 = il*10+im
                     CLOSE(fn1) ! Qeven_ini_lm.dat

                     IF (igrid == 0) THEN
                        fn2 = 100+il*10+im
                        CLOSE(fn2) ! Qodd_ini_lm.dat
                     END IF

                  END DO loop_m3

               END DO loop_l3

            END IF test_myproc3

         END IF check_timestep

      END IF test_Cauchy


c     ------------------------------------------------------------------
c
c     5. Extract waveforms at detectors
c
c     ------------------------------------------------------------------

c     Cannot use the conformal equation for ADM mass now
      do_ADMmass(2) = 0

      ndet = num_detectors

      test_detectors: IF (ndet > 0) THEN

         IF (ndet > 9) THEN
            call CCTK_WARN(0,"Too many detectors in Extract")
         END IF

         detector(1) = detector1
         detector(2) = detector2
         detector(3) = detector3
         detector(4) = detector4
         detector(5) = detector5
         detector(6) = detector6
         detector(7) = detector7
         detector(8) = detector8
         detector(9) = detector9

         DO idet = ndet+1, max_detectors
            detector(idet) = 0.0D0
         END DO

         DO idet = 1, ndet
            IF (detector(idet) > r_max) THEN
               IF (openfile == 1 .OR. cctk_iteration == it) THEN
                  call CCTK_WARN(8,"Removing detectors outside coordinate range")
               END IF
               ndet =idet-1
               GOTO 404
            ENDIF
         END DO
 404     CONTINUE

         IF (verbose == 1) THEN
            IF (openfile == 1) THEN
               WRITE(*,*) "Extracting at ",ndet," detectors "
               DO idet=1,ndet
                  WRITE(*,*) "  r = ",detector(idet)
               ENDDO
            END IF
         END IF

         detector_loop: DO idet = 1,ndet

            radius = detector(idet)

            IF (verbose == 1) THEN
              WRITE(*,*) "Calling extract at radius ... ",radius
            END IF

            CALL D3_extract(cctkGH,conformal_state,do_ADMmass,do_momentum,do_spin,
     &           igrid,orig,myproc,interpolation_operator,interpolation_order,
     &           Nt,Np,nx,ny,nz,all_modes,
     &           lmode,mmode,x_1d,y_1d,z_1d,
     &           Dx,Dy,Dz,Psi_power,Psi,g00,gxx,gxy,gxz,gyy,gyz,gzz,kxx,kxy,kxz,
     &           kyy,kyz,kzz,radius,ADMmass,momentum,spin,mass,rsch,
     &           Qodd,Qeven,temp3d,dtaudt)

            IF (verbose == 1) THEN
                   WRITE(*,*) "Detector ",idet," ..."
                   WRITE(*,*) " Sch radius/mass =",rsch,mass
            END IF

c           Output to files

            test_myproc4: IF (myproc == 0) THEN

               DO ioutput=1,number_timecoords

                  do_output = .false.
                  timestring = "ERR"
                  time = 0D0
                  IF (ioutput == 1 .AND. use_coordinate) THEN
                     timestring = ".tl"
                     do_output = .true.
                     time = cctk_time
                  ELSEIF (ioutput == 2 .AND. use_proper) THEN
                     timestring = ".ul"
                     do_output = .true.
                     time = proper_time(idet) + (cctk_time-last_time)*dtaudt
                     proper_time(idet) = time
                     if (idet == ndet) last_time = cctk_time
                  ENDIF

                  IF (do_output) THEN

c                   Output extracted radius and mass
                     rschfile = "rsch_"//"R"//
     &                    CHAR(idet+48)//timestring
                     massfile = "mass_"//"R"//
     &                    CHAR(idet+48)//timestring
                     ADMmassfile1 = "ADMmass_"//"R"//
     &                    CHAR(idet+48)//timestring
                     ADMmassfile2 = "ADMmassc_"//"R"//
     &                    CHAR(idet+48)//timestring
                     momentumfile1 = "momentum_x_"//"R"//
     &                    CHAR(idet+48)//timestring
                     momentumfile2 = "momentum_y_"//"R"//
     &                    CHAR(idet+48)//timestring
                     momentumfile3 = "momentum_z_"//"R"//
     &                    CHAR(idet+48)//timestring
                     spinfile1 = "spin_x_"//"R"//
     &                    CHAR(idet+48)//timestring
                     spinfile2 = "spin_y_"//"R"//
     &                    CHAR(idet+48)//timestring
                     spinfile3 = "spin_z_"//"R"//
     &                    CHAR(idet+48)//timestring

                     call Extract_write(cctkGH,openfile,rschfile,time,rsch)
                     call Extract_write(cctkGH,openfile,massfile,time,mass)
                     IF (do_ADMmass(1) == 1) THEN
                        call Extract_write(cctkGH,openfile,ADMmassfile1,time,ADMmass(1))
                     END IF
                     IF (do_ADMmass(2) == 1) THEN
                        call Extract_write(cctkGH,openfile,ADMmassfile2,time,ADMmass(2))
                     END IF
                     IF (do_momentum == 1) THEN
                        call Extract_write(cctkGH,openfile,momentumfile1,time,momentum(1))
                        call Extract_write(cctkGH,openfile,momentumfile2,time,momentum(2))
                        call Extract_write(cctkGH,openfile,momentumfile3,time,momentum(3))
                     END IF
                     IF (do_spin == 1) THEN
                        call Extract_write(cctkGH,openfile,spinfile1,cctk_time,spin(1))
                        call Extract_write(cctkGH,openfile,spinfile2,cctk_time,spin(2))
                        call Extract_write(cctkGH,openfile,spinfile3,cctk_time,spin(3))
                     END IF


c                    Output gauge invariant variables

                     loop_l5: DO il = lmin,lmax,lstep

                       IF (all_modes == 0) THEN
                          mmin = mmode ; mmax = mmode
                       ELSE
                          mmin = 0 ; mmax = il
                       END IF

                       loop_m5: DO im = mmin,mmax,mstep

                         out1 = "Qeven_"//"R"
     &                        //CHAR(idet+48)//"_"//CHAR(il+48)//
     &                        CHAR(im+48)//timestring
                         out2 = "Qodd_"//"R"
     &                        //CHAR(idet+48)//"_"//CHAR(il+48)//
     &                        CHAR(im+48)//timestring

c                        Write even parity waveforms
                         call Extract_write2(cctkGH,openfile,out1,time,Qeven(:,il,im))

c                        Only write odd parity waveforms if full grid
                         IF (igrid == 0) THEN
                            call Extract_write2(cctkGH,openfile,out2,time,Qodd(:,il,im))
                         END IF

                      END DO loop_m5

                   END DO loop_l5

                ENDIF
             END DO
          END IF test_myproc4

       END DO detector_loop

      END IF test_detectors

      DEALLOCATE(Qodd,Qeven)

      openfile = 0

      END SUBROUTINE Extract


      SUBROUTINE Extract_open(cctkGH,openfile,filename,filehandle)

      implicit none

      DECLARE_CCTK_PARAMETERS

      integer filehandle,openfile
      CCTK_INT nchar
      character*(*) filename
      CCTK_POINTER cctkGH

      character*200 fullname
      character*200 filestr

c     Get the output directory sorted out
      call CCTK_FortranString(nchar,out_dir,filestr)
      fullname = filestr(1:nchar)//"/"//filename

      if (openfile==1) then
         OPEN(UNIT= filehandle,FILE=fullname,STATUS="unknown")
         write(filehandle,101) char(34),filename
 101     format(a1,a14)
         call Extract_Advertise(cctkGH,fullname)
      else
         OPEN(UNIT=filehandle,FILE=fullname,STATUS="old",
     .        POSITION="append")
      end if


      END SUBROUTINE Extract_open


      SUBROUTINE Extract_write(cctkGH, openfile,filename,value1,value2)

      implicit none

      DECLARE_CCTK_PARAMETERS

      CCTK_POINTER cctkGH
      CCTK_REAL value1,value2
      integer openfile
      CCTK_INT nchar
      character*(*) filename

      character*200 fullname
      character*200 filestr

c     Get the output directory sorted out
      call CCTK_FortranString(nchar,out_dir,filestr)
      fullname = filestr(1:nchar)//"/"//filename

      if (openfile==1) then
         OPEN(UNIT= 444,FILE=fullname,STATUS="unknown")
         write(444,101) char(34),filename
 101     format(a1,a14)
         call Extract_Advertise(cctkGH,fullname)
      else
         OPEN(UNIT=444,FILE=fullname,STATUS="old",POSITION="append")
      end if

      write(444,*) value1,value2

      close(444)

      END SUBROUTINE Extract_write


      SUBROUTINE Extract_write2(cctkGH, openfile,filename,value1,value2)

      implicit none

      DECLARE_CCTK_PARAMETERS

      CCTK_POINTER cctkGH
      CCTK_REAL value1
      CCTK_REAL value2(2)
      integer openfile
      CCTK_INT nchar
      character*(*) filename

      character*200 fullname
      character*200 filestr

c     Get the output directory sorted out
      call CCTK_FortranString(nchar,out_dir,filestr)
      fullname = filestr(1:nchar)//"/"//filename

      if (openfile==1) then
         OPEN(UNIT= 444,FILE=fullname,STATUS="unknown")
         write(444,101) char(34),filename
 101     format(a1,a14)
         call Extract_Advertise(cctkGH,fullname)
      else
         OPEN(UNIT=444,FILE=fullname,STATUS="old",POSITION="append")
      end if

      write(444,*) value1,value2(1),value2(2)

      close(444)

      END SUBROUTINE Extract_write2
