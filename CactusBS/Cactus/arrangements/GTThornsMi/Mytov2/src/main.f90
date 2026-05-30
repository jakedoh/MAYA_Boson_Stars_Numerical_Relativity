
  program main

  use global_numbers

  implicit none

  integer Nxx,Ntt
  integer every_0Dt, every_1Dt

  real(kind=8) dt_temp


! --->   GET PARAMETERS   <---

    Namelist /TOV_Input/ nvars, &
                      xmin, xmax, &
                      res_num, &
                      g_cells, &
                      Nxx, courant, Ntt, &
                      boundary_type, &
                      amplitude, sigma, x0, &
                      every_0Dt, every_1Dt, &
                      diss_coef, gamma, floor, limiter, &
                      num_method, &
                      alpha_c, rho0_c, K_poly,EoS

! Dafault values
     res_num         = 1
     g_cells         = 0
     nvars           = 3
     xmin            = -1.0
     xmax            = 1.0
     Nxx             = 1000
     Ntt             = 2000
     courant         = 0.25
     boundary_type   = 2
     amplitude       = 1.0
     sigma           = 0.1
     x0              = 0.
     every_0Dt       = 100
     every_1Dt       = 100
     diss_coef       = 0.0


    open (3, file='iii.par', status = 'old' )
    read (3, nml = TOV_Input)
    close(3)

! Some numbers

   zero  = 0.0D0
   third = 1.0D0/3.0D0
   half  = 0.5D0
   one   = 1.0D0
   two   = 2.0D0
   pii = 4.0d0*atan(1.0d0)

    print *, 'pi=',pii
    print *, 'res_num=', res_num         
    print *, 'g_cells=', g_cells         
    print *, 'nvars=', nvars           
    print *, 'xmin=', xmin            
    print *, 'xmax=', xmax            
    print *, 'Nxx=', Nxx             
    print *, 'Ntt=', Ntt             
    print *, 'courant=', courant         
    print *, 'boundary_type=', boundary_type   
    print *, 'amplitude=', amplitude       
    print *, 'sigma=', sigma           
    print *, 'x0=', x0              
    print *, 'every_0Dt=', every_0Dt       
    print *, 'every_1Dt=', every_1Dt       
    print *, 'diss_coef=', diss_coef       


! --->   EVOLUTION   <---


      Nx = 2**(res_num-1)*Nxx
      Nt = 2**(res_num-1)*Ntt
      every_0D = 2**(res_num-1)*every_0Dt
      every_1D = 2**(res_num-1)*every_1Dt

      call evolve( )


! --->   END   <---

  print *
  print *, 'PROGRAM TOV HAS FINISHED'
  print *

  end program main

