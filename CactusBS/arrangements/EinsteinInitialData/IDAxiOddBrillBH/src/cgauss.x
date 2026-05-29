        gtil = amp*(exp(-((-eta0 + eta(i,j,k))**2/sigma**2)) + exp(-((eta0 +
     &   eta(i,j,k))**2/sigma**2)))
        dngtil = amp*((-2.00000000000000d0*(-eta0 + eta(i,j,k))*exp(-((-eta0
     &   + eta(i,j,k))**2/sigma**2)))/sigma**2 - (2.00000000000000d0*(eta0 +
     &   eta(i,j,k))*exp(-((eta0 + eta(i,j,k))**2/sigma**2)))/sigma**2)
        dnngtil = amp*((4.0000000000000d0*(-eta0 + eta(i,j,k))**2*exp(-((-et
     &  a0 + eta(i,j,k))**2/sigma**2)))/sigma**4 - (2.00000000000000d0*exp(-
     &  ((-eta0 + eta(i,j,k))**2/sigma**2)))/sigma**2 + (4.0000000000000d0*(
     &  eta0 + eta(i,j,k))**2*exp(-((eta0 + eta(i,j,k))**2/sigma**2)))/sigma**4 
     &  - (2.00000000000000d0*exp(-((eta0 + eta(i,j,k))**2/sigma**2)))/sigma
     &  **2)
        dnnngtil = amp*((-8.0000000000000d0*(-eta0 + eta(i,j,k))**3*exp(-((-
     &  eta0 + eta(i,j,k))**2/sigma**2)))/sigma**6 + (1.20000000000000d1*(-e
     &  ta0 + eta(i,j,k))*exp(-((-eta0 + eta(i,j,k))**2/sigma**2)))/sigma**4 - (
     &  8.0000000000000d0*(eta0 + eta(i,j,k))**3*exp(-((eta0 + eta(i,j,k))**2/si
     &  gma**2)))/sigma**6 + (1.20000000000000d1*(eta0 + eta(i,j,k))*exp(-((
     &  eta0 + eta(i,j,k))**2/sigma**2)))/sigma**4)
        dnnnngtil = amp*((1.60000000000000d1*(-eta0 + eta(i,j,k))**4*exp(-((
     &  -eta0 + eta(i,j,k))**2/sigma**2)))/sigma**8 - (4.8000000000000d1*(-e
     &  ta0 + eta(i,j,k))**2*exp(-((-eta0 + eta(i,j,k))**2/sigma**2)))/sigma**6 
     &  + (1.20000000000000d1*exp(-((-eta0 + eta(i,j,k))**2/sigma**2)))/sigm
     &  a**4 + (1.60000000000000d1*(eta0 + eta(i,j,k))**4*exp(-((eta0 + eta(i,j,k)
     &  )**2/sigma**2)))/sigma**8 - (4.8000000000000d1*(eta0 + eta(i,j,k))
     &  **2*exp(-((eta0 + eta(i,j,k))**2/sigma**2)))/sigma**6 + (1.200000000
     &  00000d1*exp(-((eta0 + eta(i,j,k))**2/sigma**2)))/sigma**4)
        dnnnnngtil = amp*((-3.2000000000000d1*(-eta0 + eta(i,j,k))**5*exp(-(
     &  (-eta0 + eta(i,j,k))**2/sigma**2)))/sigma**10 + (1.60000000000000d2*
     &  (-eta0 + eta(i,j,k))**3*exp(-((-eta0 + eta(i,j,k))**2/sigma**2)))/sigma*
     &  *8 - (1.20000000000000d2*(-eta0 + eta(i,j,k))*exp(-((-eta0 + eta(i,j,k))
     &  **2/sigma**2)))/sigma**6 - (3.2000000000000d1*(eta0 + eta(i,j,k))**5
     &  *exp(-((eta0 + eta(i,j,k))**2/sigma**2)))/sigma**10 + (1.60000000000
     &  000d2*(eta0 + eta(i,j,k))**3*exp(-((eta0 + eta(i,j,k))**2/sigma**2)))/si
     &  gma**8 - (1.20000000000000d2*(eta0 + eta(i,j,k))*exp(-((eta0 + eta(i,j,k)
     &  )**2/sigma**2)))/sigma**6)
        dnnnnnngtil = amp*((6.4000000000000d1*(-eta0 + eta(i,j,k))**6*exp(-(
     &  (-eta0 + eta(i,j,k))**2/sigma**2)))/sigma**12 - (4.8000000000000d2*(
     &  -eta0 + eta(i,j,k))**4*exp(-((-eta0 + eta(i,j,k))**2/sigma**2)))/sigma**
     &  10 + (7.2000000000000d2*(-eta0 + eta(i,j,k))**2*exp(-((-eta0 + eta(i,j,k)
     &  )**2/sigma**2)))/sigma**8 - (1.20000000000000d2*exp(-((-eta0 + 
     &  eta(i,j,k))**2/sigma**2)))/sigma**6 + (6.4000000000000d1*(eta0 + eta(i,j,k)
     &  )**6*exp(-((eta0 + eta(i,j,k))**2/sigma**2)))/sigma**12 - (4.8000
     &  000000000d2*(eta0 + eta(i,j,k))**4*exp(-((eta0 + eta(i,j,k))**2/sigma**2
     &  )))/sigma**10 + (7.2000000000000d2*(eta0 + eta(i,j,k))**2*exp(-((eta
     &  0 + eta(i,j,k))**2/sigma**2)))/sigma**8 - (1.20000000000000d2*exp(-(
     &  (eta0 + eta(i,j,k))**2/sigma**2)))/sigma**6)
        dnnnnnnngtil = amp*((-1.28000000000000d2*(-eta0 + eta(i,j,k))**7*exp
     &  (-((-eta0 + eta(i,j,k))**2/sigma**2)))/sigma**14 + (1.34400000000000
     &  d3*(-eta0 + eta(i,j,k))**5*exp(-((-eta0 + eta(i,j,k))**2/sigma**2)))/sig
     &  ma**12 - (3.3600000000000d3*(-eta0 + eta(i,j,k))**3*exp(-((-eta0 + e
     &  ta(i,j,k))**2/sigma**2)))/sigma**10 + (1.68000000000000d3*(-eta0 + e
     &  ta(i,j,k))*exp(-((-eta0 + eta(i,j,k))**2/sigma**2)))/sigma**8 - (1.28000
     &  000000000d2*(eta0 + eta(i,j,k))**7*exp(-((eta0 + eta(i,j,k))**2/sigma**2
     &  )))/sigma**14 + (1.34400000000000d3*(eta0 + eta(i,j,k))**5*exp(-((et
     &  a0 + eta(i,j,k))**2/sigma**2)))/sigma**12 - (3.3600000000000d3*(eta0
     &   + eta(i,j,k))**3*exp(-((eta0 + eta(i,j,k))**2/sigma**2)))/sigma**10 + (
     &  1.68000000000000d3*(eta0 + eta(i,j,k))*exp(-((eta0 + eta(i,j,k))**2/sigm
     &  a**2)))/sigma**8)
