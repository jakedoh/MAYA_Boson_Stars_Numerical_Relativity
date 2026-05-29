        gtil = amp*(exp(-((etagrd(i) - eta0)**2/sigma**2)) + exp(-((etagrd(i) + eta0
     &  )**2/sigma**2)))
        dngtil = amp*((-2.00000000000000d0*(etagrd(i) - eta0)*exp(-((etagrd(i) - eta
     &  0)**2/sigma**2)))/sigma**2 - (2.00000000000000d0*(etagrd(i) + eta0)*ex
     &  p(-((etagrd(i) + eta0)**2/sigma**2)))/sigma**2)
