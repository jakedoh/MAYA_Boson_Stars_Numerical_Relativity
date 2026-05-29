        mass = sqrt(2.0 + sqrt(4.0 + byJ**2
     &  ))
        a = byJ/mass
        rBL = 5.0e-1*exp(eta(i,j,k))*(1.0 + (1.0*(a
     &   + mass)*exp(-eta(i,j,k)))/sqrt(-a**2 + mass**2))*sqrt(-a**2 + mass**2)
     &  *(1.0 + (1.*(a + mass)*exp(-eta(i,j,k)))/sqrt(a**2 + m
     &  ass**2))
