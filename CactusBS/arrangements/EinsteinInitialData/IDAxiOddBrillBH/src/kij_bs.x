        o1 = cos(phi(i,j,k))
        o2 = cos(q(i,j,k))
        o3 = o2**2
        o4 = -3.00000000000000d0*eta(i,j,k)
        o5 = exp(o4)
        o6 = sin(phi(i,j,k))
        o7 = sin(q(i,j,k))
        o8 = -1.00000000000000d0 + n
        o9 = o7**o8
        o10 = 2.00000000000000d0 + n
        o11 = o7**2
        o12 = -(o10*o11)
        o13 = 1.00000000000000d0 + n + o12
        o14 = o1**2
        o15 = o6**2
        o16 = o7**n
        o17 = -2.00000000000000d0 + n
        o18 = o7**o17
        kxx(i,j,k) = -2.00000000000000d0*gtil*o1*o13*o5*o6*o9 + 2.000000
     &  00000000d0*dngtil*o1*o3*o5*o6*o9
        kxy(i,j,k) = gtil*o13*o14*o5*o9 - gtil*o13*o15*o5*o9 - dngtil*o1
     &  4*o3*o5*o9 + dngtil*o15*o3*o5*o9
        kxz(i,j,k) = -(dngtil*o16*o2*o5*o6) - gtil*o13*o18*o2*o5*o6
        kyy(i,j,k) = 2.00000000000000d0*gtil*o1*o13*o5*o6*o9 - 2.0000000
     &  0000000d0*dngtil*o1*o3*o5*o6*o9
        kyz(i,j,k) = dngtil*o1*o16*o2*o5 + gtil*o1*o13*o18*o2*o5
        kzz(i,j,k) = 0
