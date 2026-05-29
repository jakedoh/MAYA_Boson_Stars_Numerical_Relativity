        o1 = a**2
        o2 = a*o1
        o3 = -2.00000000000000d0*mass*rBL
        o4 = rBL**2
        o5 = o1 + o3 + o4
        o6 = sqrt(o5)
        o7 = cos(phi(i,j,k))
        o8 = cos(q(i,j,k))
        o9 = o8**2
        o10 = o1*o9
        o11 = o10 + o4
        o12 = o11**2
        o13 = 1/o12
        o14 = -3.00000000000000d0*eta(i,j,k)
        o15 = exp(o14)
        o16 = sin(phi(i,j,k))
        o17 = sin(q(i,j,k))
        o18 = o17**2
        o19 = o1 + o4
        o20 = 2.00000000000000d0*o19*o4
        o21 = -o1
        o22 = o21 + o4
        o23 = o11*o22
        o24 = o20 + o23
        o25 = o7**2
        o26 = o16**2
        o27 = o17*o18
        kxx(i,j,k) = -2.00000000000000d0*a*mass*o13*o15*o16*o18*o24*o7 +
     &   4.0000000000000d0*mass*o13*o15*o16*o18*o2*o6*o7*o9*rBL
        kxy(i,j,k) = a*mass*o13*o15*o18*o24*o25 - a*mass*o13*o15*o18*o24
     &  *o26 - 2.00000000000000d0*mass*o13*o15*o18*o2*o25*o6*o9*rBL + 2.
     &  00000000000000d0*mass*o13*o15*o18*o2*o26*o6*o9*rBL
        kxz(i,j,k) = -(a*mass*o13*o15*o16*o17*o24*o8) - 2.00000000000000
     &  d0*mass*o13*o15*o16*o2*o27*o6*o8*rBL
        kyy(i,j,k) = 2.00000000000000d0*a*mass*o13*o15*o16*o18*o24*o7 - 
     &  4.0000000000000d0*mass*o13*o15*o16*o18*o2*o6*o7*o9*rBL
        kyz(i,j,k) = a*mass*o13*o15*o17*o24*o7*o8 + 2.00000000000000d0*m
     &  ass*o13*o15*o2*o27*o6*o7*o8*rBL
        kzz(i,j,k) = 0
