        o1 = -xi
        o2 = cos(qgrd(j))
        o3 = xi*o2
        o4 = 1.0 + o1 + o3
        o5 = -eta0_me
        o6 = etagrd(i) + o5
        o7 = o6**2
        o8 = sigma_me**2
        o9 = o8**2
        o10 = 1/o9
        o11 = 1/o8
        o12 = -(o11*o7)
        o13 = exp(o12)
        o14 = eta0_me + etagrd(i)
        o15 = o14**2
        o16 = -(o11*o15)
        o17 = exp(o16)
        o18 = sin(qgrd(j))
        o19 = o18**n_me
        o20 = o2**2
        o21 = o13 + o17
        qfetaeta(i,j) = amp_me*o19*o4*(-2.0*o11*o13 - 2.0*o11*o17 + 4.0*o10*o15*o
     &  17 + 4.0*o10*o13*o7)
        qfqq(i,j) = -(amp_me*xi*n_me*o19*o2*o21) - amp_me*xi*(1.0 + n_me)*o19*o2*o21 -
     &   amp_me*n_me*o19*o21*o4 + amp_me*(-1.0 + n_me)*n_me*o18**(-2 + n)*o20*o21*o
     &  4
