        o1 = dq**2
        o2 = 1/o1
        o3 = 1/dq
        o4 = tan(qgrd(j))
        o5 = 1/o4
        o6 = deta**2
        o7 = 1/o6
        o8 = sin(qgrd(j))
        o9 = o8**2
        o10 = 1/o9
        o11 = cos(qgrd(j))
        o12 = o11**2
        o13 = etagrd(i)**2
        o14 = sigma**2
        o15 = 1/o14
        o16 = -(o13*o15)
        o17 = -2.00000000000000d0*etagrd(i)*eta0*o15
        o18 = eta0**2
        o19 = -(o15*o18)
        o20 = o16 + o17 + o19
        o21 = exp(o20)
        i22 = -2 + n
c       The next lines causes a floating point exception on Alphas,
c       if o22 is 0.0 but o8 is negative.
c       This is because of y**x = exp(x*ln(y)).
c       Changing the type of the exponent to integer fixes the problem.
c        o22 = -2.00000000000000d0 + n
c        o23 = o8**o22
        o23 = o8**i22
        o24 = n**2
        o25 = 2.00000000000000d0*etagrd(i)*eta0*o15
        o26 = o16 + o19 + o25
        o27 = exp(o26)
        o28 = o8**n
        o29 = o14**2
        o30 = 1/o29
        o31 = o4**2
        o32 = 1/o31
        o33 = 5.0000000000000d-1*etagrd(i)
        o34 = cosh(o33)
        Cn(i,j) = o2 + 5.0000000000000d-1*o3*o5
        Cs(i,j) = o2 - 5.0000000000000d-1*o3*o5
        Cw(i,j) = o7
        Cc(i,j) = -1.25000000000000d-1 - 1.25000000000000d-1*o10 - 2.000
     &  00000000000d0*o2 - 2.50000000000000d-1*amp*n*o12*o21*o23 + 2.500
     &  00000000000d-1*amp*o12*o21*o23*o24 - 2.50000000000000d-1*amp*n*o
     &  12*o23*o27 + 2.50000000000000d-1*amp*o12*o23*o24*o27 - 2.5000000
     &  0000000d-1*amp*n*o21*o28 - 5.0000000000000d-1*amp*o15*o21*o28 - 
     &  2.50000000000000d-1*amp*n*o27*o28 - 5.0000000000000d-1*amp*o15*o
     &  27*o28 + 2.00000000000000d0*amp*etagrd(i)*eta0*o21*o28*o30 + amp*o13*o
     &  21*o28*o30 + amp*o18*o21*o28*o30 - 2.00000000000000d0*amp*etagrd(i)*et
     &  a0*o27*o28*o30 + amp*o13*o27*o28*o30 + amp*o18*o27*o28*o30 + 1.2
     &  5000000000000d-1*o32 - 2.00000000000000d0*o7
        Ce(i,j) = o7
        Rhs(i,j) = -2.50000000000000d-1*o34 + 2.50000000000000d-1*o10*o3
     &  4 + 5.0000000000000d-1*amp*n*o12*o21*o23*o34 - 5.0000000000000d-
     &  1*amp*o12*o21*o23*o24*o34 + 5.0000000000000d-1*amp*n*o12*o23*o27
     &  *o34 - 5.0000000000000d-1*amp*o12*o23*o24*o27*o34 + 5.0000000000
     &  000d-1*amp*n*o21*o28*o34 + amp*o15*o21*o28*o34 + 5.0000000000000
     &  d-1*amp*n*o27*o28*o34 + amp*o15*o27*o28*o34 - 4.0000000000000d0*
     &  amp*etagrd(i)*eta0*o21*o28*o30*o34 - 2.00000000000000d0*amp*o13*o21*o2
     &  8*o30*o34 - 2.00000000000000d0*amp*o18*o21*o28*o30*o34 + 4.00000
     &  00000000d0*amp*etagrd(i)*eta0*o27*o28*o30*o34 - 2.00000000000000d0*amp
     &  *o13*o27*o28*o30*o34 - 2.00000000000000d0*amp*o18*o27*o28*o30*o3
     &  4 - 2.50000000000000d-1*o32*o34

