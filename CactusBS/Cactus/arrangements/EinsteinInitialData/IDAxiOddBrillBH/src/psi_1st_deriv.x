        o1 = 5.0000000000000d-1*eta(i,j,k)
        o2 = exp(o1)
        o3 = psi2d(i,j,k)
        o4 = 1/o3
        o5 = cos(phi(i,j,k))
        o6 = cos(q(i,j,k))
        o7 = dqpsi2d(i,j,k)
        o8 = -1.50000000000000d0*eta(i,j,k)
        o9 = exp(o8)
        o10 = detapsi2d(i,j,k)
        o11 = sin(q(i,j,k))
        o12 = sin(phi(i,j,k))
        psix(i,j,k) = o2*o4*(o10*o11*o5*o9 - 5.0000000000000d-1*o11*o3*o
     &  5*o9 + o5*o6*o7*o9)
        psiy(i,j,k) = o2*o4*(o10*o11*o12*o9 - 5.0000000000000d-1*o11*o12
     &  *o3*o9 + o12*o6*o7*o9)
        psiz(i,j,k) = o2*o4*(o10*o6*o9 - 5.0000000000000d-1*o3*o6*o9 - o
     &  11*o7*o9)
