        o1 = 5.0000000000000d-1*eta(i,j,k)
        o2 = exp(o1)
        o3 = psi3d(i,j,k)
        o4 = 1/o3
        o5 = cos(phi(i,j,k))
        o6 = cos(q(i,j,k))
        o7 = dqpsi3d(i,j,k)
        o8 = -1.50000000000000d0*eta(i,j,k)
        o9 = exp(o8)
        o10 = dphipsi3d(i,j,k)
        o11 = sin(phi(i,j,k))
        o12 = sin(q(i,j,k))
        o13 = 1/o12
        o14 = detapsi3d(i,j,k)
        psix(i,j,k) = o2*o4*(-(o10*o11*o13*o9) + o12*o14*o5*o9 - 5.00000
     &  00000000d-1*o12*o3*o5*o9 + o5*o6*o7*o9)
        psiy(i,j,k) = o2*o4*(o11*o12*o14*o9 - 5.0000000000000d-1*o11*o12
     &  *o3*o9 + 1.00000000000000d0*o10*o13*o5*o9 + o11*o6*o7*o9)
        psiz(i,j,k) = o2*o4*(o14*o6*o9 - 5.0000000000000d-1*o3*o6*o9 - o
     &  12*o7*o9)
