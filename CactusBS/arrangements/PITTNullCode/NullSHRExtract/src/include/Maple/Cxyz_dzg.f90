      t0 = -6*z*M*x**2/sqrt(x**2+y**2+z**2)**5
                              Cxyz_dzg(1,1)=t0

      t0 = -6*M*y*z/sqrt(x**2+y**2+z**2)**5*x
                              Cxyz_dzg(1,2)=t0

      t0 = 2*M*x*(x**2+y**2-2*z**2)/sqrt(x**2+y**2+z**2)**5
                              Cxyz_dzg(1,3)=t0

      t0 = 4*M*z/(x**2+y**2+z**2)**2*x
                              Cxyz_dzg(1,4)=t0

      t0 = -6*z*M*y**2/sqrt(x**2+y**2+z**2)**5
                              Cxyz_dzg(2,2)=t0

      t0 = 2*M*y*(x**2+y**2-2*z**2)/sqrt(x**2+y**2+z**2)**5
                              Cxyz_dzg(2,3)=t0

      t0 = 4*M*z/(x**2+y**2+z**2)**2*y
                              Cxyz_dzg(2,4)=t0

      t0 = 2*z*M*(2*x**2+2*y**2-z**2)/sqrt(x**2+y**2+z**2)**5
                              Cxyz_dzg(3,3)=t0

      t0 = -2*M*(x**2+y**2-z**2)/(x**2+y**2+z**2)**2
                              Cxyz_dzg(3,4)=t0

      t0 = -2*M*z/sqrt(x**2+y**2+z**2)**3
                              Cxyz_dzg(4,4)=t0

