	qf(i,j,k) = amp*(1.00000000000000d0+c*cos(phigrd(k))**2)*
     $	(exp(-((etagrd(i)-eta0)**2/sigma**2))+exp(-((etagrd(i)+eta0)**2
     $	/sigma**2)))*sin(qgrd(j))**n

	qfetaeta(i,j,k) = amp*(1.00000000000000d0+c*cos(phigrd(k))
     $	**2)*((4.0000000000000d0*(etagrd(i)-eta0)**2*exp(-((etagrd(i)-eta0)
     $	**2/sigma**2)))/sigma**4-(2.00000000000000d0*exp(-((etagrd(i)
     $	-eta0)**2/sigma**2)))/sigma**2+(4.0000000000000d0*(etagrd(i)+
     $	eta0)**2*exp(-((etagrd(i)+eta0)**2/sigma**2)))/sigma**4-(2.00
     $	000000000000d0*exp(-((etagrd(i)+eta0)**2/sigma**2)))/sigma**2
     $	)*sin(qgrd(j))**n

	qfqq(i,j,k) = amp*(-1.00000000000000d0+n)*n*(1.00000000
     $	000000d0+c*cos(phigrd(k))**2)*cos(qgrd(j))**2*(exp(-((etagrd(i)-eta0)**2/
     $	sigma**2))+exp(-((etagrd(i)+eta0)**2/sigma**2)))*sin(qgrd(j))**(-2
     $	+n)-amp*n*(1.00000000000000d0+c*cos(phigrd(k))
     $	**2)*(exp(-((etagrd(i)-eta0)**2/sigma**2))+exp(-((etagrd(i)+eta0)**2
     $	/sigma**2)))*sin(qgrd(j))**n

	qfphi(i,j,k) = -2.00000000000000d0*amp*c*cos(phigrd(k))*(exp(-
     $	((etagrd(i)-eta0)**2/sigma**2))+exp(-((etagrd(i)+eta0)**2/sigma**2))
     $	)*sin(phigrd(k))*sin(qgrd(j))**n

	qfphiphi(i,j,k) = -2.00000000000000d0*amp*c*cos(phigrd(k))**2*
     $	(exp(-((etagrd(i)-eta0)**2/sigma**2))+exp(-((etagrd(i)+eta0)**2/
     $	sigma**2)))*sin(qgrd(j))**n+2.00000000000000d0*amp*c*(exp(-((
     $	etagrd(i)-eta0)**2/sigma**2))+exp(-((etagrd(i)+eta0)**2/sigma**2)))*
     $	sin(phigrd(k))**2*sin(qgrd(j))**n

