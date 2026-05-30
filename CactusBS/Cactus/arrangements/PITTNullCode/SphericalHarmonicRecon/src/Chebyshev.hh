#define CHEBYSHEV_NMAX 100
namespace recon_Chebyshev
{

/*
   This code will be pretty innefficient. Especially if we only want < 10 modes. But we only call it once so ...  Note nmax is the maximum value of n you are requesting, so th arrays need to hold nmax+1 elements.
*/
void ChebyshevU(int nmax, double x, double U[], double dxU[]);
}
