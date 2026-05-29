#include "CompactObjects.hh"
#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include <VectorAlgebra.hh>
#define MAX_BOWEN_SOURCES 2

std::vector<Source*> sources;

extern "C"
void bowenid_getSourceInfo(double xx, double yy, double zz, double* jhat, double *psi, double *psiguess, double* Aij)
{
	*psi = 1.;
	*psiguess = 1.;
	VectorAlgebra::LinearOperator Aop(VectorAlgebra::zero, VectorAlgebra::zero, VectorAlgebra::zero);

	//CCTK_VInfo(CCTK_THORNSTRING, "Aop initial info at (%g, %g, %g): Axx = %g, Axy = %g, Axz = %g, Ayy = %g, Ayz = %g, Azz = %g", xx, yy, zz, Aop.x.x, Aop.x.y, Aop.z.x, Aop.y.y, Aop.z.y, Aop.z.z);

	for(std::vector<Source*>::iterator it = sources.begin(); it != sources.end(); ++it)
	{
		*psi += (*it)->getNonFlatPsi(xx, yy, zz);
		*psiguess += (*it)->getInitialPsiGuess(xx, yy, zz);
		Aop += (*it)->getAij(xx, yy, zz);

//		VectorAlgebra::LinearOperator dA = (*it)->getAij(xx, yy, zz);
//		CCTK_VInfo(CCTK_THORNSTRING, "dA loop info at (%g, %g, %g): Axx = %g, Axy = %g, Axz = %g, Ayy = %g, Ayz = %g, Azz = %g", xx, yy, zz, dA.x.x, dA.x.y, dA.z.x, dA.y.y, dA.z.y, dA.z.z);

//		Aop += dA;
	}

	Aij[0] = Aop.x.x;
	Aij[1] = Aop.x.y;
	Aij[2] = Aop.x.z;
	Aij[3] = Aop.y.y;
	Aij[4] = Aop.y.z;
	Aij[5] = Aop.z.z;
	//CCTK_VInfo(CCTK_THORNSTRING, "Aop final info at (%g, %g, %g): Axx = %g, Axy = %g, Axz = %g, Ayy = %g, Ayz = %g, Azz = %g", xx, yy, zz, Aop.x.x, Aop.x.y, Aop.z.x, Aop.y.y, Aop.z.y, Aop.z.z);

	/*
	Aij[0] = 0.;
	Aij[1] = 0.;
	Aij[2] = 0.;
	Aij[3] = 0.;
	Aij[4] = 0.;
	Aij[5] = 0.;
	*/
}


extern "C" 
CCTK_REAL bowenid_getExternalPsiGuess(CCTK_INT idx, CCTK_REAL x, CCTK_REAL y, CCTK_REAL z)
{
	CCTK_REAL res = 0.;
	for(auto it = sources.begin(); it != sources.end(); ++it)
	{
		res += (std::distance(sources.begin(), it) == idx) ? 0. : (*it)->getInitialPsiGuess(x, y, z);
	}
        //CCTK_VInfo(CCTK_THORNSTRING,"calculating external psi guess");
	return res;
}


extern "C"
CCTK_REAL bowenid_getAlpha(CCTK_REAL x, CCTK_REAL y, CCTK_REAL z)
{
        CCTK_REAL res = 1.;
        for(auto it = sources.begin(); it != sources.end(); ++it)
        {
                res -= (*it)->getconfAlpha(x, y, z);
        }
        //CCTK_VInfo(CCTK_THORNSTRING,"calculating conformal lapse");
        return res;
}
             

extern "C"
CCTK_REAL bowenid_getGradPhisqr(CCTK_REAL x, CCTK_REAL y, CCTK_REAL z)
{       
        CCTK_REAL res = 0.0;
        for(auto it = sources.begin(); it != sources.end(); ++it)
        {       
                res += (*it)->getGradPhisqr(x, y, z);
        }
        return res;
}
        

extern "C"
void bowenid_getInitialPhi(double xx, double yy, double zz, double* Phi)
{
        VectorAlgebra::Vector Phivec = VectorAlgebra::zero;
        for(std::vector<Source*>::iterator it = sources.begin(); it != sources.end(); ++it)
        {
           Phivec += (*it)->getInitialPhi(xx, yy, zz);
        }

        Phi[0] = Phivec.x;
        Phi[1] = Phivec.y;
        Phi[2] = Phivec.z;
}

extern "C"
void bowenid_getInitialConformalPi(double xx, double yy, double zz, double* confPi)
{
        VectorAlgebra::Vector confPivec = VectorAlgebra::zero;
        for(std::vector<Source*>::iterator it = sources.begin(); it != sources.end(); ++it)
        {
           confPivec += (*it)->getInitialConformalPi(xx, yy, zz);
        }

        confPi[0] = confPivec.x;
        confPi[1] = confPivec.y;
        confPi[2] = confPivec.z;
}

extern "C"
double bowenid_getSourceMass(int index)
{
	return sources[index]->getMass();
        //CCTK_VInfo(CCTK_THORNSTRING,"calculating source mass");
}

//extern "C"
//double bowenid_getRestMass(int index)
//{
//	return sources[index]->getRestMass();
//        CCTK_VInfo(CCTK_THORNSTRING,"calculating rest mass");
//}

extern "C"
void bowenid_constructSources(CCTK_ARGUMENTS)
{
	bowenid_reconstructSources(0);
        CCTK_VInfo(CCTK_THORNSTRING,"initialized");
}

extern "C"
void bowenid_reconstructSources(CCTK_INT resolve)
{
	DECLARE_CCTK_PARAMETERS;
	static int need_solve = 1;

	if(resolve == 0 && need_solve == 0)
		return;

	CCTK_VInfo(CCTK_THORNSTRING, "Deleting all previous constructed sources");
	for(std::vector<Source*>::iterator it = sources.begin(); it != sources.end(); ++it)
	{
		delete *it;
	}
	sources.clear();
	for(int i = 0; i < MAX_BOWEN_SOURCES; i++)
	{
		if(CCTK_Equals(compact_object[i], "black hole"))
		{
			CCTK_VInfo(CCTK_THORNSTRING, "Creating puncture black hole for compact object #%d", i);
			sources.push_back(new BlackHole(i));
		}
		else if(CCTK_Equals(compact_object[i], "boson star"))
		{
			CCTK_VInfo(CCTK_THORNSTRING, "Creating boson star for compact object #%d", i);
			TOVStar *star = new TOVStar(i);
			sources.push_back(star);
//			bowen_tov_mass[i] = star->getMW();
		}
			
	}
	need_solve = 0;
}
