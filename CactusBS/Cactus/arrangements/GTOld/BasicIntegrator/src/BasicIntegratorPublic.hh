#ifndef BASIC_INTEGRATOR_HDR
#define BASIC_INTEGRATOR_HDR

#include <functional>
#include <vector>
#include "cctk.h"

namespace BasicIntegrator
{
	using CartFunc = std::function<CCTK_REAL(CCTK_REAL, CCTK_REAL, CCTK_REAL)>;
	struct CoordSpec
	{
		CCTK_REAL min, max;
		CCTK_INT npoints;
		CoordSpec(CCTK_REAL mn, CCTK_REAL mx, CCTK_INT np) : min(mn), max(mx), npoints(np) {}
	};

	struct IntegralParams
	{
		CoordSpec csx, csy, csz;
		IntegralParams(CoordSpec c1, CoordSpec c2, CoordSpec c3) : csx(c1), csy(c2), csz(c3) {}
	};

	struct RealTuple
	{
		CCTK_REAL x, y, z;
		RealTuple(CCTK_REAL inx, CCTK_REAL iny, CCTK_REAL inz) : x(inx), y(iny), z(inz) {}
		RealTuple() : x(0), y(0), z(0) {}
	};

	CCTK_REAL integrate(IntegralParams pars, CartFunc ifunc);
	CCTK_REAL computeWeight(CoordSpec cs);
	CCTK_REAL getValueFromIndex(IntegralParams pars, CartFunc ifunc, CCTK_INT idx);
}

#endif
