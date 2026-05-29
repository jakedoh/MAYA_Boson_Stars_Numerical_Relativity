#include <numeric>
#include <algorithm>
#include "BasicIntegratorPublic.hh"

namespace BasicIntegrator
{
	CCTK_REAL integrate(IntegralParams pars, CartFunc ifunc)
	{
		const auto add = [](CCTK_REAL a, CCTK_REAL b) { return a+b; };


		const auto index_mapper = [pars, &ifunc](CCTK_INT in)
		{
			return getValueFromIndex(pars, ifunc, in);
		};

		const auto len = pars.csx.npoints*pars.csy.npoints*pars.csz.npoints;
		std::vector<CCTK_INT> start_int(len);
		std::vector<CCTK_REAL> start_real(len);

		std::iota(start_int.begin(), start_int.end(), 0);
		std::transform(start_int.begin(), start_int.end(), start_real.begin(), index_mapper);
		const auto unweighted_sum = std::accumulate(start_real.begin(), start_real.end(), 0., add);

		return unweighted_sum*computeWeight(pars.csx)*computeWeight(pars.csy)*computeWeight(pars.csz);
	}

	CCTK_REAL computeWeight(CoordSpec cs)
	{
		//doing npoints - 1 ensures that sampling is symmetric
		//npoints-1 then represents the number of intervals
		return (cs.max-cs.min)/(cs.npoints);
	}

	CCTK_REAL getValueFromIndex(IntegralParams pars, CartFunc ifunc, CCTK_INT idx)
	{
		const auto iz = idx/(pars.csx.npoints*pars.csy.npoints);
		const auto iy = idx/pars.csx.npoints-iz*pars.csy.npoints;
		const auto ix = idx-(iz*pars.csy.npoints+iy)*pars.csx.npoints;

		const auto rx = pars.csx.min+(ix+.5)*computeWeight(pars.csx);
		const auto ry = pars.csy.min+(iy+.5)*computeWeight(pars.csy);
		const auto rz = pars.csz.min+(iz+.5)*computeWeight(pars.csz);
		return ifunc(rx, ry, rz);
	}
}
