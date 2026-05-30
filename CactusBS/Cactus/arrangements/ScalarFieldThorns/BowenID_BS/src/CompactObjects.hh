#ifndef BOWEN_COMPACT_CLASSES
#define BOWEN_COMPACT_CLASSES

#include "cctk_Parameters.h"
#include <VectorAlgebra.hh>
#include "Utilities.hh"
#include <vector>
#include <cmath>
#define SQR(x) ((x)*(x))
#define CUBE(x) ((x)*(x)*(x))
#define QUAD(x) (SQR(x)*SQR(x))

struct tovdata_t
{
    double psi0, theta0, amp0,  psi1, theta1, amp1, rhoH, dtphi, riso, qq, cc, q2, q4;
    //psi0: zeroth derivative of conformal factor
    //theta0: zeroth derivative of psi*alpha
    //amp0: zeroth derivative of scalar field amplitude
    //psi1, theta1, amp1: first derivatives
    //riso: r isotropic (to avoid confusion with r_schwarzchild)
    //rhoH: rho_H 
    //dtphi: omega*amp0
    
    double mass;

    tovdata_t();
    tovdata_t(const tovdata_t& a);
    tovdata_t& operator=(const tovdata_t& a);
    tovdata_t& operator*=(double a);
    tovdata_t operator*(double a) const;
    tovdata_t& operator/=(double a);
    tovdata_t operator/(double a) const;
    tovdata_t& operator+=(const tovdata_t &a);
    tovdata_t operator+(const tovdata_t &a) const;
};

tovdata_t operator*(double a, const tovdata_t& b);

class Source
{
	public:
                virtual double getconfAlpha(double xx, double yy, double zz) = 0;
                virtual VectorAlgebra::Vector getInitialPhi(double xx, double yy, double zz) = 0;
		virtual VectorAlgebra::Vector getInitialConformalPi(double xx, double yy, double zz) = 0;
                virtual double getGradPhisqr(double xx, double yy, double zz) = 0;
		virtual double getNonFlatPsi(double xx, double yy, double zz) = 0;
		virtual double getInitialPsiGuess(double xx, double yy, double zz) = 0;
		virtual double getMass() { return 0.; };


		//Aij is *added* to the existing matrix
		virtual VectorAlgebra::LinearOperator getAij(double xx, double yy, double zz) = 0;

		virtual ~Source() = 0;
};

class CompactSource : public Source
{
	protected:
		VectorAlgebra::Vector position;
		int index;
	public:
		VectorAlgebra::Vector getRelativeVector(double xx, double yy, double zz) { VectorAlgebra::Vector input(xx, yy, zz); return input-position; };
		double getRelativeDistance(double xx, double yy, double zz) { return smoothR(getRelativeVector(xx, yy, zz).getMagnitude()); };

		CompactSource(int index);
		virtual ~CompactSource() = 0;
		
};

class EmptySource : public CompactSource
{
	public:
                double getconfAlpha(double xx, double yy, double zz) { return 0.0;} ;
                VectorAlgebra::Vector getInitialPhi(double xx, double yy, double zz) {return VectorAlgebra::zero; };
		VectorAlgebra::Vector getInitialConformalPi(double xx, double yy, double zz) {return VectorAlgebra::zero; };
                double getGradPhisqr(double xx, double yy, double zz) {return 0.0; };
		double getNonFlatPsi(double xx, double yy, double zz) { return 0.0;} ;
		double getInitialPsiGuess(double xx, double yy, double zz) { return 0.0; };

		double getAdmMass() { return 0.0; };

		//Aij is *added* to the existing matrix
		void addAij(double xx, double yy, double zz, double Aij[3][3]) {};
		~EmptySource() {};
		EmptySource(int index) : CompactSource(index) {};
};

class BowenSource : public CompactSource
{
	protected:
		VectorAlgebra::Vector momentum;
		VectorAlgebra::Vector spin;
		double mass;
		virtual double getMomentumDist(double xx, double yy, double zz) = 0;
		virtual double getSpinDist(double rr) = 0;

		virtual double getMomentumWeight(double rr) = 0; //Q
		virtual double getPointSourceDeviation(double rr) = 0; //C
		virtual double getSpinWeight(double rr) = 0; //N
		virtual double getNonFlatPsiOfR(double rr) = 0;
		
		double simpsonIntegrate(double (BowenSource::*dFdr)(double), double rmin, double rmid, double rmax, double fmin, double fmid, double fmax, int depth);
		double radialIntegrate(double (BowenSource::*dFdr)(double), double rmin, double rmax);
		double radius;

	public:
		virtual ~BowenSource() = 0;
		BowenSource(int idx);
		double getNonFlatPsi(double xx, double yy, double zz);
		VectorAlgebra::LinearOperator getAij(double xx, double yy, double zz);
		double getMass() { return mass; };
};


class BlackHole : public BowenSource
{
	public:
		BlackHole(int idx);
		double getInitialPsiGuess(double xx, double yy, double zz) { return getNonFlatPsi(xx, yy, zz); };
		~BlackHole() {};
		double getAdmMass();
                double getconfAlpha(double xx, double yy, double zz) { return 0.; };
                //double  getconfAlpha(double xx, double yy, double zz) { DECLARE_CCTK_PARAMETERS; double rr = std::sqrt(SQR(xx)+SQR(yy)+SQR(zz)); if(add_conformal_lapse_bh) { return 1. - pow(1. + mass/2/rr + 0.1*0.1/mass/mass/rr, -lapse_exponent); } else { return 0.; } }
                VectorAlgebra::Vector getInitialPhi(double xx, double yy, double zz) {return VectorAlgebra::zero; };
		VectorAlgebra::Vector getInitialConformalPi(double xx, double yy, double zz) {return VectorAlgebra::zero; };
                double getGradPhisqr(double xx, double yy, double zz) {return 0.0; };
	protected:
		double getMomentumDist(double xx, double yy, double zz) { return 0.; }; //sigma_p = 0 outside the puncture
		double getSpinDist(double rr) { return 0.; }; //sigma_j = 0 outside the puncture
		double getMomentumWeight(double rr) { return 1.; }; //Q
		double getPointSourceDeviation(double rr) { return 0.; }; //C
		double getSpinWeight(double rr) { return 1.; };  //N
		double getNonFlatPsiOfR(double rr) { return mass/2./rr; };
};

class SphericalHydroSource : public BowenSource
{
	public:
		SphericalHydroSource(int idx) : BowenSource(idx) {};
		virtual ~SphericalHydroSource() = 0;
                double getconfAlpha(double xx, double yy, double zz);
                VectorAlgebra::Vector getInitialPhi(double xx, double yy, double zz);
		VectorAlgebra::Vector getInitialConformalPi(double xx, double yy, double zz);
                double getGradPhisqr(double xx, double yy, double zz);
		double getInitialPsiGuess(double xx, double yy, double zz) { return 0.;} ;

	protected:
		double mbar;
                double C; 
                double psi_rescale_factor;
		//double ibar;
		virtual double getamp0(double rr) = 0;
                virtual double getamp1(double rr) = 0;
                virtual double gettheta0(double rr) = 0;
                virtual double gettheta1(double rr) = 0;
                virtual double getpsi1(double rr) = 0;
                virtual double getdtphi(double rr) = 0;
		double getMomentumDist(double xx, double yy, double zz); //sigma_p
		double getSpinDist(double rr) { return /* getMomentumDist(double xx, double yy, double zz)*mbar/ibar*/ 0.; }; //sigma_j
};

class TOVStar : public SphericalHydroSource
{
	public:
		TOVStar(int index);
		~TOVStar() { rk4data.clear(); };
		
	protected:
		int inside_points;
		double dr;
		std::vector<tovdata_t> rk4data;
		
		double getamp0(double rr) { return interpolateData(rr).amp0; };
                double getamp1(double rr) { return interpolateData(rr).amp1; };
		double getNonFlatPsiOfR(double rr) { return interpolateData(rr).psi0 - 1.; };
                double gettheta0(double rr) {return interpolateData(rr).theta0; };
                double gettheta1(double rr) {return interpolateData(rr).theta1; };
                double getpsi1(double rr) {return interpolateData(rr).psi1; };
                double getdtphi(double rr) { return interpolateData(rr).dtphi; };
		double getMomentumWeight(double rr);
		double getPointSourceDeviation(double rr);
		double getSpinWeight(double rr);
		tovdata_t interpolateData(double rr);
		tovdata_t getRHS(double rr, double w, tovdata_t xx);
		tovdata_t tov_rkstep(double rr, double dr, double w, tovdata_t in);

		
};

#endif
