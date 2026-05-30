#include "cctk_Parameters.h"
#include "CompactObjects.hh"
#include <sstream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <OutputTable.hh>
#include "cctk.h"
#define M_TO_KM  1.47695518286//1.48
#define MSUN_DENS_TO_CGI 6.17327129839e17  //6.1694e17 // Convert rest-mass density to CGI for M=Msun

//physical frequency of scalar field
double omega_rescaled;

//Pure Virtual Destructors; required
Source::~Source() {}
CompactSource::~CompactSource() {}
BowenSource::~BowenSource() {}
SphericalHydroSource::~SphericalHydroSource() {}

//Position of compact object
CompactSource::CompactSource(int idx)
{
	DECLARE_CCTK_PARAMETERS;
	index = idx;
	position.x = object_rx[index];
	position.y = object_ry[index];
	position.z = object_rz[index];
}

//Nontrivial portion of conformal factor
double BowenSource::getNonFlatPsi(double xx, double yy, double zz)
{
	//double rr = getRelativeDistance(xx, yy, zz);
	VectorAlgebra::Vector rr = getRelativeVector(xx, yy, zz);
        double rmag = smoothR(rr.getMagnitude());
	double psi = getNonFlatPsiOfR(rmag);
//	CCTK_VInfo(CCTK_THORNSTRING, "psi information at (%g, %g, %g) (r = %g): psi = %g", xx, yy, zz, rr, psi);
	return psi;
}

//Momenta and spin of compact object
BowenSource::BowenSource(int idx) : CompactSource(idx)
{
	DECLARE_CCTK_PARAMETERS;
	momentum.x = object_Px[index];
	momentum.y = object_Py[index];
	momentum.z = object_Pz[index];
	spin.x = object_Jx[index];
	spin.y = object_Jy[index];
	spin.z = object_Jz[index];
}

//A_ij
VectorAlgebra::LinearOperator BowenSource::getAij(double xx, double yy, double zz)
{
	DECLARE_CCTK_PARAMETERS;
	VectorAlgebra::Vector rr = getRelativeVector(xx, yy, zz);
	double rmag = smoothR(rr.getMagnitude());
        
	// l^i
	VectorAlgebra::Vector direction = rr/rmag;
        
	// eps^ijk l_i S_k
	VectorAlgebra::Vector lxs = VectorAlgebra::crossProduct(direction, spin);
        
        // P^i l^j + P^j l^i
	VectorAlgebra::LinearOperator psym(momentum.x*direction + direction.x*momentum, momentum.y*direction + direction.y*momentum, momentum.z*direction + direction.z*momentum);

        // l^i l^j
	VectorAlgebra::LinearOperator project(direction.x*direction, direction.y*direction, direction.z*direction);

	VectorAlgebra::LinearOperator lxssym(lxs.x*direction + direction.x*lxs, lxs.y*direction+direction.y*lxs, lxs.z*direction+direction.z*lxs);

	// P^i l_i
	double ldotp = VectorAlgebra::dotProduct(direction, momentum);
	
        //Equation
	double Q = getMomentumWeight(rmag);

        //Equation 
	double C = getPointSourceDeviation(rmag);

        //Equation
	double nn = getSpinWeight(rmag);
        
        //Equation
	VectorAlgebra::LinearOperator Aij_mom = 1.5*Q/SQR(rmag)*(psym+ ldotp*(project-VectorAlgebra::identity)) +3*C/QUAD(rmag)*(psym + ldotp*(VectorAlgebra::identity - 5.*project));

        //Equation
	VectorAlgebra::LinearOperator Aij_spin = - 3.*nn/CUBE(rmag)*(lxssym);

	if(momentum.getMagnitude() < 1e-10)
	{
		Aij_mom = 0.*VectorAlgebra::identity;
	}
	if(spin.getMagnitude() < 1e-10)
	{
		Aij_spin = 0.*VectorAlgebra::identity;
	}

	VectorAlgebra::LinearOperator Aij = Aij_mom + Aij_spin;
	
        /*
	CCTK_VInfo(CCTK_THORNSTRING, "Aij info at (%g, %g, %g) (r = %g): q = %g, c = %g, n = %g, Axx = %g, Axy = %g, Axz = %g, Ayy = %g, Ayz = %g, Azz = %g", xx, yy, zz, rmag, qq, cc, nn, Aij.x.x, Aij.x.y, Aij.z.x, Aij.y.y, Aij.z.y, Aij.z.z);
	CCTK_VInfo(CCTK_THORNSTRING, "n = (%g, %g, %g), P = (%g, %g, %g), J = (%g, %g, %g)", direction.x, direction.y, direction.z, momentum.x, momentum.y, momentum.z, spin.x, spin.y, spin.z);
	CCTK_VInfo(CCTK_THORNSTRING, "n cross J = (%g, %g, %g), (n cross J) otimes n: (xx, xy, xz, yx, yy, yz, zx, zy, zz) = (%g, %g, %g, %g, %g, %g, %g, %g, %g)", lxs.x, lxs.y, lxs.z, lxssym.x.x, lxssym.y.x, lxssym.z.x, lxssym.x.y, lxssym.y.y, lxssym.z.y, lxssym.x.z, lxssym.y.z, lxssym.z.z);
	CCTK_VInfo(CCTK_THORNSTRING, "p otimes n: (xx, xy, xz, yx, yy, yz, zx, zy, zz) = (%g, %g, %g, %g, %g, %g, %g, %g, %g)", psym.x.x, psym.y.x, psym.z.x, psym.x.y, psym.y.y, psym.z.y, psym.x.z, psym.y.z, psym.z.z);
	CCTK_VInfo(CCTK_THORNSTRING, "projector: (xx, xy, xz, yx, yy, yz, zx, zy, zz) = (%g, %g, %g, %g, %g, %g, %g, %g, %g)", project.x.x, project.y.x, project.z.x, project.x.y, project.y.y, project.z.y, project.x.z, project.y.z, project.z.z);
	*/

	return Aij;
}

BlackHole::BlackHole(int idx) : BowenSource(idx)
{
	DECLARE_CCTK_PARAMETERS;
	mass = bh_bare_mass[index];
	radius = 0.;
}

double BowenSource::radialIntegrate(double (BowenSource::*dFdr)(double), double rmin, double rmax)
{
	double fmin = (this->*dFdr)(rmin);
	double fmax = (this->*dFdr)(rmax);

	double rmid = (rmin+rmax)/2.0;
	double fmid = (this->*dFdr)(rmid);

	if(std::isnan(fmin) || std::isnan(fmid) || std::isnan(fmax))
	{
		CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING, "radial integration failed from nans: f min, mid, max = (%g, %g, %g) at (%g, %g, %g)", fmin, fmid, fmax, rmin, rmid, rmax);
	}
	return simpsonIntegrate(dFdr, rmin, rmid, rmax, fmin, fmid, fmax, 0);
}

double BowenSource::simpsonIntegrate(double (BowenSource::*dFdr)(double), double rmin, double rmid, double rmax, double fmin, double fmid, double fmax, int depth)
{
	DECLARE_CCTK_PARAMETERS;
	double S1 = (fmin+4.*fmid+fmax)/6.*(rmax-rmin);

	double rq1 = (rmin+rmid)/2.0;
	double fq1 = (this->*dFdr)(rq1);
	double S21 = (fmin+4*fq1+fmid)/6.*(rmid-rmin);


	double rq3 = (rmax+rmid)/2.0;
	double fq3 = (this->*dFdr)(rq3);
	double S22 = (fmid+4*fq3+fmax)/6.*(rmax-rmid);

	double S2 = S21 + S22;

	double delta = (S2-S1)/15.;
	if(fabs(delta) < max_quadrature_error*pow(.5,depth) || depth > 3)
	{
		return S2 + delta;
	}
	else
	{
		return simpsonIntegrate(dFdr, rmin, rq1, rmid, fmin, fq1, fmid, depth+1) + simpsonIntegrate(dFdr, rmid, rq3, rmax, fmid, fq3, fmax, depth+1);
	}
}

/*
double BowenSource::radialIntegrate(double (BowenSource::*dFdr)(double), double rmin, double rmax)
{
	DECLARE_CCTK_PARAMETERS;
	//adaptive simpson's rule (tentative; can implement other quadratures and use this function as a wrapper)	
	//TODO: a more precise integrator (e.g. Gauss-Kronrod?)
	double ffmin = (this->*dFdr)(rmin);
	double ffmax = (this->*dFdr)(rmax);

	double rmid = (rmin+rmax)/2.0;
	double ffmid = (this->*dFdr)(rmid);
	
	double rqtr1 = (rmin+rmid)/2.0;
	double ffqtr1 = (this->*dFdr)(rqtr1);

	double rqtr3 = (rmax+rmid)/2.0;
	double ffqtr3 = (this->*dFdr)(rqtr3);

	double S1 = simpson_weighted_integral(ffmin, ffmid, ffmax)*(rmax - rmin);
	double S2 = simpson_weighted_integral(ffmin, ffqtr1, ffmid)*(rmid-rmin) + simpson_weighted_integral(ffmid, ffqtr3, ffmax)*(rmax-rmid);

	bool too_small = ((rmax-rmin) < minimum_dr);
	if(too_small)
		CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING, "Simpison integration reached minimum size; results may be invalid!");

	double delta = (S2-S1)/15.;
	if(fabs(delta) < max_quadrature_error || too_small)
	{
		return S2 - delta;
	}
	else
	{
		return radialIntegrate(dFdr, rmin, rmid) + radialIntegrate(dFdr, rmid, rmax);
	}	
}
*/

//Calculating (1-conformal lapse), where the conformal lapse is set to (lapse of stationary star)/(psi of stationary star)**6
double SphericalHydroSource::getconfAlpha(double xx, double yy, double zz)
{
        DECLARE_CCTK_PARAMETERS;
	double theta0, psi, alpha, confalpha, oneminusconfalpha;
        
	//rr= getRelativeDistance(xx,yy,zz);
	VectorAlgebra::Vector rr = getRelativeVector(xx, yy, zz);
        double rmag = smoothR(rr.getMagnitude());
        
	theta0 = gettheta0(rmag);
        psi = 1.+getNonFlatPsiOfR(rmag);
        alpha = theta0/psi;
        
	confalpha = alpha/CUBE(psi)/CUBE(psi);	
	oneminusconfalpha = 1. - confalpha;

        //we use one minus confalpha to generalize to binary case	
	return oneminusconfalpha;
}

//Get initial complex scalar field at t=0
VectorAlgebra::Vector SphericalHydroSource::getInitialPhi(double xx, double yy, double zz)
{
        DECLARE_CCTK_PARAMETERS;
        double rmag, varphi, kdotx;

        VectorAlgebra::Vector rr = getRelativeVector(xx, yy, zz);
        rmag = smoothR(rr.getMagnitude());

        varphi = getamp0(rmag);

        VectorAlgebra::Vector k = momentum/mbar;
	kdotx = VectorAlgebra::dotProduct(k, rr);

        //For convenience, we will store the initial scalar field in a vector
        VectorAlgebra::Vector Phi;        
	Phi.x = varphi*std::cos(-kdotx);  //Re(Phi)
	Phi.y = varphi*std::sin(-kdotx);  //Im(Phi)
	Phi.z = 0.;

        return Phi;
}

//Calculate conformally rescaled Pi = -n^a \partial_a Phi
VectorAlgebra::Vector SphericalHydroSource::getInitialConformalPi(double xx, double yy, double zz)
{
        DECLARE_CCTK_PARAMETERS;
        double rmag, Phi_Re, Phi_Im, Theta, psi, alpha, confalpha;

        VectorAlgebra::Vector rr = getRelativeVector(xx, yy, zz);
        rmag = smoothR(rr.getMagnitude());

        VectorAlgebra::Vector Phi = getInitialPhi(xx, yy, zz);
        Phi_Re = Phi.x;
        Phi_Im = Phi.y;

        Theta = gettheta0(rmag);
        psi = 1.+getNonFlatPsiOfR(rmag);
        alpha = Theta/psi;
        confalpha = alpha/CUBE(psi)/CUBE(psi);

        VectorAlgebra::Vector confPi;
        confPi.x = omega_rescaled/confalpha*Phi_Im; //Re(Pi)
        confPi.y = -omega_rescaled/confalpha*Phi_Re; //Im(Pi)
        confPi.z = 0.;

        return confPi;
}

//Calculate gradient of scalar field amplitude
double SphericalHydroSource::getGradPhisqr(double xx, double yy, double zz)
{
        DECLARE_CCTK_PARAMETERS;
	double rmag, varphi, varphiprime, gradvarphisqr, ksqr, GradPhisqr;
        
	VectorAlgebra::Vector rr = getRelativeVector(xx, yy, zz);
        rmag = smoothR(rr.getMagnitude());
        VectorAlgebra::Vector direction = rr/rmag;
        VectorAlgebra::Vector k = momentum/mbar;
        
        varphi = getamp0(rmag);
	varphiprime = getamp1(rmag);
        
	VectorAlgebra::Vector gradvarphi = direction*varphiprime;	
        
        gradvarphisqr = VectorAlgebra::dotProduct(gradvarphi, gradvarphi);
        ksqr = VectorAlgebra::dotProduct(k, k);
 
        GradPhisqr = gradvarphisqr + ksqr*SQR(varphi);

	return GradPhisqr;
}


//Perform numerical integration of Equations (70)-(72) 
TOVStar::TOVStar(int idx) : SphericalHydroSource(idx)
{
	DECLARE_CCTK_PARAMETERS;
	dr = tov_dr[index];

 	std::stringstream ofile;
	ofile << out_dir << "/tov_internaldata_" << index << ".asc";

	OutputTable otable(ofile, 15, std::ios_base::out);
	if(testsuite)
	{
		otable << "#1:riso" << "2:amp0" << "3:amp1" << "4:psi" << "5:theta" << "6:rhoH" << "7:dtphi";
		otable << std::endl;
	}

	//Central density of boson star and associated eigenfrequency
	double A_c = BS_central_amplitude[index];
        double w = omega[index];

        //Asymptotic values of scalar field amplitude obtained via shooting method
	double amp0_inf = BS_asymptotic_amplitude[index];     

        //Initial conditions at r=0 for radial integration
	tovdata_t temp;
        temp.psi0 = 1.;  // will be adjusted later on
        temp.psi1 = 0.;                           
        temp.theta0 = 1.; // will be adjusted later on
        temp.theta1 = 0.;      
        temp.amp0 = A_c;
        temp.amp1 = 0.;                     
        temp.riso = 0.;
        
        //These are not dynamical variables, they will be calculated after the integration
        temp.rhoH = 0.;
        temp.dtphi = 0.;                           

	rk4data.push_back(temp);
           int i = 0;
              do
              {
	         if(i == tov_numpoints)
	            CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING, "Integration of hydrodynamic equations failed for TOV initial data: domain too small, still inside star after %d points, with amplitude %f\n. Initial data: A=%f\n, A'=%f\n, theta=%f\n, theta'=%f\n, Phi=%f\n, Phi'=%f\n, omega=%f\n", tov_numpoints, rk4data.back().amp0, rk4data[0].amp0, rk4data[0].amp1, rk4data[0].theta0, rk4data[0].theta1, rk4data[0].psi0, rk4data[0].psi1, w);

	         rk4data.push_back(tov_rkstep(i*dr, dr, w, rk4data[i]));
              } while(rk4data[++i].amp0 > amp0_inf);

        inside_points = i;
        int last_point = i-1;
        
	double rad_unscaled, m, psi_rescale_factor, theta_rescale_factor;
        
	rad_unscaled =rk4data[last_point].riso;
 
        psi_rescale_factor = ((rk4data[last_point]).psi0+rk4data[last_point].psi1*rad_unscaled);
        theta_rescale_factor = ((rk4data[last_point]).theta0+rk4data[last_point].theta1*rad_unscaled);
 
	omega_rescaled = w/theta_rescale_factor*psi_rescale_factor;
        
        //Equation
	mbar = rk4data[last_point].qq;
	
	otable << std::scientific;
        for(int l = 0; l < inside_points; l++)
           {           
              //Calculating rhoH and dtphi at each point
              double V = pow(rk4data[l].amp0, 2.0);
              rk4data[l].rhoH = 0.5 * ( pow(rk4data[l].psi0*rk4data[l].amp0*w/rk4data[l].theta0,2.0) + pow(rk4data[l].amp1/rk4data[l].psi0/rk4data[l].psi0,2.0) + V  );
              rk4data[l].dtphi = omega_rescaled*rk4data[l].amp0;


              rk4data[l].psi0 /= psi_rescale_factor;
              rk4data[l].psi1 /= CUBE(psi_rescale_factor);
              rk4data[l].riso *= SQR(psi_rescale_factor);
              rk4data[l].theta0 /= theta_rescale_factor;
              rk4data[l].theta1 /= (SQR(psi_rescale_factor)*theta_rescale_factor);
              rk4data[l].amp1 /= SQR(psi_rescale_factor); 
	      rk4data[l].qq /= mbar;
	      rk4data[l].cc *= QUAD(psi_rescale_factor)/mbar;
           }
    
        //asymptotic constants for scalar field 
        double k = std::sqrt(1. - SQR(omega_rescaled/rk4data[last_point].theta0*rk4data[last_point].psi0));
        C = rk4data[last_point].amp0*rk4data[last_point].riso*std::exp(k*rk4data[last_point].riso);
          
        for(int l = 0; l < inside_points; l++)
           {
           if(testsuite)
              {
                 otable << rk4data[l].riso << rk4data[l].amp0 << rk4data[l].amp1 << rk4data[l].psi0 << rk4data[l].theta0 << rk4data[l].rhoH << rk4data[l].dtphi;
                 otable << std::endl;
              } 
           }  

        dr *= SQR(psi_rescale_factor);
    
        radius = rk4data[last_point].riso;
        mass = 2.*radius*(rk4data[last_point].psi0-1.);

        if(verbose >= 1)
        {
           CCTK_INFO("Information on solved Boson star:");
		CCTK_VInfo(CCTK_THORNSTRING, "Boson star #%d info: central_amplitude = %g, omega = %g, mbar = %g", index, A_c, omega_rescaled, mbar);
        }
}

tovdata_t TOVStar::getRHS(double riso, double f, tovdata_t x)
{
        DECLARE_CCTK_PARAMETERS;
        double mass_scalar = mu;	
        tovdata_t v;
       
        double arg0 = x.psi0/x.theta0*f*x.amp0;
        double arg1 = x.amp1/x.psi0/x.psi0;
        double arg2 = x.psi0/x.theta0*f;

        double V = pow(mass_scalar*x.amp0, 2.0);
        double dVdphi = pow(mass_scalar, 2.0);
	double rhoH = 0.5 * ( pow(arg0,2.0) + pow(arg1,2.0) + V );
	double S = -0.5 * ( pow(arg1,2.0) + 3.0*(-pow(arg0,2.0) + V)  );

	v.psi1 = -2.0*M_PI*pow(x.psi0,5.0)*rhoH;
	v.theta1 = 2.0*M_PI*x.theta0*pow(x.psi0,4.0)*(rhoH+2.0*S);
        v.amp1 = -(x.theta1/x.theta0 + x.psi1/x.psi0)*x.amp1 + (dVdphi-pow(arg2,2.0))*pow(x.psi0,4.0)*x.amp0;
	v.psi0 = x.psi1;
	v.theta0 = x.theta1;
        v.amp0 = x.amp1;
	v.riso = 1.0;

        //will be computed after
        v.rhoH = 0.0;
        v.dtphi = 0.0;
       
        double confalpha = x.theta0/QUAD(x.psi0)/CUBE(x.psi0);

        double sig1 = f*SQR(x.amp0)/confalpha;
	v.qq = 4.0*M_PI*sig1*SQR(riso);
        v.cc = 2.0*M_PI/3.0*sig1*QUAD(riso);

        //used for sources of the form 3 (P^j l_j) l^i sig2(r)
        double sig2 = 0.;
        v.q2 = 4.0*M_PI*sig2*SQR(riso);
        v.q4 = 4.0*M_PI*sig2*QUAD(riso);	

	if(riso == 0.0)
	{
	   v.psi1 += 4.0*M_PI*rhoH/3.0;
	   v.theta1 += -4.0*M_PI*(rhoH+2.0*S)/3.0;
           v.amp1 += -2.0/3.0*( -(x.theta1/x.theta0 + x.psi1/x.psi0)*x.amp1 + (dVdphi-pow(arg2,2.0))*pow(x.psi0,4.0)*x.amp0 );
           //there probably ought to be some conditions here for v.qq, v.cc; see conditions for v.psi1, v.theta1
	}
	else
	{
           v.theta1 += -2.0*x.theta1/riso;
           v.psi1 += -2.0*x.psi1/riso;
           v.amp1 += -2.0*x.amp1/riso;
	}
	return v;
}

tovdata_t TOVStar::tov_rkstep(double x, double h, double ff, tovdata_t y_in)
{
        tovdata_t k1, k2, k3, k4, y_out;
        tovdata_t tmp;

        // k1
        k1 = getRHS(x, ff, y_in)*h;

        // k2
        tmp = y_in + k1/2.;
        k2 = getRHS(x + h/2, ff, tmp)*h;

        // k3
        tmp = y_in + k2/2.;
        k3 = getRHS(x + h/2, ff, tmp)*h;

        // k4
        tmp = y_in + k3;
        k4 = getRHS(x + h, ff, tmp)*h;

        // combine them
        y_out = y_in + 1/6. * (k1 + 2*k2 + 2*k3 + k4);
        return y_out;
}

tovdata_t TOVStar::interpolateData(double distance)
{
        DECLARE_CCTK_PARAMETERS;

        tovdata_t interp;
    
        int last_point = inside_points - 1;
        int i_low1, i_high1, i_low2, i_high2;
        i_low1 = (int) (distance/dr);
        i_high1 = i_low1 + 1;
	i_low2 = i_low1-1;
	i_high2 = i_high1+1;

	if(i_low2 < 0)
	{
	   i_low2 = 0;
	   i_low1 = 1;
	   i_high1 = 2;
	   i_high2 = 3;
	}
	if(i_high2 == inside_points)
	{
	   i_high2 = inside_points-1;
	   i_high1 = i_high2-1;
	   i_low1 = i_high1-1;
	   i_low2 = i_low1-1;	
	}
        if(i_high2 <= last_point)
        {
	   // box desired value such that riso[i_min] <= riso < riso[i_min+1]
	   //
	   /*
	   double driso, delta;

	   driso = distance - rk4data[i_min].riso; // distance from desired value
	   delta = rk4data[i_max].riso - rk4data[i_min].riso;  // step in riso 
	   double minweight, midweight, maxweight;
	   maxweight = (driso/delta);
	   minweight = 1. - maxweight;
	   interp = minweight*rk4data[i_min]+ maxweight*rk4data[i_max];
	   */
	   double c_low2, c_low1, c_high1, c_high2;
	   double x_low2, x_low1, x_high1, x_high2;

	   x_low2 = rk4data[i_low2].riso;
	   x_low1 = rk4data[i_low1].riso;
	   x_high1 = rk4data[i_high1].riso;
	   x_high2 = rk4data[i_high2].riso;

	   c_low2 = (distance-x_low1)*(distance-x_high1)*(distance-x_high2)/((x_low2-x_low1)*(x_low2-x_high1)*(x_low2-x_high2));
	   c_low1 = (distance-x_low2)*(distance-x_high1)*(distance-x_high2)/((x_low1-x_low2)*(x_low1-x_high1)*(x_low1-x_high2));
	   c_high1 = (distance-x_low2)*(distance-x_low1)*(distance-x_high2)/((x_high1-x_low2)*(x_high1-x_low1)*(x_high1-x_high2));
	   c_high2 = (distance-x_low2)*(distance-x_high1)*(distance-x_low1)/((x_high2-x_low2)*(x_high2-x_high1)*(x_high2-x_low1));
		
	   interp = c_low2*rk4data[i_low2]+c_low1*rk4data[i_low1]+c_high1*rk4data[i_high1]+c_high2*rk4data[i_high2];
	}
        else
	{
           double mass_scalar, w, k, V_inf;
	   
           interp = rk4data[last_point];          
               
           interp.psi0 = 1. + (interp.psi0 - 1.)*(radius/distance);
           interp.psi1 = (1. - interp.psi0)/distance;
           interp.theta0 = 1. - (1. - interp.theta0)*(radius/distance);
           interp.theta1 = (1. - interp.theta0)/distance;

           k = std::sqrt(1. - SQR(omega_rescaled/interp.theta0*interp.psi0));
           interp.amp0 = C*std::exp(-k*distance)/distance; //amp0_inf;
           interp.amp1 = -interp.amp0*(1/distance+k); //amp1_inf;
           interp.dtphi = w*interp.amp0; //w*amp0_inf;

           V_inf = pow(mass_scalar*interp.amp0, 2.);
           interp.rhoH = 0.5 * ( pow(interp.psi0*interp.amp0*w/interp.theta0,2.0) + pow(interp.amp1/interp.psi0/interp.psi0,2.0) + V_inf  );
           interp.riso = distance;
	}

	if(std::isnan(interp.amp0) || std::isnan(interp.amp1) || std::isnan(interp.rhoH) || std::isnan(interp.psi0) || std::isnan(interp.theta0) || std::isnan(interp.riso) || std::isnan(interp.rhoH))
	{
	   CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING, "interpolation of tov data returned nan: A = %g, A'=%g, rhoH = %g, psi = %g, theta = %g, riso = %g", interp.amp0, interp.amp1, interp.rhoH, interp.psi0, interp.theta0, interp.riso);
	}
        return interp;
}


double TOVStar::getMomentumWeight(double rr)
{
        DECLARE_CCTK_PARAMETERS;
        double Q = interpolateData(rr).qq + interpolateData(rr).q2;
    
	return Q;
}


double TOVStar::getPointSourceDeviation(double rr)
{
	DECLARE_CCTK_PARAMETERS;
        double C = interpolateData(rr).cc - 0.5*interpolateData(rr).q4;

        return C;
}

double TOVStar::getSpinWeight(double rr)
{
	/*
	double realrad = rr > radius ? radius : rr;
	return radialIntegrate(&BowenSource::getSpinWeightDeriv, 0., realrad);
	*/
	//return interpolateData(rr).nn;
        return 0.;
}


//NOTICE: this is consistent with the RHS routine, but not used there
double SphericalHydroSource::getMomentumDist(double xx, double yy, double zz)
{
	DECLARE_CCTK_PARAMETERS;
        double rr= getRelativeDistance(xx,yy,zz);
        double theta0 = gettheta0(rr);
        double psi = 1.+getNonFlatPsiOfR(rr);
        double alpha = theta0/psi;
	double sigma = pow(psi,6.0)/alpha*std::pow(getamp1(rr), 2.0)/mbar;
        return sigma;
} //"sigma_p"


tovdata_t::tovdata_t() {psi0=theta0=amp0=psi1=theta1=amp1=riso=rhoH=dtphi=qq=cc=q2=q4=0.0;}

tovdata_t::tovdata_t(const tovdata_t& a) {amp0 = a.amp0; amp1 = a.amp1;psi0 = a.psi0;theta0 = a.theta0;psi1 = a.psi1;theta1 = a.theta1; riso = a.riso; rhoH = a.rhoH; dtphi = a.dtphi; qq = a.qq; cc = a.cc; q2=a.q2; q4=a.q4; }

tovdata_t& tovdata_t::operator=(const tovdata_t& a) { amp0 = a.amp0; amp1 = a.amp1; psi0 = a.psi0; theta0 = a.theta0; psi1 = a.psi1; theta1 = a.theta1; riso = a.riso; rhoH = a.rhoH; dtphi = a.dtphi; qq = a.qq; cc = a.cc; q2=a.q2; q4=a.q4; return *this; }

tovdata_t& tovdata_t::operator*=(double a) {amp0 *= a; amp1 *= a; psi0 *= a; psi1 *= a; theta0 *= a;theta1*=a; riso *= a; rhoH *= a; dtphi *= a; qq *= a; cc *= a; q2 *= a; q4 *= a; return *this;}

tovdata_t tovdata_t::operator*(double a) const {return tovdata_t(*this).operator*=(a);}

tovdata_t& tovdata_t::operator/=(double a) {amp0 /= a;amp1 /= a;psi0 /= a;theta0 /= a;psi1 /= a; theta1 /= a; riso /= a; rhoH /= a; dtphi /= a; qq /= a; cc /= a; q2 /= a; q4 /= a; return *this;}

tovdata_t tovdata_t::operator/(double a) const {return tovdata_t(*this).operator/=(a);}

tovdata_t& tovdata_t::operator+=(const tovdata_t &a) {amp0 += a.amp0;amp1 += a.amp1;psi0 += a.psi0;theta0 += a.theta0;psi1 += a.psi1; theta1 += a.theta1; riso += a.riso; rhoH += a.rhoH; dtphi += a.dtphi; qq += a.qq; cc += a.cc; q2 += a.q2; q4 += a.q4; return *this;}

tovdata_t tovdata_t::operator+(const tovdata_t &a) const {return tovdata_t(*this).operator+=(a);}

tovdata_t operator*(double a, const tovdata_t& b) {return tovdata_t(b).operator*=(a); }
