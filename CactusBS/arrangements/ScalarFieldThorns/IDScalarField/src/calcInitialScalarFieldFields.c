#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"
#include "math.h"

#define SQR(x) ((x)*(x))

/*
 * provide implementations for aliased function
 */

void idscalarfield_calcInitialScalarFields(CCTK_REAL x, CCTK_REAL y, CCTK_REAL z, CCTK_REAL *phi, CCTK_REAL *grad_phi, CCTK_REAL *pi, CCTK_REAL *VV)
{
  DECLARE_CCTK_PARAMETERS;
  const double PI = 3.1415926535897932384626433832795;
  const double ReY21CONST = -0.5*sqrt( 15 / (2*PI) );

  double ret_phi = 0.0, ret_grad_phi[3] = {0,0,0}, ret_pi = 0.0, ret_VV = 0.0, Psi0m2, rp, rm;

  // calculate renormalization factor 
  //rp = sqrt( (x-xp)*(x-xp) + y*y + z*z );
  //rm = sqrt( (x-xm)*(x-xm) + y*y + z*z );
  //Psi0m2 = SQR( 4.0 * rp * rm / ( 4.0*rp*rm + 2*mp*rm + 2*mm*rp ) );

  /* we compute phi(x) = Product(phi0[i]*exp(-0.5*(r_i-r0[0])^2/chi[i]^2), i=0..num_clouds-1)
   * where r_i = sqrt((x-x0[i])^2+(y-y0[i])^2+(z-z0[i])^2)
   * and its gradient wrt. (x,y,z)
   * pi is set to zero unconditionally
   */

  if ( CCTK_EQUALS(field_type, "Bubbles") ) {
    for (int i=0; i < num_clouds; i++)
    {
      double xx, yy, zz, chi2;
      double r2, rr, curr_phi;
      double ReY21, Tank;

      xx = x-centre_x[i];
      yy = y-centre_y[i];
      zz = z-centre_z[i];

      r2 = xx*xx+yy*yy+zz*zz;
      rr = sqrt(r2) - r0[i];
      chi2 = SQR( rr / sigma0[i] );

      if ( CCTK_EQUALS(scalar_profile, "Gaussian" ) )
      {
        curr_phi = phi0[i] * exp(-0.5*chi2);

        /*ret_phi += curr_phi*Psi0m2;
        ret_grad_phi[0] += -curr_phi*rr/SQR(sigma0[i])*xx/sqrt(r2)*Psi0m2 + xx*Psi0m3_over_r3*curr_phi;
        ret_grad_phi[1] += -curr_phi*rr/SQR(sigma0[i])*yy/sqrt(r2)*Psi0m2 + yy*Psi0m3_over_r3*curr_phi;
        ret_grad_phi[2] += -curr_phi*rr/SQR(sigma0[i])*zz/sqrt(r2)*Psi0m2 + zz*Psi0m3_over_r3*curr_phi;
        */
      
        if ( sinusoidal )
          ret_phi += curr_phi*cos(kk*rr+phase);
        else
          ret_phi += curr_phi;
          ret_grad_phi[0] += -curr_phi*rr/SQR(sigma0[i])*xx/sqrt(r2);
          ret_grad_phi[1] += -curr_phi*rr/SQR(sigma0[i])*yy/sqrt(r2);
          ret_grad_phi[2] += -curr_phi*rr/SQR(sigma0[i])*zz/sqrt(r2);
          ret_pi += 0.0;

        if ( rcutoff )
        {
           if ( fabs(rr) > rcutoff )
           {
               ret_phi = 0.0;
               ret_grad_phi[0] = 0.0;
               ret_grad_phi[1] = 0.0;
               ret_grad_phi[2] = 0.0;
           } 
        }
      }

      else if ( CCTK_EQUALS(scalar_profile, "Tanh" ) )
      {
        ret_phi += phi0[i] * tanh( rr/sigma0[i] );
        ret_grad_phi[0] += phi0[i]/sigma0[i]*(1-SQR(tanh(rr/sigma0[i])))*xx/sqrt(r2);
        ret_grad_phi[1] += phi0[i]/sigma0[i]*(1-SQR(tanh(rr/sigma0[i])))*yy/sqrt(r2);
        ret_grad_phi[2] += phi0[i]/sigma0[i]*(1-SQR(tanh(rr/sigma0[i])))*zz/sqrt(r2);
        ret_pi += 0.0;
      }

      else if ( CCTK_EQUALS(scalar_profile, "Tanh10" ) )
      {
        ret_phi += 0.5*phi0[i]*(1-tanh( rr/sigma0[i] ) );
        ret_grad_phi[0] += -0.5*phi0[i]*(1-SQR(tanh(rr/sigma0[i])))*xx/sqrt(r2)/sigma0[i];
        ret_grad_phi[1] += -0.5*phi0[i]*(1-SQR(tanh(rr/sigma0[i])))*yy/sqrt(r2)/sigma0[i];
        ret_grad_phi[2] += -0.5*phi0[i]*(1-SQR(tanh(rr/sigma0[i])))*zz/sqrt(r2)/sigma0[i];
        ret_pi += 0.0;
      }
          
      else if ( CCTK_EQUALS(scalar_profile, "Y21Tanh" ) )
      {
        ret_phi += phi0[i] * (1 - tanh( rr/sigma0[i] ));
        ret_grad_phi[0] += -phi0[i]/sigma0[i]*(1-SQR(tanh(rr/sigma0[i])))*xx/sqrt(r2);
        ret_grad_phi[1] += -phi0[i]/sigma0[i]*(1-SQR(tanh(rr/sigma0[i])))*yy/sqrt(r2);
        ret_grad_phi[2] += -phi0[i]/sigma0[i]*(1-SQR(tanh(rr/sigma0[i])))*zz/sqrt(r2);
        ret_pi += 0.0;
        if ( r2 != 0 )
        {
          ReY21 = ReY21CONST * xx*zz/r2;
          Tank  = phi0[i] * (1 - tanh( rr/sigma0[i] ));
          ret_phi += eps * Tank * ReY21;
          ret_grad_phi[0] += -eps * phi0[i]/sigma0[i]*(1-SQR(tanh(rr/sigma0[i])))*xx/sqrt(r2) * ReY21 + 
                              eps * ReY21CONST * Tank * (zz/r2 * ( 1 - 2*SQR(xx)/r2 ) );
          ret_grad_phi[1] += -eps * phi0[i]/sigma0[i]*(1-SQR(tanh(rr/sigma0[i])))*yy/sqrt(r2) * ReY21 +
                              eps * ReY21CONST * Tank * (-2 * xx * yy * zz/SQR(r2) );
          ret_grad_phi[2] += -eps * phi0[i]/sigma0[i]*(1-SQR(tanh(rr/sigma0[i])))*zz/sqrt(r2) * ReY21 +
                              eps * ReY21CONST * Tank * (xx/r2 * ( 1 - 2*SQR(zz)/r2 ) );
          ret_pi += 0.0;
        }
      }
      
      else if ( CCTK_EQUALS(scalar_profile, "perturbation_sin"))
      {
        double k;
        k = (2.*PI*N_sin)/domain_length; 

        ret_phi += phiBG * (1  + epsilon*(sin(k*xx) * sin(k*yy) * sin(k*zz)));
        ret_grad_phi[0] += phiBG * epsilon * k * cos(k * xx) * sin(k * yy) * sin(k * zz);
        ret_grad_phi[1] += phiBG * epsilon * k * sin(k * xx) * cos(k * yy) * sin(k * zz);
        ret_grad_phi[2] += phiBG * epsilon * k * sin(k * xx) * sin(k * yy) * cos(k * zz);
      }

      else if  (CCTK_EQUALS(scalar_profile, "perturbation_cos"))
      {
        double k;
        k = (2.*PI*N_sin)/domain_length; 

        ret_phi += phiBG * (1  + epsilon*(cos(k*xx) * cos(k*yy) * cos(k*zz)));
        ret_grad_phi[0] += -phiBG * epsilon * k * sin(k * xx) * cos(k * yy) * cos(k * zz);
        ret_grad_phi[1] += -phiBG * epsilon * k * cos(k * xx) * sin(k * yy) * cos(k * zz);
        ret_grad_phi[2] += -phiBG * epsilon * k * cos(k * xx) * cos(k * yy) * sin(k * zz);

      }

      else if (CCTK_EQUALS(scalar_profile, "fermi_sin")){
        double k;
        k = (2.*PI*N_sin)/domain_length; 
        r2 = xx*xx+yy*yy+zz*zz;

        ret_phi +=  phiBG * (1  + epsilon*(sin(k*xx) * sin(k*yy) * sin(k*zz)) * (pow(exp((sqrt(r2) - rcutoff)/alpha) + 1.0 , -1 ) ) );

        ret_grad_phi[0] += phiBG * epsilon * (k * cos(k * xx) * sin(k * yy) * sin(k * zz) * (pow(exp((sqrt(r2) - rcutoff)/alpha) + 1.0 , -1))
          + sin(k*xx) * sin(k*yy) * sin(k*zz) * (-1 * pow((exp((sqrt(r2) - rcutoff)/alpha) + 1 ),-2)) * (xx/alpha)*sqrt(r2));

        ret_grad_phi[1] += phiBG * epsilon * (k * sin(k * xx) * cos(k * yy) * sin(k * zz) * (pow(exp((sqrt(r2) - rcutoff)/alpha)+1.0,-1)) 
          + sin(k*xx) * sin(k*yy) * sin(k*zz) * (-1 * pow(exp((sqrt(r2) - rcutoff)/alpha) + 1 , -2)) * (yy/alpha)*sqrt(r2));

        ret_grad_phi[2] += phiBG * epsilon * (k * sin(k * xx) * sin(k * yy) * cos(k * zz) * (pow(exp((sqrt(r2) - rcutoff)/alpha)+1.0,-1))
          + sin(k*xx) * sin(k*yy) * sin(k*zz) * (-1 * pow(exp((sqrt(r2) - rcutoff)/alpha) + 1, -2)) * (zz/alpha)*sqrt(r2));

      }

       else if (CCTK_EQUALS(scalar_profile, "fermi_cos")){
        double k;
        k = (2.*PI*N_sin)/domain_length; 
        r2 = xx*xx+yy*yy+zz*zz;

        ret_phi +=  phiBG * (1  + epsilon*(cos(k*xx) * cos(k*yy) * cos(k*zz)) * (pow(exp((sqrt(r2) - rcutoff)/alpha) + 1.0 , -1 ) ) );

        ret_grad_phi[0] += phiBG * epsilon * (k * -sin(k * xx) * cos(k * yy) * cos(k * zz) * (pow(exp((sqrt(r2) - rcutoff)/alpha) + 1.0 , -1))
          + cos(k*xx) * cos(k*yy) * cos(k*zz) * (-1 * pow((exp((sqrt(r2) - rcutoff)/alpha) + 1 ),-2)) * (xx/alpha)*sqrt(r2));

        ret_grad_phi[1] += phiBG * epsilon * (k * cos(k * xx) * -sin(k * yy) * cos(k * zz) * (pow(exp((sqrt(r2) - rcutoff)/alpha)+1.0,-1)) 
          + cos(k*xx) * cos(k*yy) * cos(k*zz) * (-1 * pow(exp((sqrt(r2) - rcutoff)/alpha) + 1 , -2)) * (yy/alpha)*sqrt(r2));

        ret_grad_phi[2] += phiBG * epsilon * (k * cos(k * xx) * cos(k * yy) * -sin(k * zz) * (pow(exp((sqrt(r2) - rcutoff)/alpha)+1.0,-1))
          + cos(k*xx) * cos(k*yy) * cos(k*zz) * (-1 * pow(exp((sqrt(r2) - rcutoff)/alpha) + 1, -2)) * (zz/alpha)*sqrt(r2));

      }

      else if (CCTK_EQUALS(scalar_profile, "sin_decaying")){
        double k;
        k = (2.*PI*N_sin)/domain_length; 
        r2 = xx*xx+yy*yy+zz*zz;

        ret_phi += phiBG * (1 + epsilon* (sin(k*xx)*sin(k*yy)*sin(k*zz) * (exp(-(r2)/ (sigma_decay*sigma_decay) )) ) );

        ret_grad_phi[0] += phiBG * epsilon * k * cos(k*xx)*sin(k*yy)*sin(k*zz) * (exp(-(r2)/ (sigma_decay*sigma_decay) )) 
          + phiBG * epsilon * sin(k*xx)*sin(k*yy)*sin(k*zz) * exp(-r2 / (sigma_decay * sigma_decay)) * (-2*xx / (sigma_decay * sigma_decay) );

        ret_grad_phi[1] += phiBG * epsilon * k * sin(k*xx)*cos(k*yy)*sin(k*zz) * (exp(-(r2)/ (sigma_decay*sigma_decay) )) 
          + phiBG * epsilon * sin(k*xx)*sin(k*yy)*sin(k*zz) * exp(-r2 / (sigma_decay * sigma_decay)) * (-2*yy / (sigma_decay * sigma_decay) );

        ret_grad_phi[2] += phiBG * epsilon * k * sin(k*xx)*sin(k*yy)*cos(k*zz) * (exp(-(r2)/ (sigma_decay*sigma_decay) )) 
          + phiBG * epsilon * sin(k*xx)*sin(k*yy)*sin(k*zz) * exp(-r2 / (sigma_decay * sigma_decay)) * (-2*zz / (sigma_decay * sigma_decay) );

      }
    } 
  }

  else if (CCTK_EQUALS(field_type, "Everywhere") ) {
    double r2, curr_phi, rr;
    r2 = x*x+y*y+z*z;
    rr = sqrt(r2);

    if (CCTK_EQUALS(scalar_profile, "Exponential") )
    {
      curr_phi = alphaE*exp(-betaE*r2)*(cos(gammaE*r2));
      //ret_phi += curr_phi;

      ret_phi = curr_phi;

      ret_grad_phi[0] = -2*alphaE*x*exp(-betaE*r2)*(betaE*cos(gammaE*r2) + gammaE*sin(gammaE*r2));
      ret_grad_phi[1] = -2*alphaE*y*exp(-betaE*r2)*(betaE*cos(gammaE*r2) + gammaE*sin(gammaE*r2));
      ret_grad_phi[2] = -2*alphaE*z*exp(-betaE*r2)*(betaE*cos(gammaE*r2) + gammaE*sin(gammaE*r2)); 
    }
    
    else if (CCTK_EQUALS(scalar_profile, "OneOverR") )
    {
      curr_phi = alphaE/(gammaE*r2 + betaE);

      ret_phi = curr_phi;

      ret_grad_phi[0] = -2*alphaE*gammaE*x/powl(gammaE*r2 + betaE, 2);
      ret_grad_phi[1] = -2*alphaE*gammaE*y/powl(gammaE*r2 + betaE, 2);
      ret_grad_phi[2] = -2*alphaE*gammaE*z/powl(gammaE*r2 + betaE, 2); 
    }

  }


    // add background phi and pi
    ret_phi += phiConst;
    ret_pi  += piConst;
    /*ret_phi += phiBG*Psi0m2;*/

    // add potential
    if (potential)
    {
      if ( CCTK_EQUALS(potential_type,"phiFourth") )
        ret_VV = Vlambda/8.0 * SQR( ret_phi*ret_phi - phiBG*phiBG ); 
        /*ret_VV = Vlambda* SQR( ret_phi*ret_phi);*/   
      else if ( CCTK_EQUALS(potential_type,"phiSixth") )
        ret_VV = Vlambda/8.0 * SQR( ret_phi*ret_phi - phiBG*phiBG ) * ret_phi*ret_phi;
      else if ( CCTK_EQUALS(potential_type,"phiEighth") )
        ret_VV = Vlambda/8.0 * SQR( ret_phi*ret_phi - phiBG*phiBG ) * SQR( ret_phi*ret_phi - Veta*Veta);
      else if ( CCTK_EQUALS(potential_type,"phiBubble") )
        ret_VV = Vlambda/8.0 * ( SQR( ret_phi*ret_phi - phiBG*phiBG ) - 4.0*Veps*phiBG*SQR(phiBG)* (ret_phi - phiBG) );
      else if ( CCTK_EQUALS(potential_type,"phiExpand") )
        ret_VV = V0 + Vlambda * SQR(ret_phi) * ( 1.0/6.0*SQR(ret_phi*ret_phi) - 0.25*(1+SQR(Veps))*SQR(Veta)*SQR(ret_phi) + 0.5*SQR(Veps)*SQR(SQR(Veta)) );
      /*ret_VV = Vlambda/8.0 * SQR( ret_phi*ret_phi - phiBG*phiBG*Psi0m2*Psi0m2 );*/
    }

    // set return values
    *phi = ret_phi;
    grad_phi[0] = ret_grad_phi[0];
    grad_phi[1] = ret_grad_phi[1];
    grad_phi[2] = ret_grad_phi[2];
    *pi = ret_pi;
    *VV = ret_VV;
}
