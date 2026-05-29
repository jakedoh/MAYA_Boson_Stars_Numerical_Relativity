#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "math.h"

#define SQR(x) ((x)*(x))

/*
 * provide implementations for aliased function
 */

void idscalarfield_calcInitialScalarFields(CCTK_REAL x, CCTK_REAL y, CCTK_REAL z, CCTK_REAL *phia, CCTK_REAL *grad_phia, CCTK_REAL *pia, CCTK_REAL *phib, CCTK_REAL *grad_phib, CCTK_REAL *pib, CCTK_REAL *VV)
{
  DECLARE_CCTK_PARAMETERS;

  const double PI = 3.1415926535897932384626433832795;
  const double ReY21CONST = -0.5*sqrt( 15 / (2*PI) );

  double ret_phia = 0.0, ret_grad_phia[3] = {0,0,0}, ret_pia = 0.0, ret_phib = 0.0, ret_grad_phib[3] = {0,0,0}, ret_pib = 0.0, ret_VV = 0.0;

  // calculate renormalization factor 
  //rp, rm , Psi0m2
  //rp = sqrt( (x-xp)*(x-xp) + y*y + z*z );
  //rm = sqrt( (x-xm)*(x-xm) + y*y + z*z );
  //Psi0m2 = SQR( 4.0 * rp * rm / ( 4.0*rp*rm + 2*mp*rm + 2*mm*rp ) );

  /* we compute phi(x) = Product(phi0[i]*exp(-0.5*(r_i-r0[0])^2/chi[i]^2), i=0..num_clouds-1)
   * where r_i = sqrt((x-x0[i])^2+(y-y0[i])^2+(z-z0[i])^2)
   * and its gradient wrt. (x,y,z)
   * pi is set to zero unconditionally
   */

  if ( CCTK_EQUALS(field_type, "Bubbles") ){
    for (int i=0; i < num_clouds; i++)
    {
      double xxa, yya, zza, xxb, yyb, zzb, chia2, chib2;
      double ra2, rra, curr_phia, rb2, rrb, curr_phib;
      double ReY21a, Tanka, ReY21b, Tankb;

      xxa = x-centre_xa[i];
      yya = y-centre_ya[i];
      zza = z-centre_za[i];

      xxb = x-centre_xb[i];
      yyb = y-centre_yb[i];
      zzb = z-centre_zb[i];

      ra2 = xxa*xxa+yya*yya+zza*zza;
      rra = sqrt(ra2) - r0a[i];
      chia2 = SQR( rra / sigma0a[i] );


      rb2 = xxb*xxb+yyb*yyb+zzb*zzb;
      rrb = sqrt(rb2) - r0b[i];
      chib2 = SQR( rrb / sigma0b[i] );
      if ( CCTK_EQUALS(scalar_profile, "constant" ) )
     {  
        curr_phia = phi0a[i];
        curr_phib = phi0b[i];
          ret_phia += curr_phia;
          ret_grad_phia[0] += 0;
          ret_grad_phia[1] += 0;
          ret_grad_phia[2] += 0;
          ret_pia += 0.0;

          ret_phib += curr_phib;
          ret_grad_phib[0] += 0;
          ret_grad_phib[1] += 0;
          ret_grad_phib[2] += 0;
          ret_pib += 0.0;
}
      if ( CCTK_EQUALS(scalar_profile, "Gaussian" ) )
      {
        curr_phia = phi0a[i] * exp(-0.5*chia2);
        curr_phib = phi0b[i] * exp(-0.5*chib2);
        /*ret_phi += curr_phi*Psi0m2;
        ret_grad_phi[0] += -curr_phi*rr/SQR(sigma0[i])*xx/sqrt(r2)*Psi0m2 + xx*Psi0m3_over_r3*curr_phi;
        ret_grad_phi[1] += -curr_phi*rr/SQR(sigma0[i])*yy/sqrt(r2)*Psi0m2 + yy*Psi0m3_over_r3*curr_phi;
        ret_grad_phi[2] += -curr_phi*rr/SQR(sigma0[i])*zz/sqrt(r2)*Psi0m2 + zz*Psi0m3_over_r3*curr_phi;
        */
        if ( sinusoidal )
         {
          ret_phia += curr_phia*cos(kk*rra+phase);
          ret_phib += curr_phib*cos(kk*rrb+phase);
          }
        else   
         {   
          ret_phia += curr_phia;
          ret_grad_phia[0] += -curr_phia*rra/SQR(sigma0a[i])*xxa/sqrt(ra2);
          ret_grad_phia[1] += -curr_phia*rra/SQR(sigma0a[i])*yya/sqrt(ra2);
          ret_grad_phia[2] += -curr_phia*rra/SQR(sigma0a[i])*zza/sqrt(ra2);
          ret_pia += 0.0;

          ret_phib += curr_phib;
          ret_grad_phib[0] += -curr_phib*rrb/SQR(sigma0b[i])*xxb/sqrt(rb2);
          ret_grad_phib[1] += -curr_phib*rrb/SQR(sigma0b[i])*yyb/sqrt(rb2);
          ret_grad_phib[2] += -curr_phib*rrb/SQR(sigma0b[i])*zzb/sqrt(rb2);
          ret_pib += 0.0;
          }


        if ( rcutoff )
        {
           if ( fabs(rra) > rcutoff || fabs(rrb) > rcutoff)
           {
               ret_phia = 0.0;
               ret_grad_phia[0] = 0.0;
               ret_grad_phia[1] = 0.0;
               ret_grad_phia[2] = 0.0;
               ret_phib = 0.0;
               ret_grad_phib[0] = 0.0;
               ret_grad_phib[1] = 0.0;
               ret_grad_phib[2] = 0.0;

           } 
        }
      }   
      else if ( CCTK_EQUALS(scalar_profile, "Tanh" ) )
      {
        ret_phia += phi0a[i] * tanh( rra/sigma0a[i] );
        ret_grad_phia[0] += phi0a[i]/sigma0a[i]*(1-SQR(tanh(rra/sigma0a[i])))*xxa/sqrt(ra2);
        ret_grad_phia[1] += phi0a[i]/sigma0a[i]*(1-SQR(tanh(rra/sigma0a[i])))*yya/sqrt(ra2);
        ret_grad_phia[2] += phi0a[i]/sigma0a[i]*(1-SQR(tanh(rra/sigma0a[i])))*zza/sqrt(ra2);
        ret_pia += 0.0;
        ret_phib += phi0b[i] * tanh( rrb/sigma0b[i] );
        ret_grad_phib[0] += phi0b[i]/sigma0b[i]*(1-SQR(tanh(rrb/sigma0b[i])))*xxb/sqrt(rb2);
        ret_grad_phib[1] += phi0b[i]/sigma0b[i]*(1-SQR(tanh(rrb/sigma0b[i])))*yyb/sqrt(rb2);
        ret_grad_phib[2] += phi0b[i]/sigma0b[i]*(1-SQR(tanh(rrb/sigma0b[i])))*zzb/sqrt(rb2);
        ret_pib += 0.0;
      }
else if ( CCTK_EQUALS(scalar_profile, "Tanh10" ) )
      {
        ret_phia += 0.5*phi0a[i]*(1-tanh( rra/sigma0a[i] ) );
        ret_grad_phia[0] += -0.5*phi0a[i]*(1-SQR(tanh(rra/sigma0a[i])))*xxa/sqrt(ra2)/sigma0a[i];
        ret_grad_phia[1] += -0.5*phi0a[i]*(1-SQR(tanh(rra/sigma0a[i])))*yya/sqrt(ra2)/sigma0a[i];
        ret_grad_phia[2] += -0.5*phi0a[i]*(1-SQR(tanh(rra/sigma0a[i])))*zza/sqrt(ra2)/sigma0a[i];
        ret_pia += 0.0;
        ret_phib += 0.5*phi0b[i]*(1-tanh( rrb/sigma0b[i] ) );
        ret_grad_phib[0] += -0.5*phi0b[i]*(1-SQR(tanh(rrb/sigma0b[i])))*xxb/sqrt(rb2)/sigma0b[i];
        ret_grad_phib[1] += -0.5*phi0b[i]*(1-SQR(tanh(rrb/sigma0b[i])))*yyb/sqrt(rb2)/sigma0b[i];
        ret_grad_phib[2] += -0.5*phi0b[i]*(1-SQR(tanh(rrb/sigma0b[i])))*zzb/sqrt(rb2)/sigma0b[i];
        ret_pib += 0.0;
      }
          
      else if ( CCTK_EQUALS(scalar_profile, "Y21Tanh" ) )
      {
        ret_phia += phi0a[i] * (1 - tanh( rra/sigma0a[i] ));
        ret_grad_phia[0] += -phi0a[i]/sigma0a[i]*(1-SQR(tanh(rra/sigma0a[i])))*xxa/sqrt(ra2);
        ret_grad_phia[1] += -phi0a[i]/sigma0a[i]*(1-SQR(tanh(rra/sigma0a[i])))*yya/sqrt(ra2);
        ret_grad_phia[2] += -phi0a[i]/sigma0a[i]*(1-SQR(tanh(rra/sigma0a[i])))*zza/sqrt(ra2);
        ret_pia += 0.0;
        ret_phib += phi0b[i] * (1 - tanh( rrb/sigma0b[i] ));
        ret_grad_phib[0] += -phi0b[i]/sigma0b[i]*(1-SQR(tanh(rrb/sigma0b[i])))*xxb/sqrt(rb2);
        ret_grad_phib[1] += -phi0b[i]/sigma0b[i]*(1-SQR(tanh(rrb/sigma0b[i])))*yyb/sqrt(rb2);
        ret_grad_phib[2] += -phi0b[i]/sigma0b[i]*(1-SQR(tanh(rrb/sigma0b[i])))*zzb/sqrt(rb2);
        ret_pib += 0.0;
        if ( ra2 != 0 || rb2 != 0)
        {
          ReY21a = ReY21CONST * xxa*zza/ra2;
          Tanka  = phi0a[i] * (1 - tanh( rra/sigma0a[i] ));
          ret_phia += eps * Tanka * ReY21a;
          ret_grad_phia[0] += -eps * phi0a[i]/sigma0a[i]*(1-SQR(tanh(rra/sigma0a[i])))*xxa/sqrt(ra2) * ReY21a + 
                              eps * ReY21CONST * Tanka * (zza/ra2 * ( 1 - 2*SQR(xxa)/ra2 ) );
          ret_grad_phia[1] += -eps * phi0a[i]/sigma0a[i]*(1-SQR(tanh(rra/sigma0a[i])))*yya/sqrt(ra2) * ReY21a +
                              eps * ReY21CONST * Tanka * (-2 * xxa * yya * zza/SQR(ra2) );
          ret_grad_phia[2] += -eps * phi0a[i]/sigma0a[i]*(1-SQR(tanh(rra/sigma0a[i])))*zza/sqrt(ra2) * ReY21a +
                              eps * ReY21CONST * Tanka * (xxa/ra2 * ( 1 - 2*SQR(zza)/ra2 ) );
          ret_pia += 0.0;

          ReY21b = ReY21CONST * xxb*zzb/rb2;
          Tankb  = phi0b[i] * (1 - tanh( rrb/sigma0b[i] ));
          ret_phib += eps * Tankb * ReY21b;
          ret_grad_phib[0] += -eps * phi0b[i]/sigma0b[i]*(1-SQR(tanh(rrb/sigma0b[i])))*xxb/sqrt(rb2) * ReY21b + 
                              eps * ReY21CONST * Tankb * (zzb/rb2 * ( 1 - 2*SQR(xxb)/rb2 ) );
          ret_grad_phib[1] += -eps * phi0b[i]/sigma0b[i]*(1-SQR(tanh(rrb/sigma0b[i])))*yyb/sqrt(rb2) * ReY21b +
                              eps * ReY21CONST * Tankb * (-2 * xxb * yyb * zzb/SQR(rb2) );
          ret_grad_phib[2] += -eps * phi0b[i]/sigma0b[i]*(1-SQR(tanh(rrb/sigma0b[i])))*zzb/sqrt(rb2) * ReY21b +
                              eps * ReY21CONST * Tankb * (xxb/rb2 * ( 1 - 2*SQR(zzb)/rb2 ) );
          ret_pib += 0.0;
        }
      }
      
      else if ( CCTK_EQUALS(scalar_profile, "perturbation_sin"))
      {
        double k;
        k = (2.*PI*N_sin)/domain_length; 

        ret_phia += phiBG * (1  + epsilon*(sin(k*xxa) * sin(k*yya) * sin(k*zza)));
        ret_grad_phia[0] += phiBG * epsilon * k * cos(k * xxa) * sin(k * yya) * sin(k * zza);
        ret_grad_phia[1] += phiBG * epsilon * k * sin(k * xxa) * cos(k * yya) * sin(k * zza);
        ret_grad_phia[2] += phiBG * epsilon * k * sin(k * xxa) * sin(k * yya) * cos(k * zza);

        ret_phib += phiBG * (1  + epsilon*(sin(k*xxb) * sin(k*yyb) * sin(k*zzb)));
        ret_grad_phib[0] += phiBG * epsilon * k * cos(k * xxb) * sin(k * yyb) * sin(k * zzb);
        ret_grad_phib[1] += phiBG * epsilon * k * sin(k * xxb) * cos(k * yyb) * sin(k * zzb);
        ret_grad_phib[2] += phiBG * epsilon * k * sin(k * xxb) * sin(k * yyb) * cos(k * zzb);
      }

      else if  (CCTK_EQUALS(scalar_profile, "perturbation_cos"))
      {
        double k;
        k = (2.*PI*N_sin)/domain_length; 

        ret_phia += phiBG * (1  + epsilon*(cos(k*xxa) * cos(k*yya) * cos(k*zza)));
        ret_grad_phia[0] += -phiBG * epsilon * k * sin(k * xxa) * cos(k * yya) * cos(k * zza);
        ret_grad_phia[1] += -phiBG * epsilon * k * cos(k * xxa) * sin(k * yya) * cos(k * zza);
        ret_grad_phia[2] += -phiBG * epsilon * k * cos(k * xxa) * cos(k * yya) * sin(k * zza);

        ret_phib += phiBG * (1  + epsilon*(cos(k*xxb) * cos(k*yyb) * cos(k*zzb)));
        ret_grad_phib[0] += -phiBG * epsilon * k * sin(k * xxb) * cos(k * yyb) * cos(k * zzb);
        ret_grad_phib[1] += -phiBG * epsilon * k * cos(k * xxb) * sin(k * yyb) * cos(k * zzb);
        ret_grad_phib[2] += -phiBG * epsilon * k * cos(k * xxb) * cos(k * yyb) * sin(k * zzb);
      }

      else if (CCTK_EQUALS(scalar_profile, "fermi_sin")){
        double k;
        k = (2.*PI*N_sin)/domain_length; 
        /*r2 = xx*xx+yy*yy+zz*zz;*/

        ret_phia +=  phiBG * (1  + epsilon*(sin(k*xxa) * sin(k*yya) * sin(k*zza)) * (pow(exp((sqrt(ra2) - rcutoff)/alpha) + 1.0 , -1 ) ) );

        ret_grad_phia[0] += phiBG * epsilon * (k * cos(k * xxa) * sin(k * yya) * sin(k * zza) * (pow(exp((sqrt(ra2) - rcutoff)/alpha) + 1.0 , -1))
          + sin(k*xxa) * sin(k*yya) * sin(k*zza) * (-1 * pow((exp((sqrt(ra2) - rcutoff)/alpha) + 1 ),-2)) * (xxa/alpha)*sqrt(ra2));

        ret_grad_phia[1] += phiBG * epsilon * (k * sin(k * xxa) * cos(k * yya) * sin(k * zza) * (pow(exp((sqrt(ra2) - rcutoff)/alpha)+1.0,-1)) 
          + sin(k*xxa) * sin(k*yya) * sin(k*zza) * (-1 * pow(exp((sqrt(ra2) - rcutoff)/alpha) + 1 , -2)) * (yya/alpha)*sqrt(ra2));

        ret_grad_phia[2] += phiBG * epsilon * (k * sin(k * xxa) * sin(k * yya) * cos(k * zza) * (pow(exp((sqrt(ra2) - rcutoff)/alpha)+1.0,-1))
          + sin(k*xxa) * sin(k*yya) * sin(k*zza) * (-1 * pow(exp((sqrt(ra2) - rcutoff)/alpha) + 1, -2)) * (zza/alpha)*sqrt(ra2));

        ret_phib +=  phiBG * (1  + epsilon*(sin(k*xxb) * sin(k*yyb) * sin(k*zzb)) * (pow(exp((sqrt(rb2) - rcutoff)/alpha) + 1.0 , -1 ) ) );

        ret_grad_phib[0] += phiBG * epsilon * (k * cos(k * xxb) * sin(k * yyb) * sin(k * zzb) * (pow(exp((sqrt(rb2) - rcutoff)/alpha) + 1.0 , -1))
          + sin(k*xxb) * sin(k*yyb) * sin(k*zzb) * (-1 * pow((exp((sqrt(rb2) - rcutoff)/alpha) + 1 ),-2)) * (xxb/alpha)*sqrt(rb2));

        ret_grad_phib[1] += phiBG * epsilon * (k * sin(k * xxb) * cos(k * yyb) * sin(k * zzb) * (pow(exp((sqrt(rb2) - rcutoff)/alpha)+1.0,-1)) 
          + sin(k*xxb) * sin(k*yyb) * sin(k*zzb) * (-1 * pow(exp((sqrt(rb2) - rcutoff)/alpha) + 1 , -2)) * (yyb/alpha)*sqrt(rb2));

        ret_grad_phib[2] += phiBG * epsilon * (k * sin(k * xxb) * sin(k * yyb) * cos(k * zzb) * (pow(exp((sqrt(rb2) - rcutoff)/alpha)+1.0,-1))
          + sin(k*xxb) * sin(k*yyb) * sin(k*zzb) * (-1 * pow(exp((sqrt(rb2) - rcutoff)/alpha) + 1, -2)) * (zzb/alpha)*sqrt(rb2));
      }

       else if (CCTK_EQUALS(scalar_profile, "fermi_cos")){
        double k;
        k = (2.*PI*N_sin)/domain_length; 
        /*r2 = xx*xx+yy*yy+zz*zz;*/

        ret_phia +=  phiBG * (1  + epsilon*(cos(k*xxa) * cos(k*yya) * cos(k*zza)) * (pow(exp((sqrt(ra2) - rcutoff)/alpha) + 1.0 , -1 ) ) );

        ret_grad_phia[0] += phiBG * epsilon * (k * -sin(k * xxa) * cos(k * yya) * cos(k * zza) * (pow(exp((sqrt(ra2) - rcutoff)/alpha) + 1.0 , -1))
          + cos(k*xxa) * cos(k*yya) * cos(k*zza) * (-1 * pow((exp((sqrt(ra2) - rcutoff)/alpha) + 1 ),-2)) * (xxa/alpha)*sqrt(ra2));

        ret_grad_phia[1] += phiBG * epsilon * (k * cos(k * xxa) * -sin(k * yya) * cos(k * zza) * (pow(exp((sqrt(ra2) - rcutoff)/alpha)+1.0,-1)) 
          + cos(k*xxa) * cos(k*yya) * cos(k*zza) * (-1 * pow(exp((sqrt(ra2) - rcutoff)/alpha) + 1 , -2)) * (yya/alpha)*sqrt(ra2));

        ret_grad_phia[2] += phiBG * epsilon * (k * cos(k * xxa) * cos(k * yya) * -sin(k * zza) * (pow(exp((sqrt(ra2) - rcutoff)/alpha)+1.0,-1))
          + cos(k*xxa) * cos(k*yya) * cos(k*zza) * (-1 * pow(exp((sqrt(ra2) - rcutoff)/alpha) + 1, -2)) * (zza/alpha)*sqrt(ra2));

        ret_phib +=  phiBG * (1  + epsilon*(cos(k*xxb) * cos(k*yyb) * cos(k*zzb)) * (pow(exp((sqrt(rb2) - rcutoff)/alpha) + 1.0 , -1 ) ) );

        ret_grad_phib[0] += phiBG * epsilon * (k * -sin(k * xxb) * cos(k * yyb) * cos(k * zzb) * (pow(exp((sqrt(rb2) - rcutoff)/alpha) + 1.0 , -1))
          + cos(k*xxb) * cos(k*yyb) * cos(k*zzb) * (-1 * pow((exp((sqrt(rb2) - rcutoff)/alpha) + 1 ),-2)) * (xxb/alpha)*sqrt(rb2));

        ret_grad_phib[1] += phiBG * epsilon * (k * cos(k * xxb) * -sin(k * yyb) * cos(k * zzb) * (pow(exp((sqrt(rb2) - rcutoff)/alpha)+1.0,-1)) 
          + cos(k*xxb) * cos(k*yyb) * cos(k*zzb) * (-1 * pow(exp((sqrt(rb2) - rcutoff)/alpha) + 1 , -2)) * (yyb/alpha)*sqrt(rb2));

        ret_grad_phib[2] += phiBG * epsilon * (k * cos(k * xxb) * cos(k * yyb) * -sin(k * zzb) * (pow(exp((sqrt(rb2) - rcutoff)/alpha)+1.0,-1))
          + cos(k*xxb) * cos(k*yyb) * cos(k*zzb) * (-1 * pow(exp((sqrt(rb2) - rcutoff)/alpha) + 1, -2)) * (zzb/alpha)*sqrt(rb2));

      }

      else if (CCTK_EQUALS(scalar_profile, "sin_decaying")){
        double k;
        k = (2.*PI*N_sin)/domain_length; 
        /*r2 = xx*xx+yy*yy+zz*zz;*/

        ret_phia += phiBG * (1 + epsilon* (sin(k*xxa)*sin(k*yya)*sin(k*zza) * (exp(-(ra2)/ (sigma_decay*sigma_decay) )) ) );

        ret_grad_phia[0] += phiBG * epsilon * k * cos(k*xxa)*sin(k*yya)*sin(k*zza) * (exp(-(ra2)/ (sigma_decay*sigma_decay) )) 
          + phiBG * epsilon * sin(k*xxa)*sin(k*yya)*sin(k*zza) * exp(-ra2 / (sigma_decay * sigma_decay)) * (-2*xxa / (sigma_decay * sigma_decay) );

        ret_grad_phia[1] += phiBG * epsilon * k * sin(k*xxa)*cos(k*yya)*sin(k*zza) * (exp(-(ra2)/ (sigma_decay*sigma_decay) )) 
          + phiBG * epsilon * sin(k*xxa)*sin(k*yya)*sin(k*zza) * exp(-ra2 / (sigma_decay * sigma_decay)) * (-2*yya / (sigma_decay * sigma_decay) );

        ret_grad_phia[2] += phiBG * epsilon * k * sin(k*xxa)*sin(k*yya)*cos(k*zza) * (exp(-(ra2)/ (sigma_decay*sigma_decay) )) 
          + phiBG * epsilon * sin(k*xxa)*sin(k*yya)*sin(k*zza) * exp(-ra2 / (sigma_decay * sigma_decay)) * (-2*zza / (sigma_decay * sigma_decay) );

        ret_phib += phiBG * (1 + epsilon* (sin(k*xxb)*sin(k*yyb)*sin(k*zzb) * (exp(-(rb2)/ (sigma_decay*sigma_decay) )) ) );

        ret_grad_phib[0] += phiBG * epsilon * k * cos(k*xxb)*sin(k*yyb)*sin(k*zzb) * (exp(-(rb2)/ (sigma_decay*sigma_decay) )) 
          + phiBG * epsilon * sin(k*xxb)*sin(k*yyb)*sin(k*zzb) * exp(-rb2 / (sigma_decay * sigma_decay)) * (-2*xxb / (sigma_decay * sigma_decay) );

        ret_grad_phib[1] += phiBG * epsilon * k * sin(k*xxb)*cos(k*yyb)*sin(k*zzb) * (exp(-(rb2)/ (sigma_decay*sigma_decay) )) 
          + phiBG * epsilon * sin(k*xxb)*sin(k*yyb)*sin(k*zzb) * exp(-rb2 / (sigma_decay * sigma_decay)) * (-2*yyb / (sigma_decay * sigma_decay) );

        ret_grad_phib[2] += phiBG * epsilon * k * sin(k*xxb)*sin(k*yyb)*cos(k*zzb) * (exp(-(rb2)/ (sigma_decay*sigma_decay) )) 
          + phiBG * epsilon * sin(k*xxb)*sin(k*yyb)*sin(k*zzb) * exp(-rb2 / (sigma_decay * sigma_decay)) * (-2*zzb / (sigma_decay * sigma_decay) );
      }
    } 
  }

  else if (CCTK_EQUALS(field_type, "Everywhere") ) {
    double r2, curr_phia, curr_phib;
    r2 = x*x+y*y+z*z;

    if (CCTK_EQUALS(scalar_profile, "Exponential") )
    {
      curr_phia = alphaE*exp(-betaE*r2)*(cos(gammaE*r2));
      //ret_phi += curr_phi;

      ret_phia = curr_phia;

      ret_grad_phia[0] = -2*alphaE*x*exp(-betaE*r2)*(betaE*cos(gammaE*r2) + gammaE*sin(gammaE*r2));
      ret_grad_phia[1] = -2*alphaE*y*exp(-betaE*r2)*(betaE*cos(gammaE*r2) + gammaE*sin(gammaE*r2));
      ret_grad_phia[2] = -2*alphaE*z*exp(-betaE*r2)*(betaE*cos(gammaE*r2) + gammaE*sin(gammaE*r2)); 

      curr_phib = alphaE*exp(-betaE*r2)*(cos(gammaE*r2));
      //ret_phi += curr_phi;

      ret_phib = curr_phib;

      ret_grad_phib[0] = -2*alphaE*x*exp(-betaE*r2)*(betaE*cos(gammaE*r2) + gammaE*sin(gammaE*r2));
      ret_grad_phib[1] = -2*alphaE*y*exp(-betaE*r2)*(betaE*cos(gammaE*r2) + gammaE*sin(gammaE*r2));
      ret_grad_phib[2] = -2*alphaE*z*exp(-betaE*r2)*(betaE*cos(gammaE*r2) + gammaE*sin(gammaE*r2)); 
    }
    
    else if (CCTK_EQUALS(scalar_profile, "OneOverR") )
    {
      curr_phia = alphaE/(gammaE*r2 + betaE);

      ret_phia = curr_phia;

      ret_grad_phia[0] = -2*alphaE*gammaE*x/powl(gammaE*r2 + betaE, 2);
      ret_grad_phia[1] = -2*alphaE*gammaE*y/powl(gammaE*r2 + betaE, 2);
      ret_grad_phia[2] = -2*alphaE*gammaE*z/powl(gammaE*r2 + betaE, 2); 

      curr_phib = alphaE/(gammaE*r2 + betaE);

      ret_phib = curr_phib;

      ret_grad_phib[0] = -2*alphaE*gammaE*x/powl(gammaE*r2 + betaE, 2);
      ret_grad_phib[1] = -2*alphaE*gammaE*y/powl(gammaE*r2 + betaE, 2);
      ret_grad_phib[2] = -2*alphaE*gammaE*z/powl(gammaE*r2 + betaE, 2);
    }

  }
 
   // add background phi and pi
    ret_phia += phiConst;
    ret_pia  += piConst;
    ret_phib += phiConst;
    ret_pib  += piConst;
    /*ret_phi += phiBG*Psi0m2;*/

    // add potential
    if (potential)
    {
      if ( CCTK_EQUALS(potential_type,"phiFourth") )
      ret_VV = STplankmass*STplankmass*STplankmass*STplankmass/(4*(0.5*STplankmass*STplankmass + 0.5*STcouplea*ret_phia*ret_phia + 0.5*STcoupleb* ret_phib*ret_phib )*(0.5*STplankmass*STplankmass + 0.5 *STcouplea*ret_phia*ret_phia + 0.5 *STcoupleb*ret_phib*ret_phib) )*(0.5 * STphiamass*STphiamass*ret_phia*ret_phia + 0.5 * STphibmass *STphibmass*ret_phib*ret_phib + lambdaphia/4.0 *ret_phia*ret_phia*ret_phia*ret_phia + lambdaphib/4.0 *ret_phib*ret_phib*ret_phib*ret_phib + gcouple/2.0 *ret_phia*ret_phia*ret_phib*ret_phib );
 
    }

    // set return values
    *phia = ret_phia;
    grad_phia[0] = ret_grad_phia[0];
    grad_phia[1] = ret_grad_phia[1];
    grad_phia[2] = ret_grad_phia[2];
    *pia = ret_pia;
    *phib = ret_phib;
    grad_phib[0] = ret_grad_phib[0];
    grad_phib[1] = ret_grad_phib[1];
    grad_phib[2] = ret_grad_phib[2];
    *pib = ret_pib;
    *VV = ret_VV;
}
