/* $Id: Baum.c,v 1.8 2005/10/20 14:01:19 herrmann Exp $ */
/*  check if we are happy with the ahfinderdirect mass */

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "math.h"
#include "stdio.h"
#include "stdlib.h"

# define M_PI           3.14159265358979323846  /* pi */



void TwoPuncturesMatter(CCTK_ARGUMENTS); // lives in TwoPunctures thorn




int Baum_write_output(char filename[2048],CCTK_INT create_file,
		      CCTK_INT iterations,
		      CCTK_REAL par_b,CCTK_REAL j_search, // j_search=J/mu
		      CCTK_REAL m_val_plus,CCTK_REAL m_val_minus,
		      CCTK_REAL admMass,
		      CCTK_REAL irr_mass1, CCTK_REAL irr_mass2,
		      CCTK_REAL E_bind, CCTK_INT converged)
{
  DECLARE_CCTK_PARAMETERS;

  FILE *file=NULL;
  
  if (verbose>2) {
    CCTK_INFO("in Baum_write_output");
  }


  if (create_file>0) {
    file=fopen(filename,"w");
    if (file==NULL) {
      CCTK_WARN(1,"couldn't write data file");
      return -1;
    }
    fprintf(file,"# Output data for J/mu=%.19g\n",j_search);
    fprintf(file,"# b:1 J/mu:2 E_bind/mu:3 m+:4 m-:5 M_1:6 M_2:7 M_ADM:8 Iterations:9 Converged:10\n");
  }
  else {
    file=fopen(filename,"a");
  }

  const CCTK_REAL mu=irr_mass1*irr_mass2/(irr_mass1+irr_mass2);
  const CCTK_REAL Emu=E_bind/mu;

  int ierr=fprintf(file,
		   "%.19g\t%.19g\t%.19g\t%.19g\t%.19g\t%.19g\t%.19g\t%.19g\t%d\t%d\n",
		   par_b,j_search,Emu,m_val_plus,m_val_minus,
		   irr_mass1,irr_mass2,admMass,
                   iterations,converged);
  if(ierr<0) {
    CCTK_WARN(1,"failed to write into file");
    return -1;
  }

  fclose(file);
  return 1;
}


int Baum_CallTwoPunctures(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  char val[2048];
  int ierr=-1;

  if (verbose>2) {
    CCTK_VInfo(CCTK_THORNSTRING,"calling TwoPuncturesMatterzypbaum: m+=%.19g m-=%.19g",
	       *m_val_plus,*m_val_minus);
  }

  // check that masses are greater than 0:
  if (*m_val_plus<0 || *m_val_minus<0) {
     CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                "masses have to be >0! m+=%.19g m-=%.19g",
                *m_val_plus,*m_val_minus);
  }

  // update TwoPunctures parameters
  // m+
  sprintf(val,"%.19g",*m_val_plus);
  /* update m+ */
  ierr=CCTK_ParameterSet("par_m_plus","TwoPuncturesMatterzypbaum",val);
  if (ierr!=0) CCTK_WARN(1,"failed to set m+ params");
  if (verbose>2) {
    CCTK_VInfo(CCTK_THORNSTRING,"set m+=%.19g",*m_val_plus);
  }
	
  // m-
  if (verbose>2) {
    CCTK_VInfo(CCTK_THORNSTRING,"new value m-=%.19g",*m_val_minus);
  }
  sprintf(val,"%.19g",*m_val_minus);
  ierr=CCTK_ParameterSet("par_m_minus","TwoPuncturesMatterzypbaum",val);
  if (ierr!=0) CCTK_WARN(1,"failed to set m- params");
  if (verbose>2) {
    CCTK_VInfo(CCTK_THORNSTRING,"set m-=%.19g",*m_val_minus);
  }

  // call twopunctures to do the solve
  CCTK_INFO("calling TwoPuncturesMatterzypbaum");
  TwoPuncturesMatter(CCTK_PASS_CTOC);
  CCTK_INFO("done calling TwoPuncturesMatterzypbaum");
  return 1;
}






int Baum_update_mplus_mminus(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTS;


  CCTK_INFO("compute new m_val_plus, m_val_minus");


  if (CCTK_Equals(search_method,"bisection")) {
    CCTK_INFO("bisection search");
    *m_val_plus= (*m1_max+ *m1_min)/2.;
    *m_val_minus=(*m2_max+ *m2_min)/2.;
  }
  else if (CCTK_Equals(search_method,"secant")) {
    CCTK_INFO("secant method");
    const CCTK_REAL diff1=*irr_mass1- *m1_search;
    const CCTK_REAL diff2=*irr_mass2- *m2_search;
    const CCTK_REAL dmp=(*last_m_plus-  *m_val_plus) *diff1/(*irr_mass1-*last_irr_mass1);
    const CCTK_REAL dmm=(*last_m_minus- *m_val_minus)*diff2/(*irr_mass2-*last_irr_mass2);
    *last_m_plus=*m_val_plus;
    *last_m_minus=*m_val_minus;
    *last_irr_mass1=*irr_mass1;
    *last_irr_mass2=*irr_mass2;
    *m_val_plus=*m_val_plus+dmp;
    *m_val_minus=*m_val_minus+dmm;
    
    CCTK_VInfo(CCTK_THORNSTRING,"last_m+=%g m+=%g dm+=%g last_M+=%g M+=%g",
	       *last_m_plus,*m_val_plus,dmp,*last_irr_mass1,*irr_mass1);
    CCTK_VInfo(CCTK_THORNSTRING,"last_m-=%g m-=%g dm-=%g last_M-=%g M-=%g",
	       *last_m_minus,*m_val_minus,dmm,*last_irr_mass2,*irr_mass2);
  }
  else {
    CCTK_WARN(0,"unknown search_method!");
  }
  
  return 1;
}



int Baum_change_npoints(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTS;

  char val[2048]; // XXX hope this is long enough
  CCTK_INT forcemedres=0;
  CCTK_INT forcehighres=0;
  CCTK_INT npt_A,npt_B,npt_phi, ierr=-1;



  CCTK_INFO("Baum_change_npoints");

  if (touch_npts>0)
  {
    if ( (fabs(*m1_search- *irr_mass1) < switch_high_acc)
	 || forcehighres>0) {      
      CCTK_INFO("reached accuracy to switch to high res solver");
      forcehighres=1; /* don't go back to lower res. solver */
      forcemedres=0;
      npt_A=96;
      npt_B=96;
      npt_phi=32;
    }
    else if ( (fabs(*m1_search- *irr_mass1) < switch_med_acc)
	      || forcemedres>0) {
      CCTK_INFO("reached accuracy to switch to medium res solver");
      forcemedres=1; /* don't go back to lower res. solver */
      npt_A=96;
      npt_B=93;
      npt_phi=32;
    }
    else {
      CCTK_INFO("still using lowres solver");
      npt_A=96;
      npt_B=96;
      npt_phi=32;
    }
    sprintf(val,"%d",npt_A);
    ierr=CCTK_ParameterSet("npoints_A","TwoPuncturesMatterzypbaum",val);
    if(ierr<0) CCTK_WARN(1,"failed to set npoints_A");
    sprintf(val,"%d",npt_B);
    ierr=CCTK_ParameterSet("npoints_B","TwoPuncturesMatterzypbaum",val);
    if(ierr<0) CCTK_WARN(1,"failed to set npoints_A");
    sprintf(val,"%d",npt_phi);
	  ierr=CCTK_ParameterSet("npoints_phi","TwoPuncturesMatterzypbaum",val);
	  if(ierr<0) CCTK_WARN(1,"failed to set npoints_A");
  } // change npoints

  return 1;
} 


int Baum_update_m_min_max(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTS;

  CCTK_INFO("Baum_update_m_min_max");

  if (*irr_mass1>=*m1_search) {
    *m1_max=*m_val_plus;
    *irr_mass1_max=*irr_mass1;
  }
  else {
    *m1_min=*m_val_plus;
    *irr_mass1_min=*irr_mass1;
  }
  if (*irr_mass2>=*m2_search) {
    *m2_max=*m_val_minus;
    *irr_mass2_max=*irr_mass2;
  }
  else {
    *m2_min=*m_val_minus;
    *irr_mass2_min=*irr_mass2;
  }
  
  CCTK_VInfo(CCTK_THORNSTRING,"next we have to search m+ in range : %.19g %.19g",
	     *m1_min,*m1_max);
  CCTK_VInfo(CCTK_THORNSTRING,"next we have to search m- in range : %.19g %.19g",
	     *m2_min,*m2_max);
  
  CCTK_VInfo(CCTK_THORNSTRING,"next we have to search mirr1 in range : %.19g %.19g",
	     *irr_mass1_min,*irr_mass1_max);      
  CCTK_VInfo(CCTK_THORNSTRING,"next we have to search mirr2 in range : %.19g %.19g",
	     *irr_mass2_min,*irr_mass2_max);      

  return 1;
}


int Baum_get_masses(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTS;

  CCTK_INFO("Baum_get_masses");

  *irr_mass1=*(CCTK_REAL *)CCTK_VarDataPtr(cctkGH,0,"TwoPuncturesMatterzypbaum::admMass_plus");
  *irr_mass2=*(CCTK_REAL *)CCTK_VarDataPtr(cctkGH,0,"TwoPuncturesMatterzypbaum::admMass_minus");

  if (verbose>1) {
    CCTK_VInfo(CCTK_THORNSTRING,
	       "ADM mass from other universes: m+: %.19g m-: %.19g",*irr_mass1,*irr_mass2);
  }

  if (use_irreducible_masses>0) {
    CCTK_INFO("searching for Christodolou masses!");
    const CCTK_REAL *s1x=(CCTK_REAL *) CCTK_ParameterGet("par_S_plus[0]", "TwoPuncturesMatterzypbaum", NULL);
    const CCTK_REAL *s1y=(CCTK_REAL *) CCTK_ParameterGet("par_S_plus[1]", "TwoPuncturesMatterzypbaum", NULL);
    const CCTK_REAL *s1z=(CCTK_REAL *) CCTK_ParameterGet("par_S_plus[2]", "TwoPuncturesMatterzypbaum", NULL);

    const CCTK_REAL *s2x=(CCTK_REAL *) CCTK_ParameterGet("par_S_minus[0]", "TwoPuncturesMatterzypbaum", NULL);
    const CCTK_REAL *s2y=(CCTK_REAL *) CCTK_ParameterGet("par_S_minus[1]", "TwoPuncturesMatterzypbaum", NULL);
    const CCTK_REAL *s2z=(CCTK_REAL *) CCTK_ParameterGet("par_S_minus[2]", "TwoPuncturesMatterzypbaum", NULL);

    const CCTK_REAL s1=sqrt(*s1x* *s1x+*s1y* *s1y+ *s1z* *s1z);
    const CCTK_REAL s2=sqrt(*s2x* *s2x+*s2y* *s2y+ *s2z* *s2z);

    // invert christodolou: M_irr=1/sqrt(2)*sqrt(M^2+(M^4-S^2)^(1/2))
    CCTK_REAL mi=*irr_mass1;
    CCTK_REAL mi_2=mi*mi;
    CCTK_REAL mi_4=mi_2*mi_2;
    *irr_mass1=1/sqrt(2)*sqrt(mi_2 +sqrt(mi_4-s1*s1));
    CCTK_VInfo(CCTK_THORNSTRING,"  mass1: mtot=%g s1=%g mirr=%g",mi,s1,*irr_mass1);
    mi=*irr_mass2;
    mi_2=mi*mi;
    mi_4=mi_2*mi_2;
    *irr_mass2=1/sqrt(2)*sqrt(mi_2 +sqrt(mi_4-s2*s2));
    CCTK_VInfo(CCTK_THORNSTRING,"  mass2: mtot=%g s2=%g mirr=%g",mi,s2,*irr_mass2);
  }

  // get total adm_mass out of TwoPunctures
  if (use_admrad>0) {
    CCTK_INFO("use finite radius for ADM Mass computation");
    *admMass=*(CCTK_REAL *)CCTK_VarDataPtr(cctkGH,0,"TwoPuncturesMatterzypbaum::admMassRad");
  }
  else {
    CCTK_INFO("use compactified infinity for ADM Mass computation");
    *admMass=*(CCTK_REAL *)CCTK_VarDataPtr(cctkGH,0,"TwoPuncturesMatterzypbaum::admMass");
  }
  if (verbose>1) {
    CCTK_VInfo(CCTK_THORNSTRING,"ADM Mass = %g",*admMass);
  }
  
  *E_bind=*admMass-*irr_mass1-*irr_mass2;
  
  if (verbose>0) {
    CCTK_VInfo(CCTK_THORNSTRING,"Binding Energy = %g",*E_bind);
  }

  return 1;
}


int Baum_check_convergence(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTS;

  char filename[2048]; // XXX hope we never need more

  CCTK_INFO("Baum_check_convergence");

  const CCTK_INT myproc=CCTK_MyProc(cctkGH);
  const CCTK_INT nprocs=CCTK_nProcs(cctkGH);

  if (nprocs>1) {
    sprintf(filename,"%s/baum.cpu%04d.dat",out_dir,myproc);
  }
  else {
    sprintf(filename,"%s/baum.dat",out_dir);
  }
  if (verbose>2) {
    CCTK_VInfo(CCTK_THORNSTRING,
	       "output filename is %s",filename);
  }

  // check if parameters are too close to last values
  if (fabs(*m1_max- *m1_min)<eps_mass_param) {
    CCTK_VInfo(CCTK_THORNSTRING,"m1 reached eps mass param!");
    *found=0;
  }
  if (fabs(*m2_max- *m2_min)<eps_mass_param) {
    CCTK_VInfo(CCTK_THORNSTRING,"m2 reached eps mass param!");
    *found=0;
  }
  
  /* both parameters have "converged", but not in the actual mass values!! */
  CCTK_INT bare_mass_conv=0;
  if (CCTK_EQUALS(search_method,"bisection")) {
    if (fabs(*m1_min- *m1_min)<eps_mass_param &&
	fabs(*m2_max- *m2_min)<eps_mass_param ) {
      bare_mass_conv=1;
    }
  }
  if (CCTK_EQUALS(search_method,"secant")) {
    if (fabs(*m_val_plus- *last_m_plus)<eps_mass_param &&
	fabs(*m_val_minus- *last_m_minus)<eps_mass_param ) {
      bare_mass_conv=1;
    }
  }

  if (bare_mass_conv) {
    CCTK_VInfo(CCTK_THORNSTRING,"both m1 and m2 reached eps mass param!");
    *found=1;
    
    CCTK_INT create_file=0;
    if (*first_ever==1) {
      create_file=1;
    }
    CCTK_INFO("we have only converged in the parameter search, not in the actual values");
    Baum_write_output(filename,create_file,
		      *iterations,
		      *par_b,*j_sequence,*m_val_plus,*m_val_minus,
		      *admMass,*irr_mass1,*irr_mass2,*E_bind,0);
  }
  
  // we have converged to final answer?
  if (fabs(*m1_search- *irr_mass1) < eps_mass &&
      fabs(*m2_search- *irr_mass2) < eps_mass ) 
  {
    if(verbose>0) {
      CCTK_INFO("DONE -> NAILED IT !!!");
      CCTK_VInfo(CCTK_THORNSTRING,
		 "reached value in %d iterations",*iterations);
      CCTK_VInfo(CCTK_THORNSTRING,
		 "m1_search-M_1=%g eps_mass=%g  m1_max-m1_min=%g eps_mass_param=%g",
		 fabs(*m1_search- *irr_mass1),eps_mass,fabs(*m1_max- *m1_min),eps_mass_param);
      CCTK_VInfo(CCTK_THORNSTRING,
		 "m2_search-M_2=%g eps_mass=%g  m2_max-m2_min=%g eps_mass_param=%g",
		 fabs(*m2_search- *irr_mass2),eps_mass,fabs(*m2_max- *m2_min),eps_mass_param);
      CCTK_INFO("Final Results");
      CCTK_VInfo(CCTK_THORNSTRING,"ADM Mass=%.19g",*admMass);
      CCTK_VInfo(CCTK_THORNSTRING,"Irr Mass 1=%.19g",*irr_mass1);
      CCTK_VInfo(CCTK_THORNSTRING,"Irr Mass 2=%.19g",*irr_mass2);
      CCTK_VInfo(CCTK_THORNSTRING,"m+=%.19g",*m_val_plus);
      CCTK_VInfo(CCTK_THORNSTRING,"m-=%.19g",*m_val_minus);
      CCTK_VInfo(CCTK_THORNSTRING,"E_b=%.19g",*E_bind);
    }
    *found=1;
    
    /* dump data */
    CCTK_INT create_file=0;
    if (*first_ever==1) {
      create_file=1;
    }
    Baum_write_output(filename,create_file,
		      *iterations,
		      *par_b,*j_sequence,*m_val_plus,*m_val_minus,
		      *admMass,*irr_mass1,*irr_mass2,*E_bind,1);
  }
  else if (*found == 0) /* this takes care if we have flagged this iteration as problematic */
  {
    if (verbose>2) {
      CCTK_INFO("need to try other ID");
    }
    *iterations=*iterations+1;
    *found=0;
    
    Baum_change_npoints(CCTK_PASS_CTOC); // update TwoPunctures npoints_A,... resolution params
  } // not converged yet

  return 1;
}
  


int Baum_setup_b_p(CCTK_REAL b_start, CCTK_REAL delta_b, CCTK_INT i, CCTK_REAL j_sequence, CCTK_REAL *par_b)
{
  DECLARE_CCTK_PARAMETERS;
  char val[2048]; // XXX hope we never need more
  CCTK_INT ierr;

  CCTK_INFO("Baum_setup_b_p");

  // separation for this member of the sequence
  *par_b=b_start +i*delta_b;
  sprintf(val,"%.19g",*par_b);
  ierr=CCTK_ParameterSet("par_b","TwoPuncturesMatterzypbaum",val);
  if (ierr!=0) CCTK_WARN(1,"failed to set par_b param");
  CCTK_VInfo(CCTK_THORNSTRING,"setting par_b=%s",val);


  // spin terms in z direction
  const CCTK_REAL *s1z=(CCTK_REAL *) CCTK_ParameterGet("par_S_plus[2]", "TwoPuncturesMatterzypbaum", NULL);
  const CCTK_REAL *s2z=(CCTK_REAL *) CCTK_ParameterGet("par_S_minus[2]", "TwoPuncturesMatterzypbaum", NULL);

  CCTK_REAL j_z=j_sequence;
  if (search_j_ang_mom) {
    j_z=j_sequence-*s1z-*s2z;
  }
  const CCTK_REAL par_p=j_z/(2.* *par_b); // J=2*b*p -> p=J/(2*b)
  sprintf(val,"%.19g",par_p);
  ierr=CCTK_ParameterSet("par_p_plus[1]","TwoPuncturesMatterzypbaum",val);
  if (ierr!=0) CCTK_WARN(1,"failed to set par_p_plus[1] params");
  CCTK_VInfo(CCTK_THORNSTRING,"setting par_p_plus[1]=%s",val);

  sprintf(val,"%.19g",-par_p);
  ierr=CCTK_ParameterSet("par_p_minus[1]","TwoPuncturesMatterzypbaum",val);
  if (ierr!=0) CCTK_WARN(1,"failed to set par_p_minus[1] params");
  CCTK_VInfo(CCTK_THORNSTRING,"setting par_p_minus[1]=%s",val);
  
  return 0;
}




int Baum_setup_initial_masses(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  CCTK_INT ierr;

  CCTK_INFO("Baum_setup_initial_masses: solve for bounding values");
  // setup of the initial range for m to search in
  *m1_max=*m1_guess+percent_range* *m1_guess;
  *m1_min=*m1_guess-percent_range* *m1_guess;
  
  *m2_max=*m2_guess+percent_range* *m2_guess;
  *m2_min=*m2_guess-percent_range* *m2_guess;
  
  CCTK_VInfo(CCTK_THORNSTRING,"trying to bracket solution in: m1_min: %.19g m1_max: %.19g",
	     *m1_min,*m1_max);
  CCTK_VInfo(CCTK_THORNSTRING,"                               m2_min: %.19g m2_max: %.19g",
	     *m2_min,*m2_max);
  
  // first we check the corner cases at m1_max and m1_min and make sure that
  // we actually bracket the solution!
  // update m+,m- and call TwoPunctures()
  *m_val_plus=*m1_max;
  *m_val_minus=*m2_max;
  ierr=Baum_CallTwoPunctures(CCTK_PASS_CTOC);
  *irr_mass1_max=*(CCTK_REAL *)CCTK_VarDataPtr(cctkGH,0,"TwoPuncturesMatterzypbaum::admMass_plus");
  *irr_mass2_max=*(CCTK_REAL *)CCTK_VarDataPtr(cctkGH,0,"TwoPuncturesMatterzypbaum::admMass_minus");

  *m_val_plus=*m1_min;
  *m_val_minus=*m2_min;  
  ierr=Baum_CallTwoPunctures(CCTK_PASS_CTOC);
  *irr_mass1_min=*(CCTK_REAL *)CCTK_VarDataPtr(cctkGH,0,"TwoPuncturesMatterzypbaum::admMass_plus");
  *irr_mass2_min=*(CCTK_REAL *)CCTK_VarDataPtr(cctkGH,0,"TwoPuncturesMatterzypbaum::admMass_minus");
  
  CCTK_VInfo(CCTK_THORNSTRING,"irr_mass1_min=%.19g irr_mass1_max=%.19g",
	     *irr_mass1_min,*irr_mass1_max);
  CCTK_VInfo(CCTK_THORNSTRING,"irr_mass2_min=%.19g irr_mass2_max=%.19g",
	     *irr_mass2_min,*irr_mass2_max);
  
  if (*m1_search> *irr_mass1_max) {
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "m1_search>irr_mass1_max! m1_search: %.19g irr_mass1_min: %.19g irr_mass1_max: %.19g. m1_min: %.19g m1_max: %.19g",
               *m1_search, *irr_mass1_min, *irr_mass1_max, *m1_min, *m1_max);
  }
  if (*m1_search< *irr_mass1_min) {
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "m1_search<irr_mass1_min! m1_search: %.19g irr_mass1_min: %.19g irr_mass1_max: %.19g. m1_min: %.19g m1_max: %.19g",
               *m1_search, *irr_mass1_min, *irr_mass1_max, *m1_min, *m1_max);
  }
  if (*m2_search> *irr_mass2_max) {
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "m2_search>irr_mass2_max! m2_search: %.19g irr_mass2_min: %.19g irr_mass2_max: %.19g. m2_min: %.19g m2_max: %.19g",
               *m2_search, *irr_mass2_min, *irr_mass2_max, *m2_min, *m2_max);
  }
  if (*m2_search< *irr_mass2_min) {
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "m2_search<irr_mass2_min! m2_search: %.19g irr_mass2_min: %.19g irr_mass2_max: %.19g. m2_min: %.19g m2_max: %.19g",
               *m2_search, *irr_mass2_min, *irr_mass2_max, *m2_min, *m2_max);
  }

  if (CCTK_EQUALS(search_method,"bisection")) {
    // intial starting values for this sequence member
    // note that this is redundant as it gets set again later
    *m_val_plus =(*m1_min+*m1_max)/2.;
    *m_val_minus=(*m2_min+*m2_max)/2.;
  }

  else if (CCTK_EQUALS(search_method,"secant")) {
    // arbitrarily choose minimum values as last ones!
    *last_irr_mass1=*irr_mass1_min;
    *last_irr_mass2=*irr_mass2_min;
    *last_m_plus = *m1_min;
    *last_m_minus= *m2_min;
    *irr_mass1=*m1_max;
    *irr_mass2=*m2_max;
    *m_val_plus =*m1_max;
    *m_val_minus=*m2_max;
  }

  if (verbose>1) {
    CCTK_VInfo(CCTK_THORNSTRING,"start search for m+=%.19g m-=%.19g",
	       *m_val_plus,*m_val_minus);
  }
  
  CCTK_INFO("done with initial solves!");

  return 0;
}


int Baum_parallel_setup(CCTK_ARGUMENTS, CCTK_REAL *delta_b, CCTK_INT *sep_on_cpu)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_INFO("Baum_parallel_setup");

  const CCTK_INT nruns=number_of_sequence_members;
  const CCTK_INT myproc=CCTK_MyProc(cctkGH);
  const CCTK_INT nprocs=CCTK_nProcs(cctkGH);

  if (sep_on_cpu ==NULL) CCTK_WARN(0,"no storage for sep_on_cpu!");

  // distribute over CPUs
  for (int i=0;i<nruns;i++) {
    sep_on_cpu[i]=i % nprocs;
    CCTK_VInfo(CCTK_THORNSTRING,"parallel distribute: run %d on CPU %d out of %d CPUs",
	       i,sep_on_cpu[i],nprocs);
  }

  // first the special case of nruns=1
  if (nruns==1) {
    CCTK_INFO("only 1 run requested, b_max is ignored, b_min is computed!");
    *delta_b=0;
  }
  else {
    *delta_b=(b_max-b_min)/(nruns-1);
  }

  CCTK_VInfo(CCTK_THORNSTRING,
	     "nruns=%d cpu=%d delta_b=%g",
	     nruns,myproc,*delta_b);
  
  return 0;
}

int Baum_setup_ang_mom(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTS;

  CCTK_INFO("in Baum_setup_ang_mom");

  *j_sequence=j_mu*eta;
  // only want to nail the mass for b_min...
  if (number_of_sequence_members==1) {
    const CCTK_REAL *p1y=(CCTK_REAL *) CCTK_ParameterGet("par_P_minus[1]", "TwoPuncturesMatterzypbaum", NULL);
    const CCTK_REAL *p2y=(CCTK_REAL *) CCTK_ParameterGet("par_P_plus[1]", "TwoPuncturesMatterzypbaum", NULL);

    *j_sequence=(-b_min* *p1y+b_min* *p2y);
    CCTK_VInfo(CCTK_THORNSTRING,"searching for fixed J=%.19g b=%.19g p1y=%.19g p2y=%.19g",*j_sequence,b_min,*p1y,*p2y);
  }
  else {
    CCTK_VInfo(CCTK_THORNSTRING,"searching for J/mu=%.19g -> J=%.19g",j_mu,*j_sequence);
  }
  return 1;
}



void Baum(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  FILE *file;
  char filename[2048];
  int ierr;
  char val[256];
  const CCTK_REAL *parval;
  int i;
  CCTK_REAL mass=0,area=0;
  // ---------------------------------------------

  CCTK_INFO("Starting Baum");

  // cpu specific b_start and delta_b
  const CCTK_INT myproc=CCTK_MyProc(cctkGH);
  CCTK_REAL delta_b=0;
  CCTK_INT *sep_on_cpu;
  const CCTK_INT nruns=number_of_sequence_members;
  sep_on_cpu=(CCTK_INT *) malloc(sizeof(CCTK_INT)*nruns);
  Baum_parallel_setup(CCTK_PASS_CTOC, &delta_b, sep_on_cpu);

  CCTK_VInfo(CCTK_THORNSTRING,"searching for eta=%g and m_search=%g",eta,m_search);
  // work out the search masses
  // eta=m1*m2/(m1+m2)^2
  // hence: m2=-((2*eta-1)*m1+sqrt(1-4*eta)*m1)/(2*eta)
  *m1_search=m_search;
  *m2_search=-((2*eta-1)*m_search+sqrt(1-4*eta)*m_search)/(2*eta);
  CCTK_VInfo(CCTK_THORNSTRING,"searching for mirr1: %.19g mirr2: %.19g",*m1_search,*m2_search);

  // get the angular momentum, again assume that m1+m2=1 and then mu=eta
  Baum_setup_ang_mom(CCTK_PASS_CTOC);

  if (par_m1_guess==-1) {*m1_guess=*m1_search;}
  else {*m1_guess=par_m1_guess;}
  if (par_m2_guess==-1) {*m2_guess=*m2_search;}
  else {*m2_guess=par_m2_guess;}
  CCTK_VInfo(CCTK_THORNSTRING,"using guess m1: %.19g m2: %.19g",*m1_guess,*m2_guess);

  // loop over different separations
  *first_ever=1;
  for (i=0;i<nruns;i++)
  {
    // check if we do this separation on this CPU
    if (sep_on_cpu[i]!=myproc) {
      CCTK_VInfo(CCTK_THORNSTRING,"skip separation number  %d on CPU %d",
		 i,myproc);
      continue;
    }
    CCTK_VInfo(CCTK_THORNSTRING,"analyse separation number %d on CPU %d",
	       i,myproc);

    // setup of b and p parameters for TwoPunctures
    Baum_setup_b_p(b_min, delta_b, i, *j_sequence, par_b);

    // setup of the initial range for m to search in
    Baum_setup_initial_masses(CCTK_PASS_CTOC);

    // loop until converged to criterium we want, i.e. is M_1=m1_search and M_2=M_1*mass_ratio
    *iterations=0;
    *found=0;
    while (*found==0)
    {
      if (verbose>1) {
	CCTK_VInfo(CCTK_THORNSTRING,"iteration=%d",*iterations);
      }


      ierr=Baum_update_mplus_mminus(CCTK_PASS_CTOC);
      if (ierr<0) {
	CCTK_WARN(1,"can't get updated m_plus, m_minus values");
      }

      // solve twopunctures
      int ierr=Baum_CallTwoPunctures(CCTK_PASS_CTOC);
      if (ierr<0) {
	CCTK_WARN(1,"solving in TwoPuncturesMatterzypbaum failed");
      }

      // update masses
      ierr=Baum_get_masses(CCTK_PASS_CTOC);

      // check convergence
      Baum_check_convergence(CCTK_PASS_CTOC);

      // update m_min and m_max for bisection search
      Baum_update_m_min_max(CCTK_PASS_CTOC);

    } // end of while loop to check for mass convergence

    *first_ever=0;
  } // end of sequence loop

  free(sep_on_cpu); //clean up
}
