#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "util_ErrorCodes.h"
#include "util_Table.h"


static const char * rcsid = "$Header: /numrelcvs/HerrmannCVS/Ylm_Decomp/src/Setup_Vars.c,v 1.14 2005/03/23 19:58:52 herrmann Exp $";

CCTK_FILEVERSION(Ylm_Decomp_Setup_Vars_c);


struct ylmopts {
  int active;
  int vi;
  int spnwght;
  int kick;
  int cmplgf;
  int imagIndex;
};



static void ylmdecomp_getopt (int const idx,
                    const char * const optstring,
                    void * const opts)
{
  struct ylmopts * ylmopts;
  int table;
  int cnt;
  int ierr;
  int spwght;
  int cmplgf;
  int imagInd;
  char imagVar[1024];


  assert (idx >= 0 && idx < CCTK_NumVars());
  assert (opts);

  ylmopts = &((struct ylmopts *)opts)[idx];

  assert (! ylmopts->active);

  if (! optstring) {
    char * fullname = CCTK_FullName (idx);
    CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                "The variable \"%s\" is ignored because it has no option specification in the parameter \"Ylm_Decomp::gridfunctions\"",
                fullname);
    free (fullname);
    return;
  }
  assert (optstring);
  table = Util_TableCreateFromString (optstring);
  if (table < 0) {
    char * fullname = CCTK_FullName (idx);
    CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                "The variable \"%s\" is ignored because it has an invalid option specification in the parameter \"Ylm_Decomp::gridfunctions\"",
                fullname);
    free (fullname);
    return;
  }
  assert (table >= 0);

  ylmopts->active = 1;
  ylmopts->vi = idx;

  cnt = Util_TableGetInt (table, &spwght , "sw");

  if (cnt < 0) {
    ylmopts->spnwght = -100;
  } else {
    ylmopts->spnwght = spwght;
  }
  
  cnt = Util_TableGetInt (table, &spwght , "kick");
  
  if (cnt < 0) {
    ylmopts->kick = 0;
  } else {
    ylmopts->kick = 1;
  }

  /* if this is a complex gf which varindex has the imag part? */
  //cnt = Util_TableGetInt (table, &cmplgf , "cmplx");
  cnt = Util_TableGetString (table, 1024, imagVar , "cmplx");

  //fprintf(stderr,"FIXME : extracting string retval=%d string=%s\n",cnt,imagVar);

  if (cnt < 0) {
    ylmopts->cmplgf = -1;
  } else {
    imagInd=CCTK_VarIndex(imagVar);
    //fprintf(stderr,"FIXME : imagIndex=%d\n",imagInd);
    ylmopts->cmplgf = 1;
    ylmopts->imagIndex = imagInd;
  }


  /* Don't transform if not output was requested */
  if (ylmopts->active) {
    if (ylmopts->spnwght==-100) {
      ylmopts->active = 0;
    }
  }

  ierr = Util_TableDestroy (table);
  assert (!ierr);
}











void Ylm_Setup_Vars (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  int nvars;
  struct ylmopts * ylmopts;

  int n,i;

  int ierr;

  if (verbose>6)
    CCTK_INFO("here we are in Setup_vars");

  if (verbose>0)
    CCTK_VInfo(CCTK_THORNSTRING,"Starting Ylm_Decomp at time %f",cctkGH->cctk_time);
  if (cctk_iteration<start_iteration || cctk_time<start_time)
  {
    if (verbose>1) {
      CCTK_VInfo(CCTK_THORNSTRING,"wrong iteration: %d T=%g start_iteration=%d start_time=%g",
		 cctk_iteration,cctk_time,start_iteration,start_time);
    }
    *do_nothing=1;
  }
  else {
    *do_nothing=0;
  }

  nvars = CCTK_NumVars();
  assert (nvars >= 0);
  ylmopts = malloc (nvars * sizeof *ylmopts);
  assert (ylmopts);
  for (n=0; n<nvars; ++n) {
    ylmopts[n].active = 0;
  }
  ierr = CCTK_TraverseString
    (gridfunctions, ylmdecomp_getopt, ylmopts, CCTK_GROUP_OR_VAR);
  assert (ierr >= 0);

  i=0;
  if (verbose>0) {
    CCTK_INFO("We will do the following variables");
  }
  for (n=0; n<nvars; ++n) {
    if (ylmopts[n].active) {
      varindices[i]=ylmopts[n].vi;
      spwght[i]=ylmopts[n].spnwght; // set spin weight
      if (spwght[i]!=-2) {
	CCTK_WARN(1,"not just spin weight -2 requested! -> performance penalty");
	*only_spwght_m2=0;
      }

      cmplgf[i]=ylmopts[n].cmplgf;// set indicator for complex gf
      kick[i]=ylmopts[n].kick; // set indicator for kick
      imagIndex[i]=ylmopts[n].imagIndex; // set gf index for imaginary part
      if (verbose>0) {
        char * var_name = CCTK_FullName (varindices[i]);
        if (cmplgf[i]<=0) {
          CCTK_VInfo(CCTK_THORNSTRING,"    %s spin-weight: %d real gf",
              var_name,spwght[i]);
        }
        else {
          char * imag_name = CCTK_FullName (imagIndex[i]);
          CCTK_VInfo(CCTK_THORNSTRING,"    %s spin-weight: %d complex gf: %d with imag part: %s",
                            var_name,spwght[i],
                            cmplgf[i],CCTK_FullName(imagIndex[i]));
          free(imag_name);
        }
        free(var_name);
      }
      // for a complex gf we put the imag index into the next entry!
      if (cmplgf[i]>0) {
        varindices[i+1]=imagIndex[i];
        i++;
      }
      i++;
    }
  }
  *nrdecompvars=i;
  assert(*nrdecompvars<max_nr_vars);

  // setup out_every_det to default value if needed
  char paramstr[2048]; // FIXME XXX let's hope that we never need more
  char default_out_every[2048]; // FIXME XXX let's hope that we never need more
  sprintf(default_out_every,"%d",out_every);

  for (int i=0;i<number_of_detectors;i++) {
    if (out_every_det[i]==-1) {
      sprintf(paramstr,"out_every_det[%d]",i);
      ierr=CCTK_ParameterSet(paramstr,"Ylm_Decomp",default_out_every);
      if (ierr<0) {
	CCTK_WARN(1,"parameter set for detector failed!");
      }
    }
    if (verbose>1) {
      CCTK_VInfo(CCTK_THORNSTRING,"out_every_det[%d]=%d",
		 i,out_every_det[i]);
    }
    
  }

  free (ylmopts);
}
