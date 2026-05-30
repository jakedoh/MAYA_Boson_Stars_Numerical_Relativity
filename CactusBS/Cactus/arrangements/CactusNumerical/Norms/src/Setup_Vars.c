#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "util_ErrorCodes.h"
#include "util_Table.h"


static const char * rcsid = "$Header$";

CCTK_FILEVERSION(Norms_Setup_Vars_c);


struct norms_opts {
  int active;
  int vi; /* variable index */
  CCTK_REAL bg; /* background value which gets subtracted */
};


static void normsgetopt (int const idx,
	                 const char * const optstring,
                         void * const opts)
{
  struct norms_opts * norms_opts;
  int table;
  int cnt;
  int ierr;
  CCTK_REAL bg;


  assert (idx >= 0 && idx < CCTK_NumVars());
  assert (opts);

  norms_opts = &((struct norms_opts *)opts)[idx];

  assert (! norms_opts->active);

  if (optstring) {
    assert (optstring);
    table = Util_TableCreateFromString (optstring);
    if (table < 0) {
      char * fullname = CCTK_FullName (idx);
      CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "The variable \"%s\" is ignored because it has an invalid option specification in the parameter \"Norms::gridfunctions...\"",
                  fullname);
      free (fullname);
      return;
    }
    assert (table >= 0);
  }

  norms_opts->active = 1;
  norms_opts->vi = idx;

  if (optstring) {
    cnt = Util_TableGetReal (table, &bg , "bg");

    if (cnt < 0) {
      norms_opts->bg = 0.; /* don't subtract any background */
    } else {
      norms_opts->bg = bg;
    }

    ierr = Util_TableDestroy (table);
    assert (!ierr);
  } else {
    norms_opts->bg = 0.; /* don't subtract any background */
  }

}





void Norms_Setup_Vars (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  int nvars;
  struct norms_opts * norms_opts_1st;
  struct norms_opts * norms_opts_2nd;

  int n,i;

  int ierr;

  /* init */
  *nr1stvars=0;
  *nr2ndvars=0;

  if (verbose>6)
    CCTK_INFO("here we are in Setup_vars");

  if (verbose>0)
    CCTK_VInfo(CCTK_THORNSTRING,"Starting Norms Computation at time %f",
               (double)cctkGH->cctk_time);

  for(i=0;i<max_nr_vars;i++) {
    varindices_1st[i]=-1;
    varindices_2nd[i]=-1;
    bgvals_1st[i]=0.;
    bgvals_2nd[i]=0.;
  }

  nvars = CCTK_NumVars();
  assert (nvars >= 0);
  norms_opts_1st = malloc (nvars * sizeof *norms_opts_1st);
  norms_opts_2nd = malloc (nvars * sizeof *norms_opts_2nd);
  assert (norms_opts_1st);
  assert (norms_opts_2nd);
  for (n=0; n<nvars; ++n) {
    norms_opts_1st[n].active = 0;
    norms_opts_2nd[n].active = 0;
  }
  ierr = CCTK_TraverseString
    (gridfunctions_1st, normsgetopt, norms_opts_1st, CCTK_GROUP_OR_VAR);
  assert (ierr >= 0);
  ierr = CCTK_TraverseString
    (gridfunctions_2nd, normsgetopt, norms_opts_2nd, CCTK_GROUP_OR_VAR);
  assert (ierr >= 0);

  if (verbose>0)
    CCTK_INFO("We will compute norms for the following variables");
  for (n=0; n<nvars; ++n) {
    if (norms_opts_1st[n].active) {
      assert(*nr1stvars<max_nr_vars);
      varindices_1st[*nr1stvars]=norms_opts_1st[n].vi;
      bgvals_1st[*nr1stvars]=norms_opts_1st[n].bg;
      if (verbose>0) {
        char * fullname = CCTK_FullName(varindices_1st[*nr1stvars]);
        CCTK_VInfo(CCTK_THORNSTRING,"    %s (1st order var) - %f",
                   fullname,
                   (double)bgvals_1st[*nr1stvars]);
        free (fullname);
      }
      *nr1stvars=*nr1stvars+1;
    }
    if (norms_opts_2nd[n].active) {
      assert(*nr2ndvars<max_nr_vars);
      varindices_2nd[*nr2ndvars]=norms_opts_2nd[n].vi;
      bgvals_2nd[*nr2ndvars]=norms_opts_2nd[n].bg;
      if (verbose>0) {
        char * fullname = CCTK_FullName(varindices_2nd[*nr2ndvars]);
        CCTK_VInfo(CCTK_THORNSTRING,"    %s (2nd order var) - %f",
                   fullname,
                   (double)bgvals_2nd[*nr2ndvars]);
        free (fullname);
      }
      *nr2ndvars=*nr2ndvars+1;
    }
  }
  if (verbose>2)
    fprintf(stderr,"   nr1stvars %d nr2ndvars %d\n",
            (int)*nr1stvars,
            (int)*nr2ndvars);

  if (verbose>4) {
    fprintf(stderr,"   The actual variable arrays look like\n");
    for (n=0;n<max_nr_vars;n++) {
      fprintf(stderr,"      i:%d (vi:%d vn:%s bg:%f) (vi:%d vn:%s bg:%f)\n",
              n,
              (int)varindices_1st[n],
              CCTK_VarName(varindices_1st[n]),
              (double)bgvals_1st[n],
              (int)varindices_2nd[n],
              CCTK_VarName(varindices_2nd[n]),
              (double)bgvals_2nd[n]);
    }
  }

  free (norms_opts_1st);
  free (norms_opts_2nd);
}
