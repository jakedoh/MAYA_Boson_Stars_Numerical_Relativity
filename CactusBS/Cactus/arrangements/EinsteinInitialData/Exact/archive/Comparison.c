/* Please god comment me ... */
/* $Header$ */

#include "cctk.h"

static int comp_initialized = 0;
static int comp_active = -1;

static int comp_metric = 0;
static int comp_curv   = 0;
static int comp_gauge  = 0;

#define f_comparisonmetric FORTRAN_NAME(comparisonmetric_,COMPARISONMETRIC,comparisonmetric)
void FMODIFIER f_comparisonmetric(FORT_ARGS_PROTO, PROTO_FORT_FIELDS,
				  double *,
				  double *,
				  double *,
				  double *,
				  double *,
				  double *);
#define f_comparisoncurvature FORTRAN_NAME(comparisoncurvature_,COMPARISONCURVATURE,comparisoncurvature)
void FMODIFIER f_comparisoncurvature(FORT_ARGS_PROTO, PROTO_FORT_FIELDS,
				  double *,
				  double *,
				  double *,
				  double *,
				  double *,
				  double *);
#define f_comparisongauge FORTRAN_NAME(comparisongauge_,COMPARISONGAUGE,comparisongauge)
void FMODIFIER f_comparisongauge(FORT_ARGS_PROTO, PROTO_FORT_FIELDS,
				  double *,
				  double *,
				  double *,
				  double *);

static int comptrips=1;

void Comparison(pGH *in) {
  void DoComparison(pGH *high, pGH *med, pGH *low,
		    int gfi, double *he, double *me, double *le);

  int i,j,k;
  double *h[6];
  double *m[6];
  double *l[6];

  int compevery = 1;
  
  pGH *high, *med, *low;

  /* Grab the approriate grid heirarchies */
  if (in->convlevel != 0) return;

  high = in;
  med = GetGHbyLevel(1,in->level, in->mglevel);
  low = GetGHbyLevel(2,in->level, in->mglevel);

  /* Only compare when grids are aligned. */
  if (med) compevery = compevery * 2;
  if (low) compevery = compevery * 2;

  if (iteration % compevery != 0) return;

  /* See if I'm active */
  if (comp_active < 0) {
    comp_active = Contains("comparison","yes");
  }
  if (!comp_active) return;

  if (!comp_initialized) {
    pGF *tGF;
    comp_initialized =1;
    for (i=0; i<high->ngridFuncs;i++) {
      tGF = high->gridFuncs[i];
      if (Contains_bounded("compfields",tGF->name)) {
	tGF->do_comparison = 1;
	if (tGF->gfno >= GFI_GXX &&
	    tGF->gfno <= GFI_GZZ)
	  comp_metric = 1;

	if (tGF->gfno >= GFI_HXX &&
	    tGF->gfno <= GFI_HZZ)
	  comp_curv = 1;

	if (tGF->gfno == GFI_ALP   ||
	    tGF->gfno == GFI_BETAX ||
	    tGF->gfno == GFI_BETAY ||
	    tGF->gfno == GFI_BETAZ)
	  comp_gauge = 1;
      }
    }
  }
  if (!(comp_metric || comp_gauge || comp_curv)) {
    comp_active = 0;
    printf ("You must compare gauge, curvature, or metric");
    return;
  }

  /* OK at this point we know we are going to do something, so
     allocate the memory for the comparisons */
  for (i=0;i<6;i++) {
    h[i] = (double *)malloc(high->npoints*sizeof(double));
    if (med)
      m[i] = (double *)malloc( med->npoints*sizeof(double));
    else
      m[i] = NULL;
    if (low) 
      l[i] = (double *)malloc( low->npoints*sizeof(double));
    else
      l[i] = NULL;
  }

  EnableGFDataStorage(high, high->gridFuncs[high->ngridFuncs-1]);
  if (med) 
    EnableGFDataStorage( med,  med->gridFuncs[med->ngridFuncs-1]);
  if (low) 
    EnableGFDataStorage( low,  low->gridFuncs[low->ngridFuncs-1]);

  if (comp_metric) {
    SetupFortranArrays(high);
    f_comparisonmetric(FORT_ARGS(high), PASS_FORT_FIELDS(high),
		       h[0],h[1],h[2],h[3],h[4],h[5]);

    if (med) {
      SetupFortranArrays(med);
      f_comparisonmetric(FORT_ARGS(med), PASS_FORT_FIELDS(med),
			 m[0],m[1],m[2],m[3],m[4],m[5]);
    }
    
    if (low) {
      SetupFortranArrays(low);
      f_comparisonmetric(FORT_ARGS(low), PASS_FORT_FIELDS(low),
			 l[0],l[1],l[2],l[3],l[4],l[5]);
    }

    for (i=GFI_GXX; i<=GFI_GZZ;i++) {
      if (high->gridFuncs[i]->do_comparison) {
	DoComparison(high,med,low,i,
		     h[i-GFI_GXX],m[i-GFI_GXX],l[i-GFI_GXX]);
      }
    }
  }

  if (comp_curv) {
    SetupFortranArrays(high);
    f_comparisoncurvature(FORT_ARGS(high), PASS_FORT_FIELDS(high),
		       h[0],h[1],h[2],h[3],h[4],h[5]);

    if (med) {
      SetupFortranArrays(med);
      f_comparisoncurvature(FORT_ARGS(med), PASS_FORT_FIELDS(med),
			    m[0],m[1],m[2],m[3],m[4],m[5]);
    } 

    if (low) {
      SetupFortranArrays(low);
      f_comparisoncurvature(FORT_ARGS(low), PASS_FORT_FIELDS(low),
			    l[0],l[1],l[2],l[3],l[4],l[5]);
    }

    for (i=GFI_HXX; i<=GFI_HZZ;i++) {
      if (high->gridFuncs[i]->do_comparison) {
	DoComparison(high,med,low,i,
		     h[i-GFI_HXX],m[i-GFI_HXX],l[i-GFI_HXX]);
      }
    }
  }

  if (comp_gauge) {
    SetupFortranArrays(high);
    f_comparisongauge(FORT_ARGS(high), PASS_FORT_FIELDS(high),
		       h[0],h[1],h[2],h[3]);

    if (med) {
      SetupFortranArrays(med);
      f_comparisongauge(FORT_ARGS(med), PASS_FORT_FIELDS(med),
			m[0],m[1],m[2],m[3]);
    }

    if (low) {
      SetupFortranArrays(low);
      f_comparisongauge(FORT_ARGS(low), PASS_FORT_FIELDS(low),
			l[0],l[1],l[2],l[3]);
    }

    if (high->gridFuncs[GFI_ALP]->do_comparison) {
      DoComparison(high,med,low,GFI_ALP,h[0],m[0],l[0]);
    }

    if (high->gridFuncs[GFI_BETAX]->do_comparison) {
      DoComparison(high,med,low,GFI_BETAX,h[1],m[1],l[1]);
    }

    if (high->gridFuncs[GFI_BETAY]->do_comparison) {
      DoComparison(high,med,low,GFI_BETAY,h[2],m[2],l[2]);
    }

    if (high->gridFuncs[GFI_BETAZ]->do_comparison) {
      DoComparison(high,med,low,GFI_BETAZ,h[3],m[3],l[3]);
    }
  }

  for (i=0;i<6;i++) {
    free(h[i]); if (m[i]) free(m[i]); if (l[i]) free(l[i]);
  }
  DisableGFDataStorage(high, high->gridFuncs[high->ngridFuncs-1]);
  if (med)
    DisableGFDataStorage( med,  med->gridFuncs[med->ngridFuncs-1]);
  if (low)
    DisableGFDataStorage( low,  low->gridFuncs[low->ngridFuncs-1]);
  comptrips ++ ;
}

void DoComparison(pGH *high, pGH *med, pGH *low, int gfi,
		  double *he, double *me, double *le) {

  pGF *ws, *cur;
  double max[3], nm1[3], nm2[3];
  double sigtop, sigbot, sig3w, sig10, sig21;

  int i;

  printf ("Comparing %s\n",high->gridFuncs[gfi]->name);

  /* Output the analytic solution if we need it */
  cur = high->gridFuncs[gfi]; ws = high->gridFuncs[high->ngridFuncs-1];
  if (cur->do_1dio) {
    ws->do_1dio = comptrips;
    for (i=0;i<high->npoints;i++) 
      ws->data[i] = he[i];
    sprintf(ws->name,"%s_exact",cur->name);
    IO_Write1D(high,ws);
    ws->do_1dio = 0;
    ws->lastio_it[1]--;
  }

  /* Diffs go into workspace */
  for (i=0;i<high->npoints;i++) 
    ws->data[i] = cur->data[i] - he[i];

  if (cur->do_1dio) {
    ws->do_1dio = comptrips;
    sprintf(ws->name,"%s_diff",cur->name);
    IO_Write1D(high,ws);
    ws->do_1dio = 0;
    ws->lastio_it[1]--;
  }

  max[0] = pGF_MaxVal(high,ws);
  nm1[0] = pGF_Norm1 (high,ws);
  nm2[0] = pGF_Norm2 (high,ws);

  if (med) {
    cur = med->gridFuncs[gfi]; ws = med->gridFuncs[med->ngridFuncs-1];
    for (i=0;i<med->npoints;i++)
      ws->data[i] = cur->data[i] - me[i];
    if (cur->do_1dio) {
      ws->do_1dio = comptrips;
      sprintf(ws->name,"%s_diff",cur->name);
      IO_Write1D(med,ws);
      ws->do_1dio = 0;    
      ws->lastio_it[1]--;
    }

    max[1] = pGF_MaxVal(med,ws);
    nm1[1] = pGF_Norm1 (med,ws);
    nm2[1] = pGF_Norm2 (med,ws);
  }

  if (low) {
    cur = low->gridFuncs[gfi]; ws = low->gridFuncs[low->ngridFuncs-1];
    for (i=0;i<low->npoints;i++)
      ws->data[i] = cur->data[i] - le[i];
    if (cur->do_1dio) {
      ws->do_1dio = comptrips;
      sprintf(ws->name,"%s_diff",cur->name);
      IO_Write1D(low,ws);
      ws->do_1dio = 0;
      ws->lastio_it[1]--;
    }
    
    max[2] = pGF_MaxVal(low,ws);
    nm1[2] = pGF_Norm1 (low,ws);
    nm2[2] = pGF_Norm2 (low,ws);
  } else {
    max[2] = 0.0; nm1[2] = 0.0; nm2[2] = 0.0;
  }

  if (low) {
    printf (" ------ :   High     Med       Low\n");
    printf ("Max Diff: %lf %lf %lf\n",
	    max[0],max[1],max[2]);
    printf ("NM1 Diff: %lf %lf %lf\n",
	  nm1[0],nm1[1],nm1[2]);
    printf ("NM2 Diff: %lf %lf %lf\n",
	    nm2[0],nm2[1],nm2[2]);
  } else {
    if (med) {
      printf (" ------ :   High     Med      \n");
      printf ("Max Diff: %lf %lf\n",
	      max[0],max[1]);
      printf ("NM1 Diff: %lf %lf\n",
	      nm1[0],nm1[1]);
      printf ("NM2 Diff: %lf %lf\n",
	      nm2[0],nm2[1]);
    }
  }


  if (low) {
    sigbot = nm2[0] - nm2[1];
    if (sigbot == 0) sig3w = 0.0;
    else sig3w = log(fabs((nm2[1]-nm2[2])/(nm2[0]-nm2[1])))/log(2.0);

    if (nm2[1] == 0) sig21 = 0.0;
    else sig21 = log(fabs(nm2[2]/nm2[1]))/log(2.0);
  } else {
    sig3w = 0.0; sig21 = 0;
  }
    
  if (nm2[0] == 0) sig10 = 0.0;
  else sig10 = log(fabs(nm2[1]/nm2[0]))/log(2.0);

  if (low)
    printf ("sigma: 3way %lf  hm %lf  ml %lf\n",
	    sig3w, sig10, sig21);
  else
    if (med)
      printf ("sigma: hm %lf\n",sig10);
  
}
 
