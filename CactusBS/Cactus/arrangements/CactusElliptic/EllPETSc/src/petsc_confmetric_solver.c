/*@@
   @file      PETScEll_conformal.c
   @date      Wed Apr  9 09:31:24 1997
   @author    Paul Walker
   @desc
   A nabla phi - M phi = N elliptic solver based around PETSc.
   For information, see the documentation for @seeroutine petscEll_conformal
   @enddesc
@@*/


#include <stdlib.h>

#include <cctk.h>
#include <cctk_Parameters.h>

#include <petsc.h>

#include "ellpetsc.h"
#include "CactusPUGH/PUGH/src/include/pugh.h"

/*#define DEBUG*/

/*Don't know what these actually mean! FIXME */
#define ELLCONV_ABSOLUTE 0
#define ELLCONV_RELATIVE 1
#define ELLCONV_EITHER   2

/* direction macros */
#define XDP   1
#define XDM   0
#define YDP   3
#define YDM   2
#define ZDP   5
#define ZDM   4

/* A few convenient macros */
#define a(i,j,k) a[i+1 + 3*(j+1 + 3*(k+1))]
#define SQR(a) ((a)*(a))

  
/* Some useful definitions to deal with the upper metric */
#define Uxx(i,j,k) uxx3[DATINDEX(pEx,(i),(j),(k))]
#define Uxy(i,j,k) uxy3[DATINDEX(pEx,(i),(j),(k))]
#define Uxz(i,j,k) uxz3[DATINDEX(pEx,(i),(j),(k))]
#define Uyy(i,j,k) uyy3[DATINDEX(pEx,(i),(j),(k))]
#define Uyz(i,j,k) uyz3[DATINDEX(pEx,(i),(j),(k))]
#define Uzz(i,j,k) uzz3[DATINDEX(pEx,(i),(j),(k))]
  

/* Place the matrix in global space so we can keep it around
   with the memory stripped out. Thanks, Barry!  */

static int trips=0;
static Mat     *A;           /* linear system matrix */
static Vec     soln, b;      /* approx solution, RHS */
static KSP     ksp;          /* linear solver context */


int petsc_confmetric_solver(cGH *GH, int *MetricPsiI, int MetricPsiISize, 
                             int FieldIndex, int MIndex, int NIndex, 
                             CCTK_REAL *AbsTol, CCTK_REAL *RelTol);
void *GetDataPtr_NextTL(cGH *GH, const char *field);


void *GetDataPtr_NextTL(cGH *GH, const char *field) {
  int index = CCTK_VarIndex(field);
  if (index<0) {
    CCTK_VWarn(CCTK_WARN_ALERT, __LINE__,__FILE__,CCTK_THORNSTRING,
               "Index for >%s< not found",field);
  }
  return((CCTK_REAL*)CCTK_VarDataPtrI(GH,0,index));
}


/* This passing matches that of in the convention in the
   elliptic registration routine see LinearElliptic.h*/

int petsc_confmetric_solver(cGH *GH, int *MetricPsiI, int MetricPsiISize, 
                             int FieldIndex, int MIndex, int NIndex, 
                             CCTK_REAL *AbsTol, CCTK_REAL *RelTol) {

  DECLARE_CCTK_PARAMETERS    /* CCTK passed parameters */


  PC     pc;            /* preconditioner context */
  int    num_A;         /* Number of A-arrays as needed ny Dan's MG */
  
  int    ierr;          /* Check the return status of petsc */
  int    retcode;       /* Check the return status of CCTK */

  double  a[27];        /* The stencil array */
  int     rank;         /* Rank of the matrix/vector */
  PetscInt its;         /* Number of iterations */
  
  /* Loopers */
  int i,j,k,l,m,n;

  /* loop limits put in for stencil_w !=1     Ed Evans 1/19/98 */
  int imin,imax,jmin,jmax,kmin,kmax;
   
  pGH *pughGH;                 /* The pugh Extension handle */
  pGExtras *pEx;
  pConnectivity *pCon;
  int myproc;                  /* out processor */

  /* Tolerances */
  CCTK_REAL rtol=0, atolerance=0, tolerance;
    
  /* Values to assemble the matrix */
  CCTK_REAL two=2.0, four=4.0, tmp, det;
  CCTK_REAL pm4, psixp,psiyp,psizp;
  CCTK_REAL dx,dy,dz;
  CCTK_REAL uxx=0,uxy=0,uxz=0,uyy=0,uyz=0,uzz=0;
  CCTK_REAL Gxxx,Gxxy,Gxxz,Gxyy,Gxyz,Gxzz; /* Christoffels */
  CCTK_REAL Gyxx,Gyxy,Gyxz,Gyyy,Gyyz,Gyzz;
  CCTK_REAL Gzxx,Gzxy,Gzxz,Gzyy,Gzyz,Gzzz;
  CCTK_REAL dxxx,dxxy,dxxz,dxyy,dxyz,dxzz;
  CCTK_REAL dyxx,dyxy,dyxz,dyyy,dyyz,dyzz;
  CCTK_REAL dzxx,dzxy,dzxz,dzyy,dzyz,dzzz;
  CCTK_REAL *values;
  
  int nxs,nys,nzs;      /* Size of the grid with stencils off... */
  int startpoint=0;     /* My starting index (per proc) */
  int endpoint;         /* One more than my end */
  PetscInt pstart, pend;     /* A check for PETSc layout */
  PetscInt pvstart, pvend;   /* A check for PETSc layout */
  int verbose;          /* Is the solver verbose */
  int debug;            /* Is the solver debug-verbose */
  int octant;           /* Apply octant BCs inside */
  int conformal=0;        /* Do we have conformal metric ? */
  int nabla_form=0;       /* Which form of the nable */
  int matnormalize;     /* Normalize the central mat value to one? */
  CCTK_REAL ac;         /* Storage for a(0,0,0) for renorm */
  int PetscTolStyle;

  /* For the upper metric form of nabla */
  CCTK_REAL *uxx3=NULL, *uxy3=NULL, *uxz3=NULL, *uyy3=NULL, *uyz3=NULL, *uzz3=NULL;

  
  /* Pointers to the data of : petsc workspace/ the GF to solve /
     metric / psi / derivs(psi) / Mlinear / Nlinear(source) */
  CCTK_REAL *wsp =NULL, *ell_field=NULL; 
  CCTK_REAL *gxx =NULL, *gxy =NULL; 
  CCTK_REAL *gxz =NULL, *gyy =NULL; 
  CCTK_REAL *gyz =NULL, *gzz =NULL; 
  CCTK_REAL *Mlin=NULL, *Nlin=NULL;   
  CCTK_REAL *psi =NULL;
  CCTK_REAL *psix=NULL, *psiy=NULL, *psiz=NULL; 

  int Mstorage=0, Nstorage=0;


  RelTol = RelTol;

  octant       = CCTK_Equals(domain,"octant");
  verbose      = CCTK_Equals(petsc_verbose,"yes")||
                 CCTK_Equals(petsc_verbose,"debug");
  debug        = CCTK_Equals(petsc_verbose,"debug");
  matnormalize = 0; 
  
  if (MetricPsiISize==7) conformal=1;
  else if (MetricPsiISize==7) conformal=0;
  else CCTK_WARN(0,"Size of the Metric must be either 7 (metric+conformal) or 6 (metric)");


  /* FIXME, the TolAbs/TolRel will be evaluated here */
  PetscTolStyle=0;
  tolerance    =AbsTol[0];

  /* Get the link to pugh Extension */
  pughGH = PUGH_pGH(GH);
  if (!pughGH) CCTK_WARN(0,"ETERNAL ERROR: Cannot find PUGH Extension Handle\n");

  /* Get the extras extension for 3D grid functions */
  pEx = pughGH->GFExtras[2];
  pCon = pughGH->Connectivity[2];

  /* Things to do on first iteration */
  if (trips==0) {
    int argc;
    char **argv; 

    if (debug) {
#ifdef DEBUG
      printf("PETSc: initial trip: %d \n",trips);
#endif
    }



    /* Get the commandline arguments */
    argc = CCTK_CommandLine(&argv);

    /* Set the PETSc communicator to set of PUGH and initialzie
       PETSc */
    PETSC_COMM_WORLD = pughGH->PUGH_COMM_WORLD;
    PetscInitialize(&argc,&argv,NULL,NULL); 

    CCTK_INFO("PETSc initialized");
  }

  trips++;

  /* Create a array of matrices A */
  /* num_A = MultiGridCount(); */
  num_A = 1;
  A=(Mat*) malloc (sizeof(Mat)*num_A);
  
  /* Get the data ptr of these GFs, They all have to be
     on the same timelevel */

  /* derive the metric data pointer from the index array. 
     Note the ordering in the metric */
  ell_field = (CCTK_REAL*) CCTK_VarDataPtrI(GH,0,FieldIndex);
  gxx = (CCTK_REAL*) CCTK_VarDataPtrI(GH, 0, MetricPsiI[0]);
  gxy = (CCTK_REAL*) CCTK_VarDataPtrI(GH, 0, MetricPsiI[1]);
  gxz = (CCTK_REAL*) CCTK_VarDataPtrI(GH, 0, MetricPsiI[2]);
  gyy = (CCTK_REAL*) CCTK_VarDataPtrI(GH, 0, MetricPsiI[3]);
  gyz = (CCTK_REAL*) CCTK_VarDataPtrI(GH, 0, MetricPsiI[4]);
  gzz = (CCTK_REAL*) CCTK_VarDataPtrI(GH, 0, MetricPsiI[5]);

  /* FIXME: evolved derivatives should go! -> Einstein specific*/  
  if (conformal) {
    psi = (CCTK_REAL*) CCTK_VarDataPtrI(GH, 0, MetricPsiI[6]);
    psix =(CCTK_REAL*) GetDataPtr_NextTL(GH,"einstein::psix");
    psiy =(CCTK_REAL*) GetDataPtr_NextTL(GH,"einstein::psiy");
    psiz =(CCTK_REAL*) GetDataPtr_NextTL(GH,"einstein::psiz");
  }
  /* if we have a negative index, this GF is not needed, 
     there for don't even look for it. when index positive,
     set flag Mstorage=1, dito for N */
  if (MIndex>=0)  { 
    Mlin = (CCTK_REAL*) CCTK_VarDataPtrI(GH,0,MIndex);
    Mstorage = 1;
  }
  if (NIndex>=0) {
    Nlin = (CCTK_REAL*) CCTK_VarDataPtrI(GH,0,NIndex);
    Nstorage = 1;
  }

  /* Get the workspace data pointer */
  wsp  = GetDataPtr_NextTL(GH,"ellpetsc::wsp");
  
  /* initialize the linear index lookup table, after it's 
     filled up (below), the -1 indicates a boundary */
  for (i=0;i<pEx->npoints;i++) wsp[i] = -1.0;

  /* Get myproc from the CCTK, get the gridspacing */
  myproc = PUGH_MyProc(GH);
  dx     = GH->cctk_delta_space[0];
  dy     = GH->cctk_delta_space[1];
  dz     = GH->cctk_delta_space[2];

  /* We fix the lower and upper boundary indices, that actually "active" in 
     the sense that they are no ghostzones: */
  imin=((pCon->neighbours[pughGH->myproc][XDM]<0) && !(octant) ? 1 : 
         GH->cctk_nghostzones[0]);
  imax=((pCon->neighbours[pughGH->myproc][XDP]<0) ? pEx->lnsize[0]-1 : 
         pEx->lnsize[0]-GH->cctk_nghostzones[0]);

  jmin=((pCon->neighbours[pughGH->myproc][YDM]<0) && !(octant) ? 1 : 
        GH->cctk_nghostzones[1]);
  jmax=((pCon->neighbours[pughGH->myproc][YDP]<0) ? pEx->lnsize[1]-1 : 
         pEx->lnsize[1]-GH->cctk_nghostzones[1]);
  
  kmin=((pCon->neighbours[pughGH->myproc][ZDM]<0) && !(octant) ? 1 : 
         GH->cctk_nghostzones[2]);
  kmax=((pCon->neighbours[pughGH->myproc][ZDP]<0) ? pEx->lnsize[2]-1 : 
         pEx->lnsize[2]-GH->cctk_nghostzones[2]);

  /* We need to get the global index of gridpoints (gps) owned on my proc. 
     For that we have to know how many gps are there before my processors
     (i<myproc), so we count them here (that's minus the ghostszones
     (we use PUGH's "rn" variables (number of gp on proc[#] in 
     direction [i]) and subtract the ghostzones if nec.) */

  for (i=0;i<myproc;i++) {
    
    nxs=pEx->rnsize[i][0];
    nxs=((pCon->neighbours[i][XDM]<0) && !(octant) ? nxs-1 : 
         nxs-GH->cctk_nghostzones[0]); 
    nxs=((pCon->neighbours[i][XDP]<0) ? nxs-1 : nxs-GH->cctk_nghostzones[0]);
    
    nys=pEx->rnsize[i][1];
    nys=((pCon->neighbours[i][YDM]<0) && !(octant) ? nys-1 : 
         nys-GH->cctk_nghostzones[1]);
    nys=((pCon->neighbours[i][YDP]<0) ? nys-1 : nys-GH->cctk_nghostzones[1]);
    
    nzs=pEx->rnsize[i][2];
    nzs=((pCon->neighbours[i][ZDM]<0) && !(octant) ? nzs-1 : 
         nzs-GH->cctk_nghostzones[2]);
    nzs=((pCon->neighbours[i][ZDP]<0) ? nzs-1 : nzs-GH->cctk_nghostzones[2]);    
    
    printf("PROC: %d  nxyzs %d %d %d \n",i,nxs,nys,nzs);
    startpoint += nxs*nys*nzs;
  }

#ifdef DEBUG 
  printf("STARTPOINT: %d \n",startpoint);
#endif 

  /* ...we save that as a start ...*/
  endpoint = startpoint;
    
  /* ...and continue running over our own region to get our endpoint */
  for (k=kmin;k<kmax;k++)
    for (j=jmin;j<jmax;j++)
      for (i=imin;i<imax;i++)
        wsp[DATINDEX(pEx,i,j,k)] = endpoint++;
  
  /* So now each point has a unique index of its own. If we
     do a sync that means each processor knows the indiced
     in the ghost zones (eg, on its neighbors) in a FD stencil
     of 1 */
  retcode = CCTK_SyncGroup(GH,"ellpetsc::petscworkspace");
  if (retcode<0) CCTK_WARN(1,"Synchronization failed\n");
  
  
  /* So woohoo. Now for each point in our ijk space, we have
     information about our row in the matrix (workspace->data[ijk])
     as well as the column for the stencil of the neighbors
     (workspace->data[DATINDEX(GH,i+1,j,k)] etc...). So onwards */

  nxs = pEx->nsize[0]-2;
  nys = pEx->nsize[1]-2;
  nzs = pEx->nsize[2]-2;

#ifdef DEBUG
  if (debug) {
    printf("nxyzs %d %d %d \n",nxs,nys,nzs);
    printf("lnxyz: %d %d %d \n",pEx->lnsize[0],pEx->lnsize[1],pEx->lnsize[2]);
  }
#endif

  /* Suck off boundaries... */
  rank = nxs*nys*nzs;
  
  
  /* Create the petsc matrix A and vectors x and b efficiently.
     
     Now for efficiency we have to be really careful here. By
     setting up the index array, we are abandoning the 19-way
     tri-diagonal system so we can lay out our matrices on the
     processors with the data. (One one proc we will still
     get 19 way tridiagonal). So we want each processor
     to own only part of the matrix/vector but we know exactly
     what it is, eg, from startpoint to endpoint-1 is on our
     local processor (that is, row-wise ; all columns are local).

     So we can use MatCreateMPIAIJ for exactly this, and
     VecCreateMPI for the vectors.
   */

  /* FIXME: perhaps only set up only on first iteration, this was a hack
     by PAUL, using some inofficical petsc code, check if this is 
     "official", yet */

  if (verbose) CCTK_INFO("Creating Matrices and Vectors ....");

  ierr = MatCreateAIJ(pughGH->PUGH_COMM_WORLD,     
                      (endpoint-startpoint),   /* # of rows */
                      (endpoint-startpoint),   /* This is confusing */
                      rank,rank,               /* Matrix size */
                      19,PETSC_NULL,           /* Diagonal  */
                      19,PETSC_NULL,           /* Off diagonal */
                      &A[0]);                  /* The output  */
  CHKERRQ(ierr);
  ierr = VecCreateMPI(pughGH->PUGH_COMM_WORLD,(endpoint-startpoint),rank,&soln);
  CHKERRQ(ierr);
  ierr = VecCreateMPI(pughGH->PUGH_COMM_WORLD,(endpoint-startpoint),rank,&b);
  CHKERRQ(ierr);
 
  /* Compare the PETSc layout to Cactus, this better be a match */
  ierr = VecGetOwnershipRange(soln,&pvstart,&pvend);
  ierr = MatGetOwnershipRange(A[0],&pstart,&pend);
#ifdef DEBUG
  printf("CAC M-Layout: %d %d \n",startpoint, endpoint);
  printf("PET M-Layout: %d %d \n",pstart, pend);
  printf("PET V-Layout: %d %d \n",pvstart,pvend);
#endif
  if (pstart != startpoint && pend != endpoint) 
    CCTK_WARN(1,"WARNING: PETSC and data layouts differ! (why?)\n");

  /* Decide on the nabla form in PETSc: */
  if (CCTK_EQUALS(petsc_nablaform,"down")) {
    if (verbose) 
      CCTK_INFO("Forming nabla with lower g and finite difference of dg \n");
    nabla_form = 2;

  } 
  else 
    if  (CCTK_EQUALS(petsc_nablaform,"up")) {
      if (verbose) 
        CCTK_INFO("Forming nabla with upper g and finite difference of dg \n");
      nabla_form = 3;
      uxx3 = (CCTK_REAL*)malloc(pEx->npoints*sizeof(CCTK_REAL));
      uxy3 = (CCTK_REAL*)malloc(pEx->npoints*sizeof(CCTK_REAL));
      uxz3 = (CCTK_REAL*)malloc(pEx->npoints*sizeof(CCTK_REAL));
      uyy3 = (CCTK_REAL*)malloc(pEx->npoints*sizeof(CCTK_REAL));
      uyz3 = (CCTK_REAL*)malloc(pEx->npoints*sizeof(CCTK_REAL));
      uzz3 = (CCTK_REAL*)malloc(pEx->npoints*sizeof(CCTK_REAL));
      for (i=0;i<pEx->npoints;i++) {
        CCTK_REAL p12;
        det= -(SQR(gxz[i])*gyy[i]) + 
          2*gxy[i]*gxz[i]*gyz[i] - 
          gxx[i]*SQR(gyz[i])  -
          SQR(gxy[i])*gzz[i] + 
          gxx[i]*gyy[i]*gzz[i];
        
        if (conformal) {
          pm4 = 1./pow(psi[i],4.0);
          p12 = pow(psi[i],12.0);
        } else {
          pm4 = 1.0;
          p12 = 1.0;
        }
      
        /*invert metric. This is the conformal upper metric. */
        uxx3[i]=(-SQR(gyz[i]) + gyy[i]*gzz[i])/det*pm4;
        uxy3[i]=(gxz[i]*gyz[i] - gxy[i]*gzz[i])/det*pm4;
        uyy3[i]=(-SQR(gxz[i]) + gxx[i]*gzz[i])/det*pm4;
        uxz3[i]=(-gxz[i]*gyy[i] + gxy[i]*gyz[i])/det*pm4;
        uyz3[i]=(gxy[i]*gxz[i] - gxx[i]*gyz[i])/det*pm4;
        uzz3[i]=(-SQR(gxy[i]) + gxx[i]*gyy[i])/det*pm4;
        
        det = det*p12;
      
        /* Rescaling for the uppermetric solver */
        if (Mstorage) Mlin[i] = Mlin[i]*sqrt(det);
        if (Nstorage) Nlin[i] = Nlin[i]*sqrt(det);

        uxx3[i]=uxx3[i]/(2.*dx*dx)*sqrt(det);
        uyy3[i]=uyy3[i]/(2.*dy*dy)*sqrt(det);
        uzz3[i]=uzz3[i]/(2.*dz*dz)*sqrt(det);
        uxy3[i]=uxy3[i]/(4.*dx*dy)*sqrt(det);
        uxz3[i]=uxz3[i]/(4.*dx*dz)*sqrt(det);
        uyz3[i]=uyz3[i]/(4.*dy*dz)*sqrt(det);
      }
    } 
    else  CCTK_WARN(0,"Don't know how to form nabla!\n");
  
  if (verbose) 
    CCTK_INFO("Creating the coefficient matrix...");
  if (verbose && matnormalize) 
    CCTK_INFO("...with diagonal renormalized to one");

  for (k=kmin;k<kmax;k++) {
    for (j=jmin;j<jmax;j++) {
      for (i=imin;i<imax;i++) {

        if (wsp[DATINDEX(pEx,i,j,k)] >= 0) {

          CCTK_REAL tdxgxx, tdxgxy, tdxgxz, tdxgyy, tdxgyz, tdxgzz;
          CCTK_REAL tdygxx, tdygxy, tdygxz, tdygyy, tdygyz, tdygzz;
          CCTK_REAL tdzgxx, tdzgxy, tdzgxz, tdzgyy, tdzgyz, tdzgzz;

          /* Set up indices */
          int ijk;              /* The data point for the array */
          int ig, jg, kg;       /* The position in global space */
          PetscInt I,J;         /* The col and row in the matrix */
          CCTK_REAL rhsval = 0;
          
          for (I=0;I<27;I++) a[I] = 0.0;
          
          /* these guys are easy */
          ijk = DATINDEX(pEx,i,j,k);               /* get linear index */
          ig  = i + GH->cctk_lbnd[0]; 
          jg  = j + GH->cctk_lbnd[1];
          kg  = k + GH->cctk_lbnd[2];
        
          /* Get the row we are working on */
          I = wsp[ijk];

          /* Setup Temporaries / Psi derivatives on psi */
          if (conformal) {
            pm4   = 1./pow(psi[ijk],4.0);
            psixp = psix[ijk];
            psiyp = psiy[ijk];
            psizp = psiz[ijk];
          } else {
            pm4   = 1.0;
            psixp = 0.0;
            psiyp = 0.0;
            psizp = 0.0;
          }
          if (nabla_form == 2) {
            /* Use finite differences of g for the d's */
            int ijkp, ijkm;
            
            /* X derivatives */
            ijkp = DATINDEX(pEx,i+1,j,k);
            ijkm = DATINDEX(pEx,i-1,j,k);
            
            tdxgxx = (gxx[ijkp] - gxx[ijkm])/(2.0*dx);
            tdxgxy = (gxy[ijkp] - gxy[ijkm])/(2.0*dx);
            tdxgxz = (gxz[ijkp] - gxz[ijkm])/(2.0*dx);
            tdxgyy = (gyy[ijkp] - gyy[ijkm])/(2.0*dx);
            tdxgyz = (gyz[ijkp] - gyz[ijkm])/(2.0*dx);
            tdxgzz = (gzz[ijkp] - gzz[ijkm])/(2.0*dx);
              
              
            /* Y derivatives */
            ijkp = DATINDEX(pEx,i,j+1,k);
            ijkm = DATINDEX(pEx,i,j-1,k);
            
            tdygxx = (gxx[ijkp] - gxx[ijkm])/(2.0*dy);
            tdygxy = (gxy[ijkp] - gxy[ijkm])/(2.0*dy);
            tdygxz = (gxz[ijkp] - gxz[ijkm])/(2.0*dy);
            tdygyy = (gyy[ijkp] - gyy[ijkm])/(2.0*dy);
            tdygyz = (gyz[ijkp] - gyz[ijkm])/(2.0*dy);
            tdygzz = (gzz[ijkp] - gzz[ijkm])/(2.0*dy);
            
            /* X derivatives */
            ijkp = DATINDEX(pEx,i,j,k+1);
            ijkm = DATINDEX(pEx,i,j,k-1);
            
            tdzgxx = (gxx[ijkp] - gxx[ijkm])/(2.0*dz);
            tdzgxy = (gxy[ijkp] - gxy[ijkm])/(2.0*dz);
            tdzgxz = (gxz[ijkp] - gxz[ijkm])/(2.0*dz);
            tdzgyy = (gyy[ijkp] - gyy[ijkm])/(2.0*dz);
            tdzgyz = (gyz[ijkp] - gyz[ijkm])/(2.0*dz);
            tdzgzz = (gzz[ijkp] - gzz[ijkm])/(2.0*dz);
            
            /* great ... so start hacking away at the coefficients.
               Form upper metric - compute determinant */
            det= -(SQR(gxz[ijk])*gyy[ijk]) + 
              2*gxy[ijk]*gxz[ijk]*gyz[ijk] - 
              gxx[ijk]*SQR(gyz[ijk])  -
              SQR(gxy[ijk])*gzz[ijk] + 
              gxx[ijk]*gyy[ijk]*gzz[ijk];
            
            /*invert metric. This is the conformal upper metric. */
            uxx=(-SQR(gyz[ijk]) + gyy[ijk]*gzz[ijk])/det;
            uxy=(gxz[ijk]*gyz[ijk] - gxy[ijk]*gzz[ijk])/det;
            uyy=(-SQR(gxz[ijk]) + gxx[ijk]*gzz[ijk])/det;
            uxz=(-gxz[ijk]*gyy[ijk] + gxy[ijk]*gyz[ijk])/det;
            uyz=(gxy[ijk]*gxz[ijk] - gxx[ijk]*gyz[ijk])/det;
            uzz=(-SQR(gxy[ijk]) + gxx[ijk]*gyy[ijk])/det;
            
            /* Coeff. Contributions from second derivative */
            
            /* X derivative */
            a(-1,0,0) =         pm4 * uxx / (dx*dx);
            a(0,0,0)  =  -two * pm4 * uxx / (dx*dx);
            a(1,0,0)  =         pm4 * uxx / (dx*dx);
            /* Y derivative */
            a(0,-1,0) =         pm4 * uyy / (dy*dy);
            a(0,0,0)  = a(0,0,0)  -  two * pm4 * uyy / (dy*dy);
            a(0,1,0)  =         pm4 * uyy / (dy*dy);
            /* Z derivative */
            a(0,0,-1) =         pm4 * uzz / (dz*dz);
            a(0,0,0)  = a(0,0,0)  -  two * pm4 * uzz / (dz*dz);
            a(0,0,1)  =         pm4 * uzz / (dz*dz);
            /* Mixed XY */
            a(1,1,0)  =   two * pm4 * uxy / (dx*dy);
            a(-1,-1,0)=   two * pm4 * uxy / (dx*dy);
            a(1,-1,0) =  -two * pm4 * uxy / (dx*dy);
            a(-1,1,0) =  -two * pm4 * uxy / (dx*dy);
            /* Mixed XZ */
            a(1,0,1)  =   two * pm4 * uxz / (dx*dz);
            a(-1,0,-1)=   two * pm4 * uxz / (dx*dz);
            a(1,0,-1) =  -two * pm4 * uxz / (dx*dz);
            a(-1,0,1) =  -two * pm4 * uxz / (dx*dz);
            /* Mixed YZ */
            a(0,1,1)  =   two * pm4 * uyz / (dz*dy);
            a(0,-1,-1)=   two * pm4 * uyz / (dz*dy);
            a(0,-1,1) =  -two * pm4 * uyz / (dz*dy);
            a(0,1,-1) =  -two * pm4 * uyz / (dz*dy);
            
            /*     Great so now form christoffels. Remember that
                   
                   G_kij = psi^4 (2 * (psi_j/psi g_ik + psi_i/psi g_jk - 
                   psi_k/psi g_ij)
                   + D_jik + D_ijk - D_kij)
                   
                   Since these are the DOWN christoffels, store them in d...
                   
                   NOTE however that since we will up these, we can drop the 
                   psi^4 here and psi^-4 from the up metric. */
            
            
            /*     These three have lots of cancelations */
            dxxx = (two * psixp * gxx[ijk] + tdxgxx);
            dxxy = (two * psiyp * gxx[ijk] + tdygxx);
            dxxz = (two * psizp * gxx[ijk] + tdzgxx);
            /* This one has a reduction of two identical terms */
            dxyy = (four * psiyp * gxy[ijk] -
                    two  * psixp * gyy[ijk] +
                    two  * tdygxy - tdxgyy);
            /* As does this one */
            dxzz = (four * psizp * gxz[ijk] -
                    two  * psixp * gzz[ijk] +
                    two  * tdzgxz - tdxgzz);
            /* And this one is completely general */
            dxyz = (two * psiyp * gxz[ijk] +
                    two * psizp * gxy[ijk] -
                    two * psixp * gyz[ijk] +
                    tdzgxy + tdygxz - tdxgyz);
            
            /* Now do it twice more without the explanations */
            dyyy = (two * psiyp * gyy[ijk] + tdygyy);
            dyxy = (two * psixp * gyy[ijk] + tdxgyy);
            dyyz = (two * psizp * gyy[ijk] + tdzgyy);
            dyxx = (four * psixp * gxy[ijk] -
                    two  * psiyp * gxx[ijk] +
                    two  * tdxgxy - tdygxx);
            dyzz = (four * psizp * gyz[ijk] -
                    two  * psiyp * gzz[ijk] +
                    two  * tdzgyz - tdygzz);
            dyxz = (two * psizp * gxy[ijk] +
                    two * psixp * gyz[ijk] -
                    two * psiyp * gxz[ijk] +
                    tdzgxy + tdxgyz - tdygxz);
            
            dzzz = (two * psizp * gzz[ijk] + tdzgzz);
            dzxz = (two * psixp * gzz[ijk] + tdxgzz);
            dzyz = (two * psiyp * gzz[ijk] + tdygzz);
            dzxx = (four * psixp * gxz[ijk] -
                    two  * psizp * gxx[ijk] +
                    two  * tdxgxz - tdzgxx);
            dzyy = (four * psiyp * gyz[ijk] -
                    two  * psizp * gyy[ijk] +
                    two  * tdygyz - tdzgyy);
            dzxy = (two * psiyp * gxz[ijk] +
                    two * psixp * gyz[ijk] -
                    two * psizp * gxy[ijk] +
                    tdxgyz + tdygxz - tdzgxy);
            
            /* And now raise the first index */
            Gxxx = uxx*dxxx + uxy*dyxx + uxz*dzxx;
            Gxxy = uxx*dxxy + uxy*dyxy + uxz*dzxy;
            Gxxz = uxx*dxxz + uxy*dyxz + uxz*dzxz;
            Gxyy = uxx*dxyy + uxy*dyyy + uxz*dzyy;
            Gxyz = uxx*dxyz + uxy*dyyz + uxz*dzyz;
            Gxzz = uxx*dxzz + uxy*dyzz + uxz*dzzz;
            
            Gyxx = uxy*dxxx + uyy*dyxx + uyz*dzxx;
            Gyxy = uxy*dxxy + uyy*dyxy + uyz*dzxy;
            Gyxz = uxy*dxxz + uyy*dyxz + uyz*dzxz;
            Gyyy = uxy*dxyy + uyy*dyyy + uyz*dzyy;
            Gyyz = uxy*dxyz + uyy*dyyz + uyz*dzyz;
            Gyzz = uxy*dxzz + uyy*dyzz + uyz*dzzz;
            
            Gzxx = uxz*dxxx + uyz*dyxx + uzz*dzxx;
            Gzxy = uxz*dxxy + uyz*dyxy + uzz*dzxy;
            Gzxz = uxz*dxxz + uyz*dyxz + uzz*dzxz;
            Gzyy = uxz*dxyy + uyz*dyyy + uzz*dzyy;
            Gzyz = uxz*dxyz + uyz*dyyz + uzz*dzyz;
            Gzzz = uxz*dxzz + uyz*dyzz + uzz*dzzz;
               

            /*     Great. So now start adding the summed contributions
                   from the first derivative. Note that these all have a
                   sign change since the term comes in - ... */
            
            /*  g^ij G ^x_ij */
            tmp = uxx * Gxxx + uyy * Gxyy + uzz * Gxzz +
              two * uxy * Gxxy +
              two * uxz * Gxxz + 
              two * uyz * Gxyz;
            a(1,0,0)  = a(1,0,0)  - pm4 * tmp / (two*dx);
            a(-1,0,0) = a(-1,0,0) + pm4 * tmp / (two*dx);
            
            /* g^ij G^y_ij */
            tmp = uxx * Gyxx + uyy * Gyyy + uzz * Gyzz +
              two * uxy * Gyxy +
              two * uxz * Gyxz + 
              two * uyz * Gyyz;
            a(0,1,0)  = a(0,1,0)  - pm4 * tmp / (two*dy);
            a(0,-1,0) = a(0,-1,0) + pm4 * tmp / (two*dy);
            
            /* g^ij G^z_ij */
            tmp = uxx * Gzxx + uyy * Gzyy + uzz * Gzzz +
              two * uxy * Gzxy +
              two * uxz * Gzxz + 
              two * uyz * Gzyz;
            a(0,0,1)  = a(0,0,1)  - pm4 * tmp / (two*dz);
            a(0,0,-1) = a(0,0,-1) + pm4 * tmp / (two*dz);
              
          } /* end if nable_form==2 */ 

          else { 
            /* nabla_form==3: Upper Metric Nabla Form */
            a(0,0,0) = -Uxx(i+1,j,k) -2.*Uxx(i,j,k) -Uxx(i-1,j,k) 
              -Uyy(i,j+1,k) -2.*Uyy(i,j,k) -Uyy(i,j-1,k)
              -Uzz(i,j,k+1) -2.*Uzz(i,j,k) -Uzz(i,j,k-1);
            
            /*$ae = uxx(i+1,j,k)+uxx(i,j,k)$*/
            a(1,0,0) = Uxx(i+1,j,k) + Uxx(i,j,k);
            /*$aw = uxx(i-1,j,k)+uxx(i,j,k)$*/
            a(-1,0,0) = Uxx(i-1,j,k)+Uxx(i,j,k);
                /*$an = uyy(i,j+1,k)+uyy(i,j,k)$*/              
            a(0,1,0) = Uyy(i,j+1,k)+Uyy(i,j,k);
            /*$as = uyy(i,j-1,k)+uyy(i,j,k)$*/              
            a(0,-1,0) = Uyy(i,j-1,k)+Uyy(i,j,k);
            /*$at = uzz(i,j,k+1)+uzz(i,j,k)$*/               
            a(0,0,1) = Uzz(i,j,k+1)+Uzz(i,j,k);               
            /*$ab = uzz(i,j,k-1)+uzz(i,j,k)$*/               
            a(0,0,-1) = Uzz(i,j,k-1)+Uzz(i,j,k);
            
            /*$ane = uxy(i,j+1,k)+uxy(i+1,j,k)$*/              
            a(1,1,0) = Uxy(i,j+1,k)+Uxy(i+1,j,k);
            /*$anw = - uxy(i-1,j,k)-uxy(i,j+1,k)$*/               
            a(-1,1,0) = - Uxy(i-1,j,k)-Uxy(i,j+1,k);
            /*$ase = - uxy(i+1,j,k)-uxy(i,j-1,k)$*/               
            a(1,-1,0) = - Uxy(i+1,j,k)-Uxy(i,j-1,k);
            /*$asw = uxy(i-1,j,k)+uxy(i,j-1,k)$*/               
            a(-1,-1,0) = Uxy(i-1,j,k)+Uxy(i,j-1,k);
            /*$ate = uxz(i,j,k+1)+uxz(i+1,j,k)$*/               
            a(1,0,1) = Uxz(i,j,k+1)+Uxz(i+1,j,k);
            /*$atw = - uxz(i-1,j,k)-uxz(i,j,k+1)$*/               
            a(-1,0,1) = - Uxz(i-1,j,k)-Uxz(i,j,k+1);
            /*$abe = - uxz(i+1,j,k)-uxz(i,j,k-1)$*/               
            a(1,0,-1) = - Uxz(i+1,j,k)-Uxz(i,j,k-1);
            /*$abw = uxz(i-1,j,k)+uxz(i,j,k-1)$*/               
            a(-1,0,-1) = Uxz(i-1,j,k)+Uxz(i,j,k-1);         
            /*$atn = uyz(i,j+1,k)+uyz(i,j,k+1)$*/               
            a(0,1,1) = Uyz(i,j+1,k)+Uyz(i,j,k+1);
            /*$ats = - uyz(i,j,k+1)-uyz(i,j-1,k)$*/               
            a(0,-1,1) = - Uyz(i,j,k+1)-Uyz(i,j-1,k);
            /*$abn = - uyz(i,j,k-1)-uyz(i,j+1,k)$*/               
            a(0,1,-1) = - Uyz(i,j,k-1)-Uyz(i,j+1,k);
            /*$asb = uyz(i,j,k-1)+uyz(i,j-1,k)$*/
            a(0,-1,-1) = Uyz(i,j,k-1)+Uyz(i,j-1,k);
          } /* end nabla_form=3 */

          /* M phi */
          if (Mstorage)
            a(0,0,0) = a(0,0,0) + Mlin[ijk];

          /* Great now set the values of the matrix. This is
             really painful due to the boundaries (here we force
             dirichlet).  */
          VecSetValues(soln,1,&I,&(ell_field[ijk]),INSERT_VALUES);
          
          /* Put in the octant boundary conditions.
           
             We do these by reflecting the -1 stencil if we are at a
             boundary at -1. That is, imagine the 1D laplace equation
             with the boundary condition a(0) = a(1). This means
             the point 1 stencil (which is our first stencil in
             the matrix) becomes
             
             a(0) + a(2) - 2 a(1) = rho(i) * deltax  ->
             a(2) - 2 a(1) + a(1) = rho(i) * deltax
             
             OK so think about this with a general stenci
             
             S_0 a_0 + S_1 a_1 + S_2 a_2 = rho_1 ->
             S_0 a_1 + S_1 a_1 + S_2 a_2 = rho_1 ->
             (S_0 + S_1) a_1 + S_2 a_2 = rho_1
             
             eg, S_0 -> 0 and S_1 -> S_1 + S_0
             
             Great, so implement this in 3D. It is a bit trickier,
             since we don't want to zero neighbors in the wrong
             direction, but not impossible.
             
            */
#ifdef DEBUG
          if (ijk==1918) {
            for (l=-1;l<=1;l++)
              for (m=-1;m<=1;m++)
                for (n=-1;n<=1;n++) printf(" (%d %d %d): %f \n",l,m,n,a(l,m,n));
            printf("\n");
          }
#endif

          if (octant)
            for (l=-1;l<=0;l++)
              for (m=-1;m<=0;m++)
                for (n=-1;n<=0;n++)
                  if (l*m*n == 0) 
                    if (wsp[DATINDEX(pEx,i+l,j+m,k+n)] < 0 &&
                        (ig == imin || jg == jmin || kg == kmin))  
                      {
                        /* We are on an inner boundary point, eg, 
                           at an inner face */
                        int ll,mm,nn;
                        /* Only zero the guys at the boundaries */
                        ll = (ig == imin ? 0 : l);
                        mm = (jg == jmin ? 0 : m);
                        nn = (kg == kmin ? 0 : n);
                        a(ll,mm,nn) += a(l,m,n);
                        a(l,m,n)     = 0.0;
                      }
      
          /* renormalize */
          if (matnormalize) {
            ac = a(0,0,0);
            for (J=0;J<27;J++) a[J] = a[J] / ac;
          } 
          else ac = 1;
          

          /* This is the new way-clever look. Note it relies
             heavily on the index array in the workspace and
             on multiple processors it will *NOT* make a
             19-way banded matrix. */

          for (l=-1;l<=1;l++)
            for (m=-1;m<=1;m++)
              for (n=-1;n<=1;n++) {
                if (l*m*n == 0) {       /* This is the 19 point ... if none are
                                         * zero, then we have no stencil here.
                                         */
                  if (wsp[DATINDEX(pEx,i+l,j+m,k+n)] < 0) {
                    /* This is a boundary. */
                    rhsval += a(l,m,n) * ell_field[DATINDEX(pEx,i+l,j+m,k+n)];
                  } 
                  else {
                    J = wsp[DATINDEX(pEx,i+l,j+m,k+n)];
                    ierr = MatSetValues(A[0],1,&I,1,&J,&(a(l,m,n)),INSERT_VALUES); 
                    CHKERRQ(ierr);
                  }
                }
              }

          if (Nstorage)
            rhsval = -rhsval - Nlin[ijk] / ac;

          ierr   = VecSetValues(b,1,&I,&rhsval,INSERT_VALUES);
          CHKERRQ(ierr);
          
        }
      }
    }
  }
  
  if (verbose) CCTK_INFO("Assembling the vectors");
  ierr = MatAssemblyBegin(A[0],MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(soln);                    CHKERRQ(ierr);
  ierr = VecAssemblyBegin(b);                       CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A[0],MAT_FINAL_ASSEMBLY);   CHKERRQ(ierr);
  ierr = VecAssemblyEnd(soln);                      CHKERRQ(ierr);
  ierr = VecAssemblyEnd(b);                         CHKERRQ(ierr);
#if PETSC_VERSION_MAJOR < 3
  ierr = MatSetOption(A[0],MAT_NO_NEW_NONZERO_LOCATIONS); CHKERRQ(ierr);
#else
  ierr = MatSetOption(A[0],MAT_NEW_NONZERO_LOCATIONS,PETSC_FALSE); CHKERRQ(ierr);
#endif

  if (nabla_form == 3) {
    if (verbose)
        printf("Freeing upper metric storage\n");
    free(uxx3);
    free(uxy3);
    free(uxz3);
    free(uyy3);
    free(uyz3);
    free(uzz3);
  }
  if (trips==0) 
  {
    PetscOptionsSetValue("-ksp_monitor","");
  }

  if (verbose) 
  {
    CCTK_INFO("CREATING KSP");
  }

  ierr = KSPCreate(pughGH->PUGH_COMM_WORLD,&ksp);  CHKERRQ(ierr);
  
  
  /* 
     Set operators. Here the matrix that defines the linear system
     also serves as the preconditioning matrix. At a later date
     we can probably optimize this using SAME_NONZERO_PATTERN
     and a static trip-though flag (eg, on the second trip
     through do something different). Also we need to think about
     using A to precondition itself. Since I don't know what this
     means, we'll leave it for now.  */

  if (verbose) CCTK_INFO("CREATING KSP OPERATOR");
  ierr = KSPSetOperators(ksp,A[0],A[0]); CHKERRQ(ierr);

  /* Set linear solver defaults for this problem. Later this
     should be replaced/modified with appropriate parsing from the 
     parameter parser. For now it is not. These defaults are
     reasonable, I hope.  */
  if (verbose) CCTK_INFO("KSPGet KSP/PC");
  ierr = KSPGetPC(ksp,&pc);   CHKERRQ(ierr);
  
 /* Get the PC Type */
  if (CCTK_Equals(petsc_PC_type,"PCJACOBI")) 
    ierr = PCSetType(pc,PCJACOBI);
  else if (CCTK_Equals(petsc_PC_type,"PCBJACOBI"))
    ierr = PCSetType(pc,PCBJACOBI);
  else if (CCTK_Equals(petsc_PC_type,"PCICC"))
    ierr = PCSetType(pc,PCICC);
  else if (CCTK_Equals(petsc_PC_type,"PCILU"))
    ierr = PCSetType(pc,PCILU);
  else if (CCTK_Equals(petsc_PC_type,"PCASM"))
    ierr = PCSetType(pc,PCASM);
  else if (CCTK_Equals(petsc_PC_type,"PCLU"))
    ierr = PCSetType(pc,PCLU);
  else if (CCTK_Equals(petsc_PC_type,"PCNONE"))
    ierr = PCSetType(pc,PCNONE);
  else {
    CCTK_WARN(1,"Don't understand petsc_PC_type. Using PCNONE\n");
    ierr = PCSetType(pc,PCNONE);
  }
  CHKERRQ(ierr);

  
  /* Now the same thing for the KSP Type */
  if (CCTK_Equals(petsc_KSP_type,"KSPBCGS")) 
    ierr = KSPSetType(ksp,KSPBCGS); 
  else if (CCTK_Equals(petsc_KSP_type,"KSPRICHARDSON"))
    ierr = KSPSetType(ksp,KSPRICHARDSON);
  else if (CCTK_Equals(petsc_KSP_type,"KSPCHEBYCHEV"))
    ierr = KSPSetType(ksp,KSPCHEBYSHEV);
  else if (CCTK_Equals(petsc_KSP_type,"KSPCHEBYSHEV"))
    ierr = KSPSetType(ksp,KSPCHEBYSHEV);
  else if (CCTK_Equals(petsc_KSP_type,"KSPCG"))
    ierr = KSPSetType(ksp,KSPCG);
  else if (CCTK_Equals(petsc_KSP_type,"KSPGMRES"))
    ierr = KSPSetType(ksp,KSPGMRES);
  else if (CCTK_Equals(petsc_KSP_type,"KSPCGS"))
    ierr = KSPSetType(ksp,KSPCGS);
  else if (CCTK_Equals(petsc_KSP_type,"KSPTFQMR"))
    ierr = KSPSetType(ksp,KSPTFQMR);
  else if (CCTK_Equals(petsc_KSP_type,"KSPTCQMR"))
    ierr = KSPSetType(ksp,KSPTCQMR);
  else if (CCTK_Equals(petsc_KSP_type,"KSPCR"))
    ierr = KSPSetType(ksp,KSPCR);
  else if (CCTK_Equals(petsc_KSP_type,"KSPLSQR"))
    ierr = KSPSetType(ksp,KSPLSQR);
  else {
    CCTK_WARN (1,"I don't understand petsc_KSP_type. Using BiCGSTAB\n");
    ierr = KSPSetType(ksp,KSPBCGS);
  }
  CHKERRQ(ierr);



 /* Set up tolerances */
  
  if (PetscTolStyle == ELLCONV_ABSOLUTE) {
    rtol = 1.0e-15; 
    atolerance = tolerance;
  } else if (PetscTolStyle == ELLCONV_RELATIVE) {
    rtol = tolerance; 
    atolerance = 1.0e-15;
  } else if (PetscTolStyle == ELLCONV_EITHER) {
    rtol = tolerance; 
    atolerance = tolerance;
  } else {
    printf("PETSC Solver: PetscTolStyle set incorrectly [%d]\n",
            PetscTolStyle);
  }
  
  ierr = KSPSetTolerances(ksp,rtol,atolerance,PETSC_DEFAULT,
                          PETSC_DEFAULT);
  CHKERRQ(ierr);


    /* We are warned in the manual that 

     The default technique for orthogonalization of the Hessenberg matrix
     in GMRES is the modified Gram-Schmidt method, which employs many
     VecDot() operations and can thus be slow in parallel. A fast approach
     is to use the unmodified Gram-Schmidt method, which can be set with

     ierr = KSPGMRESSetOrthogonalization(KSP ksp, 
        KSPGMRESUnmodifiedGramSchmidtOrthogonalization); 

     or the options database command
     -ksp_gmres_unmodifiedgramschmidt. Note that this algorithm is
     numerically unstable, but may deliver much better speed
     performance. One can also use unmodifed Gram-Schmidt with
     iterative refinement, by setting the orthogonalization routine,
     KSPGMRESIROrthog(), by using the command line option
     -ksp_gmres_irorthog.

     So this my not be a smart thing to do. But lets put it here 
     so we can turn it on or off later.

   */
  /*$ierr = KSPGMRESSetOrthogonalization(ksp, 
                KSPGMRESUnmodifiedGramSchmidtOrthogonalization); 
  CHKERRQ(ierr);$*/

  /* We also learn that  ...

     By default, KSP assumes an initial guess of zero by zeroing the
     initial value for the solution vector that is given. To use a
     nonzero initial guess, the user must call

     ierr = KSPSetInitialGuessNonzero(KSP ksp); 

   */
  if (verbose) CCTK_INFO("KSPSetInitialGuess\n");
#if PETSC_VERSION_MAJOR < 2 || (PETSC_VERSION_MAJOR == 2 && PETSC_VERSION_MINOR <= 0 && PETSC_VERSION_SUBMINOR < 5)
  ierr = KSPSetInitialGuessNonzero(ksp); CHKERRQ(ierr);
#else
  ierr = KSPSetInitialGuessNonzero(ksp, PETSC_TRUE); CHKERRQ(ierr);
#endif
  
  /* 
    Set runtime options. Since we don't use PETSC runtime options
    but rather use the CACTUS parameter parser, we do this before
    we parse things out. But that may not be such a good idea.
    So lets put it here.
  */
  ierr = KSPSetFromOptions(ksp);   CHKERRQ(ierr);
    

  /* OK so finally we are able to ... */

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                      Solve the linear system
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  if (verbose) CCTK_INFO("KSP solve...");
  ierr = KSPSolve(ksp,b,soln); CHKERRQ(ierr); 
  ierr = KSPGetIterationNumber(ksp,&its); CHKERRQ(ierr); 
  
  /* Here we can form a "res = Ax  - b" and find its norm to get
     an idea of how well we've converged. Put this in later
     but we can do it using MatMultAdd() after we flip the sign
     on b with VecScale() (into res, a new vector, then find the
     norm of res with VecNorm()
  */
  
  /* Now since we now have local layout, we can just get the values
     from the soln vector now matter how many processors we are
     on (eg, the VecCreateMPI has made the solution local) */

  VecGetArray(soln,&values); 
  for (k=kmin;k<kmax;k++) {
    for (j=jmin;j<jmax;j++) {
      for (i=imin;i<imax;i++) {
        /* Set up indices */
        int ijk;                /* The data point for the array */
        PetscInt I;             /* The col in the matrix */

        /* these guys are easy */
        ijk = DATINDEX(pEx,i,j,k);
        if (wsp[ijk] >= 0) {
         
          /* Now this one. "Fortran-order" the matrix. But remember
             we have stripped off the ghost zones. Hence ig-1 and nxs...
          */
          I =wsp[ijk]-startpoint;
          ell_field[ijk] = values[I];
        }
      }
    }
  }
  
  VecRestoreArray(soln,&values);   

  /* Sync var in CCTK speak. */ 
  CCTK_SyncGroup(GH,"ellpetsc::petscworkspace");
  
  /* And finally free up the matrix memory   FIXME 
  ierr = MatReleaseValuesMemory(A); CHKERRQ(ierr); */

  /* This code is not used anymore */
  if (verbose) CCTK_INFO("Destroying Matrices");
  ierr = KSPDestroy(&ksp);  CHKERRQ(ierr);
  ierr = VecDestroy(&soln); CHKERRQ(ierr);
  ierr = VecDestroy(&b);    CHKERRQ(ierr);  
  ierr = MatDestroy(A);     CHKERRQ(ierr);
  free(A);
    
  /*$PetscFinalize();$*/

  if (verbose) CCTK_INFO("LEAVING ELLPETSC");
  
  return (0);
}

