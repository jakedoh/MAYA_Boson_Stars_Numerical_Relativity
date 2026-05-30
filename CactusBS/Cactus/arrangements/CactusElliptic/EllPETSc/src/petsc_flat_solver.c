 /*@@
   @file      petsc_flat_solver.c
   @date      
   @author    Gerd Lanfermann
   @desc 
   
   @enddesc 
 @@*/

#include <stdlib.h>

#include <cctk.h>
#include <cctk_CommandLine.h>
#include <cctk_Parameters.h>

#include <petsc.h>

#include "ellpetsc.h"
#include "CactusPUGH/PUGH/src/include/pugh.h"

/* #define DEBUG */

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

/* A convenient macros for the stencil molecule*/
#define a(i,j,k) a[i+1 + 3*(j+1 + 3*(k+1))]
  

static int trips=0;
static Mat     *A;           /* linear system matrix */
static Vec     soln, b;      /* approx solution, RHS */
static KSP     ksp;          /* linear solver context */

int petsc_flat(cGH *GH, int FieldIndex, int MIndex, int NIndex, 
                int *AbsTol, int *RelTol);

int petsc_flat(cGH *GH, int FieldIndex, int MIndex, int NIndex, 
                int *AbsTol, int *RelTol) 
{

  DECLARE_CCTK_PARAMETERS    /* CCTK passed parameters */

  PC     pc;            /* preconditioner context */
  int    num_A;         /* Number of A-arrays as needed ny Dan's MG */
  

  int    ierr;          /* Check the return status of petsc */
  int    retcode;       /* Check the return status of CCTK */

  int    rank;         /* Rank of the matrix/vector */
  PetscInt its;        /* Number of iterations */
  
  CCTK_REAL a[27];      /* The stencil array */
  CCTK_REAL ac;         /* Storage for a(0,0,0) for renorm */
  
  /* Loopers */
  int i,j,k,l,m,n;

  /* loop limits put in for stencil_w !=1     Ed Evans 1/19/98 */
  int imin,imax,jmin,jmax,kmin,kmax;
   
  pGH *pughGH;          /* The pugh Extension handle */
  pGExtras *pEx;
  pConnectivity *pCon;

  /* Tolerances */
  CCTK_REAL rtol=0, atolerance=0, tolerance;

  /* Values to assemble the matrix */
  CCTK_REAL *values;

  int myproc; 
  int nxs,nys,nzs;      /* Size of the grid with stencils off... */
  int startpoint=0;     /* My starting index (per proc) */
  int endpoint;         /* One more than my end */
  PetscInt pstart, pend;     /* A check for PETSc layout */
  PetscInt pvstart, pvend;   /* A check for PETSc layout */
  int verbose;          /* Is the solver verbose */
  int debug;            /* Is the solver debug-verbose */
  int octant;           /* Apply octant BCs inside */
  int matnormalize;     /* Normalize the central mat value to one? */
  PetscInt I,J;              /* The PETSc col and row in the matrix */
            
  int PetscTolStyle;

  
  /* Pointers to the data of : petsc workspace/ the GF to solve /
     Mlinear / Nlinear(source) */
  CCTK_REAL *wsp =NULL, *ell_field=NULL; 
  CCTK_REAL *Mlin=NULL, *Nlin=NULL;   

  /* flags to signal if storage is ON or OFF (used/not-used) for M/N*/
  int Mstorage=0, Nstorage=0;

  RelTol = RelTol;

  octant       = CCTK_Equals(domain,"octant");
  verbose      = CCTK_Equals(petsc_verbose,"yes")||
                 CCTK_Equals(petsc_verbose,"debug");
  debug        = CCTK_Equals(petsc_verbose,"debug");
  matnormalize = 0; 
  

  /* FIXME, the TolAbs/TolRel will be evaluated here */
  PetscTolStyle=0;
  tolerance    =AbsTol[0];

  /* Get the link to pugh Extension */
  pughGH = (pGH*)GH->extensions[CCTK_GHExtensionHandle("PUGH")];
  if (!pughGH) 
  {
    CCTK_WARN(0,"ETERNAL ERROR: Cannot find PUGH Extension Handle");
  }

  /* Get the extras extension for 3D grid functions */
  pEx = pughGH->GFExtras[2];
  pCon = pughGH->Connectivity[2];

  /* Things to do on first iteration */
  if (trips==0) 
  {
    int argc;
    char **argv; 

#ifdef DEBUG
    if (debug)
      printf("PETSc: initial trip: %d \n",trips);
#endif

    argc = CCTK_CommandLine(&argv);
    PETSC_COMM_WORLD = pughGH->PUGH_COMM_WORLD;
    PetscInitialize(&argc,&argv,NULL,NULL); 

    CCTK_INFO("PETSc initialized");
  }

  trips++;

  /* Allocate storage for system matrix, will have more of them later */
  num_A = 1;
  A=(Mat*) malloc (sizeof(Mat)*num_A);


  /* Get the data pointers */
  
  /* if we have a negative index, this GF is not needed, 
     there for don't even look for it. when index positive,
     set flag Mstorage=1, dito for N */
  if (MIndex>=0)  
  { 
    Mlin = (CCTK_REAL*) CCTK_VarDataPtrI(GH,0,MIndex);
    Mstorage = 1;
  }
  if (NIndex>=0) 
  {
    Nlin = (CCTK_REAL*) CCTK_VarDataPtrI(GH,0,NIndex);
    Nstorage = 1;
  }
  ell_field = (CCTK_REAL*) CCTK_VarDataPtrI(GH,0,FieldIndex);
  wsp       = (CCTK_REAL*) CCTK_VarDataPtrI(GH,0,CCTK_VarIndex("ellpetsc::wsp"));
  

  /* initialize the linear index lookup table with -1; 
     this will indicate a boudary later */
  for (i=0;i<pEx->npoints;i++) 
  {
    wsp[i] = -1.0;
  }

  /* Get myproc from the CCTK, get the gridspacing */
  myproc = PUGH_MyProc(GH);

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

  for (i=0;i<myproc;i++) 
  {
    
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

#ifdef DEBUG  
  printf("ENDPOINT:  %d \n",endpoint);
#endif 

  /* So now each point has a unique index of its own. If we
     do a sync that means each processor knows the indiced
     in the ghost zones (eg, on its neighbors) in a FD stencil of 1 */
  retcode = CCTK_SyncGroup(GH,"ellpetsc::petscworkspace");
  if (retcode<0) 
  {
    CCTK_WARN(1,"Synchronization failed");
  }
  
  /* So woohoo. Now for each point in our ijk space, we have
     information about our row in the matrix (workspace->data[ijk])
     as well as the column for the stencil of the neighbors
     (workspace->data[DATINDEX(GH,i+1,j,k)] etc...). So onwards */

  nxs = pEx->nsize[0]-2;
  nys = pEx->nsize[1]-2;
  nzs = pEx->nsize[2]-2;

#ifdef DEBUG
  if (debug) 
  {
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
  if (debug) 
  { 
    printf("CAC M-Layout: %d %d \n", (int)startpoint, (int)endpoint);
    printf("PET M-Layout: %d %d \n", (int)pstart, (int)pend);
    printf("PET V-Layout: %d %d \n", (int)pvstart,(int)pvend);
  }
  if (pstart != startpoint && pend != endpoint) 
    CCTK_WARN(1,"WARNING: PETSC and data layouts differ! (why?)");
         

  /* Initialize the stencil array */
  for (I=0;I<27;I++) a[I] = 0.0;
  a( 0, 0, 0) =-6.0;
  a( 1, 0, 0) = 1.0;  /* ae */
  a(-1, 0, 0) = 1.0;  /* aw */
  a( 0, 1, 0) = 1.0;  /* an */       
  a( 0,-1, 0) = 1.0;  /* as */      
  a( 0, 0, 1) = 1.0;  /* at */              
  a( 0, 0,-1) = 1.0;  /* ab */


  for (k=kmin;k<kmax;k++) 
  {
    for (j=jmin;j<jmax;j++) 
    {
      for (i=imin;i<imax;i++) 
      {

        if (wsp[DATINDEX(pEx,i,j,k)] >= 0) 
        {

          /* Set up indices */
          int ijk;              /* The lineax index for the array */
          int ig, jg, kg;       /* The position in global space */
          CCTK_REAL rhsval = 0;

          /* local linear index / global 3d index */
          ijk = DATINDEX(pEx,i,j,k);    
          ig  = i + GH->cctk_lbnd[0]; 
          jg  = j + GH->cctk_lbnd[1];
          kg  = k + GH->cctk_lbnd[2];
        
          /* Get the row we are working on */
          I = wsp[ijk];
       
          /* M phi */
          if (Mstorage) 
          {
            a(0,0,0) = a(0,0,0) + Mlin[ijk];
          }

          /* Now set the values of the matrix. Here we force dirichlet. */
          VecSetValues(soln,1,&I,&(ell_field[ijk]),INSERT_VALUES);
          
#ifdef DEBUG
          if (ijk==1918) 
          {
            for (l=-1;l<=1;l++)
            {
              for (m=-1;m<=1;m++)
              {
                for (n=-1;n<=1;n++) 
                {
                  printf(" (%d %d %d): %f \n",l,m,n,a(l,m,n));
                }
              }
            }
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
          if (matnormalize) 
          {
            ac = a(0,0,0);
            for (J=0;J<27;J++) a[J] = a[J] / ac;
          } 
          else 
            ac = 1;
          

          /* This is the new way-clever look. Note it relies
             heavily on the index array in the workspace and
             on multiple processors it will *NOT* make a
             19-way banded matrix. */

          for (l=-1;l<=1;l++)
          {
            for (m=-1;m<=1;m++)
            {
              for (n=-1;n<=1;n++) 
              {
                if (l*m*n == 0) 
                {       /* This is the 19 point ... if none are
                         * zero, then we have no stencil here.
                         */
                  if (wsp[DATINDEX(pEx,i+l,j+m,k+n)] < 0) 
                  {
                    /* This is a boundary. */
                    rhsval += a(l,m,n) * ell_field[DATINDEX(pEx,i+l,j+m,k+n)];
                  } 
                  else 
                  {
                    J = wsp[DATINDEX(pEx,i+l,j+m,k+n)];
                    ierr = MatSetValues(A[0],1,&I,1,&J,&(a(l,m,n)),INSERT_VALUES); 
                    CHKERRQ(ierr);
                  }
                }
              }
            }
          }
          if (Nstorage) 
          {
            rhsval = -rhsval - Nlin[ijk] / ac;
          }

          ierr   = VecSetValues(b,1,&I,&rhsval,INSERT_VALUES);
          CHKERRQ(ierr);
          
          /* reset the central value*/
          a( 0, 0, 0) = -6.0;
  
        }   /* if (wsp>=0) */
      }   /*for i */
    }   /* for j */
  }   /* for k */
  
  if (verbose) 
  {
    CCTK_INFO ("Assembling the vectors");
  }
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
  
  if (trips==0) 
  {
#if PETSC_VERSION_MAJOR < 2 || (PETSC_VERSION_MAJOR == 2 && PETSC_VERSION_MINOR == 0)    
    OptionsSetValue("-ksp_monitor","");
#else  
    PetscOptionsSetValue("-ksp_monitor","");
#endif
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

  if (verbose) 
  {
    CCTK_INFO("CREATING KSP OPERATOR");
  }
  ierr = KSPSetOperators(ksp,A[0],A[0]); CHKERRQ(ierr);

  /* Set linear solver defaults for this problem. Later this
     should be replaced/modified with appropriate parsing from the 
     parameter parser. For now it is not. These defaults are
     reasonable, I hope.  */
  if (verbose) 
  {
    CCTK_INFO("KSPGet KSP/PC");
  }
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
  else 
  {
    CCTK_WARN(1,"Don't understand petsc_PC_type. Using PCNONE");
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
  else 
  {
    CCTK_WARN (1,"I don't understand petsc_KSP_type. Using BiCGSTAB");
    ierr = KSPSetType(ksp,KSPBCGS);
  }
  CHKERRQ(ierr);



 /* Set up tolerances */
  
  if (PetscTolStyle == ELLCONV_ABSOLUTE) 
  {
    rtol = 1.0e-15; 
    atolerance = tolerance;
  } 
  else if (PetscTolStyle == ELLCONV_RELATIVE) 
  {
    rtol = tolerance; 
    atolerance = 1.0e-15;
  } 
  else if (PetscTolStyle == ELLCONV_EITHER) 
  {
    rtol = tolerance; 
    atolerance = tolerance;
  } 
  else 
  {
    printf("PETSC Solver: PetscTolStyle set incorrectly [%d]\n",
            PetscTolStyle);
  }
  
  ierr = KSPSetTolerances(ksp, rtol, atolerance, PETSC_DEFAULT, PETSC_DEFAULT);
  CHKERRQ(ierr);

  if (verbose) 
  {
    CCTK_INFO("KSPSetInitialGuess");
  }
#if PETSC_VERSION_MAJOR < 2 || (PETSC_VERSION_MAJOR == 2 && PETSC_VERSION_MINOR <= 0)
  ierr = KSPSetInitialGuessNonzero(ksp); 
#else
  ierr = KSPSetInitialGuessNonzero(ksp, PETSC_TRUE); 
#endif
  CHKERRQ(ierr);
  
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
  if (verbose) 
  {
    CCTK_INFO("KSP solve...");
  }

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
  for (k=kmin;k<kmax;k++) 
  {
    for (j=jmin;j<jmax;j++) 
    {
      for (i=imin;i<imax;i++) 
      {
        /* Set up indices */
        int ijk;                /* The data point for the array */

        /* these guys are easy */
        ijk = DATINDEX(pEx,i,j,k);
        if (wsp[ijk] >= 0) 
        {
          /* Now this one. "Fortran-order" the matrix. But remember
             we have stripped off the ghost zones. Hence ig-1 and nxs */
          I = wsp[ijk]-startpoint;
          ell_field[ijk] = values[I];
        }
      }
    }
  }
  
  VecRestoreArray(soln,&values);   

  CCTK_SyncGroup(GH,"ellpetsc::petscworkspace");
  
  /* And finally free up the matrix memory */  
  /*$ierr = MatReleaseValuesMemory(A); CHKERRQ(ierr);$*/ 

  /* This code is not used anymore */
  if (verbose) 
  {
    CCTK_INFO("Destroying Matrices");
  }
  ierr = KSPDestroy(&ksp);  CHKERRQ(ierr);
  ierr = VecDestroy(&soln); CHKERRQ(ierr);
  ierr = VecDestroy(&b);    CHKERRQ(ierr);  
  ierr = MatDestroy(A);     CHKERRQ(ierr);
  free(A);
    
  /*$PetscFinalize();$*/

  if (verbose) 
  {
    CCTK_INFO("LEAVING ELLPETSC");
  }
  
  return (0);
}
