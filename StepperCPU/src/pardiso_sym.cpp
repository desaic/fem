/* -------------------------------------------------------------------- */
/*      Example program to show the use of the "PARDISO" routine        */
/*      on symmetric linear systems                                     */
/* -------------------------------------------------------------------- */
/* -------------------------------------------------------------------- */
/*      This program can be downloaded from the following site:         */
/*      http://www.pardiso-project.org                                  */
/*                                                                      */
/*  (C) Olaf Schenk, Institute of Computational Science                 */
/*      Universita della Svizzera italiana, Lugano, Switzerland.        */
/*      Email: olaf.schenk@usi.ch                                       */
/* -------------------------------------------------------------------- */


// C++ compatible

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "SparseLin.hpp"
#include <omp.h>

#ifdef PARDISO_AVAILABLE

/* PARDISO prototype. */
extern "C" void pardisoinit (void   *, int    *,   int *, int *, double *, int *);
extern "C" void pardiso     (void   *, int    *,   int *, int *,    int *, int *, 
                  double *, int    *,    int *, int *,   int *, int *,
                     int *, double *, double *, int *, double *);
extern "C" void pardiso_chkmatrix  (int *, int *, double *, int *, int *, int *);
extern "C" void pardiso_chkvec     (int *, int *, double *, int *);
extern "C" void pardiso_printstats (int *, int *, double *, int *, int *, int *,
                           double *, int *);

/* Internal solver memory pointer pt,                  */
/* 32-bit: int pt[64]; 64-bit: long int pt[64]         */
/* or void *pt[64] should be OK on both architectures  */
void    *pt[64];

/* Pardiso control parameters. */
int      iparm[64];
double   dparm[64];
int      maxfct, mnum, phase, error, msglvl, solver;

/* Number of processors. */
int      num_procs;

/* Auxiliary variables. */
char    *var;
int      mtype = -2;        /* Real symmetric matrix */

void sparseInit()
{
  error = 0;
  solver = 0; /* use sparse direct solver */
  iparm[34] = 1;
  pardisoinit (pt,  &mtype, &solver, iparm, dparm, &error);

  if (error != 0)
  {
      if (error == -10 )
         printf("No license file found \n");
      if (error == -11 )
         printf("License is expired \n");
      if (error == -12 )
         printf("Wrong username or hostname \n");
       return ;
  }
//  else{
//      printf("[PARDISO]: License check was successful ... \n");
//  }
  /* Numbers of processors, value of OMP_NUM_THREADS */
  var = getenv("OMP_NUM_THREADS");
  if(var != NULL){
      sscanf( var, "%d", &num_procs );
  }
  iparm[2]  = num_procs;

  maxfct = 1;		/* Maximum number of numerical factorizations.  */
  mnum   = 1;         /* Which factorization to use. */

  msglvl = 0;//1;         /* Print statistical information  */
  error  = 0;         /* Initialize error flag */
}

int
sparseSolve(int * ia, int * ja, double * a, int n, double* x, double * b)
{
    int      nnz = ia[n];
    /* RHS and solution vectors. */
    int      nrhs = 1;          /* Number of right hand sides. */
    int      i;
    double   ddum;              /* Double dummy */
    int      idum;              /* Integer dummy. */
/* -------------------------------------------------------------------- */
/* ..  Convert matrix from 0-based C-notation to Fortran 1-based        */
/*     notation.                                                        */
/* -------------------------------------------------------------------- */
    for (i = 0; i < n + 1; i++) {
      ia[i] += 1;
    }
    for (i = 0; i < nnz; i++) {
      ja[i] += 1;
    }
/* -------------------------------------------------------------------- */
/*  .. pardiso_chk_matrix(...)                                          */
/*     Checks the consistency of the given matrix.                      */
/*     Use this functionality only for debugging purposes               */
/* -------------------------------------------------------------------- */
    
//    pardiso_chkmatrix  (&mtype, &n, a, ia, ja, &error);
//    if (error != 0) {
//        printf("\nERROR in consistency of matrix: %d", error);
//    }

/* -------------------------------------------------------------------- */
/* ..  pardiso_chkvec(...)                                              */
/*     Checks the given vectors for infinite and NaN values             */
/*     Input parameters (see PARDISO user manual for a description):    */
/*     Use this functionality only for debugging purposes               */
/* -------------------------------------------------------------------- */

//    pardiso_chkvec (&n, &nrhs, b, &error);
//    if (error != 0) {
//        printf("\nERROR  in right hand side: %d", error);
//    }

/* -------------------------------------------------------------------- */
/* .. pardiso_printstats(...)                                           */
/*    prints information on the matrix to STDOUT.                       */
/*    Use this functionality only for debugging purposes                */
/* -------------------------------------------------------------------- */

//    pardiso_printstats (&mtype, &n, a, ia, ja, &nrhs, b, &error);
//    if (error != 0) {
//        printf("\nERROR right hand side: %d", error);
//    }

/* -------------------------------------------------------------------- */
/* ..  Reordering and Symbolic Factorization.  This step also allocates */
/*     all memory that is necessary for the factorization.              */
/* -------------------------------------------------------------------- */
    phase = 11; 

    pardiso (pt, &maxfct, &mnum, &mtype, &phase,
	     &n, a, ia, ja, &idum, &nrhs,
             iparm, &msglvl, &ddum, &ddum, &error, dparm);
    
    if (error != 0) {
        printf("\nERROR during symbolic factorization: %d", error);
    }
//    printf("\nReordering completed ... ");
//    printf("\nNumber of nonzeros in factors  = %d", iparm[17]);
//    printf("\nNumber of factorization MFLOPS = %d", iparm[18]);
   
/* -------------------------------------------------------------------- */
/* ..  Numerical factorization.                                         */
/* -------------------------------------------------------------------- */    
    phase = 22;
    iparm[32] = 0; /* compute determinant */

    pardiso (pt, &maxfct, &mnum, &mtype, &phase,
             &n, a, ia, ja, &idum, &nrhs,
             iparm, &msglvl, &ddum, &ddum, &error,  dparm);
   
    if (error != 0) {
        printf("\nERROR during numerical factorization: %d", error);
    }
//    printf("\nFactorization completed ...\n ");

/* -------------------------------------------------------------------- */    
/* ..  Back substitution and iterative refinement.                      */
/* -------------------------------------------------------------------- */    
    phase = 33;

    iparm[7] = 1;       /* Max numbers of iterative refinement steps. */
   
    pardiso (pt, &maxfct, &mnum, &mtype, &phase,
             &n, a, ia, ja, &idum, &nrhs,
             iparm, &msglvl, b, x, &error,  dparm);
   
    if (error != 0) {
        printf("\nERROR during solution: %d", error);
    }

    /*double  norm_b = 0.0;
    for (int i = 0; i < n; ++i)
      norm_b  += b[i] * b[i];
    norm_b = sqrt(norm_b );
    printf("norm(b) = %f \n", norm_b );
    double  norm_x = 0.0;
    for (i = 0; i < n; ++i)
      norm_x  += x[i] * x[i];
    norm_x = sqrt(norm_x );
    printf("norm(x) = %f \n", norm_x );*/ 

//    printf("\nSolve completed ... ");
//    printf("\nThe solution of the system is: ");
//    for (i = 0; i < n; i++) {
//        printf("\n x [%d] = % f", i, x[i] );
//    }
//    printf ("\n");

/* -------------------------------------------------------------------- */    
/* ..  Termination and release of memory.                               */
/* -------------------------------------------------------------------- */    
    phase = -1;                 /* Release internal memory. */
    
    pardiso (pt, &maxfct, &mnum, &mtype, &phase,
             &n, &ddum, ia, ja, &idum, &nrhs,
             iparm, &msglvl, &ddum, &ddum, &error,  dparm);

    return 0;
}

#else
void sparseInit()
{
  printf("\nERROR pardiso is not available!\n");
}

int
sparseSolve(int * ia, int * ja, double * a, int n, double* x, double * b)
{
  printf("\nERROR pardiso is not available!\n");
  return -1;
}

#endif
