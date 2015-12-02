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
#include "pardiso_sym.hpp"
#include "Timer.hpp"

#ifdef PARDISO_AVAILABLE

#ifdef __linux__
#include <omp.h>

/* PARDISO prototype. */
extern "C" void pardisoinit (void   *, int    *,   int *, int *, double *, int *);
extern "C" void pardiso     (void   *, int    *,   int *, int *,    int *, int *, 
                  double *, int    *,    int *, int *,   int *, int *,
                     int *, double *, double *, int *, double *);
extern "C" void pardiso_chkmatrix  (int *, int *, double *, int *, int *, int *);
extern "C" void pardiso_chkvec     (int *, int *, double *, int *);
extern "C" void pardiso_printstats (int *, int *, double *, int *, int *, int *,
                           double *, int *);

void pardisoInit(PardisoState * state)
{
  int error = 0, solver = 0; /* use sparse direct solver */
  pardisoinit (state->pt,  &state->mtype, &solver, state->iparm, state->dparm, &error);
  int num_procs=1;
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
  char * var = getenv("OMP_NUM_THREADS");
  if(var != NULL){
      sscanf( var, "%d", &num_procs );
  }
  state->iparm[2]  = num_procs;
  printf("%d omp threads\n", state->iparm[2]);
}

int pardisoSymbolicFactorize(int * ia, int * ja, int n, PardisoState *state)
{
  Timer timer;
  timer.startWall();
  int error = 0;

  /* RHS and solution vectors. */
  int      nrhs = 1;          /* Number of right hand sides. */

  double   ddum;              /* Double dummy */
  int      idum;              /* Integer dummy. */
  int phase=0;

    phase = 11;

    pardiso (state->pt, &state->maxfct, &state->mnum, &state->mtype, &phase,
       &n, &ddum, ia, ja, &idum, &nrhs,
       state->iparm, &state->msglvl, &ddum, &ddum, &error, state->dparm);

    if (error != 0) {
        printf("\nERROR during symbolic factorization: %d", error);
    }
//    printf("\nReordering completed ... ");
//    printf("\nNumber of nonzeros in factors  = %d", iparm[17]);
//    printf("\nNumber of factorization MFLOPS = %d", iparm[18]);

    timer.endWall();
    printf("up to symbolic factorization: %f \n",timer.getSecondsWall());

}

int pardisoNumericalFactorize(int *ia, int *ja, double *val, int n, PardisoState *state)
{
  int error = 0;
  /* RHS and solution vectors. */
  int      nrhs = 1;          /* Number of right hand sides. */

  double   ddum;              /* Double dummy */
  int      idum;              /* Integer dummy. */
  int phase=22;

  pardiso (state->pt, &state->maxfct, &state->mnum, &state->mtype, &phase,
           &n, val, ia, ja, &idum, &nrhs,
           state->iparm, &state->msglvl, &ddum, &ddum, &error,  state->dparm);

  if (error != 0) {
    printf("\nERROR during numerical factorization: %d", error);
  }

  return error;
}

int
pardisoBackSubstitute(int * ia, int * ja, double* val, int n, double* x, double * b, PardisoState *state)
{
  int error = 0;
  /* RHS and solution vectors. */
  int      nrhs = 1;          /* Number of right hand sides. */
  int      idum;              /* Integer dummy. */
  int phase = 33;

  state->iparm[7] = 1;       /* Max numbers of iterative refinement steps. */

  pardiso (state->pt, &state->maxfct, &state->mnum, &state->mtype, &phase,
           &n, val, ia, ja, &idum, &nrhs,
           state->iparm, &state->msglvl, b, x, &error,  state->dparm);
  if (error != 0) {
    printf("\nERROR during solution: %d", error);
  }
  return 0;
}

int
pardisoSolve(int * ia, int * ja, double * a, int n, double* x, double * b,
             PardisoState *state)
{
//  Timer timer;
//  timer.startWall();
  pardisoNumericalFactorize(ia, ja, a, n, state);
//  timer.endWall();
//  printf("numerical factorization: %f\n",timer.getSecondsWall());

//  timer.startWall();
  pardisoBackSubstitute(ia,ja,a,n,x,b,state);
//  printf("back substitution: %f \n", timer.getSecondsWall());
//  timer.endWall();
}

void pardisoFree(int * ia, int * ja, int n, PardisoState *state)
{
  int phase = 0;
  int error = 0;

  /* RHS and solution vectors. */
  int      nrhs = 1;          /* Number of right hand sides. */

  double   ddum;              /* Double dummy */
  int      idum;              /* Integer dummy. */

  /* -------------------------------------------------------------------- */
  /* ..  Termination and release of memory.                               */
  /* -------------------------------------------------------------------- */


  phase = -1;                 /* Release internal memory. */

  pardiso (state->pt, &state->maxfct, &state->mnum, &state->mtype, &phase,
           &n, &ddum, ia, ja, &idum, &nrhs,
           state->iparm, &state->msglvl, &ddum, &ddum, &error,  state->dparm);

}

#endif

//=======================linux end==============

#ifdef __APPLE__
#include <mkl.h>

/* Internal solver memory pointer pt, */
/* 32-bit: int pt[64]; 64-bit: long int pt[64] */
/* or void *pt[64] should be OK on both architectures */
void *pt[64];
/* Pardiso control parameters. */
MKL_INT iparm[64];


void pardisoInit(){
MKL_INT i;
for ( i = 0; i < 64; i++ )
    {
        iparm[i] = 0;
    }
    iparm[0] = 1;         /* No solver default */
    iparm[1] = 2;         /* Fill-in reordering from METIS */
    iparm[3] = 0;         /* No iterative-direct algorithm */
    iparm[4] = 0;         /* No user fill-in reducing permutation */
    iparm[5] = 0;         /* Write solution into x */
    iparm[6] = 0;         /* Not in use */
    iparm[7] = 2;         /* Max numbers of iterative refinement steps */
    iparm[8] = 0;         /* Not in use */
    iparm[9] = 13;        /* Perturb the pivot elements with 1E-13 */
    iparm[10] = 1;        /* Use nonsymmetric permutation and scaling MPS */
    iparm[11] = 0;        /* Not in use */
    iparm[12] = 0;        /* Maximum weighted matching algorithm is switched-off (default for symmetric). Try iparm[12] = 1 in case of inappropriate accuracy */
    iparm[13] = 0;        /* Output: Number of perturbed pivots */
    iparm[14] = 0;        /* Not in use */
    iparm[15] = 0;        /* Not in use */
    iparm[16] = 0;        /* Not in use */
    iparm[17] = -1;       /* Output: Number of nonzeros in the factor LU */
    iparm[18] = -1;       /* Output: Mflops for LU factorization */
    iparm[19] = 0;        /* Output: Numbers of CG Iterations */
    
/* -------------------------------------------------------------------- */
/* .. Initialize the internal solver memory pointer. This is only */
/* necessary for the FIRST call of the PARDISO solver. */
/* -------------------------------------------------------------------- */
    for ( i = 0; i < 64; i++ )
    {
        pt[i] = 0;
    }
}

int
pardisoSolve(int * ia, int * ja, double * a, int n, double* x, double * b)
{
  MKL_INT mtype = -2;       /* Real symmetric matrix */
  /* RHS and solution vectors. */
  MKL_INT nrhs = 1;     /* Number of right hand sides. */

  MKL_INT maxfct, mnum, phase, error, msglvl;

  MKL_INT i;
  double ddum;          /* Double dummy */
  MKL_INT idum;         /* Integer dummy. */
  MKL_INT      nnz = ia[n];
/* -------------------------------------------------------------------- */
/* ..  Convert matrix from 0-based C-notation to Fortran 1-based        */
/*     notation.                                                        */
/* -------------------------------------------------------------------- */
    for (i = 0; i < n+1; i++) {
        ia[i] += 1;
    }
    for (i = 0; i < nnz; i++) {
        ja[i] += 1;
    }

  maxfct = 1;           /* Maximum number of numerical factorizations. */
  mnum = 1;         /* Which factorization to use. */
  msglvl = 0;           /* Print statistical information in file */
  error = 0;            /* Initialize error flag */
  phase = 11;
  PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
           &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
//   if ( error != 0 )
//   {
//       printf ("\nERROR during symbolic factorization: %d", error);
//       exit (1);
//   }
//   printf ("\nReordering completed ... ");
//   printf ("\nNumber of nonzeros in factors = %d", iparm[17]);
//   printf ("\nNumber of factorization MFLOPS = %d", iparm[18]);
/* -------------------------------------------------------------------- */
/* .. Numerical factorization. */
/* -------------------------------------------------------------------- */
  phase = 22;
  PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
           &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
//   if ( error != 0 )
//   {
//       printf ("\nERROR during numerical factorization: %d", error);
//       exit (2);
//   }
//   printf ("\nFactorization completed ... ");
/* -------------------------------------------------------------------- */
/* .. Back substitution and iterative refinement. */
/* -------------------------------------------------------------------- */
  phase = 33;
  iparm[7] = 2;         /* Max numbers of iterative refinement steps. */
  
  PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
           &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, b, x, &error);
  // if ( error != 0 )
//   {
//       printf ("\nERROR during solution: %d", error);
//       exit (3);
//   }
//   printf ("\nSolve completed ... ");
//   printf ("\nThe solution of the system is: ");
//   for ( i = 0; i < n; i++ )
//   {
//       printf ("\n x [%d] = % f", i, x[i]);
//   }
//   printf ("\n");
/* -------------------------------------------------------------------- */
/* .. Termination and release of memory. */
/* -------------------------------------------------------------------- */
  phase = -1;           /* Release internal memory. */
  PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
           &n, &ddum, ia, ja, &idum, &nrhs,
           iparm, &msglvl, &ddum, &ddum, &error);
  return 0;

}


#endif



#else

//========================no pardiso============
void pardisoInit()
{
  printf("\nERROR pardiso is not available!\n");
}

int
pardisoSolve(int * ia, int * ja, double * a, int n, double* x, double * b)
{
  printf("\nERROR pardiso is not available!\n");
  return -1;
}

#endif
