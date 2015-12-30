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

#include <omp.h>

/* PARDISO prototype. */
extern "C" void pardisoinit(void   *, int    *, int *, int *, double *, int *);
extern "C" void pardiso     (void   *, int    *,   int *, int *,    int *, int *, 
                  double *, int    *,    int *, int *,   int *, int *,
                     int *, double *, double *, int *, double *);
extern "C" void pardiso_chkmatrix  (int *, int *, double *, int *, int *, int *);
extern "C" void pardiso_chkvec     (int *, int *, double *, int *);
extern "C" void pardiso_printstats (int *, int *, double *, int *, int *, int *,
                           double *, int *);

void pardisoInit(PardisoState * state)
{
  state->mtype = -2;
  int error = 0, solver = 0; /* use sparse direct solver */
  pardisoinit(state->pt, &state->mtype, &solver, state->iparm, state->dparm, &error);
  state->maxfct = 1;		/* Maximum number of numerical factorizations.  */
  state->mnum = 1;         /* Which factorization to use. */
  state->msglvl = 0;//1;         /* Print statistical information  */

  int num_procs = 1;
  if (error != 0)
  {
    if (error == -10)
      printf("No license file found \n");
    if (error == -11)
      printf("License is expired \n");
    if (error == -12)
      printf("Wrong username or hostname \n");
    return;
  }
  //  else{
  //      printf("[PARDISO]: License check was successful ... \n");
  //  }
  /* Numbers of processors, value of OMP_NUM_THREADS */
  char * var = getenv("OMP_NUM_THREADS");
  if (var != NULL){
    sscanf(var, "%d", &num_procs);
  }
  state->iparm[2] = num_procs;
  printf("%d omp threads\n", state->iparm[2]);
}

int pardisoSymbolicFactorize(int * ia, int * ja, int n, PardisoState *state)
{
  //Timer timer;
  //timer.startWall();
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

    //timer.endWall();
    //printf("up to symbolic factorization: %f \n",timer.getSecondsWall());
    return error;
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
  int status = pardisoNumericalFactorize(ia, ja, a, n, state);
//  timer.endWall();
//  printf("numerical factorization: %f\n",timer.getSecondsWall());

//  timer.startWall();
  pardisoBackSubstitute(ia,ja,a,n,x,b,state);
//  printf("back substitution: %f \n", timer.getSecondsWall());
//  timer.endWall();
  return status;
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

//=======================linux end==============

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
