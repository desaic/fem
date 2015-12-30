#ifndef SPARSELIN_HPP
#define SPARSELIN_HPP

void sparseInit();
int sparseSolve(int * ia, int * ja, double* val,int n, double* x, double * b);

#endif
