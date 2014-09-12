#include "MatrixXd.hpp"

MatrixXd::MatrixXd():M(0),mm(0),nn(0){}

MatrixXd::MatrixXd(int _m, int _n):M(0), mm(_m), nn(_n){allocate(_m,_n);}

MatrixXd::MatrixXd(const MatrixXd & ma):M(0),mm(ma.mm),nn(ma.nn){
  allocate(mm,nn);
  copyArray(ma.M, M, mm*nn);
}

void MatrixXd::fill(double val){
  int nEntry=mm*nn;
  for(int ii = 0; ii<nEntry; ii++){
    M[ii] = val;
  }
}

MatrixXd&MatrixXd::operator=(const MatrixXd & ma){
  if( mm * nn != ma.mm * ma.nn && M!=0 ){
    delete[]M;
    M=0;
  }
  mm=ma.mm;
  nn=ma.nn;
  if(M==0){
    M=new double[mm*nn];
  }
  copyArray(ma.M,M,mm*nn);
  return *this;
}

void MatrixXd::allocate(int _m, int _n){
  if( mm * nn != _m * _n && M!=0){
    delete[]M;
    M=0;
  }
  nn = _n;
  mm = _m;
  if(M==0){
    M=new double[mm*nn];
  }
}

void MatrixXd::resize(int _m, int _n)
{
  if(M==0){
    allocate(_m,_n);
    return;
  }
  if(mm==_m && nn==_n){
    return;
  }

  double * Mold = M;
  M = new double[ _m*_n ];
  for(int ii = 0; ii<_m; ii++){
    for(int jj= 0; jj<_n; jj++){
      if(ii<mm && jj<nn){
        M[ii*_n + jj] = Mold[ii*nn + jj];
      }else{
        M[ii*_n + jj] = 0;
      }
    }
  }
  mm=_m;
  nn=_n;
  delete [] Mold;
}

int MatrixXd::get1dIndex(int ii, int jj)const{
  return ii*nn + jj;
}

double & MatrixXd::operator()(int ii , int jj){
  return  M[get1dIndex(ii,jj)];
}

double MatrixXd::operator()(int ii , int jj)const{
  return  M[get1dIndex(ii,jj)];
}

MatrixXd::~MatrixXd(){
  if(M!=0){
    delete []M;
  }
}

void MatrixXd::print(std::ostream & stream){
  for(int ii = 0;ii<mm;ii++){
    for(int jj = 0;jj<nn;jj++){
      stream<<operator ()(ii,jj)<<" ";
    }
    stream<<"\n";
  }
  stream<<"\n";
}

