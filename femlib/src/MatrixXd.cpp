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

