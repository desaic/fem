#ifndef MATRIXXD_HPP
#define MATRIXXD_HPP
#include <iostream>

template <typename T>
void copyArray(const T* src, T*dst, int size)
{
  for(int ii = 0;ii<size;ii++){
    dst[ii] = src[ii];
  }
}

struct MatrixXd
{
  MatrixXd():M(0),mm(0),nn(0){}
  
  MatrixXd(const MatrixXd & ma):mm(ma.mm),nn(ma.nn){
    allocate(mm,nn);
    copyArray(ma.M, M, mm*nn);
  }
  
  MatrixXd&operator=(const MatrixXd & ma){
    if( mm * nn != ma.mm * ma.nn && M!=0 ){
      delete[]M;
      M=0;
    }
    if(M==0){
      M=new double[mm*nn];
    }
    mm=ma.mm;
    nn=ma.nn;
    copyArray(ma.M,M,mm*nn);
  }

  ///@brief array of matrix entries
  double * M;
  
  ///@brief number of rows.
  int mm;
  ///@brief number of columns
  int nn;
  
  ///@param _m number of rows.
  ///@param _n number of columns.
  void allocate(int _m, int _n){
    if( mm * nn != _m * _n && M!=0){
      delete[]M;
      M=0;
    }
    if(M==0){
      M=new double[mm*nn];
    }
    nn = _n;
    mm = _m;
  }
  
  virtual int get1dIndex(int ii, int jj)const{
    return ii*nn + jj;
  }

  ///@param ii row number
  ///@param jj column number
  double & operator()(int ii , int jj){
    return  M[get1dIndex(ii,jj)];
  }
  
  double operator()(int ii , int jj)const{
    return  M[get1dIndex(ii,jj)];
  }

  ~MatrixXd(){
    if(M!=0){delete []M;}
  }

  void print(std::ostream & stream){
    for(int ii = 0;ii<mm;ii++){
      for(int jj = 0;jj<nn;jj++){
        stream<<operator ()(ii,jj)<<" ";
      }
      stream<<"\n";
    }
    stream<<"\n";
  }
};

template<typename T1, typename T2>
void copySubMat(const T1 & src, T2 & dst, int si1, int sj1, int si2, int sj2,
                int m, int n)
{
  for(int ii = 0;ii<m;ii++){
    for(int jj = 0;jj<n;jj++){
      dst(si2 + ii, sj2 + jj) = src(si1+ii, sj1+jj);
    }
  }
}

#endif
