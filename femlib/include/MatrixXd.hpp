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
  MatrixXd();
  
  MatrixXd(int _m, int _n);

  MatrixXd(const MatrixXd & ma);
  
  void fill(double val);

  MatrixXd&operator=(const MatrixXd & ma);

  ///@brief array of matrix entries
  double * M;
  
  ///@brief number of rows.
  int mm;
  ///@brief number of columns
  int nn;
  
  ///@param _m number of rows.
  ///@param _n number of columns.
  void allocate(int _m, int _n);

  ///@brief preserves values if possible.
  ///Add 0s if necessary.
  void resize(int _m, int _n);
  
  virtual int get1dIndex(int ii, int jj)const;

  ///@param ii row number
  ///@param jj column number
  double & operator()(int ii , int jj);
  
  double operator()(int ii , int jj)const;

  ~MatrixXd();

  void print(std::ostream & stream);
};

template<typename T1, typename T2>
void addSubMat(const T1 & src, T2 & dst, int si1, int sj1, int si2, int sj2,
                int m, int n)
{
  for(int ii = 0;ii<m;ii++){
    for(int jj = 0;jj<n;jj++){
      dst(si2 + ii, sj2 + jj) += src(si1+ii, sj1+jj);
    }
  }
}

template<typename T1, typename T2>
void copyMat(const T1 & src, T2& dst, int m, int n)
{
  for(int ii = 0;ii<m;ii++){
    for(int jj = 0;jj<n;jj++){
      dst(ii,jj) = src(ii,jj);
    }
  }
}

#endif
