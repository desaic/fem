#ifndef UTIL_HPP
#define UITL_HPP

#include <vector>

template<class T>
inline T clamp(T a, T lower, T upper)
{
   if(a<lower) return lower;
   else if(a>upper) return upper;
   else return a;
}


template<typename T>
std::vector<T> operator+(const std::vector<T> & a, const std::vector<T> & b)
{
  std::vector<T> sum(a.size());
#ifdef DEBUG
  if(a.size()!=b.size()){
    std::cout<<"adding vectors of different sizes\n";
  }
#endif
  for(unsigned int ii =0;ii<a.size();ii++){
    sum[ii] = a[ii] + b[ii];
  }
  return sum;
}

template<typename T>
std::vector<T> operator*(float a, const std::vector<T> & b)
{
  std::vector<T> prod(b.size());
  for(unsigned int ii =0;ii<b.size();ii++){
    prod[ii] = a * b[ii];
  }
  return prod;
}

#endif