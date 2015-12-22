#include "ArrayUtil.hpp"
void BBox(const std::vector<Eigen::Vector3f >& v,
  Eigen::Vector3f & mn, Eigen::Vector3f & mx)
{
  mn = v[0];
  mx = v[0];
  for(unsigned int ii = 1 ;ii<v.size();ii++){
    for(int dim = 0 ; dim<3;dim++){
      if(v[ii][dim]<mn[dim]){
        mn[dim] = v[ii][dim];
      }
      if(v[ii][dim]>mx[dim]){
        mx[dim] = v[ii][dim];
      }
    }
  }
}
