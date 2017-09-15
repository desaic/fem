/// \file
#ifndef MICROSTRUCTURE_PARAM_2D_HPP
#define MICROSTRUCTURE_PARAM_2D_HPP

#include <iostream>
#include <vector>
#include <math.h>
#include <Eigen/Dense>

struct Rectangle2D
{
  Rectangle2D() :rot(0){
    center[0] = 0;
    center[1] = 0;
    length[0] = 0;
    length[1] = 0;
  }
  float center[2];
  float length[2];
  float rot;
};

struct Ellipse2D
{
  Ellipse2D() :rot(0){
    center[0] = 0;
    center[1] = 0;
    r[0] = 0;
    r[1] = 0;
  }
  float center[2];
  float r[2];
  float rot;
};


std::vector<int> linearToGridIdx(int linearIdx, const std::vector<int> & arrSize);

int linearIdx(int x, int y, const std::vector<int> & arrSize);

int gridToLinearIdx(const std::vector<int> & coord, const std::vector<int> & arrSize);

bool inside(const std::vector<int> & coord,
  const std::vector<int> & gridSize);

bool ptInRect(const float * pt, const float * len, int dim);

bool ptInEllipse(const float * pt, const float * r, int dim);

/// \param r. center and length are specified relative to a unit square.
/// The array represents a unit square with its origin at lower left corner and side lengths scaled to 1.
/// Area outside the unit square is ignored.
/// \param arr represents a 2D array.
/// \param arrSize x and y size of array.
/// \param val Area covered by the rectangle is written this value.
template <typename T>
void drawRectangle(const Rectangle2D & r, std::vector<T> & arr,
  const std::vector<int> & arrSize, T val)
{
  int dim = 2;
  //lower and upper bounds of a loop.
  int l[2];
  int u[2];
  float maxLen = std::sqrt(r.length[0] * r.length[0] + r.length[1] * r.length[1]);
  for (int i = 0; i < 2; i++){
    l[i] = (int)((r.center[i] - 0.5 * maxLen) * arrSize[i] + 0.5 );
    u[i] = (int)((r.center[i] + 0.5 * maxLen) * arrSize[i] );
  }
  Eigen::Rotation2D<float> R(r.rot);
  Eigen::Vector2f recCenter;
  recCenter << r.center[0] , r.center[1];
  R = R.inverse();
  std::vector<int> coord(2);
  for (int i = l[0]; i < u[0]; i++){
    for (int j = l[1]; j < u[1]; j++){
      
      coord[0] = i;
      coord[1] = j;
      if (!inside(coord, arrSize)){
        continue;
      }
      Eigen::Vector2f pixel;
      pixel << (i + 0.5f)/arrSize[0], (j + 0.5f)/arrSize[1];
      pixel -= recCenter;
      pixel = R * pixel;
      if (!ptInRect(pixel.data(), r.length, 2)){
        continue;
      }

      int linearIdx = gridToLinearIdx(coord, arrSize);
      arr[linearIdx] = val;
    }
  }
}

template <typename T>
void drawEllipse(const Ellipse2D & e, std::vector<T> & arr,
  const std::vector<int> & arrSize, T val)
{
  int dim = 2;
  //lower and upper bounds of a loop.
  int l[2];
  int u[2];
  float maxLen = std::max(e.r[0], e.r[1]);
  for (int i = 0; i < 2; i++){
    l[i] = (int)((e.center[i] - maxLen) * arrSize[i] + 0.5);
    u[i] = (int)((e.center[i] + maxLen) * arrSize[i]);
  }
  Eigen::Rotation2D<float> R(e.rot);
  Eigen::Vector2f eCenter;
  eCenter << e.center[0], e.center[1];
  R = R.inverse();
  std::vector<int> coord(2);
  for (int i = l[0]; i < u[0]; i++){
    for (int j = l[1]; j < u[1]; j++){

      coord[0] = i;
      coord[1] = j;
      if (!inside(coord, arrSize)){
        continue;
      }
      Eigen::Vector2f pixel;
      pixel << (i + 0.5f) / arrSize[0], (j + 0.5f) / arrSize[1];
      pixel -= eCenter;
      pixel = R * pixel;
      if (!ptInEllipse(pixel.data(), e.r, 2)){
        continue;
      }

      int linearIdx = gridToLinearIdx(coord, arrSize);
      arr[linearIdx] = val;
    }
  }
}

//copy the lower left corner with x>=y to all 8 triangles.
template <typename T>
void make2DCubic(std::vector<T> & arr, const std::vector<int> & arrSize)
{
  int mid[2] = { arrSize[0] / 2, arrSize[1] / 2 };
  for (int i = 0; i < arrSize[0]; i++){
    for (int j = 0; j < arrSize[1]; j++){
      int x = i;
      if (x >= mid[0]){
        x = arrSize[0] - x - 1;
      }
      int y = j;
      if (y>mid[1]){
        y = arrSize[1] - y - 1;
      }
      if (y > x){
        int tmp = y;
        y = x;
        x = tmp;
      }
      std::vector<int> dstCoord(2), srcCoord(2);
      dstCoord[0] = i;
      dstCoord[1] = j;
      srcCoord[0] = x;
      srcCoord[1] = y;
      int dstIdx = gridToLinearIdx(dstCoord, arrSize);
      int srcIdx = gridToLinearIdx(srcCoord, arrSize);
      if (srcIdx == dstIdx){
        continue;
      }
      arr[dstIdx] = arr[srcIdx];
    }
  }
}

template <typename T>
void saveStructure2D(const std::vector<T> & arr,
  const std::vector<int> & arrSize, std::ostream & out)
{
  std::vector<int> coord(2);
  out << arrSize[0] << " " << arrSize[1] << "\n";;
  for (int i = 0; i < arrSize[0]; i++){
    for (int j = 0; j < arrSize[1]; j++){
      coord[0] = i;
      coord[1] = j;
      int linearIdx = gridToLinearIdx(coord, arrSize);
      out << arr[linearIdx] << " ";
    }
    out << "\n";
  }
}

template <typename T>
void loadStructure2D(std::vector<T> & arr,
  std::vector<int> & arrSize, std::istream & in)
{
  arrSize.resize(2);
  in >> arrSize[0] >> arrSize[1];
  arr.resize(arrSize[0] * arrSize[1]);
  std::vector<int> coord(2);
  for (int i = 0; i < arrSize[0]; i++){
    for (int j = 0; j < arrSize[1]; j++){
      coord[0] = i;
      coord[1] = j;
      int linearIdx = gridToLinearIdx(coord, arrSize);
      in >> arr[linearIdx];
    }
  }
}

template <typename T>
void resize2D(const std::vector<T> & s, const std::vector<int> & size0,
  std::vector<T> & s1, const std::vector<int> & size1)
{
  int nEle = size1[0] * size1[1];
  s1.resize(nEle);
  for (int i = 0; i < size1[0]; i++){
    for (int j = 0; j < size1[1]; j++){
      int i0 = (int)((i + 0.5) / size1[0] * size0[0]);
      int j0 = (int)((j + 0.5) / size1[1] * size0[1]);
      int l0 = linearIdx(i0, j0, size0);
      int l1 = linearIdx(i, j, size1);
      s1[l1] = s[l0];
    }
  }
}

void test2DParam();
#endif
