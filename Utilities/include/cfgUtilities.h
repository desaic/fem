#ifndef cfgUtilities_h
#define cfgUtilities_h

#include <iostream>
#include <numeric>
#include <set>
#include <vector>
#include <map>
#include <fstream>

namespace cfgUtil
{
  /*-------------------------------------------------------------------------------------*/
  // Container manipulations
  /*-------------------------------------------------------------------------------------*/
  template<class T1, class T2>
  std::vector<T2> convertVec(const std::vector<T1> &iVec);

  template<class T1, class T2>
  std::vector<std::vector<T2> > convertVec(const std::vector<std::vector<T1> > &iVec);

  template<class T1, class T2>
  T2 * createArray(const std::vector<T1> &iVec);

  template<class T>
  std::vector<T> getSubVector(const std::vector<T> &iVec, const std::vector<int> &iIndices);

  template<class T>
  std::vector<T> getSubVector(const std::vector<T> &iVec, int iDim, const std::vector<int> &iIndices);

  template<class T>
  std::vector<std::vector<T> > getSubVector(const std::vector<std::vector<T> > &iVec, const std::vector<int> &iIndices);

  template<class T>
  std::vector<T> toStdVector(const std::set<T> &iSet);

  template<class T>
  std::set<T> toStdSet(const std::vector<T> &iVec);

  template<class T>
  std::map<std::vector<T>, int> computeVector2IndexMap(const std::vector<std::vector<T> > &iVec);

  /*-------------------------------------------------------------------------------------*/
  // numerics
  /*-------------------------------------------------------------------------------------*/
  template<class T>
  T average(const std::vector<T> &iVec);

  template<class T, int Dim>
  void getBarycenter(const std::vector<T> &iVec, T oBarycenter[Dim]);

  //add iValue to every element in iVec
  template<class T>
  std::vector<T> add(const std::vector<T> &iVec, const T &iValue); 

  // return iVec1+iVec2
  template<class T>
  std::vector<T> add(const std::vector<T> &iVec1, const std::vector<T> &iVec2); 

  //substract iValue to every element in iVec
  template<class T>
  std::vector<T> sub(const std::vector<T> &iVec, const T &iValue); 

  // return iVec1-iVec2
  template<class T>
  std::vector<T> sub(const std::vector<T> &iVec1, const std::vector<T> &iVec2); 

  template<class T>
  std::vector<T> mult(const std::vector<T> &iVec, const T &iValue); 

  template<class T, int Dim>
  std::vector<T> mult(const std::vector<T> &iVec, const T iValue[Dim]); 

  template<class T>
  T innerProd(const std::vector<T> &iVec1, const std::vector<T> &iVec2);

  template<class T>
  void fitLinearFunction(const std::vector<T> &iX, const std::vector<T> &iY, T &oSlope);

  /*-------------------------------------------------------------------------------------*/
  // I/O operations
  /*-------------------------------------------------------------------------------------*/

  // matlab
  template<class T>
  void writeVector2File(const std::vector<int> &iVector, const std::string &iFileName);

  template<class T>
  void serialize(std::ostream &ioStream, const std::vector<T> &iVec, const std::string &iString="");

  template<class T, int n>
  void serialize(std::ostream &ioStream, const std::vector<T> iVec[n], const std::string &iString="");

  template<class T>
  void serialize(std::ostream &ioStream, const std::vector<std::vector<T> > &iVec, const std::string &iString="");

  template<class T>
  void serialize(std::ostream &ioStream, const T *iArray, int iArraySize, const std::string &iString="");

  template<class T>
  void deSerialize(std::istream &iStream, std::vector<T> &oVec, const std::string &iString="");

  template<class T, int n>
  void deSerialize(std::istream &iStream, std::vector<T> oVec[n], const std::string &iString="");

  template<class T>
  void deSerialize(std::istream &iStream, std::vector<std::vector<T> > &oVec, const std::string &iString="");

  template<class T>
  void deSerialize(std::istream &iStream, int iArraySize, T *oArray, const std::string &iString="");

  template<class T>
  bool writeBinary(const std::string &iFileName, const std::vector<T> &iVec);

  template<class T>
  bool readBinary(const std::string &iFileName, std::vector<T> &oVec);

  template<class T>
  bool writeBinary(const std::string &iFileName, const std::vector<std::vector<T> > &iVec);

  template<class T>
  bool readBinary(const std::string &iFileName, std::vector<std::vector<T> > &oVec);

  /*=====================================================================================*/
  // Implementations
  /*=====================================================================================*/

  template<class T1, class T2>
  std::vector<T2> convertVec(const std::vector<T1> &iVec)
  {
    std::vector<T2> newVec;
    newVec.reserve(iVec.size());
    size_t ielem, nelem=iVec.size();
    for (ielem=0; ielem<nelem; ielem++)
    {
      newVec.push_back(iVec[ielem]);
    }
    return newVec;
  }

  template<class T1, class T2>
  std::vector<std::vector<T2> > convertVec(const std::vector<std::vector<T1> > &iVec)
  {
    std::vector<std::vector<T2> > newVec;
    newVec.reserve(iVec.size());
    size_t ivec, nvec=iVec.size();
    for (ivec=0; ivec<nvec; ivec++)
    {
      newVec.push_back(convertVec<T1,T2>(iVec[ivec]));
    }
    return newVec;
  }

  template<class T1, class T2>
  T2 * createArray(const std::vector<T1> &iVec)
  {
    T2 * newArray = new T2[iVec.size()];
    size_t ielem, nelem=iVec.size();
    for (ielem=0; ielem<nelem; ielem++)
    {
      newArray[ielem] = iVec[ielem];
    }
    return newArray;
  }

  template<class T>
  std::vector<T> getSubVector(const std::vector<T> &iVec, const std::vector<int> &iIndices)
  {
    std::vector<T> subVector;
    size_t ielem, nelem=iIndices.size();
    for (ielem=0; ielem<nelem; ielem++)
    {
      int ind = iIndices[ielem];
      assert(ind>=0 && ind<iVec.size());
      subVector.push_back(iVec[ind]);
    }
    return subVector;
  }

  template<class T>
  std::vector<T> getSubVector(const std::vector<T> &iVec, int iDim, const std::vector<int> &iIndices)
  {
    assert(iVec.size()%iDim==0 && iDim>0);

    std::vector<T> subVector;
    size_t ielem, nelem=iIndices.size();
    for (ielem=0; ielem<nelem; ielem++)
    {
      int ind = iIndices[ielem];
      assert(ind>=0 && ind<iVec.size()/iDim);
      int icoord;
      for (icoord=0; icoord<iDim; icoord++)
      {
        subVector.push_back(iVec[iDim*ind+icoord]);
      }
    }
    return subVector;
  }

  template<class T>
  std::vector<std::vector<T> > getSubVector(const std::vector<std::vector<T> > &iVec, const std::vector<int> &iIndices)
  {
    size_t ielem, nelem=iIndices.size();
    std::vector<std::vector<T> > subVector(nelem);
    for (ielem=0; ielem<nelem; ielem++)
    {
      int ind = iIndices[ielem];
      assert(ind>=0 && ind<iVec.size());
      subVector[ielem].insert(subVector[ielem].end(), iVec[ind].begin(), iVec[ind].end());
    }
    return subVector;
  }

  template<class T>
  std::vector<T> toStdVector(const std::set<T> &iSet)
  {
    std::vector<T> vec;
    typename std::set<T>::const_iterator it, it_end=iSet.end();
    for (it=iSet.begin(); it!=it_end; it++)
    {
      vec.push_back(*it);
    }
    return vec;
  }

  template<class T>
  std::set<T> toStdSet(const std::vector<T> &iVec)
  {
    std::set<T> Set;
    typename std::vector<T>::const_iterator it, it_end=iVec.end();
    for (it=iVec.begin(); it!=it_end; it++)
    {
      Set.insert(*it);
    }
    return Set;
  }

  template<class T>
  std::map<std::vector<T>, int> computeVector2IndexMap(const std::vector<std::vector<T> > &iVec)
  {
    std::map<std::vector<T>, int> mapVec2Index;
    int ielem, nelem=(int)iVec.size();
    for (ielem=0; ielem<nelem; ielem++)
    {
      mapVec2Index[iVec[ielem]] = ielem;
    }
    return mapVec2Index;
  }

  template<class T>
  T average(const std::vector<T> &iVec)
  {
    T avg = std::accumulate(iVec.begin(), iVec.end(), (T)0);
    if (iVec.size()>0)
      avg /= iVec.size();

	return avg;
  }

  template<class T, int Dim>
  void getBarycenter(const std::vector<T> &iVec,  T oBarycenter[Dim])
  {
    assert(iVec.size()%Dim==0 && Dim>0);
    int ipoint, npoint=(int)iVec.size()/Dim;
    int icoord;
    for (icoord=0; icoord<Dim; icoord++)
    {
      oBarycenter[icoord]=0;
      for (ipoint=0; ipoint<npoint; ipoint++)
      {
        oBarycenter[icoord] += iVec[Dim*ipoint+icoord];
      }
    }
    if (npoint>0)
    {
      for (icoord=0; icoord<Dim; icoord++)
      {
        oBarycenter[icoord] /= npoint;
      }
    }
  }

  template<class T>
  std::vector<T> add(const std::vector<T> &iVec, const T &iValue)
  {
    std::vector<T> vec = iVec;
    size_t ielem, nelem=iVec.size();
    for (ielem=0; ielem<nelem; ielem++)
    {
      vec[ielem] += iValue;
    }
    return vec;
  }

  template<class T>
  std::vector<T> add(const std::vector<T> &iVec1, const std::vector<T> &iVec2)
  {
    assert(iVec1.size()==iVec2.size());
    std::vector<T> vec = iVec1;
    size_t ielem, nelem=iVec1.size();
    for (ielem=0; ielem<nelem; ielem++)
    {
      vec[ielem] += iVec2[ielem];
    }
    return vec;
  }

  template<class T>
  std::vector<T> sub(const std::vector<T> &iVec, const T &iValue)
  {
    std::vector<T> vec = iVec;
    size_t ielem, nelem=iVec.size();
    for (ielem=0; ielem<nelem; ielem++)
    {
      vec[ielem] -= iValue;
    }
    return vec;
  }

  template<class T>
  std::vector<T> sub(const std::vector<T> &iVec1, const std::vector<T> &iVec2)
  {
    assert(iVec1.size()==iVec2.size());
    std::vector<T> vec = iVec1;
    size_t ielem, nelem=iVec1.size();
    for (ielem=0; ielem<nelem; ielem++)
    {
      vec[ielem] -= iVec2[ielem];
    }
    return vec;
  }


  template<class T>
  std::vector<T> mult(const std::vector<T> &iVec, const T &iValue)
  {
    std::vector<T> vec = iVec;
    size_t ielem, nelem=iVec.size();
    for (ielem=0; ielem<nelem; ielem++)
    {
      vec[ielem] *= iValue;
    }
    return vec;
  }

  template<class T, int Dim>
  std::vector<T> mult(const std::vector<T> &iVec, const T iValue[Dim])
  {
    assert(iVec.size()%Dim==0 && Dim>0);

    std::vector<T> vec = iVec;
    int ipoint, npoint=(int)iVec.size()/Dim;
    int icoord;
    for (icoord=0; icoord<Dim; icoord++)
    {
      for (ipoint=0; ipoint<npoint; ipoint++)
      {
        vec[Dim*ipoint+icoord] *= iValue[Dim];
      }
    }
    return vec;
  }

  template<class T>
  T innerProd(const std::vector<T> &iVec1, const std::vector<T> &iVec2)
  {
    T val = std::inner_product(iVec1.begin(),iVec1.end(), iVec2.begin(), (T)0);
    return val;
  }

  template<class T>
  void fitLinearFunction(const std::vector<T> &iX, const std::vector<T> &iY, T &oSlope)
  {
    T avg_X = average<T>(iX);
    T avg_Y = average<T>(iY);
    std::vector<T> X_minus_Xavg = sub(iX, avg_X);
    std::vector<T> Y_minus_Yavg = sub(iY, avg_Y);

    oSlope = innerProd<T>(X_minus_Xavg, Y_minus_Yavg) / innerProd<T>(X_minus_Xavg, X_minus_Xavg);
  }

  template<class T>
  void writeVector2File(const std::vector<T> &iVector, const std::string &iFileName)
  {
    std::ofstream stream(iFileName);
    assert(stream.is_open());

    size_t ielem, nelem=iVector.size();
    for (ielem=0; ielem<nelem; ielem++)
    {
      stream << iVector[ielem] << " ";
    }
    stream << std::endl;
  }

  template<class T>
  void serialize(std::ostream &ioStream, const std::vector<T> &iVec, const std::string &iString)
  {
    size_t ielem, nelem=iVec.size();
    ioStream << nelem << " ";
    for (ielem=0; ielem<nelem; ielem++)
    {
      ioStream << iVec[ielem] << " "; 
    }
    ioStream << std::endl;
  }

  template<class T, int n>
  void serialize(std::ostream &oStream, const std::vector<T> iVec[n], const std::string &iString)
  {
    size_t ielem, nelem=iVec.size();
    oStream << nelem << " ";
    for (ielem=0; ielem<nelem; ielem++)
    {
      int  icol;
      for (icol=0; icol<n; icol++)
      {
        oStream << iVec[ielem][icol] << " "; 
      }
      oStream << std::endl;
    }
  }

  template<class T>
  void serialize(std::ostream &oStream, const std::vector<std::vector<T> > &iVec, const std::string &iString)
  {
    size_t ielem, nelem=iVec.size();
    oStream << nelem << " " << std::endl;
    for (ielem=0; ielem<nelem; ielem++)
    {
      serializeVec(oStream, iVec[ielem], "");
    }
  }

  template<class T>
  void serialize(std::ostream &oStream, const T *iArray, int iArraySize, const std::string &iString)
  {
    int ielem;
    oStream << iArraySize << " ";
    for (ielem=0; ielem<iArraySize; ielem++)
    {
      oStream << iArray[ielem] << " "; 
    }
    oStream << std::endl;
  }

  template<class T>
  void deSerialize(std::istream &iStream, std::vector<T> &oVec, const std::string &iString)
  {
    if (iString!="")
    {
      std::string dummy;
      iStream >> dummy;
    }
    size_t ielem, nelem=0;
    iStream >> nelem;
	  oVec.resize(nelem);
    for (ielem=0; ielem<nelem; ielem++)
    {
      iStream >> oVec[ielem]; 
    }
  }

  template<class T, int n>
  void deSerialize(std::istream &iStream, std::vector<T> oVec[n], const std::string &iString)
  {
    if (iString!="")
    {
      std::string dummy;
      iStream >> dummy;
    }
    size_t ielem, nelem=0;
    iStream >> nelem;
    for (ielem=0; ielem<nelem; ielem++)
    {
      T elemArray[n];
      deSerialize(iStream, n, elemArray);
      int ival;
      for (ival=0; ival<n; ival++)
      {
        oVec[ival].push_back(elemArray[ival]);
      }
    }
  }

  template<class T>
  void deSerialize(std::istream &iStream, std::vector<std::vector<T> > &oVec, const std::string &iString)
  {
    if (iString!="")
    {
      std::string dummy;
      iStream >> dummy;
    }
    size_t ielem, nelem=0;
    iStream >> nelem;
    oVec.resize(nelem);
    for (ielem=0; ielem<nelem; ielem++)
    {
      deSerialize(iStream, oVec[ielem]);
    }
  }

  template<class T>
  void deSerialize(std::istream &iStream, int iArraySize, T *oArray, const std::string &iString)
  {
    if (iString!="")
    {
      std::string dummy;
      iStream >> dummy;
    }
    int ielem=0;
    for (ielem=0; ielem<iArraySize; ielem++)
    {
      iStream >> oArray[ielem]; 
    }
  }

  template<class T>
  bool writeBinary(const std::string &iFileName, const std::vector<T> &iVec)
  {
    std::fstream file(iFileName, std::ios::out | std::ios::binary);
    int nelem = (int)iVec.size();
    file.write((char*)&nelem, sizeof(int));
    file.write((char*)&iVec[0], nelem*sizeof(T));
    file.close();
    return true;
  }

  template<class T>
  bool readBinary(const std::string &iFileName, std::vector<T> &oVec)
  {
    std::ifstream file (iFileName, std::ios::in | std::ios::binary);
    if (!file.is_open())
      return false;

    int nelem;
    file.read((char*)&nelem, sizeof(int));
    oVec.resize(nelem);
    file.read((char*)&oVec[0], nelem*sizeof(T));
    file.close();
    return true;
  }

  template<class T>
  bool writeBinary(const std::string &iFileName, const std::vector<std::vector<T> > &iVec)
  {
    std::fstream file(iFileName, std::ios::out | std::ios::binary);
    int nvec = (int)iVec.size();
    file.write((char*)&nvec, sizeof(int));
    int ivec;
    for (ivec=0; ivec<nvec; ivec++)
    {
      int nelem = (int)iVec[ivec].size();
      file.write((char*)&nelem, sizeof(int));
      file.write((char*)&iVec[ivec][0], nelem*sizeof(T));
    }
    file.close();
    return true;
  }

  template<class T>
  bool readBinary(const std::string &iFileName, std::vector<std::vector<T> > &oVec)
  {
    std::ifstream file (iFileName, std::ios::in | std::ios::binary);
    if (!file.is_open())
      return false;

    int nvec;
    file.read((char*)&nvec, sizeof(int));
    oVec.resize(nvec);
    int ivec;
    for (ivec=0; ivec<nvec; ivec++)
    {
      int nelem;
      file.read((char*)&nelem, sizeof(int));
      oVec[ivec].resize(nelem);
      file.read((char*)&oVec[ivec][0], nelem*sizeof(T));
      bool toto = true;
    }
    file.close();
    return true;
  }

} //namespace cfgUtil

#endif 




