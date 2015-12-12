#include "cfgMaterialUtilities.h"

#include <iostream>
#include <fstream>
#include <assert.h>

#include "ElementRegGrid.hpp"
#include "Element.hpp"
#include "ElementRegGrid2D.h"
#include "Element2D.h"

#include "qhullUtilities.h"

#include <Qhull.h>
#include <QhullFacetList.h>
#include <QhullVertex.h>
#include <QhullVertexSet.h>
using namespace orgQhull;

#include <cfgUtilities.h>
using namespace cfgUtil;

#include "DistanceTool.h"

#include <list>

bool cfgMaterialUtilities::readMaterialCombinations(const std::string iFileName, std::vector<std::vector<int> > &oMaterials)
{
  std::ifstream stream(iFileName);
  if (!stream.is_open())
  {
    return false;
  }

  deSerialize<int>(stream, oMaterials, "Materials");
  return true;
}


bool cfgMaterialUtilities::writeMaterialCombinations(const std::string iFileName, const std::vector<std::vector<int> > &iMaterials)
{
  std::ofstream stream(iFileName);
  if (!stream.is_open())
  {
    std::cout << "cfgMaterialUtilities::writeMaterialCombinations: invalid fileName" << std::endl;
    return false;
  }
   
  stream << "Materials" << std::endl;
  int ivec=0, nvec=(int)iMaterials.size();
  stream << nvec << " " << std::endl;
  for (ivec=0; ivec<nvec; ivec++)
  {
    const std::vector<int> &iMaterialVec = iMaterials[ivec];
    int imat=0, nmat=(int)iMaterialVec.size();
    stream << nmat << " ";
    for (imat=0; imat<nmat; imat++)
    {
      stream << iMaterialVec[imat] << " ";
    }
    stream << std::endl;
  }
  return true;
}

bool cfgMaterialUtilities::readData(const std::string &iFileName,  Vector3f &oForceAxis, std::vector<int> &oMaterials, std::vector<float> &oStresses, std::vector<float> &oStrains)
{
  std::ifstream stream(iFileName);
  if (!stream.is_open())
  {
    return false;
  }

  std::string dummy;
  // Force Axis
  stream >> dummy; 
  stream >> oForceAxis[0] >>  oForceAxis[1] >> oForceAxis[2];
  
  // Materials
  deSerialize<int>(stream, oMaterials,"Materials");

  // Stress strain
  std::vector<float> StrainStress[2];
  deSerialize<float,2>(stream, StrainStress, "StressStrain");
  oStrains = StrainStress[0];
  oStresses = StrainStress[1];

  return true;
}

bool cfgMaterialUtilities::writeData(const std::string &iFileName, const Vector3f &iForceAxis, const std::vector<int> &iMaterials, const std::vector<float> &iStresses, const std::vector<float> &iStrains)
{
  std::ofstream stream(iFileName);
  if (!stream.is_open())
  {
    return false;
  }

  stream << "ForceAxis" << std::endl; 
  stream << iForceAxis[0] << " " << iForceAxis[1] << " " << iForceAxis[2]  << std::endl;
  
  stream << "Materials" << std::endl;
  int imat=0, nmat=(int)iMaterials.size();
  stream << nmat << " ";
  for (imat=0; imat<nmat; imat++)
  {
    stream << iMaterials[imat] << " ";
  }
  stream << std::endl;

  stream << "StressStrain" << std::endl;
  int isample, nsample=(int)iStresses.size();
  stream << nsample+1 << " " << std::endl;
  stream << 0 << " " << 0 << std::endl; 
  for (isample=0; isample<nsample; isample++)
  {
    stream << iStrains[isample] << " " << iStresses[isample]  << std::endl; 
  }
  return false;
}

bool cfgMaterialUtilities::writeData(const std::string &iFileName, const Vector3f &iForceAxis, const std::vector<std::vector<int> > &iMaterials, const std::vector<std::vector<float> > &iStresses, const std::vector<std::vector<float> > &iStrains)
{
  std::ofstream stream(iFileName);
  if (!stream.is_open())
  {
    return false;
  }
  stream << iForceAxis[0] << " " << iForceAxis[1] << " " << iForceAxis[2]  << std::endl;

  int ncomb = (int)iMaterials.size();
  assert(iStresses.size()==ncomb && iStrains.size()==ncomb);
  stream << ncomb << " " << std::endl;

  int icomb;
  for (icomb=0; icomb<ncomb; icomb++)
  {
    const std::vector<int> &materials = iMaterials[icomb];
    const std::vector<float> &stresses = iStresses[icomb];
    const std::vector<float> &strains = iStrains[icomb];

    cfgUtil::serialize<int>(stream, materials, "");

    int isample=0, nsample=(int)stresses.size();
    stream << nsample << std::endl;
    for (isample=0; isample<nsample; isample++)
    {
      stream << strains[isample] << " " << stresses[isample]  << std::endl; 
    }
  }
  return true;
}

bool cfgMaterialUtilities::writeData(const std::string &iFileName, const Vector3f &iForceAxis, const std::vector<std::vector<int> > &iMaterials, const std::vector<std::vector<float> > &iStresses, 
                                     const std::vector<std::vector<std::vector<float> > > &ix, int iN[3])
{
  std::ofstream stream(iFileName);
  if (!stream.is_open())
  {
    return false;
  }
  stream << iN[0] << " " << iN[1] << " " << iN[2]  << std::endl;

  stream << iForceAxis[0] << " " << iForceAxis[1] << " " << iForceAxis[2]  << std::endl;

  int ncomb = (int)iMaterials.size();
  stream << ncomb << " " << std::endl;

  int icomb;
  for (icomb=0; icomb<ncomb; icomb++)
  {
    const std::vector<int> &materials = iMaterials[icomb];
    const std::vector<float> &stresses = iStresses[icomb];
    const std::vector<std::vector<float> > &x = ix[icomb];

    cfgUtil::serialize<int>(stream, materials, "");

    int isample=0, nsample=(int)stresses.size();
    stream << nsample << std::endl;
    for (isample=0; isample<nsample; isample++)
    {
      stream << stresses[isample] << " "; 
      cfgUtil::serialize<float>(stream, x[isample], "");
    }
  }
  return true;
}

bool cfgMaterialUtilities::writeDataBinary(const std::string &iFileName, const Vector3f &iForceAxis, const std::vector<std::vector<int> > &iMaterials, const std::vector<std::vector<float> > &iStresses, 
                                           const std::vector<std::vector<std::vector<float> > > &ix, int iN[3])
{
  std::fstream file(iFileName, std::ios::out | std::ios::binary);
  if (!file.is_open())
  {
    return false;
  }

  file.write((char*)&iN[0], 3*sizeof(int));
  file.write((char*)&iForceAxis[0], 3*sizeof(float));
  int ncomb = (int)iMaterials.size();
  file.write((char*)&ncomb, sizeof(int));

  int icomb;
  for (icomb=0; icomb<ncomb; icomb++)
  {
    const std::vector<int> &materials = iMaterials[icomb];
    const std::vector<float> &stresses = iStresses[icomb];
    const std::vector<std::vector<float> > &x = ix[icomb];

    int nmat = (int)materials.size();
    file.write((char*)&nmat, sizeof(int));
    file.write((char*)&materials[0], nmat*sizeof(int));

    int isample=0, nsample=(int)stresses.size();
    file.write((char*)&nsample, sizeof(int));
    file.write((char*)&stresses[isample], nsample*sizeof(float));

    for (isample=0; isample<nsample; isample++)
    {
      int nx = (int)x[isample].size();
      file.write((char*)&nx, sizeof(int));
      file.write((char*)&x[isample][0], nx*sizeof(float));
    }
  }
  file.close();
  return true;
}

bool cfgMaterialUtilities::concatenateData(const std::string &iInputBaseFileName, int iNbFiles, std::string &iOutputFileName)
{
  std::string fileName = iInputBaseFileName + "_" + std::to_string(0) + ".txt";
  Vector3f forceAxis;
  std::vector<int> materials;
  std::vector<float> stresses, strains;
  bool ResOk = readData(fileName,  forceAxis, materials, stresses, strains);
  if (!ResOk)
    return false;

  std::ofstream stream(iOutputFileName);
  if (!stream.is_open())
  {
    return false;
  }
  stream << forceAxis[0] << " " << forceAxis[1] << " " << forceAxis[2]  << std::endl;
  stream << iNbFiles << " " << std::endl;

  int ifile;
  for (ifile=0; ifile<iNbFiles && ResOk; ifile++)
  {
    std::string fileName = iInputBaseFileName + "_" + std::to_string(ifile) + ".txt";
    Vector3f forceAxis;
    std::vector<int> materials;
    std::vector<float> stresses, strains;

    bool ResOk = readData(fileName,  forceAxis, materials, stresses, strains);
    if (ResOk)
    {
       cfgUtil::serialize<int>(stream, materials, "");
       int isample=0, nsample=(int)stresses.size();
       stream << nsample << std::endl;
       for (isample=0; isample<nsample; isample++)
       {
         stream << strains[isample] << " " << stresses[isample]  << std::endl; 
       }
    }
  }
  return ResOk;
}


bool cfgMaterialUtilities::readData(const std::string &iFileName,  Vector3f &oForceAxis, std::vector<std::vector<int> > &oMaterials, std::vector<std::vector<float> > &oStresses, std::vector<std::vector<float> > &oStrains)
{
  oMaterials.clear();
  oStresses.clear();
  oStrains.clear();

  std::ifstream stream(iFileName);
  if (!stream.is_open())
  {
    return false;
  }
  // Force Axis
  stream >> oForceAxis[0] >>  oForceAxis[1] >> oForceAxis[2];
  int nsample;
  stream >> nsample;

  oMaterials.resize(nsample);
  oStresses.resize(nsample);
  oStrains.resize(nsample);

  int isample=0;
  for (isample=0; isample<nsample; isample++)
  {
    // Materials
    deSerialize<int>(stream, oMaterials[isample],"");
    // Stress strain
    std::vector<float> StrainStress[2];
    deSerialize<float,2>(stream, StrainStress, "");
    oStrains[isample] = StrainStress[0];
    oStresses[isample] = StrainStress[1];
  }
  return true;
}

bool cfgMaterialUtilities::readData(const std::string &iFileName,  Vector3f &oForceAxis, std::vector<std::vector<int> > &oMaterials, std::vector<std::vector<float> > &oStresses, std::vector<std::vector<std::vector<float> > > &ox, int oN[3])
{
  oMaterials.clear();
  oStresses.clear();
  ox.clear();

  std::ifstream stream(iFileName);
  if (!stream.is_open())
  {
    return false;
  }
  // structure size
  stream >> oN[0] >> oN[1] >> oN[2];

  // Force Axis
  stream >> oForceAxis[0] >>  oForceAxis[1] >> oForceAxis[2];
  int ncomb;
  stream >> ncomb;
  std::cout << "Nb comb = " << ncomb << std::endl;

  oMaterials.resize(ncomb);
  oStresses.resize(ncomb);
  ox.resize(ncomb);

  int icomb=0;
  for (icomb=0; icomb<ncomb; icomb++)
  {
    // Materials
    deSerialize<int>(stream, oMaterials[icomb],"");
    // Stress strain

    int nsample;
    stream >> nsample;
    ox[icomb].resize(nsample);
    oStresses[icomb].resize(nsample);
    int isample;
    for (isample=0; isample<nsample; isample++)
    {
      stream >> oStresses[icomb][isample];
      deSerialize<float>(stream, ox[icomb][isample], "");
    }

    if (icomb%1000==0)
    {
      std::cout << " " << icomb;
    }
  }
  return true;
}

bool cfgMaterialUtilities::readDataBinary(const std::string &iFileName,  Vector3f &oForceAxis, std::vector<std::vector<int> > &oMaterials, std::vector<std::vector<float> > &oStresses, std::vector<std::vector<std::vector<float> > > &ox, int oN[3])
{
  oMaterials.clear();
  oStresses.clear();
  ox.clear();

  std::ifstream file (iFileName, std::ios::in | std::ios::binary);
  if (!file.is_open())
  {
    return false;
  }

  // structure size
  file.read((char*)&oN[0], 3*sizeof(int));

  // Force Axis
  file.read((char*)&oForceAxis[0], 3*sizeof(float));

  int ncomb;
  file.read((char*)&ncomb, sizeof(int));

  oMaterials.resize(ncomb);
  oStresses.resize(ncomb);
  ox.resize(ncomb);

  int icomb=0;
  for (icomb=0; icomb<ncomb; icomb++)
  {
    // Materials
    int nmat;
    file.read((char*)&nmat, sizeof(int));
    oMaterials[icomb].resize(nmat);
    file.read((char*)&oMaterials[icomb][0], nmat*sizeof(int));

    // Stress strain
    int nsample;
    file.read((char*)&nsample, sizeof(int));
    ox[icomb].resize(nsample);
    oStresses[icomb].resize(nsample);
    file.read((char*)&oStresses[icomb][0], nsample*sizeof(float));
    int isample;
    for (isample=0; isample<nsample; isample++)
    {
      int nx;
      file.read((char*)&nx, sizeof(int));
      ox[icomb][isample].resize(nx);
      file.read((char*)&ox[icomb][isample][0], nx*sizeof(float));
    }

    if (icomb%1000==0)
    {
      std::cout << " " << icomb;
    }
  }
  return true;
}

bool cfgMaterialUtilities::computeMaterialParameters(const std::string &iMaterialFile, const std::string iStressStrainFilesDirectories[2], const std::string iStressStrainBaseFileName, int iNbFiles,
                                                    std::vector<cfgScalar> &oPhysicalParameters, std::vector<std::vector<int> > &oBaseMaterialStructures, std::vector<std::vector<int> > &oMaterialAssignmentsOneCell)
{
  oBaseMaterialStructures.clear();
  readMaterialCombinations(iMaterialFile,   oBaseMaterialStructures);

  std::vector<std::vector<int> > &baseMaterialStructures = oBaseMaterialStructures;

  oPhysicalParameters.clear(); //YoungModulus, density;
  int  ncomb=iNbFiles;
  oMaterialAssignmentsOneCell.clear();
  oMaterialAssignmentsOneCell.resize(ncomb);

  std::vector<std::vector<int> > & materials = oMaterialAssignmentsOneCell; 
  std::vector<std::vector<cfgScalar> > stresses[2], strains[2];

  int iaxis, naxis = 2;
  for (iaxis=0; iaxis<naxis; iaxis++)
  {
    stresses[iaxis].resize(ncomb);
    strains[iaxis].resize(ncomb);
  }
  int icomb;
  for (icomb=0; icomb<ncomb; icomb++)
  {
    int iaxis;
    for (iaxis=0; iaxis<naxis; iaxis++)
    {
      std::string sampleFile = iStressStrainFilesDirectories[iaxis] + iStressStrainBaseFileName + "_" + std::to_string(icomb) + ".txt"; 
      Vector3f forceAxis;
      std::vector<float> stresses_f, strains_f;
      readData(sampleFile, forceAxis, materials[icomb], stresses_f, strains_f);
      stresses[iaxis][icomb] = convertVec<float, cfgScalar>(stresses_f);
      strains[iaxis][icomb] = convertVec<float, cfgScalar>(strains_f);
    }
  }
  computeMaterialParameters(materials, baseMaterialStructures, stresses, strains, oPhysicalParameters);
  return true;
}

void cfgMaterialUtilities::computeStrain(int n[2], const std::vector<std::vector<std::vector<cfgScalar> > > iX[2], std::vector<std::vector<cfgScalar> > oStrains[2][2])
{
  int icomb=0, ncomb=(int)iX[0].size();
  for (int iaxis=0; iaxis<2; iaxis++)
  {
    for (int jaxis=0; jaxis<2; jaxis++)
    {
      oStrains[iaxis][jaxis].resize(ncomb);
      for (icomb=0; icomb<ncomb; icomb++)
      {
        int isample=0, nsample=(int)iX[iaxis][icomb].size();
        for (isample=0; isample<nsample; isample++)
        {
          cfgScalar strain = computeStrain(iX[iaxis][icomb][isample], n[0], n[1], jaxis);
          oStrains[iaxis][jaxis][icomb].push_back(strain);
        }
      }
    }
  }
}

void cfgMaterialUtilities::computeStrain3D(int n[3], const std::vector<std::vector<std::vector<cfgScalar> > > iX[3], std::vector<std::vector<cfgScalar> > oStrains[3][3])
{
  int icomb=0, ncomb=(int)iX[0].size();
  for (int iaxis=0; iaxis<3; iaxis++)
  {
    for (int jaxis=0; jaxis<3; jaxis++)
    {
      oStrains[iaxis][jaxis].resize(ncomb);
      for (icomb=0; icomb<ncomb; icomb++)
      {
        int isample=0, nsample=(int)iX[iaxis][icomb].size();
        for (isample=0; isample<nsample; isample++)
        {
          cfgScalar strain = computeStrain(iX[iaxis][icomb][isample], n[0], n[1], n[2], jaxis);
          oStrains[iaxis][jaxis][icomb].push_back(strain);
        }
      }
    }
  }
}

bool cfgMaterialUtilities::computeMaterialParametersFromDeformations(const std::string &iMaterialFile, const std::string iStressDeformationFileNames[2], bool iBinary,
                                                                     std::vector<cfgScalar> &oPhysicalParameters, std::vector<std::vector<int> > &oBaseMaterialStructures, std::vector<std::vector<int> > &oMaterialAssignmentsOneCell)
{
  std::cout << "Reading " << iMaterialFile << "...";

  oBaseMaterialStructures.clear();
  readMaterialCombinations(iMaterialFile,   oBaseMaterialStructures);

  std::vector<std::vector<int> > &baseMaterialStructures = oBaseMaterialStructures;
  std::vector<std::vector<int> > materials; 
  std::vector<std::vector<cfgScalar> > stresses[2], strains[2][2];
  
  bool ResOk = true;
  for (int iaxis=0; iaxis<2 && ResOk; iaxis++)
  {
    Vector3f forceAxis;
    std::vector<std::vector<std::vector<float> > > x;
    std::vector<std::vector<float> > stresses_f;
    int n[3];
    if (iBinary)
    {
      ResOk = readDataBinary(iStressDeformationFileNames[iaxis],  forceAxis, oMaterialAssignmentsOneCell, stresses_f, x, n);
    }
    else
    {
      ResOk = readData(iStressDeformationFileNames[iaxis],  forceAxis, oMaterialAssignmentsOneCell, stresses_f, x, n);
    }
    stresses[iaxis] = convertVec<float, cfgScalar>(stresses_f);
    if (ResOk)
    {
      int icomb=0, ncomb=(int)x.size();
      int jaxis;
      for (jaxis=0; jaxis<2; jaxis++)
      {
        strains[iaxis][jaxis].resize(ncomb);
        for (icomb=0; icomb<ncomb; icomb++)
        {
          int isample=0, nsample=(int)x[icomb].size();
          for (isample=0; isample<nsample; isample++)
          {
            cfgScalar strain = computeStrain(convertVec<float, cfgScalar>(x[icomb][isample]), n[0], n[1], jaxis);
            strains[iaxis][jaxis][icomb].push_back(strain);
          }
        }
      }
    }
  }
  std::cout << " OK" << std::endl;

  if (ResOk)
  {
    computeMaterialParameters(oMaterialAssignmentsOneCell, baseMaterialStructures, stresses, strains, oPhysicalParameters);
  }
  return ResOk;
}

bool cfgMaterialUtilities::computeMaterialParameters(const std::string &iMaterialFile, const std::string iStressStrainFileNames[2],
                                                    std::vector<cfgScalar> &oPhysicalParameters, std::vector<std::vector<int> > &oBaseMaterialStructures, std::vector<std::vector<int> > &oMaterialAssignmentsOneCell)
{
  oBaseMaterialStructures.clear();
  readMaterialCombinations(iMaterialFile,   oBaseMaterialStructures);

  std::vector<std::vector<int> > &baseMaterialStructures = oBaseMaterialStructures;
  std::vector<std::vector<int> > &materials = oMaterialAssignmentsOneCell; 
  std::vector<std::vector<cfgScalar> > stresses[2], strains[2];

  bool ResOk = true;
  for (int iaxis=0; iaxis<2 && ResOk; iaxis++)
  {
    Vector3f forceAxis;
    std::vector<std::vector<float> > stresses_f, strains_f;

    ResOk = readData(iStressStrainFileNames[iaxis],  forceAxis, materials, stresses_f, strains_f);
    stresses[iaxis] = convertVec<float, cfgScalar>(stresses_f);
    strains[iaxis] = convertVec<float, cfgScalar>(strains_f);
  }

  int ncomb=(int)materials.size();
  oPhysicalParameters.clear(); //YoungModulus, density;
  oMaterialAssignmentsOneCell.clear();
  oMaterialAssignmentsOneCell.resize(ncomb);

  if (ResOk)
  {
    computeMaterialParameters(materials, baseMaterialStructures, stresses, strains, oPhysicalParameters);
  }
  return ResOk;
}

// iStresses = [sigma_xx sigma_yy sigma_zz sigma_xy sigma_yz sigma_zx]
cfgScalar cfgMaterialUtilities::computeVonMisesStress(std::vector<cfgScalar> &iStresses)
{
  assert(iStresses.size()==6);

  cfgScalar diffSigma[3];
  for (int iaxis=0; iaxis<3; iaxis++)
  {
    cfgScalar val = (iStresses[iaxis]-iStresses[(iaxis+1)%3]);
    diffSigma[iaxis] = val*val;
  }
  cfgScalar vonMisesStress = sqrt(0.5*( (diffSigma[0] + diffSigma[1] + diffSigma[2] + 6*(iStresses[3]+iStresses[4]+iStresses[5]) ) ) );
  return vonMisesStress;
}

//iStresses = [sigma_xx sigma_yy sigma_xy]
cfgScalar cfgMaterialUtilities::computeVonMisesStress2D(std::vector<cfgScalar> &iStresses)
{
  cfgScalar val = sqrt(iStresses[0]*iStresses[0] - iStresses[0]*iStresses[1] + iStresses[1]*iStresses[1] + 3*iStresses[2]);
  return val;
}

cfgScalar cfgMaterialUtilities::computeYoungModulus(cfgScalar iStrain, cfgScalar iStress)
{
  cfgScalar YoungModulus=iStress/iStrain;
  return YoungModulus;
}

cfgScalar cfgMaterialUtilities::computeYoungModulus(const std::vector<cfgScalar> &iStrains, const std::vector<cfgScalar>  &iStresses)
{
  cfgScalar YoungModulus=0;
  if (iStrains.size()==1)
  {
    std::vector<cfgScalar> strains(1, 0), stresses(1, 0);
    strains.push_back(iStrains[0]);
    stresses.push_back(iStresses[0]);
    fitLinearFunction<cfgScalar>(strains, stresses, YoungModulus);
  }
  else
  {
    fitLinearFunction<cfgScalar>(iStrains, iStresses, YoungModulus);
  }
  return YoungModulus;
}

cfgScalar cfgMaterialUtilities::computePoissonRatio(const std::vector<cfgScalar> &iStrainsAxialDir, const std::vector<cfgScalar> &iStrainsTransversalDir)
{
  cfgScalar poissonRatio = 0;
  if (iStrainsAxialDir.size()>0 && iStrainsTransversalDir.size()>0)
  {
    poissonRatio = - iStrainsTransversalDir[0]/iStrainsAxialDir[0];

    /*cfgScalar lambda1 = 1 + iStrainsAxialDir[0];
    cfgScalar lambda2 = 1 + iStrainsTransversalDir[0];

    cfgScalar Vol = lambda1*lambda2;
    cfgScalar a = Vol/lambda1;
    poissonRatio = (1-sqrt(a))/(lambda1-1);*/ 
  }
  return poissonRatio;
}

cfgScalar cfgMaterialUtilities::computePoissonRatio(cfgScalar &iStrainsAxialDir, cfgScalar&iStrainsTransversalDir)
{
  cfgScalar poissonRatio = - iStrainsTransversalDir/iStrainsAxialDir;
  return poissonRatio;
}

cfgScalar cfgMaterialUtilities::computeMaterialRatio(const std::vector<int> &iMaterialAssignment, const std::vector<std::vector<int> > &iBaseMaterialStructures)
{
  cfgScalar r=0;
  int imat, nmat=(int)iMaterialAssignment.size();
  for (imat=0; imat<nmat; imat++)
  {
    int matIndex = iMaterialAssignment[imat];
    const std::vector<int> & baseMaterialStructure = iBaseMaterialStructures[matIndex];
    int icell, ncell=(int)baseMaterialStructure.size();
    for (icell=0; icell<ncell; icell++)
    {
      r += baseMaterialStructure[icell];
    }
  }
  int ncell = (int)iBaseMaterialStructures[0].size()*nmat;
  r /= (cfgScalar)ncell;
  return r;
}

bool cfgMaterialUtilities::computeMaterialParameters(const std::vector<std::vector<int> > &iMaterials, const std::vector<std::vector<int> > &iBaseMaterialStructures,
                                                     const std::vector<std::vector<cfgScalar> > iStresses[2], const std::vector<std::vector<cfgScalar> > iStrains[2],
                                                     std::vector<cfgScalar> &oPhysicalParameters)
{
  int naxis = 2;
  int icomb, ncomb = (int)iMaterials.size();
  for (icomb=0; icomb<ncomb; icomb++)
  {
    cfgScalar r = computeMaterialRatio(iMaterials[icomb], iBaseMaterialStructures);
    oPhysicalParameters.push_back(r);

    int iaxis;
    for (iaxis=0; iaxis<naxis; iaxis++)
    {
      cfgScalar Y = computeYoungModulus(iStrains[iaxis][icomb], iStresses[iaxis][icomb]);
      oPhysicalParameters.push_back(Y);
    }
  }
  return true;
}

bool cfgMaterialUtilities::computeMaterialParameters(const std::vector<std::vector<int> > &iMaterials, const std::vector<std::vector<int> > &iBaseMaterialStructures,
                                                     const std::vector<std::vector<cfgScalar> > iStresses[2], const std::vector<std::vector<cfgScalar> > iStrains[2][2],
                                                     std::vector<cfgScalar> &oPhysicalParameters)
{
  oPhysicalParameters.clear();

  int iVersion=2; // 0: Density Yx Yy; 1: Density Nux Nuy 2: Density Yx Yy, Nux Nuy

  int naxis = 2;
  int icomb, ncomb = (int)iMaterials.size();
  for (icomb=0; icomb<ncomb; icomb++)
  {
    cfgScalar r = computeMaterialRatio(iMaterials[icomb], iBaseMaterialStructures);
    oPhysicalParameters.push_back(r);

    if (iVersion==0 || iVersion==2)
    {
      int iaxis;
      for (iaxis=0; iaxis<naxis; iaxis++)
      {
        cfgScalar Y = computeYoungModulus(iStrains[iaxis][iaxis][icomb], iStresses[iaxis][icomb]);
        oPhysicalParameters.push_back(Y);
      }
    }
    if (iVersion==1 || iVersion==2)
    {
      int iaxis;
      for (iaxis=0; iaxis<naxis; iaxis++)
      {
        cfgScalar poissonRatio = computePoissonRatio(iStrains[iaxis][iaxis][icomb], iStrains[iaxis][1-iaxis][icomb]); 
        oPhysicalParameters.push_back(poissonRatio);
      }
    }
  }
  return true;
}

bool cfgMaterialUtilities::computeMaterialParameters(const std::vector<std::vector<int> > &iMaterials, const std::vector<std::vector<int> > &iBaseMaterialStructures,
                                                     const std::vector<std::vector<cfgScalar> > iStresses[3], const std::vector<std::vector<cfgScalar> > iStrains[3][3],
                                                     std::vector<cfgScalar> &oPhysicalParameters)
{
  oPhysicalParameters.clear();

  int iVersion=2; // 0: Density Yx Yy; 1: Density Nux Nuy, 2: Density Yx Yy Nux Nuy

  int naxis = 3;
  int icomb, ncomb = (int)iMaterials.size();
  for (icomb=0; icomb<ncomb; icomb++)
  {
    cfgScalar r = computeMaterialRatio(iMaterials[icomb], iBaseMaterialStructures);
    oPhysicalParameters.push_back(r);

    if (iVersion==0 || iVersion==2)
    {
      int iaxis;
      for (iaxis=0; iaxis<naxis; iaxis++)
      {
        cfgScalar Y = computeYoungModulus(iStrains[iaxis][iaxis][icomb], iStresses[iaxis][icomb]);
        oPhysicalParameters.push_back(Y);
      }
    }
    if (iVersion==1 || iVersion==2)
    {
      cfgScalar poissonRatio_xy = computePoissonRatio(iStrains[0][0][icomb], iStrains[0][1][icomb]); 
      oPhysicalParameters.push_back(poissonRatio_xy);

      cfgScalar poissonRatio_xz = computePoissonRatio(iStrains[0][0][icomb], iStrains[0][2][icomb]); 
      oPhysicalParameters.push_back(poissonRatio_xz);

      cfgScalar poissonRatio_yz = computePoissonRatio(iStrains[1][1][icomb], iStrains[1][2][icomb]); 
      oPhysicalParameters.push_back(poissonRatio_yz);
    }
  }
  return true;
}

bool cfgMaterialUtilities::getMaterialAssignment(int iDim, std::string iMaterialFile, const std::vector<int> &iMaterialCombination, int iBlockRep, int iNbSubdiv, std::vector<int> &oCellMaterials)
{
  std::vector<std::vector<int> > baseMaterialStructures;
  bool resOk = readMaterialCombinations(iMaterialFile, baseMaterialStructures);
  if (!resOk)
    return false;

  int baseMatStructureSize =  (int)pow(baseMaterialStructures[0].size(), 1.f/iDim);

  int N[3] = {baseMatStructureSize, baseMatStructureSize, baseMatStructureSize};
  std::vector<int> cellMaterials;
  if (iDim==3)
  {
    getMaterialAssignment(2, 2, 2, iMaterialCombination, N[0], N[1], N[2], baseMaterialStructures, iBlockRep, iBlockRep, iBlockRep, iNbSubdiv, oCellMaterials);
  }
  else
  {
    getMaterialAssignment(2, 2, iMaterialCombination, N[0], N[1], baseMaterialStructures, iBlockRep, iBlockRep, iNbSubdiv, oCellMaterials);
  }
  return true;
}

void cfgMaterialUtilities::getMaterialAssignment(int nx, int ny, int nz, const std::vector<int> &iMaterialCombIndices, int nX, int nY, int nZ, const std::vector<std::vector<int> > &iBaseCellMaterials, 
                                                 int repX, int repY, int repZ, int iNbSubdivisions, 
                                                 std::vector<int> &oMaterials)
{
  oMaterials.clear();
  oMaterials.resize(nx*ny*nz*nX*nY*nZ*repX*repY*repZ*iNbSubdivisions*iNbSubdivisions*iNbSubdivisions, -1);

  int l, m, n;
  for (l=0; l<nx; l++)
  {
    for (m=0; m<ny; m++)
    {
      for (n=0; n<nz; n++)
      {
        int matCombinationIndex = iMaterialCombIndices[l*ny*nz + m*ny + n];
        const std::vector<int> &matCombination = iBaseCellMaterials[matCombinationIndex];

        int ii,jj,kk;
        for (ii=0; ii<nX; ii++)
        {
          for (jj=0; jj<nY; jj++)
          {
            for (kk=0; kk<nZ; kk++)
            {       
              int indMat = matCombination[ii*nY*nZ + jj*nZ + kk];

              int i,j,k;
              for (i=0; i<repX; i++)
              {
                for (j=0; j<repY; j++)
                {
                  for (k=0; k<repZ; k++)
                  {
                    int ix, iy, iz;
                    for (ix=0; ix<iNbSubdivisions; ix++)
                    {
                      for (iy=0; iy<iNbSubdivisions; iy++)
                      {
                        for (iz=0; iz<iNbSubdivisions; iz++)
                        {
                          int indx = i*(nX*nx*iNbSubdivisions) + nX*l*iNbSubdivisions + ii*iNbSubdivisions + ix;
                          int indy = j*(nY*ny*iNbSubdivisions) + nY*m*iNbSubdivisions + jj*iNbSubdivisions + iy;
                          int indz = k*(nZ*nz*iNbSubdivisions) + nZ*n*iNbSubdivisions + kk*iNbSubdivisions + iz;

                          int elementIndex = indx*nY*nZ*repY*repZ*ny*nz*iNbSubdivisions*iNbSubdivisions + indy * nZ*repZ*nz*iNbSubdivisions + indz;
                          oMaterials[elementIndex] = indMat;
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}

void cfgMaterialUtilities::repMaterialAssignment(int nx, int ny, const std::vector<int> &iMaterialCombIndices, int repX, int repY, int iNbSubdivisions, std::vector<int> &oMaterials)
{
  oMaterials.clear();
  oMaterials.resize(nx*ny*repX*repY*iNbSubdivisions*iNbSubdivisions, -1);

  int m, n;
  for (m=0; m<nx; m++)
  {
    for (n=0; n<ny; n++)
    {
      int matCombinationIndex = iMaterialCombIndices[getGridToVectorIndex(m, n, nx, ny)];

      int i,j;
      for (i=0; i<repX; i++)
      {
        for (j=0; j<repY; j++)
        {
          int k,l;
          for (k=0; k<iNbSubdivisions; k++)
          {
            for (l=0; l<iNbSubdivisions; l++)
            {
              int indx = i*iNbSubdivisions*nx + m*iNbSubdivisions + k;
              int indy = j*iNbSubdivisions*ny + n*iNbSubdivisions + l;
              int elementIndex = getGridToVectorIndex(indx, indy, repX*nx*iNbSubdivisions, repY*ny*iNbSubdivisions);
              oMaterials[elementIndex] = matCombinationIndex;
            }
          }
        }
      }
    }
  }
}

void cfgMaterialUtilities::repMaterialAssignment(int nx, int ny, int nz, const std::vector<int> &iMaterialCombIndices, int repX, int repY, int repZ, int iNbSubdivisions, std::vector<int> &oMaterials)
{
  oMaterials.clear();
  oMaterials.resize(nx*ny*nz*repX*repY*repZ*iNbSubdivisions*iNbSubdivisions*iNbSubdivisions, -1);

  int m, n, p;
  for (m=0; m<nx; m++)
  {
    for (n=0; n<ny; n++)
    {
      for (p=0; p<nz; p++)
      {
        int matCombinationIndex = iMaterialCombIndices[getGridToVectorIndex(m, n, p, nx, ny, nz)];
        int i,j,k;
        for (i=0; i<repX; i++)
        {
          for (j=0; j<repY; j++)
          {
            for (k=0; k<repZ; k++)
            {
              int kk,ll,mm;
              for (kk=0; kk<iNbSubdivisions; kk++)
              {
                for (ll=0; ll<iNbSubdivisions; ll++)
                {
                  for (mm=0; mm<iNbSubdivisions; mm++)
                  {
                    int indx = i*iNbSubdivisions*nx + m*iNbSubdivisions + kk;
                    int indy = j*iNbSubdivisions*ny + n*iNbSubdivisions + ll;
                    int indz = k*iNbSubdivisions*nz + p*iNbSubdivisions + mm;
                    int elementIndex = getGridToVectorIndex(indx, indy, indz, repX*nx*iNbSubdivisions, repY*ny*iNbSubdivisions, repZ*nz*iNbSubdivisions);
                    oMaterials[elementIndex] = matCombinationIndex;
                  }
                }
              }
            }
          }
        }
      }
    }
  }
} 

void cfgMaterialUtilities::getMaterialAssignment(int nx, int ny, const std::vector<int> &iMaterialCombIndices, int nX, int nY, const std::vector<std::vector<int> > &iBaseCellMaterials, int repX, int repY, 
                                                 int iNbSubdivisions, std::vector<int> &oMaterials)
{
  oMaterials.clear();
  oMaterials.resize(nx*ny*nX*nY*repX*repY*iNbSubdivisions*iNbSubdivisions, -1);

  int l, m;
  for (l=0; l<nx; l++)
  {
    for (m=0; m<ny; m++)
    {
      int matCombinationIndex = iMaterialCombIndices[l*ny + m];
      const std::vector<int> &matCombination = iBaseCellMaterials[matCombinationIndex];

      int ii,jj;
      for (ii=0; ii<nX; ii++)
      {
        for (jj=0; jj<nY; jj++)
        {      
          int indMat = matCombination[ii*nY + jj];

          int i,j;
          for (i=0; i<repX; i++)
          {
            for (j=0; j<repY; j++)
            {
              int ix, iy;
              for (ix=0; ix<iNbSubdivisions; ix++)
              {
                for (iy=0; iy<iNbSubdivisions; iy++)
                {
                  int indx = i*(nX*nx*iNbSubdivisions) + nX*l*iNbSubdivisions + ii*iNbSubdivisions + ix;
                  int indy = j*(nY*ny*iNbSubdivisions) + nY*m*iNbSubdivisions + jj*iNbSubdivisions + iy;

                  int elementIndex = indx*nY*repY*ny*iNbSubdivisions + indy;
                  oMaterials[elementIndex] = indMat;
                }
              }
            }
          }
        }
      }
    }
  }
}

int cfgMaterialUtilities::getGridToVectorIndex(int i, int j, int nx, int ny)
{
  int index =  i*ny+ j;
  return index;
}

int cfgMaterialUtilities::getGridToVectorIndex(int i, int j, int k, int nx, int ny, int nz)
{
  int index =  i*ny*nz + j*nz +k;
  return index;
}

void cfgMaterialUtilities::getVectorIndexToGrid(int iIndex, int nx, int ny, int &oIndex_i, int &oIndex_j)
{
  oIndex_i = iIndex/ ny;
  oIndex_j = iIndex % ny;
}

void cfgMaterialUtilities::getSideElements(int iSide, int nx, int ny, std::vector<int> &oElementIndices)
{
  getSideVertices(iSide, nx-1, ny-1, oElementIndices);
}

//0: Left, 1: Right, 2: Bottom, 3:Top
void cfgMaterialUtilities::getSideVertices(int iSide, int nx, int ny, std::vector<int> &oVertexIndices)
{
  oVertexIndices.clear();

  int nxvertex = nx+1, nyvertex = ny+1;

   int ii, jj;
   if (iSide==0) //left
   {    
     for(jj=0; jj<nyvertex; jj++)
     {
       int indVertex = getGridToVectorIndex(0, jj, nxvertex, nyvertex);
       oVertexIndices.push_back(indVertex);
     }
   }
   else if (iSide==1) //right
   {
     for(jj=0; jj<nyvertex; jj++)
     {
       int indVertex = getGridToVectorIndex(nxvertex-1, jj, nxvertex, nyvertex);
       oVertexIndices.push_back(indVertex);
     }
   }
   else if (iSide==2) //bottom
   {
     for(ii=0; ii<nxvertex; ii++)
     {
       int indVertex = getGridToVectorIndex(ii, 0, nxvertex, nyvertex);
       oVertexIndices.push_back(indVertex);
     }
   }
   else if (iSide==3) //top
   {
     for(ii=0; ii<nxvertex; ii++)
     {
       int indVertex = getGridToVectorIndex(ii, nyvertex-1, nxvertex, nyvertex);   
       oVertexIndices.push_back(indVertex);
     }
   }
}

//0: Left, 1: Right, 2: Bottom, 3:Top, 4: Back, 5: Front
void cfgMaterialUtilities::getSideVertices(int iSide, int nx, int ny, int nz, std::vector<int> &oVertexIndices)
{
  oVertexIndices.clear();

  int nxvertex = nx+1, nyvertex = ny+1, nzvertex = nz+1;

   int ii, jj, kk;
   if (iSide==0) //left
   {    
     for(jj=0; jj<nyvertex; jj++)
     {
       for (kk=0; kk<nzvertex; kk++)
       {
        int indVertex = getGridToVectorIndex(0, jj, kk, nxvertex, nyvertex, nzvertex);
        oVertexIndices.push_back(indVertex);
       }
     }
   }
   else if (iSide==1) //right
   {
     for(jj=0; jj<nyvertex; jj++)
     {
       for (kk=0; kk<nzvertex; kk++)
       {
         int indVertex = getGridToVectorIndex(nxvertex-1, jj, kk, nxvertex, nyvertex, nzvertex);
         oVertexIndices.push_back(indVertex);
       }
     }
   }
   else if (iSide==2) //bottom
   {
     for(ii=0; ii<nxvertex; ii++)
     {
       for (kk=0; kk<nzvertex; kk++)
       {
         int indVertex = getGridToVectorIndex(ii, 0, kk, nxvertex, nyvertex, nzvertex);
         oVertexIndices.push_back(indVertex);
       }
     }
   }
   else if (iSide==3) //top
   {
     for(ii=0; ii<nxvertex; ii++)
     {
       for (kk=0; kk<nzvertex; kk++)
       {
         int indVertex = getGridToVectorIndex(ii, nyvertex-1, kk, nxvertex, nyvertex, nzvertex);   
         oVertexIndices.push_back(indVertex);
       }
     }
   }
   else if (iSide==4) // back
   {
     for(ii=0; ii<nxvertex; ii++)
     {
       for (jj=0; jj<nyvertex; jj++)
       {
         int indVertex = getGridToVectorIndex(ii, jj, 0, nxvertex, nyvertex, nzvertex);   
         oVertexIndices.push_back(indVertex);
       }
     }
   }
   else if (iSide==5) // front
   {
     for(ii=0; ii<nxvertex; ii++)
     {
       for (jj=0; jj<nyvertex; jj++)
       {
         int indVertex = getGridToVectorIndex(ii, jj, nzvertex-1, nxvertex, nyvertex, nzvertex);   
         oVertexIndices.push_back(indVertex);
       }
     }
   }
}

void cfgMaterialUtilities::updateMaterialSubstructure(std::vector<int> &ioMaterialAsssignment, int nx, int ny, const std::vector<int> &iSubMaterialStructure, int Nx, int Ny, int iSubMatStructLocation)
{
  int Cx, Cy;
  getVectorIndexToGrid(iSubMatStructLocation, nx, ny, Cx, Cy);
  for (int i=0; i<Nx; i++)
  {
    for (int j=0; j<Ny; j++)
    {
      int x = (Cx+i)%nx;
      int y = (Cy+j)%ny;
      int index = getGridToVectorIndex(x, y, nx, ny);
      int oldMat = ioMaterialAsssignment[index];

      int IndexSubstructure = getGridToVectorIndex(i, j, Nx, Ny);
      int newMat = iSubMaterialStructure[IndexSubstructure];
      ioMaterialAsssignment[index] = newMat;
    }
  }
}

bool cfgMaterialUtilities::isStructureManifold(int nx, int ny, int nz, const std::vector<int> &iMaterials, bool isStructuredMirrored, int nX, int nY, int nZ)
{
  bool isManifold = true;

  int nxComb = nx/nX;
  int nyComb = ny/nY;
  int nzComb = nz/nZ;

  int startIndex = (isStructuredMirrored? 1: 0);

  int icomb, jcomb, kcomb;
  for (icomb=0; icomb<nxComb && isManifold; icomb++)
  {
    for (jcomb=0; jcomb<nyComb && isManifold; jcomb++)
    {
      for (kcomb=0; kcomb<nzComb && isManifold; kcomb++)
      {
        int icell;
        for (icell=0; icell<nX && isManifold; icell++)
        {
          int x = nX*icomb + icell;
          int y = nY*jcomb;
          int z = nZ*kcomb;
          if (x>=startIndex && y>=startIndex && z>=startIndex)
          {
            int indMat111 = getGridToVectorIndex(x, y, z, nx, ny, nz);
            int indMat121 = getGridToVectorIndex(x, (y+ny-1)%ny, z, nx, ny, nz);
            int indMat211 = getGridToVectorIndex((x+nx-1)%nx, y, z, nx, ny, nz);
            int indMat221 = getGridToVectorIndex((x+nx-1)%nx, (y+ny-1)%ny, z, nx, ny, nz);
            int indMat112 = getGridToVectorIndex(x, y, (z+nz-1)%nz, nx, ny, nz);
            int indMat122 = getGridToVectorIndex(x, (y+ny-1)%ny, (z+nz-1)%nz, nx, ny, nz);
            int indMat212 = getGridToVectorIndex((x+nx-1)%nx, y, (z+nz-1)%nz, nx, ny, nz);
            int indMat222 = getGridToVectorIndex((x+nx-1)%nx, (y+ny-1)%ny, (z+nz-1)%nz, nx, ny, nz);

            int mat111 = iMaterials[indMat111];
            int mat121 = iMaterials[indMat121];
            int mat211 = iMaterials[indMat211];
            int mat221 = iMaterials[indMat221];
            int mat112 = iMaterials[indMat112];
            int mat122 = iMaterials[indMat122];
            int mat212 = iMaterials[indMat212];
            int mat222 = iMaterials[indMat222];

            isManifold = mat111!=mat221 || mat121!=mat211 || mat121==mat111;
            isManifold &= (mat112!=mat222 || mat122!=mat212 || mat122==mat112);
            isManifold &= (mat111!=mat122 || mat121!=mat112 || mat111==mat112);
            isManifold &= (mat211!=mat222 || mat212!=mat221 || mat211==mat221);
            isManifold &= (mat111!=mat212 || mat211!=mat112 || mat111==mat112);
            isManifold &= (mat121!=mat222 || mat221!=mat122 || mat121==mat122);
            isManifold &= (mat121!=mat222 || mat221!=mat122 || mat121==mat122);
            isManifold &= (mat111!=mat222 || mat211!=mat122 || mat221!=mat112 || mat121!=mat212 || (mat111==mat211 && mat111==mat112 && mat111==mat212));
          }
        }

        int jcell;
        for (jcell=1; jcell<nY && isManifold; jcell++)
        {
          int x = nX*icomb;
          int y = nY*jcomb + jcell;
          int z = nZ*kcomb;
          if (x>=startIndex && y>=startIndex && z>=startIndex)
          {
            int indMat111 = getGridToVectorIndex(x, y, z, nx, ny, nz);
            int indMat121 = getGridToVectorIndex(x, (y+ny-1)%ny, z, nx, ny, nz);
            int indMat211 = getGridToVectorIndex((x+nx-1)%nx, y, z, nx, ny, nz);
            int indMat221 = getGridToVectorIndex((x+nx-1)%nx, (y+ny-1)%ny, z, nx, ny, nz);
            int indMat112 = getGridToVectorIndex(x, y, (z+nz-1)%nz, nx, ny, nz);
            int indMat122 = getGridToVectorIndex(x, (y+ny-1)%ny, (z+nz-1)%nz, nx, ny, nz);
            int indMat212 = getGridToVectorIndex((x+nx-1)%nx, y, (z+nz-1)%nz, nx, ny, nz);
            int indMat222 = getGridToVectorIndex((x+nx-1)%nx, (y+ny-1)%ny, (z+nz-1)%nz, nx, ny, nz);

            int mat111 = iMaterials[indMat111];
            int mat121 = iMaterials[indMat121];
            int mat211 = iMaterials[indMat211];
            int mat221 = iMaterials[indMat221];
            int mat112 = iMaterials[indMat112];
            int mat122 = iMaterials[indMat122];
            int mat212 = iMaterials[indMat212];
            int mat222 = iMaterials[indMat222];

            isManifold = mat111!=mat221 || mat121!=mat211 || mat121==mat111;
            isManifold &= (mat112!=mat222 || mat122!=mat212 || mat122==mat112);
            isManifold &= (mat111!=mat122 || mat121!=mat112 || mat111==mat112);
            isManifold &= (mat211!=mat222 || mat212!=mat221 || mat211==mat221);
            isManifold &= (mat111!=mat212 || mat211!=mat112 || mat111==mat112);
            isManifold &= (mat121!=mat222 || mat221!=mat122 || mat121==mat122);
            isManifold &= (mat121!=mat222 || mat221!=mat122 || mat121==mat122);
            isManifold &= (mat111!=mat222 || mat211!=mat122 || mat221!=mat112 || mat121!=mat212 || (mat111==mat211 && mat111==mat112 && mat111==mat212));
          }
        }

        int kcell;
        for (kcell=1; kcell<nZ && isManifold; kcell++)
        {
          int x = nX*icomb;
          int y = nY*jcomb;
          int z = nZ*kcomb + kcell;
          if (x>=startIndex && y>=startIndex && z>=startIndex)
          {
            int indMat111 = getGridToVectorIndex(x, y, z, nx, ny, nz);
            int indMat121 = getGridToVectorIndex(x, (y+ny-1)%ny, z, nx, ny, nz);
            int indMat211 = getGridToVectorIndex((x+nx-1)%nx, y, z, nx, ny, nz);
            int indMat221 = getGridToVectorIndex((x+nx-1)%nx, (y+ny-1)%ny, z, nx, ny, nz);
            int indMat112 = getGridToVectorIndex(x, y, (z+nz-1)%nz, nx, ny, nz);
            int indMat122 = getGridToVectorIndex(x, (y+ny-1)%ny, (z+nz-1)%nz, nx, ny, nz);
            int indMat212 = getGridToVectorIndex((x+nx-1)%nx, y, (z+nz-1)%nz, nx, ny, nz);
            int indMat222 = getGridToVectorIndex((x+nx-1)%nx, (y+ny-1)%ny, (z+nz-1)%nz, nx, ny, nz);

            int mat111 = iMaterials[indMat111];
            int mat121 = iMaterials[indMat121];
            int mat211 = iMaterials[indMat211];
            int mat221 = iMaterials[indMat221];
            int mat112 = iMaterials[indMat112];
            int mat122 = iMaterials[indMat122];
            int mat212 = iMaterials[indMat212];
            int mat222 = iMaterials[indMat222];

            isManifold = mat111!=mat221 || mat121!=mat211 || mat121==mat111;
            isManifold &= (mat112!=mat222 || mat122!=mat212 || mat122==mat112);
            isManifold &= (mat111!=mat122 || mat121!=mat112 || mat111==mat112);
            isManifold &= (mat211!=mat222 || mat212!=mat221 || mat211==mat221);
            isManifold &= (mat111!=mat212 || mat211!=mat112 || mat111==mat112);
            isManifold &= (mat121!=mat222 || mat221!=mat122 || mat121==mat122);
            isManifold &= (mat121!=mat222 || mat221!=mat122 || mat121==mat122);
            isManifold &= (mat111!=mat222 || mat211!=mat122 || mat221!=mat112 || mat121!=mat212 || (mat111==mat211 && mat111==mat112 && mat111==mat212));
          }
        }
      }
    }
  }
  return isManifold;
}

bool cfgMaterialUtilities::isStructureManifold(int nx, int ny, const std::vector<int> &iMaterials, int nX, int nY, bool isStructuredMirrored, int iNbCellsToCheck)
{
  bool isManifold = true;

  int nxComb = nx/nX;
  int nyComb = ny/nY;

  int startIndex = (isStructuredMirrored? 1: 0);

  int icomb, jcomb;
  for (icomb=0; icomb<nxComb && isManifold; icomb++)
  {
    for (jcomb=0; jcomb<nyComb && isManifold; jcomb++)
    {
      int indComb = icomb*nxComb + jcomb;

      //if (iNbCellsToCheck>=indComb+2)
      {
        int icell;
        for (icell=0; icell<nX && isManifold; icell++)
        {
          int x = nX*icomb + icell;
          int y = nY*jcomb;
          if (x>=startIndex && y>=startIndex)
          {
            int indMat11 = getGridToVectorIndex(x, y, nx, ny);
            int indMat12 = getGridToVectorIndex(x, (y+ny-1)%ny, nx, ny);
            int indMat21 = getGridToVectorIndex((x+nx-1)%nx, y, nx, ny);
            int indMat22 = getGridToVectorIndex((x+nx-1)%nx, (y+ny-1)%ny, nx, ny);

            int mat11 = iMaterials[indMat11];
            int mat12 = iMaterials[indMat12];
            int mat21 = iMaterials[indMat21];
            int mat22 = iMaterials[indMat22];

            isManifold = mat11!=mat22 || mat12!=mat21 || mat12==mat11;
          }
        }
      }

      //if (iNbCellsToCheck>=indComb+2)
      {
        int jcell;
        for (jcell=1; jcell<nY && isManifold; jcell++)
        {
          int x = nX*icomb;
          int y = nY*jcomb + jcell;
          if (x>=startIndex && y>=startIndex)
          {
            int indMat11 = getGridToVectorIndex(x, y, nx, ny);
            int indMat12 = getGridToVectorIndex(x, (y+ny-1)%ny, nx, ny);
            int indMat21 = getGridToVectorIndex((x+nx-1)%nx, y, nx, ny);
            int indMat22 = getGridToVectorIndex((x+nx-1)%nx, (y+ny-1)%ny, nx, ny);

            int mat11 = iMaterials[indMat11];
            int mat12 = iMaterials[indMat12];
            int mat21 = iMaterials[indMat21];
            int mat22 = iMaterials[indMat22];

            isManifold = mat11!=mat22 || mat12!=mat21 || mat12==mat11;
          }
        }
      }
    }
  }
  return isManifold;
}

bool cfgMaterialUtilities::isStructureManifold(int nx, int ny, const std::vector<int> &iMaterials, bool isStructuredMirrored)
{
  bool isManifold = true;
  int startIndex = (isStructuredMirrored? 1: 0);

  for (int x=startIndex; x<nx && isManifold; x++)
  {
    for (int y=startIndex; y<ny && isManifold; y++)
    {
      int indMat11 = getGridToVectorIndex(x, y, nx, ny);
      int indMat12 = getGridToVectorIndex(x, (y+ny-1)%ny, nx, ny);
      int indMat21 = getGridToVectorIndex((x+nx-1)%nx, y, nx, ny);
      int indMat22 = getGridToVectorIndex((x+nx-1)%nx, (y+ny-1)%ny, nx, ny);

      int mat11 = iMaterials[indMat11];
      int mat12 = iMaterials[indMat12];
      int mat21 = iMaterials[indMat21];
      int mat22 = iMaterials[indMat22];

      isManifold = mat11!=mat22 || mat12!=mat21 || mat12==mat11;
    }
  } 
  return isManifold;
}

void cfgMaterialUtilities::getLayer(int Nx, int Ny, const std::vector<int> &iStructureElements, int iIndex, int iDim, std::vector<int> &oNewLayerElements)
{
  oNewLayerElements.clear();

  int min[2] = {0, 0};
  int max[2] = {Nx, Ny};
  min[iDim]=iIndex;
  max[iDim]=iIndex+1;

  for (int i=min[0]; i<max[0]; i++)
  {
    for (int j=min[1]; j<max[1]; j++)
    {
      int indMat = getGridToVectorIndex(i, j, Nx, Ny);
      oNewLayerElements.push_back(iStructureElements[indMat]);
    }
  }
}

void cfgMaterialUtilities::getLayer(int Nx, int Ny, int Nz, const std::vector<int> &iStructureElements, int iIndex, int iDim, std::vector<int> &oNewLayerElements)
{
  oNewLayerElements.clear();

  int min[3] = {0, 0};
  int max[3] = {Nx, Ny, Nz};
  min[iDim]=iIndex;
  max[iDim]=iIndex+1;

  for (int i=min[0]; i<max[0]; i++)
  {
    for (int j=min[1]; j<max[1]; j++)
    {
      for (int k=min[2]; k<max[2]; k++)
      {
        int indMat = getGridToVectorIndex(i, j, k, Nx, Ny ,Nz);
        oNewLayerElements.push_back(iStructureElements[indMat]);
      }
    }
  }
}

void cfgMaterialUtilities::insertLayer(int Nx, int Ny, const std::vector<int> &iStructureElements, std::vector<int> &iNewLayerElements, int iIndex, int iDim, std::vector<int> &oNewStructureElements)
{
  int n[2] = {Nx, Ny};
  n[iDim]++;

  oNewStructureElements.resize(n[0]*n[1]);

  int indInsertion[2] = {Nx, Ny};
  indInsertion[iDim] = iIndex;

  for (int i=0; i<n[0]; i++)
  {
    for (int j=0; j<n[1]; j++)
    {
      int mat;
      if (iDim==0 && i==iIndex)
      {
        mat = iNewLayerElements[j];
      }
      else if (iDim==1 && j==iIndex)
      {
        mat = iNewLayerElements[i];
      }
      else
      {
        int X = (i<indInsertion[0]? i: i-1);
        int Y = (j<indInsertion[1]? j: j-1);
        int indMat = getGridToVectorIndex(X, Y, Nx, Ny);
        mat = iStructureElements[indMat];
      }
      int indMatNew = getGridToVectorIndex(i, j, n[0], n[1]);
      oNewStructureElements[indMatNew] = mat;
    }
  }
}

void cfgMaterialUtilities::mirrorStructure(int Nx, int Ny, const std::vector<int> &iStructureElements, std::vector<int> &oNewStructureElements)
{
  int n[2] = {2*Nx, 2*Ny};
  oNewStructureElements.resize(n[0]*n[1]);

  for (int i=0; i<Nx; i++)
  {
    for (int j=0; j<Ny; j++)
    {
      int indMat = getGridToVectorIndex(i, j, Nx, Ny);
      int mat = iStructureElements[indMat];

      int indMat1 = getGridToVectorIndex(i, j, n[0], n[1]);
      oNewStructureElements[indMat1] = mat;

      int indMat2 = getGridToVectorIndex(i, Ny-1-j + Ny, n[0], n[1]);
      oNewStructureElements[indMat2] = mat;

      int indMat3 = getGridToVectorIndex(Nx-1-i + Nx, j, n[0], n[1]);
      oNewStructureElements[indMat3] = mat;

      int indMat4 = getGridToVectorIndex(Nx-1-i + Nx, Ny-1-j + Ny, n[0], n[1]);
      oNewStructureElements[indMat4] = mat;
    }
  }
}

void cfgMaterialUtilities::mirrorStructure(int Nx, int Ny, int Nz ,const std::vector<int> &iStructureElements, std::vector<int> &oNewStructureElements)
{
  int n[3] = {2*Nx, 2*Ny, 2*Nz};
  oNewStructureElements.resize(n[0]*n[1]*n[2]);

  for (int i=0; i<Nx; i++)
  {
    for (int j=0; j<Ny; j++)
    {
      for (int k=0; k<Nz; k++)
      {
        int indMat = getGridToVectorIndex(i, j, k, Nx, Ny, Nz);
        int mat = iStructureElements[indMat];

        int indMat1 = getGridToVectorIndex(i, j, k, n[0], n[1], n[2]);
        oNewStructureElements[indMat1] = mat;

        int indMat2 = getGridToVectorIndex(i, Ny-1-j + Ny, k, n[0], n[1], n[2]);
        oNewStructureElements[indMat2] = mat;

        int indMat3 = getGridToVectorIndex(Nx-1-i + Nx, j, k, n[0], n[1], n[2]);
        oNewStructureElements[indMat3] = mat;

        int indMat4 = getGridToVectorIndex(Nx-1-i + Nx, Ny-1-j + Ny, k, n[0], n[1], n[2]);
        oNewStructureElements[indMat4] = mat;

        int indMat5 = getGridToVectorIndex(i, j, Nz-1-k + Nz, n[0], n[1], n[2]);
        oNewStructureElements[indMat5] = mat;

        int indMat6 = getGridToVectorIndex(i, Ny-1-j + Ny, Nz-1-k + Nz, n[0], n[1], n[2]);
        oNewStructureElements[indMat6] = mat;

        int indMat7 = getGridToVectorIndex(Nx-1-i + Nx, j, Nz-1-k + Nz, n[0], n[1], n[2]);
        oNewStructureElements[indMat7] = mat;

        int indMat8 = getGridToVectorIndex(Nx-1-i + Nx, Ny-1-j + Ny, Nz-1-k + Nz, n[0], n[1], n[2]);
        oNewStructureElements[indMat8] = mat;
      }
    }
  }
}

int cfgMaterialUtilities::getDiagonalStructureNumberOfElements(int N, int iDim)
{
  int nelem=0;
  if (iDim==2)
  {
    nelem=(N*(N+1))/2;
  }
  else
  {
    for (int i=0; i<N; i++)
    {
      for (int j=0; j<=i; j++)
      {
        for (int k=0; k<=j; k++)
        {
          nelem++;
        }
      }
    }
  }
  return nelem;
}

void cfgMaterialUtilities::mirrorStructureAlongDiagonal(int Nx, int Ny, const std::vector<int> &iStructureElements, std::vector<int> &oNewStructureElements)
{
  int n[2] = {Nx, Ny};
  oNewStructureElements.resize(n[0]*n[1]);

  assert(n[0]==n[1]);
  assert(iStructureElements.size()==(n[0]*n[1]-n[0])/2+n[0]);

  int index=0;
  for (int i=0; i<Nx; i++)
  {
    for (int j=0; j<=i; j++)
    {
      int mat = iStructureElements[index++];

      int indMat1 = getGridToVectorIndex(i, j, n[0], n[1]);
      oNewStructureElements[indMat1] = mat;
      if (i!=j)
      {
        int indMat2 = getGridToVectorIndex(j, i, n[0], n[1]);
        oNewStructureElements[indMat2] = mat;
      }
    }
  }
}

void cfgMaterialUtilities::mirrorStructureAlongDiagonal(int Nx, int Ny, int Nz, const std::vector<int> &iStructureElements, std::vector<int> &oNewStructureElements)
{
  int n[3] = {Nx, Ny, Nz};
  oNewStructureElements.resize(n[0]*n[1]*n[2]);

  assert(n[0]==n[1] && n[0]==n[2]);
  //assert(iStructureElements.size()==(n[0]*n[1]-n[0])/2+n[0]);

  int index=0;
  for (int i=0; i<Nx; i++)
  {
    for (int j=0; j<=i; j++)
    {
      for (int k=0; k<=j; k++)
      {
        int mat = iStructureElements[index++];

        int indMat1 = getGridToVectorIndex(i, j, k, n[0], n[1], n[2]);
        oNewStructureElements[indMat1] = mat;
        if (k!=j)
        {
          int indMat2 = getGridToVectorIndex(i, k, j, n[0], n[1], n[2]);
          oNewStructureElements[indMat2] = mat;
        }
        if (i!=j)
        {
          int indMat3 = getGridToVectorIndex(j, i, k, n[0], n[1], n[2]);
          oNewStructureElements[indMat3] = mat;
          if (k!=i)
          {
            int indMat4 = getGridToVectorIndex(j, k, i, n[0], n[1], n[2]);
            oNewStructureElements[indMat4] = mat;
          }
        }
        if (k!=i)
        {
          int indMat5 = getGridToVectorIndex(k, j, i, n[0], n[1], n[2]);
          oNewStructureElements[indMat5] = mat;
          if (j!=i)
          {
            int indMat6 = getGridToVectorIndex(k, i, j, n[0], n[1], n[2]);
            oNewStructureElements[indMat6] = mat;
          }
        }   
      }
    }
  }
}


void cfgMaterialUtilities::getQuarter(int Nx, int Ny, const std::vector<int> &iStructureElements, std::vector<int> &oNewStructureElements)
{
  int n[3] = {Nx/2, Ny/2};
  oNewStructureElements.clear();
  oNewStructureElements.reserve(n[0]*n[1]);

  for (int i=0; i<n[0]; i++)
  {
    for (int j=0; j<n[1]; j++)
    {
      int indMat = getGridToVectorIndex(i, j, Nx, Ny);
      int mat = iStructureElements[indMat];
      oNewStructureElements.push_back(mat);
    }
  }
}

void cfgMaterialUtilities::getQuarter(int Nx, int Ny, int Nz, const std::vector<int> &iStructureElements, std::vector<int> &oNewStructureElements)
{
  int n[3] = {Nx/2, Ny/2, Nz/2};
  oNewStructureElements.clear();
  oNewStructureElements.reserve(n[0]*n[1]*n[2]);

  for (int i=0; i<n[0]; i++)
  {
    for (int j=0; j<n[1]; j++)
    {
      for (int k=0; k<n[2]; k++)
      {
        int indMat = getGridToVectorIndex(i, j, k, Nx, Ny, Nz);
        int mat = iStructureElements[indMat];
        oNewStructureElements.push_back(mat);
      }
    }
  }
}

void cfgMaterialUtilities::getTriangularStructure(int Nx, int Ny, const std::vector<int> &iStructureElements, std::vector<int> &oNewStructureElements)
{
  oNewStructureElements.clear();

  for (int i=0; i<Nx/2; i++)
  {
    for (int j=0; j<=i; j++)
    {
      int indMat = getGridToVectorIndex(i, j, Nx, Ny);
      int mat = iStructureElements[indMat];
      oNewStructureElements.push_back(mat);
    }
  }
}

void cfgMaterialUtilities::getTetrahedralStructure(int Nx, int Ny, int Nz, const std::vector<int> &iStructureElements, std::vector<int> &oNewStructureElements)
{
  oNewStructureElements.clear();

  for (int i=0; i<Nx/2; i++)
  {
    for (int j=0; j<=i; j++)
    {
      for (int k=0; k<=j; k++)
      {
        int indMat = getGridToVectorIndex(i, j, k, Nx, Ny, Nz);
        int mat = iStructureElements[indMat];
        oNewStructureElements.push_back(mat);
      }
    }
  }
}

cfgScalar cfgMaterialUtilities::computeStrain(const std::vector<cfgScalar> &ix, int nx, int ny, int iAxis)
{
  Vector2S Barycenters[2];

  int iside;
  for (iside=0; iside<2; iside++)
  {
    Barycenters[iside] = Vector2S::Zero();

    std::vector<int> vertexIndices;
    getSideVertices(2*iAxis+iside, nx, ny, vertexIndices);

    int ivertex=0, nvertex=(int)vertexIndices.size();
    for (ivertex=0; ivertex<nvertex; ivertex++)
    {
      int indVertex = vertexIndices[ivertex];
      Barycenters[iside] += getVector2S(indVertex, ix);
    }
    Barycenters[iside] /= (cfgScalar)nvertex;
  }
  float strain = (Barycenters[0]-Barycenters[1]).norm()-1;
  return strain;
}

cfgScalar cfgMaterialUtilities::computeStrain(const std::vector<cfgScalar> &ix, int nx, int ny, int nz, int iAxis)
{
  Vector3S Barycenters[2];

  int iside;
  for (iside=0; iside<2; iside++)
  {
    Barycenters[iside] = Vector3S::Zero();

    std::vector<int> vertexIndices;
    getSideVertices(2*iAxis+iside, nx, ny, nz, vertexIndices);

    int ivertex=0, nvertex=(int)vertexIndices.size();
    for (ivertex=0; ivertex<nvertex; ivertex++)
    {
      int indVertex = vertexIndices[ivertex];
      Barycenters[iside] += getVector3S(indVertex, ix);
    }
    Barycenters[iside] /= (cfgScalar)nvertex;
  }
  float strain = (Barycenters[0]-Barycenters[1]).norm()-1;
  return strain;
}

void cfgMaterialUtilities::getClosestPoints(const std::vector<float> &iPoints, int iDim, std::vector<int> &iRefPointIndices, float iRange, std::vector<int> &oPoints)
{
  oPoints.clear();
  std::vector<float> pointCloud = getSubVector(iPoints, iDim, iRefPointIndices);
  std::vector<float> distances;
  computeDistancesToPointCloud(iPoints, pointCloud, distances);
  int ipoint=0, npoint=(int)distances.size();
  for (ipoint=0; ipoint<npoint; ipoint++)
  {
    if (distances[ipoint] < iRange)
    {
      oPoints.push_back(ipoint);
    }
  }
}

void cfgMaterialUtilities::getFurthestPointsGreedy(int iOutputNbPoints, const std::vector<cfgScalar> &iPoints, int iDim, std::vector<int> &oPointIndices)
{
  oPointIndices.clear();

  std::set<int> points;
  int ipoint, npoint=(int)iPoints.size()/iDim;
  for (ipoint=1; ipoint<npoint; ipoint++)
  {
    points.insert(ipoint);
  }
  int pointIndex = 0;
  oPointIndices.push_back(pointIndex);
 
  for (ipoint=1; ipoint<iOutputNbPoints && points.size()>0; ipoint++)
  {
    cfgScalar maxDist = 0;
    int furthestPoint = *points.begin();
    std::set<int>::iterator it, it_end=points.end();
    for (it=points.begin(); it!=it_end; it++)
    {
      int indVertex = *it;
      Vector3f p = getVector3f(indVertex, iPoints);

      cfgScalar currentMinDist = FLT_MAX;
      int jpoint, njpoint=(int)oPointIndices.size();
      for (jpoint=0; jpoint<njpoint; jpoint++)
      {
        int indPoint2 = oPointIndices[jpoint];
        Vector3f q = getVector3f(indPoint2, iPoints);
        cfgScalar dist = (p-q).absSquared();
        if (dist < currentMinDist)
        {
           currentMinDist = dist;
        }
      }
      if (currentMinDist > maxDist)
      {
        maxDist = currentMinDist;
        furthestPoint = indVertex;
      }
    }
    oPointIndices.push_back(furthestPoint);
    points.erase(furthestPoint);
  }
}

void cfgMaterialUtilities::getBoundingBox(const std::vector<cfgScalar> &iPoints, int iDim, std::vector<cfgScalar> oBox[2])
{
  assert(iPoints.size()%iDim==0);

  oBox[0].resize(iDim, FLT_MAX);
  oBox[1].resize(iDim, -FLT_MAX);

  int ipoint, npoint=(int)iPoints.size()/iDim;
  for (ipoint=0; ipoint<npoint; ipoint++)
  {
    for (int icoord=0; icoord<iDim; icoord++)
    {
      cfgScalar val = iPoints[iDim*ipoint+icoord];
      if (val < oBox[0][icoord])
      {
        oBox[0][icoord] = val;
      }
      if (val > oBox[1][icoord])
      {
        oBox[1][icoord] = val;
      }
    }
  }
}

bool cfgMaterialUtilities::rescaleData(std::vector<cfgScalar> &ioPoints, int iDim,  const std::vector<cfgScalar> &iTargetBoxLengths)
{
  assert(iTargetBoxLengths.size()==iDim);

  cfgScalar eps = (cfgScalar)1.e-6;
  std::vector<cfgScalar> box[2];
  getBoundingBox(ioPoints, iDim, box);
  std::vector<cfgScalar> scaleFactors;
  int icoord;
  bool degenerated = true;
  for (icoord=0; icoord<iDim; icoord++)
  {
    double length = box[1][icoord]-box[0][icoord];
    double target_length = iTargetBoxLengths[icoord];
    double scale = (length>eps? target_length/length: 0);
    scaleFactors.push_back(scale);
    degenerated &= scale==0;
  }
  int ipoint, npoint=(int)ioPoints.size()/iDim;
  for (ipoint=0; ipoint<npoint; ipoint++)
  {
    for (icoord=0; icoord<iDim; icoord++)
    {
      cfgScalar val = ioPoints[iDim*ipoint+icoord];
      val = scaleFactors[icoord]*(val-box[0][icoord]);
      ioPoints[iDim*ipoint+icoord] = val;
    }
  }
  return !degenerated;
}

void cfgMaterialUtilities::getKMeans(int iNbIterations, int iNbClusters, const std::vector<cfgScalar> &iPoints, int iDim, std::vector<std::vector<int> > &oClusters, std::vector<int> *oCenters)
{
  oClusters.clear();

  std::vector<int> centers;
  int nbClusters = std::min(iNbClusters, (int)iPoints.size()/iDim);

  getFurthestPointsGreedy(nbClusters, iPoints, iDim, centers);
  std::vector<double> centerPositions = convertVec<cfgScalar,double>(getSubVector(iPoints, iDim, centers));

  std::vector<std::vector<int> > clusters;
  int iter;
  for (iter=0; iter<iNbIterations; iter++)
  {
    clusters.clear();
    clusters.resize(nbClusters);

    DistanceTool distanceTool(centerPositions);

    int ipoint, npoint=(int)iPoints.size()/iDim;
    for (ipoint=0; ipoint<npoint; ipoint++)
    {
      double p[3] = {iPoints[3*ipoint], iPoints[3*ipoint+1], iPoints[3*ipoint+2]};
      int closestCenter = distanceTool.getClosestPointIndex(p);
      clusters[closestCenter].push_back(ipoint);
    }
    centerPositions.clear();
    int icluster;
    for (icluster=0; icluster<nbClusters; icluster++)
    {
      std::vector<double> clusterPoints = convertVec<cfgScalar,double>(getSubVector(iPoints, iDim, clusters[icluster]));
      double center[3];
      cfgUtil::getBarycenter<double, 3>(clusterPoints, center);
      int icoord;
      for (icoord=0; icoord<3; icoord++)
      {
        centerPositions.push_back(center[icoord]);
      }
    }
  }
  int icluster;
  for (icluster=0; icluster<nbClusters; icluster++)
  {
    if (clusters[icluster].size()>0)
    {
      oClusters.push_back(clusters[icluster]);
    }
  }

  if (oCenters)
  {
    oCenters->clear();
    std::vector<double> points = convertVec<cfgScalar,double>(iPoints);
    DistanceTool distanceTool(points);
    int icluster;
    for (icluster=0; icluster<nbClusters; icluster++)
    {
      if (clusters[icluster].size()>0)
      {
        double p[3] = {centerPositions[3*icluster], centerPositions[3*icluster+1], centerPositions[3*icluster+2]};
        int closestPoint = distanceTool.getClosestPointIndex(p);
        oCenters->push_back(closestPoint);
      }
    }
  }
}

void cfgMaterialUtilities::computeConvexHull(const std::vector<float> &iPoints, int iDim, std::vector<int> &oConvexHullVertices)
{
  assert(iPoints.size()%iDim==0);

  oConvexHullVertices.clear();
  coordT * points = cfgUtil::createArray<float,coordT>(iPoints);
  int numpoints = (int)iPoints.size()/iDim;
  boolT newIsMalloc = 0;

  /*qh_init_A (stdin, stdout, stderr, 0, NULL);
  int exitcode = setjmp (qh errexit);
  if (!exitcode)
  {
    qh_init_B(points, numpoints, iDim, newIsMalloc);
    qh_qhull();
    qh_check_output();
    print_summary();
  }*/ 

  std::string rboxCommand2, qhullCommand2;
  Qhull qhull(rboxCommand2.c_str(), iDim, numpoints, points, qhullCommand2.c_str());

  QhullVertexList vertices = qhull.vertexList();
  QhullLinkedList<QhullVertex>::const_iterator v_it, v_end=vertices.end();
  for (v_it=vertices.begin(); v_it!=v_end; ++v_it)
  {
    QhullPoint p = (*v_it).point();
    oConvexHullVertices.push_back(p.id());

    /*const coordT *coord = p.coordinates();
    int icoord;
    for (icoord=0; icoord<iDim; icoord++)
    {
      oConvexHullVertices.push_back(coord[icoord]);
    }*/ 
  }

  /*qhull.outputQhull();
  if(qhull.hasQhullMessage()){
    std::cout << "\nResults of qhull\n" << qhull.qhullMessage();
    qhull.clearQhullMessage();
  }

  QhullFacetList facets= qhull.facetList();
  std::cout << "\nFacets created by Qhull::runQhull()\n" << facets;

  QhullVertexList vertices = qhull.vertexList();
  std::cout << "\nVertices created by Qhull::runQhull()\n" << vertices; */ 
}

void cfgMaterialUtilities::computeDelaundayTriangulation(const std::vector<float> &iPoints, int iDim, std::vector<int> &oFaceVertices, std::vector<int> &oBoundaryVertices, std::vector<int> &oBoundaryFaces, std::vector<float> *oDistancesToBoundary)
{
  assert(iPoints.size()%iDim==0);

  oFaceVertices.clear();
  coordT * points = cfgUtil::createArray<float,coordT>(iPoints);
  int numpoints = (int)iPoints.size()/iDim;
  boolT newIsMalloc = 0;

  std::string rboxCommand2, qhullCommand2="d QJ";
  Qhull qhull(rboxCommand2.c_str(), iDim, numpoints, points, qhullCommand2.c_str());

  std::vector<int> allFaceVertices;
  qhullUtilities::getFacetsVertices(qhull.facetList(), allFaceVertices);
  //oFaceVertices = allFaceVertices;

  Vector3f meanLengths = getMeanEdgeLengthPerAxis(iPoints, allFaceVertices, 4);
  Vector3f medianLengths = getMedianEdgeLengthPerAxis(iPoints, allFaceVertices, 4);

  float meanEdgeLength = getMeanEdgeLength(iPoints, allFaceVertices, 4);
  float minEdgeLength = getMinEdgeLength(iPoints, allFaceVertices, 4);
  //qhullUtilities::getBoundaryVertices(qhull, qhull.facetList(), 1.65, *oBoundaryVertices);

  std::vector<QhullFacet> smallFacets, largeFacets;
  //float maxEdgeSize = 1.8; //3D
  float maxEdgeSize = 2.5; //2D 
  //qhullUtilities::sortFaces(qhull, qhull.facetList(), maxEdgeSize, smallFacets, largeFacets);

  maxEdgeSize = 3; //10
  std::vector<float> scale;
  scale.push_back(1.f/meanLengths[0]);
  scale.push_back(1.f/meanLengths[1]);
  scale.push_back(1.f/meanLengths[2]);
  qhullUtilities::sortFaces(qhull, qhull.facetList(), scale, maxEdgeSize, smallFacets, largeFacets);

  qhullUtilities::getBoundaryVertices(smallFacets, oBoundaryVertices, oBoundaryFaces);

  oFaceVertices.clear();
  qhullUtilities::getFacetsVertices(smallFacets, oFaceVertices);

  /*oFaceVertices.clear();
  qhullUtilities::filterFaces(qhull, qhull.facetList(), 1.5, oFaceVertices);
  *oBoundaryVertices = oFaceVertices;*/ 

  if (oDistancesToBoundary)
  {
    oDistancesToBoundary->clear();

    QhullPoints points = qhull.points();
    int ipoint=0, npoint=(int)points.size();
    for (ipoint=0; ipoint<npoint; ipoint++)
    {
      QhullPoint p = points.at(ipoint);
      float dist = qhullUtilities::distToClosestFacet(p, smallFacets);
      oDistancesToBoundary->push_back(dist);
    }
  }

  if (0)
  {
    std::ofstream stream("Summary.txt");
    qhull.outputQhull();
    if(qhull.hasQhullMessage()){
      stream << "\nResults of qhull\n" << qhull.qhullMessage();
      qhull.clearQhullMessage();
    }

    QhullFacetList facets= qhull.facetList();
    stream << "\nFacets created by Qhull::runQhull()\n" << facets;

    QhullVertexList vertices = qhull.vertexList();
    stream << "\nVertices created by Qhull::runQhull()\n" << vertices;
  }
}

Vector3f cfgMaterialUtilities::getMeanEdgeLengthPerAxis(const std::vector<float> &iX,  const std::vector<int> &iIndexArray, int iDim)
{
  assert(iIndexArray.size()%iDim==0);
  assert(iDim==4);

  cfgScalar eps = (cfgScalar)1.e-3;

  int tetfaces[4][3] = { {0, 1, 2}, {0, 1, 3}, {1, 2, 3}, {2, 3, 0} };

  Vector3f MeanLength(0.f,0.f,0.f);
  int nedge=0;
  int itet=0, ntet=(int)iIndexArray.size()/iDim;
  for (itet=0; itet<ntet; itet++)
  {
    int iface;
    for (iface=0; iface<4; iface++)
    {
      int ivertex=0;
      for (ivertex=0; ivertex<3; ivertex++)
      {
        int indVertex1 = iIndexArray[iDim*itet+tetfaces[iface][ivertex]];
        int indVertex2 = iIndexArray[iDim*itet+tetfaces[iface][(ivertex+1)%3]];

        Vector3f p1 = getVector3f(indVertex1, iX);
        Vector3f p2 = getVector3f(indVertex2, iX);
        int icoord=0;
        bool added = false;
        for (icoord=0; icoord<3; icoord++)
        {
          float length = fabs(p2[icoord]-p1[icoord]);
          if (length > eps)
          {
            MeanLength[icoord] += length;
            added = true;
          }
        }
        if (added)
          nedge++;
      }
    }
  }
  if (nedge>0)
    MeanLength /= (float)nedge;

  return MeanLength;
}

Vector3f cfgMaterialUtilities::getMedianEdgeLengthPerAxis(const std::vector<float> &iX,  const std::vector<int> &iIndexArray, int iDim)
{
  assert(iIndexArray.size()%iDim==0);
  assert(iDim==4);

  cfgScalar eps = (cfgScalar)1.e-3;

  int tetfaces[4][3] = { {0, 1, 2}, {0, 1, 3}, {1, 2, 3}, {2, 3, 0} };

  std::vector<float> MeanLength[3];
  int itet=0, ntet=(int)iIndexArray.size()/iDim;
  for (itet=0; itet<ntet; itet++)
  {
    int iface;
    for (iface=0; iface<4; iface++)
    {
      int ivertex=0;
      for (ivertex=0; ivertex<3; ivertex++)
      {
        int indVertex1 = iIndexArray[iDim*itet+tetfaces[iface][ivertex]];
        int indVertex2 = iIndexArray[iDim*itet+tetfaces[iface][(ivertex+1)%3]];

        Vector3f p1 = getVector3f(indVertex1, iX);
        Vector3f p2 = getVector3f(indVertex2, iX);
        int icoord=0;
        for (icoord=0; icoord<3; icoord++)
        {
          float length = fabs(p2[icoord]-p1[icoord]);
          if (length > eps)
          {
            MeanLength[icoord].push_back(length);
          }
        }
      }
    }
  }
  std::sort(MeanLength[0].begin(), MeanLength[0].end());
  std::sort(MeanLength[1].begin(), MeanLength[1].end());
  std::sort(MeanLength[2].begin(), MeanLength[2].end());

  int indMedian = (int)MeanLength[0].size()/2;
  Vector3f median;
  median[0] = MeanLength[0][indMedian];
  median[1] = MeanLength[1][indMedian];
  median[2] = MeanLength[2][indMedian];

  return median;
}


float cfgMaterialUtilities::getMeanEdgeLength(const std::vector<float> &iX,  const std::vector<int> &iIndexArray, int iDim)
{
  assert(iIndexArray.size()%iDim==0);
  assert(iDim==4);

  int tetfaces[4][3] = { {0, 1, 2}, {0, 1, 3}, {1, 2, 3}, {2, 3, 0} };

  float MeanEdgeLength=0;
  int nedge=0;
  int itet=0, ntet=(int)iIndexArray.size()/iDim;
  for (itet=0; itet<ntet; itet++)
  {
    int iface;
    for (iface=0; iface<4; iface++)
    {
      int ivertex=0;
      for (ivertex=0; ivertex<3; ivertex++)
      {
        int indVertex1 = iIndexArray[iDim*itet+tetfaces[iface][ivertex]];
        int indVertex2 = iIndexArray[iDim*itet+tetfaces[iface][(ivertex+1)%3]];

        Vector3f p1 = getVector3f(indVertex1, iX);
        Vector3f p2 = getVector3f(indVertex2, iX);
        float length = (p2-p1).abs();
        MeanEdgeLength += length;
        nedge++;
      }
    }
  }
  if (nedge>0)
    MeanEdgeLength /= nedge;

  return MeanEdgeLength;
}

float cfgMaterialUtilities::getMinEdgeLength(const std::vector<float> &iX,  const std::vector<int> &iIndexArray, int iDim)
{
  assert(iIndexArray.size()%iDim==0);

  int tetfaces[4][3] = { {0, 1, 2}, {0, 1, 3}, {1, 2, 3}, {2, 3, 0} };

  float SqMinEdgeLength=FLT_MAX;
  int itet=0, ntet=(int)iIndexArray.size()/iDim;
  for (itet=0; itet<ntet; itet++)
  {
    int iface;
    for (iface=0; iface<4; iface++)
    {
      int ivertex=0;
      for (ivertex=0; ivertex<3; ivertex++)
      {
        int indVertex1 = iIndexArray[iDim*itet+tetfaces[iface][ivertex]];
        int indVertex2 = iIndexArray[iDim*itet+tetfaces[iface][(ivertex+1)%3]];

        Vector3f p1 = getVector3f(indVertex1, iX);
        Vector3f p2 = getVector3f(indVertex2, iX);
        float SqLength = (p2-p1).absSquared();
        if (SqLength < SqMinEdgeLength)
        {
          SqMinEdgeLength = SqLength;
        }
      }
    }
  }
  return sqrt(SqMinEdgeLength);
}

void cfgMaterialUtilities::getEdgesFromTriFaceIndexArray(const std::vector<int> &iTriIndexArray, std::vector<int> &oEdgeIndexArray)
{
  std::vector<int> edgeIndices;
  int ifacet=0, nfacet=(int)iTriIndexArray.size()/3;
  for (ifacet=0; ifacet<nfacet; ifacet++)
  {
    int ipoint=0;
    for (ipoint=0; ipoint<3; ipoint++)
    {
      int indPoint1 = iTriIndexArray[3*ifacet+ipoint];
      oEdgeIndexArray.push_back(indPoint1);

      int indPoint2 = iTriIndexArray[3*ifacet+(ipoint+1)%3];
      oEdgeIndexArray.push_back(indPoint2);
    }
  }
}

void cfgMaterialUtilities::getEdgesFromTetFaceIndexArray(const std::vector<int> &iTetIndexArray, std::vector<int> &oEdgeIndexArray)
{
  oEdgeIndexArray.clear();

  int tetfaces[4][3] = { {0, 1, 2}, {0, 1, 3}, {1, 2, 3}, {2, 3, 0} };

  std::vector<int> edgeIndices;
  int ifacet=0, nfacet=(int)iTetIndexArray.size()/4;
  for (ifacet=0; ifacet<nfacet; ifacet++)
  {
    int iface=0;
    for (iface=0; iface<4; iface++)
    {
      int ipoint=0;
      for (ipoint=0; ipoint<4; ipoint++)
      {
        int indVertex1 = tetfaces[iface][ipoint%3];
        int indPoint1 = iTetIndexArray[4*ifacet+indVertex1];
        oEdgeIndexArray.push_back(indPoint1);

        int indVertex2 = tetfaces[iface][(ipoint+1)%3];
        int indPoint2 = iTetIndexArray[4*ifacet+indVertex2];
        oEdgeIndexArray.push_back(indPoint2);
      }
    }
  }
}
void cfgMaterialUtilities::fromLamesParametersToYoungModulusPoissonRatio(cfgScalar iLambda, cfgScalar iMu, cfgScalar &oYoungModulus, cfgScalar &oPoissonRatio)
{
  oYoungModulus = iMu*(3*iLambda+2*iMu)/(iLambda+iMu);
  oPoissonRatio = iLambda/(2*(iLambda+iMu));
}

void cfgMaterialUtilities::fromYoungModulusPoissonRatioToLamesParameters(cfgScalar iYoungModulus, cfgScalar iPoissonRatio, cfgScalar &oLambda, cfgScalar &oMu)
{
  oLambda = iYoungModulus*iPoissonRatio/((1+iPoissonRatio)*(1-2*iPoissonRatio));
  oMu = iYoungModulus/(2*(1+iPoissonRatio));
}

Vector3f cfgMaterialUtilities::getVector3f(int indVertex, const std::vector<float> &iPoints)
{
  assert(iPoints.size()%3==0);
  int indPoint = 3*indVertex;
  return Vector3f(iPoints[indPoint], iPoints[indPoint+1], iPoints[indPoint+2]);
}

Vector2f cfgMaterialUtilities::getVector2f(int indVertex, const std::vector<float> &iPoints)
{
  assert(iPoints.size()%2==0);
  int indPoint = 2*indVertex;
  return Vector2f(iPoints[indPoint], iPoints[indPoint+1]);
}

Vector2S cfgMaterialUtilities::getVector2S(int indVertex, const std::vector<cfgScalar> &iPoints)
{
  assert(iPoints.size()%2==0);
  int indPoint = 2*indVertex;
  return Vector2S(iPoints[indPoint], iPoints[indPoint+1]);
}

Vector3S cfgMaterialUtilities::getVector3S(int indVertex, const std::vector<cfgScalar> &iPoints)
{
  assert(iPoints.size()%3==0);
  int indPoint = 3*indVertex;
  return Vector3S(iPoints[indPoint], iPoints[indPoint+1], iPoints[indPoint+2]);
}

Vector3d cfgMaterialUtilities::getVector3d(int indVertex, const std::vector<double> &iPoints)
{
  assert(iPoints.size()%3==0);
  int indPoint = 3*indVertex;
  return Vector3d(iPoints[indPoint], iPoints[indPoint+1], iPoints[indPoint+2]);
}

Vector2d cfgMaterialUtilities::getVector2d(int indVertex, const std::vector<double> &iPoints)
{
  assert(iPoints.size()%2==0);
  int indPoint = 2*indVertex;
  return Vector2d(iPoints[indPoint], iPoints[indPoint+1]);
}

std::vector<float> cfgMaterialUtilities::toVectorFloat(const std::vector<Vector3f> &iPoints)
 {
   std::vector<float> vec;
   vec.reserve(iPoints.size()*3);
   int ipoint=0, npoint=(int)iPoints.size();
   for (ipoint=0; ipoint<npoint; ipoint++)
   {
     const Vector3f & p = iPoints[ipoint];
     int icoord=0;
     for (icoord=0; icoord<3; icoord++)
     {
       vec.push_back(p[icoord]);
     }
   }
   return vec;
 }

 std::vector<float> cfgMaterialUtilities::toVectorFloat(const std::vector<Vector2f> &iPoints)
 {
   std::vector<float> vec;
   vec.reserve(iPoints.size()*2);
   int ipoint=0, npoint=(int)iPoints.size();
   for (ipoint=0; ipoint<npoint; ipoint++)
   {
     const Vector2f & p = iPoints[ipoint];
     int icoord=0;
     for (icoord=0; icoord<2; icoord++)
     {
       vec.push_back(p[icoord]);
     }
   }
   return vec;
 }

  std::vector<cfgScalar> cfgMaterialUtilities::toVectorScalar(const std::vector<Vector3S> &iPoints)
 {
   std::vector<cfgScalar> vec;
   vec.reserve(iPoints.size()*3);
   int ipoint=0, npoint=(int)iPoints.size();
   for (ipoint=0; ipoint<npoint; ipoint++)
   {
     const Vector3S & p = iPoints[ipoint];
     int icoord=0;
     for (icoord=0; icoord<3; icoord++)
     {
       vec.push_back(p[icoord]);
     }
   }
   return vec;
 }

 std::vector<cfgScalar> cfgMaterialUtilities::toVectorScalar(const std::vector<Vector2S> &iPoints)
 {
   std::vector<cfgScalar> vec;
   vec.reserve(iPoints.size()*2);
   int ipoint=0, npoint=(int)iPoints.size();
   for (ipoint=0; ipoint<npoint; ipoint++)
   {
     const Vector2S & p = iPoints[ipoint];
     int icoord=0;
     for (icoord=0; icoord<2; icoord++)
     {
       vec.push_back(p[icoord]);
     }
   }
   return vec;
 }

 std::vector<Vector2f> cfgMaterialUtilities::toVector2f(const std::vector<float> &iPoints)
 {
   std::vector<Vector2f> vec;
   assert(iPoints.size()%2==0);
   int ipoint=0, npoint=(int)iPoints.size()/2;
   for (ipoint=0; ipoint<npoint; ipoint++)
   {
     vec.push_back(Vector2f(iPoints[2*ipoint], iPoints[2*ipoint+1]));
   }
   return vec;
 }

 std::vector<Vector3f> cfgMaterialUtilities::toVector3f(const std::vector<float> &iPoints)
 {
   std::vector<Vector3f> vec;
   assert(iPoints.size()%3==0);
   int ipoint=0, npoint=(int)iPoints.size()/3;
   for (ipoint=0; ipoint<npoint; ipoint++)
   {
     vec.push_back(Vector3f(iPoints[3*ipoint], iPoints[3*ipoint+1], iPoints[3*ipoint+2]));
   }
   return vec;
 }

 std::vector<Vector2S> cfgMaterialUtilities::toVector2S(const std::vector<cfgScalar> &iPoints)
 {
   std::vector<Vector2S> vec;
   assert(iPoints.size()%2==0);
   int ipoint=0, npoint=(int)iPoints.size()/2;
   for (ipoint=0; ipoint<npoint; ipoint++)
   {
     vec.push_back(Vector2S(iPoints[2*ipoint], iPoints[2*ipoint+1]));
   }
   return vec;
 }

 MatrixXS cfgMaterialUtilities::toMatrixScalar(const MatrixEXd &iMatrix)
 {
   int nrow = iMatrix.rows();
   int ncol = iMatrix.cols();

   MatrixXS mat(nrow, ncol);
   for (int irow=0; irow<nrow; irow++)
   {
     for (int icol=0; icol<ncol; icol++)
     {
       mat(irow, icol) = (cfgScalar)iMatrix(irow, icol);
     }
   }
   return mat;
 }

 MatrixEXd cfgMaterialUtilities::toMatrixDouble(const MatrixXS &iMatrix)
 {
   int nrow = iMatrix.rows();
   int ncol = iMatrix.cols();

   MatrixEXd mat(nrow, ncol);
   for (int irow=0; irow<nrow; irow++)
   {
     for (int icol=0; icol<ncol; icol++)
     {
       mat(irow, icol) = (cfgScalar)iMatrix(irow, icol);
     }
   }
   return mat;
 }

void cfgMaterialUtilities::sampleMesh(const std::vector<float> &iPoints, const std::vector<int> &iTriIndexArray, int iNbSamplesPerFace, std::vector<float> &oPoints)
{
  int ifacet=0, nfacet=(int)iTriIndexArray.size()/3;
  for (ifacet=0; ifacet<nfacet; ifacet++)
  {
    Vector3f p[3];
    float a[3];
    int ipoint=0;
    for (ipoint=0; ipoint<3; ipoint++)
    {
      int indPoint = iTriIndexArray[3*ifacet+ipoint];
      p[ipoint] =  getVector3f(indPoint, iPoints);
    }
    int isample=0;
    for (isample=0; isample<iNbSamplesPerFace; isample++)
    {
      float sumA = 0;
      for (ipoint=0; ipoint<3; ipoint++)
      {
        a[ipoint] = (float)rand()/float(RAND_MAX) + 1.e-6f;
        sumA += a[ipoint];
      }

      Vector3f newPoint(0,0,0);
      for (ipoint=0; ipoint<3; ipoint++)
      {
        newPoint +=  a[ipoint]*p[ipoint];
      }
      newPoint /= sumA;
      int icoord;
      for (icoord=0; icoord<3; icoord++)
      {
        oPoints.push_back(newPoint[icoord]);
      }
    }
  }
}

void cfgMaterialUtilities::computeDistancesToPointCloud(const std::vector<float> &iPoints, const std::vector<float> &iPointCloud, std::vector<float> &oDistances)
{
  oDistances.clear();

  std::vector<double> cloudPositions = convertVec<float,double>(iPointCloud);
  std::vector<double> pointsToTest = convertVec<float,double>(iPoints);

  DistanceTool distanceTool(cloudPositions);
  int ipoint=0, npoint=(int)pointsToTest.size()/3;
  for (ipoint=0; ipoint<npoint; ipoint++)
  {
    double ClosestPoint[3];
    distanceTool.getClosestPoint(&pointsToTest[3*ipoint], ClosestPoint);
    float Dist = (Vector3f(ClosestPoint[0],ClosestPoint[1],ClosestPoint[2])- Vector3f(pointsToTest[3*ipoint], pointsToTest[3*ipoint+1], pointsToTest[3*ipoint+2])).abs();
    oDistances.push_back(Dist);
  }
}

void cfgMaterialUtilities::clusterByBoundary(int nx, int ny, const std::vector<std::vector<int> > &icellMaterials, int iNbMat, std::vector<std::vector<int> > &oClusters)
{
  int nboundary = 2*(nx-1) + 2*(ny-1);
  int ncluster = (int)pow(iNbMat, nboundary);
 
  std::vector<int> boundaryIndices;
  getBoundaryCellIndices(nx, ny, boundaryIndices);

  std::vector<std::vector<int> > clusters(ncluster);
  int icomb, ncomb=(int)icellMaterials.size();
  for (icomb=0; icomb<ncomb; icomb++)
  {
    int signature = getMaterialSignature(nx, ny, icellMaterials[icomb], boundaryIndices);
    clusters[signature].push_back(icomb);
  }

   oClusters.reserve(ncluster);
   int icluster;
   for (icluster=0; icluster<ncluster; icluster++)
   {
     if (clusters[icluster].size()>0)
     {
       oClusters.push_back(clusters[icluster]);
     }
   }
}

void cfgMaterialUtilities::getBoundaryCellIndices(int nx, int ny, std::vector<int> &oCellIndices)
{
  assert(nx>1 && ny>1);

  oCellIndices.clear();
  int iside;
  for (iside=0; iside<2; iside++)
  {
    std::vector<int> elemIndices;
    getSideElements(iside, nx, ny, elemIndices);
    int ielem, nelem=(int)elemIndices.size();
    for (ielem=0; ielem<nelem; ielem++)
    {
      oCellIndices.push_back(elemIndices[ielem]);
    }
  }
  for (iside=2; iside<4; iside++)
  {
    std::vector<int> elemIndices;
    getSideElements(iside, nx, ny, elemIndices);
    int ielem, nelem=(int)elemIndices.size();
    for (ielem=1; ielem<nelem-1; ielem++)
    {
      oCellIndices.push_back(elemIndices[ielem]);
    }
  }
}

int cfgMaterialUtilities::getMaterialSignature(int nx, int ny, const std::vector<int> &icellMaterials, const std::vector<int> &iIndicesToUse)
{
  int signature = 0;
  int icell, ncell=(int)iIndicesToUse.size();
  assert(ncell>0);

  signature = icellMaterials[iIndicesToUse[0]];
  int factor = 2;
  for (icell=1; icell<ncell; icell++)
  {
    int mat = icellMaterials[iIndicesToUse[icell]];
    signature += factor*mat;
    factor *= 2;
  }
  return signature;
}

std::vector<int> cfgMaterialUtilities::genIncrementalSequence(int iMin, int iMax, int iStep)
{
  std::vector<int> values;
  for (int ival=iMin; ival<=iMax; ival+=iStep)
  {
    values.push_back(ival);
  }
  return values;
}






