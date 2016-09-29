#include "MaterialParametersView.h"

#include <vtkRenderView.h>
#include <vtkContextView.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>

#include <vtkContextScene.h>

#include <cfgChartXYZ.h>
#include <vtkFloatArray.h>
#include <vtkTable.h>
#include <cfgPlotPoints3D.h>
#include <vtkChart.h>
#include <cfgPlotLine3D.h>
#include <cfgPlotSurface.h>

#include <vtkPointPicker.h>
#include <vtkPolyData.h>
#include <vtkIdTypeArray.h>
#include <vtkPen.h>

#include <vtkContext2D.h>
#include <vtkContext3D.h>
#include <vtkAxis.h>

#include "exProject.h"

#include "ExplorerUtilities.h"
using namespace ExplorerUtilities;

#include <cfgMaterialUtilities.h>
using namespace cfgMaterialUtilities;

#include <cfgUtilities.h>
using namespace cfgUtil;

#include "QVTKWidgetUtilities.h"
using namespace QVTKWidgetUtilities;

//#include "ScoringFunction.h"
#include "Resampler.h"
#include "DistanceField.h"

MaterialParametersView::MaterialParametersView()
  :QVTKWidget()
{
  m_samplerRadius = 0;
  m_samplerNbPoints = 0;

  m_smallDotSize = 7; //5;
  m_largeDotSize = 15; //10;
}

MaterialParametersView::~MaterialParametersView()
{
  SAFE_DELETE_VEC(m_resamplers);
}

std::vector<cfgScalar> computeScores(const std::vector<cfgScalar> &iPoints, int iDim)
{
  std::vector<cfgScalar> scores;

  int ipoint, npoint=(int)iPoints.size()/iDim;
  for (ipoint=0; ipoint<npoint; ipoint++)
  {
    Vector3S P = getVector3S(ipoint, iPoints);
    cfgScalar score = (P[0]+P[1]+P[2]);
    scores.push_back(score);
  }
  return scores;
}

void MaterialParametersView::rescaleParameters(int iIndex, cfgScalar iScalingFactor, std::vector<cfgScalar> &ioPoints, int iDim)
{
  assert(ioPoints.size()%iDim==0);
  int ipoint=0, npoint=(int)ioPoints.size()/iDim;
  for (ipoint=0; ipoint<npoint; ipoint++)
  {
    ioPoints[iDim*ipoint + iIndex] *= iScalingFactor;
  }
}

cfgScalar MaterialParametersView::getMaxValue(int iIndex, const std::vector<cfgScalar> &iPoints, int iDim)
{
  cfgScalar maxValue = -FLT_MAX;

  assert(iPoints.size()%iDim==0);
  int ipoint=0, npoint=(int)iPoints.size()/iDim;
  for (ipoint=0; ipoint<npoint; ipoint++)
  {
    if (iPoints[iDim*ipoint + iIndex]>maxValue)
    {
      maxValue = iPoints[iDim*ipoint + iIndex];
    }
  }
  return maxValue;
}

void MaterialParametersView::convertToLogValues(int iIndex, const std::vector<cfgScalar> &iPoints, int iDim, std::vector<cfgScalar> &oPoints, cfgScalar iEpsilon)
{
  assert(iPoints.size()%iDim==0);

  oPoints = iPoints;
  int ipoint=0, npoint=(int)iPoints.size()/iDim;
  for (ipoint=0; ipoint<npoint; ipoint++)
  {
    int coord = iDim*ipoint + iIndex;
    cfgScalar value = iPoints[coord];
    oPoints[coord] = log(value+iEpsilon);
  }
}

int MaterialParametersView::getParameterDim()
{
  exProject::MicrostructureType type = m_project->getType();

  int paramdim = 0;
  if (type == exProject::CubicType)
  {
    paramdim = 4;
  }
  else if (type == exProject::OrthotropicType)
  {
    if (m_project->getDim()==2)
    {
      paramdim = 5;
    }
    else
    {
      paramdim = 10;
    }
  }
  return paramdim;
}

void MaterialParametersView::computeDistancesToBoundary(const std::vector<std::vector<cfgScalar> > &iParameters, int iParamDim)
{
  int paramdim = getParameterDim(); 
  int nlevel=(int)iParameters.size();

  m_distancesToBoundary.resize(nlevel);
  m_resamplers.resize(nlevel);

  for (int ilevel=0; ilevel<nlevel; ilevel++)
  {
    if (!m_resamplers[ilevel])
    {
      DistanceField distanceField(iParamDim);
      std::vector<cfgScalar> derivatives;
      m_distancesToBoundary[ilevel] = distanceField.computeDistances(iParameters[ilevel], &derivatives);

       //derivatives = cfgUtil::mult<cfgScalar>(derivatives, 0.2);
      //std::vector<cfgScalar> newPoints = cfgUtil::add<cfgScalar>(iParameters, derivatives);

      Resampler * resampler = new Resampler();
      resampler->setPoints(&iParameters[ilevel], iParamDim);
      resampler->setScores(&(m_distancesToBoundary[ilevel]));
      resampler->init();
      m_resamplers[ilevel] = resampler;
    }
  }
}

std::string MaterialParametersView::parameterTypeToString(const exProject::MaterialParameterType &iType)
{
  std::string str;
  switch (iType)
  {
  case exProject::EXType:
    str = "EX";
    break;
  case exProject::EYType:
    str = "EY";
    break;
  case exProject::EZType:
    str = "EZ";
    break;
  case exProject::NuXYType:
    str = "NuXY";
    break;
  case exProject::NuXZType:
    str = "NuXZ";
    break;
  case exProject::NuYZType:
    str = "NuYZ";
    break;
  case exProject::MuXYType:
    str = "MuXY";
    break;
  case exProject::MuXZType:
    str = "MuXZ";
    break;
    case exProject::MuYZType:
    str = "MuYZ";
    break;
  case exProject::DensityType:
    str = "Density";
    break;
  case exProject::StrengthType:
    str = "Strength";
    break;
  }
  return str;
}

void MaterialParametersView::getValidParameterTypes(exProject::MicrostructureType iMicrostructureType, int iMicrostructureDimension, std::vector<exProject::MaterialParameterType> &oValidParameterTypes, std::vector<int> &oDimToRescale)
{
  oValidParameterTypes.clear();

  if (iMicrostructureType == exProject::CubicType)
  {
    oValidParameterTypes.push_back(exProject::DensityType);
    oValidParameterTypes.push_back(exProject::EXType);
    oValidParameterTypes.push_back(exProject::NuXYType);
    oValidParameterTypes.push_back(exProject::MuXYType);
    oValidParameterTypes.push_back(exProject::StrengthType);
  
    oDimToRescale.push_back(1);
  }
  else if (iMicrostructureType == exProject::OrthotropicType)
  {
    if (iMicrostructureDimension==2)
    {
      oValidParameterTypes.push_back(exProject::DensityType);
      oValidParameterTypes.push_back(exProject::EXType);
      oValidParameterTypes.push_back(exProject::EYType);
      oValidParameterTypes.push_back(exProject::NuXYType);
      oValidParameterTypes.push_back(exProject::MuXYType);
      oValidParameterTypes.push_back(exProject::StrengthType);
 
      oDimToRescale.push_back(1);
      oDimToRescale.push_back(2);
    }
    else
    {
      oValidParameterTypes.push_back(exProject::DensityType);
      oValidParameterTypes.push_back(exProject::EXType);
      oValidParameterTypes.push_back(exProject::EYType);
      oValidParameterTypes.push_back(exProject::EZType);
      oValidParameterTypes.push_back(exProject::NuXYType);
      oValidParameterTypes.push_back(exProject::NuXZType);
      oValidParameterTypes.push_back(exProject::NuYZType);
      oValidParameterTypes.push_back(exProject::MuXYType);
      oValidParameterTypes.push_back(exProject::MuXZType);
      oValidParameterTypes.push_back(exProject::MuYZType);
      oValidParameterTypes.push_back(exProject::StrengthType);

      oDimToRescale.push_back(1);
      oDimToRescale.push_back(2);
      oDimToRescale.push_back(3);
    }
  }
}

void MaterialParametersView::updatePlots()
{
  srand(0);

  Q_ASSERT(m_project);

  const std::vector<std::vector<cfgScalar> > & physicalParametersPerLevelInit = m_project->getPhysicalParameters();
  int nlevel=(int)physicalParametersPerLevelInit.size();
  int paramdim = getParameterDim();

  std::vector<std::vector<cfgScalar> > physicalParametersPerLevel = physicalParametersPerLevelInit;

  // Add additional properties
  // -------------------------
  bool addStrengthParameters = false;
  if (addStrengthParameters)
  {
    for (int ilevel=0; ilevel<nlevel; ilevel++)
    {
      int npoints = (int)physicalParametersPerLevel[ilevel].size()/paramdim;

      std::vector<cfgScalar> newParameters;
      int newParamDim = paramdim;
      const std::vector<std::vector<cfgScalar> > & strengths = m_project->getUltimateStrengths();
      if (strengths[ilevel].size()>0)
      {
        newParamDim ++;
        for (int ipoint=0; ipoint<npoints; ipoint++)
        {
          for (int iparam=0; iparam<paramdim; iparam++)
          {
            newParameters.push_back(physicalParametersPerLevel[ilevel][paramdim*ipoint+iparam]);
          }
          cfgScalar strength = strengths[ilevel][ipoint];
          if (strength > 1)
          {
            strength = 1;
          }
          newParameters.push_back(strength);
        }
        physicalParametersPerLevel[ilevel] = newParameters;
        paramdim = newParamDim;
      }
    }
  }

  exProject::MicrostructureType microstructureType = m_project->getType();
  int microstructureDim = m_project->getDim();
  std::vector<int> dimToRescale;
  std::vector<exProject::MaterialParameterType> validTypes;
  getValidParameterTypes(microstructureType, microstructureDim, validTypes, dimToRescale);

  // Rescale parameters
  // ------------------
  std::vector<cfgScalar> maxValues(paramdim, -FLT_MAX);
  for (int idim=0; idim<paramdim; idim++)
  {
    for (int ilevel=0; ilevel<nlevel; ilevel++)
    {
      cfgScalar maxValue = getMaxValue(idim, physicalParametersPerLevel[ilevel], paramdim);
      if (maxValue > maxValues[idim])
      {
        maxValues[idim] = maxValue;
      }
    }
  }
  int ilevel=0;
  bool useLogScale = true;
  for (ilevel=0; ilevel<nlevel; ilevel++)
  {
    if (useLogScale)
    {
      if (0)
      {
        for (int idim=0; idim<paramdim; idim++)
        {
          cfgScalar maxValue = maxValues[idim];
          rescaleParameters(idim, 1./maxValue, physicalParametersPerLevel[ilevel], paramdim);
          std::vector<cfgScalar> newPoints;
          convertToLogValues(idim, physicalParametersPerLevel[ilevel], paramdim, newPoints, 1.e-6);
          physicalParametersPerLevel[ilevel] = newPoints;
        } 
      }
      else
      {
        cfgScalar scalingFactor_E = 1.f/290.91;
        for (int idim=0; idim<dimToRescale.size(); idim++)
        {
          std::vector<cfgScalar> newPoints;
          convertToLogValues(dimToRescale[idim], physicalParametersPerLevel[ilevel], paramdim, newPoints);
          physicalParametersPerLevel[ilevel] = newPoints;
          float maxE = getMaxValue(dimToRescale[idim], physicalParametersPerLevel[ilevel], paramdim);
          if (maxE > 2)
          {
            rescaleParameters(dimToRescale[idim], scalingFactor_E, physicalParametersPerLevel[ilevel], paramdim);
          }
        }
      }
    }
  }

  // Compute distances to boundary
  // -----------------------------
  bool computeDistances = true;
  if (computeDistances)
  {
    computeDistancesToBoundary(physicalParametersPerLevel, paramdim);
    m_project->setDistancesToBoundary(m_distancesToBoundary);
  }

  std::vector<cfgScalar> lengths(paramdim, 1);
  for (ilevel=0; ilevel<nlevel; ilevel++)
  {
    bool resOk = rescaleData(physicalParametersPerLevel[ilevel], paramdim, lengths);
  }
  m_project->setPhysicalParameters(physicalParametersPerLevel);

  // Plot gamuts
  // -----------
  std::vector<std::string> labels;
  for (int itype=0; itype<(int)validTypes.size(); itype++)
  {
    labels.push_back(parameterTypeToString(validTypes[itype]));
  }
  labels.push_back("Level");

  std::vector<vtkVector3i> matColors(nlevel, vtkVector3i(0, 0, 255));
  if (nlevel>1)
  {
    for (ilevel=0; ilevel<nlevel; ilevel++)
    {
      QColor colHSV = QColor::fromHsv(ilevel*255/(nlevel-1), 255, 255);
      QColor colRGB = colHSV.toRgb();
      matColors[ilevel] = vtkVector3i(colRGB.red(), colRGB.green(), colRGB.blue());
    }
  }

  std::set<exProject::MaterialParameterType> setValidTypes = toStdSet<exProject::MaterialParameterType>(validTypes);

  m_plotsPerLevel.clear();
  m_plotsPerLevel.resize(nlevel);
  m_tables.clear();

  for (ilevel=0; ilevel<nlevel; ilevel++)
  {
    int npoints = (int)physicalParametersPerLevel[ilevel].size()/paramdim;
    std::cout << "npoints = " << npoints << std::endl;

    std::vector<int> levels(npoints, ilevel);
    std::vector<vtkVector3i> colors(npoints, matColors[ilevel]);

    if (computeDistances)
    {
      getColorsForDistanceVisualization(m_distancesToBoundary[ilevel], colors);
    }
    vtkSmartPointer<vtkTable> table =  createTable(physicalParametersPerLevel[ilevel], levels, paramdim, 1, labels);
    m_tables.push_back(table);

    exProject::MaterialParameterType paramType1 = m_project->getParameterToVisualize(0);
    exProject::MaterialParameterType paramType2 = m_project->getParameterToVisualize(1);
    exProject::MaterialParameterType paramType3 = m_project->getParameterToVisualize(2);

    if (setValidTypes.count(paramType1)>0 && setValidTypes.count(paramType2)>0 && setValidTypes.count(paramType3)>0)
    {
      vtkSmartPointer<cfgPlotPoints3D> plot = createPointPlot3D(table, parameterTypeToString(paramType1), parameterTypeToString(paramType2), parameterTypeToString(paramType3), colors, m_smallDotSize);
      m_plotsPerLevel[ilevel].push_back(plot);
    }

    if (0)
    {
      int nTargetParticules = 100;
      cfgScalar minRadius = 0.1;
      std::vector<int> newparticules3;
      Resampler resampler;
      //resampler.resampleBoundary(minRadius, paramdim, physicalParametersPerLevel[ilevel], distances, nTargetParticules, nTargetParticules);
      int nTargetParticules2 = 10000;
      resampler.resampleUsingNormalDistribution(m_distancesToBoundary[ilevel], nTargetParticules2, newparticules3);

      //std::vector<int> newparticules3tmp = newparticules3;
      //newparticules3.clear();
      //int ind = 0;
      //newparticules3.push_back(newparticules3tmp[ind]);

      vtkVector3i green(0, 255, 0);
      vtkSmartPointer<vtkTable> table3 =  createTable(physicalParametersPerLevel[ilevel], levels, paramdim, 1, labels, &newparticules3);
      vtkSmartPointer<cfgPlotPoints3D> plot = createPointPlot3D(table3, parameterTypeToString(paramType1), parameterTypeToString(paramType2), parameterTypeToString(paramType3), green, m_largeDotSize);
      //m_plotsPerLevel[ilevel].push_back(plot3);

      /*vtkVector3i magenta(255, 0, 255);
      vtkSmartPointer<vtkTable> table4 =  createTable(newPoints, levels, paramdim, 1, labels, &newparticules3);
      vtkSmartPointer<cfgPlotPoints3D> plot4 = createPointPlot3D(table4, "Y1", "Nu1", "Density", magenta, 10);
      //m_plotsPerLevel[ilevel].push_back(plot4); */ 
    }
    if (0)
    {
      int nTargetParticules = 1000;
      cfgScalar minRadius = 0.1;
      std::vector<int> newparticules;
      Resampler resampler;
      resampler.resampleBoundary(minRadius, paramdim, physicalParametersPerLevel[ilevel], m_distancesToBoundary[ilevel], nTargetParticules, newparticules);
      //cfgUtil::writeBinary<int>("boundaryPoints.bin", newparticules);

      vtkVector3i green(0, 255, 0);
      vtkSmartPointer<vtkTable> table2 =  createTable(physicalParametersPerLevel[ilevel], levels, paramdim, 1, labels, &newparticules);
      vtkSmartPointer<cfgPlotPoints3D> plot = createPointPlot3D(table2, parameterTypeToString(paramType1), parameterTypeToString(paramType2), parameterTypeToString(paramType3), green, m_largeDotSize);
      m_plotsPerLevel[ilevel].push_back(plot);
    }
    /*if (0)
    {
      ScoringFunction scoringFunction(paramdim);
      //scoringFunction.setDensityRadius(0.083);
      bool useDistanceField = true;
      scoringFunction.setUseDistanceField(useDistanceField);
      std::vector<cfgScalar> scores = scoringFunction.computeScores(physicalParametersPerLevel[ilevel]);
      cfgScalar maxValue = log(*max_element(scores.begin(), scores.end())+1);
      std::vector<vtkVector3i> colors;
      for (int ipoint=0; ipoint<npoints; ipoint++)
      {
        vtkVector3i matColor;
        cfgScalar w = log(scores[ipoint]+1)/maxValue;
        matColor = vtkVector3i(w*255 + (1-w)*255, (1-w)*255, (1-w)*255);
        colors.push_back(matColor);
      }
      std::vector<int> newparticules3;
      newparticules3= genIncrementalSequence(0, npoints-1);
      vtkSmartPointer<vtkTable> table3 =  createTable(physicalParametersPerLevel[ilevel], levels, paramdim, 1, labels, &newparticules3);
      vtkSmartPointer<cfgPlotPoints3D> plot3 = createPointPlot3D(table3, "Y1", "Nu1", "Density", colors, m_smallDotSize);
      m_plotsPerLevel[ilevel].push_back(plot3);

      cfgScalar minRadius = 0; //0.05;
      int nTargetParticules = 300;
      std::vector<int> newparticules;
      Resampler resampler;
      resampler.resample(scores, nTargetParticules, newparticules);
      resampler.resample(minRadius, paramdim, physicalParametersPerLevel[ilevel], scores, nTargetParticules, newparticules);

      int ind = newparticules[0];
      //newparticules.clear();
      //newparticules.push_back(ind);

      vtkVector3i red(255, 0, 0);
      vtkVector3i green(0, 255, 0);
      vtkSmartPointer<vtkTable> table2 =  createTable(physicalParametersPerLevel[ilevel], levels, paramdim, 1, labels, &newparticules);
      vtkSmartPointer<cfgPlotPoints3D> plot2 = createPointPlot3D(table2, "Y1", "Nu1", "Density", green, m_largeDotSize);
      //m_plotsPerLevel[ilevel].push_back(plot2);
    }*/ 

    //const std::vector<std::vector<cfgScalar> > & thermalExpansionCoeffs = m_project->getThermalExpansionCoeffs();
    /*if (0) //strengths[ilevel].size()>0)
    {
      cfgScalar coeffMin = *std::min_element(strengths[ilevel].begin(), strengths[ilevel].end());
      cfgScalar coeffMax = *std::max_element(strengths[ilevel].begin(), strengths[ilevel].end());
      cfgScalar maxVal = 1; //20//50;
      cfgScalar minVal = 0;
      if (coeffMax>maxVal)
        coeffMax = maxVal;
      for (int ipoint=0; ipoint<npoints; ipoint++)
      {
        vtkVector3i matColor;
        cfgScalar val = strengths[ilevel][ipoint];

        if (ipoint % 100==0)
          std::cout << "ipoint " << ipoint << std::endl;
        int n = m_project->getLevels()[ilevel];
        //std::vector<int> matAssignment = m_project->getMaterialAssignments()[ilevel][ipoint];
        //bool modified=false;
        //bool resOk = filterOutNonConnectedComponents(n, n, n, matAssignment, modified);
        //if (!resOk)
        //{
        //  val = coeffMin;
        //}

        if (val > 1.01)
        {
          matColor = vtkVector3i(0,0,255);
        }
        else
        {
          if (val > coeffMax)
            val = coeffMax;
          if (val <0)
            val = 0;
          //cfgScalar w = (val-coeffMin)/(coeffMax-coeffMin);
          if (val < 0)
          {
            cfgScalar w = log(-val+1)/log(-coeffMin+1);
            matColor = vtkVector3i((1-w)*255, (1-w)*255, w*255 + (1-w)*255);
          }
          else
          {
            cfgScalar w = log(val+1)/log(coeffMax+1);
            //cfgScalar w = (val-minVal)/(coeffMax-minVal);
            matColor = vtkVector3i(w*255 + (1-w)*255, (1-w)*255, (1-w)*255);
          }
        }
        colors.push_back(matColor);
      }
    }*/ 
  }

  /*std::vector<cfgScalar> lengths(paramdim, 1);
  for (ilevel=0; ilevel<nlevel; ilevel++)
  {
    bool resOk = rescaleData(physicalParametersPerLevel[ilevel], paramdim, lengths);
  }
  m_project->setPhysicalParameters(physicalParametersPerLevel);*/ 

  updateChart();
}


void MaterialParametersView::setProject(QPointer<exProject> iProject)
{
  m_project = iProject;
  Q_ASSERT(m_project);

  connect(m_project, SIGNAL(levelVisibilityModified()), this, SLOT(onLevelVisibilityModified()));
  connect(m_project, SIGNAL(levelsModified()), this, SLOT(onLevelsModified()));
  connect(m_project, SIGNAL(pickedReducedPointModified()), this, SLOT(onReducedPointModified()));
  
  GetRenderWindow();

  m_vtkView = vtkContextView::New();
  m_vtkView->SetInteractor(GetInteractor());
  SetRenderWindow(m_vtkView->GetRenderWindow());

  vtkRenderer * renderer = m_vtkView->GetRenderer();
  renderer->SetBackground(1,1,1);

  int width = this->width();
  int height = this->height();
  int size = std::min(width, height);
  int x = width/2 - size/2;
  int y = height/2 - size/2;

  vtkSmartPointer<cfgChartXYZ> chart = vtkSmartPointer<cfgChartXYZ>::New();
  chart->SetGeometry(vtkRectf(x, y, size, size));
  for (int i=0; i<3; i++)
  {
    chart->GetAxis(i)->SetPrecision(4);
  }

  m_vtkView->GetScene()->AddItem(chart.GetPointer()); 
  m_chart = chart;
  updateChart();

  m_vtkView->GetRenderWindow()->SetMultiSamples(10);
  m_vtkView->GetRenderWindow()->Render();
}

void MaterialParametersView::updateChart()
{
  m_chart->ClearPlots();

  m_plots2Levels.clear();
  m_plots2PlotIndex.clear();

  int ilevel, nlevel=(int)m_plotsPerLevel.size();
  for (ilevel=0; ilevel<nlevel; ilevel++)
  {
    bool levelVisiblity = m_project->getLevelVisibility(ilevel);
    if (levelVisiblity==true)
    {
      int iplot, nplot=(int)m_plotsPerLevel[ilevel].size();
      for (iplot=0; iplot<nplot; iplot++)
      {
        m_chart->AddPlot(m_plotsPerLevel[ilevel][iplot]);
        m_plots2Levels.push_back(ilevel);
        m_plots2PlotIndex.push_back(iplot);
      }
    }
  }
  if (m_auxiliaryPlot.Get()!=NULL)
  {
    m_chart->AddPlot(m_auxiliaryPlot);
  }
  update();
}

int MaterialParametersView::getStructureIndex(int iPointIndex)
{
  int shift=0;
  int ilevel=0, nlevel=(int)m_lastPointIndices.size();
  for (ilevel=0; ilevel<nlevel; ++ilevel)
  {
    int lastPointIndex = m_lastPointIndices[ilevel];
    if (iPointIndex > lastPointIndex)
    {
      shift = m_lastPointIndices[ilevel]+1;
    }
    else
    {
      break;
    }
  }
  int indStructure = iPointIndex-shift;
  return indStructure;
}

int MaterialParametersView::getParameterIndex(int iStructureIndex, int iLevel)
{
  int matIndex = iStructureIndex;
  if (iLevel>0)
  {
    matIndex += m_lastPointIndices[iLevel-1]+1;
  }
  return matIndex;
}

void MaterialParametersView::displayPoint(int iPointIndex, int iLevel, const vtkVector3i &iColor)
{  
  int indPlot = 0;
  cfgPlotPoints3D * plot = (cfgPlotPoints3D*)(m_plotsPerLevel[iLevel][indPlot].Get());

  std::vector<int> rowIndices;
  rowIndices.push_back(iPointIndex);

  vtkSmartPointer<vtkTable> table = QVTKWidgetUtilities::createTable(plot->getTable(), rowIndices);
  m_auxiliaryPlot = createPointPlot3D(table, plot->GetXAxisLabel(), plot->GetYAxisLabel(), plot->GetZAxisLabel(), iColor, m_largeDotSize);

  updateChart();
  update();
}

void MaterialParametersView::changePointsColor(const std::vector<int> &iPointIndices, int iIndPlot, const vtkVector3i &iColor)
{
  int npreviousPoints = (int)m_savedPointIndices.size();
  for (int ipoint=0; ipoint<npreviousPoints; ipoint++)
  {
    assert(m_savedPlot);
    m_savedPlot->setColor(m_savedPointIndices[ipoint], m_savedColors[ipoint]);
  }
  m_savedPointIndices.clear();
  m_savedColors.clear();
  m_savedPlot = NULL;

  int indPlot = m_plots2PlotIndex[iIndPlot];
  int indLevel = m_plots2Levels[iIndPlot];
  cfgPlotPoints3D * plot = (cfgPlotPoints3D*)(m_plotsPerLevel[indLevel][indPlot].Get());
  m_savedPlot = plot;

  int npoint=(int)iPointIndices.size();
  for (int ipoint=0; ipoint<npoint; ipoint++)
  {
    int indPoint = iPointIndices[ipoint];
    m_savedPointIndices.push_back(indPoint);
    m_savedColors.push_back(plot->getColor(indPoint));
    plot->setColor(indPoint, iColor);
  }
  m_chart->Update();
  update();
}

void MaterialParametersView::changePointsColor(const std::vector<int> &iPointIndices, int iIndPlot, const std::vector<vtkVector3i> &iColors)
{
  int npreviousPoints = (int)m_savedPointIndices.size();
  for (int ipoint=0; ipoint<npreviousPoints; ipoint++)
  {
    assert(m_savedPlot);
    m_savedPlot->setColor(m_savedPointIndices[ipoint], m_savedColors[ipoint]);
  }
  m_savedPointIndices.clear();
  m_savedColors.clear();
  m_savedPlot = NULL;

  int indPlot = m_plots2PlotIndex[iIndPlot];
  int indLevel = m_plots2Levels[iIndPlot];
  cfgPlotPoints3D * plot = (cfgPlotPoints3D*)(m_plotsPerLevel[indLevel][indPlot].Get());
  m_savedPlot = plot;

  int npoint=(int)iPointIndices.size();
  for (int ipoint=0; ipoint<npoint; ipoint++)
  {
    int indPoint = iPointIndices[ipoint];
    m_savedPointIndices.push_back(indPoint);
    m_savedColors.push_back(plot->getColor(indPoint));
    plot->setColor(indPoint, iColors[ipoint]);
  }
  m_chart->Update();
  update();
}

void MaterialParametersView::mouseReleaseEvent(QMouseEvent * iMouseEvent)
{  
  QVTKWidget::mouseReleaseEvent(iMouseEvent);
  int pickedPointIndex = m_chart->pickedPointIndex();
  int pickedPlotIndex = m_chart->pickedPlotIndex();
  if (pickedPointIndex>=0)
  {
    int indLevel = m_plots2Levels[pickedPlotIndex];
    m_project->setPickedStructure(pickedPointIndex, indLevel);

    std::cout << "ActiveTool = " << m_project->getActiveTool() << std::endl;

    if (m_project->getActiveTool()==ex::RegionSelectionToolType)
    {
      int paramDim = getParameterDim();
      cfgScalar maxRadius = m_samplerRadius;
      cfgScalar ratioDist = 0.2;
      int nbPoints = m_samplerNbPoints;
      std::vector<int> pointIndices;
      m_resamplers[pickedPlotIndex]->resampleBoundary(pickedPointIndex, 0, maxRadius, nbPoints, pointIndices);
      //m_resamplers[pickedPlotIndex]->resampleBoundary(pickedPointIndex, 0, maxRadius, ratioDist, pointIndices);
      m_project->setSelectedPoints(pointIndices, pickedPlotIndex);

      if (pointIndices.size()>0)
      {
        std::string fileNameIndices = "..//..//Output//microstructureIndices.txt";
        cfgUtil::writeVector2File<int>(pointIndices, fileNameIndices); 

        const std::vector<std::vector<cfgScalar> > & parameters = m_project->getPhysicalParameters();
        std::vector<cfgScalar> selectedMicrostructuresParameters = getSubVector(parameters[indLevel], paramDim, pointIndices);
        std::string fileNameParameters = "..//..//Output//microstructureParameters.txt";
        cfgUtil::writeVector2File<cfgScalar>(selectedMicrostructuresParameters, fileNameParameters); 

        const std::vector<std::vector<std::vector<int> > > & matAssignments = m_project->getMaterialAssignments();
        std::vector<int> selectedMicrostructuresMatAssignments = toStdVector(getSubVector(matAssignments[indLevel], pointIndices));
        std::string fileNameMaterials = "..//..//Output//microstructureMaterials.txt";
        cfgUtil::writeVector2File<int>(selectedMicrostructuresMatAssignments, fileNameMaterials); 
      }
      vtkVector3i color(0, 255, 0);
      changePointsColor(pointIndices,pickedPlotIndex, color);
    }
    const std::vector<std::vector<cfgScalar> > & tensors = m_project->getElasticityTensors();
    if (tensors.size()>0 && tensors[indLevel].size()>0)
    {
      int indStart = 6*pickedPointIndex;
      const std::vector<cfgScalar> & t = tensors[indLevel];
      std::cout << "Elasticity tensor: " << std::endl;
      std::cout << t[indStart] << std::endl;
      std::cout << t[indStart+1] << " " << t[indStart+2] << std::endl;
      std::cout << t[indStart+3] << " " << t[indStart+4] << " " << t[indStart+5] << std::endl;
    }
  }

  /*if (pickedPointIndex>=0)
  {
    //int pickedStructureLevel = m_tablePoint3D->GetValue(pickedPointIndex, 3).ToInt();
    //int pickedStructureIndex = getStructureIndex(pickedPointIndex);

    int pickedPlotIndex = m_chart->pickedPlotIndex();
    cout << "pickedPlotIndex " << pickedPlotIndex << std::endl;
    int pickedStructureLevel = m_plots2Levels[pickedPlotIndex];
    cout << "pickedStructureLevel " << pickedStructureLevel << std::endl;
    int indPlot = m_plots2PlotIndex[pickedPlotIndex];
    cout << "indPlot " << indPlot << std::endl;
    int pickedStructureIndex = m_plotPointIndices[pickedStructureLevel][indPlot][pickedPointIndex];
    pickedStructureIndex = getStructureIndex(pickedStructureIndex);
    cout << "picked structure : " << pickedStructureIndex << std::endl;

    cout << m_plotPointIndices.size() << std::endl;
    cout << m_plotPointIndices[pickedStructureLevel].size() << std::endl;
    cout << m_plotPointIndices[pickedStructureLevel][indPlot].size() << std::endl;

    std::vector<int> pickedStructurebaseMaterials;
    m_project->setPickedStructure(pickedStructureIndex, pickedStructureLevel, pickedStructurebaseMaterials);
    cout << "picked structure : " << pickedStructureIndex << std::endl;

    int prevLevelIndex = pickedStructureLevel-1;
    
    const std::vector<int> & materials = m_materialAssignements[pickedStructureLevel][pickedStructureIndex];
    int imat=0, nmat=(int)materials.size();
    for (imat=0; imat<nmat; imat++)
    {
      int indMat = materials[imat];
      const std::vector<int> & baseMaterials = m_baseMaterials[pickedStructureLevel][indMat];
      int baseMaterialIndex = m_baseMaterials2Indices[pickedStructureLevel][baseMaterials];
      int pointIndex = getParameterIndex(baseMaterialIndex, prevLevelIndex);
      pickedStructurebaseMaterials.push_back(pointIndex);
    }
    //highlighPoints(pickedStructurebaseMaterials);
    m_project->setPickedStructure(pickedStructureIndex, pickedStructureLevel, pickedStructurebaseMaterials); 
  } */ 
}

void MaterialParametersView::highlighPoints(const std::vector<int> &iPointIndices)
{
  /*if (m_plotHighLightedPoints.Get() == NULL)
  {
    std::string labels[4]= {"Density", "Y1", "Y2", "Level"};
    m_tableHighLightedPoints =  createTable(m_physicalParameters, m_levels, 3, 1, labels, &iPointIndices);
    m_plotHighLightedPoints = createPointPlot3D(m_tableHighLightedPoints, "Y1", "Y2", "Density", vtkVector3i(0,0,0), m_largeDotSize);
    m_chart->AddPlot(m_plotHighLightedPoints);
  }
  else
  {
    std::string labels[4]= {"Density", "Y1", "Y2", "Level"};
    m_tableHighLightedPoints = createTable(m_physicalParameters, m_levels, 3, 1, labels, &iPointIndices);
    m_plotHighLightedPoints->SetInputData(m_tableHighLightedPoints, "Y1", "Y2", "Density");
  }*/ 
}

void MaterialParametersView::onReducedPointModified()
{
  int indPoint = m_project->getPickedStructureIndex();
  if (indPoint>=0)
  {
    int indLevel = m_project->getPickedStructureLevel();

    vtkVector3i col(255, 0, 255);
    displayPoint(indPoint, indLevel, col);
  }
}

void MaterialParametersView::onLevelVisibilityModified()
{
  updateChart();
}

void MaterialParametersView::onLevelsModified()
{
  updatePlots();
}

void MaterialParametersView::onParamsToVisualizeModified()
{
  updatePlots();
}

void MaterialParametersView::onRadiusValueChanged(double iValue)
{
  m_samplerRadius = iValue;
}

void MaterialParametersView::onNbPointsValueChanged(int iValue)
{
  m_samplerNbPoints = iValue;
}



