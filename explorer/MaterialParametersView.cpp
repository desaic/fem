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

#include "exProject.h"

#include <cfgMaterialUtilities.h>
using namespace cfgMaterialUtilities;

#include <cfgUtilities.h>
using namespace cfgUtil;

#include "ScoringFunction.h"
#include "Resampler.h"

#include "DistanceField.h"

MaterialParametersView::MaterialParametersView()
  :QVTKWidget()
{
}

MaterialParametersView::~MaterialParametersView()
{
}

vtkSmartPointer<vtkTable> MaterialParametersView::createTable(const std::vector<float> &iValues1, int idim1, const std::vector<std::string> &iLabels, const std::vector<int> *iPointIndices)
{
  vtkSmartPointer<vtkTable> table = vtkSmartPointer<vtkTable>::New();

  assert(idim1>=3);
  int npoint = (int)(iPointIndices? iPointIndices->size(): iValues1.size()/idim1);
  int icol=0;
  for (icol=0; icol<idim1; icol++)
  {
    table->AddColumn(vtkSmartPointer<vtkFloatArray>::New());
    table->GetColumn(icol)->SetName(iLabels[icol].c_str());
  }
  table->SetNumberOfRows(npoint);

  assert(iValues1.size()%idim1==0);
  int ipoint=0;
  for (ipoint=0; ipoint<npoint; ipoint++)
  {
    int indPoint = (iPointIndices? (*iPointIndices)[ipoint]: ipoint);
    for (icol=0; icol<idim1; icol++)
    {
      table->SetValue(ipoint, icol,  iValues1[idim1*indPoint + icol]);
    }
  }
  return table;
}

vtkSmartPointer<vtkTable> MaterialParametersView::createTable(const std::vector<float> &iValues1, const std::vector<int> &iValues2, int idim1, int idim2, const std::vector<std::string> &iLabels, const std::vector<int> *iPointIndices)
{
  vtkSmartPointer<vtkTable> table = vtkSmartPointer<vtkTable>::New();

  assert(idim1>=3);
  int npoint = (int)(iPointIndices? iPointIndices->size(): iValues1.size()/idim1);
  int icol=0;
  for (icol=0; icol<idim1; icol++)
  {
    table->AddColumn(vtkSmartPointer<vtkFloatArray>::New());
    table->GetColumn(icol)->SetName(iLabels[icol].c_str());
  }
  for (icol=0; icol<idim2; icol++)
  {
    table->AddColumn(vtkSmartPointer<vtkIntArray>::New());
    table->GetColumn(icol+idim1)->SetName(iLabels[icol+idim1].c_str());
  }
  table->SetNumberOfRows(npoint);

  assert(iValues1.size()%idim1==0 && iValues2.size()%idim2==0 && (iValues2.size()>0?  iValues1.size()/idim1==iValues2.size()/idim2: true));
  int ipoint=0;
  for (ipoint=0; ipoint<npoint; ipoint++)
  {
    int indPoint = (iPointIndices? (*iPointIndices)[ipoint]: ipoint);
    for (icol=0; icol<idim1; icol++)
    {
      table->SetValue(ipoint, icol,  iValues1[idim1*indPoint + icol]);
    }
    for (icol=0; icol<idim2;  icol++)
    {
      table->SetValue(ipoint, icol+idim1,  iValues2[idim2*indPoint + icol]);
    }
  }
  return table;
}

vtkSmartPointer<cfgPlotPoints3D> MaterialParametersView::createPointPlot3D(vtkSmartPointer<vtkTable> &iTable, const std::string &iLabel1, const std::string &iLabel2, const std::string &iLabel3, std::vector<vtkVector3i> colors, float iWidth)
{
  vtkSmartPointer<cfgPlotPoints3D> plot = vtkSmartPointer<cfgPlotPoints3D>::New();
  plot->SetInputData(iTable, iLabel1, iLabel2, iLabel3);
  plot->SetColors(colors);
  plot->GetPen()->SetWidth(iWidth);
  return plot;
}

vtkSmartPointer<cfgPlotPoints3D> MaterialParametersView::createPointPlot3D(vtkSmartPointer<vtkTable> &iTable, const std::string &iLabel1, const std::string &iLabel2, const std::string &iLabel3, vtkVector3i &iColor, float iWidth)
{
  vtkSmartPointer<cfgPlotPoints3D> plot = vtkSmartPointer<cfgPlotPoints3D>::New();
  plot->SetInputData(iTable, iLabel1, iLabel2, iLabel3);
  plot->GetPen()->SetColor(iColor[0], iColor[1], iColor[2]);
  plot->GetPen()->SetWidth(iWidth);
  return plot;
}

vtkSmartPointer<cfgPlotLine3D> MaterialParametersView::createLinePlot3D(vtkSmartPointer<vtkTable> &iTable, const std::string &iLabel1, const std::string &iLabel2, const std::string &iLabel3, vtkVector3i &iColor, float iWidth)
{
  vtkSmartPointer<cfgPlotLine3D> plot = vtkSmartPointer<cfgPlotLine3D>::New();
  plot->SetInputData(iTable, iLabel1, iLabel2, iLabel3);
  plot->GetPen()->SetColor(iColor[0], iColor[1], iColor[2]);
  plot->GetPen()->SetWidth(iWidth);
  return plot;
}

vtkSmartPointer<cfgPlotSurface> MaterialParametersView::createSurfacePlot3D(vtkSmartPointer<vtkTable> &iTable, const std::string &iLabel1, const std::string &iLabel2, const std::string &iLabel3, vtkVector3i &iColor, int iAlpha)
{
  vtkSmartPointer<cfgPlotSurface> plot = vtkSmartPointer<cfgPlotSurface>::New();
  plot->SetInputData(iTable, iLabel1, iLabel2, iLabel3);
  plot->GetPen()->SetColor(iColor[0], iColor[1], iColor[2]);
  plot->GetPen()->SetOpacity(iAlpha);
  return plot;
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

void MaterialParametersView::updatePlots()
{
  Q_ASSERT(m_project);

  const std::vector<std::vector<cfgScalar> > & physicalParametersPerLevel = m_project->getPhysicalParameters();
  int ilevel, nlevel=(int)physicalParametersPerLevel.size();

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

  m_plotsPerLevel.clear();
  m_plotsPerLevel.resize(nlevel);

  //int paramdim = m_project->getDim()+1;
  //int paramdim = 2*m_project->getDim()+1;
  int paramdim = 0;
  bool cubic = true;
  bool orthotropic = true;

  std::vector<std::string> labels;
  if (paramdim==3)
  {
    labels.push_back("Density");
    labels.push_back("Y1");
    labels.push_back("Y2");
    labels.push_back("Level");
  }
  else if (cubic)
  {
    paramdim = 4;
    labels.push_back("Density");
    labels.push_back("Y1");
    labels.push_back("Nu1");
    labels.push_back("Mu1");
    labels.push_back("Level");
  }
  else if (orthotropic)
  {
    if (m_project->getDim()==2)
    {
      paramdim = 5;
      labels.push_back("Density");
      labels.push_back("Y1");
      labels.push_back("Y2");
      labels.push_back("Nu1");
      labels.push_back("Mu1");
      labels.push_back("Level");
    }
    else
    {
      paramdim = 10;
      labels.push_back("Density");
      labels.push_back("Y1");
      labels.push_back("Y2");
      labels.push_back("Y3");
      labels.push_back("Nu1"); //12
      labels.push_back("Nu2"); //13
      labels.push_back("Nu3"); //23
      labels.push_back("Mu1");
      labels.push_back("Mu2");
      labels.push_back("Mu3");
      labels.push_back("Level");
    }
  }
  else if (paramdim==5 && cubic)
  {
    labels.push_back("Density");
    labels.push_back("Y1");
    labels.push_back("Nu1");
    labels.push_back("Mu1");
    labels.push_back("Level");
  }
  else if (paramdim==7)
  {
    labels.push_back("Density");
    labels.push_back("Y1");
    labels.push_back("Y2");
    labels.push_back("Y3");
    labels.push_back("Nu1");
    labels.push_back("Nu2");
    labels.push_back("Nu3");
    labels.push_back("Level");
  }
  for (ilevel=0; ilevel<nlevel; ilevel++)
  {
    vtkVector3i col = matColors[ilevel];

    int npoints = (int)physicalParametersPerLevel[ilevel].size()/paramdim;
    std::vector<int> levels(npoints, ilevel);
    vtkSmartPointer<vtkTable> table =  createTable(physicalParametersPerLevel[ilevel], levels, paramdim, 1, labels);

    if (0)
    {
      DistanceField distanceField(paramdim);
      std::vector<cfgScalar> distances = distanceField.computeDistances(physicalParametersPerLevel[ilevel]);

      cfgScalar minValue = FLT_MAX;
      cfgScalar maxValue = -FLT_MAX;
      for (int ipoint=0; ipoint<npoints; ipoint++)
      {
        cfgScalar dist = distances[ipoint];
        if (dist < minValue)
          minValue = dist;
        if (dist > maxValue)
          maxValue = dist;
      }
      std::cout << "min dist = " << minValue << ", max dist = " << maxValue << std::endl;

      std::vector<vtkVector3i> colors;
      for (int ipoint=0; ipoint<npoints; ipoint++)
      {
        vtkVector3i matColor;
        if (distances[ipoint] < 0)
        {
          cfgScalar w = distances[ipoint]/minValue;
          matColor = vtkVector3i((1-w)*255, (1-w)*255, w*255 + (1-w)*255);
        }
        else if (distances[ipoint] > 0)
        {
          cfgScalar w = distances[ipoint]/maxValue;
          matColor = vtkVector3i(w*255 + (1-w)*255, (1-w)*255, (1-w)*255);
        }
        else
        {
           matColor = vtkVector3i(255, 0, 0);
        }
        colors.push_back(matColor);
      }
      vtkSmartPointer<cfgPlotPoints3D> plot = createPointPlot3D(table, "Y1", "Nu1", "Density", colors, 10);
      //m_plotsPerLevel[ilevel].push_back(plot);
    }

    if (0)
    {
      ScoringFunction scoringFunction(paramdim);
      //scoringFunction.setDensityRadius(0.083);
      std::vector<cfgScalar> scores = scoringFunction.computeScores(physicalParametersPerLevel[ilevel]);

      int nTargetParticules = 300;
      std::vector<int> newparticules;
      Resampler resampler;
      resampler.resample(scores, nTargetParticules, newparticules);

      vtkVector3i red(255, 0, 0);
      vtkSmartPointer<vtkTable> table2 =  createTable(physicalParametersPerLevel[ilevel], levels, paramdim, 1, labels, &newparticules);
      vtkSmartPointer<cfgPlotPoints3D> plot2 = createPointPlot3D(table2, "Y1", "Y2", "Density", red, 15);
      //m_plotsPerLevel[ilevel].push_back(plot2);
    }
   
    std::vector<vtkVector3i> colors;
    for (int ipoint=0; ipoint<npoints; ipoint++)
    {
      cfgScalar w = 1 ;//scores[ipoint];
      vtkVector3i matColor(w*col[0] + (1-w)*255, w*col[1] + (1-w)*255, w*col[2] + (1-w)*255);
      colors.push_back(matColor);
    }
    //vtkSmartPointer<cfgPlotPoints3D> plot = createPointPlot3D(table, "Y1", "Y2", "Density", colors, 10);
    //vtkSmartPointer<cfgPlotPoints3D> plot = createPointPlot3D(table, "Nu1", "Nu2", "Density", colors, 10);
    vtkSmartPointer<cfgPlotPoints3D> plot = createPointPlot3D(table, "Y1", "Nu1", "Density", colors, 10);
    m_plotsPerLevel[ilevel].push_back(plot);
  }
  updateChart();
}

void MaterialParametersView::setProject(QPointer<exProject> iProject)
{
  m_project = iProject;
  Q_ASSERT(m_project);

  connect(m_project, SIGNAL(levelVisibilityModified()), this, SLOT(onLevelVisibilityModified()));
  connect(m_project, SIGNAL(levelsModified()), this, SLOT(onLevelsModified()));

#if 0
  // read data
  int maxlevel = m_project->getMaxLevel();
  int nlevel = maxlevel+1;
  std::vector<std::vector<cfgScalar> > physicalParametersPerLevel(nlevel);
  std::vector<cfgScalar> physicalParameters;

  m_baseMaterials.clear();
  m_baseMaterials.resize(nlevel);

  m_baseMaterials2Indices.clear();
  m_baseMaterials2Indices.resize(nlevel);

  m_materialAssignements.clear();
  m_materialAssignements.resize(nlevel);

  m_plotsPerLevel.resize(nlevel);
  m_plotPointIndices.resize(nlevel);

  int nparam = 3;

  bool writeBinary=false;
  bool readBinary=!writeBinary;
  int nextralevel =  13 ;//14; //8;

  int ilevel=0;
  for (ilevel=1; ilevel<nlevel; ilevel++)
  {
    bool ResOk = true;
    std::string fileRootName = m_project->getFileDirectory() + "level" + std::to_string(ilevel) + "_";
    std::string fileExtension = ".bin";
    if (readBinary && ilevel<nlevel)
    {
      ResOk = cfgUtil::readBinary<float>(fileRootName + "params" + fileExtension, physicalParametersPerLevel[ilevel]);
      if (ResOk)
        ResOk = cfgUtil::readBinary<int>(fileRootName + "baseMat" + fileExtension, m_baseMaterials[ilevel]);
      if (ResOk)
      ResOk = cfgUtil::readBinary<int>(fileRootName + "matAssignments" + fileExtension, m_materialAssignements[ilevel]);
    }
    else
    {
      ResOk = m_project->getMaterialParameters(ilevel, physicalParametersPerLevel[ilevel], m_baseMaterials[ilevel], m_materialAssignements[ilevel]);
      if (ResOk /*&& writeBinary*/ )
      {
        cfgUtil::writeBinary<float>(fileRootName + "params" + fileExtension, physicalParametersPerLevel[ilevel]);
        cfgUtil::writeBinary<int>(fileRootName + "baseMat" + fileExtension, m_baseMaterials[ilevel]);
        cfgUtil::writeBinary<int>(fileRootName + "matAssignments" + fileExtension, m_materialAssignements[ilevel]);
      }
    }
    if (!ResOk)
    {
      std::cout << "Cannot read file for level " << ilevel << std::endl;
    }
    m_baseMaterials2Indices[ilevel] = computeVector2IndexMap(m_baseMaterials[ilevel]);
  }
 
  m_baseMaterials[0] = m_baseMaterials[1];
  std::vector<int> mat1, mat2;
  mat1.push_back(0);
  mat2.push_back(1);
  m_materialAssignements[0].push_back(mat1);
  m_materialAssignements[0].push_back(mat2);
  int iparam=0;
  for (iparam=0; iparam<nparam; iparam++)
  {
    physicalParametersPerLevel[0].push_back(physicalParametersPerLevel[1][iparam]);
  }
  int shift = physicalParametersPerLevel[1].size()-nparam;
  for (iparam=0; iparam<nparam; iparam++)
  {
    physicalParametersPerLevel[0].push_back(physicalParametersPerLevel[1][iparam+shift]);
  }

  /*m_materialAssignements[1] = m_materialAssignements[0];
  physicalParametersPerLevel[1] = physicalParametersPerLevel[0];
  std::string fileRootName = m_project->getFileDirectory() + "level" + std::to_string(1) + "_";
  std::string fileExtension = ".bin";
  cfgUtil::writeBinary<float>(fileRootName + "params" + fileExtension, physicalParametersPerLevel[1]);
  cfgUtil::writeBinary<int>(fileRootName + "baseMat" + fileExtension, m_baseMaterials[1]);
  cfgUtil::writeBinary<int>(fileRootName + "matAssignments" + fileExtension, m_materialAssignements[1]);*/ 

  for (ilevel=0; ilevel<nlevel; ilevel++)
  {
    physicalParameters.insert(physicalParameters.end(), physicalParametersPerLevel[ilevel].begin(), physicalParametersPerLevel[ilevel].end());
  }
  m_physicalParameters = physicalParameters;

  std::vector<vtkVector3i> colors;
  std::vector<int> levels;

  std::vector<vtkVector3i> matColors;
  /*matColors.push_back(vtkVector3i(0, 0, 0));
  matColors.push_back(vtkVector3i(255, 0, 0));
  matColors.push_back(vtkVector3i(0, 0, 255));
  matColors.push_back(vtkVector3i(255, 0, 255)); 
  matColors.push_back(vtkVector3i(255, 255, 0)); */ 
  for (ilevel=0; ilevel<nlevel; ilevel++)
  {
    QColor colHSV = QColor::fromHsv(ilevel*255/(nlevel-1), 255, 255);
    QColor colRGB = colHSV.toRgb();

    matColors.push_back(vtkVector3i(colRGB.red(), colRGB.green(), colRGB.blue()));
    //matColors.push_back(vtkVector3i(0, 0, ilevel*255/nlevel));
  }
  
  int nrow=0;
  for (ilevel=0; ilevel<nlevel; ilevel++)
  {
    nrow += (int)physicalParametersPerLevel[ilevel].size()/nparam;
    colors.resize(nrow, matColors[ilevel]);
    levels.resize(nrow, ilevel);
    m_lastPointIndices.push_back(nrow-1);
  }
  m_levels = levels;

  std::string labels[4]= {"Density", "Y1", "Y2", "Level"};
  vtkSmartPointer<vtkTable> table =  createTable(physicalParameters, levels, 3, 1, labels);
  vtkSmartPointer<cfgPlotPoints3D> points3D = createPointPlot3D(table, "Y1", "Y2", "Density", colors);
  m_tablePoint3D = table;

  for (ilevel=0; ilevel<nlevel; ilevel++)
  {
    bool computeTriangulation = false; //(ilevel>1);
    std::vector<int> convexHull, boundaryVertices, boundaryFaces;
    if (computeTriangulation && ilevel<nlevel-nextralevel)
    {
      //computeBoundaryPoints(physicalParametersPerLevel[0], 3, convexHull);
      //computeDelaundayTriangulation(physicalParametersPerLevel[1], 3, convexHull, boundaryVertices, boundaryFaces);

      /*computeConvexHull(physicalParametersPerLevel[indLevelBoundary], nparam, boundaryVertices);
      std::vector<int> newBoundaryVertices;
      getClosestPoints(physicalParametersPerLevel[indLevelBoundary], 3, boundaryVertices, 1.e-3, newBoundaryVertices);*/ 

      std::vector<cfgScalar> parameterPoints;

      parameterPoints.insert(parameterPoints.end(), physicalParameters.begin(), physicalParameters.begin()+ 3*m_lastPointIndices[ilevel]+3);
      std::cout << "Computing delaunay triangulation... " << std::endl;
      //computeDelaundayTriangulation(parameterPoints, 3, convexHull, boundaryVertices, boundaryFaces);

      computeDelaundayTriangulation(physicalParametersPerLevel[ilevel], 3, convexHull, boundaryVertices, boundaryFaces);
      boundaryVertices = add(boundaryVertices, m_lastPointIndices[ilevel-1]+1); 
      boundaryFaces = add(boundaryFaces, m_lastPointIndices[ilevel-1]+1); 
      std::cout << "# params = " << physicalParametersPerLevel[ilevel].size() << ", # boundary points = " << boundaryVertices.size() << std::endl;
      //writeVector2File(physicalParametersPerLevel[indLevelBoundary], m_project->getFileDirectory() + "params.m");
    }
    else if (ilevel>0)
    {
      int ipoint, npoint = physicalParametersPerLevel[ilevel].size()/3;
      for (ipoint=0; ipoint<npoint; ipoint++)
      {
        boundaryVertices.push_back(ipoint);
      }
      boundaryVertices = add(boundaryVertices, m_lastPointIndices[ilevel-1]+1); 
    }
    bool computeDistToBoundary = false;
    if (computeDistToBoundary)
    {
      int indLevel = maxlevel;
      std::vector<cfgScalar> sampledBoundary;
      sampleMesh(physicalParameters, boundaryFaces, 1000, sampledBoundary);
      std::vector<cfgScalar> distToBoundary;
      computeDistancesToPointCloud(physicalParametersPerLevel[indLevel], sampledBoundary, distToBoundary);
      cfgScalar maxDist = *std::max_element(distToBoundary.begin(), distToBoundary.end()); 
      distToBoundary = cfgUtil::mult<cfgScalar>(distToBoundary, 1./maxDist);

      shift = m_lastPointIndices[indLevel-1]+1;
      int ipoint=0, npoint=m_lastPointIndices[indLevel]-m_lastPointIndices[indLevel-1];
      for (ipoint=0; ipoint<npoint; ipoint++)
      {
        cfgScalar w = 1-distToBoundary[ipoint];
        vtkVector3i matColor(w*matColors[indLevel][0] + (1-w)*255, w*matColors[indLevel][1] + (1-w)*255, w*matColors[indLevel][2] + (1-w)*255);
        colors[ipoint+shift] = matColor;
      }
    }
    bool displayClusters = false; //ilevel>0;
    if (displayClusters/* && ilevel==5*/)
    {
      int npoints = 100;

      std::vector<cfgScalar> points = physicalParametersPerLevel[ilevel];
      std::vector<cfgScalar> lengths(3, 1);
      rescaleData(points, 3, lengths);

      int nbIter = 10;
      std::vector<std::vector<int> > clusters;
      std::vector<int> pointIndices; 

      getKMeans(nbIter, npoints, points, 3, clusters, &pointIndices);

      int icluster, ncluster=(int)clusters.size();
      for (icluster=0; icluster<ncluster; icluster++)
      {
        std::vector<int> pointIndices = clusters[icluster];
        pointIndices = add(pointIndices, m_lastPointIndices[ilevel-1]+1); 

        vtkSmartPointer<vtkTable> table =  createTable(physicalParameters, levels, 3, 1, labels, &pointIndices);
        vtkVector3i col = matColors[icluster%matColors.size()];
        vtkSmartPointer<cfgPlotPoints3D> plot = createPointPlot3D(table, "Y1", "Y2", "Density", col, 15);
        m_plotsPerLevel[ilevel].push_back(plot);
        m_plotPointIndices[ilevel].push_back(pointIndices);
      }

      if (0)
      {
        //getFurthestPointsGreedy(npoints, points, 3, pointIndices);
        pointIndices = add(pointIndices, m_lastPointIndices[ilevel-1]+1); 

        vtkSmartPointer<vtkTable> table =  createTable(physicalParameters, levels, 3, 1, labels, &pointIndices);
        vtkVector3i col = vtkVector3i(255,0,255);
        vtkSmartPointer<cfgPlotPoints3D> plot = createPointPlot3D(table, "Y1", "Y2", "Density", col, 15);
        m_plotsPerLevel[ilevel].push_back(plot);
        m_plotPointIndices[ilevel].push_back(pointIndices);
      }
    }

    bool displayInsidePoints = false;
    if (displayInsidePoints)
    {

    }

    if (boundaryVertices.size()>0)
    {
      vtkSmartPointer<vtkTable> table =  createTable(physicalParameters, levels, 3, 1, labels, &boundaryVertices);
      vtkVector3i col = matColors[ilevel];
      //vtkVector3i col = vtkVector3i(0,255,0);
      vtkSmartPointer<cfgPlotPoints3D> boundaryPlot = createPointPlot3D(table, "Y1", "Y2", "Density", col, 10);
      m_plotsPerLevel[ilevel].push_back(boundaryPlot);
      m_plotPointIndices[ilevel].push_back(boundaryVertices);
    } 
    if (convexHull.size()>0)
    {
      std::vector<int> edgeIndices;
      //getEdgesFromTetFaceIndexArray(convexHull, edgeIndices);

      getEdgesFromTriFaceIndexArray(boundaryFaces, edgeIndices);

      vtkSmartPointer<vtkTable> edgesTable =  createTable(physicalParameters, levels, 3, 1, labels, &edgeIndices);
      vtkSmartPointer<cfgPlotLine3D> edgesPlot = createLinePlot3D(edgesTable, "Y1", "Y2", "Density", vtkVector3i(0,0,0));
      m_plotsPerLevel[ilevel].push_back(edgesPlot);
      m_plotPointIndices[ilevel].push_back(edgeIndices);

      vtkSmartPointer<vtkTable> facesTable =  createTable(physicalParameters, levels, 3, 1, labels, &boundaryFaces);
      vtkVector3i col = matColors[ilevel];
      //vtkVector3i col = vtkVector3i(255,0,0);
      vtkSmartPointer<cfgPlotSurface> surfacePlot = createSurfacePlot3D(facesTable, "Y1", "Y2", "Density", col, 10);
      m_plotsPerLevel[ilevel].push_back(surfacePlot);
      m_plotPointIndices[ilevel].push_back(boundaryFaces);
    }
    /*if (sampledBoundary.size()>0)
    {
      vtkSmartPointer<vtkTable> sampledBoundaryTable =  createTable(sampledBoundary, 3, labels, NULL);
      vtkSmartPointer<cfgPlotPoints3D> sampledBoundaryPlot = createPointPlot3D(sampledBoundaryTable, "Y1", "Y2", "Density", vtkVector3i(0,0,0), 5);
      m_plotsPerLevel[ilevel].push_back(sampledBoundaryPlot);
    }*/ 
  } 
#endif 

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

void MaterialParametersView::mouseReleaseEvent(QMouseEvent * iMouseEvent)
{  
  QVTKWidget::mouseReleaseEvent(iMouseEvent);
  int pickedPointIndex = m_chart->pickedPointIndex();
  int pickedPlotIndex = m_chart->pickedPlotIndex();
  if (pickedPointIndex>=0)
  {
    int indLevel = m_plots2Levels[pickedPlotIndex];
    m_project->setPickedStructure(pickedPointIndex, indLevel);
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
    m_plotHighLightedPoints = createPointPlot3D(m_tableHighLightedPoints, "Y1", "Y2", "Density", vtkVector3i(0,0,0), 20);
    m_chart->AddPlot(m_plotHighLightedPoints);
  }
  else
  {
    std::string labels[4]= {"Density", "Y1", "Y2", "Level"};
    m_tableHighLightedPoints = createTable(m_physicalParameters, m_levels, 3, 1, labels, &iPointIndices);
    m_plotHighLightedPoints->SetInputData(m_tableHighLightedPoints, "Y1", "Y2", "Density");
  }*/ 
}

void MaterialParametersView::onLevelVisibilityModified()
{
  updateChart();
}

void MaterialParametersView::onLevelsModified()
{
  updatePlots();
}

