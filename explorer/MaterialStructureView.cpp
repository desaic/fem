
#include "MaterialStructureView.h"

#include <QKeyEvent>

#include <vtkContextView.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>

#include <vtkPoints.h>
#include <vtkTriangle.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkCellArray.h>
#include <vtkPolygon.h>
#include <vtkProperty.h>
#include <vtkCellData.h>

#include <ElementMesh.hpp>
#include <Element.hpp>
#include <array>

#include "exProject.h"

#include "ExplorerUtilities.h"
using namespace ExplorerUtilities;

#include "cfgMaterialUtilities.h"
using namespace cfgMaterialUtilities;

MaterialStructureView::MaterialStructureView()
  :QVTKWidget()
{
}

MaterialStructureView::~MaterialStructureView()
{
}

void MaterialStructureView::setProject(QPointer<exProject> iProject)
{
  m_project = iProject;

  connect(m_project, SIGNAL(pickedStructureModified()), this, SLOT(onPickedStructureModified()));

  Q_ASSERT(m_project);

  vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
  vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
  renderWindow->AddRenderer(renderer);
  SetRenderWindow(renderWindow);

  renderer->SetBackground(1, 1, 1);
  m_renderer = renderer;

  renderWindow->Render();
  m_renderWindow = renderWindow;
}

void MaterialStructureView::onPickedStructureModified()
{
  int pickedStructureIndex = m_project->getPickedStructureIndex();
  cout << "MaterialStructureView::onPickedStructureModified()   pickedStructureIndex = " << pickedStructureIndex << std::endl;  
  if (pickedStructureIndex>=0)
  {
    int pickedStructureLevel = m_project->getPickedStructureLevel();
    cout << "MaterialStructureView::onPickedStructureModified()   updateMesh" << std::endl;  
    updateMesh(pickedStructureIndex, pickedStructureLevel);
  }
  update();
}

void MaterialStructureView::updateMesh(int iCombIndex, int iLevel)
{
  //QSharedPointer<ElementMesh> elementMesh = m_project->computeElementMesh(iCombIndex, iLevel);
  QSharedPointer<ElementMesh> elementMesh = m_project->computeElementMeshIncr(iCombIndex, iLevel);
  vtkSmartPointer<vtkPolyData> polyData = createVtkPolyData(elementMesh.data());

  vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  mapper->SetInputData(polyData);

  if (m_actor)
  {
    m_renderer->RemoveActor(m_actor);
  }

  m_actor = vtkSmartPointer<vtkActor>::New();
  m_actor->SetMapper(mapper);
  m_renderer->AddActor(m_actor);
  if (m_renderWindow)
  {
    m_renderWindow->Render();
  }
}

vtkSmartPointer<vtkPolyData> MaterialStructureView::createVtkPolyData(const ElementMesh * iElementMesh)
{
  assert(iElementMesh);

  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray> faces = vtkSmartPointer<vtkCellArray>::New();

  unsigned char colMat[2][4] =  {{255, 255, 255, 120}, {0, 0, 255, 255}};

  vtkSmartPointer<vtkUnsignedCharArray> colors = vtkSmartPointer<vtkUnsignedCharArray>::New();
  colors->SetNumberOfComponents(4);
  colors->SetName("Colors");

  int externalFaceIndices[3][2][4] = {
    { {0,1,3,2}, {5,4,6,7} },
    { {1,0,4,5}, {2,3,7,6} }, 
    { {4,0,2,6}, {1,5,7,3} } };

  int indPoint = 0;
  int ielement=0, nelement=(int)iElementMesh->e.size();
  for(ielement=0; ielement<nelement; ielement++)
  {
    Element * ele = iElementMesh->e[ielement];
    const std::vector<int> & nodes = ele->getNodeIndices();
    int indMat = iElementMesh->me[ielement];
    if (indMat==1)
    {
      for (int iaxis=0; iaxis<3; iaxis++)
      {
        for (int iside=0; iside<2; iside++)
        {
          vtkSmartPointer<vtkPolygon> polygon = vtkSmartPointer<vtkPolygon>::New();
          polygon->GetPointIds()->SetNumberOfIds(4); //make a quad
          for (int ivertex=0; ivertex<4; ivertex++)
          {
            int ind = externalFaceIndices[iaxis][iside][ivertex];
            int indNode = nodes[ind];

            Vector3f v = iElementMesh->X[indNode];
            points->InsertNextPoint (v[0], v[1], v[2]);
            polygon->GetPointIds()->SetId(ivertex, indPoint++);
          }
          faces->InsertNextCell(polygon);
          colors->InsertNextTupleValue(colMat[indMat]);
        }
      }
    }

    /*std::vector<std::array<int,2> > edges = ele->getEdges();
    int indPoint=0;
    for(unsigned int ii = 0;ii<edges.size();ii++)
    {
      int indV1 = (*ele)[edges[ii][0]];
      int indV2 = (*ele)[edges[ii][1]];

      Vector3f v1 = iElementMesh->X[indV1];
      Vector3f v2 = iElementMesh->X[indV2];

      points->InsertNextPoint (v1[0], v1[1], v1[2]);
      points->InsertNextPoint (v2[0], v2[1], v2[2]);
      points->InsertNextPoint (v1[0], v1[1], v1[2]);

      vtkSmartPointer<vtkTriangle> triangle = vtkSmartPointer<vtkTriangle>::New();
      triangle->GetPointIds()->SetId(0, indPoint++);
      triangle->GetPointIds()->SetId(1, indPoint++);
      triangle->GetPointIds()->SetId(2, indPoint++);

      triangles->InsertNextCell(triangle);
    }*/ 
  } 

  vtkSmartPointer<vtkPolyData> polyData =  vtkSmartPointer<vtkPolyData>::New();
  polyData->SetPoints(points);
  polyData->SetPolys(faces); 
  polyData->GetCellData()->SetScalars(colors);

  return polyData;
}

void MaterialStructureView::keyPressEvent(QKeyEvent *pEvent)
{
  switch (pEvent->key()) 
  {
  case Qt::Key_S:
    saveStructure();
    break;
  default:
    QVTKWidget::keyPressEvent(pEvent);
  }
}

void MaterialStructureView::saveStructure()
{
  int pickedStructureIndex = m_project->getPickedStructureIndex();
  if (pickedStructureIndex>=0)
  {
    int pickedStructureLevel = m_project->getPickedStructureLevel();
    const std::vector<std::vector<std::vector<int> > > & materialAssignments = m_project->getMaterialAssignments();
    const std::vector<int> & matAssignment = materialAssignments[pickedStructureLevel][pickedStructureIndex];
    const std::vector<int> & levels = m_project->getLevels();

    int n = levels[pickedStructureLevel];
    std::string fileName = m_project->getFileDirectory() + "MaterialAssignment_" + std::to_string(n) + "_" + std::to_string(pickedStructureIndex) + ".txt";

    saveMicrostructure(fileName, n, n, matAssignment);
  }
}

