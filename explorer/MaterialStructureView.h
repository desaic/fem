#ifndef MaterialStructureView_h
#define MaterialStructureView_h

#include <QVTKWidget.h> 
#include <QPointer>

#include <vtkSmartPointer.h>

class vtkPolyData;
class vtkPolyDataMapper;
class ElementMesh;
class exProject;

class vtkActor;
class vtkRenderer;
class vtkRenderWindow;

class MaterialStructureView: public QVTKWidget
{
  Q_OBJECT
public:
  MaterialStructureView();
  ~MaterialStructureView();

  void setProject(QPointer<exProject> iProject);

private:
  vtkSmartPointer<vtkPolyData> createVtkPolyData(const ElementMesh * iElementMesh);
  void updateMesh(int iCombIndex, int iLevel);

  private slots:
    void onPickedStructureModified();

private:
  QPointer<exProject> m_project;

  vtkSmartPointer<vtkPolyDataMapper> m_mapper;
  vtkSmartPointer<vtkActor> m_actor;
  vtkSmartPointer<vtkRenderer> m_renderer;
  vtkSmartPointer<vtkRenderWindow> m_renderWindow;
};

#endif 

