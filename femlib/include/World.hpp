#ifndef WORLDHPP
#define WORLDHPP
#include <vector>
class ElementMesh;
class ElementMesh2D;
class Stepper;
class World{
public:
  World() :stepper(0), u(0), fe(0){}
  std::vector<ElementMesh*> em;
  std::vector<ElementMesh2D*> em2d;
  Stepper * stepper;
  //displacements for em[0]. Initialized to NULL.
  std::vector< std::vector<double> > * u;
  std::vector< std::vector<double> > * fe;
};
#endif