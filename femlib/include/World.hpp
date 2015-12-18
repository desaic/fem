#ifndef WORLDHPP
#define WORLDHPP
#include <vector>
class ElementMesh;
class ElementMesh2D;
class Stepper;
class World{
public:
  World() :stepper(0){}
  std::vector<ElementMesh*> em;
  std::vector<ElementMesh2D*> em2d;
  Stepper * stepper;
};
#endif