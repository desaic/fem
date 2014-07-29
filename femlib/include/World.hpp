#ifndef WORLDHPP
#define WORLDHPP
#include <vector>
class ElementMesh;
class World{
public:
  std::vector<ElementMesh*> em;
};
#endif