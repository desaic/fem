#include "Render.hpp"
#include "World.hpp"
#include "ElementRegGrid.hpp"
int main(int argc, char* argv[])
{
  int nx = 2,ny=3,nz=4;
  ElementRegGrid * em = new ElementRegGrid(nx,ny,nz);
  World * world = new World();
  world->em.push_back(em);
  Render render;
  render.init(world);
  render.loop();
	return 0;
}

