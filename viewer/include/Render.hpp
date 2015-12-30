/*
 * Render.hpp
 *
 *  Created on: Aug 20, 2013
 *      Author: desaic
 */

#ifndef RENDER_HPP_
#define RENDER_HPP_
#include <vector>
#include "Camera.hpp"
#include <Eigen/Dense>
struct GLFWwindow;
class ElementHex;
class World;
class Element;
class ElementMesh;
class Element2D;
class ElementMesh2D;
class Stepper;
class Render
{
public:
  Render();
  void init(World * world);
  int loop();
  void draw();
  void drawEle(int eidx, ElementMesh * eMesh);
  void drawEleMesh(ElementMesh * eMesh);
  void drawEle2D(int eidx, ElementMesh2D * eMesh);
  void drawEleMesh2D(ElementMesh2D * eMesh);
  void moveCamera(float dt);
  Stepper * getStepper();
  virtual ~Render();
  bool anim;
  Camera camera;
  GLFWwindow* window;
  
  std::vector<Eigen::Vector3f> matColor;

  //how fast to rotate in x and y axis
  float xRotSpeed, yRotSpeed;
  float camSpeed;
  void toggleForce();
  int forceIdx;
private:
  ///@brief Render does not own this pointer.
  World * world;
};

#endif /* RENDER_HPP_ */
