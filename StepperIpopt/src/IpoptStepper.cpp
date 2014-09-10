#include "TimeStepper/IpoptStepper.hpp"
#include "TimeStepper/IpoptInterface.hpp"
#include "World/World.hpp"
#include "IpIpoptApplication.hpp"
#include <iostream>
class ElementMesh;
using namespace Ipopt;
void IpoptStepper::Step(World * world)
{
  SmartPtr<IpoptInterface> problem = new IpoptInterface();
  
  SmartPtr<IpoptApplication> app = IpoptApplicationFactory();
  app->Options()->SetNumericValue("tol", 0.5);
  app->Options()->SetStringValue("mu_strategy", "adaptive");
//  app->Options()->SetStringValue("hessian_approximation","limited-memory");
  app->Options()->SetIntegerValue("max_iter", NSteps);
//steps each mesh separately
  for(size_t ii =0 ;ii<world->element.size();ii++){
 //   char outfile[32];
//    snprintf(outfile,32,"out%lu",ii);
 //   app->Options()->SetStringValue("output_file", outfile);
    ElementMesh * ele = world->element[ii];
    problem->ele = ele;
    ApplicationReturnStatus status;
    status = app->Initialize();
    if (status != Solve_Succeeded) {
      std::cout << std::endl << std::endl << "*** Error during initialization!" << std::endl;
      return;
    }
    // Ask Ipopt to solve the problem
    status = app->OptimizeTNLP(problem);

    if (status == Solve_Succeeded) {
      std::cout << std::endl << std::endl << ii << ": *** The problem solved!" << std::endl;
    }
    else {
      std::cout << std::endl << std::endl << ii << ":*** The problem FAILED!" << std::endl;
    }
  }

}

IpoptStepper::IpoptStepper():NSteps(400)
{}

IpoptStepper::~IpoptStepper()
{}
