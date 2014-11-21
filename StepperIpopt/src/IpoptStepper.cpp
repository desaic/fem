#include "Stepper.hpp"
#include "IpoptStepper.hpp"
#include "IpoptInterface.hpp"
#include "ElementMesh.hpp"
#include "IpIpoptApplication.hpp"
#include <iostream>

class ElementMesh;
using namespace Ipopt;

void IpoptStepper::step(ElementMesh * mesh)
{
  SmartPtr<IpoptInterface> problem = new IpoptInterface();
  
  SmartPtr<IpoptApplication> app = IpoptApplicationFactory();
  app->Options()->SetNumericValue("tol", 0.5);
  app->Options()->SetStringValue("mu_strategy", "adaptive");
//  app->Options()->SetStringValue("hessian_approximation","limited-memory");
  app->Options()->SetIntegerValue("max_iter", NSteps);
//steps each mesh separately

 //   char outfile[32];
//    snprintf(outfile,32,"out%lu",ii);
 //   app->Options()->SetStringValue("output_file", outfile);
    problem->ele = mesh;
    ApplicationReturnStatus status;
    status = app->Initialize();
    if (status != Solve_Succeeded) {
      std::cout << std::endl << std::endl << "*** Error during initialization!" << std::endl;
      return;
    }
    // Ask Ipopt to solve the problem
    status = app->OptimizeTNLP(problem);

    if (status == Solve_Succeeded) {
      std::cout << std::endl << std::endl << ": *** The problem solved!" << std::endl;
    }
    else {
      std::cout << std::endl << std::endl << ":*** The problem FAILED!" << std::endl;
    }


}

IpoptStepper::IpoptStepper():NSteps(400)
{}

IpoptStepper::~IpoptStepper()
{}
