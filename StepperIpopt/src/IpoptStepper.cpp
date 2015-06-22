
#include <iostream>
#include "IpoptStepper.hpp"
class ElementMesh;
using namespace Ipopt;

void IpoptStepper::init(ElementMesh * _m)
{
  std::cout<<"Init ipopt solver\n";
  m = _m;
  problem = new IpoptInterface();

  app = IpoptApplicationFactory();
  app->Options()->SetNumericValue("tol", 0.5);
  app->Options()->SetStringValue("mu_strategy", "adaptive");
//  app->Options()->SetStringValue("hessian_approximation","limited-memory");
  app->Options()->SetIntegerValue("max_iter", NSteps);
//steps each mesh separately

 //   char outfile[32];
//    snprintf(outfile,32,"out%lu",ii);
 //   app->Options()->SetStringValue("output_file", outfile);
    problem->ele = m;
    ApplicationReturnStatus status;
    status = app->Initialize();
    if (status != Solve_Succeeded) {
      std::cout << std::endl << std::endl << "*** Error during initialization!" << std::endl;
      return;
    }

    std::cout<<"end init ipoptstepper\n";
}

int IpoptStepper::oneStep()
{
    // Ask Ipopt to solve the problem
    ApplicationReturnStatus status = app->OptimizeTNLP(problem);

    //just need to call once for statics
    return -1;
    if (status == Solve_Succeeded) {
      std::cout << std::endl << std::endl << ": *** The problem solved!" << std::endl;
      return 0;
    }
    else {
      std::cout << std::endl << std::endl << ":*** The problem FAILED!" << std::endl;
      return -1;
    }
}

IpoptStepper::IpoptStepper():NSteps(400)
{}

IpoptStepper::~IpoptStepper()
{}
