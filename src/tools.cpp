#include "tools.h"
#include <iostream>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
   * TODO: Calculate the RMSE here.
   */
  VectorXd rmse(4);
  rmse << 0,0,0,0;

  // check the validity of the following inputs:
  //  * the estimation vector size should not be zero
  //  * the estimation vector size should equal ground truth vector size
  if (estimations.size() != ground_truth.size()){
      std::cout << "input error" << std::endl;
      return rmse;
  }

  // accumulate squared residuals
  for (unsigned int i=0; i < estimations.size(); ++i) {
    // ... your code here
    rmse(0) += pow((estimations[i](0) - ground_truth[i](0)),2);
    rmse(1) += pow((estimations[i](1) - ground_truth[i](1)),2);
    rmse(2) += pow((estimations[i](2) - ground_truth[i](2)),2);
    rmse(3) += pow((estimations[i](3) - ground_truth[i](3)),2);
    
  }

  // calculate the mean
  rmse(0) = rmse(0)/estimations.size();
  rmse(1) = rmse(1)/estimations.size();
  rmse(2) = rmse(2)/estimations.size();
  rmse(3) = rmse(3)/estimations.size();
  

  // calculate the squared root
  rmse(0) = pow(rmse(0),1.0/2.0); 
  rmse(1) = pow(rmse(1),1.0/2.0); 
  rmse(2) = pow(rmse(2),1.0/2.0); 
  rmse(3) = pow(rmse(3),1.0/2.0); 

  return rmse;

}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
   * TODO:
   * Calculate a Jacobian here.
   */
   Eigen::MatrixXd Hj;
   Hj = MatrixXd(3,4);

   float c1 = x_state(0)*x_state(0) + x_state(1)*x_state(1);
   float c2 = sqrt(c1);
   float c3 = (c1*c2);
   Hj << (x_state(0)/c2), (x_state(1)/c2), 0, 0,
      -(x_state(1)/c1), (x_state(0)/c1), 0, 0,
      x_state(1)*(x_state(2)*x_state(1) - x_state(3)*x_state(0))/c3, x_state(0)*(x_state(3)*x_state(0) - x_state(2)*x_state(1))/c3, x_state(0)/c2, x_state(1)/c2;

   return Hj;
}
