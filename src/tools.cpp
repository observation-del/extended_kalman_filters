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
  if (estimations.size() != ground_truth.size()
      || estimations.size() == 0) {
    std::cout << "Invalid estimation or ground_truth data" << std::endl;
    return rmse;
  }

  // accumulate squared residuals
  for (unsigned int i=0; i < estimations.size(); ++i) {

    VectorXd residual = estimations[i] - ground_truth[i];

    // coefficient-wise multiplication
    residual = residual.array()*residual.array();
    rmse += residual;
  }

  // calculate the mean
  rmse = rmse/estimations.size();

  // calculate the squared root
  rmse = rmse.array().sqrt();

  // return the result
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
