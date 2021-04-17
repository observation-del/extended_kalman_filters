#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

/* 
 * Please note that the Eigen library does not initialize 
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
   * TODO: predict the state
   */
  x_ = F_*x_;
  P_ = F_ * P_ * F_.transpose();

}

void KalmanFilter::Update(const VectorXd &z) {
  /**
   * TODO: update the state by using Kalman Filter equations
   */
  static Eigen::MatrixXd I;
  I = MatrixXd(4,4);
  I << 1,0,0,0,
    0,1,0,0,
    0,0,1,0,
    0,0,0,1;

  Eigen::VectorXd y_laser;
  y_laser = VectorXd(4);
  y_laser = z - H_*x_;

  Eigen::MatrixXd S_laser;
  S_laser = MatrixXd(4, 4);
  S_laser = H_ * P_ * H_.transpose() + R_;
  
  Eigen::MatrixXd K_laser;
  K_laser = MatrixXd(4,4);
  K_laser = P_ * H_.transpose() * S_laser.inverse();

  x_ = x_ + K_laser*y_laser;
  
  P_ = (I - K_laser*H_)*P_;


}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
   * TODO: update the state by using Extended Kalman Filter equations
   */

  static Eigen::MatrixXd I;
  I = MatrixXd(4,4);
  I << 1,0,0,0,
    0,1,0,0,
    0,0,1,0,
    0,0,0,1;

  Eigen::VectorXd y_radar;
  y_radar = VectorXd(3);
  y_radar = z - H_*x_;

  Eigen::MatrixXd S_radar;
  S_radar = MatrixXd(3, 3);
  S_radar = H_ * P_ * H_.transpose() + R_;
  
  Eigen::MatrixXd K_radar;
  K_radar = MatrixXd(4,3);
  K_radar = P_ * H_.transpose() * S_radar.inverse();

  x_ = x_ + K_radar*y_radar;
  
  P_ = (I - K_radar*H_)*P_;

}
