#include "FusionEKF.h"
#include <iostream>
#include "Eigen/Dense"
#include "tools.h"
#include <cstdio>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;
using std::vector;

/**
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;

  /**
   * TODO: Finish initializing the FusionEKF.
   * TODO: Set the process and measurement noises
   */

  // the initial transition matrix F_
  ekf_.F_ = MatrixXd(4, 4);
  ekf_.F_ << 1, 0, 1, 0,
            0, 1, 0, 1,
            0, 0, 1, 0,
            0, 0, 0, 1;

  // set the process covariance matrix Q
  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.Q_ <<  0, 0, 0, 0,
         0, 0, 0, 0,
         0, 0, 0, 0,
         0, 0, 0, 0;

  // set the measurement matrix H for laser
  H_laser_ << 1, 0, 0, 0,
          0, 1, 0, 0;
  
  // set the measurement matrix H
  Hj_ << 1, 1, 0, 0,
      1, 1, 0, 0,
      1, 1, 1, 1;
  

}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  /**
   * Initialization
   */
  printf("FusionEKF.cpp line: %d\n", __LINE__);
  if (!is_initialized_) {
    /**
     * TODO: Initialize the state ekf_.x_ with the first measurement.
     * TODO: Create the covariance matrix.
     * You'll need to convert radar from polar to cartesian coordinates.
     */

    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    // state covariance matrix P
    ekf_.P_ = MatrixXd(4, 4);
    ekf_.P_ << 1, 0, 0, 0,
            0, 1, 0, 0,
            0, 0, 1000, 0,
            0, 0, 0, 1000;

    Eigen::MatrixXd I;
    I = MatrixXd(4,4);
    I << 1,0,0,0,
      0,1,0,0,
      0,0,1,0,
      0,0,0,1;

    printf("FusionEKF.cpp line: %d\n", __LINE__);

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      // TODO: Convert radar from polar to cartesian coordinates 
      //         and initialize state.
      // pre-compute a set of terms to avoid repeated calculation

      float c1 = ekf_.x_(0)*ekf_.x_(0) + ekf_.x_(1)*ekf_.x_(1);
      
      // check division by zero
      if (fabs(c1) < 0.0001) {
        cout << "CalculateJacobian () - Error - Division by Zero" << endl;
        return;
      }

      else{
        printf("FusionEKF.cpp line: %d\n", __LINE__);
        // compute the Jacobian matrix
        Hj_ = tools.CalculateJacobian(ekf_.x_); 

        Eigen::VectorXd y_radar;
        y_radar = VectorXd(3);
        y_radar = measurement_pack.raw_measurements_ - Hj_*ekf_.x_;

        Eigen::MatrixXd S_radar;
        S_radar = MatrixXd(3, 3);
        S_radar = Hj_ * ekf_.P_ * Hj_.transpose() + R_radar_;
        
        Eigen::MatrixXd K_radar;
        K_radar = MatrixXd(4,3);
        K_radar = ekf_.P_ * Hj_.transpose() * S_radar.inverse();

        ekf_.x_ = ekf_.x_ + K_radar*y_radar;
        
        ekf_.P_ = (I - K_radar*Hj_)*ekf_.P_;

      }
      printf("FusionEKF.cpp line: %d\n", __LINE__);
    }

    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      printf("FusionEKF.cpp line: %d\n", __LINE__);
      // TODO: Initialize state.
      Eigen::VectorXd y_laser;
      y_laser = VectorXd(4);
      y_laser = measurement_pack.raw_measurements_ - H_laser_*ekf_.x_;
      printf("FusionEKF.cpp line: %d\n", __LINE__);
      Eigen::MatrixXd S_laser;
      S_laser = MatrixXd(4, 4);
      S_laser = H_laser_ * ekf_.P_ * H_laser_.transpose() + R_laser_;
      printf("FusionEKF.cpp line: %d\n", __LINE__);
      Eigen::MatrixXd K_laser;
      K_laser = MatrixXd(4,4);
      K_laser = ekf_.P_ * H_laser_.transpose() * S_laser.inverse();
      printf("FusionEKF.cpp line: %d\n", __LINE__);
      ekf_.x_ = ekf_.x_ + K_laser*y_laser;
      printf("FusionEKF.cpp line: %d\n", __LINE__);
      ekf_.P_ = (I - K_laser*H_laser_)*ekf_.P_;
      printf("FusionEKF.cpp line: %d\n", __LINE__);
    }

    previous_timestamp_ = measurement_pack.timestamp_;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /**
   * Prediction
   */

  /**
   * TODO: Update the state transition matrix F according to the new elapsed time.
   * Time is measured in seconds.
   * TODO: Update the process noise covariance matrix.
   * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */

  long long dt = measurement_pack.timestamp_ - previous_timestamp_;
  long long dt2 = dt*dt;
  long long dt3 = dt*dt2;
  long long dt4 = dt2*dt2;
  long long noise_ax = 9;
  long long noise_ay = 9; 
  long long noise_ax2 = noise_ax*noise_ax;
  long long noise_ay2 = noise_ay*noise_ay;

  ekf_.F_ << 1, 0, measurement_pack.timestamp_ - previous_timestamp_, 0,
        0, 1, 0, measurement_pack.timestamp_ - previous_timestamp_,
        0, 0, 1, 0,
        0, 0, 0, 1;

  
  ekf_.Q_ << dt4*noise_ax2/4, 0, dt3*noise_ax2/2, 0,
            0, dt4*noise_ay2/4, 0, dt3*noise_ay2/2,
            dt3*noise_ax2/2, 0, dt2*noise_ax2, 0,
            0, dt3*noise_ay2/2, 0, dt2*noise_ay2;

  ekf_.Predict();

  /**
   * Update
   */

  /**
   * TODO:
   * - Use the sensor type to perform the update step.
   * - Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // TODO: Radar updates

    float c1 = ekf_.x_(0)*ekf_.x_(0) + ekf_.x_(1)*ekf_.x_(1);

    // check division by zero
    if (fabs(c1) < 0.0001) {
      cout << "CalculateJacobian () - Error - Division by Zero" << endl;
      return;
    }

    // compute the Jacobian matrix
    Hj_ = tools.CalculateJacobian(ekf_.x_); 

    ekf_.H_ = Hj_;
    ekf_.R_ = R_radar_;

    ekf_.UpdateEKF(measurement_pack.raw_measurements_);


  } else {
    // TODO: Laser updates
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;

    ekf_.Update(measurement_pack.raw_measurements_);

  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
