#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
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
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */

  ekf_.P_ = MatrixXd(4, 4);
  ekf_.P_ << 1000, 0, 0, 0,
          0, 1000, 0, 0,
          0, 0, 1000, 0,
          0, 0, 0, 1000;
 

}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    cout<< "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      const VectorXd& m = measurement_pack.raw_measurements_;
      double rho = m[0];
      double phi = m[1];
      double px = rho / (sqrt(1 + pow(tan(phi), 2)));
      double py = sqrt(pow(rho, 2) - pow(px, 2));
      ekf_.x_(0) = px;
      ekf_.x_(1) = py;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      const VectorXd& m = measurement_pack.raw_measurements_;
      ekf_.x_(0) = m(0);
      ekf_.x_(1) = m(1);
    }

    // done initializing, no need to predict or update
    previous_timestamp_ = measurement_pack.timestamp_;
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
  long long dt_ll = measurement_pack.timestamp_ - previous_timestamp_;
  double dt = dt_ll / 1e6;
  ekf_.F_ = MatrixXd(4, 4);
  ekf_.F_.setIdentity();
  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;
  
  ekf_.Q_ = MatrixXd(4, 4);
  const double noise_ax = 9;
  const double noise_ay = 9; 
  const double dt2 = dt * dt;
  const double dt3_h = dt2 * dt / 2.0;
  const double dt4_q = dt2 * dt2 / 4.0;
  ekf_.Q_ << dt4_q * noise_ax, 0, dt3_h * noise_ax, 0,
         0, dt4_q * noise_ay, 0, dt3_h * noise_ay,
        dt3_h * noise_ax, 0, dt2 * noise_ax, 0,
        0, dt3_h * noise_ay, 0, dt2 * noise_ay;

  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
    ekf_.R_ = R_radar_;
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else {
    // Laser updates
    ekf_.H_ = MatrixXd(2, 4);
    ekf_.H_ << 1, 0, 0, 0,
           0, 1, 0, 0;
  
    ekf_.R_ = R_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
  previous_timestamp_ = measurement_pack.timestamp_;
}
