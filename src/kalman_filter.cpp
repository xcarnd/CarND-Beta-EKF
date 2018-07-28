#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

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
  TODO:
    * predict the state
  */
  x_ = F_ * x_;
  P_ = F_ * P_ * F_.transpose() + Q_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  VectorXd z_p(3);
  const double rho = sqrt(pow(x_(0), 2) + pow(x_(1), 2));
  double phi = atan2(x_(1), x_(0));
  z_p << rho,
       phi,
       (x_(0) * x_(2) + x_(1) * x_(3)) / rho;
  while (z_p(1) < M_PI) z_p(1) += 2 * M_PI;
  while (z_p(1) > M_PI) z_p(1) -= 2 * M_PI;
  
  VectorXd y = z - z_p;
  MatrixXd S = H_ * P_ * H_.transpose() + R_;
  MatrixXd K = P_ * H_.transpose() * S.inverse();
  
  x_ = x_ + K * y;
  P_ = P_ - K * H_ * P_;
}

void KalmanFilter::Update(const VectorXd &z) {
  VectorXd y = z - H_ * x_;
  MatrixXd S = H_ * P_ * H_.transpose() + R_;
  MatrixXd K = P_ * H_.transpose() * S.inverse();
  x_ = x_ + K * y;
  P_ = P_ - K * H_ * P_;
}
