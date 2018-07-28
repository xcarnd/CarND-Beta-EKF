#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
  VectorXd sum_error(4);
  sum_error.fill(0.0);
  for (size_t i = 0; i < estimations.size(); ++i) {
    VectorXd error = estimations[i] - ground_truth[i];
    sum_error = sum_error.array() + error.array() * error.array();
  }
  sum_error /= estimations.size();
  return sum_error.array().sqrt();
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
  TODO:
    * Calculate a Jacobian here.
  */
  const double px = x_state(0);
  const double py = x_state(1);
  const double vx = x_state(2);
  const double vy = x_state(3);
  const double k2 = pow(px, 2) + pow(py, 2);
  const double k1 = sqrt(k2);
  const double k3 = k1 * k2;
  const double km = vx * py - vy * px;
  
  MatrixXd Hj(3, 4);
  Hj << px / k1, py / k1, 0, 0,
      -py / k2, px / k2, 0, 0,
      py * km / k3, px * km / k3, px / k1, py / k1;
  return Hj;
}
