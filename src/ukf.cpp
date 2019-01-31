#include "ukf.h"
#include "Eigen/Dense"
#include <cmath>

#define EPS 1e-8

using std::cos;
using std::sin;
using std::pow;
using std::sqrt;

using Eigen::MatrixXd;
using Eigen::VectorXd;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1.0;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.5;

  /**
   * DO NOT MODIFY measurement noise values below.
   * These are provided by the sensor manufacturer.
   */

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  /**
   * End DO NOT MODIFY section for measurement noise values
   */

  /**
   * TODO: Complete the initialization. See ukf.h for other member properties.
   * Hint: one or more values initialized above might be wildly off...
   */

  // State dimension
  n_x_ = 5;

  // Augmented state dimension
  n_aug_ = 7;

  // Sigma point spreading parameter
  lambda_ = 3 - n_aug_;

  // Number of sigma points
  n_sig_ = 2 * n_aug_ + 1;

  // Sigma points weights
  weights_ = VectorXd(n_sig_);
  weights_(0) = lambda_ / (lambda_ + n_aug_);
  for (int i = 1; i < n_sig_; i++) {
    weights_(i) = 1 / (lambda_ + n_aug_) * 0.5;
  }
}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  float delta_t =  (meas_package.timestamp_ - time_us_) * 1e-6;
  time_us_ = meas_package.timestamp_;

  /**
   * 1. Initialization
   */
  if (!is_initialized_) {
    x_.fill(0.0);
    P_.fill(0.0);
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      /**
       * Convert radar from polar to cartesian coordinates and initialize state
       */

      float rho = meas_package.raw_measurements_[0];
      float phi = meas_package.raw_measurements_[1];

      x_(0) = rho * cos(phi);
      x_(1) = rho * sin(phi);

      P_(0, 0) = pow(std_radr_ * cos(phi), 2)
               + pow(rho       * sin(phi), 2);

      P_(1, 1) = pow(std_radr_ * sin(phi), 2)
               + pow(rho       * cos(phi), 2);
    } else {
      /**
       * Initialize state
       */

      x_(0) = meas_package.raw_measurements_[0];
      x_(1) = meas_package.raw_measurements_[1];

      P_(0, 0) = pow(std_laspx_, 2);
      P_(1, 1) = pow(std_laspy_, 2);
    }

    /**
     * Initialize  state covariance matrix
     */
    P_(2, 2) = 0.5;
    P_(3, 3) = 0.5;
    P_(4 ,4) = 0.5;

    is_initialized_ = true;
    return;
  }

  /**
   * 2. Prediction
   */
  Prediction(delta_t);

  /**
   * 3. Update
   */
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    UpdateRadar(meas_package);
  } else {
    UpdateLidar(meas_package);
  }
}

void UKF::Prediction(double delta_t) {

  /**
   * Create augmented mean vector
   */
  VectorXd x_aug = VectorXd::Zero(7);
  x_aug.head(5) = x_;

  /**
   * Create augmented state covariance
   */
  MatrixXd P_aug = MatrixXd::Zero(7, 7);
  P_aug.topLeftCorner(5, 5) = P_;
  P_aug(5, 5) = pow(std_a_, 2);
  P_aug(6, 6) = pow(std_yawdd_, 2);

  /**
   * Create square root matrix
   */
  MatrixXd A = P_aug.llt().matrixL();
  double coeff = sqrt(lambda_ + n_aug_);

  /**
   * Create sigma points
   */
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  Xsig_aug << VectorXd::Zero(7), coeff * A, -coeff * A;
  Xsig_aug.colwise() += x_aug;

  /**
   * Predict sigma points
   */
  Xsig_pred_ = MatrixXd::Zero(n_x_, n_sig_);
  for (int i = 0; i < n_sig_; i++) {
    VectorXd sig_point = Xsig_aug.col(i);
    VectorXd sig_point_pred = sig_point.head(5);

    double v, phi, phi_d, a, phi_dd;

    v      = sig_point(2);
    phi    = sig_point(3);
    phi_d  = sig_point(4);
    a      = sig_point(5);
    phi_dd = sig_point(6);

    if (fabs(phi_d) < EPS) {
      // avoid division by zero
      sig_point_pred(0) += v * cos(phi) * delta_t;
      sig_point_pred(1) += v * sin(phi) * delta_t;
      sig_point_pred(2) += 0;
      sig_point_pred(3) += phi * delta_t;
      sig_point_pred(4) += 0;
    } else {
      sig_point_pred(0) += v/phi_d * ( sin(phi + phi_d * delta_t) - sin(phi)),
      sig_point_pred(1) += v/phi_d * (-cos(phi + phi_d * delta_t) + cos(phi)),
      sig_point_pred(2) += 0;
      sig_point_pred(3) += phi_d * delta_t;
      sig_point_pred(4) += 0;
    }

    /**
     * Apply process noise
     */
    sig_point_pred(0) += 0.5 * pow(delta_t, 2) * cos(phi) * a;
    sig_point_pred(1) += 0.5 * pow(delta_t, 2) * sin(phi) * a;
    sig_point_pred(2) += delta_t * a;
    sig_point_pred(3) += 0.5 * pow(delta_t, 2) * phi_dd;
    sig_point_pred(4) += delta_t * phi_dd;

    Xsig_pred_.col(i) = sig_point_pred;
  }

  /**
   * Calculate predicted mean
   */
  x_.fill(0.0);
  for (int i = 0; i < n_sig_; i++) {
    x_ += weights_(i) * Xsig_pred_.col(i);
  }

  /**
   * Calculate predicted state covariance
   */
  P_.fill(0.0);
  for (int i = 0; i < n_sig_; i++) {
    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    // angle normalization
    while (x_diff(3) >  M_PI) x_diff(3) -= 2*M_PI;
    while (x_diff(3) < -M_PI) x_diff(3) += 2*M_PI;

    P_ += weights_(i) * x_diff * x_diff.transpose();
  }
}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
  VectorXd z(2);
  z << meas_package.raw_measurements_[0],
       meas_package.raw_measurements_[1];

  /**
   * Create measurement matrix
   */
  MatrixXd H = MatrixXd::Zero(2, 5);
  H(0, 0) = 1;
  H(1, 1) = 1;

  /*
   * Create measurement covatiance matrix
   */
  MatrixXd R = MatrixXd::Zero(2, 2);
  R(0, 0) = std_laspx_;
  R(1, 1) = std_laspy_;


  MatrixXd I = MatrixXd::Identity(5, 5);

  /*
   * Update state using regular Kalman Filter equations
   */
  VectorXd y = z - H * x_;
  MatrixXd S = H * P_ * H.transpose() + R;
  MatrixXd K = P_ * H.transpose() * S.inverse();

  x_ = x_ + K * y;
  P_ = (I - K * H) * P_;
}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  int n_z = 3;
  VectorXd z = VectorXd(n_z);
  z << meas_package.raw_measurements_[0],
       meas_package.raw_measurements_[1],
       meas_package.raw_measurements_[2];

  // transform sigma points into measurement space
  MatrixXd Zsig = MatrixXd(n_z, n_sig_);
  for (int i = 0; i < n_sig_; i++) {
    double p_x, p_y, v, phi;

    p_x   = Xsig_pred_(0, i);
    p_y   = Xsig_pred_(1, i);
    v     = Xsig_pred_(2, i);
    phi   = Xsig_pred_(3, i);

    Zsig(0, i) = sqrt(pow(p_x, 2) + pow(p_y, 2));
    Zsig(1, i) = atan2(p_y, p_x);
    Zsig(2, i) = (p_x*cos(phi) + p_y*sin(phi)) * v / Zsig(0, i);
  }

  // calculate mean predicted measurement
  VectorXd z_pred = VectorXd::Zero(n_z);
  for (int i = 0; i < n_sig_; i++) {
    z_pred += weights_(i) * Zsig.col(i);
  }

  // calculate innovation covariance matrix S
  MatrixXd S = MatrixXd::Zero(n_z, n_z);
  for (int i = 0; i < n_sig_; i++) {
    VectorXd z_diff = Zsig.col(i) - z_pred;

    // angle normalization
    while (z_diff(1)>  M_PI) z_diff(1) -= 2*M_PI;
    while (z_diff(1)< -M_PI) z_diff(1) += 2*M_PI;

    S += weights_(i) * z_diff * z_diff.transpose();
  }

  // add measurement noise to covariance matrix S
  S(0, 0) += pow(std_radr_,   2);
  S(1, 1) += pow(std_radphi_, 2);
  S(2, 2) += pow(std_radrd_,  2);

  // calculate cross correlation matrix
  MatrixXd Tc = MatrixXd::Zero(n_x_, n_z);
  for (int i = 0; i < n_sig_; i++) {
    // residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    // angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    // angle normalization
    while (x_diff(3) >  M_PI) x_diff(3) -= 2*M_PI;
    while (x_diff(3) < -M_PI) x_diff(3) += 2*M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  // calculate Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  // update state mean and covariance matrix
  VectorXd z_diff = z - z_pred;

  // angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  x_ += K * z_diff;
  P_ -= K * S * K.transpose();
}
