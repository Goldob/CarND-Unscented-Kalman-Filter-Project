#include "ukf.h"
#include "Eigen/Dense"
#include <cmath>

#define EPS 1e-8

#define NORMALIZE_ANGLE(a) (a = std::fmod(a + M_PI, 2*M_PI) - M_PI)

using std::cos;
using std::sin;
using std::abs;
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
  std_a_ = 30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 30;

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
  float delta_t =  (time_us_ - meas_package.timestamp_) * 1e-6;
  time_us_ = meas_package.timestamp_;

  /**
   * 1. Initialization
   */
  if (!is_initialized_) {
    x_.fill(0.0);
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      /**
       * Convert radar from polar to cartesian coordinates and initialize state
       */

      float rho = meas_package.raw_measurements_[0];
      float phi = meas_package.raw_measurements_[1];

      x_(0) = rho * cos(phi);
      x_(1) = rho * sin(phi);
    } else {
      /**
       * Initialize state
       */

      x_(0) = meas_package.raw_measurements_[0];
      x_(1) = meas_package.raw_measurements_[1];
    }

    /**
     * Initialize  state covariance matrix
     */
    P_.fill(0.0);
    P_(0, 0) = 10;
    P_(1, 1) = 10;
    P_(2, 2) = 100;
    P_(3, 3) = 100;
    P_(4 ,4) = 100;

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
  MatrixXd Xsig_pred = MatrixXd::Zero(n_x_, n_sig_);
  for (int i = 0; i < n_sig_; i++) {
    VectorXd sig_point = Xsig_aug.col(i);
    VectorXd sig_point_pred = sig_point.head(5);

    double v, phi, phi_d, a, phi_dd;

    v      = sig_point(2);
    phi    = sig_point(3);
    phi_d  = sig_point(4);
    a      = sig_point(5);
    phi_dd = sig_point(6);

    if (abs(phi_d) < EPS) {
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

    NORMALIZE_ANGLE(sig_point_pred(3));
    NORMALIZE_ANGLE(sig_point_pred(4));

    Xsig_pred.col(i) = sig_point_pred;
  }

  /**
   * Calculate predicted mean
   */
  x_.fill(0.0);
  for (int i = 0; i < n_sig_; i++) {
    x_ += weights_(i) * Xsig_pred.col(i);
  }

  NORMALIZE_ANGLE(x_(3));
  NORMALIZE_ANGLE(x_(4));

  /**
   * Calculate predicted state covariance
   */
  P_.fill(0.0);
  for (int i = 0; i < n_sig_; i++) {
    VectorXd diff = Xsig_pred.col(i) - x_;
    P_ += weights_(i) * diff * diff.transpose();
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
  /**
   * TODO: Complete this function! Use radar data to update the belief
   * about the object's position. Modify the state vector, x_, and
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */
}
