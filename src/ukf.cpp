#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

#define MIN 0.001

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

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
  std_a_ = 3.0;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.5;

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

  /* my code starts below this line */

  // number of dimensions
  n_x_ = 5;
  // number of augmented dimensions
  n_aug_ = 7;
  //lambda value
  lambda_ = 3 - n_aug_;
  // avoid additional computation
  num_sigma_points = 2 * n_aug_ + 1;

  // set weights
  weights_ = VectorXd(num_sigma_points);
  weights_(0) = lambda_ / (lambda_ + n_aug_);

  float w_2_plus = 1 / (2 * (lambda_ + n_aug_));

  for(int j = 1; j < num_sigma_points; j++) {
    weights_(j) = w_2_plus;
  }

  // don't initialize until we process our first measurement
  is_initialized_ = false;
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {

  if(!is_initialized_) {
    // use our first measurement to initialize our kalman filter
    P_ << 1, 0, 0, 0, 0,
          0, 1, 0, 0, 0,
          0, 0, 1, 0, 0,
          0, 0, 0, 1, 0,
          0, 0, 0, 0, 1;

    previous_timestamp_ = meas_package.timestamp_;

    if(meas_package.sensor_type_ == MeasurementPackage::LASER) {
      x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0, 0;
    }

    if(meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      VectorXd cartesian = ConvertPolar(meas_package.raw_measurements_);
      x_ << cartesian[0], cartesian[1], meas_package.raw_measurements_[2], 0, 0;
    }

    is_initialized_ = true;

    return;
  }

  long current_timestamp = meas_package.timestamp_;
  float delta_t = (current_timestamp - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = meas_package.timestamp_;

  Prediction(delta_t);

  // Make our prediction and then based on measurement type update object representation
  if(use_laser_ && meas_package.sensor_type_ == MeasurementPackage::LASER) {
    UpdateLidar(meas_package);
  }

  if(use_radar_ && meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    UpdateRadar(meas_package);
  }
}


/**
 * Converts a state vector from polar coordinates to Cartesian coordinates
 * @param x_state the state vector containing polar coordinates
 */
VectorXd UKF::ConvertPolar(const VectorXd& x_state) {
  float rho     = x_state[0];
  float phi     = x_state[1];
  float rho_dot = x_state[2];

  float px      = rho * cos(phi);
  float py      = rho * sin(phi);
  float vx      = rho_dot * cos(phi);
  float vy      = rho_dot * sin(phi);

  VectorXd cartesian(4);
  cartesian << px,
               py,
               vx,
               vy;

  return cartesian;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  // generate sigma points and use them to make predictions

  // Generate Augmented Sigma Points

  // augmented mean
  VectorXd x_aug = VectorXd(n_aug_);
  x_aug.fill(0.0);

  // augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  P_aug.fill(0.0);

  // sigma point mtx
  MatrixXd Xsig_aug = MatrixXd(n_aug_, num_sigma_points);
  Xsig_aug.fill(0.0);

  // augmented mean
  x_aug << x_, 0, 0;

  // augmented covariance mtx
  P_aug.topLeftCorner(n_x_, n_x_) << P_;
  P_aug.bottomRightCorner(2, 2) << (std_a_ * std_a_), 0, 0, (std_yawdd_ * std_yawdd_);

  // create square root matrix
  MatrixXd A = P_aug.llt().matrixL();

  //create augmented sigma points
  double sqrt_lambda = sqrt(lambda_ + n_aug_);
  MatrixXd lambda_times_A = sqrt_lambda * A;

  MatrixXd x_aug_spread(n_aug_, n_aug_);
  x_aug_spread << x_aug, x_aug, x_aug, x_aug, x_aug, x_aug, x_aug;

  // compute our augmented sigma points
  MatrixXd first_term  = x_aug;
  MatrixXd second_term = x_aug_spread + lambda_times_A;
  MatrixXd third_term  = x_aug_spread - lambda_times_A;

  Xsig_aug << first_term, second_term, third_term;


  // Sigma Point Prediction

  // create matrix with predicted sigma points as columns
  MatrixXd Xsig_pred = MatrixXd(n_x_, num_sigma_points);

  //predict sigma points
  //avoid division by zero
  //write predicted sigma points into right column
  for(int i = 0; i < (num_sigma_points); i++) {
    float px              = Xsig_aug(0,i);
    float py              = Xsig_aug(1,i);
    float v               = Xsig_aug(2,i);
    float psi             = Xsig_aug(3,i);
    float psi_d           = Xsig_aug(4,i);
    float nu_a            = Xsig_aug(5,i);
    float nu_psi_dd       = Xsig_aug(6,i);

    float delta_t_squared = delta_t * delta_t;
    float sin_psi         = sin(psi);
    float cos_psi         = cos(psi);

    VectorXd pred(n_x_);

    if(fabs(psi_d) < MIN) {
      pred << px + v * cos_psi * delta_t + 0.5 * delta_t_squared * cos_psi * nu_a,
              py + v * sin_psi * delta_t + 0.5 * delta_t_squared * sin_psi * nu_a,
              v  + 0 + delta_t * nu_a,
              psi + psi_d  * delta_t + 0.5 * delta_t_squared * nu_psi_dd,
              psi_d  + 0 + delta_t * nu_psi_dd;

    } else {
      pred << px + (v / psi_d ) * ( sin(psi + psi_d * delta_t) - sin_psi) + 0.5 * delta_t_squared * cos_psi * nu_a,
              py + (v / psi_d ) * (-cos(psi + psi_d * delta_t) + cos_psi) + 0.5 * delta_t_squared * sin_psi * nu_a,
              v + 0 + delta_t * nu_a,
              psi + psi_d  * delta_t + 0.5 * delta_t_squared * nu_psi_dd,
              psi_d  + 0 + delta_t * nu_psi_dd;
    }
    Xsig_pred.col(i) = pred;
  }
  Xsig_pred_ = Xsig_pred;

  // predicted state
  VectorXd x = VectorXd(n_x_);
  x.fill(0.0);

  // covariance mtx for prediction
  MatrixXd P = MatrixXd(n_x_, n_x_);
  P.fill(0.0);

  // predict state mean
  for(int i = 0; i < num_sigma_points; i++) {
    x = x + weights_(i) * Xsig_pred_.col(i);
  }

  // predict state covariance mtx
  for(int i = 0; i < num_sigma_points; i++) {
    VectorXd x_diff = Xsig_pred_.col(i) - x;
    // atan2 handles angle normalization for us too.
    x_diff(3) = atan2(sin(x_diff(3)), cos(x_diff(3))) ;

    P = P + weights_(i) * x_diff * x_diff.transpose();
  }

  x_ = x;
  P_ = P;
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {

  //set measurement dimension, lidar can measure px, py
  int n_z = 2;

  // create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, num_sigma_points);
  Zsig.fill(0.0);

  // mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);

  // measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z, n_z);
  S.fill(0.0);

  // FIXME an interesting problem happens when
  // these for loops are combined.
  // The RMSE in position jumps up by .2
  // I leave them separate for now due to time but
  // I'm definitely curious what's going on

  // transform sigma points into measurement space
  for(int i = 0; i < num_sigma_points; i++) {
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);

    Zsig(0,i) = p_x;
    Zsig(1,i) = p_y;
  }

  // mean prediction
  for(int i = 0; i < num_sigma_points; i++) {
    z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  // measurement covariance mtx
  for(int i = 0; i < num_sigma_points; i++) {
    VectorXd z_diff = Zsig.col(i) - z_pred;
    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  S(0,0) += std_laspx_ * std_laspx_;
  S(1,1) += std_laspy_ * std_laspy_;

  // cross correlation mtx
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);

  for(int i = 0; i < num_sigma_points; i++) {
    VectorXd z_diff = Zsig.col(i) - z_pred;
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    // atan2 handles this for us
    x_diff(3) = atan2(sin(x_diff(3)), cos(x_diff(3))) ;
    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  // calculate Kalman gain K;
  MatrixXd kalman_gain = Tc * S.inverse();

  // set z to raw measurements
  VectorXd z = meas_package.raw_measurements_;
  VectorXd z_diff = z - z_pred;

  // update state mean and covariance matrix
  x_ = x_ + kalman_gain * z_diff;
  P_ = P_ - kalman_gain *  S * kalman_gain.transpose();
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  //set measurement dimension, radar can measure r, phi, and r_dot
  int n_z = 3;

  // TODO update Radar and Update Lidar do repeat
  // some small tasks that could be reduced
  // to make the code cleaner. ignoring for time.

  // matrix for sigma points
  MatrixXd Zsig = MatrixXd(n_z, num_sigma_points);
  Zsig.fill(0.0);

  // mean prediction
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);

  // measurement covariance mtx
  MatrixXd S = MatrixXd(n_z, n_z);
  S.fill(0.0);

  // transform sigma points into measurement space
  for(int i = 0; i < num_sigma_points; i++) {
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v = Xsig_pred_(2,i);
    double psi = Xsig_pred_(3,i);

    Zsig(0,i) = sqrt(p_x * p_x + p_y * p_y);
    Zsig(1,i) = atan2(p_y , p_x);
    Zsig(2,i) = (p_x * cos(psi) * v + p_y * sin(psi) * v) / Zsig(0,i);
  }

  // mean predicted measurement
  for(int i = 0; i < num_sigma_points; i++) {
    z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  // measurement covariance mtx
  for(int i = 0; i < num_sigma_points; i++) {
    VectorXd z_diff = Zsig.col(i) - z_pred;
    z_diff(1) = atan2(sin(z_diff(1)), cos(z_diff(1))) ;
    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  S(0,0) += std_radr_ * std_radr_;
  S(1,1) += std_radphi_ * std_radphi_;
  S(2,2) += std_radrd_ * std_radrd_;

  // TODO this could also be put into a separate function
  // ignoring because late submission

  // cross correlation mtx
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);

  for(int i = 0; i < num_sigma_points; i++) {
    VectorXd z_diff = Zsig.col(i) - z_pred;
    z_diff(1) = atan2(sin(z_diff(1)), cos(z_diff(1))) ;
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    x_diff(3) = atan2(sin(x_diff(3)), cos(x_diff(3))) ;
    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  // kalman gain;
  MatrixXd kalman_gain = Tc * S.inverse();

  // real measured values
  VectorXd z = meas_package.raw_measurements_;
  VectorXd z_diff = z - z_pred;

  z_diff(1) = atan2(sin(z_diff(1)), cos(z_diff(1))) ;

  // update state and covariance mtx
  x_ = x_ + kalman_gain * z_diff;
  P_ = P_ - kalman_gain *  S * kalman_gain.transpose();
}

