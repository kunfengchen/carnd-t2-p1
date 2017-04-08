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
  DONE:
    * Finish initializing the FusionEKF.
  */
  //create a 4D state vector, we don't know yet the values of the x state
  ekf_.x_ = VectorXd(4);

  //state covariance matrix P
  ekf_.P_ = MatrixXd(4, 4);
  ekf_.P_ << 1, 0, 0, 0,
          0, 1, 0, 0,
          0, 0, 1000, 0,
          0, 0, 0, 1000;


  //measurement covariance
  ekf_.R_ = R_laser_;

  ekf_.R_radar_ = R_radar_;

  //measurement matrix
  H_laser_ << 1, 0, 0, 0,
          0, 1, 0, 0;
  ekf_.H_ = H_laser_;

  //the initial transition matrix F_
  ekf_.F_ = MatrixXd(4, 4);
  ekf_.F_ << 1, 0, 1, 0,
          0, 1, 0, 1,
          0, 0, 1, 0,
          0, 0, 0, 1;

  ekf_.Q_ = MatrixXd(4, 4);

  //set the acceleration noise components
  noise_ax = 9;
  noise_ay = 9;

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
    DONE:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    cout << "EKF Init: " << endl << ekf_.x_ << endl;
    ekf_.x_ = VectorXd(4);

    ekf_.x_ << 1, 1, 1, 1;

    previous_timestamp_ = measurement_pack.timestamp_;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      double ro = measurement_pack.raw_measurements_[0];
      double theta = measurement_pack.raw_measurements_[1];
      double x = ro * cos(theta);
      double y = ro * sin(theta);
      if (x == 0 or y ==0) {
        cout << "SKIP the zeros !!!" << endl;
      }
      ekf_.x_ << x, y, 0, 0;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
    }

    if (ekf_.x_.norm() > 1e-4) {
      // done initializing, no need to predict or update
      is_initialized_ = true;
    } else {
      cout << "input is too SMALL. Not inited." << endl;
    }

    // KFC print the output
    cout << "x_ = " << endl << ekf_.x_ << endl;
    cout << "P_ = " << endl << ekf_.P_ << endl;
    cout << "F_ = " << endl << ekf_.F_ << endl;
    cout << "Q_ = " << endl << ekf_.Q_ << endl;
    cout << "H_ = " << endl << ekf_.H_ << endl;
    cout << "R_ = " << endl << ekf_.R_ << endl;
    cout << "R_Radar_ = " << endl << ekf_.R_radar_ << endl;

    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   DONE:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
   */
  //compute the time elapsed between the current and previous measurements
  cout << " timestamps " << measurement_pack.timestamp_<< " pre = " << previous_timestamp_ << endl;
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
  cout << "dt = " << dt << endl;
  previous_timestamp_ = measurement_pack.timestamp_;

  if ( dt < 0.001 ) {
    cout << "dt is too SMALL!!!" << endl;
  }
  double dt_2 = dt * dt;
  double dt_3 = dt_2 * dt;
  double dt_4 = dt_3 * dt;

  //Modify the F matrix so that the time is integrated
  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;

  //set the process covariance matrix Q
  ekf_.Q_ <<  dt_4/4*noise_ax, 0, dt_3/2*noise_ax, 0,
          0, dt_4/4*noise_ay, 0, dt_3/2*noise_ay,
          dt_3/2*noise_ax, 0, dt_2*noise_ax, 0,
          0, dt_3/2*noise_ay, 0, dt_2*noise_ay;

  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   DONE:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else {
    // Laser updates
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << endl << ekf_.x_ << endl;
  cout << "P_ = " << endl << ekf_.P_ << endl;
  cout << "F_ = " << endl << ekf_.F_ << endl;
  cout << "Q_ = " << endl << ekf_.Q_ << endl;
  cout << "H_ = " << endl << ekf_.H_ << endl;
  cout << "R_ = " << endl << ekf_.R_ << endl;
  cout << "R_Radar_ = " << endl << ekf_.R_radar_ << endl;
}
