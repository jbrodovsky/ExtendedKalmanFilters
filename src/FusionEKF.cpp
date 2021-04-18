#include "FusionEKF.h"
#include <iostream>
#include "Eigen/Dense"
#include "tools.h"
#include <math.h>

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
  //measurement covariance matrix - laser
  R_laser_ = MatrixXd(2, 2);
  R_laser_ << 0.0225, 0,
              0, 0.0225;
  //measurement covariance matrix - radar
  R_radar_ = MatrixXd(3, 3);
  R_radar_ << 0.09, 0,      0,
              0,    0.0009, 0,
              0,    0,      0.09;
  //measurement matrix - laser
  H_laser_ = MatrixXd(2, 4);
  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;
  //measurement matrix - radar, need to use jacobian
  Hj_ = MatrixXd(3, 4);
  Hj_ << 0, 0, 0, 0,
         0, 0, 0, 0,
         0, 0, 0, 0;
  /**
   * TODO: Finish initializing the FusionEKF.
   * TODO: Set the process and measurement noises
   */
  
  // Derived from lesson 25.13
  ekf_.x_ = VectorXd(4);
  ekf_.x_ << 0, 0, 0, 0;
  
  ekf_.P_ = MatrixXd(4,4);
  ekf_.P_ << 1, 0, 0, 0,
             0, 1, 0, 0,
             0, 0, 1000, 0,
             0, 0, 0, 1000;
  
  ekf_.F_ = MatrixXd(4,4);
  ekf_.F_ << 1, 0, 1, 0,
             0, 1, 0, 1,
             0, 0, 1, 0,
             0, 0, 0, 1;
  ekf_.Q_ = MatrixXd(4,4);
  ekf_.Q_ << 0, 0, 0, 0,
             0, 0, 0, 0,
  		     0, 0, 0, 0,
  			 0, 0, 0, 0;
  
  // set the acceleration noise components
  noise_ax = 9;
  noise_ay = 9;
}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  /**
   * Initialization
   */
  if (!is_initialized_) {
    /**
     * TODO: Initialize the state ekf_.x_ with the first measurement.
     * TODO: Create the covariance matrix.
     * You'll need to convert radar from polar to cartesian coordinates.
     */

    // first measurement
    cout << "EKF: " << endl;
	previous_timestamp_ = measurement_pack.timestamp_;
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      // TODO: Convert radar from polar to cartesian coordinates 
      //         and initialize state.
      cout << "INITIALIZING from RADAR measurement: "<<measurement_pack.raw_measurements_<< endl;
      ekf_.x_(0) = measurement_pack.raw_measurements_[0]*sin(measurement_pack.raw_measurements_[1]);
      ekf_.x_(1) = measurement_pack.raw_measurements_[0]*cos(measurement_pack.raw_measurements_[1]);
      ekf_.x_(2) = 0;//measurement_pack.raw_measurements_[2]*sin(measurement_pack.raw_measurements_[1]);
      ekf_.x_(3) = 0;//measurement_pack.raw_measurements_[2]*cos(measurement_pack.raw_measurements_[1]);
      cout << "Initial state: "<<ekf_.x_<<endl;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      // TODO: Initialize state.
      cout << "INITIALIZING from LASER measurement: "<<measurement_pack.raw_measurements_<< endl;
      ekf_.x_(0) =  measurement_pack.raw_measurements_[0];
      ekf_.x_(1) =  measurement_pack.raw_measurements_[1];
      ekf_.x_(2) =  0;//measurement_pack.raw_measurements_[2];
      ekf_.x_(3) =  0;//measurement_pack.raw_measurements_[3];
      cout << "Initial state: "<<ekf_.x_<<endl;
    }

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
  //cout<<"Input time: "<<measurement_pack.timestamp_<<endl;
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;
  
  // update the F matrix to include the time change for velocity updates
  ekf_.F_(0,2) = dt;
  ekf_.F_(1,3) = dt;
  
  //cout<<"Updated F matrix: "<<ekf_.F_<<endl;
  
  // Set the process covariance matrix Q
  float dt_2 = dt * dt;		//dt squared
  float dt_3 = dt_2 * dt;	//dt cubed
  float dt_4 = dt_3 * dt;	//dt^4
   
  ekf_.Q_ << dt_4/4*noise_ax, 0,               dt_3/2*noise_ax, 0,
             0,               dt_4/4*noise_ay, 0,               dt_3/2*noise_ay,
             dt_3/2*noise_ax, 0,               dt_2*noise_ax,   0,
             0,               dt_3/2*noise_ay, 0,               dt_2*noise_ay;
  //Execute EKF prediction
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
    ekf_.R_ = R_radar_;
    ekf_.H_ = tools.CalculateJacobian(ekf_.x_); //Hj
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } 
  else {
    // TODO: Laser updates
    //cout<<"Updating based on laser measurement"<<endl;
    ekf_.R_ = R_laser_;
    ekf_.H_ = H_laser_;
	ekf_.Update(measurement_pack.raw_measurements_);
    //cout<<"Laser update successfull"<<endl;
  }

  // print the output
  //cout << "x_ = " << ekf_.x_ << endl;
  //cout << "P_ = " << ekf_.P_ << endl;
}
