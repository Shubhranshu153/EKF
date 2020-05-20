#include "FusionEKF.h"
#include <iostream>
#include "Eigen/Dense"
#include "tools.h"

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
  P = MatrixXd(4, 4);
  /*Tune the values to optimize errors*/
  P << 1, 0, 0, 0,
  0, 1, 0, 0,
  0, 0, 1000, 0,
  0, 0, 0, 1000;

 
  H_laser_ << 1, 0, 0, 0,
  0, 1, 0, 0;
  
  // state transition matrix (initially Δt is 0)
  F=MatrixXd(4, 4);
  F << 1, 0, 0, 0,
  0, 1, 0, 0,
  0, 0, 1, 0,
  0, 0, 0, 1;

  // process covariance matrix (initially Δt is 0, hence Q consists of 0's;
  MatrixXd Q(4, 4);
  Q<< 0, 0, 0, 0,
  0, 0, 0, 0,
  0, 0, 0, 0,
  0, 0, 0, 0;

  noise_ax=9.0;
  noise_ay=9.0;

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

  
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;
    previous_timestamp_=measurement_pack.timestamp_;
    Eigen::VectorXd x_=VectorXd(4);
    
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
     
     
      double rho     = measurement_pack.raw_measurements_[0];
      double phi     = measurement_pack.raw_measurements_[1];
      double px      = rho     * cos(phi);
      double py      = rho     * sin(phi);

      // initial state in case the first measurement comes from radar sensor
      x_ << px,py,0,0; // We dont know phi_dot so we cant predict velocity in cartisan 
      ekf_.Init(x_, P, F, H_laser_, R_radar_, Q);

    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
     
       x_ << measurement_pack.raw_measurements_[0],measurement_pack.raw_measurements_[1],0.0,0.0;    // we have no info about the velocity for lidar measurement
       ekf_.Init(x_, P, F, H_laser_, R_laser_, Q);

    }

    is_initialized_ = true;
    return;
  }
 
  
  double dt = ( measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds

  previous_timestamp_ = measurement_pack.timestamp_;
 

   ekf_.F_(0, 2) = dt;
   ekf_.F_(1, 3) = dt;
 
  
  double dt_2=dt*dt;
  double dt_3=dt*dt*dt;
  double dt_4=dt*dt*dt*dt;
 

  // update the process covariance matrix
  ekf_.Q_ = MatrixXd(4, 4);

  ekf_.Q_ <<  (dt_4/4)*noise_ax,              0.0, (dt_3/2)*noise_ax,              0.0,
                           0.0, (dt_4/4)*noise_ay,              0.0, (dt_3/2)*noise_ay,
              (dt_3/2)*noise_ax,              0.0,   dt_2*noise_ax,              0.0,
                           0.0, (dt_3/2)*noise_ay,              0.0,   dt_2*noise_ay;


  ekf_.Predict();



  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {

	
    ekf_.R_= R_radar_;
     ekf_.UpdateEKF(measurement_pack.raw_measurements_);

  } else {
   
 
     ekf_.R_= R_laser_;
	 ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
