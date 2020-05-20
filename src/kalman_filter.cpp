#include "kalman_filter.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;


/* 
 * Please note that the Eigen library does not initialize 
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

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
VectorXd CartToPolar(const VectorXd &v_cart){
  const double & px = v_cart(0);
  const double & py = v_cart(1);
  const double & vx = v_cart(2);
  const double & vy = v_cart(3);

  double rho, phi, rho_dot;
  rho = sqrt(px*px + py*py);
  phi = atan2(py, px);

  // protection from division by zero
  if (rho < 0.000001) {
    rho = 0.000001;
  }

  rho_dot = (px * vx + py * vy) / rho;

  VectorXd v_polar = VectorXd(3);
  v_polar << rho, phi, rho_dot;

  return v_polar;

}

void KalmanFilter::Predict() {
 
  x_ = F_ * x_; // u is zero vector
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) 
{
    I_= MatrixXd::Identity(4, 4);
    VectorXd y;
   
    y = z - H_*x_;
   
  
   MatrixXd Ht_ = H_.transpose();
   MatrixXd psi = H_ * P_ * Ht_ + R_;
   MatrixXd psi_inv = psi.inverse();
   MatrixXd K =  P_ * Ht_ * psi_inv;
  
  // new state
  x_ = x_ + (K * y);
  P_ = (I_ - K * H_) * P_;
    
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  
  
   I_= MatrixXd::Identity(4, 4);
   MatrixXd Hj_ = tools_.CalculateJacobian(x_);
   VectorXd y;
   VectorXd x_polar =CartToPolar(x_);
   y = z - x_polar;
    // Round up the angle between -pi to pi
    // This makes the function discontinous though
    while(y(1) > M_PI){
      y(1) -= 2 * M_PI;
    }

    while(y(1) < -M_PI){
      y(1) += 2 * M_PI;
    }
  
   MatrixXd Hjt_ = Hj_.transpose();
   MatrixXd psi = Hj_ * P_ * Hjt_ + R_; 
   MatrixXd psi_inv = psi.inverse(); 
   MatrixXd K =  P_ * Hjt_ * psi_inv;
    // new state
   x_ = x_ + (K * y);
   P_ = (I_ - K * Hj_) * P_;

  
}
