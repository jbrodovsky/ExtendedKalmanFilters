#include "kalman_filter.h"
#include <math.h>
#include "tools.h"
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

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in, MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
   * TODO: predict the state
   */
  
  /*
  Error occuring in prediction, appears to be matrix dimensions
  */
 // std::cout<<"In KF predict"<<std::endl;
  x_ = F_*x_;					//prime state
  //std::cout<<"xp: "<<x_<<std::endl;
  //std::cout<<"Transition: "<<std::endl;
  //std::cout<<F_<<std::endl;
  
  P_ = F_*P_*F_.transpose() + Q_;	//prime covariance
  //std::cout<<"Pp: "<<P_<<std::endl;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
   * TODO: update the state by using Kalman Filter equations
   */
  VectorXd y = z - H_*x_;
  MatrixXd S = H_*P_*H_.transpose() + R_;
  MatrixXd K = P_*H_.transpose()*S.inverse();
  x_ = x_ + (K*y);
  P_ = (MatrixXd::Identity(x_.size(),x_.size()) - K*H_)*P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
   * TODO: update the state by using Extended Kalman Filter equations
   */
  
  // Measurement Matrix H is toggled to Hj by a conditional block in FusionEKF.ProcessMeasurement()
  
  //Calculate h(x)
  // use h function in 25.18
  VectorXd h = VectorXd(3);
  h << sqrt(pow(x_(0), 2) + pow(x_(1), 2)), 
       atan2(x_(1), x_(0)), 
       (x_(0)*x_(2) + x_(1)*x_(3))/(sqrt(pow(x_(0),2) + pow(x_(1),2)));
  //cout<<"h'(x) function intialized!"<<endl;
  //cout<<h<<endl;
  //cout<<"Measured angle: "<<z(1)<<endl;
  //cout<<"Calculated angle: "<<h(1)<<endl;
  VectorXd y = z - h;
  cout<<"angle error: "<<y(1)<<endl;
  
  /*
  // Constrain angle to [-pi, pi)
  float a = fmod(y(1) + 3.14159265359, 2*3.14159265359);
  if (a < 0){
    cout<<"ANGLE NEEDS TO BE NORMAILIZED!!!"<<endl;
    y(1) += 2*3.14159265359;
  }*/
  
  if(y(1)>3.14159265359){y(1) = y(1) - 2*3.14159265359;}
  else if(y(1)<-3.14159265359){y(1) = y(1) + 2*3.14159265359;}
  
  MatrixXd S = H_*P_*H_.transpose() + R_;
  MatrixXd K = P_*H_.transpose()*S.inverse();
  x_ = x_ + (K*y);
  P_ = (MatrixXd::Identity(x_.size(),x_.size()) - K*H_)*P_;
}
