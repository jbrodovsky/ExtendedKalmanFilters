#include "tools.h"
#include <iostream>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;
using std::cout;
using std::endl;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations, const vector<VectorXd> &ground_truth) {
  /**
   * TODO: Calculate the RMSE here.
   */
  // Based off of lesson 25.24
  //std::cout << "Im calculating RMSE!" << std::endl;
  
  VectorXd rmse(4);
  rmse << 0,0,0,0;
  if (estimations.size() != ground_truth.size() || estimations.size() == 0){
    std::cout << "Invalid data" << std::endl;
    return rmse;
  }
  
  for (unsigned int i = 0; i < estimations.size(); ++i){
    VectorXd residual = estimations[i] - ground_truth[i];
    residual = residual.array() * residual.array();		// Coefficient-wise multiplication
    rmse += residual;
  }
  
  rmse = rmse/estimations.size();		// mean
  rmse = rmse.array().sqrt();			// square root
  return rmse;
  
}

//MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {		//original
MatrixXd Tools::CalculateJacobian(VectorXd& x_state){
/**
   * TODO:
   * Calculate a Jacobian here.
   */
  // from Lesson 25.20
  
  MatrixXd Hj(3,4);
  // recover state parameters
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);

  // pre-compute a set of terms to avoid repeated calculation
  float c1 = px*px+py*py;
  float c2 = sqrt(c1);
  float c3 = (c1*c2);

  // check division by zero
  if (fabs(c1) < 0.0001) {
    cout << "CalculateJacobian () - Error - Division by Zero" << endl;
    return Hj;
  }

  // compute the Jacobian matrix
  Hj << (px/c2), (py/c2), 0, 0,
      -(py/c1), (px/c1), 0, 0,
      py*(vx*py - vy*px)/c3, px*(px*vy - py*vx)/c3, px/c2, py/c2;
  //cout<<"Hj: "<<endl<<Hj<<endl;
  return Hj;
  
}
