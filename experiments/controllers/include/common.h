#include <ros/ros.h>
#include <sensor_msgs/Image.h>
#include <sensor_msgs/JointState.h>
#include <geometry_msgs/Pose.h>
#include <std_msgs/Float64MultiArray.h>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Geometry>
#include <visp3/core/vpImage.h>
#include <visp3/core/vpConfig.h>
#include <visp3/core/vpDisplay.h>
#include <visp3/gui/vpDisplayX.h>
#include <visp3/core/vpCameraParameters.h>
#include <visp3/core/vpPixelMeterConversion.h>
#include <visp3/detection/vpDetectorAprilTag.h>
#include <chrono>
#include <thread>
#include <iostream>
#include <fstream>
#include <sstream>
#include <ctime>

using namespace std;
using namespace Eigen;


const int N = 5;

MatrixXf Gu(3,12);
MatrixXf Go(3,12);
MatrixXf Sn(3*N,3);
MatrixXf Snn(N,N);
MatrixXf Gamma(N,N);
Vector3f g(0.0,0.0,9.81);

float m = 0.5;
float I1 = 0.1;
float I2 = 0.1;
float I3 = 0.1;
float mu = 0.1;
DiagonalMatrix<float,3> I;

Vector3f r1( 0.075,-0.075,-0.075);
Vector3f r2( 0.075, 0.075,-0.075);
Vector3f r3(-0.075, 0.075,-0.075);
Vector3f r4(-0.075,-0.075,-0.075);
MatrixXf r(4,3);

float dt = 0.01;
float vt = 1.0;
float vr = 1.0;
float fc = 1e5;
float epsilon = 0.01;
float kappa_u = 1.0;
float alpha_u = 300.0;
float kappa_o = 0.25;
float alpha_o = 150.0;
float rho_u = 2.0;
float beta_u = 50.0;
float rho_o = 2.0;
float beta_o = 50.0;
float gamma_v = 0.05;
float kappa_f = 0.5;
float rho_f = 2.0;
float kappa_pbc = 1.0;
float rho_pbc = 5.0;
float gamma_pbc = 0.5;

Matrix3f skewMat (const Vector3f& w);
Vector3f skewVec (const Matrix3f& M);
MatrixXf kroneckerProduct (const MatrixXf& A, const MatrixXf& B);

void init () {
  I.diagonal() << I1, I2, I3;
  r << r1.transpose(), r2.transpose(), r3.transpose(), r4.transpose();
  Gamma = MatrixXf::Identity(N,N);
  for (size_t n=1; n<N; ++n)
	Gamma(n,n-1) = -1.0;
  Sn.setZero();
  Snn.setZero();
  for (size_t n=0; n<N; ++n) {
	MatrixXf sn = MatrixXf::Zero(N,1);
	for (size_t i=0; i<=n; ++i)
		sn(i) = 1.0;
	Sn += kroneckerProduct(sn,Matrix3f::Identity());
	Snn += sn*sn.transpose();
  }
  for (size_t i=0; i<4; ++i) {
	Gu.block(0,3*i,3,3) = Matrix3f::Identity();
	Go.block(0,3*i,3,3) = skewMat(r.row(i));
  }
}


Matrix3f skewMat (const Vector3f& w) {
    Matrix3f W;
    W << 0.0, -w(2), w(1), w(2), 0.0, -w(0), -w(1), w(0), 0.0;
    return W;
}

Vector3f skewVec (const Matrix3f& M) {
    Vector3f m;
    m << M(2,1), M(0,2), M(1,0);
    return m;
}

MatrixXf kroneckerProduct (const MatrixXf& A, const MatrixXf& B) {
    MatrixXf result(A.rows()*B.rows(), A.cols()*B.cols());
    for (size_t i = 0; i < A.rows(); ++i) {
        for (size_t j=0; j<A.cols(); ++j)
            result.block(i*B.rows(), j*B.cols(), B.rows(), B.cols()) = A(i,j)*B;
    }
    return result;
}
