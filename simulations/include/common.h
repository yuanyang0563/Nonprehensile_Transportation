#include <ros/ros.h>
#include <sensor_msgs/Image.h>
#include <sensor_msgs/JointState.h>
#include <std_msgs/Float64MultiArray.h>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
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

MatrixXd Gu(3,12);
MatrixXd Go(3,12);
MatrixXd Sn(3*N,3);
MatrixXd Snn(N,N);
MatrixXd Gamma(N,N);
Vector3d g(0.0,0.0,9.81);

double m = 0.5;
double I1 = 0.1;
double I2 = 0.1;
double I3 = 0.1;
double mu = 0.1;
DiagonalMatrix<double,3> I;

Vector3d r1( 0.075,-0.075,-0.075);
Vector3d r2( 0.075, 0.075,-0.075);
Vector3d r3(-0.075, 0.075,-0.075);
Vector3d r4(-0.075,-0.075,-0.075);
MatrixXd r(4,3);

double dt = 0.033;
double vt = 1.0;
double vr = 1.0;
double fc = 1e5;
double epsilon = 0.01;
double kappa_u = 1.0;
double alpha_u = 300.0;
double kappa_o = 1.0;
double alpha_o = 300.0;
double rho_u = 0.1;
double beta_u = 5.0;
double rho_o = 0.1;
double beta_o = 5.0;
double gamma_v = 0.05;
double kappa_f = 0.1;
double rho_f = 3.0;

Matrix3d skewMat (const Vector3d& w);
Vector3d skewVec (const Matrix3d& M);
MatrixXd kroneckerProduct (const MatrixXd& A, const MatrixXd& B);

void init () {
  I.diagonal() << I1, I2, I3;
  r << r1.transpose(), r2.transpose(), r3.transpose(), r4.transpose();
  Gamma = MatrixXd::Identity(N,N);
  for (int n=1; n<N; ++n)
	Gamma(n,n-1) = -1.0;
  Sn.setZero();
  Snn.setZero();
  for (int n=0; n<N; ++n) {
	MatrixXd sn = MatrixXd::Zero(N,1);
	for (int i=0; i<=n; ++i)
		sn(i) = 1.0;
	Sn += kroneckerProduct(sn,Matrix3d::Identity());
	Snn += sn*sn.transpose();
  }
  for (int i=0; i<4; ++i) {
	Gu.block(0,3*i,3,3) = Matrix3d::Identity();
	Go.block(0,3*i,3,3) = skewMat(r.row(i));
  }
}


Matrix3d skewMat (const Vector3d& w) {
    Matrix3d W;
    W << 0.0, -w(2), w(1), w(2), 0.0, -w(0), -w(1), w(0), 0.0;
    return W;
}

Vector3d skewVec (const Matrix3d& M) {
    Vector3d m;
    m << M(2,1), M(0,2), M(1,0);
    return m;
}

MatrixXd kroneckerProduct (const MatrixXd& A, const MatrixXd& B) {
    MatrixXd result(A.rows()*B.rows(), A.cols()*B.cols());
    for (int i = 0; i < A.rows(); ++i) {
        for (int j=0; j<A.cols(); ++j)
            result.block(i * B.rows(), j * B.cols(), B.rows(), B.cols()) = A(i, j) * B;
    }
    return result;
}
