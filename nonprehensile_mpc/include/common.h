#include <ros/ros.h>
#include <sensor_msgs/Image.h>
#include <std_msgs/Float32MultiArray.h>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>

using namespace std;
using namespace Eigen;

extern MatrixXd Gu;
extern MatrixXd Go;
extern MatrixXd Sn;
extern MatrixXd Snn;
extern MatrixXd Gamma;
extern Vector3d g;

extern double m;
extern double I1;
extern double I2;
extern double I3;
extern double mu;
extern DiagonalMatrix<double, 3> I;

extern Vector3d r1;
extern Vector3d r2;
extern Vector3d r3;
extern Vector3d r4;
extern MatrixXd r;

extern const int N;
extern double dt;
extern double t_max;
extern double vt;
extern double vr;
extern double fc;
extern double epsilon;
extern double kappa_u;
extern double alpha_u;
extern double kappa_o;
extern double alpha_o;
extern double rho_u;
extern double rho_o;
extern double uof_ub[];
extern double uof_lb[];

void init ();

Matrix3d skewMat (const Vector3d& w);
Vector3d skewVec (const Matrix3d& M);
MatrixXd kroneckerProduct (const MatrixXd& A, const MatrixXd& B);
