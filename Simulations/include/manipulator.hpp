#include "common.h"

using namespace std;
using namespace Eigen;

class manipulator {

  public:
  
	VectorXd q, dq, bu_d, bo_d, eta, lambda, h;
	Vector3d base, upsilon, omega;
	Vector3d x, xeo, xd;
	Matrix3d R, Reo, Rd;
	MatrixXd J, Au_d, Ao_d, Hu, Ho, Hf;
	int jointHandles[6];
	
	manipulator() {
		q = VectorXd(6);
		dq = VectorXd(6);
		upsilon.setZero();
		omega.setZero();
		J = MatrixXd(6,6);
		Hu = MatrixXd(6*N,3*N);
		Ho = MatrixXd(6*N,3*N);
		Hf = MatrixXd(6*N,12*N);
		eta = VectorXd::Zero(3*N);
		lambda = VectorXd::Zero(3*N);
		h = VectorXd::Zero(6*N);
		Au_d = MatrixXd(3*N,3*N);
		bu_d = VectorXd(3*N);
		Ao_d = MatrixXd(3*N,3*N);
		bo_d = VectorXd(3*N);
	}

	void get_pose_jacobian () {
		MatrixXd a(6,1), alpha(6,1), d(6,1);
		a << 0.0, -0.24365, -0.21325, 0.0, 0.0, 0.0;
		alpha << M_PI/2.0, 0.0, 0.0, M_PI/2.0, -M_PI/2.0, 0.0;
		d << 0.1519, 0.0, 0.0, 0.11235, 0.08535, 0.0819;
		Matrix4d A;
		MatrixXd T = Matrix4d::Identity();
		Matrix<double, 3, 7> o, z;
		o.col(0) << 0.0, 0.0, 0.0;
		z.col(0) << 0.0, 0.0, 1.0;
		for (int i=0; i<6; ++i) {
			A << cos(q(i)), -sin(q(i))*cos(alpha(i)),  sin(q(i))*sin(alpha(i)), a(i)*cos(q(i)),
	     	     	     sin(q(i)),  cos(q(i))*cos(alpha(i)), -cos(q(i))*sin(alpha(i)), a(i)*sin(q(i)),
	     	             0.0,        sin(alpha(i)),            cos(alpha(i)),           d(i),
	     	             0.0,        0.0,                      0.0,                     1.0;
	     		T *= A;
			o.col(i+1) << T.block(0,3,3,1);
			z.col(i+1) << T.block(0,2,3,1);
		}
		x = T.block(0,3,3,1)+base;
		R = T.block(0,0,3,3);
		for (int i=0; i<6; ++i)
			J.col(i)<< z.col(i).cross(o.col(6)-o.col(i)), R.transpose()*z.col(i);
	}
	
	void move_one_step () {
		dq << upsilon, omega;
		dq = J.transpose()*(J*J.transpose()).inverse()*dq;
		q += dt*dq;
	}
	
	void set_tar_pars () {
		Au_d = alpha_u*MatrixXd::Identity(3*N,3*N)+kappa_u*kroneckerProduct(Snn,Matrix3d::Identity());
		Ao_d = alpha_o*MatrixXd::Identity(3*N,3*N);
		Hu.setZero();
		Hu.block(0,0,3*N,3*N) = m*kroneckerProduct(Gamma,Matrix3d::Identity());
	}
	
	void update_tar_pars () {
		bu_d = kappa_u*Sn*(x-xd);
		bo_d = kappa_o*Sn*skewVec(Rd.transpose()*R);
		eta.head(3) = upsilon-R*xeo.cross(omega);
		lambda.head(3) = I*Reo.transpose()*omega;
		h.head(3*N) = m*eta-dt*m*kroneckerProduct(MatrixXd::Ones(N,1),R*skewMat(omega)*skewMat(omega)*xeo+g);
		h.tail(3*N) = lambda-dt*kroneckerProduct(MatrixXd::Ones(N,1),skewMat(Reo.transpose()*omega)*I*Reo.transpose()*omega);
		Ho << m*kroneckerProduct(Gamma,R*skewMat(xeo)), -kroneckerProduct(Gamma,I*Reo.transpose());
		Hf << kroneckerProduct(MatrixXd::Identity(N,N),R*Reo*Gu), kroneckerProduct(MatrixXd::Identity(N,N),Go);
	}
	
	inline double cost () {
		return (xd-x).norm()-(Rd*R.transpose()).trace()+3.0;
	}

};

