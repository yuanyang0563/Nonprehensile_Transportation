#include "common.h"

using namespace std;
using namespace Eigen;

class manipulator {

  public:
  
  	VectorXd a, alpha, d;
	VectorXd q, dq;
	VectorXd bu_d, bo_d, eta, lambda, h;
	MatrixXd J, Au_d, Ao_d, Hu, Ho, Hf;
	Vector3d upsilon, omega;
	Vector3d x, xb, xeo, xd;
	Matrix3d R, Rb, Reo, Rd;
	int jointHandles[6];
	
	manipulator() {
		a = VectorXd(6);
		alpha = VectorXd(6);
		d = VectorXd(6);
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
		x = Rb*T.block(0,3,3,1)+xb;
		R = Rb*T.block(0,0,3,3);
		for (int i=0; i<6; ++i)
			J.col(i)<< Rb*z.col(i).cross(o.col(6)-o.col(i)), R.transpose()*Rb*z.col(i);
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

