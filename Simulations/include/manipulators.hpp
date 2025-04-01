#include "manipulator.hpp"

using namespace std;
using namespace Eigen;

class manipulators {

  public:
  
	manipulator left;
	manipulator right;
	
	MatrixXd A_d, A_s, A_v;
	VectorXd b_d, b_s, b_v;
	
	manipulators (string arm_left, string arm_right) : left(arm_left), right(arm_right) {
		A_d = MatrixXd::Zero(12*N,12*N);
		b_d = VectorXd::Zero(12*N);
		A_s = MatrixXd::Zero(12*N,12*N);
		b_s = VectorXd::Zero(12*N);
		A_v = MatrixXd::Zero(12*N,12*N);
		b_v = VectorXd::Zero(12*N);
		left.display->setWindowPosition(0,0);
		right.display->setWindowPosition(0,640);
	}

	void get_pose_jacobian () {
		left.get_pose_jacobian();
		right.get_pose_jacobian();
	}
	
	void move_one_step () {
		left.move_one_step();
		right.move_one_step();
	}
	
	void set_tar_pars () {
		left.set_tar_pars();
		right.set_tar_pars();
		A_d.block(0,0,3*N,3*N) = left.Au_d;
		A_d.block(3*N,3*N,3*N,3*N) = right.Au_d;
		A_d.block(6*N,6*N,3*N,3*N) = left.Ao_d;
		A_d.block(9*N,9*N,3*N,3*N) = right.Ao_d;
	}
	
	void update_tar_pars () {
		left.update_tar_pars();
		right.update_tar_pars();
		b_d << left.bu_d, right.bu_d, left.bo_d, right.bo_d;
	}
	
	void set_syn_pars () {
		A_s.block(0,0,3*N,3*N) = rho_u*kroneckerProduct(Snn,Matrix3d::Identity());
		A_s.block(0,3*N,3*N,3*N) = -rho_u*kroneckerProduct(Snn,left.Rd*right.Rd.transpose());
		A_s.block(3*N,0,3*N,3*N) = -rho_u*kroneckerProduct(Snn,right.Rd*left.Rd.transpose());
		A_s.block(3*N,3*N,3*N,3*N) = rho_u*kroneckerProduct(Snn,Matrix3d::Identity());
	}

	void update_syn_pars () {
		b_s.segment(0,3*N) = rho_u*Sn*(left.x-left.xd-left.Rd*right.Rd.transpose()*(right.x-right.xd));
		b_s.segment(3*N,3*N) = rho_u*Sn*(right.x-right.xd-right.Rd*left.Rd.transpose()*(left.x-left.xd));
		b_s.segment(6*N,3*N) = rho_o*Sn*skewVec(right.R.transpose()*right.Rd*left.Rd.transpose()*left.R);
		b_s.segment(9*N,3*N) = rho_o*Sn*skewVec(left.R.transpose()*left.Rd*right.Rd.transpose()*right.R);
		A_s.block(6*N,9*N,3*N,3*N) = 0.5*rho_o*kroneckerProduct(Snn,right.R.transpose()*right.Rd*left.Rd.transpose()*left.R-(right.R.transpose()*right.Rd*left.Rd.transpose()*left.R).trace()*Matrix3d::Identity());
		A_s.block(9*N,6*N,3*N,3*N) = 0.5*rho_o*kroneckerProduct(Snn,left.R.transpose()*left.Rd*right.Rd.transpose()*right.R-(left.R.transpose()*left.Rd*right.Rd.transpose()*right.R).trace()*Matrix3d::Identity());
	}
	
	void update_vis_pars () {
		A_v.setZero();
		b_v.setZero();
		MatrixXd A, B, C;
		for (int i=0; i<4; ++i) {
			A = left.R*left.Rec*left.L[i].transpose();
			B = skewMat(left.Reo*left.p[i]+left.xeo)*left.Rec*left.L[i].transpose();
			C = skewMat(right.Reo*left.p[i]+right.xeo)*right.R.transpose()*left.R*left.Rec*left.L[i].transpose();
			MatrixXd Lt_left(12*N,2);
			Lt_left << -Sn*A, Sn*A, -Sn*B, Sn*C;
			A_v += gamma_v*Lt_left*Lt_left.transpose();
			b_v += gamma_v*Lt_left*(left.zeta[i]-left.zeta_d[i]);
			A = right.R*right.Rec*right.L[i].transpose();
			B = skewMat(left.Reo*right.p[i]+left.xeo)*left.R.transpose()*right.R*right.Rec*right.L[i].transpose();
			C = skewMat(right.Reo*right.p[i]+right.xeo)*right.Rec*right.L[i].transpose();
			MatrixXd Lt_right(12*N,2);
			Lt_right << Sn*A, -Sn*A, Sn*B, -Sn*C;
			A_v += gamma_v*Lt_right*Lt_right.transpose();
			b_v += gamma_v*Lt_right*(right.zeta[i]-right.zeta_d[i]);
		}
	}
	
	inline double cost () {
		return left.cost()+right.cost();
	}

};





