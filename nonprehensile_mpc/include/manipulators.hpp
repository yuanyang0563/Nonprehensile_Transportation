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
		A_d.block(0,0,3*N,3*N) = left.A_d.block(0,0,3*N,3*N);
		A_d.block(3*N,3*N,3*N,3*N) = right.A_d.block(0,0,3*N,3*N);
		A_d.block(6*N,6*N,3*N,3*N) = left.A_d.block(3*N,3*N,3*N,3*N);
		A_d.block(9*N,9*N,3*N,3*N) = right.A_d.block(3*N,3*N,3*N,3*N);
	}
	
	void update_tar_pars () {
		left.update_tar_pars();
		right.update_tar_pars();
		b_d << left.b_d.head(3*N), right.b_d.head(3*N), left.b_d.tail(3*N), right.b_d.tail(3*N);
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
	
	void set_vis_pars () {
		left.set_vis_pars();
		right.set_vis_pars();
	}
	
	void update_vis_pars () {
		left.update_vis_pars();
		right.update_vis_pars();
		A_v.block(0,0,3*N,3*N) = 5.0*kroneckerProduct(Snn,left.Lu.transpose()*left.Lu);
		A_v.block(0,6*N,3*N,3*N) = 5.0*kroneckerProduct(Snn,left.Lu.transpose()*left.Lo);
		A_v.block(3*N,3*N,3*N,3*N) = 5.0*kroneckerProduct(Snn,right.Lu.transpose()*right.Lu);
		A_v.block(3*N,9*N,3*N,3*N) = 5.0*kroneckerProduct(Snn,right.Lu.transpose()*right.Lo);
		A_v.block(6*N,0,3*N,3*N) = 5.0*kroneckerProduct(Snn,left.Lo.transpose()*left.Lu);
		A_v.block(6*N,6*N,3*N,3*N) = 5.0*kroneckerProduct(Snn,left.Lo.transpose()*left.Lo);
		A_v.block(9*N,3*N,3*N,3*N) = 5.0*kroneckerProduct(Snn,right.Lo.transpose()*right.Lu);
		A_v.block(9*N,9*N,3*N,3*N) = 5.0*kroneckerProduct(Snn,right.Lo.transpose()*right.Lo);
		b_v.segment(0,3*N) = Sn*left.Lu.transpose()*(left.zeta-left.zeta_d);
		b_v.segment(3*N,3*N) = Sn*right.Lu.transpose()*(right.zeta-right.zeta_d);
		b_v.segment(6*N,3*N) = Sn*left.Lo.transpose()*(left.zeta-left.zeta_d);
		b_v.segment(9*N,3*N) = Sn*right.Lo.transpose()*(right.zeta-right.zeta_d);
	}
	
	inline double cost () {
		return left.cost()+right.cost();
	}

};





