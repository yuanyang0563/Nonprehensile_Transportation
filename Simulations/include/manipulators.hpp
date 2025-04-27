#include "manipulator.hpp"

using namespace std;
using namespace Eigen;

class manipulators {

  public:
  
	manipulator left;
	manipulator right;
	
	MatrixXd A_d, A_s, A_v;
	VectorXd b_d, b_s, b_v;
	
	Vector3d x0_lr, x0_rl, x_lr, x_rl;
	Matrix3d R0_lr, R0_rl, R_lr, R_rl;
	
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
		A_d.block(0*N,0*N,3*N,3*N) = left.A_d.block(0*N,0*N,3*N,3*N);
		A_d.block(3*N,3*N,3*N,3*N) = right.A_d.block(0*N,0*N,3*N,3*N);
		A_d.block(6*N,6*N,3*N,3*N) = left.A_d.block(3*N,3*N,3*N,3*N);
		A_d.block(9*N,9*N,3*N,3*N) = right.A_d.block(3*N,3*N,3*N,3*N);
	}
	
	void update_tar_pars () {
		left.update_tar_pars();
		right.update_tar_pars();
		b_d << left.b_d.head(3*N), right.b_d.head(3*N), left.b_d.tail(3*N), right.b_d.tail(3*N);
	}

	void set_syn_pars () {
		left.get_pose_jacobian();
		right.get_pose_jacobian();
		left.x0 = left.x;
		left.R0 = left.R;
		right.x0 = right.x;
		right.R0 = right.R;
		x0_lr = left.R0.transpose()*(right.x0-left.x0);
		x0_rl = right.R0.transpose()*(left.x0-right.x0);
		R0_lr = left.R0.transpose()*right.R0;
		R0_rl = right.R0.transpose()*left.R0;
		A_s.block(0*N,0*N,3*N,3*N)  =  rho_u*kroneckerProduct(Snn,4.0*Matrix3d::Identity());
		A_s.block(0*N,3*N,3*N,3*N)  = -rho_u*kroneckerProduct(Snn,4.0*Matrix3d::Identity());
		A_s.block(3*N,0*N,3*N,3*N)  = -rho_u*kroneckerProduct(Snn,4.0*Matrix3d::Identity());
		A_s.block(3*N,3*N,3*N,3*N)  =  rho_u*kroneckerProduct(Snn,4.0*Matrix3d::Identity());
		A_s.block(6*N,6*N,3*N,3*N)  = -rho_u*kroneckerProduct(Snn,skewMat(x0_lr)*skewMat(x0_lr));
		A_s.block(9*N,9*N,3*N,3*N)  = -rho_u*kroneckerProduct(Snn,skewMat(x0_rl)*skewMat(x0_rl));
		A_s.block(0*N,0*N,3*N,3*N) +=  beta_u*MatrixXd::Identity(3*N,3*N);
		A_s.block(3*N,3*N,3*N,3*N) +=  beta_u*MatrixXd::Identity(3*N,3*N);
		A_s.block(6*N,6*N,3*N,3*N) +=  beta_o*MatrixXd::Identity(3*N,3*N);
		A_s.block(9*N,9*N,3*N,3*N) +=  beta_o*MatrixXd::Identity(3*N,3*N);
	}

	void update_syn_pars () {
		x_lr = left.R.transpose()*(right.x-left.x);
		x_rl = right.R.transpose()*(left.x-right.x);
		R_lr = left.R.transpose()*right.R;
		R_rl = right.R.transpose()*left.R;
		A_s.block(0*N,6*N,3*N,3*N)  = -rho_u*kroneckerProduct(Snn,2.0*left.R*skewMat(x0_lr));
		A_s.block(0*N,9*N,3*N,3*N)  =  rho_u*kroneckerProduct(Snn,2.0*right.R*skewMat(x0_rl));
		A_s.block(3*N,6*N,3*N,3*N)  =  rho_u*kroneckerProduct(Snn,2.0*left.R*skewMat(x0_lr));
		A_s.block(3*N,9*N,3*N,3*N)  = -rho_u*kroneckerProduct(Snn,2.0*right.R*skewMat(x0_rl));
		A_s.block(6*N,0*N,3*N,3*N)  =  rho_u*kroneckerProduct(Snn,2.0*skewMat(x0_lr)*left.R.transpose());
		A_s.block(6*N,3*N,3*N,3*N)  = -rho_u*kroneckerProduct(Snn,2.0*skewMat(x0_lr)*left.R.transpose());
		A_s.block(9*N,0*N,3*N,3*N)  = -rho_u*kroneckerProduct(Snn,2.0*skewMat(x0_rl)*right.R.transpose());
		A_s.block(9*N,3*N,3*N,3*N)  =  rho_u*kroneckerProduct(Snn,2.0*skewMat(x0_rl)*right.R.transpose());
		A_s.block(6*N,9*N,3*N,3*N)  =  rho_u*kroneckerProduct(Snn,skewMat(x0_lr)*R_lr*skewMat(x0_rl));
		A_s.block(6*N,9*N,3*N,3*N) +=  rho_o*kroneckerProduct(Snn,R_lr*R0_lr.transpose()*R_lr-R_lr*(R0_lr.transpose()*R_lr).trace());
		A_s.block(9*N,6*N,3*N,3*N)  =  rho_u*kroneckerProduct(Snn,skewMat(x0_rl)*R_rl*skewMat(x0_lr));
		A_s.block(9*N,6*N,3*N,3*N) +=  rho_o*kroneckerProduct(Snn,R_rl*R0_rl.transpose()*R_rl-R_rl*(R0_rl.transpose()*R_rl).trace());
		b_s.segment(0*N,3*N)  = rho_u*Sn*(2.0*left.R*(x0_lr-x_lr)-2.0*right.R*(x0_rl-x_rl));
		b_s.segment(3*N,3*N)  = rho_u*Sn*(2.0*right.R*(x0_rl-x_rl)-2.0*left.R*(x0_lr-x_lr));
		b_s.segment(6*N,3*N)  = rho_u*Sn*(skewMat(x0_lr)*R_lr*(2.0*x_rl-x0_rl));
		b_s.segment(6*N,3*N) += rho_o*Sn*(2.0*skewVec(R0_lr*R_lr.transpose()));
		b_s.segment(9*N,3*N)  = rho_u*Sn*(skewMat(x0_rl)*R_rl*(2.0*x_lr-x0_lr));
		b_s.segment(9*N,3*N) += rho_o*Sn*(2.0*skewVec(R0_rl*R_rl.transpose()));
	}
	
	void set_vis_pars () {
		left.set_vis_pars();
		right.set_vis_pars();
	}
	
	void update_vis_pars () {
		left.update_vis_pars();
		right.update_vis_pars();
		A_v.block(0*N,0*N,3*N,3*N) = 5.0*gamma_v*kroneckerProduct(Snn,left.Lu.transpose()*left.Lu);
		A_v.block(0*N,6*N,3*N,3*N) = 5.0*gamma_v*kroneckerProduct(Snn,left.Lu.transpose()*left.Lo);
		A_v.block(3*N,3*N,3*N,3*N) = 5.0*gamma_v*kroneckerProduct(Snn,right.Lu.transpose()*right.Lu);
		A_v.block(3*N,9*N,3*N,3*N) = 5.0*gamma_v*kroneckerProduct(Snn,right.Lu.transpose()*right.Lo);
		A_v.block(6*N,0*N,3*N,3*N) = 5.0*gamma_v*kroneckerProduct(Snn,left.Lo.transpose()*left.Lu);
		A_v.block(6*N,6*N,3*N,3*N) = 5.0*gamma_v*kroneckerProduct(Snn,left.Lo.transpose()*left.Lo);
		A_v.block(9*N,3*N,3*N,3*N) = 5.0*gamma_v*kroneckerProduct(Snn,right.Lo.transpose()*right.Lu);
		A_v.block(9*N,9*N,3*N,3*N) = 5.0*gamma_v*kroneckerProduct(Snn,right.Lo.transpose()*right.Lo);
		b_v.segment(0*N,3*N) = gamma_v*Sn*left.Lu.transpose()*(left.zeta-left.zeta_d);
		b_v.segment(3*N,3*N) = gamma_v*Sn*right.Lu.transpose()*(right.zeta-right.zeta_d);
		b_v.segment(6*N,3*N) = gamma_v*Sn*left.Lo.transpose()*(left.zeta-left.zeta_d);
		b_v.segment(9*N,3*N) = gamma_v*Sn*right.Lo.transpose()*(right.zeta-right.zeta_d);
	}
	
	inline double cost () {
		return left.cost()+right.cost();
	}

};
