#include "manipulator.hpp"

class manipulator_dual {

  public:
  
	manipulator left;
	manipulator right;
	
	MatrixXf A_obj, A_ls, A_rs, A_cs;
	VectorXf b_obj, b_ls, b_rs;
	
	Vector3f x_lr, x0_lr, x_rl, x0_rl;
	Matrix3f R_lr, R0_lr, R_rl, R0_rl;
	
	float mode;
	
	double uof_ub[24*N], uof_lb[24*N];
	
	bool getJoints, getFeatures, setSynpars, isInitialized;
	
	manipulator_dual (string arm_left, string arm_right, float control_mode) : left(arm_left), right(arm_right), mode(control_mode) {
		A_obj = MatrixXf::Zero(24*N,24*N);
		b_obj = VectorXf::Zero(24*N);
		A_ls = MatrixXf::Zero(6*N,6*N);
		b_ls = VectorXf::Zero(6*N);
		A_rs = MatrixXf::Zero(6*N,6*N);
		b_rs = VectorXf::Zero(6*N);
		A_cs = MatrixXf::Zero(6*N,6*N);
		getJoints = false;
		getFeatures = false;
		setSynpars = false;
		isInitialized = false;
		if (mode==0.0) {
			gamma_v = 5.0;
			kappa_f = 0.0;
		}
	}

	void get_pose_jacobian () {
		left.get_pose_jacobian();
		right.get_pose_jacobian();
		if (!getJoints)
			getJoints = left.getJoints && right.getJoints;
		else {
			if (!setSynpars) {
				set_syn_pars();
				setSynpars = !setSynpars;
			}
			update_syn_pars();
		}
	}
	
	void move_one_step () {
		left.move_one_step();
		right.move_one_step();
	}
	
	void stop_moving () {
		left.upsilon.setZero();
		left.omega.setZero();
		right.upsilon.setZero();
		right.omega.setZero();
		move_one_step();
	}
	
	void set_tar_pars () {
		left.set_tar_pars();
		right.set_tar_pars();
	}
	
	void update_tar_pars () {
		left.update_tar_pars();
		right.update_tar_pars();
	}

	void set_syn_pars () {
		R0_lr = left.R0.transpose()*right.R0;
		R0_rl = right.R0.transpose()*left.R0;
		x0_lr = left.R0.transpose()*(right.x0-left.x0);
		x0_rl = right.R0.transpose()*(left.x0-right.x0);
		A_ls.block(0*N,0*N,3*N,3*N)  =  rho_u*kroneckerProduct(Snn,4.0*Matrix3f::Identity());
		A_ls.block(0*N,0*N,3*N,3*N) +=  beta_u*MatrixXf::Identity(3*N,3*N);
		A_ls.block(3*N,3*N,3*N,3*N)  = -rho_u*kroneckerProduct(Snn,skewMat(x0_lr)*skewMat(x0_lr));
		A_ls.block(3*N,3*N,3*N,3*N) +=  beta_o*MatrixXf::Identity(3*N,3*N);
		A_rs.block(0*N,0*N,3*N,3*N)  =  rho_u*kroneckerProduct(Snn,4.0*Matrix3f::Identity());
		A_rs.block(0*N,0*N,3*N,3*N) +=  beta_u*MatrixXf::Identity(3*N,3*N);
		A_rs.block(3*N,3*N,3*N,3*N)  = -rho_u*kroneckerProduct(Snn,skewMat(x0_rl)*skewMat(x0_rl));
		A_rs.block(3*N,3*N,3*N,3*N) +=  beta_o*MatrixXf::Identity(3*N,3*N);
		A_cs.block(0*N,0*N,3*N,3*N)  = -rho_u*kroneckerProduct(Snn,4.0*Matrix3f::Identity());
	}

	void update_syn_pars () {
		R_lr = left.R.transpose()*right.R;
		R_rl = right.R.transpose()*left.R;
		x_lr = left.R.transpose()*(right.x-left.x);
		x_rl = right.R.transpose()*(left.x-right.x);
		A_ls.block(0*N,3*N,3*N,3*N)  = -rho_u*kroneckerProduct(Snn,2.0*left.R*skewMat(x0_lr));
		A_ls.block(3*N,0*N,3*N,3*N)  =  A_ls.block(0*N,3*N,3*N,3*N).transpose();
		A_rs.block(0*N,3*N,3*N,3*N)  = -rho_u*kroneckerProduct(Snn,2.0*right.R*skewMat(x0_rl));
		A_rs.block(3*N,0*N,3*N,3*N)  =  A_rs.block(0*N,3*N,3*N,3*N).transpose();
		A_cs.block(0*N,3*N,3*N,3*N)  = -A_rs.block(0*N,3*N,3*N,3*N);
		A_cs.block(3*N,0*N,3*N,3*N)  = -A_ls.block(3*N,0*N,3*N,3*N);
		A_cs.block(3*N,3*N,3*N,3*N)  =  rho_u*kroneckerProduct(Snn,skewMat(x0_lr)*R_lr*skewMat(x0_rl));
		A_cs.block(3*N,3*N,3*N,3*N) +=  2.0*rho_o*kroneckerProduct(Snn,R_lr*R0_lr.transpose()*R_lr-R_lr*(R0_lr.transpose()*R_lr).trace());
		b_ls.head(3*N) =  rho_u*Sn*(4.0*(left.x-right.x)+2.0*(left.R*x0_lr-right.R*x0_rl));
		b_ls.tail(3*N) = -rho_u*Sn*(2.0*skewMat(x0_lr)*x_lr+skewMat(x0_lr)*R_lr*x0_rl)+2.0*rho_o*Sn*skewVec(R0_lr*R_lr.transpose());
		b_rs.head(3*N) =  rho_u*Sn*(4.0*(right.x-left.x)+2.0*(right.R*x0_rl-left.R*x0_lr));
		b_rs.tail(3*N) = -rho_u*Sn*(2.0*skewMat(x0_rl)*x_rl+skewMat(x0_rl)*R_rl*x0_lr)+2.0*rho_o*Sn*skewVec(R0_rl*R_rl.transpose());
	}
	
	void set_vis_pars () {
		left.set_vis_pars();
		right.set_vis_pars();
	}
	
	void update_vis_pars () {
		left.update_vis_pars();
		right.update_vis_pars();
		if (!getFeatures)
			getFeatures = left.getFeatures && right.getFeatures;
	}
	
	void set_cst_pars () {
		left.set_cst_pars();
		right.set_cst_pars();
	}
	
	void update_cst_pars () {
		left.update_cst_pars();
		right.update_cst_pars();
	}
	
	void set_opt_pars () {
		for (size_t i=0; i<24*N; ++i) {
			if (i<6*N)
				uof_ub[i] = vt;
			else if (i<12*N)
				uof_ub[i] = vr;
			else
				uof_ub[i] = fc;
			uof_lb[i] =-uof_ub[i];
		}
		for (size_t n=0; n<N; ++n) {
			for (size_t i=0; i<4; ++i)
				uof_lb[12*N+12*n+3*i+2] = epsilon;
		}
		set_tar_pars();
		set_vis_pars();
		set_cst_pars();
		A_obj.block(12*N,12*N,12*N,12*N) = kappa_f*MatrixXf::Identity(12*N,12*N);
	}
	
	void update_opt_pars () {
		MatrixXf A_l = A_ls;
		MatrixXf A_r = A_rs;
		VectorXf b_l = b_ls;
		VectorXf b_r = b_rs;
		update_tar_pars();
		if (mode!=0.0) {
			A_l += left.A_d;
			A_r += right.A_d;
			b_l += left.b_d;
			b_r += right.b_d;
		}
		update_vis_pars();
		if (mode!=1.0) {
			A_l += left.A_v;
			A_r += right.A_v;
			b_l += left.b_v;
			b_r += right.b_v;
		}
		A_obj.block(0*N,0*N,6*N,6*N) =  A_l;
		A_obj.block(0*N,6*N,6*N,6*N) = -A_cs;
		A_obj.block(6*N,0*N,6*N,6*N) = -A_cs.transpose();
		A_obj.block(6*N,6*N,6*N,6*N) =  A_r;
		b_obj.segment(0*N,6*N) = b_l;
		b_obj.segment(6*N,6*N) = b_r;
		update_cst_pars();
	}
	
	void store_data (float t_duration) {
		left.store_data(t_duration);
		right.store_data(t_duration);
	}

};





