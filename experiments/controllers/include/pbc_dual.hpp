#include "manipulator.hpp"

class manipulator_dual {

  public:
  
  	manipulator left;
  	manipulator right;

  	float mode;
  	
  	bool getJoints, getFeatures;
  	
  	manipulator_dual (string arm_left, string arm_right, float control_mode) : left(arm_left), right(arm_right), mode(control_mode) {
		getJoints = false;
		getFeatures = false;
		if (mode==0.0)
			gamma_pbc = 5.0;
  	}
  	
  	void get_pose_jacobian () {
  		left.get_pose_jacobian();
  		right.get_pose_jacobian();
  		if (!getJoints)
  			getJoints = left.getJoints && right.getJoints;
  	}
  	
  	void move_one_step () {
  		left.move_one_step();
  		right.move_one_step();
  	}
  	
  	void stop_moving () {
  		left.stop_moving();
  		right.stop_moving();
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
	
	void get_vel_inputs () {
		Vector3f upsilon_ld = kappa_pbc*(left.xd-left.x);
		Vector3f upsilon_rd = kappa_pbc*(right.xd-right.x);
		Vector3f omega_ld = kappa_pbc*left.R.transpose()*skewVec(left.Rd*left.R.transpose());
		Vector3f omega_rd = kappa_pbc*right.R.transpose()*skewVec(right.Rd*right.R.transpose());
		Vector3f upsilon_ls = rho_pbc*(right.x-left.x-0.5*(left.R*left.R0.transpose()+right.R*right.R0.transpose())*(right.x0-left.x0));
		Vector3f upsilon_rs = rho_pbc*(left.x-right.x-0.5*(right.R*right.R0.transpose()+left.R*left.R0.transpose())*(left.x0-right.x0));
		Vector3f omega_ls = rho_pbc*left.R.transpose()*skewVec(right.R*right.R0.transpose()*left.R0*left.R.transpose());
		Vector3f omega_rs = rho_pbc*right.R.transpose()*skewVec(left.R*left.R0.transpose()*right.R0*right.R.transpose());
		VectorXf twist_lv = gamma_pbc*(left.Lm.transpose()*left.Lm).inverse()*left.Lm.transpose()*(left.zeta_d-left.zeta);
		VectorXf twist_rv = gamma_pbc*(right.Lm.transpose()*right.Lm).inverse()*right.Lm.transpose()*(right.zeta_d-right.zeta);
		left.upsilon = upsilon_ls;
		left.omega = omega_ls;
		right.upsilon = upsilon_rs;
		right.omega = omega_rs;
		if (mode!=0) {
			left.upsilon += upsilon_ld;
			left.omega += omega_ld;
			right.upsilon += upsilon_rd;
			right.omega += omega_rd;
		}
		if (mode!=1) {
			left.upsilon += twist_lv.head(3);
			left.omega += twist_lv.tail(3);
			right.upsilon += twist_rv.head(3);
			right.omega += twist_rv.tail(3);
		}
		left.upsilon = left.upsilon.array().min(vt).max(-vt);
		left.omega = left.omega.array().min(vr).max(-vr);
		right.upsilon = right.upsilon.array().min(vt).max(-vt);
		right.omega = right.omega.array().min(vr).max(-vr);
	}
  	
  	void store_data (float t_duration) {
		left.store_data(t_duration);
		right.store_data(t_duration);
	}

};
