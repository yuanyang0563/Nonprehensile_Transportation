#include "manipulator.hpp"

struct Arm {
	Vector3f x, x0;
	Matrix3f R, R0;
	VectorXf f;
	Arm () {
		f = VectorXf(12);
	}
};

class manipulator_side : public manipulator {

  public:
  	ros::Publisher pub_peer;
  	ros::Subscriber sub_peer;
  	std_msgs::Float64MultiArray pofo_msg;
  	
  	Arm self, peer;
  	bool getPeer;
  	
  	MatrixXf A_s, A_obj;
  	VectorXf b_s, b_obj;
  	
  	Vector3f x_sp, x0_sp, x0_ps;
  	Matrix3f R_sp, R0_sp;
  	
  	float mode;
  	
  	double uof_ub[18*N], uof_lb[18*N];

	manipulator_side (string name_self, string name_peer, float control_mode) : manipulator (name_self), mode(control_mode) {
		A_s = MatrixXf::Zero(6*N,6*N);
		b_s = VectorXf::Zero(6*N);
		A_obj = MatrixXf::Zero(18*N,18*N);
		b_obj = VectorXf::Zero(18*N);
		pub_peer = nh.advertise<std_msgs::Float64MultiArray>(name_self+"/pofo",1);
		sub_peer = nh.subscribe(name_peer+"/pofo",1,&manipulator_side::peer_callback,this);
		pofo_msg.data.resize(19);
		getPeer = false;
		if (mode==0.0) {
			gamma_v = 5.0;
			kappa_f = 0.0;
		}
	}
	
	void get_pose_jacobian () override {
		manipulator::get_pose_jacobian();
		self.x = x;
		self.R = R;
		self.x0 = x0;
		self.R0 = R0;
	}
	
	void move_one_step () override {
		pofo_msg.data[0] = x(0);
		pofo_msg.data[1] = x(1);
		pofo_msg.data[2] = x(2);
		Quaternionf quaternion(R);
		pofo_msg.data[3] = quaternion.w();
		pofo_msg.data[4] = quaternion.x();
		pofo_msg.data[5] = quaternion.y();
		pofo_msg.data[6] = quaternion.z();
		for (size_t i=0; i<12; ++i)
			pofo_msg.data[7+i] = f(i);
		pub_peer.publish(pofo_msg);
		manipulator::move_one_step();
	}
	
	void peer_callback (const std_msgs::Float64MultiArray::ConstPtr& msg) {
		peer.x << msg->data[0], msg->data[1], msg->data[2];
		Quaternionf quaternion (msg->data[3], msg->data[4], msg->data[5], msg->data[6]);
		peer.R << quaternion.toRotationMatrix();
		for (size_t i=0; i<12; ++i)
			peer.f(i) = msg->data[7+i];
		if (!getPeer) {
			peer.x0 = peer.x;
			peer.R0 = peer.R;
			set_syn_pars();
			getPeer = !getPeer;
		} else {
			update_syn_pars();
		}
	}
	
	void set_syn_pars () {
		R0_sp = self.R0.transpose()*peer.R0;
		x0_sp = self.R0.transpose()*(peer.x0-self.x0);
		x0_ps = peer.R0.transpose()*(self.x0-peer.x0);
		A_s.block(0*N,0*N,3*N,3*N)  =  rho_u*kroneckerProduct(Snn,4.0*Matrix3f::Identity());
		A_s.block(0*N,0*N,3*N,3*N) +=  beta_u*MatrixXf::Identity(3*N,3*N);
		A_s.block(3*N,3*N,3*N,3*N)  = -rho_u*kroneckerProduct(Snn,skewMat(x0_sp)*skewMat(x0_sp));
		A_s.block(3*N,3*N,3*N,3*N) +=  beta_o*MatrixXf::Identity(3*N,3*N);
	}
	
	void update_syn_pars () {
		R_sp = self.R.transpose()*peer.R;
		x_sp = self.R.transpose()*(peer.x-self.x);
		A_s.block(0*N,3*N,3*N,3*N) = -rho_u*kroneckerProduct(Snn,2.0*self.R*skewMat(x0_sp));
		A_s.block(3*N,0*N,3*N,3*N) =  A_s.block(0*N,3*N,3*N,3*N).transpose();
		b_s.head(3*N) =  rho_u*Sn*(4.0*(self.x-peer.x)+2.0*(self.R*x0_sp-peer.R*x0_ps));
		b_s.tail(3*N) = -rho_u*Sn*(2.0*skewMat(x0_sp)*x_sp+skewMat(x0_sp)*R_sp*x0_ps)+2.0*rho_o*Sn*skewVec(R0_sp*R_sp.transpose());
	}
	
	void set_opt_pars () {
		for (size_t i=0; i<18*N; ++i) {
			if (i<3*N)
				uof_ub[i] = vt;
			else if (i<6*N)
				uof_ub[i] = vr;
			else
				uof_ub[i] = fc;
			uof_lb[i] =-uof_ub[i];
		}
		for (size_t n=0; n<N; ++n) {
			for (size_t i=0; i<4; ++i)
				uof_lb[6*N+12*n+3*i+2] = epsilon;
		}
		set_tar_pars();
		set_vis_pars();
		set_cst_pars();
		A_obj.block(6*N,6*N,12*N,12*N) = (kappa_f+rho_f)*MatrixXf::Identity(12*N,12*N);
	}
	
	void update_opt_pars () {
		MatrixXf A = A_s;
		VectorXf b = b_s;
		update_tar_pars();
		if (mode!=0) {
			A += A_d;
			b += b_d;
		}
		update_vis_pars();
		if (mode!=1) {
			A += A_v;
			b += b_v;
		}
		A_obj.block(0*N,0*N, 6*N, 6*N) = A;
		b_obj.segment(0*N, 6*N) = b;
		b_obj.segment(6*N,12*N) = -rho_f*kroneckerProduct(MatrixXf::Ones(N,1),dt*MatrixXf(peer.f));
		update_cst_pars();
	}

};
