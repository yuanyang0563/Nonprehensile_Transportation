#include "manipulator.hpp"

struct Arm {
	Vector3f x, x0;
	Matrix3f R, R0;
};

class manipulator_side : public manipulator {

  public:
  	ros::Publisher pub_peer;
  	ros::Subscriber sub_peer;
  	geometry_msgs::Pose pose_msg;
  	
  	Arm self, peer;
  	
  	float mode;
  	
  	bool getPeer;

	manipulator_side (string name_self, string name_peer, float control_mode) : manipulator (name_self), mode(control_mode) {
		pub_peer = nh.advertise<geometry_msgs::Pose>(name_self+"/pose",1);
		sub_peer = nh.subscribe(name_peer+"/pose",1,&manipulator_side::peer_callback,this);
		getPeer = false;
		if (mode==0.0)
			gamma_pbc = 5.0;
	}
	
	void get_pose_jacobian () override {
		manipulator::get_pose_jacobian();
		if (getJoints) {
			if (!getPeer) {
				self.x0 = x0;
				self.R0 = R0;
			}
			self.x = x;
			self.R = R;
		}
	}
	
	void move_one_step () override {
		pose_msg.position.x = x(0);
		pose_msg.position.y = x(1);
		pose_msg.position.z = x(2);
		Quaternionf quaternion(R);
		pose_msg.orientation.w = quaternion.w();
		pose_msg.orientation.x = quaternion.x();
		pose_msg.orientation.y = quaternion.y();
		pose_msg.orientation.z = quaternion.z();
		if (getJoints && getFeatures)
			pub_peer.publish(pose_msg);
		manipulator::move_one_step();
	}
	
	void peer_callback (const geometry_msgs::Pose::ConstPtr& msg) {
		peer.x << msg->position.x, msg->position.y, msg->position.z;
		Quaternionf quaternion (msg->orientation.w, msg->orientation.x, msg->orientation.y, msg->orientation.z);
		peer.R << quaternion.toRotationMatrix();
		if (!getPeer) {
			peer.x0 = peer.x;
			peer.R0 = peer.R;
			getPeer = !getPeer;
		}
	}
	
	void get_vel_input () {
		Vector3f upsilon_d = kappa_pbc*(xd-x);
		Vector3f omega_d = kappa_pbc*R.transpose()*skewVec(Rd*R.transpose());
		Vector3f upsilon_s = rho_pbc*(peer.x-self.x-0.5*(self.R*self.R0.transpose()+peer.R*peer.R0.transpose())*(peer.x0-self.x0));
		Vector3f omega_s = rho_pbc*self.R.transpose()*skewVec(peer.R*peer.R0.transpose()*self.R0*self.R.transpose());
		VectorXf twist_v = gamma_pbc*(Lm.transpose()*Lm).inverse()*Lm.transpose()*(zeta_d-zeta);
		upsilon = upsilon_s;
		omega = omega_s;
		if (mode!=0) {
			upsilon += upsilon_d;
			omega += omega_d;
		}
		if (mode!=1) {
			upsilon += twist_v.head(3);
			omega += twist_v.tail(3);
		}
		upsilon = upsilon.array().min(vt).max(-vt);
		omega = omega.array().min(vr).max(-vr);
	}

};
