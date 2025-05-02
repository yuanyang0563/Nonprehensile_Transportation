#include "common.h"

class manipulator {

  protected:
  	ros::NodeHandle nh;
	ros::Publisher pub_jv;
	ros::Subscriber sub_js;
	ros::Subscriber sub_if;
	std_msgs::Float64MultiArray jv_msg;
	vector<float*> data;
	stringstream file_name;
	bool file_init;

  public:
  	VectorXf a, alpha, d;
	VectorXf q, dq, f;
	VectorXf b_d, b_v, eta, lambda, h;
	MatrixXf A_d, A_v, J, Hu, Ho, Hf, Tu, To, L, Lu, Lo, Lm;
	Vector3f upsilon, omega;
	Vector3f xb, xeo, xec, xd, x, x0;
	Matrix3f Rb, Reo, Rec, Rd, R, R0;
	bool getJoints, getFeatures;
	VectorXf zeta, zeta_d;
	MatrixXf zeta_h;
	
	manipulator (string name) {
		a = VectorXf(6);
		alpha = VectorXf(6);
		d = VectorXf(6);
		q = VectorXf(6);
		dq = VectorXf(6);
		f = VectorXf::Zero(12);
		upsilon = Vector3f::Zero();
		omega = Vector3f::Zero();
		J = MatrixXf(6,6);
		Hu = MatrixXf::Zero(6*N,3*N);
		Ho = MatrixXf::Zero(6*N,3*N);
		Hf = MatrixXf::Zero(6*N,12*N);
		eta = VectorXf::Zero(3*N);
		lambda = VectorXf::Zero(3*N);
		h = VectorXf::Zero(6*N);
		A_d = MatrixXf::Zero(6*N,6*N);
		b_d = VectorXf::Zero(6*N);
		A_v = MatrixXf::Zero(6*N,6*N);
		b_v = VectorXf::Zero(6*N);
		pub_jv = nh.advertise<std_msgs::Float64MultiArray>(name+"/joint_group_vel_controller/command",1);
		jv_msg.data.resize(6);
		sub_js = nh.subscribe(name+"/joint_states",1,&manipulator::joint_callback,this);
		getJoints = false;
		sub_if = nh.subscribe(name+"/image_features",1,&manipulator::feature_callback,this);
		getFeatures = false;
		Tu = MatrixXf::Zero(6,3);
		To = MatrixXf::Zero(6,3);
		L = MatrixXf::Zero(8,6);
		Lu = MatrixXf::Zero(8,3);
		Lo = MatrixXf::Zero(8,3);
		Lm = MatrixXf::Zero(8,6);
		zeta = VectorXf(8);
		zeta_d = VectorXf(8);
		zeta_h = MatrixXf(8,8);
		data.resize(32);
		for (size_t i=0; i<6; ++i)
			data[i] = &q(i);
		for (size_t i=6; i<12; ++i)
			data[i] = &dq(i-6);
		for (size_t i=12; i<20; ++i)
			data[i] = &zeta(i-12);
		for (size_t i=20; i<32; ++i)
			data[i] = &f(i-20);
		file_name << "../data/" << ros::this_node::getName() << "_" << name << "_" << time(0) << ".txt";
		file_init = false;
	}
	
	virtual void get_pose_jacobian () {
		Matrix4f A;
		MatrixXf T = Matrix4f::Identity();
		Matrix<float, 3, 7> o, z;
		o.col(0) << 0.0, 0.0, 0.0;
		z.col(0) << 0.0, 0.0, 1.0;
		for (size_t i=0; i<6; ++i) {
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
		for (size_t i=0; i<6; ++i)
			J.col(i)<< Rb*z.col(i).cross(o.col(6)-o.col(i)), R.transpose()*Rb*z.col(i);
	}
	
	void feature_callback(const std_msgs::Float64MultiArray::ConstPtr& msg) {
		zeta << msg->data[0], msg->data[1], msg->data[2], msg->data[3], msg->data[4], msg->data[5], msg->data[6], msg->data[7];
		if (getFeatures) {
			for (size_t i=7; i>0; --i)
				zeta_h.col(i) = zeta_h.col(i-1);
			zeta_h.col(0) = zeta;
		} else {
			zeta_d = zeta;
			for (size_t i=0; i<8; ++i)
				zeta_h.col(i) = zeta;
			getFeatures = !getFeatures;
		}
		zeta = 0.125*zeta_h.rowwise().sum();
            	for (size_t i=0; i<4; ++i) {
            		L.row(2*i+0) << -10.0,   0.0, 10.0*zeta(2*i+0), zeta(2*i+0)*zeta(2*i+1), -1.0-zeta(2*i+0)*zeta(2*i+0),  zeta(2*i+1);
            		L.row(2*i+1) <<   0.0, -10.0, 10.0*zeta(2*i+1), 1.0+zeta(2*i+1)*zeta(2*i+1), -zeta(2*i+0)*zeta(2*i+1), -zeta(2*i+0);
            	}
	}
	
	void joint_callback (const sensor_msgs::JointState::ConstPtr& msg) {
		q << msg->position[2], msg->position[1], msg->position[0], msg->position[3], msg->position[4], msg->position[5];
		if (!getJoints) {
			get_pose_jacobian();
			x0 = x;
			R0 = R;
			getJoints = !getJoints;
		}
	}
	
	virtual void move_one_step () {
		dq << upsilon, omega;
		dq = J.transpose()*(J*J.transpose()).inverse()*dq;
		for (size_t i=0; i<6; i++)
			jv_msg.data[i] = dq(i);
		pub_jv.publish(jv_msg);
		ros::spinOnce();
	}
	
	void set_tar_pars () {
		A_d.block(0*N,0*N,3*N,3*N) = kappa_u*kroneckerProduct(Snn,Matrix3f::Identity())+alpha_u*MatrixXf::Identity(3*N,3*N);
		A_d.block(3*N,3*N,3*N,3*N) = alpha_o*MatrixXf::Identity(3*N,3*N);
	}
	
	void update_tar_pars () {
		b_d.head(3*N) = kappa_u*Sn*(x-xd);
		b_d.tail(3*N) = kappa_o*Sn*skewVec(Rd.transpose()*R);
	}
	
	void set_vis_pars () {
		To << -Rec.transpose()*skewMat(xec), Rec.transpose();
	}
	
	void update_vis_pars () {
		Tu.block(0,0,3,3) = (R*Rec).transpose();
		Lu = L*Tu;
		Lo = L*To;
		Lm << Lu, Lo;
		A_v.block(0*N,0*N,3*N,3*N) = 5.0*gamma_v*kroneckerProduct(Snn,Lu.transpose()*Lu);
		A_v.block(0*N,3*N,3*N,3*N) = 5.0*gamma_v*kroneckerProduct(Snn,Lu.transpose()*Lo);
		A_v.block(3*N,0*N,3*N,3*N) = 5.0*gamma_v*kroneckerProduct(Snn,Lo.transpose()*Lu);
		A_v.block(3*N,3*N,3*N,3*N) = 5.0*gamma_v*kroneckerProduct(Snn,Lo.transpose()*Lo);
		b_v.head(3*N) = gamma_v*Sn*Lu.transpose()*(zeta-zeta_d);
		b_v.tail(3*N) = gamma_v*Sn*Lo.transpose()*(zeta-zeta_d);
	}
	
	void set_cst_pars () {
		Hu.block(0*N,0*N,3*N, 3*N) =  kroneckerProduct(Gamma,m*Matrix3f::Identity());
		Ho.block(3*N,0*N,3*N, 3*N) = -kroneckerProduct(Gamma,I*Reo.transpose());
		Hf.block(3*N,0*N,3*N,12*N) =  kroneckerProduct(MatrixXf::Identity(N,N),Go);
	}
	
	void update_cst_pars () {
		eta.head(3) = upsilon-R*xeo.cross(omega);
		lambda.head(3) = I*Reo.transpose()*omega;
		h.head(3*N) = m*eta-m*dt*kroneckerProduct(MatrixXf::Ones(N,1),R*skewMat(omega)*skewMat(omega)*xeo+g);
		h.tail(3*N) = lambda-dt*kroneckerProduct(MatrixXf::Ones(N,1),skewMat(Reo.transpose()*omega)*I*Reo.transpose()*omega);
		Ho.block(0*N,0*N,3*N, 3*N) = kroneckerProduct(Gamma,m*R*skewMat(xeo));
		Hf.block(0*N,0*N,3*N,12*N) = kroneckerProduct(MatrixXf::Identity(N,N),R*Reo*Gu);
	}
	
	void store_data (float t_duration) {
		if (getJoints && getFeatures) {
			ofstream data_stream;
			if (!file_init) {
				data_stream.open(file_name.str());
				file_init = !file_init;
			} else
				data_stream.open(file_name.str(),ios_base::app);
			data_stream << setiosflags(ios::fixed) << setprecision(2) << ros::Time::now().toSec();
			data_stream << ", " << setprecision(4) << t_duration;
			vector<float*>::iterator it;
			for (it=data.begin(); it!=data.end(); ++it)
				data_stream << ", " << setprecision(3) << **it;
			data_stream << endl;
			data_stream.close();
		}
	}

};

