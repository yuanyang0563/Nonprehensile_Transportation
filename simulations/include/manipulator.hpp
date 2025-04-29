#include "common.h"

class manipulator {

  protected:
  	ros::NodeHandle nh;
	ros::Publisher pub_joint;
	ros::Subscriber sub_image;
	sensor_msgs::JointState joint_msg;
	vpImage<unsigned char> image;
	vpCameraParameters cam;
	vector<vector<vpImagePoint>> tagsCorners;
	vector<double*> data;
	stringstream file_name;
	bool file_init;

  public:
  	VectorXd a, alpha, d;
	VectorXd q, dq, f;
	VectorXd b_d, b_v, eta, lambda, h;
	MatrixXd A_d, A_v, J, Hu, Ho, Hf, Tu, To, L, Lu, Lo, Lm;
	Vector3d upsilon, omega;
	Vector3d xb, xeo, xec, xd, x, x0;
	Matrix3d Rb, Reo, Rec, Rd, R, R0;
	bool getJoints;
	bool getImage;
	vpDisplay *display;
	vpDetectorAprilTag detector;
	VectorXd zeta, zeta_d;
	
	manipulator (string name) {
		a = VectorXd(6);
		alpha = VectorXd(6);
		d = VectorXd(6);
		q = VectorXd(6);
		dq = VectorXd(6);
		f = VectorXd::Zero(12);
		upsilon = Vector3d::Zero();
		omega = Vector3d::Zero();
		J = MatrixXd(6,6);
		Hu = MatrixXd::Zero(6*N,3*N);
		Ho = MatrixXd::Zero(6*N,3*N);
		Hf = MatrixXd::Zero(6*N,12*N);
		eta = VectorXd::Zero(3*N);
		lambda = VectorXd::Zero(3*N);
		h = VectorXd::Zero(6*N);
		A_d = MatrixXd::Zero(6*N,6*N);
		b_d = VectorXd::Zero(6*N);
		A_v = MatrixXd::Zero(6*N,6*N);
		b_v = VectorXd::Zero(6*N);
		pub_joint = nh.advertise<sensor_msgs::JointState>(name+"/joint_position",1);
		joint_msg.position.resize(6);
		getJoints = false;
		sub_image = nh.subscribe(name+"/image",1,&manipulator::image_callback,this);
		image.resize(480,480);
		getImage = false;
		display = new vpDisplayX(image);
		vpDisplay::setTitle(image, name+"_image");
		cam.initPersProjWithoutDistortion(480,480,240,240);
		Tu = MatrixXd::Zero(6,3);
		To = MatrixXd::Zero(6,3);
		L = MatrixXd::Zero(8,6);
		Lu = MatrixXd::Zero(8,3);
		Lo = MatrixXd::Zero(8,3);
		Lm = MatrixXd::Zero(8,6);
		zeta = VectorXd(8);
		zeta_d = VectorXd(8);
		data.resize(32);
		for (int i=0; i<6; ++i)
			data[i] = &q(i);
		for (int i=6; i<12; ++i)
			data[i] = &dq(i-6);
		for (int i=12; i<20; ++i)
			data[i] = &zeta(i-12);
		for (int i=20; i<32; ++i)
			data[i] = &f(i-20);
		file_name << "../data/" << ros::this_node::getName() << "_" << name << "_" << time(0) << ".txt";
		file_init = false;
	}
	
	~manipulator () {
		delete display;
	}
	
	virtual void get_pose_jacobian () {
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
		if (!getJoints) {
			x0 = x;
			R0 = R;
			getJoints = !getJoints;
		}
		for (int i=0; i<6; ++i)
			J.col(i)<< Rb*z.col(i).cross(o.col(6)-o.col(i)), R.transpose()*Rb*z.col(i);
	}
	
	void image_callback(const sensor_msgs::Image::ConstPtr& msg) {
        	try {
            		for (int i=0; i<msg->height; ++i) {
                		for (int j=0; j<msg->width; ++j)
                			image[msg->height-i-1][j] = msg->data[i*msg->step+3*j];
            		}
            		detector.detect(image);
            		tagsCorners = detector.getTagsCorners();
            		if (tagsCorners[0].size()==4) {
				for (int i=0; i<4; ++i) {
            				vpPixelMeterConversion::convertPoint(cam,tagsCorners[0][i],zeta(2*i+0),zeta(2*i+1));
            				L.row(2*i+0) << -10.0,   0.0, 10.0*zeta(2*i+0), zeta(2*i+0)*zeta(2*i+1), -1.0-zeta(2*i+0)*zeta(2*i+0),  zeta(2*i+1);
            				L.row(2*i+1) <<   0.0, -10.0, 10.0*zeta(2*i+1), 1.0+zeta(2*i+1)*zeta(2*i+1), -zeta(2*i+0)*zeta(2*i+1), -zeta(2*i+0);
            			}
            			if (!getImage)
            				zeta_d = zeta;
            			getImage = true;
            		} else {
            			getImage = false;
            		}
			vpDisplay::display(image);
			vpDisplay::flush(image);
        	} catch (const vpException &e) {
	    		cerr << "Catch an exception: " << e.getMessage() << endl;
        	}
	}
	
	virtual void move_one_step () {
		dq << upsilon, omega;
		dq = J.transpose()*(J*J.transpose()).inverse()*dq;
		q += dt*dq;
		for (int i=0; i<6; i++)
			joint_msg.position[i] = q(i);
		pub_joint.publish(joint_msg);
		ros::spinOnce();
	}
	
	void set_tar_pars () {
		A_d.block(0*N,0*N,3*N,3*N) = kappa_u*kroneckerProduct(Snn,Matrix3d::Identity())+alpha_u*MatrixXd::Identity(3*N,3*N);
		A_d.block(3*N,3*N,3*N,3*N) = alpha_o*MatrixXd::Identity(3*N,3*N);
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
		Hu.block(0*N,0*N,3*N, 3*N) =  kroneckerProduct(Gamma,m*Matrix3d::Identity());
		Ho.block(3*N,0*N,3*N, 3*N) = -kroneckerProduct(Gamma,I*Reo.transpose());
		Hf.block(3*N,0*N,3*N,12*N) =  kroneckerProduct(MatrixXd::Identity(N,N),Go);
	}
	
	void update_cst_pars () {
		eta.head(3) = upsilon-R*xeo.cross(omega);
		lambda.head(3) = I*Reo.transpose()*omega;
		h.head(3*N) = m*eta-m*dt*kroneckerProduct(MatrixXd::Ones(N,1),R*skewMat(omega)*skewMat(omega)*xeo+g);
		h.tail(3*N) = lambda-dt*kroneckerProduct(MatrixXd::Ones(N,1),skewMat(Reo.transpose()*omega)*I*Reo.transpose()*omega);
		Ho.block(0*N,0*N,3*N, 3*N) = kroneckerProduct(Gamma,m*R*skewMat(xeo));
		Hf.block(0*N,0*N,3*N,12*N) = kroneckerProduct(MatrixXd::Identity(N,N),R*Reo*Gu);
	}
	
	inline double cost () {
		return (xd-x).norm()-(Rd*R.transpose()).trace()+3.0;
	}
	
	void store_data (double t_duration) {
		if (getImage) {
			ofstream data_stream;
			if (!file_init) {
				data_stream.open(file_name.str());
				file_init = !file_init;
			} else
				data_stream.open(file_name.str(),ios_base::app);
			data_stream << setiosflags(ios::fixed) << setprecision(2) << ros::Time::now().toSec();
			data_stream << ", " << setprecision(3) << t_duration;
			vector<double*>::iterator it;
			for (it=data.begin(); it!=data.end(); ++it)
				data_stream << ", " << setprecision(3) << **it;
			data_stream << endl;
			data_stream.close();
		}
	}

};

