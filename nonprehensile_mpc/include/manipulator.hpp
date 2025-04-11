#include "common.h"

using namespace std;
using namespace Eigen;

class manipulator {

  private:
  	ros::NodeHandle nh;
	ros::Publisher pub;
	ros::Subscriber sub;
	sensor_msgs::JointState msg;
	vpImage<unsigned char> image;
	vpCameraParameters cam;
	vector<vector<vpImagePoint>> tagsCorners;

  public:
  	VectorXd a, alpha, d;
	VectorXd q, dq;
	VectorXd b_d, b_v, eta, lambda, h;
	MatrixXd A_d, A_v, J, Hu, Ho, Hf, Tu, To, L, Lu, Lo, Lm;
	Vector3d upsilon, omega;
	Vector3d x, xb, xeo, xec, xd, x0;
	Matrix3d R, Rb, Reo, Rec, Rd, R0;
	bool getImage;
	vpDisplay *display;
	vpDetectorAprilTag detector;
	VectorXd zeta, zeta_d;
	
	manipulator(string arm) {
		a = VectorXd(6);
		alpha = VectorXd(6);
		d = VectorXd(6);
		q = VectorXd(6);
		dq = VectorXd(6);
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
		pub = nh.advertise<sensor_msgs::JointState>(arm+"/joint_position",1);
		msg.position.resize(6);
		sub = nh.subscribe(arm+"/image",1,&manipulator::image_callback,this);
		image.resize(480,480);
		getImage = false;
		display = new vpDisplayX(image);
		vpDisplay::setTitle(image, arm+"_image");
		cam.initPersProjWithoutDistortion(480,480,240,240);
		Tu = MatrixXd::Zero(6,3);
		To = MatrixXd::Zero(6,3);
		L = MatrixXd::Zero(8,6);
		Lu = MatrixXd::Zero(8,3);
		Lo = MatrixXd::Zero(8,3);
		Lm = MatrixXd::Zero(8,6);
		zeta = VectorXd(8);
		zeta_d = VectorXd(8);
	}
	
	~manipulator() {
		delete display;
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
	
	void image_callback(const sensor_msgs::ImageConstPtr& msg) {
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
	
	void move_one_step () {
		dq << upsilon, omega;
		dq = J.transpose()*(J*J.transpose()).inverse()*dq;
		q += dt*dq;
		for (int i=0; i<6; i++)
			msg.position[i] = q(i);
		pub.publish(msg);
		ros::spinOnce();
	}
	
	void set_tar_pars () {
		A_d.block(0*N,0*N,3*N,3*N) = kappa_u*kroneckerProduct(Snn,Matrix3d::Identity())+alpha_u*MatrixXd::Identity(3*N,3*N);
		A_d.block(3*N,3*N,3*N,3*N) = alpha_o*MatrixXd::Identity(3*N,3*N);
		Hu.block(0,0,3*N,3*N) = m*kroneckerProduct(Gamma,Matrix3d::Identity());
	}
	
	void update_tar_pars () {
		b_d.head(3*N) = kappa_u*Sn*(x-xd);
		b_d.tail(3*N) = kappa_o*Sn*skewVec(Rd.transpose()*R);
		eta.head(3) = upsilon-R*xeo.cross(omega);
		lambda.head(3) = I*Reo.transpose()*omega;
		h.head(3*N) = m*eta-dt*m*kroneckerProduct(MatrixXd::Ones(N,1),R*skewMat(omega)*skewMat(omega)*xeo+g);
		h.tail(3*N) = lambda-dt*kroneckerProduct(MatrixXd::Ones(N,1),skewMat(Reo.transpose()*omega)*I*Reo.transpose()*omega);
		Ho << m*kroneckerProduct(Gamma,R*skewMat(xeo)), -kroneckerProduct(Gamma,I*Reo.transpose());
		Hf << kroneckerProduct(MatrixXd::Identity(N,N),R*Reo*Gu), kroneckerProduct(MatrixXd::Identity(N,N),Go);
	}
	
	void set_vis_pars () {
		To << -Rec.transpose()*skewMat(xec), Rec.transpose();
	}
	
	void update_vis_pars () {
		A_v.setZero();
		b_v.setZero();
		Tu.block(0,0,3,3) = (R*Rec).transpose();
		Lu = L*Tu;
		Lo = L*To;
		Lm << Lu, Lo;
		A_v.block(0*N,0*N,3*N,3*N) = 5.0*kroneckerProduct(Snn,Lu.transpose()*Lu);
		A_v.block(0*N,3*N,3*N,3*N) = 5.0*kroneckerProduct(Snn,Lu.transpose()*Lo);
		A_v.block(3*N,0*N,3*N,3*N) = 5.0*kroneckerProduct(Snn,Lo.transpose()*Lu);
		A_v.block(3*N,3*N,3*N,3*N) = 5.0*kroneckerProduct(Snn,Lo.transpose()*Lo);
		b_v.head(3*N) = Sn*Lu.transpose()*(zeta-zeta_d);
		b_v.tail(3*N) = Sn*Lo.transpose()*(zeta-zeta_d);
	}
	
	inline double cost () {
		return (xd-x).norm()-(Rd*R.transpose()).trace()+3.0;
	}

};

