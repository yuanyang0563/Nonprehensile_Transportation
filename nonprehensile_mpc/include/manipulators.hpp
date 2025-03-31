#include "manipulator.hpp"

using namespace std;
using namespace Eigen;

class manipulators {

  public:
  
	manipulator left;
	manipulator right;
	
	MatrixXd Au_d, Ao_d, Au_s, Ao_s;
	VectorXd bu_d, bo_d, bu_s, bo_s;
	
	manipulators (string arm_left, string arm_right) : left(arm_left), right(arm_right) {
		Au_d = MatrixXd(6*N,6*N);
		Ao_d = MatrixXd(6*N,6*N);
		bu_d = VectorXd(6*N);
		bo_d = VectorXd(6*N);
		Au_s = MatrixXd(6*N,6*N);
		Ao_s = MatrixXd(6*N,6*N);
		bu_s = VectorXd(6*N);
		bo_s = VectorXd(6*N);
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
		bu_d << left.bu_d, right.bu_d;
		bo_d << left.bo_d, right.bo_d;
	}
	
	void set_tar_pars () {
		Au_d.setZero();
		Ao_d.setZero();
		left.set_tar_pars();
		right.set_tar_pars();
		Au_d.block(0,0,3*N,3*N) = left.Au_d;
		Ao_d.block(0,0,3*N,3*N) = left.Ao_d;
		Au_d.block(3*N,3*N,3*N,3*N) = right.Au_d;
		Ao_d.block(3*N,3*N,3*N,3*N) = right.Ao_d;
	}
	
	void update_tar_pars () {
		left.update_tar_pars();
		right.update_tar_pars();
	}
	
	void set_syn_pars () {
		Au_s.block(0,0,3*N,3*N) = kroneckerProduct(Snn,Matrix3d::Identity());
		Au_s.block(0,3*N,3*N,3*N) = -kroneckerProduct(Snn,left.Rd*right.Rd.transpose());
		Au_s.block(3*N,0,3*N,3*N) = -kroneckerProduct(Snn,right.Rd*left.Rd.transpose());
		Au_s.block(3*N,3*N,3*N,3*N) = kroneckerProduct(Snn,Matrix3d::Identity());
		Au_s *= rho_u;
		Ao_s.setZero();
	}

	void update_syn_pars () {
		bu_s.head(3*N) = rho_u*Sn*(left.x-left.xd-left.Rd*right.Rd.transpose()*(right.x-right.xd));
		bu_s.tail(3*N) = rho_u*Sn*(right.x-right.xd-right.Rd*left.Rd.transpose()*(left.x-left.xd));
		bo_s.head(3*N) = rho_o*Sn*skewVec(right.R.transpose()*right.Rd*left.Rd.transpose()*left.R);
		bo_s.tail(3*N) = rho_o*Sn*skewVec(left.R.transpose()*left.Rd*right.Rd.transpose()*right.R);
		Ao_s.block(0,3*N,3*N,3*N) = 0.5*rho_o*kroneckerProduct(Snn,right.R.transpose()*right.Rd*left.Rd.transpose()*left.R-(right.R.transpose()*right.Rd*left.Rd.transpose()*left.R).trace()*Matrix3d::Identity());
		Ao_s.block(3*N,0,3*N,3*N) = 0.5*rho_o*kroneckerProduct(Snn,left.R.transpose()*left.Rd*right.Rd.transpose()*right.R-(left.R.transpose()*left.Rd*right.Rd.transpose()*right.R).trace()*Matrix3d::Identity());
	}
	
	inline double cost () {
		return left.cost()+right.cost();
	}

};





