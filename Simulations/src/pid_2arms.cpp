#include "manipulators.hpp"

using namespace std;
using namespace Eigen;

int main(int argc, char *argv[])
{
	init();
	// set ROS
	ros::init(argc,argv,"pid_controller");
	// create the arms object
	manipulators arms("arm_l","arm_r");
	// set DH parameters
	arms.left.a << 0.0, -0.24365, -0.21325, 0.0, 0.0, 0.0;
	arms.left.alpha << M_PI/2.0, 0.0, 0.0, M_PI/2.0, -M_PI/2.0, 0.0;
	arms.left.d << 0.1519, 0.0, 0.0, 0.11235, 0.08535, 0.0819;
	arms.right.a << 0.0, -0.24365, -0.21325, 0.0, 0.0, 0.0;
	arms.right.alpha << M_PI/2.0, 0.0, 0.0, M_PI/2.0, -M_PI/2.0, 0.0;
	arms.right.d << 0.1519, 0.0, 0.0, 0.11235, 0.08535, 0.0819;
	// set the initial joint coordinates
	arms.left.q  <<  -M_PI, -1.0*M_PI/3.0,  M_PI/4.0, -5.0*M_PI/12.0, -M_PI/2.0,  M_PI/4.0;
	arms.right.q << -M_PI, -2.0*M_PI/3.0, -M_PI/4.0, -7.0*M_PI/12.0,  M_PI/2.0, -M_PI/4.0;
	arms.left.xeo  << 0.0, -0.2405, 0.025;
	arms.left.Reo  <<  1.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, -1.0;
	arms.right.xeo << 0.0, -0.2405, 0.025;
	arms.right.Reo << -1.0, 0.0, 0.0, 0.0,  1.0, 0.0, 0.0, 0.0, -1.0;
	arms.left.xec  << 0.045, -0.02, 0.01;
	arms.left.Rec  << 1.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 1.0, 0.0;
	arms.right.xec << 0.045, -0.02, 0.01;
	arms.right.Rec << 1.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 1.0, 0.0;
	arms.left.xb  <<  0, -0.44, 0.0;
	arms.left.Rb  << cos(3.0*M_PI/4.0), -sin(3.0*M_PI/4.0), 0.0, sin(3.0*M_PI/4.0), cos(3.0*M_PI/4.0), 0.0, 0.0, 0.0, 1.0;
	arms.right.xb << 0, 0.44, 0.0;
	arms.right.Rb << cos(M_PI/4.0), -sin(M_PI/4.0), 0.0, sin(M_PI/4.0), cos(M_PI/4.0), 0.0, 0.0, 0.0, 1.0;
	arms.left.xd  <<  0.35, -0.2275, 0.20;
	arms.left.Rd  <<  1.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, -1.0;
	arms.right.xd <<  0.35,  0.2275, 0.20;
	arms.right.Rd << -1.0, 0.0, 0.0, 0.0,  1.0, 0.0, 0.0, 0.0, -1.0;
	double mode;
	if (argc<2) {
		cout << "Please set the control mode:" << endl;
		cout << "0 for visual servoing, 1 for target reaching." << endl;
		return 1;
	} else {
		mode = stod(argv[1]);
		if (mode==0.0 || mode==1.0)
			ros::param::set("/control_mode",mode);
		else {
			cout << "Wrong mode." << endl;
			return 1;
		}
	}
	arms.get_pose_jacobian();
	arms.left.x0 = arms.left.x;
	arms.left.R0 = arms.left.R;
	arms.right.x0 = arms.right.x;
	arms.right.R0 = arms.right.R;
	arms.left.set_vis_pars();
	arms.right.set_vis_pars();
	while (ros::ok()) {
		// record the starting time of each round
		auto t_start = chrono::high_resolution_clock::now();
		// get the current robot poses and Jacobians
		arms.get_pose_jacobian();
		// get the pose error and stops the simulation when the error is small enough
		if (arms.cost()<1e-5) {
			cout << "Target pose reached!" << endl;
			break;
		}
		// compute controls and drive the arms
		if (arms.left.getImage && arms.right.getImage) {
			arms.left.update_vis_pars();
			arms.right.update_vis_pars();
			VectorXd tar_ul  = 1.0*(arms.left.xd-arms.left.x);
			VectorXd tar_ur  = 1.0*(arms.right.xd-arms.right.x);
			VectorXd tar_ol  = 1.0*skewVec(arms.left.Rd*arms.left.R.transpose());
			VectorXd tar_or  = 1.0*skewVec(arms.right.Rd*arms.right.R.transpose());
			VectorXd syn_ul  = 5.0*(arms.right.x-arms.left.x-0.5*(arms.left.R*arms.left.R0.transpose()+arms.right.R*arms.right.R0.transpose())*(arms.right.x0-arms.left.x0));
			VectorXd syn_ur  = 5.0*(arms.left.x-arms.right.x-0.5*(arms.left.R*arms.left.R0.transpose()+arms.right.R*arms.right.R0.transpose())*(arms.left.x0-arms.right.x0));
			VectorXd syn_ol  = 5.0*skewVec(arms.right.R*arms.right.R0.transpose()*arms.left.R0*arms.left.R.transpose());
			VectorXd syn_or  = 5.0*skewVec(arms.left.R*arms.left.R0.transpose()*arms.right.R0*arms.right.R.transpose());
			VectorXd vis_vl  = 5.0*(arms.left.Lm.transpose()*arms.left.Lm).inverse()*arms.left.Lm.transpose()*(arms.left.zeta_d-arms.left.zeta);
			VectorXd vis_vr  = 5.0*(arms.right.Lm.transpose()*arms.right.Lm).inverse()*arms.right.Lm.transpose()*(arms.right.zeta_d-arms.right.zeta);
			arms.left.upsilon  = mode*tar_ul+(1.0-mode)*vis_vl.head(3)+syn_ul;
			arms.left.omega  = arms.left.R.transpose()*(mode*tar_ol+syn_ol)+(1.0-mode)*vis_vl.tail(3);
			arms.right.upsilon = mode*tar_ur+(1.0-mode)*vis_vr.head(3)+syn_ur;
			arms.right.omega = arms.right.R.transpose()*(mode*tar_or+syn_or)+(1.0-mode)*vis_vr.tail(3);
		}
		arms.move_one_step();
		// count the time spent in solving the control per round and the maximum time
		auto t_stop = chrono::high_resolution_clock::now();
    		auto t_duration = chrono::duration<double>(t_stop-t_start);
    		t_max = max(t_max,t_duration.count());
    		// fix the time duration of each round
    		if (t_duration.count()<dt)
    			this_thread::sleep_for(chrono::duration<double>(dt-t_duration.count()));
	}
	cout << "Maximum solution time: " << t_max << " s." << endl;
	
	return 0;
}

