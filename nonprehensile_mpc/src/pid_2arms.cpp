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
	arms.left.q <<  -M_PI, -1.0*M_PI/3.0,  M_PI/4.0, -5.0*M_PI/12.0, -M_PI/2.0,  M_PI/4.0;
	arms.right.q << -M_PI, -2.0*M_PI/3.0, -M_PI/4.0, -7.0*M_PI/12.0,  M_PI/2.0, -M_PI/4.0;
	arms.left.xeo << 0.0, -0.2405, 0.025;
	arms.left.Reo <<  1.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, -1.0;
	arms.right.xeo << 0.0, -0.2405, 0.025;
	arms.right.Reo << -1.0, 0.0, 0.0, 0.0,  1.0, 0.0, 0.0, 0.0, -1.0;
	arms.left.xec << 0.045, -0.02, 0.01;
	arms.left.Rec << 1.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 1.0, 0.0;
	arms.right.xec << 0.045, -0.02, 0.01;
	arms.right.Rec << 1.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 1.0, 0.0;
	arms.left.xb <<  0, -0.44, 0.0;
	arms.left.Rb << cos(3.0*M_PI/4.0), -sin(3.0*M_PI/4.0), 0.0, sin(3.0*M_PI/4.0), cos(3.0*M_PI/4.0), 0.0, 0.0, 0.0, 1.0;
	arms.right.xb << 0, 0.44, 0.0;
	arms.right.Rb << cos(M_PI/4.0), -sin(M_PI/4.0), 0.0, sin(M_PI/4.0), cos(M_PI/4.0), 0.0, 0.0, 0.0, 1.0;
	arms.left.xd <<  0.35, -0.2275, 0.20;
	arms.left.Rd << 1.0,  0.0, 0.0,  0.0, -1.0, 0.0, 0.0, 0.0, -1.0;
	arms.right.xd << 0.35,  0.2275, 0.20;
	arms.right.Rd << -1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, -1.0;
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
	arms.left.set_vis_pars();
	arms.right.set_vis_pars();
	while (ros::ok()) {
		// record the starting time of each round
		auto t_start = chrono::high_resolution_clock::now();
		// get the current robot poses and Jacobians
		arms.get_pose_jacobian();
		// get the pose error and stops the simulation when the error is small enough
		if (arms.cost()<1e-5)
			break;
		// compute controls and drive the arms
		if (arms.left.getImage && arms.right.getImage) {
			arms.left.update_vis_pars();
			arms.right.update_vis_pars();
			MatrixXd L_left(8,6), L_right(8,6);
			for (int i=0; i<4; ++i) {
				L_left.row(2*i+0) = arms.left.L[i].row(0);
				L_left.row(2*i+1) = arms.left.L[i].row(1);
				L_right.row(2*i+0) = arms.right.L[i].row(0);
				L_right.row(2*i+1) = arms.right.L[i].row(1);
			}
			L_left = L_left*arms.left.Tv;
			L_right = L_right*arms.right.Tv;
			VectorXd Zeta_l(8), Zeta_ld(8), Zeta_r(8), Zeta_rd(8);
			Zeta_l << arms.left.zeta[0], arms.left.zeta[1], arms.left.zeta[2], arms.left.zeta[3];
			Zeta_ld << arms.left.zeta_d[0], arms.left.zeta_d[1], arms.left.zeta_d[2], arms.left.zeta_d[3];
			Zeta_r << arms.right.zeta[0], arms.right.zeta[1], arms.right.zeta[2], arms.right.zeta[3];
			Zeta_rd << arms.right.zeta_d[0], arms.right.zeta_d[1], arms.right.zeta_d[2], arms.right.zeta_d[3];
			VectorXd vel_left = 5.0*(L_left.transpose()*L_left).inverse()*L_left.transpose()*(Zeta_ld-Zeta_l);
			VectorXd vel_right = 5.0*(L_right.transpose()*L_right).inverse()*L_right.transpose()*(Zeta_rd-Zeta_r);
			arms.left.upsilon = mode*(arms.left.xd-arms.left.x)+(1.0-mode)*vel_left.head(3);
			arms.left.omega = mode*skewVec(arms.left.Rd*arms.left.R.transpose())+(1.0-mode)*vel_left.tail(3);
			arms.right.upsilon = mode*(arms.right.xd-arms.right.x)+(1.0-mode)*vel_right.head(3);
			arms.right.omega = mode*skewVec(arms.right.Rd*arms.right.R.transpose())+(1.0-mode)*vel_right.tail(3);
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
	cout << "Target pose reached!" << endl;
	cout << "Maximum solution time: " << t_max << " s." << endl;
	
	return 0;
}

