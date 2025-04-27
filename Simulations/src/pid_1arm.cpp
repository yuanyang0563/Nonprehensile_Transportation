#include "manipulator.hpp"

using namespace std;
using namespace Eigen;

int main(int argc, char *argv[])
{
	init();
	// set ROS
	ros::init(argc,argv,"pid_controller");
	// create the arm object
	manipulator arm("ur3");
	// set DH parameters
	arm.a << 0.0, -0.24365, -0.21325, 0.0, 0.0, 0.0;
	arm.alpha << M_PI/2.0, 0.0, 0.0, M_PI/2.0, -M_PI/2.0, 0.0;
	arm.d << 0.1519, 0.0, 0.0, 0.11235, 0.08535, 0.0819;
	// set the initial joint coordinates
	arm.q << -M_PI/4.0, -2.0*M_PI/3.0, -M_PI/4.0, -7.0*M_PI/12.0, M_PI/2.0, -M_PI/4.0;
	arm.xec << 0.045, -0.02, 0.01;
	arm.Rec << 1.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 1.0, 0.0;
	arm.xb << 0.0, 0.0, 0.0;
	arm.Rb << 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0;
	arm.xd << 0.4, 0.25, 0.20;
	arm.Rd << 0.0, -1.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, -1.0;
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
	arm.set_vis_pars();
	while (ros::ok()) {
		// record the starting time of each round
		auto t_start = chrono::high_resolution_clock::now();
		// get the current robot poses and Jacobians
		arm.get_pose_jacobian();
		// get the pose error and stops the simulation when the error is small enough
		if (arm.cost()<1e-5) {
			cout << "Target pose reached!" << endl;
			break;
		}
		// compute control and drive the arm
		if (arm.getImage) {
			arm.update_vis_pars();
			VectorXd vel = 3.0*(arm.Lm.transpose()*arm.Lm).inverse()*arm.Lm.transpose()*(arm.zeta_d-arm.zeta);
			arm.upsilon = mode*(arm.xd-arm.x)+(1.0-mode)*vel.head(3);
			arm.omega = mode*arm.R.transpose()*skewVec(arm.Rd*arm.R.transpose())+(1.0-mode)*vel.tail(3);
		}
		arm.move_one_step();
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

