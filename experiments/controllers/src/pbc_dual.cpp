#include "pbc_dual.hpp"

int main(int argc, char *argv[])
{
	init();
	// set ROS
	ros::init(argc,argv,"pbc_dual");
	// set control mode
	float mode;
	if (argc<2) {
		cout << "Please set the control mode:" << endl;
		cout << "0 for visual servoing." << endl;
		cout << "1 for target reaching." << endl;
		cout << "2 for combined control." << endl;
		return 1;
	} else {
		mode = stod(argv[1]);
		if (mode!=0.0 && mode!=1.0 && mode!=2.0) {
			ROS_ERROR("Wrong mode.");
			return 1;
		}
	}
	// create the arms object
	manipulator_dual arms("arm_l","arm_r",mode);
	// set the DH parameters of arms
	arms.left.a << 0.0, -0.24365, -0.21325, 0.0, 0.0, 0.0;
	arms.left.alpha << M_PI/2.0, 0.0, 0.0, M_PI/2.0, -M_PI/2.0, 0.0;
	arms.left.d << 0.15190, 0.0, 0.0, 0.11235, 0.08535, 0.0819;
	arms.right.a << 0.0, -0.24355, -0.21320, 0.0, 0.0, 0.0;
	arms.right.alpha << M_PI/2.0, 0.0, 0.0, M_PI/2.0, -M_PI/2.0, 0.0;
	arms.right.d << 0.15185, 0.0, 0.0, 0.13105, 0.08535, 0.0921;
	// set the poses of camera frames relative to tne end-effector frames of arms
	arms.left.xec << 0.035, 0.050, -0.075;
	arms.left.Rec << 0.0, -1.0, 0.0, 0.0, 0.0, -1.0, 1.0, 0.0, 0.0;
	arms.right.xec << -0.035, -0.050, -0.075;
	arms.right.Rec << 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0;
	// set the poses of arm bases
	arms.left.xb << 0.0, -0.44, 0.0;
	arms.left.Rb << cos(3.0*M_PI/4.0), -sin(3.0*M_PI/4.0), 0.0, sin(3.0*M_PI/4.0), cos(3.0*M_PI/4.0), 0.0, 0.0, 0.0, 1.0;
	arms.right.xb << 0.0, 0.44, 0.0;
	arms.right.Rb << cos(M_PI/4.0), -sin(M_PI/4.0), 0.0, sin(M_PI/4.0), cos(M_PI/4.0), 0.0, 0.0, 0.0, 1.0;
	// set the desired poses of arms
	arms.left.xd << 0.1392, -0.1616, 0.15;
	arms.left.Rd << cos(-M_PI/4.0), sin(-M_PI/4.0), 0.0, sin(-M_PI/4.0), -cos(-M_PI/4.0), 0.0, 0.0, 0.0, -1.0;
	arms.right.xd << 0.4008, 0.1000, 0.15;
	arms.right.Rd << cos(-M_PI/4.0), sin(-M_PI/4.0), 0.0, sin(-M_PI/4.0), -cos(-M_PI/4.0), 0.0, 0.0, 0.0, -1.0;
	// set parameters for visual servoing
	arms.set_vis_pars();
	while (ros::ok()) {
		// record the starting time of each round
		auto t_start = chrono::high_resolution_clock::now();
		// get the current robot poses and Jacobians
		arms.get_pose_jacobian();
		// update parameters for visual servoing
		arms.update_vis_pars();
		if (arms.getJoints && arms.getFeatures)
			arms.get_vel_inputs();
		arms.move_one_step();
		// count the time spent and store data
		auto t_stop = chrono::high_resolution_clock::now();
    		auto t_duration = chrono::duration<float>(t_stop-t_start);
    		arms.store_data(t_duration.count());
    		// fix the time duration of each round
    		if (t_duration.count()<dt)
    			this_thread::sleep_for(chrono::duration<float>(dt-t_duration.count()));
	}
	arms.stop_moving();
	
	return 0;
}

