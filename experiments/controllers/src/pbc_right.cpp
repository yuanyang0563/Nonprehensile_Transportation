#include "pbc_side.hpp"

int main(int argc, char *argv[])
{
	init();
	// set ROS
	ros::init(argc,argv,"pbc_right");
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
	// create the arm object
	manipulator_side arm("arm_r","arm_l",mode);
	// set the DH parameters of the arm
	arm.a << 0.0, -0.24355, -0.21320, 0.0, 0.0, 0.0;
	arm.alpha << M_PI/2.0, 0.0, 0.0, M_PI/2.0, -M_PI/2.0, 0.0;
	arm.d << 0.15185, 0.0, 0.0, 0.13105, 0.08535, 0.0921;
	// set the pose of the camera frame relative to tne end-effector frame of the arm
	arm.xec << -0.035, -0.050, -0.075;
	arm.Rec << 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0;
	// set the pose of the arm base
	arm.xb << 0.0, 0.44, 0.0;
	arm.Rb << cos(M_PI/4.0), -sin(M_PI/4.0), 0.0, sin(M_PI/4.0), cos(M_PI/4.0), 0.0, 0.0, 0.0, 1.0;
	// set the desired pose of the arm
	arm.xd << 0.4008, 0.1000, 0.15;
	arm.Rd << cos(-M_PI/4.0), sin(-M_PI/4.0), 0.0, sin(-M_PI/4.0), -cos(-M_PI/4.0), 0.0, 0.0, 0.0, -1.0;
	// set parameters for visual servoing
	arm.set_vis_pars();
	while (ros::ok()) {
		// record the starting time of each round
		auto t_start = chrono::high_resolution_clock::now();
		// get the current robot pose and Jacobian
		arm.get_pose_jacobian();
		// update parameters for visual servoing
		arm.update_vis_pars();
		if (arm.getJoints && arm.getFeatures && arm.getPeer)
			arm.get_vel_input();
		arm.move_one_step();
		// count the time spent and store data
		auto t_stop = chrono::high_resolution_clock::now();
    		auto t_duration = chrono::duration<float>(t_stop-t_start);
    		arm.store_data(t_duration.count());
    		// fix the time duration of each round
    		if (t_duration.count()<dt)
    			this_thread::sleep_for(chrono::duration<float>(dt-t_duration.count()));
	}
	arm.stop_moving();
	
	return 0;
}

