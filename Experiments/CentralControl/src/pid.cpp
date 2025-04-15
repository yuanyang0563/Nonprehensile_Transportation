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
	arms.left.a  << 0.0, -0.24365, -0.21325, 0.0, 0.0, 0.0;
	arms.left.alpha  << M_PI/2.0, 0.0, 0.0, M_PI/2.0, -M_PI/2.0, 0.0;
	arms.left.d  << 0.15190, 0.0, 0.0, 0.11235, 0.08535, 0.0819;
	arms.right.a << 0.0, -0.24355, -0.21320, 0.0, 0.0, 0.0;
	arms.right.alpha << M_PI/2.0, 0.0, 0.0, M_PI/2.0, -M_PI/2.0, 0.0;
	arms.right.d << 0.15185, 0.0, 0.0, 0.13105, 0.08535, 0.0921;
	// set configuration parameters
	arms.left.xeo  <<  0.0, -0.20, -0.05;
	arms.left.Reo  <<  1.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, -1.0;
	arms.right.xeo <<  0.0,  0.20, -0.05;
	arms.right.Reo <<  1.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, -1.0;
	arms.left.xec  <<  0.035,  0.050, -0.075;
	arms.left.Rec  <<  0.0, -1.0, 0.0, 0.0, 0.0, -1.0, 1.0, 0.0, 0.0;
	arms.right.xec << -0.035, -0.050, -0.075;
	arms.right.Rec <<  0.0,  1.0, 0.0, 0.0, 0.0,  1.0, 1.0, 0.0, 0.0;
	arms.left.xb   <<  0.0, -0.44, 0.0;
	arms.left.Rb   <<  cos(3.0*M_PI/4.0), -sin(3.0*M_PI/4.0), 0.0, sin(3.0*M_PI/4.0), cos(3.0*M_PI/4.0), 0.0, 0.0, 0.0, 1.0;
	arms.right.xb  <<  0.0,  0.44, 0.0;
	arms.right.Rb  <<  cos(1.0*M_PI/4.0), -sin(1.0*M_PI/4.0), 0.0, sin(1.0*M_PI/4.0), cos(1.0*M_PI/4.0), 0.0, 0.0, 0.0, 1.0;
	arms.left.xd   <<  0.4308, -0.2616, 0.15;
	arms.left.Rd   <<  cos(M_PI/4.0), sin(M_PI/4.0), 0.0, sin(M_PI/4.0), -cos(M_PI/4.0), 0.0, 0.0, 0.0, -1.0;
	arms.right.xd  <<  0.1692,  0.0000, 0.15;
	arms.right.Rd  <<  cos(M_PI/4.0), sin(M_PI/4.0), 0.0, sin(M_PI/4.0), -cos(M_PI/4.0), 0.0, 0.0, 0.0, -1.0;
	double mode;
	if (argc<2) {
		cout << "Please set the control mode:" << endl;
		cout << "0 for visual servoing." << endl;
		cout << "1 for target reaching." << endl;
		cout << "2 for combined control." << endl;
		return 1;
	} else {
		mode = stod(argv[1]);
		if (mode!=0.0 && mode!=1.0 && mode!=2.0) {
			cout << "Wrong mode." << endl;
			return 1;
		}
	}
	arms.set_vis_pars();
	while (ros::ok()) {
		// record the starting time of each round
		auto t_start = chrono::high_resolution_clock::now();
		// get the current robot poses and Jacobians
		arms.get_pose_jacobian();
		// compute controls and drive the arms
		if (arms.getJoints) {
			if (!arms.getPaths) {
				arms.plan_paths(1);
			} else {
				if (arms.cost()<1e-4) {
					if (arms.waypoint==1) {
						cout << "Target pose reached!" << endl;
						break;
					}
					if (arms.left.getFeature && arms.right.getFeature)
						arms.set_waypoints();
				}
				arms.left.upsilon   = 2.5*(arms.right.x-arms.left.x+arms.left.R*arms.left.R0.transpose()*(arms.left.x0-arms.right.x0));
				arms.left.upsilon  += 2.5*(arms.right.x-arms.left.x+arms.right.R*arms.right.R0.transpose()*(arms.left.x0-arms.right.x0));
				arms.left.omega     = 2.5*arms.left.R.transpose()*skewVec(arms.right.R*arms.right.R0.transpose()*arms.left.R0*arms.left.R.transpose());
				arms.right.upsilon  = 2.5*(arms.left.x-arms.right.x+arms.left.R*arms.left.R0.transpose()*(arms.right.x0-arms.left.x0));
				arms.right.upsilon += 2.5*(arms.left.x-arms.right.x+arms.right.R*arms.right.R0.transpose()*(arms.right.x0-arms.left.x0));
				arms.right.omega    = 2.5*arms.right.R.transpose()*skewVec(arms.left.R*arms.left.R0.transpose()*arms.right.R0*arms.right.R.transpose());
				if (mode==0.0 || mode==2.0) {
					arms.update_vis_pars();
					VectorXd vis_vl = 0.5*(arms.left.Lm.transpose()*arms.left.Lm).inverse()*arms.left.Lm.transpose()*(arms.left.zeta_d-arms.left.zeta);
					VectorXd vis_vr = 0.5*(arms.right.Lm.transpose()*arms.right.Lm).inverse()*arms.right.Lm.transpose()*(arms.right.zeta_d-arms.right.zeta);
					arms.left.upsilon += vis_vl.head(3);
					arms.left.omega += vis_vl.tail(3);
					arms.right.upsilon += vis_vr.head(3);
					arms.right.omega += vis_vr.tail(3);
				}
				if (mode==1.0 || mode==2.0) {
					arms.left.upsilon  += 0.1*(arms.left.xt-arms.left.x);
					arms.left.omega  += 0.1*arms.left.R.transpose()*skewVec(arms.left.Rt*arms.left.R.transpose());
					arms.right.upsilon  += 0.1*(arms.right.xt-arms.right.x);
					arms.right.omega  += 0.1*arms.right.R.transpose()*skewVec(arms.right.Rt*arms.right.R.transpose());
				}
			}
		}
		arms.control_forward();
		// count the time spent in solving the control per round and the maximum time
		auto t_stop = chrono::high_resolution_clock::now();
    		auto t_duration = chrono::duration<double>(t_stop-t_start);
    		t_max = max(t_max,t_duration.count());
    		// fix the time duration of each round
    		if (t_duration.count()<dt)
    			this_thread::sleep_for(chrono::duration<double>(dt-t_duration.count()));
	}
	arms.stop_moving();
	cout << "Maximum solution time: " << t_max << " s." << endl;
	
	return 0;
}

