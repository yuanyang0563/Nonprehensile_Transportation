#include "manipulators.hpp"
#include "remote_api/RemoteAPIClient.h"
#include <chrono>
#include <thread>

using namespace std;
using namespace Eigen;

int main(int argc, char *argv[])
{
	init();
	manipulators arms;
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
	arms.left.xb <<  0, -0.44, 0.0;
	arms.left.Rb << cos(3.0*M_PI/4.0), -sin(3.0*M_PI/4.0), 0.0, sin(3.0*M_PI/4.0), cos(3.0*M_PI/4.0), 0.0, 0.0, 0.0, 1.0;
	arms.right.xb << 0, 0.44, 0.0;
	arms.right.Rb << cos(M_PI/4.0), -sin(M_PI/4.0), 0.0, sin(M_PI/4.0), cos(M_PI/4.0), 0.0, 0.0, 0.0, 1.0;
	arms.left.xd <<  0.35, -0.2275, 0.20;
	arms.left.Rd << 1.0,  0.0, 0.0,  0.0, -1.0, 0.0, 0.0, 0.0, -1.0;
	arms.right.xd << 0.35, 0.2275, 0.20;
	arms.right.Rd << -1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, -1.0;
	// initialize the remote API to interface with CoppeliaSim
	RemoteAPIClient client;
    	auto sim = client.getObject().sim();
    	for (int i=0; i<6; ++i) {
    		arms.left.jointHandles[i] = sim.getObject("/joint_l"+to_string(i+1));
    		arms.right.jointHandles[i] = sim.getObject("/joint_r"+to_string(i+1));
    		sim.setJointPosition(arms.left.jointHandles[i],arms.left.q(i));
    		sim.setJointPosition(arms.right.jointHandles[i],arms.right.q(i));
    	}
    	sim.startSimulation();
	while (true) {
		// record the starting time of each round
		auto t_start = chrono::high_resolution_clock::now();
		// get the current robot poses and Jacobians
		arms.get_pose_jacobian();
		// get the pose error and stops the simulation when the error is small enough
		if (arms.cost()<1e-5)
			break;
		arms.left.upsilon = arms.left.xd-arms.left.x;
		arms.left.omega = skewVec(arms.left.Rd*arms.left.R.transpose());
		arms.right.upsilon = arms.right.xd-arms.right.x;
		arms.right.omega = skewVec(arms.right.Rd*arms.right.R.transpose());
		arms.move_one_step();
		// count the time spent in solving the control per round and the maximum time
		auto t_stop = chrono::high_resolution_clock::now();
    		auto t_duration = chrono::duration<double>(t_stop-t_start);
    		t_max = max(t_max,t_duration.count());
    		// fix the time duration of each round
    		if (t_duration.count()<dt)
    			this_thread::sleep_for(chrono::duration<double>(dt-t_duration.count()));
    		// update the position of the simulated robot in CoppeliaSim
    		for (int i=0; i<6; ++i) {
			sim.setJointPosition(arms.left.jointHandles[i],arms.left.q(i));
			sim.setJointPosition(arms.right.jointHandles[i],arms.right.q(i));
		}
	}
	sim.stopSimulation();
	cout << "Target pose reached!" << endl;
	cout << "Maximum solution time: " << t_max << " s." << endl;
	
	return 0;
}

