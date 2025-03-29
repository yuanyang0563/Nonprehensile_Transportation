#include <gurobi_c++.h>
#include "manipulator.hpp"
#include <chrono>
#include <thread>

using namespace std;
using namespace Eigen;

int main(int argc, char *argv[])
{
	init();
	// create the arm object
	manipulator arm;
	// set DH parameters
	arm.a << 0.0, -0.24365, -0.21325, 0.0, 0.0, 0.0;
	arm.alpha << M_PI/2.0, 0.0, 0.0, M_PI/2.0, -M_PI/2.0, 0.0;
	arm.d << 0.1519, 0.0, 0.0, 0.11235, 0.08535, 0.0819;
	// set the initial joint coordinates
	arm.q << -M_PI/4.0, -2.0*M_PI/3.0, -M_PI/4.0, -7.0*M_PI/12.0, M_PI/2.0, -M_PI/4.0;
	arm.xeo << 0.0, -0.125, 0.025;
	arm.Reo << -1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, -1.0;
	arm.xb << 0.0, 0.0, 0.0;
	arm.Rb << 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0;
	arm.xd << 0.4, 0.25, 0.20;
	arm.Rd << 0.0, -1.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, -1.0;
    	// set parameters for the objective functions
	arm.set_tar_pars();
	// set upper and lower bounds of linear and angular velocities
	double uof_ub[18*N], uof_lb[18*N];
	for (int i=0; i<18*N; ++i) {
		if (i<3*N)
			uof_ub[i] = vt;
		else if (i<6*N)
			uof_ub[i] = vr;
		else
			uof_ub[i] = fc;
		uof_lb[i] =-uof_ub[i];
	}
	// set ROS
	ros::init(argc,argv,"mpc_controller");
	ros::NodeHandle nh;
	ros::Publisher pub = nh.advertise<std_msgs::Float32MultiArray>("joint_position",1);
	//ros::Subscriber sub = nh.subscribe("/joint_position",1,&manipulator::position_callback,this);
	// call Gurobi to perform optimal control
	GRBEnv env = GRBEnv();
	env.set(GRB_IntParam_OutputFlag, 0);
	while (ros::ok()) {
		// record the starting time of each round
		auto t_start = chrono::high_resolution_clock::now();
		// get the current robot poses and Jacobians
		arm.get_pose_jacobian();
		// get the pose error and stops the simulation when the error is small enough
		if (arm.cost()<1e-5)
			break;
		// update parameters for the optimal control problem
		arm.update_tar_pars();
		try {
			// create a gurobi model and add optimization variables uof with lower and upper bounds
			GRBModel model = GRBModel(env);
			GRBVar *uof = model.addVars(uof_lb,uof_ub,NULL,NULL,NULL,18*N);
			// set the objective function
			GRBQuadExpr obj = 0.0;
			for (int i=0; i<3*N; ++i) {
				for (int j=0; j<3*N; ++j)
					obj += dt*uof[i]*arm.Au_d(i,j)*uof[j]+dt*uof[3*N+i]*arm.Ao_d(i,j)*uof[3*N+j];
				obj += 2.0*arm.bu_d(i)*uof[i]+2.0*arm.bo_d(i)*uof[3*N+i];
			}
			model.setObjective(obj);
			// add motion constraints on the transported object
			for (int i=0; i<6*N; ++i) {
				GRBLinExpr cst = 0.0;
				for (int j=0; j<3*N; ++j)
					cst += arm.Hu(i,j)*uof[j];
				for (int j=3*N; j<6*N; ++j)
					cst -= arm.Ho(i,j-3*N)*uof[j];
				for (int j=6*N; j<18*N; ++j)
					cst -= dt*arm.Hf(i,j-6*N)*uof[j];
				model.addConstr(cst==arm.h(i));
			}
			// add constraints on the contact forces between the object and the tray
			for (int n=0; n<N; ++n) {
				for (int i=0; i<4; ++i) {
					int ind = 6*N+12*n+3*i;
					GRBQuadExpr cstc = uof[ind+0]*uof[ind+0]+uof[ind+1]*uof[ind+1]-mu*mu*uof[ind+2]*uof[ind+2];
					model.addQConstr(cstc<=0.0);
					GRBLinExpr cstn = uof[ind+2];
					model.addConstr(cstn>=epsilon);
				}
			}
			// solve and set the Cartesian velocity control inputs
			model.optimize();
			arm.upsilon << uof[0].get(GRB_DoubleAttr_X), uof[1].get(GRB_DoubleAttr_X), uof[2].get(GRB_DoubleAttr_X);
			arm.omega << uof[3*N].get(GRB_DoubleAttr_X), uof[3*N+1].get(GRB_DoubleAttr_X), uof[3*N+2].get(GRB_DoubleAttr_X);
		} catch(GRBException e) {
			cout << "No solution. Retry." << endl;
			// reduce the prediction horizon (N) or increase the weight (alpha) in the objective also make sense.
		} catch(...) {
			cout << "Exception during optimization." << endl;
		}
		arm.move_one_step();
		// publish joint position
    		std_msgs::Float32MultiArray msg;
		msg.data.resize(6);
		for (int i=0; i<6; i++)
			msg.data[i] = arm.q(i);
		pub.publish(msg);
		ros::spinOnce();
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

