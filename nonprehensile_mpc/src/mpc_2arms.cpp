#include <gurobi_c++.h>
#include "manipulators.hpp"
#include <chrono>
#include <thread>

using namespace std;
using namespace Eigen;

int main(int argc, char *argv[])
{
	init();
	// create the arms object
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
	// set parameters for the objective functions
	arms.set_tar_pars();
	arms.set_syn_pars();
	// set upper and lower bounds of linear and angular velocities
	double uof_ub[24*N], uof_lb[24*N];
	for (int i=0; i<24*N; ++i) {
		if (i<6*N)
			uof_ub[i] = vt;
		else if (i<12*N)
			uof_ub[i] = vr;
		else
			uof_ub[i] = fc;
		uof_lb[i] =-uof_ub[i];
	}
	// set ROS
	ros::init(argc,argv,"pid_controller");
	ros::NodeHandle nh;
	ros::Publisher pub = nh.advertise<std_msgs::Float32MultiArray>("joint_position",1);
	//ros::Subscriber sub = nh.subscribe("/joint_position",1,&manipulator::position_callback,this);
	// call Gurobi to perform optimal control
	GRBEnv env = GRBEnv();
	env.set(GRB_IntParam_OutputFlag, 0);
	while (true) {
		// record the starting time of each round
		auto t_start = chrono::high_resolution_clock::now();
		// get the current robot poses and Jacobians
		arms.get_pose_jacobian();
		// get the pose error and stops the simulation when the error is small enough
		if (arms.cost()<1e-5)
			break;
		// update parameters for the optimal control problem
		arms.update_tar_pars();
		arms.update_syn_pars();
		try {
			// create a gurobi model and add optimization variables uof with lower and upper bounds
			GRBModel model = GRBModel(env);
			GRBVar *uof = model.addVars(uof_lb,uof_ub,NULL,NULL,NULL,24*N);
			// set the objective function
			GRBQuadExpr obj = 0.0;
			for (int i=0; i<6*N; ++i) {
				for (int j=0; j<6*N; ++j)
					obj += dt*uof[i]*(arms.Au_d(i,j)+arms.Au_s(i,j))*uof[j]+dt*uof[6*N+i]*(arms.Ao_d(i,j)+arms.Ao_s(i,j))*uof[6*N+j];
				obj += 2.0*(arms.bu_d(i)+arms.bu_s(i))*uof[i]+2.0*(arms.bo_d(i)+arms.bo_s(i))*uof[6*N+i];
			}
			model.setObjective(obj);
			// add motion constraints on the transported object
			for (int i=0; i<6*N; ++i) {
				GRBLinExpr cst = 0.0;
				for (int j=0; j<3*N; ++j)
					cst += arms.left.Hu(i,j-0*N)*uof[j];
				for (int j=3*N; j<6*N; ++j)
					cst += arms.right.Hu(i,j-3*N)*uof[j];
				for (int j=6*N; j<9*N; ++j)
					cst -= arms.left.Ho(i,j-6*N)*uof[j];
				for (int j=9*N; j<12*N; ++j)
					cst -= arms.right.Ho(i,j-9*N)*uof[j];
				for (int j=12*N; j<24*N; ++j) {
					cst -= dt*arms.left.Hf(i,j-12*N)*uof[j];
					cst -= dt*arms.right.Hf(i,j-12*N)*uof[j];
				}
				model.addConstr(cst==arms.left.h(i)+arms.right.h(i));
			}
			// add constraints on the contact forces between the object and the tray
			for (int n=0; n<N; ++n) {
				for (int i=0; i<4; ++i) {
					int ind = 12*N+12*n+3*i;
					GRBQuadExpr cstc = uof[ind+0]*uof[ind+0]+uof[ind+1]*uof[ind+1]-mu*mu*uof[ind+2]*uof[ind+2];
					model.addQConstr(cstc<=0.0);
					GRBLinExpr cstn = uof[ind+2];
					model.addConstr(cstn>=epsilon);
				}
			}
			// solve and set the Cartesian velocity control inputs
			model.optimize();
			arms.left.upsilon << uof[0].get(GRB_DoubleAttr_X), uof[1].get(GRB_DoubleAttr_X), uof[2].get(GRB_DoubleAttr_X);
			arms.right.upsilon << uof[3*N].get(GRB_DoubleAttr_X), uof[3*N+1].get(GRB_DoubleAttr_X), uof[3*N+2].get(GRB_DoubleAttr_X);
			arms.left.omega << uof[6*N].get(GRB_DoubleAttr_X), uof[6*N+1].get(GRB_DoubleAttr_X), uof[6*N+2].get(GRB_DoubleAttr_X);
			arms.right.omega << uof[9*N].get(GRB_DoubleAttr_X), uof[9*N+1].get(GRB_DoubleAttr_X), uof[9*N+2].get(GRB_DoubleAttr_X);
		} catch(GRBException e) {
			cout << "No solution. Retry." << endl;
			// reduce the prediction horizon (N) or increase the damping (alpha) in the objective.
		} catch(...) {
			cout << "Exception during optimization." << endl;
		}
		arms.move_one_step();
		// publish joint position
    		std_msgs::Float32MultiArray msg;
		msg.data.resize(12);
		for (int i=0; i<6; ++i) {
			msg.data[i] = arms.left.q(i);
			msg.data[i+6] = arms.right.q(i);
		}
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

