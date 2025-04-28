#include <gurobi_c++.h>
#include "manipulator_side.hpp"

using namespace std;
using namespace Eigen;

int main(int argc, char *argv[])
{
	init();
	// set ROS
	ros::init(argc,argv,"mpc_left");
	// set control mode
	double mode;
	if (argc<2) {
		cout << "Please set the control mode:" << endl;
		cout << "0 for visual servoing." << endl;
		cout << "1 for target reaching." << endl;
		cout << "2 for combined control." << endl;
		return 1;
	} else {
		mode = stod(argv[1]);
		if (mode==0.0 || mode==1.0 || mode==2.0)
			ros::param::set("/control_mode",mode);
		else {
			ROS_ERROR("Wrong mode.");
			return 1;
		}
	}
	// create the arm object
	manipulator_side arm("arm_l","arm_r",mode);
	// set the DH parameters of the arm
	arm.a << 0.0, -0.24365, -0.21325, 0.0, 0.0, 0.0;
	arm.alpha << M_PI/2.0, 0.0, 0.0, M_PI/2.0, -M_PI/2.0, 0.0;
	arm.d << 0.1519, 0.0, 0.0, 0.11235, 0.08535, 0.0819;
	// set the initial joint coordinates of the arm
	arm.q << -M_PI, -1.0*M_PI/3.0, M_PI/4.0, -5.0*M_PI/12.0, -M_PI/2.0,  M_PI/4.0;
	// set the object pose relative to the end-effector frame of the arm
	arm.xeo << 0.0, -0.2405, 0.025;
	arm.Reo << 1.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, -1.0;
	// set the pose of camera frame relative to tne end-effector frame of the arm
	arm.xec << 0.045, -0.02, 0.01;
	arm.Rec << 1.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 1.0, 0.0;
	// set the pose of the arm base
	arm.xb << 0, -0.44, 0.0;
	arm.Rb << cos(3.0*M_PI/4.0), -sin(3.0*M_PI/4.0), 0.0, sin(3.0*M_PI/4.0), cos(3.0*M_PI/4.0), 0.0, 0.0, 0.0, 1.0;
	// set the desired pose of the arm
	arm.xd << 0.35, -0.2275, 0.20;
	arm.Rd << 1.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, -1.0;
	// set parameters for Gurobi to solve the problem
	arm.set_opt_pars();
	// call Gurobi to perform optimal control
	GRBEnv env = GRBEnv();
	env.set(GRB_IntParam_OutputFlag, 0);
	while (ros::ok()) {
		// record the starting time of each round
		auto t_start = chrono::high_resolution_clock::now();
		// get the current robot pose and Jacobian
		arm.get_pose_jacobian();
		// get the pose error and stops the simulation when the error is small enough
		if (arm.cost()<1e-5) {
			cout << "Target pose reached!" << endl;
			break;
		}
		// update parameters for Gurobi to the solve the problem
		arm.update_opt_pars();
		if (arm.getImage && arm.getPeer) {
			try {
				// create a gurobi model and add optimization variables uof with lower and upper bounds
				GRBModel model = GRBModel(env);
				GRBVar *uof = model.addVars(arm.uof_lb,arm.uof_ub,NULL,NULL,NULL,18*N);
				// set the objective function
				GRBQuadExpr obj = 0.0;
				for (int i=0; i<18*N; ++i) {
					for (int j=0; j<18*N; ++j)
						obj += dt*uof[i]*arm.A_obj(i,j)*uof[j];
					obj += 2.0*arm.b_obj(i)*uof[i];
				}
				model.setObjective(obj);
				// add motion constraints on the transported object
				for (int i=0; i<6*N; ++i) {
					GRBLinExpr cst = 0.0;
					for (int j=0; j<3*N; ++j)
						cst += arm.Hu(i,j-0*N)*uof[j];
					for (int j=3*N; j<6*N; ++j)
						cst -= arm.Ho(i,j-3*N)*uof[j];
					for (int j=6*N; j<18*N; ++j) {
						cst -= dt*arm.Hf(i,j-6*N)*uof[j];
					}
					model.addConstr(cst==arm.h(i));
				}
				// add constraints on the contact forces between the object and the tray
				for (int n=0; n<N; ++n) {
					for (int i=0; i<4; ++i) {
						int ind = 6*N+12*n+3*i;
						GRBQuadExpr cstc = uof[ind+0]*uof[ind+0]+uof[ind+1]*uof[ind+1]-mu*mu*uof[ind+2]*uof[ind+2];
						model.addQConstr(cstc<=0.0);
					}
				}
				// solve and set the Cartesian velocity control inputs
				model.optimize();
				for (int i=0; i<3; ++i) {
					arm.upsilon(i) = uof[i].get(GRB_DoubleAttr_X);
					arm.omega(i) = uof[3*N+i].get(GRB_DoubleAttr_X);
				}
				for (int i=0; i<12; ++i)
					arm.f(i) = uof[6*N+i].get(GRB_DoubleAttr_X);
			} catch(GRBException e) {
				cout << "No solution. Retry." << endl;
				// reduce the prediction horizon (N) or increase the damping (alpha) in the objective.
			} catch(...) {
				cout << "Exception during optimization." << endl;
			}
		}
		arm.move_one_step();
		// count the time spent and store data
		auto t_stop = chrono::high_resolution_clock::now();
    		auto t_duration = chrono::duration<double>(t_stop-t_start);
    		arm.store_data(t_duration.count());
    		// fix the time duration of each round
    		if (t_duration.count()<dt)
    			this_thread::sleep_for(chrono::duration<double>(dt-t_duration.count()));
	}
	
	return 0;
}

