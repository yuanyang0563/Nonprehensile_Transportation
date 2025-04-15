#include <gurobi_c++.h>
#include "manipulators.hpp"

using namespace std;
using namespace Eigen;

int main(int argc, char *argv[])
{
	init();
	// set ROS
	ros::init(argc,argv,"mpc_controller");
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
	arms.set_tar_pars();
	arms.set_vis_pars();
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
	// call Gurobi to perform optimal control
	GRBEnv env = GRBEnv();
	env.set(GRB_IntParam_OutputFlag, 0);
	while (ros::ok()) {
		// record the starting time of each round
		auto t_start = chrono::high_resolution_clock::now();
		// get the current robot poses and Jacobians
		arms.get_pose_jacobian();
		// update parameters for the optimal control problem
		if (arms.getJoints) {
			if (!arms.getPaths) {
				arms.plan_paths(1);
				arms.set_syn_pars();
			} else {
				if (arms.cost()<1e-4) {
					if (arms.waypoint==1) {
						cout << "Target pose reached!" << endl;
						break;
					}
					if (arms.left.getFeature && arms.right.getFeature)
						arms.set_waypoints();
				}
				arms.update_tar_pars();
				arms.update_syn_pars();
				arms.update_vis_pars();
				MatrixXd A_obj = arms.A_s;
				VectorXd b_obj = arms.b_s;
				if (mode==0.0 || mode==2.0) {
					A_obj += arms.A_v;
					b_obj += arms.b_v;
				}
				if (mode==1.0 || mode==2.0) {
					A_obj += arms.A_d;
					b_obj += arms.b_d;
				}
				try {
					// create a gurobi model and add optimization variables uof with lower and upper bounds
					GRBModel model = GRBModel(env);
					GRBVar *uof = model.addVars(uof_lb,uof_ub,NULL,NULL,NULL,24*N);
					// set the objective function
					GRBQuadExpr obj = 0.0;
					for (int i=0; i<12*N; ++i) {
						for (int j=0; j<12*N; ++j)
							obj += dt*uof[i]*A_obj(i,j)*uof[j];
						obj += 2.0*b_obj(i)*uof[i];
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
					// reduce the prediction horizon or increase the damping in the objective.
					arms.left.upsilon.setZero();
					arms.right.upsilon.setZero();
					arms.left.omega.setZero();
					arms.right.omega.setZero();
				} catch(...) {
					cout << "Exception during optimization." << endl;
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

