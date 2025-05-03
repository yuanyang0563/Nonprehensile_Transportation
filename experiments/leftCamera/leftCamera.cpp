#include "camera.hpp"

int main (int argc, char **argv)
{
	ros::init(argc,argv,"cam_l");
	camera realsense("arm_l");
	try {
		realsense.getFeature();
	} catch (const vpException &e) {
		cerr << "Catch an exception: " << e.getMessage() << endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}
