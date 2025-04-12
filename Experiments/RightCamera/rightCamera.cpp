#include "camera.hpp"

using namespace std;

int main (int argc, char **argv)
{
	ros::init(argc,argv,"cam_r");
	camera realsense("arm_r");
	try {
		realsense.getImage();
	} catch (const vpException &e) {
		cerr << "Catch an exception: " << e.getMessage() << endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}
