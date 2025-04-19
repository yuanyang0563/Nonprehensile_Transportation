#include <ros/ros.h>
#include <std_msgs/Float64MultiArray.h>
#include <visp3/core/vpConfig.h>
#include <visp3/core/vpDisplay.h>
#include <visp3/gui/vpDisplayX.h>
#include <visp3/sensor/vpRealSense2.h>
#include <visp3/core/vpImageConvert.h>
#include <visp3/core/vpCameraParameters.h>
#include <visp3/core/vpPixelMeterConversion.h>
#include <visp3/detection/vpDetectorAprilTag.h>

using namespace std;

class camera {

  private:
  	ros::NodeHandle nh;
  	ros::Publisher pub;
  	ros::Rate* loop_rate;
  	std_msgs::Float64MultiArray msg;
  	vpDisplay *display;
  	vpCameraParameters cam;
  	vpImage<vpRGBa> image_color;
  	vpImage<unsigned char> image_grey;
  	vector<vector<vpImagePoint>> current_corners, initial_corners;
  	vpDetectorAprilTag detector;
  	vpRealSense2 g;
  	rs2::config config;
  	bool getImage;
  	
  public:
  	camera (string name) {
  		pub = nh.advertise<std_msgs::Float64MultiArray>(name+"/image_features",1);
  		loop_rate = new ros::Rate(30);
  		msg.data.resize(8);
  		image_color.resize(1080,1920);
  		display = new vpDisplayX(image_color);
  		vpDisplay::setTitle(image_color,name+"_image");
  		config.disable_stream(RS2_STREAM_DEPTH);
  		config.disable_stream(RS2_STREAM_INFRARED);
  		config.enable_stream(RS2_STREAM_COLOR,1920,1080,RS2_FORMAT_RGBA8,30);
  		g.open(config);
  		g.acquire(image_color);
  		cam = g.getCameraParameters(RS2_STREAM_COLOR,vpCameraParameters::perspectiveProjWithDistortion);
  		getImage = false;
  	}
  	
  	~camera() {
  		delete display;
  	}
  	
  	void getFeature () {
  		while (ros::ok()) {
  			g.acquire(image_color);
  			vpDisplay::display(image_color);
  			vpImageConvert::convert(image_color,image_grey);
  			bool status = detector.detect(image_grey);
  			if (status) {
  				current_corners = detector.getTagsCorners();
  				if (!getImage) {
  					initial_corners = current_corners;
  					getImage = !getImage;
  				} 	
  				for (int i=0; i<4; ++i) {
  					vpDisplay::displayCross(image_color,initial_corners[0][i],15,vpColor::blue,3);
  					vpDisplay::displayCross(image_color,current_corners[0][i],15,vpColor::red,3);
  					vpPixelMeterConversion::convertPoint(cam,current_corners[0][i],msg.data[2*i+0],msg.data[2*i+1]);
  				}
  				pub.publish(msg);
  				
  			}
  			vpDisplay::flush(image_color);
  			if (vpDisplay::getClick(image_color,false))
  				break;
  			ros::spinOnce();
  			loop_rate->sleep();
  		}
  	}

};
