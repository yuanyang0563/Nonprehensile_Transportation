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
  	vector<vector<vpImagePoint>> tagsCorners;
  	vpDetectorAprilTag detector;
  	vpRealSense2 g;
  	rs2::config config;
  	
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
  	}
  	
  	~camera() {
  		delete display;
  	}
  	
  	void getFeature () {
  		while (ros::ok()) {
  			g.acquire(image_color);
  			vpImageConvert::convert(image_color,image_grey);
  			detector.detect(image_grey);
  			tagsCorners = detector.getTagsCorners();
  			for (int i=0; i<4; ++i)
  				vpPixelMeterConversion::convertPoint(cam,tagsCorners[0][i],msg.data[2*i+0],msg.data[2*i+1]);
  			pub.publish(msg);
  			vpDisplay::display(image_color);
  			vpDisplay::flush(image_color);
  			if (vpDisplay::getClick(image_color,false))
  				break;
  			ros::spinOnce();
  			loop_rate->sleep();
  		}
  	}

};
