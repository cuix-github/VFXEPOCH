/*******************************************************************************
    VFXEPOCH - Physically based simulation VFX

    Copyright (c) 2016 Snow Tsui <trevor.miscellaneous@gmail.com>

    All rights reserved. Use of this source code is governed by
    the MIT license as written in the LICENSE file.
*******************************************************************************/

// This utility is usually used for testing grids within several digits.
// WARNING! Higher dimensions could violate the display format.

#include "Helpers.h"
#include "VisualizerHelpers.h"
#include "VisualizerHelpers_OpenGL.h"

#include <OpenEXR/ImfRgbaFile.h>
#include <OpenEXR/ImfStringAttribute.h>
#include <OpenEXR/ImfMatrixAttribute.h>
#include <OpenEXR/ImfArray.h>
#include <OpenEXR/ImfNamespace.h>

// For visualization
using namespace VFXEpoch::Gluvi;
using namespace VFXEpoch::OpenGL_Utility;
namespace IMF = OPENEXR_IMF_NAMESPACE;
using namespace IMF;
using namespace IMATH_NAMESPACE;

#include <string>
#include <sstream>
#include <cctype>
#include <algorithm>

/***************************** For Visualization ******************************/
unsigned int win_width = 720;
unsigned int win_height = 720;
int total_frames = 240;

PanZoom2D cam(-0.5, 0.5, -0.5, 0.5);
bool process_cmd_params(int argc, char* argv[]);
bool init_solver_params();
void display();
void mouse();
void drag(int x, int y);
void timer(int arg);

// Visualization & Output Data
bool preview = true;
bool outputParticles = false;
bool outputAlembic = false;
/***************************** For Visualization ******************************/

EulerGAS2D::Parameters params;
VFXEpoch::Solvers::Euler_Fluid2D_Base *solver = NULL;
VFXEpoch::Solvers::EulerGAS2D *gas_solver = NULL;
const float source = 1.0f;
VFXEpoch::Vector2Df c0(0.5,0.5), c1(0.7,0.5), c2(0.3,0.35), c3(0.5,0.7);
float rad0 = 0.4,  rad1 = 0.1,  rad2 = 0.1,   rad3 = 0.1;
const int tracers = 100000;

using namespace Helpers;

/****************************** For Dbuge ******************************/
float circle_phi(const VFXEpoch::Vector2Df& pos, const VFXEpoch::Vector2Df& center, float radius) {
	return (VFXEpoch::Dist2D(pos, center) - radius);	
}

float boundary_phi(const VFXEpoch::Vector2Df& position) {
   float phi0 = -circle_phi(position, c0, rad0);
   float phi1 = circle_phi(position, c1, rad1);
   float phi2 = circle_phi(position, c2, rad2);
   float phi3 = circle_phi(position, c3, rad3);
   return std::min(std::min(phi0,phi1), std::min(phi2,phi3));
   return phi0;
}
/****************************** For Dbuge ******************************/

//write tracer particle data to a file so that can be visulized by matlab
void write_file(double *pos_x, double *pos_y, int num, int frame)
{
	char filename[256];
	sprintf(filename,"../../outputs/sims/Particle_data%04d.bin",frame);
	float *data;
	data = new float[num*4];
	for (int i=0;i<num;i++)
	{
		data[i*4+0] = pos_x[i];
		data[i*4+1] = pos_y[i];
		data[i*4+2] = 0;
		data[i*4+3] = 1;
	}
	FILE *f;
	f = fopen(filename,"wb");
	fwrite(&num, sizeof(int), 1, f);
	fwrite(data,sizeof(float),4*num,f);

	delete[]data;
	fclose(f);
}

bool process_cmd_params(int argc, char* argv[])
{
	// Simulation Resolution is specified
	if (3 == argc) {
		std::string resx = argv[1];
		std::string resy = argv[2];
		std::string::size_type sz;
		params.dimension[0] = std::stoi(resx, &sz);
		params.dimension[1] = std::stoi(resy, &sz);
		std::cout << "Dimension: " <<  params.dimension[0] << " x " << params.dimension[1] << std::endl;
	}
	// Window Size is specified
	else if (5 == argc) {
		std::string resx = argv[1];
		std::string resy = argv[2];
		std::string window_width = argv[3];
		std::string window_height = argv[4];
		std::string::size_type sz;
		params.dimension[0] = std::stoi(resx, &sz);
		params.dimension[1] = std::stoi(resy, &sz);
		win_width = std::stoi(window_width, &sz);
		win_height = std::stoi(window_height, &sz);
	}
	// Total Frame is specified
	else if (6 == argc) {
		std::string resx = argv[1];
		std::string resy = argv[2];
		std::string window_width = argv[3];
		std::string window_height = argv[4];
		std::string frames = argv[5];
		std::string::size_type sz;
		params.dimension[0] = std::stoi(resx, &sz);
		params.dimension[1] = std::stoi(resy, &sz);
		win_width = std::stoi(window_width, &sz);
		win_height = std::stoi(window_height, &sz);
		total_frames = std::stoi(frames, &sz);
	}
	// Preview on/off
	else if (7 == argc) {
		std::string resx = argv[1];
		std::string resy = argv[2];
		std::string window_width = argv[3];
		std::string window_height = argv[4];
		std::string frames = argv[5];
		std::stringstream preview_switch(argv[6]);
		std::string::size_type sz;
		params.dimension[0] = std::stoi(resx, &sz);
		params.dimension[1] = std::stoi(resy, &sz);
		win_width = std::stoi(window_width, &sz);
		win_height = std::stoi(window_height, &sz);
		total_frames = std::stoi(frames, &sz);
		if (!(preview_switch >> std::boolalpha >> preview)) {
			return false;
		}
	}
	else if (1 == argc) {
		params.dimension[0] = 32;
		params.dimension[1] = 32;
		return true;
	}
	else
		return false;

	return true;
}

bool init_solver_params()
{
	params.h = 1.f / (float)params.dimension[0];
	params.dt = 0.005f;
	params.buoyancy_alpha = 0.1;
	params.buoyancy_beta = 0.3;
	params.vort_conf_eps = 0.55;
	params.density_source = 20;
	params.external_force_strength = 10;
	params.max_iterations = 300;
	params.min_tolerance = 1e-5;
	params.diff = 0.01;
	params.visc = 0.01;
	return true;
}

void display()
{
	glColor3f(1.0, 1.0, 1.0);
	OpenGL_Utility::draw_grid2d(VFXEpoch::Vector2Dd(-0.5, -0.5), params.h, params.dimension[0], params.dimension[1]);
}

void mouse(int button, int state, int x, int y)
{
	if(state==GLUT_UP)
		std::cout << "Mouse Clicked" << std::endl;
    else if(button==GLUT_LEFT_BUTTON)
        std::cout << "Mouse Left Button Clicked" << std::endl;
    else if(button==GLUT_MIDDLE_BUTTON)
        std::cout << "Mouse Middle Button Clicked" << std::endl;
    else if(button==GLUT_RIGHT_BUTTON)
        std::cout << "Mouse Right Button Clicked" << std::endl;
}

void timer(int arg)
{
	if (arg > 0) {
		int i = (total_frames - arg) + 1;
		cout << "****************** Frame " << i << " ******************" << endl;
		gas_solver->step();
		cout << "**************** Step " << i << " done ****************" << endl;
		cout << endl;
		glutPostRedisplay();
		glutTimerFunc(1000 / 24, timer, arg - 1);
	}
	else {
		std::cout << total_frames << " steps of simulation is done" << std::endl;
		exit(-1);
	}
}

int main(int argc, char** argv)
{
	if (!process_cmd_params(argc, argv)) {
		std::cout << "ERROR: Input arguments do not meet the requirements:" << std::endl;
		std::cout << "./program resolutionX(int) resolutionY(int) window_width(int) window_height(int) total_frames(int) preview_swtich(true/false)" << std::endl;
		return -1;
	}
	init_solver_params();

	solver = new VFXEpoch::Solvers::EulerGAS2D();
	if (!solver) {
		return -1;
	}
	gas_solver = dynamic_cast<VFXEpoch::Solvers::EulerGAS2D*>(solver);
	if(gas_solver){
		if(gas_solver = dynamic_cast<VFXEpoch::Solvers::EulerGAS2D*>(solver)){
			cout << "Successfully initialize gas solver" << endl;
		}
		else{
			cout << "Failed dynamic cast fluid base pointer to an instance" << endl;
			return -1;
		}
	}

	if (!preview) {
		
		std::cout << '\n' << "VFXEpoch Lib - Ubuntu v16.04 Unit Test" << '\n' << '\n';
		gas_solver->set_user_params(params);
		params.clear();
		params = gas_solver->get_user_params();

		cout << "Simulation User Parameters:" << endl;
		cout << params << endl;

		/************************************ Test solver functions ***********************************/
		bool isInit = gas_solver->init(params);	if(!isInit) return -1;
		gas_solver->set_source_location(params.dimension[0] / 2, params.dimension[1] / 2);
		gas_solver->set_external_force_location(VFXEpoch::VECTOR_COMPONENTS::Y, params.dimension[0]/2, params.dimension[1]/2);
		gas_solver->set_static_boundary(boundary_phi);

		for(int i = 0; i != total_frames; i++){
			cout << "****************** Frame " << i << " ******************" << endl;
			gas_solver->step();
			cout << "**************** Step " << i << " done ****************" << endl;
			cout << endl;
		}
		/************************************ Test solver functions ***********************************/
	}

	else {
		gas_solver->set_user_params(params);
		params.clear();
		params = gas_solver->get_user_params();

		cout << "Simulation User Parameters:" << endl;
		cout << params << endl;

		cout << "Simulation User Parameters:" << endl;
		cout << params << endl;

		bool isInit = gas_solver->init(params);	if(!isInit) return -1;
		gas_solver->set_source_location(params.dimension[0] / 2, params.dimension[1] / 2);
		gas_solver->set_external_force_location(VFXEpoch::VECTOR_COMPONENTS::Y, params.dimension[0]/2, params.dimension[1]/2);
		gas_solver->set_static_boundary(boundary_phi);

		Gluvi::init("VFXEPOCH - Example - Smoke", &argc, argv, win_width, win_height);
		Gluvi::camera = &cam;
		Gluvi::userDisplayFunc = display;
		Gluvi::userMouseFunc = mouse;
		glClearColor(0, 0, 0, 1);
		glutTimerFunc(1000 / 24, timer, total_frames);
		Gluvi::run();
	}

	gas_solver->close();
	if(solver) 
		delete solver;
	solver = NULL;

	return 0;
}