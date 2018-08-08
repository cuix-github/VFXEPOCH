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

#include <X11/Xlib.h>
#include <X11/Xutil.h>

int X_error_handler(Display *d, XErrorEvent *e)
{
        char msg[80];
        XGetErrorText(d, e->error_code, msg, sizeof(msg));

        fprintf(stderr, "Error %d (%s): request %d.%d\n",
                        e->error_code, msg, e->request_code, 
                        e->minor_code);
}

/***************************** For Visualization ******************************/
unsigned int win_width = 720;
unsigned int win_height = 720;
float sim_width = 100.0f;
int total_frames = 240;

// cam(bottom, height, left, right)
float pan_zoom_cam_bottom = -0.1f;
float pan_zoom_cam_height = 1.2f;
float pan_zoom_cam_left = -0.1f;
float pan_zoom_cam_right = pan_zoom_cam_left + (pan_zoom_cam_height * win_width) / win_height;
PanZoom2D cam(pan_zoom_cam_bottom, pan_zoom_cam_bottom + pan_zoom_cam_height, pan_zoom_cam_left, pan_zoom_cam_right);

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
VFXEpoch::Vector2Df c0(50, 50), c1(0.7,0.5), c2(0.3,0.35), c3(0.5,0.7);
VFXEpoch::Vector2Df o0(0.0, 0.0);
float rad0 = 40.0f,  rad1 = 0.1,  rad2 = 0.1,   rad3 = 0.1;
const int tracers = 100000;

using namespace Helpers;

/****************************** For Dbuge ******************************/
float circle_phi(const VFXEpoch::Vector2Df& pos, const VFXEpoch::Vector2Df& center, float radius) {
	return (VFXEpoch::Dist2D(pos, center) - radius);	
}

float boundary_phi(const VFXEpoch::Vector2Df& position) {
   float phi0 = -circle_phi(position, c0, rad0);
   // double phi1 = circle_phi(position, c1, rad1);
   // double phi2 = circle_phi(position, c2, rad2);
   // double phi3 = circle_phi(position, c3, rad3);
   // return std::min(std::min(phi0,phi1), std::min(phi2,phi3));
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
		total_frames = -1;
	}
	// Total Frame specificed
	else if (4 == argc) {
		std::string resx = argv[1];
		std::string resy = argv[2];
		std::string frames = argv[3];
		std::string::size_type sz;
		params.dimension[0] = std::stoi(resx, &sz);
		params.dimension[1] = std::stoi(resy, &sz);
		total_frames = std::stoi(frames, &sz);
	}
	// Window Size specified
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
		preview = true;
		total_frames = -1;
	}
	// Windows Size & Total Frame specified
	else if (6 == argc) {
		std::string resx = argv[1];
		std::string resy = argv[2];
		std::string frames = argv[3];
		std::string window_width = argv[4];
		std::string window_height = argv[5];
		std::string::size_type sz;
		params.dimension[0] = std::stoi(resx, &sz);
		params.dimension[1] = std::stoi(resy, &sz);
		win_width = std::stoi(window_width, &sz);
		win_height = std::stoi(window_height, &sz);
		preview = true;
		total_frames = std::stoi(frames, &sz);
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

bool init_particles(int num_particles)
{
	
}


bool init_solver_params()
{
	params.origin = o0;
	
	// params.dimension = VFXEpoch::Vector2Di(32, 32);
	// Since druing debuggin, we change the dimension very often, so we move this param
	// to be passed from command line
	
	params.h = sim_width / (float)params.dimension[0];
	params.dt = 0.1f;
	params.buoyancy_alpha = 0.1;
	params.buoyancy_beta = 0.3;
	params.vort_conf_eps = 0.55;
	params.density_source = 20;
	params.external_force_strength = 50;
	params.use_gravity = true;
	params.max_iterations = 300;
	params.min_tolerance = 1e-5;
	params.diff = 0.01;
	params.visc = 0.01;
	return true;
}

vector<VFXEpoch::Vector2Df> particles;
void display()
{
	glPushMatrix();

	float glScale = 1.0f / (params.h* params.dimension[0]);
	glScaled(glScale, glScale, glScale);

	OpenGL_Utility::draw_grid2d(o0, params.h, params.dimension[0], params.dimension[1], VFXEpoch::Vector3Df(0.5, 0.5, 0.5));

	OpenGL_Utility::draw_circle2d(c0, rad0, 100, VFXEpoch::Vector3Df(0.0, 1.0, 0.0));
	
	OpenGL_Utility::draw_particles2d(gas_solver->get_particles(), 10, VFXEpoch::Vector3Df(1.0, 0.0, 0.0), true);

	OpenGL_Utility::draw_arrows(gas_solver, 0.1f, VFXEpoch::Vector3Df(0.0, 1.0, 0.0));

	glPopMatrix();
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
	else if (total_frames == -1) {
		int i = (total_frames - arg) + 1;
		cout << "****************** Frame " << i << " ******************" << endl;
		gas_solver->step();
		cout << "**************** Step " << i << " done ****************" << endl;
		cout << endl;
		glutPostRedisplay();
		glutTimerFunc(1000 / 24, timer, 0);
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

	// Silent run for the simulation
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

	// Toggle the OpenGL window to preview the simulation
	else {
		gas_solver->set_user_params(params);
		params.clear();
		params = gas_solver->get_user_params();

		cout << "Simulation User Parameters:" << endl;
		cout << params << endl;

		bool isInit = gas_solver->init(params);	if(!isInit) return -1;
		gas_solver->set_source_location(params.dimension[0] / 2 - 10, params.dimension[1] / 2);
		gas_solver->set_source_location(params.dimension[0] / 2 + 1, params.dimension[1] / 2);
		gas_solver->set_source_location(params.dimension[0] / 2 + 2, params.dimension[1] / 2);
		gas_solver->set_source_location(params.dimension[0] / 2 + 3, params.dimension[1] / 2);
		gas_solver->set_external_force_location(VFXEpoch::VECTOR_COMPONENTS::Y, params.dimension[0]/2, params.dimension[1]/2);
		gas_solver->set_static_boundary(boundary_phi);

		VFXEpoch::Particle2Df p;
		p.pos = VFXEpoch::Vector2Df(42.079145937118049,  47.019514491033981);
		gas_solver->add_particles(p);
		p.pos = VFXEpoch::Vector2Df(33.079145937118049,  52.019514491033981);
		gas_solver->add_particles(p);
		p.pos = VFXEpoch::Vector2Df(55.079145937118049,  25.019514491033981);
		gas_solver->add_particles(p);
		p.pos = VFXEpoch::Vector2Df(30.079145937118049,  45.019514491033981);
		gas_solver->add_particles(p);
		p.pos = VFXEpoch::Vector2Df(23.079145937118049,  66.019514491033981);
		gas_solver->add_particles(p);
		p.pos = VFXEpoch::Vector2Df(11.079145937118049,  21.019514491033981);
		gas_solver->add_particles(p);
		p.pos = VFXEpoch::Vector2Df(7.079145937118049,  10.019514491033981);
		gas_solver->add_particles(p);
		p.pos = VFXEpoch::Vector2Df(9.079145937118049,  27.019514491033981);
		gas_solver->add_particles(p);

		Gluvi::init("VFXEPOCH - Smoke Prototype", &argc, argv, win_width, win_height);
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

	XSetErrorHandler(X_error_handler);

	return 0;
}