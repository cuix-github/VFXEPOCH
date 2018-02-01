/*******************************************************************************
    VFXEPOCH - Physically based simulation VFX

    Copyright (c) 2016 Snow Tsui <trevor.miscellaneous@gmail.com>

    All rights reserved. Use of this source code is governed by
    the MIT license as written in the LICENSE file.
*******************************************************************************/

// This utility is usually used for testing grids within several digits.
// WARNING! Higher dimensions could violate the display format.

#include "Helpers.h"

const int Nx = 1024;
const int Ny = 1024;
const double source = 1.0f;
const double h = 1.f / Nx;
VFXEpoch::Vector2Dd c0(0.5,0.5), c1(0.7,0.5), c2(0.3,0.35), c3(0.5,0.7);
double rad0 = 0.4,  rad1 = 0.1,  rad2 = 0.1,   rad3 = 0.1;

using namespace Helpers;

/****************************** For Dbuge ******************************/
double circle_phi(const VFXEpoch::Vector2Dd& pos, const VFXEpoch::Vector2Dd& center, double radius){
	return (VFXEpoch::Dist2D(pos, center) - radius);	
}

double boundary_phi(const VFXEpoch::Vector2Dd& position) {
   double phi0 = -circle_phi(position, c0, rad0);
   // float phi1 = circle_phi(position, c1, rad1);
   // float phi2 = circle_phi(position, c2, rad2);
   // float phi3 = circle_phi(position, c3, rad3);
   // return std::min(std::min(phi0,phi1), std::min(phi2,phi3));
   return phi0;
}
/****************************** For Dbuge ******************************/

int main(int argc, char** argv)
{
	std::cout << '\n' << "VFXEpoch Lib - Ubuntu v16.04 Unit Test" << '\n' << '\n';
	VFXEpoch::Solvers::Euler_Fluid2D_Base *solver = new VFXEpoch::Solvers::EulerGAS2D();
	VFXEpoch::Solvers::EulerGAS2D* gas_solver = dynamic_cast<VFXEpoch::Solvers::EulerGAS2D*>(solver);

	if(gas_solver){
		if(gas_solver = dynamic_cast<VFXEpoch::Solvers::EulerGAS2D*>(solver)){
			cout << "Successfully initialize gas solver" << endl;
		}
		else{
			cout << "Failed dynamic cast fluid base pointer to a instance" << endl;
		}
		cout << endl;
	}
	else cout << "Euler GAS solver has been set up." << '\n' << '\n';

	EulerGAS2D::Parameters params;
	params.dimension = VFXEpoch::Vector2Di(Nx, Ny);
	params.h = h;
	params.dt = 0.005f;
	params.buoyancy_alpha = 0.1;
	params.buoyancy_beta = 0.3;
	params.vort_conf_eps = 0.55;
	params.max_iterations = 30;
	params.density_source = 20;
	params.external_force_strength = 10;
	params.max_iterations = 300;
	params.min_tolerance = 1e-5;
	params.diff = 0.01;
	params.visc = 0.01;

	gas_solver->set_user_params(params);
	params.clear();
	params = gas_solver->get_user_params();

	cout << "Simulation User Parameters:" << endl;
	cout << params << endl;

	/************************************ Test solver functions ************************************/
	bool isInit = gas_solver->init(params);	if(!isInit) return -1;
	gas_solver->set_source_location(Nx/2, Ny/2);
	gas_solver->set_external_force_location(VFXEpoch::VECTOR_COMPONENTS::Y, Nx/2, Ny/2);
	gas_solver->set_static_boundary(boundary_phi);
	
	// TODO: A small size simulation for unit test
	int total_frames = 300;
	for(int i = 0; i != total_frames; i++){
		cout << "****************** Frame " << i << " ******************" << endl;
		gas_solver->step();
		cout << "**************** Step " << i << " done ****************" << endl;
		cout << endl;
	}

	gas_solver->close();
	if(solver) delete solver;
	solver = NULL;
	/************************************ Test solver functions ************************************/

	return 0;
}