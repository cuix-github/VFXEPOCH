/*******************************************************************************
    VFXEPOCH - Physically based simulation VFX

    Copyright (c) 2016 Snow Tsui <trevor.miscellaneous@gmail.com>

    All rights reserved. Use of this source code is governed by
    the MIT license as written in the LICENSE file.
*******************************************************************************/

// This utility is usually used for testing grids within several digits.
// Higher dimensions could violate the display format.

#include "Helpers.h"

const int Nx = 9;
const int Ny = 9;
const float source = 1.0f;
const float h = 1.f / Nx;
VFXEpoch::Vector2Df c0(0.5,0.5), c1(0.7,0.5), c2(0.3,0.35), c3(0.5,0.7);
float rad0 = 0.4,  rad1 = 0.1,  rad2 = 0.1,   rad3 = 0.1;

using namespace Helpers;

/****************************** For Dbuge ******************************/
float circle_phi(const VFXEpoch::Vector2Df& pos, const VFXEpoch::Vector2Df& center, float radius){
	return (VFXEpoch::Dist2D(pos, center) - radius);	
}

float boundary_phi(const VFXEpoch::Vector2Df& position) {
   float phi0 = -circle_phi(position, c0, rad0);
   //float phi1 = circle_phi(position, c1, rad1);
   //float phi2 = circle_phi(position, c2, rad2);
   //float phi3 = circle_phi(position, c3, rad3);
   //return std::min(std::min(phi0,phi1), std::min(phi2,phi3));
   return phi0;
}
/****************************** For Dbuge ******************************/

int main(int argc, char** argv)
{
	std::cout << '\n' << "VFXEpoch Lib - Ubuntu v16.04 Unit Test" << '\n' << '\n';

	VFXEpoch::Grid2DfScalarField gridf(Nx, Ny, h, h);
	VFXEpoch::Grid2DVector2DfField gridv(Nx, Ny, h, h);
	VFXEpoch::Grid2DCellTypes domain_mask(Nx, Ny);
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
	Helpers::randomInitScalarField(gridf, -1.0, 1.0);
	Helpers::randomInitVectorField(gridv, -1.0, 1.0);
	Helpers::randomInitCellStatus(domain_mask);

	cout << "2D Scalar field access test:" << endl;
	displayScalarField(gridf);
	cout << endl;

	cout << "2D Vector field access test:" << endl;
	displayVectorField(gridv);
	cout << endl;

	cout << "Domain status (F->Fluid, B->Boundary):" << endl;
	displayCellStatus(domain_mask);
	cout << endl;

	int i = VFXEpoch::RandomI(0, Ny - 1);
	int j = VFXEpoch::RandomI(0, Nx - 1);
	std::cout << "The value of gridf at position" << "(" << i << ", " << j << ") is " << gridf(i, j) << '\n';
	std::cout << "The value of gridv at position" << "(" << i << ", " << j << ") is " << "Vector2Df(" << gridv(i, j).m_x << ", " << gridv(i, j).m_y << ")" << '\n';
	if(VFXEpoch::BOUNDARY_MASK::NOTHING == domain_mask(i, j))
	std::cout << "The value of domain_mask at position" << "(" << i << ", " << j << ") is Fluid" << '\n';
	else
	std::cout << "The value of domain_mask at position" << "(" << i << ", " << j << ") is Boundary" << '\n';
	cout << endl;

	EulerGAS2D::Parameters params;
	params.dimension = VFXEpoch::Vector2Di(Nx, Ny);
	params.h = h;
	params.dt = 0.005;
	params.buoyancy_alpha = 0.1;
	params.buoyancy_beta = 0.3;
	params.vort_conf_eps = 0.55;
	params.max_iterations = 30;
	params.density_source = 20;
	params.external_force_strength = 10;
	params.diff = 0.01;
	params.visc = 0.01;

	// Error cast.
	gas_solver->set_user_params(params);
	params.clear();
	params = gas_solver->get_user_params();

	cout << "Simulation User Parameters:" << endl;
	cout << params;

	/************************************ Test Vector functions ************************************/
	VFXEpoch::Vector2Df p0f(0.5f, 0.0f), p1f(0.0f, 10.0f);
	VFXEpoch::Vector2Df p2f = p0f / 0.1f;
	cout << endl << "The distance between p0 & p1 is: " << VFXEpoch::Dist2D(p0f, p1f) << endl;

	VFXEpoch::Vector2Dd p0d(0.9, 0.0), p1d(0.0, 5.0);
	cout << endl << "The distance between p0 & p1 is: " << VFXEpoch::Dist2D(p0d, p1d) << endl;
	/************************************ Test Vector functions ************************************/

	/************************************ Test solver functions ************************************/
	bool isInit = gas_solver->init(params);
	if(!isInit) return -1;
	gas_solver->set_source_location(2, 2);
	gas_solver->set_static_boundary(boundary_phi);
	gas_solver->set_external_force_location(VFXEpoch::VECTOR_COMPONENTS::Y, 4, 4);
	
	// TODO: A small size simulation for unit test
	int total_frames = 300;
	for(int i = 0; i != total_frames; i++){
		cout << endl << "Solver starts working for frame " << i;
		gas_solver->step();
	}

	gas_solver->close();
	if(solver) delete solver;
	solver = NULL;
	/************************************ Test solver functions ************************************/

	return 0;
}
