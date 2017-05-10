/*******************************************************************************
    VFXEPOCH - Physically based simulation VFX

    Copyright (c) 2016 Snow Tsui <trevor.miscellaneous@gmail.com>

    All rights reserved. Use of this source code is governed by
    the MIT license as written in the LICENSE file.
*******************************************************************************/

// This utility is usually used for testing grids within several digits.
// Higher dimensions could violate the display format.

#include "Helpers.h"

const int Nx = 5;
const int Ny = 4;
const float source = 1.0f;
const float spacex = 1.0f / (Nx + 1);
const float spcaey = 1.0f / (Ny + 1);

using namespace Helpers;

int main(int argc, char** argv)
{
	std::cout << '\n' << "VFXEpoch Lib - Ubuntu v14.04 Unit Test" << '\n' << '\n';

	VFXEpoch::Grid2DfScalarField gridf(Nx, Ny);
	VFXEpoch::Grid2DVector2DfField gridv(Nx, Ny);
	VFXEpoch::Grid2DCellTypes domain_mask(Nx, Ny);
	VFXEpoch::Solvers::Euler_Fluid2D_Base *solver = new VFXEpoch::Solvers::EulerGAS2D();
	VFXEpoch::Solvers::EulerGAS2D* gas_solver = dynamic_cast<VFXEpoch::Solvers::EulerGAS2D*>(solver);

	// TODO: Dynamic cast check
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

	cout << "Domain status (F->NOTHING, B->SOMTHING):" << endl;
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
	params.dt = 0.01;
	params.buoyancy_alpha = 0.1;
	params.buoyancy_beta = 0.3;
	params.vort_conf_eps = 0.55;
	params.max_iterations = 30;

	// Error cast.
	gas_solver->set_user_params(params);
	params.clear();
	params = gas_solver->get_user_params();

	cout << "Simulation User Parameters:" << endl;
	cout << params;

	gas_solver->close();
	if(solver) delete solver;
	solver = NULL;

	return 0;
}
