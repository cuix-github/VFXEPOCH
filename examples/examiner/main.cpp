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
const int Ny = 8;
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
	if(!solver) return -1;
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
	std::cout << "The value at gridv at position" << "(" << i << ", " << j << ") is " << "Vector2Df(" << gridv(i, j).m_x << ", " << gridv(i, j).m_y << ")" << '\n';
	if(VFXEpoch::BOUNDARY_MASK::NOTHING == domain_mask(i, j))
	std::cout << "The value of domain_mask at position" << "(" << i << ", " << j << ") is Fluid" << '\n';
	else
	std::cout << "The value of domain_mask at position" << "(" << i << ", " << j << ") is Boundary" << '\n';
	cout << endl;

	EulerGAS2D::Parameters params;
	params._dimension = VFXEpoch::Vector2Di(Nx, Ny);
	params._dt = 0.01;
	params._buoyancy_alpha = 0.1;
	params._buoyancy_beta = 0.3;
	params._vort_conf_eps = 0.55;
	params._max_iterations = 30;

	// Error cast.
	solver->set_user_params(params);

	if(solver) delete solver;
	solver = NULL;

	return 0;
}
