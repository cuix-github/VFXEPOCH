/*******************************************************************************
    VFXEPOCH - Physically based simulation VFX

    Copyright (c) 2016 Snow Tsui <trevor.miscellaneous@gmail.com>

    All rights reserved. Use of this source code is governed by
    the MIT license as written in the LICENSE file.
*******************************************************************************/
#include "utl/Common/Helpers.h"

const int Nx = 4;
const int Ny = 10;
const float source = 1.0f;
const float spacex = 1.0f / (Nx + 1);
const float spcaey = 1.0f / (Ny + 1);

using namespace Helpers;

int main(int argc, char** argv)
{
	std::cout << '\n' << "VFXEpoch Lib - Ubuntu v14.04 Unit Test" << '\n' << '\n';

	VFXEpoch::Grid2DfScalarField gridf(Nx, Ny);
	VFXEpoch::Grid2DVector2DfField gridv(Nx, Ny);
	Helpers::randomInitScalarField(gridf, -1.0, 1.0);
	Helpers::randomInitVectorField(gridv, -1.0, 1.0);

	cout << "2D Scalar field access test:" << endl;
	displayScalarField(gridf);
	cout << endl;

	cout << "2D Vector field access test:" << endl;
	displayVectorField(gridv);
	cout << endl;

	std::cout << '\n' << "The value of gridf at position" << "(2, 3) is " << gridf(2, 3) << '\n';
	std::cout << '\n' << "The value at gridv at position" << "(2, 1) is VFXEpoch::Vector2Df(" << gridv(2, 1).m_x << ", " << gridv(2, 1).m_y << ")" << '\n';
	cout << endl;

	return 0;
}
