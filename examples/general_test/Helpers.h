/*******************************************************************************
    VFXEPOCH - Physically based simulation VFX

    Copyright (c) 2016 Snow Tsui <trevor.miscellaneous@gmail.com>

    All rights reserved. Use of this source code is governed by
    the MIT license as written in the LICENSE file.
*******************************************************************************/
#ifndef _HELPERS_H_
#define _HELPERS_H_
#include <iostream>
#include <iomanip>
#include <time.h>

#ifdef __cplusplus__
#include <cstdlib>
#else
#include <stdlib.h>
#endif

using namespace std;
using namespace VFXEpoch;

#ifdef _WIN32
//TODO: Visual Studio pre-compile symbols
#include <Windows.h>
#endif

namespace Helper
{
	typedef struct _sim_params
	{
	public:
		_sim_params()
		{
			nx = ny = 32;
			dt = 0.01f;
			src = 100.0f;
			diff = 0.0f;
			user_force = 5.0f;
			padding1 = padding2 = 0.0f;
		}

		friend ostream& operator<<(ostream& os, const _sim_params& simParams)
		{
			os << std::setprecision(4) << setiosflags(ios::fixed);
			os << "Dimension = (" << simParams.nx << ", " << simParams.ny << ")" << endl;
			os << "Time step = " << simParams.dt << endl;
			os << "Source    = " << simParams.src << endl;
			os << "Diffuse   = " << simParams.diff << endl;
			os << "Force     = " << simParams.user_force << endl;
			return os;
		}

		int nx, ny;
		float dt;
		float src;
		float diff;
		float user_force;
		float padding1, padding2;
	}SimulationParameters;

	void displayField(int row, int col, float* field);
	void displayGrid(int row, int col, VFXEpoch::Grid2DfScalarField field);
	void displayVectorField(int row, int col, VFXEpoch::Grid2DVector2DfField field);
	void randomInitScalarField(VFXEpoch::Grid2DfScalarField& field, float min, float max);
	void randomInitVectorField(VFXEpoch::Grid2DVector2DfField& field, float min, float max);
}


#endif
