/*******************************************************************************
    VFXEPOCH - Physically based simulation VFX

    Copyright (c) 2016 Snow Tsui <trevor.miscellaneous@gmail.com>

    All rights reserved. Use of this source code is governed by
    the MIT license as written in the LICENSE file.
*******************************************************************************/
#ifndef _HELPERS_H_
#define _HELPERS_H_
#include "utl/UTL_General.h"
#include "utl/PCGSolver/pcg_solver.h"
#include "fluids/euler/SIM_EulerGAS.h"
#include "Eigen/Dense"
#include <iostream>
#include <iomanip>
#include <time.h>

#ifdef _WIN32
#include <Windows.h>
#endif

#ifdef __cplusplus__
#include <cstdlib>
#else
#include <stdlib.h>
#endif

using namespace std;
using namespace VFXEpoch;
using namespace VFXEpoch::Solvers;

namespace Helpers
{
	void displayFieldf(int row, int col, float* field);
	void displayScalarField(VFXEpoch::Grid2DfScalarField field);
	void displayVectorField(VFXEpoch::Grid2DVector2DfField field);
	void displayCellStatus(VFXEpoch::Grid2DCellTypes field);
	void randomInitScalarField(VFXEpoch::Grid2DfScalarField& field, float min, float max);
	void randomInitVectorField(VFXEpoch::Grid2DVector2DfField& field, float min, float max);
	void randomInitCellStatus(VFXEpoch::Grid2DCellTypes &field);
}


#endif
