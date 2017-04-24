/*******************************************************************************
    VFXEPOCH - Physically based simulation VFX

    Copyright (c) 2016 Snow Tsui <trevor.miscellaneous@gmail.com>

    All rights reserved. Use of this source code is governed by
    the MIT license as written in the LICENSE file.
*******************************************************************************/
#ifndef _HELPERS_H_
#define _HELPERS_H_
#include "fluids/euler/SIM_Gas.h"
#include "fluids/euler/SIM_Mac.h"
#include "utl/UTL_General.h"
#include "utl/PCGSolver/pcg_solver.h"
#include <iostream>
#include <iomanip>
#include <time.h>
#include <GL/glut.h>

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

namespace Helper
{
	typedef struct _sim_params
	{
	public:
		_sim_params()
		{
			nx = ny = 32;
			dt = 0.01f;
			src = 1.0f;
			src_rate = 1.0f;
			diff = 0.01f;
			visc = 0.0f;
			user_force = 5.0f;
			streamer_len = 1.0f;
			num_particles = 10.0f;
			heat_source = 10.0f;
			particle_life_span_rev = 1.0f;
			vort_conf_eps = 0.0f;
			linear_solver_iterations = 30;
		}

		friend ostream& operator<<(ostream& os, const _sim_params& simParams)
		{
			os << std::setprecision(4) << setiosflags(ios::fixed);
			os << "Dimension			= " << simParams.nx << " x " << simParams.ny << endl;
			os << "Time step			= " << simParams.dt << endl;
			os << "Source				= " << simParams.src << endl;
			os << "Source rate			= " << simParams.src_rate << endl;
			os << "Diffuse				= " << simParams.diff << endl;
			os << "Heat source			= " << simParams.heat_source << endl;
			os << "Force				= " << simParams.user_force << endl;
			os << "Vorticity Confinement Epsilon	= " << simParams.vort_conf_eps << endl;
			os << "Streamer length			= " << simParams.streamer_len << endl;
			os << "Iterations in solver		= " << simParams.linear_solver_iterations << endl;
			return os;
		}

		int nx, ny;
		float dt;
		float src;
		float diff;
		float visc;
		float user_force;
		float streamer_len;
		float src_rate;
		float particle_life_span_rev;
		float heat_source;
		float vort_conf_eps;
		int linear_solver_iterations;
		int num_particles;
	}SimulationParameters;

	void displayFieldf(int row, int col, float* field);
	void displayScalarField(int row, int col, VFXEpoch::Grid2DfScalarField field);
	void displayVectorField(int row, int col, VFXEpoch::Grid2DVector2DfField field);
	void randomInitScalarField(VFXEpoch::Grid2DfScalarField& field, float min, float max);
	void randomInitVectorField(VFXEpoch::Grid2DVector2DfField& field, float min, float max);
}


#endif
