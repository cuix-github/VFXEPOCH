/*******************************************************************************
    VFXEPOCH - Physically based simulation VFX

    Copyright (c) 2016 Snow Tsui <trevor.miscellaneous@gmail.com>

    All rights reserved. Use of this source code is governed by
    the MIT license as written in the LICENSE file.
	
      v6  v2   v5
        \  |  /
         \ | /
          \|/
      v3---v0---v1
          /|\
         / | \
        /  |  \
      v7  v4   v8
*******************************************************************************/
#ifndef _SIM_LBM_H_
#define _SIM_LBM_H_
#include "../../utl/UTL_Grid.h"
#include "../../utl/UTL_Matrix.h"
#include "../../utl/UTL_Vector.h"
#include "../../utl/UTL_General.h"
#include "../../utl/UTL_LinearSolvers.h"

#include <math.h>

using namespace std;

#define LBMIDX(k, i, j) (((k * i)) + j)

namespace VFXEpoch
{
	namespace Solvers
	{
		/*
		* Lattice Boltzmann Method (LBM)
		*/
		struct LBM2DParameters
		{
			float tau, rho;
			const float auxFactor1 = 4.f / 9.f;
			const float auxFactor2 = 1.f / 9.f;
			const float auxFaCTOR3 = 1.F / 36.f;
			const VFXEpoch::Vector2Di e0 = VFXEpoch::Vector2Di(0, 0);
			const VFXEpoch::Vector2Di e1 = VFXEpoch::Vector2Di(1, 0);
			const VFXEpoch::Vector2Di e2 = VFXEpoch::Vector2Di(0, 1);
			const VFXEpoch::Vector2Di e3 = VFXEpoch::Vector2Di(-1, 0);
			const VFXEpoch::Vector2Di e4 = VFXEpoch::Vector2Di(0, -1);
			const VFXEpoch::Vector2Di e5 = VFXEpoch::Vector2Di(1, 1);
			const VFXEpoch::Vector2Di e6 = VFXEpoch::Vector2Di(-1, 1);
			const VFXEpoch::Vector2Di e7 = VFXEpoch::Vector2Di(-1, -1);
			const VFXEpoch::Vector2Di e8 = VFXEpoch::Vector2Di(1, -1);
		};

		class LBM2D
		{
		public:
			LBM2D();
			LBM2D(int resx, int resy);
			LBM2D(const LBM2D& source);
			LBM2D& operator=(const LBM2D& source);
			~LBM2D();

			// TODO: Will change to private and implement interfaces
		public:
			bool _initialize(int x, int y);
			void _set_sim_params(float rho, float tau);
			void _stream();
			void _collide();
			void _process_boundary();
			void _set_solid_at_cell(BOUNDARY_MASK flag, int x, int y);
			void _bounce_back();
			void _clear();

		private:
			int resolutionX, fieldX;
			int resolutionY, fieldY;
			VFXEpoch::Grid2DVector2DfField vel;
			VFXEpoch::Grid2DfScalarField mag_vel;
			VFXEpoch::Grid2D<BOUNDARY_MASK> solid_mask;
			VFXEpoch::Grid2DfScalarField v0, v1, v2, v3, v4, v5, v6, v7, v8;
			VFXEpoch::Grid2DfScalarField auxv0, auxv1, auxv2, auxv3, auxv4, auxv5, auxv6, auxv7, auxv8;
			LBM2DParameters params;
		};

		// TODO: Implementation of 3D Lattice Boltzmann Method
		// Follow by Nils Thuerey
		class LBM3D
		{
		public:
			LBM3D();
			LBM3D(const LBM3D& source);
			LBM3D& operator=(const LBM3D& source);
			~LBM3D();

		private:
			bool _initialize();
			void _clear();
		};
	}
}

#endif
