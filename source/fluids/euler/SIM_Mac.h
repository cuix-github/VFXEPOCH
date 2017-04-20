/*******************************************************************************
    VFXEPOCH - Physically based simulation VFX

    Copyright (c) 2016 Snow Tsui <trevor.miscellaneous@gmail.com>

    All rights reserved. Use of this source code is governed by
    the MIT license as written in the LICENSE file.
*******************************************************************************/
#ifndef _MAC_H_
#define _MAC_H_

#include "../../utl/UTL_Analysis.h"

namespace VFXEpoch
{
	class Mac2D
	{
	public:
		enum class FIELD2D
		{
			CURL,
			U_VEL,
			V_VEL,
			DENSITY,
			PRESSURE
		};
	private:
		VFXEpoch::Grid2DfScalarField curlField;
		VFXEpoch::Grid2DfScalarField curls;
		VFXEpoch::Grid2DfScalarField vel_u;
		VFXEpoch::Grid2DfScalarField vel_u1;
		VFXEpoch::Grid2DfScalarField vel_v;
		VFXEpoch::Grid2DdScalarField vel_v1;
		VFXEpoch::Grid2DfScalarField denField;
		VFXEpoch::Grid2DfScalarField densities;
		VFXEpoch::Grid2DfScalarField pressure;
		VFXEpoch::Grid2DfScalarField pressures;

		int dimenssion;
		float dt;
		float kappa;

	public:
		Mac2D();
		Mac2D(int N, float _dt, float _kappa);
		Mac2D(const Mac2D& source);
		Mac2D& operator=(const Mac2D& source);
		~Mac2D();

	public:
		bool InitGasSimParameters(int N, float _dt, float _kappa);
		void SetField(FIELD2D field, VFXEpoch::Grid2DfScalarField fieldData);
		void GetFeatures();
		VFXEpoch::Grid2DfScalarField GetField(FIELD2D field);
		void Reset();

	private:
		bool _init_gas_sim_parameters(int N, float _dt, float _kappa);
		void _set_field(FIELD2D field, VFXEpoch::Grid2DfScalarField fieldData);
		VFXEpoch::Grid2DfScalarField _get_field(FIELD2D field);
		void _features_calculation();
		void _clear();

		// Microsolvers
	private:
		void _get_curls();
	};
}

#endif
