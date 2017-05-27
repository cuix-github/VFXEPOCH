/*******************************************************************************
    VFXEPOCH - Physically based simulation VFX

    Copyright (c) 2016 Snow Tsui <trevor.miscellaneous@gmail.com>

    All rights reserved. Use of this source code is governed by
    the MIT license as written in the LICENSE file.
*******************************************************************************/
#ifndef _UTL_ANALYSIS_H_
#define _UTL_ANALYSIS_H_

#include "UTL_Grid.h"
#include <vector>

using namespace std;

namespace VFXEpoch
{
	namespace Analysis
	{
		void
		computeCurl_uniform(VFXEpoch::Grid2DfScalarField& dest, VFXEpoch::Grid2DVector2DfField ref);

		void
		computeCurl_uniform_Stokes(VFXEpoch::Grid2DfScalarField& dest, VFXEpoch::Grid2DVector2DfField ref);

		void
		computeCurl_uniform_LS(VFXEpoch::Grid2DfScalarField& dest, VFXEpoch::Grid2DVector2DfField ref);

		void
		computeCurl_uniform_Richardson(VFXEpoch::Grid2DfScalarField& dest, VFXEpoch::Grid2DVector2DfField ref);

		void
		computeCurl_uniform(VFXEpoch::Grid3DVector3DfField& dest, VFXEpoch::Grid3DVector3DfField ref);

		void
		computeCurl_mac(VFXEpoch::Grid2DfScalarField& dest,
						VFXEpoch::Grid2DfScalarField u, VFXEpoch::Grid2DfScalarField v);

		void
		computeCurl_mac(VFXEpoch::Grid3DVector3DfField& dest,
						VFXEpoch::Grid2DfScalarField u, VFXEpoch::Grid2DfScalarField v, VFXEpoch::Grid2DfScalarField w);

		void
		computeGradient_uniform(VFXEpoch::Grid2DVector2DfField& dest, VFXEpoch::Grid2DfScalarField ref);

		void
		computeGradient_mac(VFXEpoch::Grid2DVector2DfField& dest, VFXEpoch::Grid2DfScalarField ref);

		void
		find_vector_from_vector_potential_2D(VFXEpoch::Grid2DVector2DfField& u, VFXEpoch::Grid2DfScalarField psi);

		void
		computeLaplace(VFXEpoch::Grid2DfScalarField& dest, VFXEpoch::Grid2DfScalarField& ref);

		// Comments for computeDivergence_uniform(...):
		// Q: Why using "float scale" parameters ?
		// A: Usually we want negative divergence field.
		// Desc:
		// Possion Equation usually has negative divergence calculation
		// at right hand side (RHS)
		void
		computeDivergence_uniform(VFXEpoch::Grid2DfScalarField& dest, VFXEpoch::Grid2DVector2DfField ref);

		void
		computeDivergence_with_weights_mac(VFXEpoch::Grid2DdScalarField& dest, VFXEpoch::Grid2DfScalarField u, VFXEpoch::Grid2DfScalarField v,
										   VFXEpoch::Grid2DfScalarField _uw, VFXEpoch::Grid2DfScalarField _vw);
		void
		computeDivergence_mac(VFXEpoch::Grid2DfScalarField& dest, VFXEpoch::Grid2DfScalarField u, VFXEpoch::Grid2DfScalarField v);
	}
}

#endif
