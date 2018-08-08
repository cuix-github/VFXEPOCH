/*******************************************************************************
    VFXEPOCH - Physically based simulation VFX

    Copyright (c) 2016 Snow Tsui <trevor.miscellaneous@gmail.com>

    All rights reserved. Use of this source code is governed by
    the MIT license as written in the LICENSE file.
*******************************************************************************/
#ifndef _UTL_LINEAR_SOLVERS_H_
#define _UTL_LINEAR_SOLVERS_H_
#include <math.h>
#include <iostream>
#include "UTL_Grid.h"
#include "UTL_Analysis.h"
#include "UTL_General.h"

using namespace std;

namespace VFXEpoch
{
	namespace LinearSolver
	{
		void
		GSSolve(VFXEpoch::Grid2DfScalarField& x, VFXEpoch::Grid2DfScalarField x0, VFXEpoch::BndConditionPerEdge b[], float coefMatrixAElement, float c, int iterations);

		void
		GSSolve(VFXEpoch::Grid2DdScalarField& x, VFXEpoch::Grid2DdScalarField x0, VFXEpoch::BndConditionPerEdge b[], float coefMatrixAElement, float c, int iterations);		

		void
		RBGSSolve(float h, VFXEpoch::Grid2DfScalarField& x, VFXEpoch::Grid2DfScalarField x0, VFXEpoch::BndConditionPerEdge b[], float coefMatrixAElement, float c);

		void
		JacobiSolve(VFXEpoch::Grid2DfScalarField& x, VFXEpoch::Grid2DfScalarField x0, VFXEpoch::BndConditionPerEdge b[], float coefMatrixAElement, float c, int iterations);

		void
		MultigridSolve_V_Cycle(float h, VFXEpoch::Grid2DfScalarField& x, VFXEpoch::Grid2DfScalarField x0, VFXEpoch::BndConditionPerEdge b[], float coefMatrixAElement, float c, int nSmooth);

		// TODO:
		// Conjugate Gradient (CG);
		// Pre-Conditioned Conjugate Gradient (PCG);
	}
}

#endif
