/*******************************************************************************
    VFXEPOCH - Physically based simulation VFX

    Copyright (c) 2016 Snow Tsui <trevor.miscellaneous@gmail.com>

    All rights reserved. Use of this source code is governed by
    the MIT license as written in the LICENSE file.
*******************************************************************************/
#include "UTL_LinearSolvers.h"

namespace VFXEpoch
{
	namespace LinearSolver
	{
		void
		GSSolve(VFXEpoch::Grid2DfScalarField& x, VFXEpoch::Grid2DfScalarField x0, VFXEpoch::BndConditionPerEdge b[], float coefMatrixAElement, float c, int iterations)	{
			for (int m = 0; m != iterations; m++)	{
				for (int i = 1; i != x.getDimY() - 1; i++)	{
					for (int j = 1; j != x.getDimX() - 1; j++)	{
						float sumFluxout = x(i + 1, j) + x(i - 1, j) + x(i, j + 1) + x(i, j - 1);
						sumFluxout *= coefMatrixAElement;
						float curApproximation = x0(i, j) + sumFluxout;
						curApproximation /= c;
						x(i, j) = curApproximation;
					}
				}
				x.setBoundaries(b[0].boundaryType, b[0].side);
				x.setBoundaries(b[1].boundaryType, b[1].side);
				x.setBoundaries(b[2].boundaryType, b[2].side);
				x.setBoundaries(b[3].boundaryType, b[3].side);
				x.setBoundariesOnCorners();
			}
		}

		void
		GSSolve(VFXEpoch::Grid2DdScalarField& x, VFXEpoch::Grid2DdScalarField x0, VFXEpoch::BndConditionPerEdge b[], float coefMatrixAElement, float c, int iterations)	{
			for (int m = 0; m != iterations; m++)	{
				for (int i = 1; i != x.getDimY() - 1; i++)	{
					for (int j = 1; j != x.getDimX() - 1; j++)	{
						double sumFluxout = x(i + 1, j) + x(i - 1, j) + x(i, j + 1) + x(i, j - 1);
						sumFluxout *= coefMatrixAElement;
						double curApproximation = x0(i, j) + sumFluxout;
						curApproximation /= c;
						x(i, j) = curApproximation;
					}
				}
				x.setBoundaries(b[0].boundaryType, b[0].side);
				x.setBoundaries(b[1].boundaryType, b[1].side);
				x.setBoundaries(b[2].boundaryType, b[2].side);
				x.setBoundaries(b[3].boundaryType, b[3].side);
				x.setBoundariesOnCorners();
			}
		}

		void
		RBGSSolve(float h, VFXEpoch::Grid2DfScalarField& x, VFXEpoch::Grid2DfScalarField x0, VFXEpoch::BndConditionPerEdge b[], float coefMatrixAElement, float c) {
			for (int rb = 0; rb != 2; rb++)	{
				for (int i = 1; i != x.getDimY() - 1; i++) {
					for (int j = 1; j != x.getDimX() - 1; j++) {
						if ((i + j) % 2 == rb)	{
							float sumFluxout = x(i + 1, j) + x(i - 1, j) + x(i, j + 1) + x(i, j - 1);
							sumFluxout *= coefMatrixAElement;
							float curApproximation = x0(i, j) + sumFluxout;
							curApproximation /= c;
							x(i, j) = curApproximation;
						}
					}
				}
				// TODO: Boundary values need to be set
				x.setBoundaries(b[0].boundaryType, b[0].side);
				x.setBoundaries(b[1].boundaryType, b[1].side);
				x.setBoundaries(b[2].boundaryType, b[2].side);
				x.setBoundaries(b[3].boundaryType, b[3].side);
				x.setBoundariesOnCorners();
			}
		}

		void
		JacobiSolve(VFXEpoch::Grid2DfScalarField& x, VFXEpoch::Grid2DfScalarField x0, VFXEpoch::BndConditionPerEdge b[], float coefMatrixAElement, float c, int iterations)	{
				// Extra memory to avoid covering previous value
				VFXEpoch::Grid2DfScalarField auxiliary(x0.getDimX(), x0.getDimY(), (float)(1 / (x0.getDimX() - 1)), (float)(1 / (x0.getDimY() - 1)));

				int xCells = x0.getDimX();
				int yCells = x0.getDimY();

				float h = 1.0f / (float)(x0.getDimY() - 1);

				for (int m = 0; m != iterations; m++)	{
					for (int i = 1; i != yCells - 1; i++)	{
						for (int j = 1; j != xCells - 1; j++)	{
							float sumFluxout = x(i + 1, j) + x(i - 1, j) + x(i, j + 1) + x(i, j - 1);
							sumFluxout *= coefMatrixAElement;
							float curApproximation = x0(i, j) * h * h + sumFluxout;
							curApproximation /= c;
							auxiliary(i, j) = curApproximation;
						}
					}

					for (int i = 1; i != yCells - 1; i++)
						for (int j = 1; j != xCells - 1; j++)
							x(i, j) = auxiliary(i, j);

					x.setBoundaries(b[0].boundaryType, b[0].side);
					x.setBoundaries(b[1].boundaryType, b[1].side);
					x.setBoundaries(b[2].boundaryType, b[2].side);
					x.setBoundaries(b[3].boundaryType, b[3].side);
					x.setBoundariesOnCorners();
				}
				auxiliary.clear();
		}

		// TODO: Fix bugs in V-Cycle Multigrid algorithm
		void
		MultigridSolve_V_Cycle(float h, VFXEpoch::Grid2DfScalarField& x, VFXEpoch::Grid2DfScalarField x0, VFXEpoch::BndConditionPerEdge b[], float coefMatrixAElement, float c, int nSmooth) {
			int L = x.getDimY() - 2;
			if (L == 1){
				x(1, 1) = 0.25f * (x(0, 1) + x(1, 0) + x(1, 2) + x(2, 1) + x0(1, 1));
				return;
			}

			for (int i = 0; i != nSmooth; i++)
				VFXEpoch::LinearSolver::RBGSSolve(h, x, x0, b, 1, 4);

			VFXEpoch::Grid2DfScalarField r(L + 2, L + 2);
			VFXEpoch::Zeros(r);
			VFXEpoch::Grid2DfScalarField laplace_x(L + 2, L + 2, h, h);
			VFXEpoch::Analysis::computeLaplace(laplace_x, x);
			r = x0 + laplace_x;

			int L2 = L / 2;
			VFXEpoch::Grid2DfScalarField R(L2 + 2, L2 + 2);
			for (int i = 1; i <= L2; i++) {
				int m = 2 * i - 1;
				for (int j = 1; j <= L2; j++) {
					int n = 2 * j - 1;
					R(i, j) = 0.25f * (r(m, n) + r(m + 1, n) + r(m, n + 1) + r(m + 1, n + 1));
				}
			}

			VFXEpoch::Grid2DfScalarField V(L2 + 2, L2 + 2);
			VFXEpoch::Zeros(V);
			float H = 2 * h;
			VFXEpoch::LinearSolver::MultigridSolve_V_Cycle(H, V, R, b, 1, 4, nSmooth);

			VFXEpoch::Grid2DfScalarField v(L + 2, L + 2);
			for (int i = 1; i <= L2; i++){
				int m = 2 * i - 1;
				for (int j = 1; j <= L2; j++){
					int n = 2 * j - 1;
					v(m, n) = v(m + 1, n) = v(m, n + 1) = v(m + 1, n + 1) = V(i, j);
				}
			}

			x += v;
			for (int i = 0; i != nSmooth; i++)
				//VFXEpoch::LinearSolver::JacobiSolve(x, x0, b, 1, 4, 30);
				VFXEpoch::LinearSolver::RBGSSolve(h, x, x0, b, 1, 4);

			x.setBoundaries(b[0].boundaryType, b[0].side);
			x.setBoundaries(b[1].boundaryType, b[1].side);
			x.setBoundaries(b[2].boundaryType, b[2].side);
			x.setBoundaries(b[3].boundaryType, b[3].side);
			x.setBoundariesOnCorners();
		}
	}
}
