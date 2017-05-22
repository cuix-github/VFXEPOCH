/*******************************************************************************
    VFXEPOCH - Physically based simulation VFX

    Copyright (c) 2016 Snow Tsui <trevor.miscellaneous@gmail.com>

    All rights reserved. Use of this source code is governed by
    the MIT license as written in the LICENSE file.
*******************************************************************************/
#include "UTL_Analysis.h"

namespace VFXEpoch
{
	namespace Analysis
	{
		// How to calculate curl in the uniform grid��
		// Localtion for storing the curl: Omega
		// 1.Calculate u component gradient from Cell6 & Cell4
		// 2.Calculate v component gradient from Cell8 & Cell2
		// 3.Subtract du from dv to get curl (dv - du).
		// Omega = dv - du;
		/*
			-------------------------
			|				| 			| 			|
			| Cell1	| Cell2 | Cell3 |
			|				|				|				|
			-------------------------
			|				|				|				|
			| Cell4	|  Curl | Cell6	|
			|				|				|				|
			-------------------------
			|				|				|				|
			| Cell7 | Cell8 | Cell9 |
			|				|				|				|
			-------------------------
		*/
		void
		computeCurl_uniform(VFXEpoch::Grid2DfScalarField& dest, VFXEpoch::Grid2DVector2DfField ref) {
			float hx = 1.0f / (dest.getDimX() - 2);
			float hy = 1.0f / (dest.getDimY() - 2);
			for (int i = 1; i != dest.getDimY() - 1; i++){
				for (int j = 1; j != dest.getDimX() - 1; j++){
					float du, dv, curl;
					du = 0.5f * (ref.getData(i + 1, j).m_x - ref.getData(i - 1, j).m_x) / hx;
					dv = 0.5f * (ref.getData(i, j + 1).m_y - ref.getData(i, j - 1).m_y) / hy;
					curl = dv - du;
					dest(i, j) = curl;
				}
			}
		}

		void
		computeCurl_uniform_Stokes(VFXEpoch::Grid2DfScalarField& dest, VFXEpoch::Grid2DVector2DfField ref){
			// TODO: Stkes theorem based vorticity calculation
			float hx = 1.0f / (dest.getDimX() - 2);
			float hy = 1.0f / (dest.getDimY() - 2);
			for (int i = 1; i != dest.getDimY() - 1; i++){
				for (int j = 1; j != dest.getDimX() - 1; j++){
					float du0, du1, dv0, dv1, coef;
					coef = 1.0f / (8.0f * (hx * hy));
					du0 = hx * (ref(i - 1, j - 1).m_x + 2 * ref(i, j - 1).m_x + ref(i + 1, j - 1).m_x);
					du1 = hx * (ref(i + 1, j + 1).m_x + 2 * ref(i, j + 1).m_x + ref(i - 1, j + 1).m_x);
					dv0 = hy * (ref(i + 1, j - 1).m_y + 2 * ref(i + 1, j).m_y + ref(i + 1, j + 1).m_y);
					dv1 = hy * (ref(i - 1, j + 1).m_y + 2 * ref(i - 1, j).m_y + ref(i - 1, j - 1).m_y);
					dest(i, j) = coef * (du0 + dv0 - du1 - dv1);
				}
			}
		}

		void
		computeCurl_uniform_LS(VFXEpoch::Grid2DfScalarField& dest, VFXEpoch::Grid2DVector2DfField ref) {
			// TODO: Least Square to get vorticity
			float hx = 1.0f / (dest.getDimX() - 2);
			float hy = 1.0f / (dest.getDimY() - 2);
			for (int i = 2; i != dest.getDimY() - 2; i++){
				for (int j = 2; j != dest.getDimX() - 2; j++){
					float dv, du, coefx, coefy;
					coefx = 1.0f / (10.0f * hx);
					coefy = 1.0f / (10.0f * hy);
					du = 2 * ref(i, j + 2).m_y + ref(i, j + 1).m_y - ref(i, j - 1).m_y - 2 * ref(i, j - 2).m_y;
					dv = 2 * ref(i + 2, j).m_x + ref(i + 1, j).m_x - ref(i - 1, j).m_x - 2 * ref(i - 2, j).m_x;
					dest(i, j) = coefx * dv - coefy * du;
				}
			}
		}

		void
		computeCurl_uniform_Richardson(VFXEpoch::Grid2DfScalarField& dest, VFXEpoch::Grid2DVector2DfField ref){
			// TODO: Richarson theorem based vorticity calculation
			float hx = 1.0f / (dest.getDimX() - 2);
			float hy = 1.0f / (dest.getDimY() - 2);
			for (int i = 2; i != dest.getDimY() - 2; i++){
				for (int j = 2; j != dest.getDimX() - 2; j++){
					float dv, du, coefx, coefy;
					coefx = 1.0f / (12.0f * hx);
					coefy = 1.0f / (12.0f * hy);
					du = -ref(i, j + 2).m_y + 8 * ref(i, j + 1).m_y - 8 * ref(i, j - 1).m_y + ref(i, j - 2).m_y;
					dv = -ref(i + 2, j).m_x + 8 * ref(i + 1, j).m_x - 8 * ref(i - 1, j).m_x + ref(i - 2, j).m_x;
					dest(i, j) = coefx * dv - coefy * du;
				}
			}
		}

		void
		computeCurl_uniform(VFXEpoch::Grid3DVector3DfField& dest, VFXEpoch::Grid3DVector3DfField ref)
		{
			// TODO: Process 3D curl calculation in the uniform grid

		}

		// How to calculate curl in the mac grid��
		// Omega = Nabla x velocity field;
		// du = u2 - u1;
		// dv = v2 - v1;
		// Omega = du - dv;
		//
		/*
		-------------------------
		|			|			|
		|			u1
		|			|			|
		----v1-----Curl----v2----
		|			|			|
		|			u2			|
		|			|			|
		-------------------------
		*/
		void
		computeCurl_mac(VFXEpoch::Grid2DfScalarField& dest, VFXEpoch::Grid2DfScalarField u, VFXEpoch::Grid2DfScalarField v) {
			for (int i = 1; i != dest.getDimY() - 1; i++){
				for (int j = 1; j != dest.getDimX() - 1; j++){
					float du, dv;
					float curl;
					du = (u(i, j) - u.getData(i - 1, j)) / u.getDy();
					dv = (v(i, j) - v.getData(i, j - 1)) / v.getDx();
					curl = dv - du;
					dest(i, j) = curl;
				}
			}
		}

		void
		computeCurl_mac(VFXEpoch::Grid3DVector3DfField& dest, VFXEpoch::Grid3DfScalarField u, VFXEpoch::Grid3DfScalarField v, VFXEpoch::Grid3DfScalarField w)
		{
			// TODO: Process 3D curl calculation in the uniform grid
		}

		void
		computeGradient_uniform(VFXEpoch::Grid2DVector2DfField& dest, VFXEpoch::Grid2DfScalarField ref) {
			// TODO: Compute gradients of a scalar field.
			// Bugs here, central difference coef 0.5f * N;
			for (int i = 1; i != dest.getDimY() - 1; i++){
				for (int j = 1; j != dest.getDimX() - 1; j++){
					VFXEpoch::Vector2Df data(0.5f * (dest.getDimY() - 2) * (ref(i + 1, j) - ref(i - 1, j)),
											 0.5f * (dest.getDimX() - 2) * (ref(i, j + 1) - ref(i, j - 1)));
					dest(i, j) = data;
				}
			}
		}

		void
		computeGradient_mac(VFXEpoch::Grid2DVector2DfField& dest, VFXEpoch::Grid2DfScalarField ref)
		{
			// TODO: Compute gradients on staggered grid.
		}

		void
		find_vector_from_vector_potential_2D(VFXEpoch::Grid2DVector2DfField& u, VFXEpoch::Grid2DfScalarField psi)
		{
			// TODO: Calculate vector potential from the given function (Psi for stream function)
			float dpsi_dy, dpsi_dx;
			int Nx = u.getDimX() - 2;
			int Ny = u.getDimY() - 2;
			u.zeroVectors();

			for (int i = 1; i != u.getDimY() - 1; i++){
				for (int j = 1; j != u.getDimX() - 1; j++){
					dpsi_dy = 0.5f * (psi(i + 1, j) - psi(i - 1, j)) * Ny;
					dpsi_dx = 0.5f * (psi(i, j + 1) - psi(i, j - 1)) * Nx;
					u(i, j).m_x = dpsi_dy;
					u(i, j).m_y = -dpsi_dx;
				}
			}
		}

		void
		computeLaplace(VFXEpoch::Grid2DfScalarField& dest, VFXEpoch::Grid2DfScalarField& ref)
		{
			for (int i = 1; i != dest.getDimY() - 1; i++){
				for (int j = 1; j != dest.getDimX() - 1; j++){
					dest(i, j) = (ref(i + 1, j) + ref(i - 1, j) + ref(i, j + 1) + ref(i, j - 1) - 4.0f * ref(i, j));
				}
			}
		}

		void
		computeDivergence_uniform(VFXEpoch::Grid2DfScalarField& dest, VFXEpoch::Grid2DVector2DfField ref) {
			int Nx = ref.getDimX() - 2;
			int Ny = ref.getDimY() - 2;
			for (int i = 1; i != dest.getDimY() - 1; i++){
				for (int j = 1; j != dest.getDimX() - 1; j++){
					float div(0.0);
					float du(0.0f), dv(0.0f);
					du = 0.5f * (ref.getData(i + 1, j).m_x - ref.getData(i - 1, j).m_x) / Nx;
					dv = 0.5f * (ref.getData(i, j + 1).m_y - ref.getData(i, j - 1).m_y) / Ny;
					div = du + dv;
					dest(i, j) = div;
				}
			}
		}

		void
		computeDivergence_with_weights_mac(VFXEpoch::Grid2DfScalarField& dest, VFXEpoch::Grid2DfScalarField u, VFXEpoch::Grid2DfScalarField v, 
										   VFXEpoch::Grid2DfScalarField _uw, VFXEpoch::Grid2DfScalarField _vw) {
			assert(u.getDimY() == _uw.getDimY() && u.getDimX() == _uw.getDimX() &&
				   v.getDimY() == _vw.getDimY() && v.getDimX() == _vw.getDimX() &&
				   u.getDimY() == v.getDimX() && u.getDimX() == v.getDimY());
			dest.clear();
			float h = u.getDx();
			LOOP_GRID2D_WITHOUT_DOMAIN_BOUNDARY(dest){
				dest(i, j) -= _uw(i, j+1) * u(i, j+1) / h;
				dest(i, j) += _uw(i, j) * u(i, j) / h;
				dest(i, j) -= _vw(i+1, j) * v(i, j) / h;
				dest(i, j) += _vw(i, j) * v(i, j) / h;
			}
		}

		void
		computeDivergence_mac(VFXEpoch::Grid2DfScalarField& dest, VFXEpoch::Grid2DfScalarField u, VFXEpoch::Grid2DfScalarField v) {
			// TODO: Compute divergence on staggered grid
		}
	}
}
