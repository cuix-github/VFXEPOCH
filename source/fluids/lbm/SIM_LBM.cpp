/*******************************************************************************
    VFXEPOCH - Physically based simulation VFX

    Copyright (c) 2016 Snow Tsui <trevor.miscellaneous@gmail.com>

    All rights reserved. Use of this source code is governed by
    the MIT license as written in the LICENSE file.
*******************************************************************************/
#include "SIM_LBM.h"

using namespace VFXEpoch;
using namespace VFXEpoch::Solvers;

LBM2D::LBM2D()
{
	_initialize(1, 1);
}

LBM2D::LBM2D(int resx, int resy)
{
	_initialize(resx, resy);
}

LBM2D::LBM2D(const LBM2D& source)
{
	if (resolutionX <= 0 || resolutionY <= 0)
		return;
	else
	{
		resolutionX = source.resolutionX; fieldX = resolutionX - 2;
		resolutionY = source.resolutionY; fieldY = resolutionY - 2;

		v0 = source.v0;
		v1 = source.v1;
		v2 = source.v2;
		v3 = source.v3;
		v4 = source.v4;
		v5 = source.v5;
		v6 = source.v6;
		v7 = source.v7;
		v8 = source.v8;

		auxv0 = source.auxv0;
		auxv1 = source.auxv1;
		auxv2 = source.auxv2;
		auxv3 = source.auxv3;
		auxv4 = source.auxv4;
		auxv5 = source.auxv5;
		auxv6 = source.auxv6;
		auxv7 = source.auxv7;
		auxv8 = source.auxv8;

		vel = source.vel;
		mag_vel = source.mag_vel;
		solid_mask = source.solid_mask;
	}
}

LBM2D&
LBM2D::operator=(const LBM2D& source)
{
	if (resolutionX <= 0 || resolutionY <= 0){
		LBM2D* nullObj = nullptr;
		return *nullObj;
	}
	else{
		resolutionX = source.resolutionX; fieldX = resolutionX - 2;
		resolutionY = source.resolutionY; fieldY = resolutionY - 2;

		v0 = source.v0;
		v1 = source.v1;
		v2 = source.v2;
		v3 = source.v3;
		v4 = source.v4;
		v5 = source.v5;
		v6 = source.v6;
		v7 = source.v7;
		v8 = source.v8;

		auxv0 = source.auxv0;
		auxv1 = source.auxv1;
		auxv2 = source.auxv2;
		auxv3 = source.auxv3;
		auxv4 = source.auxv4;
		auxv5 = source.auxv5;
		auxv6 = source.auxv6;
		auxv7 = source.auxv7;
		auxv8 = source.auxv8;

		vel = source.vel;
		mag_vel = source.mag_vel;
		solid_mask = source.solid_mask;
	}
}

LBM2D::~LBM2D()
{
	_clear();
}

bool
LBM2D::_initialize(int x, int y)
{
	resolutionX = x; fieldX = resolutionX - 2;
	resolutionY = y; fieldY = resolutionY - 2;

	v0.ResetDimension(x, y);
	v1.ResetDimension(x, y);
	v2.ResetDimension(x, y);
	v3.ResetDimension(x, y);
	v4.ResetDimension(x, y);
	v5.ResetDimension(x, y);
	v6.ResetDimension(x, y);
	v7.ResetDimension(x, y);
	v8.ResetDimension(x, y);

	auxv0.ResetDimension(x, y);
	auxv1.ResetDimension(x, y);
	auxv2.ResetDimension(x, y);
	auxv3.ResetDimension(x, y);
	auxv4.ResetDimension(x, y);
	auxv5.ResetDimension(x, y);
	auxv6.ResetDimension(x, y);
	auxv7.ResetDimension(x, y);
	auxv8.ResetDimension(x, y);

	vel.ResetDimension(x, y);
	mag_vel.ResetDimension(x, y);
	solid_mask.ResetDimension(x, y);

	// TODO: Other operations if faild return false

	return true;
}

void
LBM2D::_set_sim_params(float rho, float tau)
{
	params.rho = rho;
	params.tau = tau;
}

void
LBM2D::_stream()
{
	if (resolutionX <= 0 || resolutionY <= 0)
		return;

	else{
		for (int i = 1; i <= fieldY; i++){
			for (int j = 1; j <= fieldX; j++){
				auxv1(i, j) = v1(i, j - 1);
				auxv2(i, j) = v2(i - 1, j);
				auxv3(i, j) = v3(i, j + 1);
				auxv4(i, j) = v4(i + 1, j);
				auxv5(i, j) = v5(i - 1, j - 1);
				auxv6(i, j) = v6(i - 1, j + 1);
				auxv7(i, j) = v7(i + 1, j + 1);
				auxv8(i, j) = v8(i + 1, j - 1);
			}
		}

		for (int i = 1; i <= fieldY; i++){
			for (int j = 1; j <= fieldX; j++){
				v1(i, j) = auxv1(i, j);
				v2(i, j) = auxv2(i, j);
				v3(i, j) = auxv3(i, j);
				v4(i, j) = auxv4(i, j);
				v5(i, j) = auxv5(i, j);
				v6(i, j) = auxv6(i, j);
				v7(i, j) = auxv7(i, j);
				v8(i, j) = auxv8(i, j);
			}
		}
	}
}

// BGK Collision algorithm
void
LBM2D::_collide()
{
	float rtau = 1.0f / params.tau, v_sqr = 0.0f;
	float rho, u, v;
	float eq0, eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8;

	float weight0 = 4.0f / 9.0f;	// i = 0
	float weight1 = 1.0f / 9.0f;	// i = 1, 2, 3, 4
	float weight2 = 1.0f / 36.0f;	// i = 5, 6, 7, 8
	float coef0 = 3.0f, coef1 = 4.5f, coef2 = 1.5f;

	for (int i = 1; i <= fieldY; i++){
		for (int j = 1; j <= fieldX; j++){

			rho = v0(i, j) + v1(i, j) + v2(i, j) + v3(i, j) + v4(i, j) + v5(i, j) + v6(i, j) + v7(i, j) + v8(i, j);

			u = (v1(i, j) * params.e1.m_x +
				 v3(i, j) * params.e3.m_x +
				 v5(i, j) * params.e5.m_x + v6(i, j) * params.e6.m_x + v7(i, j) * params.e7.m_x + v8(i, j) * params.e8.m_x) / rho;

			v = (v2(i, j) * params.e2.m_y + v4(i, j) * params.e4.m_y +
				 v5(i, j) * params.e5.m_y + v6(i, j) * params.e6.m_y + v7(i, j) * params.e7.m_y + v8(i, j) * params.e8.m_y) / rho;

			vel(i, j) = VFXEpoch::Vector2Df(u, v);

			// TODO: Update equilibrium function.
			mag_vel(i, j) = sqrt(pow(u, 2) + pow(v, 2));
			v_sqr = pow(u, 2) + pow(v, 2);

			// TODO: Evaluate local equilibrium Bhatnagar-Gross-Krook(BGK) collision
			// Tips for understanding, try to simplify the nested equation.
			eq0 = rho * weight0 * (1.0f - coef2 * v_sqr);
			eq1 = rho * weight1 * (1.0f + coef0 * u + coef1 * pow(u, 2) - coef2 * v_sqr);
			eq2 = rho * weight1 * (1.0f + coef0 * v + coef1 * pow(v, 2) - coef2 * v_sqr);
			eq3 = rho * weight1 * (1.0f - coef0 * u + coef1 * pow(u, 2) - coef2 * v_sqr);
			eq4 = rho * weight1 * (1.0f - coef0 * v + coef1 * pow(v, 2) - coef2 * v_sqr);
			eq5 = rho * weight2 * (1.0f + coef0 * (u + v) + coef1 * pow(u + v, 2) - coef2 * v_sqr);
			eq6 = rho * weight2 * (1.0f + coef0 * (-u + v) + coef1 * pow(-u + v, 2) - coef2 * v_sqr);
			eq7 = rho * weight2 * (1.0f + coef0 * (-u + -v) + coef1 * pow(-u + -v, 2) - coef2 * v_sqr);
			eq8 = rho * weight2 * (1.0f + coef0 * (u - v) + coef1 * pow(u - v, 2) - coef2 * v_sqr);

			// TODO: Local equilibrium update
			v0(i, j) = (1.0f - rtau) * v0(i, j) + eq0;
			v1(i, j) = (1.0f - rtau) * v1(i, j) + eq1;
			v2(i, j) = (1.0f - rtau) * v2(i, j) + eq2;
			v3(i, j) = (1.0f - rtau) * v3(i, j) + eq3;
			v4(i, j) = (1.0f - rtau) * v4(i, j) + eq4;
			v5(i, j) = (1.0f - rtau) * v5(i, j) + eq5;
			v6(i, j) = (1.0f - rtau) * v6(i, j) + eq6;
			v7(i, j) = (1.0f - rtau) * v7(i, j) + eq7;
			v8(i, j) = (1.0f - rtau) * v8(i, j) + eq8;
		}
	}
}

void
LBM2D::_process_boundary()
{

}

/*
 * Flag means "Something" or "Nothing" in the specific cell.
 * This will be pointing to a particular case that sometimes
 * "Nothing" could mean air or other medium.
 */
void
LBM2D::_set_solid_at_cell(BOUNDARY_MASK flag, int x, int y)
{
	assert(x >= 0 && x < solid_mask.getDimX() && y >= 0 && y < solid_mask.getDimY());
	solid_mask(x, y) = flag;
}

void
LBM2D::_bounce_back()
{
	if (resolutionX == 0 || resolutionY == 0) return;

	float v1_prev, v2_prev, v3_prev, v4_prev, v5_prev, v6_prev, v7_prev, v8_prev;

	for (int i = 1; i <= fieldY; i++){
		for (int j = 1; j <= fieldX; j++){
			if (solid_mask(i, j) == VFXEpoch::BOUNDARY_MASK::SOMETHING)
			{
				// v0 doesn't need to be proceesed it will stay on the original position.
				v1_prev = v1(i, j); v2_prev = v2(i, j); v3_prev = v3(i, j);	v4_prev = v4(i, j);
				v5_prev = v5(i, j); v6_prev = v6(i, j);	v7_prev = v7(i, j);	v8_prev = v8(i, j);

				// if the vector --> goes to the solid wall then we just turn back its direction
				// as <--
				v1(i, j) = v3_prev;	v2(i, j) = v4_prev;	v3(i, j) = v1_prev;	v4(i, j) = v2_prev;
				v5(i, j) = v7_prev;	v6(i, j) = v8_prev;	v7(i, j) = v5_prev;	v8(i, j) = v6_prev;
			}
			else
			{
				// TODO��
			}
		}
	}
}

void
LBM2D::_clear()
{
	vel.clear();
	mag_vel.clear();
	solid_mask.clear();

	v0.clear();
	v1.clear();
	v2.clear();
	v3.clear();
	v4.clear();
	v5.clear();
	v6.clear();
	v7.clear();
	v8.clear();

	auxv0.clear();
	auxv1.clear();
	auxv2.clear();
	auxv3.clear();
	auxv4.clear();
	auxv5.clear();
	auxv6.clear();
	auxv7.clear();
}
