/*******************************************************************************
    VFXEPOCH - Physically based simulation VFX

    Copyright (c) 2016 Snow Tsui <trevor.miscellaneous@gmail.com>

    All rights reserved. Use of this source code is governed by
    the MIT license as written in the LICENSE file.
*******************************************************************************/
#ifndef _UTL_GENERAL_H_
#define _UTL_GENERAL_H_

#include <stdlib.h>
#include <iostream>
#include <random>
#include "UTL_Grid.h"

#define REAL double

using namespace std;

namespace VFXEpoch
{
	/*
	* 2D/3D Particle
	* Particles are used for tracking the motion of the Euler grid
	*/
	typedef struct _particle_2d{
		VFXEpoch::Vector2Df pos;
		VFXEpoch::Vector2Df vel;
		VFXEpoch::Vector3Df color;
		float padding;
	}Particle2D;

	typedef struct _particle_3d{
		VFXEpoch::Vector3Df pos;
		VFXEpoch::Vector3Df vel;
		VFXEpoch::Vector3Df color;
		float padding[7];
	}Particles3D;

	enum class VORT_METHODS
	{
		STOKES,
		LEAST_SQUARE,
		RICHARDSON_EXTRAPOLATION
	};

	enum class VECTOR_COMPONENTS
	{
		X, Y, Z
	};

	// Helpers
	float RandomF(float min, float max);
	int RandomI(int min, int max);
	float Lerp(float t, float x0, float x1);
	float Bilerp(float t, float s, float x0, float x1, float y0, float y1);
	void ExtractComponents(VFXEpoch::Grid2DfScalarField& component, VFXEpoch::Grid2DVector2DfField vectorField, VECTOR_COMPONENTS axis);
	void ExtractComponents(VFXEpoch::Grid2DfScalarField& component, VFXEpoch::Grid3DVector3DfField vectorField, VECTOR_COMPONENTS axis);
	void InsertComponents(VFXEpoch::Grid2DfScalarField component, VFXEpoch::Grid2DVector2DfField& vectorField, VECTOR_COMPONENTS axis);
	void InsertComponents(VFXEpoch::Grid2DfScalarField component, VFXEpoch::Grid3DVector3DfField& vectorField, VECTOR_COMPONENTS axis);
	float InterpolateGrid(int N, float x, float y, VFXEpoch::Grid2DfScalarField& field);
	void Zeros(VFXEpoch::Grid2DfScalarField& field);
	void Zeros(VFXEpoch::Grid2DVector2DfField& field);
	void Zeros(VFXEpoch::Grid2DiScalarField& field);
	void Zeros(VFXEpoch::Grid2DVector2DiField& field);

	// Template Helpers
	template <class T>
	void Swap(T& a, T& b)
	{
		T c(a);
		a = b;
		b = c;
	}

	// Get barycentric
	template<class T>
	inline void get_barycentric(T x, int& i, T& f, int i_low, int i_high)
	{
		T s = std::floor(x);
		i = (int)s;
		if (i<i_low){
			i = i_low;
			f = 0;
		}
		else if (i>i_high - 2){
			i = i_high - 2;
			f = 1;
		}
		else
			f = (T)(x - s);
	}

	template<class T>
	const T& _max(const T& a, const T& b)
	{
		return (a < b) ? b : a;
	}

	template<class T>
	const T& _min(const T& a, const T& b)
	{
		return (b < a) ? b : a;
	}

	template <typename T>
	T clamp(const T& n, const T& lower, const T& upper) {
		return VFXEpoch::_max(lower, VFXEpoch::_min(n, upper));
	}
}

#endif