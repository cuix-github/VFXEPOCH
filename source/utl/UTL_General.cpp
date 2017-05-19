/*******************************************************************************
    VFXEPOCH - Physically based simulation VFX

    Copyright (c) 2016 Snow Tsui <trevor.miscellaneous@gmail.com>

    All rights reserved. Use of this source code is governed by
    the MIT license as written in the LICENSE file.
*******************************************************************************/
#include "UTL_General.h"

using namespace VFXEpoch;

float
VFXEpoch::RandomF(float min, float max) {
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<float> range(min, max);
	return range(gen);
}

int
VFXEpoch::RandomI(int min, int max)
{
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<int> range(min, max);
	return range(gen);
}

float
VFXEpoch::Lerp(float t, float x0, float x1) {
	return (1.0f - t) * x0 + t * x1;
}

float
VFXEpoch::Bilerp(float t, float s, float x0, float x1, float y0, float y1) {
	return VFXEpoch::Lerp(s, VFXEpoch::Lerp(t, x0, x1), VFXEpoch::Lerp(t, y0, y1));
}

void
VFXEpoch::ExtractComponents(VFXEpoch::Grid2DfScalarField& component, VFXEpoch::Grid2DVector2DfField vectorField, VECTOR_COMPONENTS axis) {
	if (component.getDimY() != vectorField.getDimY() ||
		component.getDimX() != vectorField.getDimX()) {
		assert(component.getDimY() == vectorField.getDimY() && component.getDimX() == vectorField.getDimX());
	}
	else {
		switch (axis)	{
		case VFXEpoch::VECTOR_COMPONENTS::X:
			for (int i = 0; i != component.getDimY(); i++){
				for (int j = 0; j != component.getDimX(); j++){
					component.setData(vectorField.getData(i, j).m_x, i, j);
				}
			}
			break;
		case VFXEpoch::VECTOR_COMPONENTS::Y:
			for (int i = 0; i != component.getDimY(); i++){
				for (int j = 0; j != component.getDimX(); j++){
					component.setData(vectorField.getData(i, j).m_y, i, j);
				}
			}
			break;
		case VFXEpoch::VECTOR_COMPONENTS::Z:
			break;
		default:
			for (int i = 0; i != component.getDimY(); i++){
				for (int j = 0; j != component.getDimX(); j++){
					component.setData(vectorField.getData(i, j).m_x, i, j);
				}
			}
			break;
		}
	}
}

void
VFXEpoch::ExtractComponents(VFXEpoch::Grid2DfScalarField& component, VFXEpoch::Grid3DVector3DfField vectorField, VECTOR_COMPONENTS axis) {
	// TODO: Coming soon
}

void
VFXEpoch::InsertComponents(VFXEpoch::Grid2DfScalarField component, VFXEpoch::Grid2DVector2DfField& vectorField, VECTOR_COMPONENTS axis)
{
	if (component.getDimY() != vectorField.getDimY() ||
		component.getDimX() != vectorField.getDimX())
	{
		assert(component.getDimY() == vectorField.getDimY() && component.getDimX() != vectorField.getDimX());
	}
	else
	{
		VFXEpoch::Vector2Df vec(0.0f, 0.0f);
		switch (axis)
		{
		case VFXEpoch::VECTOR_COMPONENTS::X:
			for (int i = 0; i != component.getDimY(); i++){
				for (int j = 0; j != component.getDimX(); j++){
					vec = vectorField.getData(i, j);
					vec.m_x = component.getData(i, j);
					vectorField.setData(vec, i, j);
				}
			}
			break;
		case VFXEpoch::VECTOR_COMPONENTS::Y:
			for (int i = 0; i != component.getDimY(); i++){
				for (int j = 0; j != component.getDimX(); j++){
					vec = vectorField.getData(i, j);
					vec.m_y = component.getData(i, j);
					vectorField.setData(vec, i, j);
				}
			}
			break;
		case VFXEpoch::VECTOR_COMPONENTS::Z:
			break;
		default:
			break;
		}
	}
}

void
VFXEpoch::InsertComponents(VFXEpoch::Grid2DfScalarField component, VFXEpoch::Grid3DVector3DfField& vectorField, VECTOR_COMPONENTS axis) {
	// TODO: Coming soon
}

float
VFXEpoch::InterpolateGrid(float x, float y, VFXEpoch::Grid2DfScalarField& field)
{
	int i, j;
	float fx, fy;

	VFXEpoch::get_barycentric(x, j, fx, 0, field.getDimX());
	VFXEpoch::get_barycentric(y, i, fy, 0, field.getDimY());
	return VFXEpoch::Bilerp(fx, fy, field(i, j), field(i + 1, j), field(i, j + 1), field(i + 1, j + 1));
}

float
VFXEpoch::InterpolateGrid(Vector2Df pos, Grid2DfScalarField& field){
	int i, j;
	float fx, fy;
	VFXEpoch::get_barycentric(pos.m_x, j, fx, 0, field.getDimX());
	VFXEpoch::get_barycentric(pos.m_y, i, fy, 0, field.getDimY());
	return VFXEpoch::Bilerp(fx, fy, field(i, j), field(i + 1, j), field(i, j + 1), field(i + 1, j + 1));
}

float 
VFXEpoch::InterpolateGradient(Vector2Df gradient, Vector2Df pos, VFXEpoch::Grid2DfScalarField& field){
	int i, j;
	float fx, fy;
	VFXEpoch::get_barycentric(pos.m_x, j, fx, 0, field.getDimX());
	VFXEpoch::get_barycentric(pos.m_y, i, fy, 0, field.getDimY());

	assert(i >= 0 && i < field.getDimY() && j >= 0 && j < field.getDimX());
	float dy0 = field(i+1, j) - field(i, j);
	float dy1 = field(i+1, j+1) - field(i, j+1);
	float dx0 = field(i, j+1) - field(i, j);
	float dx1 = field(i+1, j+1) - field(i+1, j);

	gradient.m_x = VFXEpoch::Lerp(fy, dx0, dx1);
	gradient.m_y = VFXEpoch::Lerp(fx, dy0, dy1);
}

float
VFXEpoch::InteralFrac(float left, float right){
	if(left < 0 && right < 0)
      return 1.0f;
   	if (left < 0 && right >= 0)
      return left / (left - right);
  	if(left >= 0 && right < 0)
      return right / (right - left);
   	else
      return 0;
}

void
VFXEpoch::Zeros(VFXEpoch::Grid2DfScalarField& field)
{
	for (int i = 0; i != field.getDimY(); i++)
	{
		for (int j = 0; j != field.getDimX(); j++)
		{
			field(i, j) = 0.0f;
		}
	}
}

void
VFXEpoch::Zeros(VFXEpoch::Grid2DVector2DfField& field)
{
	for (int i = 0; i != field.getDimY(); i++)
	{
		for (int j = 0; j != field.getDimX(); j++)
		{
			field(i, j) = VFXEpoch::Vector2Df(0.0f, 0.0f);
		}
	}
}

void
VFXEpoch::Zeros(VFXEpoch::Grid2DiScalarField& field)
{
	for (int i = 0; i != field.getDimY(); i++){
		for (int j = 0; j != field.getDimX(); j++){
			field(i, j) = 0;
		}
	}
}

void
VFXEpoch::Zeros(VFXEpoch::Grid2DVector2DiField& field)
{
	for (int i = 0; i != field.getDimY(); i++){
		for (int j = 0; j != field.getDimX(); j++){
			field(i, j) = VFXEpoch::Vector2Di(0, 0);
		}
	}
}
