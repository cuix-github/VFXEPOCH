/*******************************************************************************
    VFXEPOCH - Physically based simulation VFX

    Copyright (c) 2016 Snow Tsui <trevor.miscellaneous@gmail.com>

    All rights reserved. Use of this source code is governed by
    the MIT license as written in the LICENSE file.
*******************************************************************************/
#include "SIM_Mac.h"

VFXEpoch::Mac2D::Mac2D(){

}

VFXEpoch::Mac2D::Mac2D(int N, float _dt, float _kappa) :dimenssion(N), dt(_dt), kappa(_kappa){
	curlField.ResetDimension(N + 1, N + 1);
	vel_u.ResetDimension(N, N + 1);
	vel_v.ResetDimension(N + 1, N);
	denField.ResetDimension(N, N);
	pressure.ResetDimension(N, N);
}

VFXEpoch::Mac2D::Mac2D(const Mac2D& source){
	dimenssion = source.dimenssion;
	dt = source.dt;
	kappa = source.kappa;

	curlField = source.curlField;
	vel_u = source.vel_u;
	vel_v = source.vel_v;
	denField = source.denField;
}

VFXEpoch::Mac2D& VFXEpoch::Mac2D::operator=(const VFXEpoch::Mac2D& source){
	dimenssion = source.dimenssion;
	dt = source.dt;
	kappa = source.kappa;

	curlField = source.curlField;
	vel_u = source.vel_u;
	vel_v = source.vel_v;
	denField = source.denField;
	return *this;
}

VFXEpoch::Mac2D::~Mac2D(){
	dimenssion = 0;
	dt = 0.0f;
	kappa = 0.0f;
	_clear();
}

bool VFXEpoch::Mac2D::InitGasSimParameters(int N, float _dt, float _kappa){
	bool result;
	result = _init_gas_sim_parameters(N, _dt, _kappa);
	if (!result)
		return false;

	return true;
}

void VFXEpoch::Mac2D::SetField(FIELD2D field, VFXEpoch::Grid2DfScalarField fieldData){
	_set_field(field, fieldData);
}

VFXEpoch::Grid2DfScalarField VFXEpoch::Mac2D::GetField(FIELD2D field){
	return _get_field(field);
}

void VFXEpoch::Mac2D::Reset(){
	dimenssion = 0;
	dt = 0.0f;
	kappa = 0.0f;
	_clear();
}

void VFXEpoch::Mac2D::GetFeatures(){
	_features_calculation();
}

bool VFXEpoch::Mac2D::_init_gas_sim_parameters(int N, float _dt, float _kappa){
	dimenssion = N;
	dt = _dt;
	kappa = _kappa;

	return true;
}

void VFXEpoch::Mac2D::_set_field(FIELD2D field, VFXEpoch::Grid2DfScalarField fieldData){
	switch (field)
	{
	case VFXEpoch::Mac2D::FIELD2D::CURL:
		curlField = fieldData;
		break;
	case VFXEpoch::Mac2D::FIELD2D::U_VEL:
		vel_u = fieldData;
		break;
	case VFXEpoch::Mac2D::FIELD2D::V_VEL:
		vel_v = fieldData;
		break;
	case VFXEpoch::Mac2D::FIELD2D::DENSITY:
		denField = fieldData;
		break;
	case VFXEpoch::Mac2D::FIELD2D::PRESSURE:
		pressure = fieldData;
		break;
	default:
		break;
	}
}

VFXEpoch::Grid2DfScalarField VFXEpoch::Mac2D::_get_field(FIELD2D field)
{
	switch (field)
	{
	case VFXEpoch::Mac2D::FIELD2D::CURL:
		return curlField;
		break;
	case VFXEpoch::Mac2D::FIELD2D::U_VEL:
		return vel_u;
		break;
	case VFXEpoch::Mac2D::FIELD2D::V_VEL:
		return vel_v;
		break;
	case VFXEpoch::Mac2D::FIELD2D::DENSITY:
		return denField;
		break;
	case VFXEpoch::Mac2D::FIELD2D::PRESSURE:
		return pressure;
		break;
	default:
		break;
	}
}

void VFXEpoch::Mac2D::_features_calculation(){
	_get_curls();
}

void VFXEpoch::Mac2D::_clear(){
	curlField.clear();
	vel_u.clear();
	vel_v.clear();
	denField.clear();
}

void VFXEpoch::Mac2D::_get_curls(){
	VFXEpoch::Analysis::computeCurl_mac(curlField, vel_u, vel_v);
}
