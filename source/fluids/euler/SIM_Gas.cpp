/*******************************************************************************
    VFXEPOCH - Physically based simulation VFX

    Copyright (c) 2016 Snow Tsui <trevor.miscellaneous@gmail.com>

    All rights reserved. Use of this source code is governed by
    the MIT license as written in the LICENSE file.
*******************************************************************************/
#include "SIM_Gas.h"

using namespace VFXEpoch::Solvers;

SL2D::SL2D(){
	dimX = dimY = 0;
	iterations = 0;
	viscosity = 0.0;
	spacingX = 0.0;
	spacingY = 0.0;
	timeStep = 0.0;
	diffuseRate = 0.0f;
}

SL2D::SL2D(int dimX, int dimY, int iterations, double timeStep, float diffuseRate, float viscosity, float sourceRate, float spacingX, float spacingY,
	VFXEpoch::Grid2DVector2DfField velField, VFXEpoch::Grid2DVector2DfField velFieldPrev,
	VFXEpoch::Grid2DfScalarField denField, VFXEpoch::Grid2DfScalarField denFieldPrev, VFXEpoch::Grid2DfScalarField pressure,
	VFXEpoch::Grid2DfScalarField divergence, VFXEpoch::Grid2DVector2DfField gradPressure, VFXEpoch::Grid2DfScalarField temperature, VFXEpoch::Grid2DfScalarField temperaturePrev){
	this->dimX = dimX;
	this->dimY = dimY;
	this->timeStep = timeStep;
	this->diffuseRate = diffuseRate;
	this->sourceRate = sourceRate;
	this->spacingX = spacingX;
	this->spacingY = spacingY;
	this->velocityField = velField; velocityFieldPrev = velFieldPrev;
	this->densityField = denField; densityFieldPrev = denFieldPrev;
	this->temperature = temperature;
	this->temperaturePrev = temperaturePrev;
	this->pressure = pressure;
	this->divergence = divergence;
	this->gradPressure = gradPressure;
}

SL2D::SL2D(const SL2D& source){
	this->dimX = source.dimX; dimY = source.dimY;
	this->timeStep = source.timeStep;
	this->diffuseRate = source.diffuseRate;
	this->sourceRate = source.sourceRate;
	this->spacingX = source.spacingX; spacingY = source.spacingY;
	this->velocityField = source.velocityField;
	this->velocityFieldPrev = source.velocityFieldPrev;
	this->densityField = source.densityField;
	this->densityFieldPrev = source.densityFieldPrev;
	this->temperature = source.temperature;
	this->temperaturePrev = source.temperaturePrev;
	this->pressure = pressure;
	this->divergence = divergence;
	this->gradPressure = gradPressure;
}

SL2D& SL2D::operator=(const SL2D& source){
	this->dimX = source.dimX; dimY = source.dimY;
	this->timeStep = source.timeStep;
	this->diffuseRate = source.diffuseRate;
	this->sourceRate = source.sourceRate;
	this->spacingX = source.spacingX; spacingY = source.spacingY;
	this->velocityField = source.velocityField;
	this->velocityFieldPrev = source.velocityFieldPrev;
	this->densityField = source.densityField;
	this->densityFieldPrev = source.densityFieldPrev;
	this->temperature = source.temperature;
	this->temperaturePrev = source.temperaturePrev;
	this->pressure = pressure;
	this->divergence = divergence;
	this->gradPressure = gradPressure;
	return *this;
}

SL2D::~SL2D(){
	_clear();
}

bool
SL2D::Initialize(VFXEpoch::Vector2Di dimension, VFXEpoch::Vector2Df spacing, int iterations, float timeStep, float diffuseRate, float viscosity, float sourceRate)
{
	bool result;
	result = _init(dimension.m_y, dimension.m_x, iterations, timeStep, diffuseRate, viscosity, sourceRate, spacing.m_x, spacing.m_y);
	if (!result)
		return false;

	return true;
}

void
SL2D::SetField(VFXEpoch::Grid2DfScalarField field, VFXEpoch::COMPUTATIONAL_SCALAR_FIELD_2D fieldName)
{
	_set_field(field, fieldName);
}

void
SL2D::SetField(VFXEpoch::Grid2DVector2DfField field, VFXEpoch::COMPUTATIONAL_VECTOR_FIELD_2D fieldName)
{
	_set_field(field, fieldName);
}

void
SL2D::SetField(VFXEpoch::Mac2D mac)
{
	_set_field(mac);
}

void
SL2D::SetFieldBoundary(VFXEpoch::BOUNDARY boundaryType, VFXEpoch::EDGES_2DSIM edge)
{
	switch (edge)
	{
	case VFXEpoch::EDGES_2DSIM::TOP:
		fieldBoundaries[0].side = edge;
		fieldBoundaries[0].boundaryType = boundaryType;
		break;
	case VFXEpoch::EDGES_2DSIM::BOTTOM:
		fieldBoundaries[1].side = edge;
		fieldBoundaries[1].boundaryType = boundaryType;
		break;
	case VFXEpoch::EDGES_2DSIM::LEFT:
		fieldBoundaries[2].side = edge;
		fieldBoundaries[2].boundaryType = boundaryType;
		break;
	case VFXEpoch::EDGES_2DSIM::RIGHT:
		fieldBoundaries[3].side = edge;
		fieldBoundaries[3].boundaryType = boundaryType;
		break;
	default:
		break;
	}
}

void
SL2D::AddSource(VFXEpoch::Grid2DfScalarField& target, VFXEpoch::Grid2DfScalarField source)
{
	_set_source(target, source);
}

void
SL2D::AddSource(VFXEpoch::Grid2DVector2DfField& target, VFXEpoch::Grid2DVector2DfField source)
{
	_set_source(target, source);
}

// Q:Is it worth to advect
// 1.Previous timestep density field?
// 2.Pressure?
// 3.Divergence?
// A:
// So far there is no proof showing that we need to
// diffuse those 3 fields. But in the future, the curl
// could be one.
void
SL2D::Diffuse(VFXEpoch::Grid2DfScalarField& target, VFXEpoch::Grid2DfScalarField source)
{
	_diffuse(target, source);
}

// Comments:
// Same question as previous one.
void
SL2D::Diffuse(VFXEpoch::Grid2DVector2DfField& target, VFXEpoch::Grid2DVector2DfField source)
{
	_diffuse(target, source);
}

void
SL2D::Advect(VFXEpoch::Grid2DfScalarField& target, VFXEpoch::Grid2DfScalarField previous_timestep_field, VFXEpoch::Grid2DVector2DfField ref)
{
	_advect(target, previous_timestep_field, ref);
}

void
SL2D::Advect(VFXEpoch::Grid2DVector2DfField& target, VFXEpoch::Grid2DVector2DfField previous_timestep_field, VFXEpoch::Grid2DVector2DfField ref)
{
	_advect(target, previous_timestep_field, ref);
}

void
SL2D::PressureSolve(VFXEpoch::Grid2DVector2DfField& target, VFXEpoch::Grid2DfScalarField pressure, VFXEpoch::Grid2DfScalarField divergence)
{
	_project(target, pressure, divergence);
}

void
SL2D::AddVortConf(VFXEpoch::Grid2DVector2DfField& target, float vort_conf_eps, VFXEpoch::VORT_METHODS method){
	_vort_conf_update(target, vort_conf_eps, method);
}

void
SL2D::Reset()
{
	densityField.zeroScalars();
	densityFieldPrev.zeroScalars();
	velocityField.zeroVectors();
	velocityField.zeroVectors();
	pressure.zeroScalars();
	divergence.zeroScalars();
}

void
SL2D::Shutdown()
{
	_clear();
}

const VFXEpoch::Grid2DfScalarField&
SL2D::getDensityField() const {
	return this->densityField;
}

const VFXEpoch::Grid2DfScalarField&
SL2D::getDensityPrevField() const {
	return this->densityFieldPrev;
}

const VFXEpoch::Grid2DVector2DfField&
SL2D::getVelocityField() const {
	return this->velocityField;
}

const VFXEpoch::Grid2DVector2DfField&
SL2D::getVelocityPrevField() const {
	return this->velocityFieldPrev;
}

const VFXEpoch::BndConditionPerEdge* SL2D::getFieldBoundaries() const
{
	return this->fieldBoundaries;
}

VFXEpoch::Grid2DfScalarField& SL2D::getDensityField() {
	return this->densityField;
}

VFXEpoch::Grid2DfScalarField& SL2D::getDensityPrevField() {
	return this->densityFieldPrev;
}

VFXEpoch::Grid2DfScalarField& SL2D::getTemperatureField() {
	return this->temperature;
}
VFXEpoch::Grid2DfScalarField& SL2D::getTemperaturePrevField() {
	return this->temperaturePrev;
}

VFXEpoch::Grid2DVector2DfField& SL2D::getVelocityField() {
	return this->velocityField;
}

VFXEpoch::Grid2DVector2DfField& SL2D::getVelocityPrevField() {
	return this->velocityFieldPrev;
}

VFXEpoch::BndConditionPerEdge* SL2D::getFieldBoundaries()
{
	return this->fieldBoundaries;
}

bool
SL2D::_init(int dimY, int dimX, int iterations, double timeStep, float diffuseRate, float viscosity, float sourceRate, float spacingX, float spacingY)
{
	this->dimX = dimX; this->dimY = dimY;
	this->timeStep = timeStep;
	this->diffuseRate = diffuseRate;
	this->sourceRate = sourceRate;
	this->spacingX = spacingX; this->spacingY = spacingY;
	this->iterations = iterations;
	this->viscosity = viscosity;

	velocityFieldPrev.ResetDimension(dimY, dimX);
	velocityField.ResetDimension(dimY, dimX);
	gradPressure.ResetDimension(dimY, dimX);
	pressure.ResetDimension(dimY, dimX);
	divergence.ResetDimension(dimY, dimX);
	densityFieldPrev.ResetDimension(dimY, dimX);
	densityField.ResetDimension(dimY, dimX);

	this->velocityFieldPrev.zeroVectors();
	this->velocityField.zeroVectors();
	this->densityFieldPrev.zeroScalars();
	this->densityField.zeroScalars();

	return true;
}

void
SL2D::_set_field(VFXEpoch::Grid2DfScalarField field, VFXEpoch::COMPUTATIONAL_SCALAR_FIELD_2D fieldName)
{
	switch (fieldName)
	{
	case VFXEpoch::COMPUTATIONAL_SCALAR_FIELD_2D::DENSITY:
		this->densityField = field;
		break;
	case VFXEpoch::COMPUTATIONAL_SCALAR_FIELD_2D::DENSITY_PREV:
		this->densityFieldPrev = field;
		break;
	case VFXEpoch::COMPUTATIONAL_SCALAR_FIELD_2D::TEMPERATURE:
		this->temperature = field;
		break;
	case VFXEpoch::COMPUTATIONAL_SCALAR_FIELD_2D::TEMPERATURE_PREV:
		this->temperaturePrev = field;
		break;
	case VFXEpoch::COMPUTATIONAL_SCALAR_FIELD_2D::PRESSURE:
		this->pressure = field;
		break;
	case VFXEpoch::COMPUTATIONAL_SCALAR_FIELD_2D::DIVERGENCE:
		this->divergence = field;
		break;
	default:
		field.clear();
		break;
	}
}

void
SL2D::_set_field(VFXEpoch::Grid2DVector2DfField field, VFXEpoch::COMPUTATIONAL_VECTOR_FIELD_2D fieldName)
{
	switch (fieldName)
	{
	case VFXEpoch::COMPUTATIONAL_VECTOR_FIELD_2D::VEL_U:
		// TODO: Will be executed when staggered grid applied
		break;
	case VFXEpoch::COMPUTATIONAL_VECTOR_FIELD_2D::VEL_V:
		// TODO: Will be executed when staggered grid applied
		break;
	case VFXEpoch::COMPUTATIONAL_VECTOR_FIELD_2D::VEL_U_PREV:
		// TODO: Will be executed when staggered grid applied
		break;
	case VFXEpoch::COMPUTATIONAL_VECTOR_FIELD_2D::VEL_V_PREV:
		// TODO: Will be executed when staggered grid applied
		break;

		// Central difference do not need to process velocity component separately
	case VFXEpoch::COMPUTATIONAL_VECTOR_FIELD_2D::VEL:
		this->velocityField = field;
		break;
	case VFXEpoch::COMPUTATIONAL_VECTOR_FIELD_2D::VEL_PREV:
		this->velocityFieldPrev = field;
		break;

		// Pressure gradient field is used for velocity correction
	case VFXEpoch::COMPUTATIONAL_VECTOR_FIELD_2D::GRAD_PRESSURE:
		this->gradPressure = field;
		break;
	default:
		break;
	}
}

void
SL2D::_set_field(VFXEpoch::Mac2D mac)
{
	staggeredGrid = mac;
}

void
SL2D::_set_source(VFXEpoch::Grid2DfScalarField& scalarField, VFXEpoch::Grid2DfScalarField source)
{
	for (int i = 0; i != scalarField.m_yCell; i++)
	{
		for (int j = 0; j != scalarField.m_xCell; j++)
			scalarField.setData(scalarField(i, j) + source(i, j) * (float)timeStep * sourceRate, i, j);
	}
}

void
SL2D::_set_source(VFXEpoch::Grid2DVector2DfField& vectorField, VFXEpoch::Grid2DVector2DfField source)
{
	for (int i = 0; i != vectorField.m_yCell; i++)
	{
		for (int j = 0; j != vectorField.m_xCell; j++)
			vectorField.setData(vectorField(i, j) + source(i, j) * (float)timeStep * sourceRate, i, j);
	}
}

void
SL2D::_set_field_boundary(VFXEpoch::Grid2DfScalarField& field)
{
	for (int i = 0; i != 4; i++)
	{
		field.setBoundaries(fieldBoundaries[i].boundaryType, fieldBoundaries[i].side);
	}

	field.setBoundariesOnCorners();
}

void
SL2D::_set_field_boundary(VFXEpoch::Grid2DVector2DfField& field)
{
	for (int i = 0; i != 4; i++)
	{
		field.setBoundaries(fieldBoundaries[i].boundaryType, fieldBoundaries[i].side);
	}

	field.setBoundariesOnCorners();
}

VFXEpoch::Vector2Df
SL2D::Runge_Kutta_II(int dimX, int dimY, VFXEpoch::Grid2DfScalarField& u, VFXEpoch::Grid2DfScalarField& v, const VFXEpoch::Vector2Df& position, float dt)
{
	// TODO: Runge-Kutta II order ODEs integration
	VFXEpoch::Vector2Df vel = this->_get_vel(dimX, dimY, position, u, v);
	vel = this->_get_vel(dimX, dimY, VFXEpoch::Vector2Df(position.m_x + 0.5f * dt * vel.m_x, position.m_y + 0.5f * dt * vel.m_y), u, v);
	return VFXEpoch::Vector2Df(position.m_x + dt * vel.m_x, position.m_y + dt * vel.m_y);
}

void
SL2D::_advect(VFXEpoch::Grid2DVector2DfField& vectorField, VFXEpoch::Grid2DVector2DfField vectorFieldOirigin, VFXEpoch::Grid2DVector2DfField referenceField)
{
	float dtx = (float)timeStep * (vectorFieldOirigin.getDimX() - 2);
	float dty = (float)timeStep * (vectorFieldOirigin.getDimY() - 2);
	float exactx, exacty;
	float xlerp, ylerp, subLerp0, subLerp1;
	int nx, ny;
	int idxi, idxj;
	nx = vectorField.getDimX() - 2;
	ny = vectorField.getDimY() - 2;
	VFXEpoch::Vector2Df data;

	for (int i = 1; i <= ny; i++){
		for (int j = 1; j <= nx; j++){
			exactx = i - dtx * referenceField(i, j).m_x;
			exacty = j - dty * referenceField(i, j).m_y;

			if (exactx < 0.5f) exactx = 0.5f;
			if (exacty < 0.5f) exacty = 0.5f;
			if (exactx > (nx + 0.5f)) exactx = nx + 0.5f;
			if (exacty > (ny + 0.5f)) exacty = ny + 0.5f;

			idxi = (int)exactx;
			idxj = (int)exacty;

			// Lerp by aixs order
			// Lerp y direction first and then take the results as
			// inputs to do the x direction lerp or vice versa.
			xlerp = exactx - idxi;
			ylerp = exacty - idxj;

			// Lerp x direction first
			subLerp0 = VFXEpoch::Lerp(ylerp, vectorFieldOirigin(idxi, idxj).m_x, vectorFieldOirigin(idxi, idxj + 1).m_x);
			subLerp1 = VFXEpoch::Lerp(ylerp, vectorFieldOirigin(idxi + 1, idxj).m_x, vectorFieldOirigin(idxi + 1, idxj + 1).m_x);
			data.m_x = VFXEpoch::Lerp(xlerp, subLerp0, subLerp1);
			subLerp0 = VFXEpoch::Lerp(ylerp, vectorFieldOirigin(idxi, idxj).m_y, vectorFieldOirigin(idxi, idxj + 1).m_y);
			subLerp1 = VFXEpoch::Lerp(ylerp, vectorFieldOirigin(idxi + 1, idxj).m_y, vectorFieldOirigin(idxi + 1, idxj + 1).m_y);
			data.m_y = VFXEpoch::Lerp(xlerp, subLerp0, subLerp1);
			vectorField.setData(data, i, j);
		}
	}
	_set_field_boundary(vectorField);
}

void
SL2D::_advect(VFXEpoch::Grid2DfScalarField& scalarField, VFXEpoch::Grid2DfScalarField scalarFieldOrigin, VFXEpoch::Grid2DVector2DfField referenceField)
{
	float dtx, dty;
	float exactx, exacty;
	float xlerp, ylerp, sublerp0, sublerp1;
	int idxi, idxj;
	int nx, ny;
	nx = scalarField.getDimX() - 2;
	ny = scalarField.getDimY() - 2;
	dtx = (float)timeStep * (float)nx;
	dty = (float)timeStep * (float)ny;

	for (int i = 1; i <= ny; i++)
	{
		for (int j = 1; j <= nx; j++)
		{
			exactx = i - dtx * referenceField(i, j).m_x;
			exacty = j - dty * referenceField(i, j).m_y;

			if (exactx < 0.5f) exactx = 0.5f;
			if (exacty < 0.5f) exacty = 0.5f;
			if (exactx > nx + 0.5f) exactx = nx + 0.5f;
			if (exacty > ny + 0.5f) exacty = ny + 0.5f;

			idxi = (int)exactx;
			idxj = (int)exacty;

			xlerp = exactx - idxi;
			ylerp = exacty - idxj;

			sublerp0 = VFXEpoch::Lerp(xlerp, scalarFieldOrigin(idxi, idxj), scalarFieldOrigin(idxi + 1, idxj));
			sublerp1 = VFXEpoch::Lerp(xlerp, scalarFieldOrigin(idxi, idxj + 1), scalarFieldOrigin(idxi + 1, idxj + 1));
			scalarField.setData(VFXEpoch::Lerp(ylerp, sublerp0, sublerp1), i, j);
		}
	}
	_set_field_boundary(scalarField);
}

void
SL2D::_advect_rk2()
{
	// TODO: RK II based advection
}

void
SL2D::_advect_particles_rk2(VFXEpoch::Grid2DfScalarField& u, VFXEpoch::Grid2DfScalarField& v, std::vector<VFXEpoch::Particle2D>& particleContainer)
{
	VFXEpoch::Vector2Df position(0.0f, 0.0f), vel(0.0f, 0.0f);
	// TODO: Move particles along with velocity field
	for (std::vector<VFXEpoch::Particle2D>::iterator ite = particleContainer.begin(); ite != particleContainer.end(); ite++){
		position.ZeroComponents();
		position = this->Runge_Kutta_II(dimX, dimY, u, v, VFXEpoch::Vector2Df(ite->pos.m_x, ite->pos.m_y), (float)this->timeStep);
		ite->pos.m_x = position.m_x; ite->pos.m_y = position.m_y;
		vel = this->_get_vel(dimX, dimY, VFXEpoch::Vector2Df(ite->pos.m_x, ite->pos.m_y), u, v);
		ite->vel = vel;
	}

	// TODO: Correct nodal-solid boundary conditions.
}

void
SL2D::_diffuse(VFXEpoch::Grid2DfScalarField& densityField, VFXEpoch::Grid2DfScalarField densityFieldOrigin)
{
	int size = (densityFieldOrigin.getDimX() - 2) * (densityFieldOrigin.getDimY() - 2);
	float coef = (float)timeStep * diffuseRate * size;
	VFXEpoch::LinearSolver::GSSolve(densityField, densityFieldOrigin, fieldBoundaries, coef, 1 + 4 * coef, iterations);
}

void
SL2D::_diffuse(VFXEpoch::Grid2DVector2DfField& vectorField, VFXEpoch::Grid2DVector2DfField vectorFieldOrigin)
{
	VFXEpoch::Grid2DfScalarField x(vectorFieldOrigin.getDimY(), vectorFieldOrigin.getDimX(), vectorFieldOrigin.getDx(), vectorFieldOrigin.getDy());
	VFXEpoch::Grid2DfScalarField xPrev(vectorFieldOrigin.getDimY(), vectorFieldOrigin.getDimX(), vectorFieldOrigin.getDx(), vectorFieldOrigin.getDy());
	VFXEpoch::Grid2DfScalarField y(vectorFieldOrigin.getDimY(), vectorFieldOrigin.getDimX(), vectorFieldOrigin.getDx(), vectorFieldOrigin.getDy());
	VFXEpoch::Grid2DfScalarField yPrev(vectorFieldOrigin.getDimY(), vectorFieldOrigin.getDimX(), vectorFieldOrigin.getDx(), vectorFieldOrigin.getDy());
	ExtractComponents(xPrev, vectorFieldOrigin, VFXEpoch::VECTOR_COMPONENTS::X);
	ExtractComponents(yPrev, vectorFieldOrigin, VFXEpoch::VECTOR_COMPONENTS::Y);
	this->_diffuse(x, xPrev);
	this->_diffuse(y, yPrev);
	InsertComponents(x, vectorField, VFXEpoch::VECTOR_COMPONENTS::X);
	InsertComponents(y, vectorField, VFXEpoch::VECTOR_COMPONENTS::Y);
}

void
SL2D::_project(VFXEpoch::Grid2DVector2DfField& field, VFXEpoch::Grid2DfScalarField pressure, VFXEpoch::Grid2DfScalarField divergence)
{
	double h = 1.0f / (pressure.getDimY() - 1);
	VFXEpoch::Analysis::computeDivergence_uniform(divergence, field);
	_set_field_boundary(divergence);
	divergence.scale(-1.0f);
	VFXEpoch::Zeros(pressure);
	_set_field_boundary(pressure);
	VFXEpoch::LinearSolver::GSSolve(pressure, divergence, this->fieldBoundaries, 1, 4, this->iterations);
	//VFXEpoch::LinearSolver::MultigridSolve_V_Cycle((float)h, pressure, divergence, this->fieldBoundaries, 1, 4, 30);
	VFXEpoch::Analysis::computeGradient_uniform(gradPressure, pressure);
	_update_vel(field, gradPressure);
}

void
SL2D::_update_vel(VFXEpoch::Grid2DVector2DfField& dest, VFXEpoch::Grid2DVector2DfField _gradPressure)
{
	dest -= _gradPressure;
	_set_field_boundary(dest);
}

void
SL2D::_get_uniform_curls(VFXEpoch::Grid2DfScalarField& w, VFXEpoch::Grid2DVector2DfField vel)
{
	VFXEpoch::Analysis::computeCurl_uniform(w, vel);
	_set_field_boundary(w);
}

void
SL2D::_get_buoyancy(VFXEpoch::Grid2DfScalarField density, VFXEpoch::Grid2DfScalarField temp, VFXEpoch::Grid2DVector2DfField& buoyancy, float alpha, float beta)
{
	float sampledDensity(0.0f), sampledTemperature(0.0f);
	for (int i = 1; i != buoyancy.getDimY() - 1; i++) {
		for (int j = 1; j != buoyancy.getDimX() - 1; j++) {
			sampledDensity = density(i, j);
			sampledTemperature = temp(i, j);
			buoyancy(i, j).m_x = 0.0f;
			buoyancy(i, j).m_y = -alpha * sampledDensity + beta * sampledTemperature;
		}
	}
}

void
SL2D::_vort_conf_update(VFXEpoch::Grid2DVector2DfField& velocity, float vort_conf_eps, VFXEpoch::VORT_METHODS method)
{
	VFXEpoch::Grid2DVector2DfField fvort(velocity.getDimX(), velocity.getDimY(), 1.0f / (velocity.getDimX() - 2), 1.0f / (velocity.getDimY() - 2));
	VFXEpoch::Grid2DfScalarField vort(velocity.getDimX(), velocity.getDimY(), 1.0f / (velocity.getDimX() - 2), 1.0f / (velocity.getDimY() - 2));
	VFXEpoch::Grid2DVector2DfField vort_grad(velocity.getDimX(), velocity.getDimY(), 1.0f / (velocity.getDimX() - 2), 1.0f / (velocity.getDimY() - 2));
	VFXEpoch::Vector3Df vort_dir(0.0f, 0.0f, 0.0f);
	VFXEpoch::Vector3Df vorticity(0.0f, 0.0f, 0.0f);
	VFXEpoch::Vector3Df fconf(0.0f, 0.0f, 0.0f);
	VFXEpoch::Zeros(fvort);
	VFXEpoch::Zeros(vort);

	// Users' choice to pickup an appropriate method for vorticity calculation
	// 1.LEAST_SQUARE:
	// Reduces the effects of fluctuations
	// 2.RICHARDSON_EXTRAPOLATION:
	// Smaller truncation error
	if (method == VFXEpoch::VORT_METHODS::STOKES){
		VFXEpoch::Analysis::computeCurl_uniform_Stokes(vort, velocity);
	}
	else if (method == VFXEpoch::VORT_METHODS::LEAST_SQUARE){
		VFXEpoch::Analysis::computeCurl_uniform_Richardson(vort, velocity);
	}
	else{
		VFXEpoch::Analysis::computeCurl_uniform_Richardson(vort, velocity);
	}

	// Get gradient of vorticity
	VFXEpoch::Analysis::computeGradient_uniform(vort_grad, vort);

	for (int i = 1; i != fvort.getDimY() - 1; i++){
		for (int j = 1; j != fvort.getDimX() - 1; j++){
			float gradu(0.0f), gradv(0.0f), gradlen(0.0f);
			gradu = vort_grad(i, j).m_x;
			gradv = vort_grad(i, j).m_y;
			gradlen = VFXEpoch::Vector2Df(gradu, gradv).length();
			gradu *= 1.0f / (gradlen + 10e-20f);
			gradv *= 1.0f / (gradlen + 10e-20f);
			vort_dir.m_x = gradu;
			vort_dir.m_y = gradv;
			vort_dir.m_z = 0.0f;
			vorticity.m_x = vorticity.m_y = 0.0f;
			vorticity.m_z = vort(i, j);
			fconf = VFXEpoch::Vector3Df::cross(vorticity, vort_dir);
			fconf.m_x *= (vort_conf_eps * (1.0f / (velocity.getDimY() - 2)));
			fconf.m_y *= (vort_conf_eps * (1.0f / (velocity.getDimY() - 2)));
			fconf.m_z *= (vort_conf_eps * (1.0f / (velocity.getDimY() - 2)));

			// Check if cross product is correct.
			if (fconf.m_z > 0 || fconf.m_z < 0){
				cout << endl << "Vorticity confinement calculation error" << endl;
				continue;
			}

			fvort(i, j).m_x = fconf.m_x;
			fvort(i, j).m_y = fconf.m_y;
		}
	}

	for (int i = 1; i != velocity.getDimY() - 2; i++){
		for (int j = 1; j != velocity.getDimX() - 2; j++){
			velocity(i, j).m_x += (fvort(i, j).m_x + fvort(i - 1, j).m_x) * 0.5f;
			velocity(i, j).m_y += (fvort(i, j).m_y + fvort(i, j - 1).m_y) * 0.5f;
		}
	}
}

VFXEpoch::Vector2Df
SL2D::_get_vel(int dimX, int dimY, const VFXEpoch::Vector2Df& position, VFXEpoch::Grid2DfScalarField& u, VFXEpoch::Grid2DfScalarField& v)
{
	float u0, v0;
	u0 = VFXEpoch::InterpolateGrid(dimX, position.m_x * this->dimX - 0.0f, position.m_y * this->dimX - 0.5f, u);
	v0 = VFXEpoch::InterpolateGrid(dimY, position.m_x * this->dimY - 0.5f, position.m_y * this->dimY - 0.0f, v);
	return VFXEpoch::Vector2Df(u0, v0);
}

void SL2D::_clear()
{
	dimX = dimY = 0;
	iterations = 0;
	timeStep = 0.0;
	viscosity = 0.0f;
	diffuseRate = 0.0f;
	sourceRate = 0.0f;
	spacingX = spacingY = 0.0f;
	this->velocityFieldPrev.clear();
	this->velocityField.clear();
	this->densityFieldPrev.clear();
	this->densityField.clear();
	this->temperature.clear();
	this->temperaturePrev.clear();
	this->pressure.clear();
	this->gradPressure.clear();
	this->divergence.clear();
}
