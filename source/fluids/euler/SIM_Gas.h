/*******************************************************************************
    VFXEPOCH - Physically based simulation VFX

    Copyright (c) 2016 Snow Tsui <trevor.miscellaneous@gmail.com>

    All rights reserved. Use of this source code is governed by
    the MIT license as written in the LICENSE file.
*******************************************************************************/
#ifndef _SIM_GAS_H_
#define _SIM_GAS_H_

#include "../../utl/UTL_Grid.h"
#include "../../utl/UTL_Matrix.h"
#include "../../utl/UTL_Vector.h"
#include "../../utl/UTL_General.h"
#include "../../utl/UTL_LinearSolvers.h"
#include "SIM_Mac.h"

namespace VFXEpoch
{
	namespace Solvers
	{
		class SL2D;
		class SL3D;

		/*
		* Reference:
		* Jos Stam
		* Stalbe Fluids
		* Proceedings of ACM SIGGRAPH 1999, 26th annual conference on
		* computer graphics & interactive techniques
		* Link: http://www.autodeskresearch.com/pdf/ns.pdf
		*
		*/
		class SL2D
		{
		private:
			int dimX;
			int dimY;
			int iterations;
			double timeStep;
			float viscosity;
			float diffuseRate;
			float sourceRate;
			float spacingX;
			float spacingY;
			VFXEpoch::BndConditionPerEdge fieldBoundaries[4];
			VFXEpoch::Grid2DVector2DfField velocityField;
			VFXEpoch::Grid2DVector2DfField velocityFieldPrev;
			VFXEpoch::Grid2DVector2DfField gradPressure;
			VFXEpoch::Grid2DfScalarField pressure;
			VFXEpoch::Grid2DfScalarField divergence;
			VFXEpoch::Grid2DfScalarField densityField;
			VFXEpoch::Grid2DfScalarField densityFieldPrev;
			VFXEpoch::Grid2DfScalarField temperature;
			VFXEpoch::Grid2DfScalarField temperaturePrev;
			std::vector<VFXEpoch::Particle2D> particles;

			// If using staggered scheme, use following variable to store fluid variables
			VFXEpoch::Mac2D staggeredGrid;

		public:
			SL2D();
			SL2D(int dimX, int dimY, int iterations, double timeStep, float diffuseRate, float viscosity, float sourceRate, float spacingX, float spacingY,
				VFXEpoch::Grid2DVector2DfField velField, VFXEpoch::Grid2DVector2DfField velFieldPrev,
				VFXEpoch::Grid2DfScalarField denField, VFXEpoch::Grid2DfScalarField denFieldPrev, VFXEpoch::Grid2DfScalarField pressure,
				VFXEpoch::Grid2DfScalarField divergence, VFXEpoch::Grid2DVector2DfField gradPressure, VFXEpoch::Grid2DfScalarField temperature, VFXEpoch::Grid2DfScalarField temperaturePrev);
			SL2D(const SL2D& source);
			SL2D& operator=(const SL2D& source);
			~SL2D();

		public:
			// TODO: Higher supervision of Stable Fluids scheme functions
			bool Initialize(VFXEpoch::Vector2Di dimension, VFXEpoch::Vector2Df spacing, int iterations, float timeStep, float diffuseRate, float viscosity, float sourceRate);
			void SetField(VFXEpoch::Grid2DfScalarField field, VFXEpoch::COMPUTATIONAL_SCALAR_FIELD_2D fieldName);
			void SetField(VFXEpoch::Grid2DVector2DfField field, VFXEpoch::COMPUTATIONAL_VECTOR_FIELD_2D fieldName);
			void SetField(VFXEpoch::Mac2D mac);
			void SetFieldBoundary(VFXEpoch::BOUNDARY boundaryType, VFXEpoch::EDGES_2DSIM edge);
			void AddSource(VFXEpoch::Grid2DfScalarField& target, VFXEpoch::Grid2DfScalarField source);
			void AddSource(VFXEpoch::Grid2DVector2DfField& target, VFXEpoch::Grid2DVector2DfField source);
			void Diffuse(VFXEpoch::Grid2DfScalarField& target, VFXEpoch::Grid2DfScalarField source);
			void Diffuse(VFXEpoch::Grid2DVector2DfField& target, VFXEpoch::Grid2DVector2DfField source);
			void Advect(VFXEpoch::Grid2DfScalarField& target, VFXEpoch::Grid2DfScalarField previous_timestep_field, VFXEpoch::Grid2DVector2DfField ref);
			void Advect(VFXEpoch::Grid2DVector2DfField& target, VFXEpoch::Grid2DVector2DfField previous_timestep_field, VFXEpoch::Grid2DVector2DfField ref);
			void PressureSolve(VFXEpoch::Grid2DVector2DfField& target, VFXEpoch::Grid2DfScalarField pressure, VFXEpoch::Grid2DfScalarField divergence);
			void AddVortConf(VFXEpoch::Grid2DVector2DfField& target, float vort_conf_eps, VORT_METHODS method);
			void Reset();
			void Shutdown();

		public:
			// TODO: Public access functions
			const VFXEpoch::Grid2DfScalarField& getDensityField() const;
			const VFXEpoch::Grid2DfScalarField& getDensityPrevField() const;
			const VFXEpoch::Grid2DVector2DfField& getVelocityField() const;
			const VFXEpoch::Grid2DVector2DfField& getVelocityPrevField() const;
			const VFXEpoch::BndConditionPerEdge* getFieldBoundaries() const;
			VFXEpoch::Grid2DfScalarField& getDensityField();
			VFXEpoch::Grid2DfScalarField& getDensityPrevField();
			VFXEpoch::Grid2DfScalarField& getTemperatureField();
			VFXEpoch::Grid2DfScalarField& getTemperaturePrevField();
			VFXEpoch::Grid2DVector2DfField& getVelocityField();
			VFXEpoch::Grid2DVector2DfField& getVelocityPrevField();
			VFXEpoch::BndConditionPerEdge* getFieldBoundaries();


			// For test, make the microsolvers as public members.
		public:
			// TODO: Micro-Solvers
			bool _init(int dimX, int dimY, int iterations, double timeStep, float diffuseRate, float viscosity, float sourceRate, float spacingX, float spacingY);
			void _set_field(VFXEpoch::Grid2DfScalarField field, VFXEpoch::COMPUTATIONAL_SCALAR_FIELD_2D fieldName);
			void _set_field(VFXEpoch::Grid2DVector2DfField field, VFXEpoch::COMPUTATIONAL_VECTOR_FIELD_2D fieldName);
			void _set_field_boundary(VFXEpoch::Grid2DfScalarField& field);
			void _set_field_boundary(VFXEpoch::Grid2DVector2DfField& field);
			void _set_field(VFXEpoch::Mac2D mac);
			void _set_source(VFXEpoch::Grid2DfScalarField& scalarField, VFXEpoch::Grid2DfScalarField source);
			void _set_source(VFXEpoch::Grid2DVector2DfField& vectorField, VFXEpoch::Grid2DVector2DfField source);

			// TODO: You may want to advect other properties, add your own
			// advect function here
			// Runge-Kutta 2nd order ODEs integrator.
			Vector2Df Runge_Kutta_II(int dimX, int dimY, VFXEpoch::Grid2DfScalarField& u, VFXEpoch::Grid2DfScalarField& v, const VFXEpoch::Vector2Df& position, float dt);
			// ADVECT
			void _advect(VFXEpoch::Grid2DVector2DfField& vectorField, VFXEpoch::Grid2DVector2DfField vectorFieldOirigin, VFXEpoch::Grid2DVector2DfField referenceField);
			void _advect(VFXEpoch::Grid2DfScalarField& scalarField, VFXEpoch::Grid2DfScalarField scalarFieldOrigin, VFXEpoch::Grid2DVector2DfField referenceField);
			void _advect_rk2();
			void _advect_particles_rk2(VFXEpoch::Grid2DfScalarField& u, VFXEpoch::Grid2DfScalarField& v, std::vector<VFXEpoch::Particle2D>& particleContainer);
			// END ADVECT

			void _diffuse(VFXEpoch::Grid2DfScalarField& densityField, VFXEpoch::Grid2DfScalarField densityFieldOrigin);
			void _diffuse(VFXEpoch::Grid2DVector2DfField& vectorField, VFXEpoch::Grid2DVector2DfField vectorFieldOrigin);
			void _project(VFXEpoch::Grid2DVector2DfField& field, VFXEpoch::Grid2DfScalarField pressure, VFXEpoch::Grid2DfScalarField divergence);
			void _update_vel(VFXEpoch::Grid2DVector2DfField& dest, Grid2DVector2DfField _gradPressure);

			void _get_uniform_curls(VFXEpoch::Grid2DfScalarField& w, VFXEpoch::Grid2DVector2DfField vel);
			void _get_buoyancy(VFXEpoch::Grid2DfScalarField density, VFXEpoch::Grid2DfScalarField temp, VFXEpoch::Grid2DVector2DfField& buoyancy, float alpha, float beta);
			void _vort_conf_update(VFXEpoch::Grid2DVector2DfField& target, float vort_conf_eps, VFXEpoch::VORT_METHODS method);

			// TODO: Used for Runge-Kutta 2nd order
			VFXEpoch::Vector2Df _get_vel(int dimX, int dimY, const VFXEpoch::Vector2Df& position, VFXEpoch::Grid2DfScalarField& u, VFXEpoch::Grid2DfScalarField& v);

			void _clear();
		};

		class SL3D
		{
		public:
			SL3D();
			SL3D(const SL3D& src);
			SL3D& operator=(const SL3D& rhs);
			~SL3D();

		public:
			//TODO: Public interfaces

		private:
			//TODO: Kernel functions
		};
	}
}

#endif
