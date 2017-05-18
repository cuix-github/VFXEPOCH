/*******************************************************************************
    VFXEPOCH - Physically based simulation VFX

    Copyright (c) 2016 Snow Tsui <trevor.miscellaneous@gmail.com>

    All rights reserved. Use of this source code is governed by
    the MIT license as written in the LICENSE file.
*******************************************************************************/

/*******************************************************************************
* Desc:
* Be aware of the boundary in the grid, so the actual length of
* size N grid should be N+2.(Imagine that you need one more cell
* on both sides.
* e.g:
* If we are building 3x3 grid for the simulation, actually we need
* 5x5 grid as the outter sides are boundaries.
*
*    N+2
* /       \
* # # # # # \
* #	* *	* #
* # * * * #  N+2
* # * * * #
* # # # # # /
* N = 3;
*******************************************************************************/
#ifndef _UTL_GRID_H_
#define _UTL_GRID_H_

#include "UTL_Matrix.h"

#define LOOP_GRID2D(grid)	for(int i=0; i != grid.getDimY(); i++) \
								for (int j=0; j != grid.getDimX(); j++)
#define LOOP_GRID3D(grid)	for(int i=0; i != grid.getDimZ(); i++) \
								for(int j=0; j != grid.getDimY(); j++) \
									for(int k=0; k != grid.getDimX(); k++)

#define IDX2D(i, j) ((i) * (m_xCell) + (j))
#define IDX3D(i, j, k) ((i) * (m_xCell * m_yCell) + (j) * (m_xCell) + (k))

namespace VFXEpoch
{
	enum class BOUNDARY
	{
		DIRICHLET,
		NEUMANN_CLOSE,
		NEUMANN_OPEN,
		STREAK,
		ROBIN,
		CAUCHY,
	};

	enum class BOUNDARY_MASK{
		NOTHING,
		SOMETHING
	};

	enum class EDGES_2DSIM
	{
		TOP, BOTTOM,
		LEFT, RIGHT
	};

	enum class FACES_3DSIM
	{
		FRONT, REAR,
		LEFT, RIGHT,
		TOP, BOTTOM
	};

	enum class FIELD_TYPE
	{
		SCALAR_FIELD,
		VECTOR_FIELD
	};

	typedef struct _boundary_state_2D
	{
		EDGES_2DSIM side;
		BOUNDARY boundaryType;
	}BoundaryState2D;

	typedef struct _boundary_state_3D
	{
		FACES_3DSIM face;
		BOUNDARY boundaryType;
	}BoundaryState3D;

	typedef struct _boundary_condition_per_edge
	{
		VFXEpoch::EDGES_2DSIM side;
		VFXEpoch::BOUNDARY boundaryType;
	}BndConditionPerEdge;

	enum class COMPUTATIONAL_SCALAR_FIELD_2D
	{
		DENSITY,
		DENSITY_PREV,
		TEMPERATURE,
		TEMPERATURE_PREV,
		PRESSURE,
		DIVERGENCE,
	};

	enum class COMPUTATIONAL_VECTOR_FIELD_2D
	{
		VEL_U,
		VEL_V,
		VEL_U_PREV,
		VEL_V_PREV,
		VEL,
		VEL_PREV,
		GRAVITY,
		GRAD_PRESSURE,
	};

	enum class COMPUTATIONAL_SCALAR_FIELD_3D
	{
		//TODO: Potential 3D scalar field
	};

	enum class COMPUTATIONAL_VECTOR_FIELD_3D
	{
		//TODO: Potential 3D vector field
	};

	enum class COMPUTATIONAL_FIELD_2D
	{
		COMPUTATIONAL_SCALAR_FIELD_2D,
		COMPUTATIONAL_VECTOR_FIELD_2D,
	};

	enum class COMPUTATIONAL_FIELD_3D
	{
		COMPUTATIONAL_SCALAR_FIELD_3D,
		COMPUTATIONAL_VECTOR_FIELD_3D,
	};

	// xCell for column loop
	// yCell for row loop
	template <class T>
	class Grid2D
	{
	public:
		int m_xCell, m_yCell;
		float dx, dy;
		std::vector<T> data;
		BoundaryState2D boundaryState[4];

	public:
		Grid2D(){ m_xCell = m_yCell = 0; dx = 0.0f; data.clear(); }
		Grid2D(int x, int y) : m_xCell(x), m_yCell(y){ data.clear(); data.resize(m_xCell * m_yCell); }
		Grid2D(int x, int y, float _dx, float _dy) : m_xCell(x), m_yCell(y), dx(_dx), dy(_dy){ data.clear(); data.resize(m_xCell * m_yCell); }
		Grid2D(const Grid2D& source){ this->m_xCell = source.m_xCell; this->m_yCell = source.m_yCell; this->dx = source.dx; this->dy = source.dy; this->data.clear(); this->data = source.data; }
		Grid2D<T>& operator=(const Grid2D<T>& source) {
			m_xCell = source.m_xCell;
			m_yCell = source.m_yCell;
			dx = source.dx;
			dy = source.dy;

			data.clear();
			data = source.data;
			return *this;
		}

		friend Grid2D<T> operator+(Grid2D<T>& a, Grid2D<T>& b) {
			if (a.m_yCell != b.m_yCell || a.m_xCell != b.m_xCell)
				assert(a.m_yCell == b.m_yCell && a.m_xCell == b.m_xCell);

			Grid2D<T> result(a);
			for (int i = 0; i != result.getDimY(); i++){
				for (int j = 0; j != result.getDimX(); j++){
					T t1 = a.getData(i, j);
					T t2 = b.getData(i, j);
					result.setData(t1 + t2, i, j);
				}
			}

			return result;
		}

		friend Grid2D<T> operator-(Grid2D<T>& a, Grid2D<T>& b) {
			if (a.m_yCell != b.m_yCell || a.m_xCell != b.m_xCell)
				assert(a.m_yCell == b.m_yCell && a.m_xCell == b.m_xCell);

			Grid2D<T> result(a);
			for (int i = 0; i != result.getDimY(); i++){
				for (int j = 0; j != result.getDimX(); j++){
					T t1 = a.getData(i, j);
					T t2 = b.getData(i, j);
					result.setData(t1 - t2, i, j);
				}
			}

			return result;
		}

		Grid2D<T>& operator+=(const VFXEpoch::Grid2D<T>& rhs) {
			assert(rhs.m_xCell == m_xCell && rhs.m_yCell == m_yCell);

			for (int i = 0; i != m_yCell; i++){
				for (int j = 0; j != m_xCell; j++){
					data[IDX2D(i, j)] += rhs(i, j);
				}
			}

			return *this;
		}

		Grid2D<T>& operator-=(const VFXEpoch::Grid2D<T>& rhs) {
			assert(rhs.m_xCell == m_xCell && rhs.m_yCell == m_yCell);

			for (int i = 0; i != m_yCell; i++){
				for (int j = 0; j != m_xCell; j++){
					data[IDX2D(i, j)] -= rhs(i, j);
				}
			}

			return *this;
		}

		// TODO: Check the logic for m_xCell & m_yCell
		const T& operator()(int i, int j) const {
			assert(i >= 0 && i <= m_yCell - 1 && j >= 0 && j <= m_xCell - 1);
			return data[IDX2D(i, j)];
		}

		// TODO: Check the logic for m_xCell & m_yCell
		T& operator()(int i, int j)	{
			assert(i >= 0 && i <= m_yCell - 1 && j >= 0 && j <= m_xCell - 1);
			return data[IDX2D(i, j)];
		}

		// TODO: Overload operator "*"

		~Grid2D(){ m_xCell = m_yCell = 0; dx = dy = 0.0f; data.clear(); }

	public:
		void zeroVectors(){
			int size = m_xCell * m_yCell;
			for (int i = 0; i != size; i++) {
				data[i].m_x = data[i].m_y = 0.0f;
			}
		}

		void zeroScalars(){
			int size = m_xCell * m_yCell;
			for (int i = 0; i != size; i++) {
				data[i] = 0.0f;
			}
		}

		void ResetDimension(int xCell, int yCell){
			data.clear();
			m_xCell = xCell;
			m_yCell = yCell;
			data.resize(xCell * yCell);
		}

		void setData(T _data, int i, int j) {
			assert(i >= 0 && i <= (m_yCell - 1) && j >= 0 && j <= (m_xCell - 1));
			data[IDX2D(i, j)] = _data;
		}

		T getData(int i, int j)	{
			assert(i >= 0 && i <= (m_yCell - 1) && j >= 0 && j <= (m_xCell - 1));
			return data[IDX2D(i, j)];
		}

		Grid2D<T> scale(T f){
			for (int i = 0; i != m_yCell; i++){
				for (int j = 0; j != m_xCell; j++){
					data[IDX2D(i, j)] *= f;
				}
			}
			return *this;
		}

		const BoundaryState2D* getBoundaryState() const{
			return boundaryState;
		}

		inline int getDimY(){
			return m_yCell;
		}

		inline int getDimX(){
			return m_xCell;
		}

		inline float getDy(){
			return dy;
		}

		inline float getDx(){
			return dx;
		}

		// TODO: Grid properties
		void setBoundaries(BOUNDARY boundaryType, EDGES_2DSIM edge)
		{
			// Boundary types:
			// -> Dirichlet
			// -> Neumann-close
			// -> Neumann-open
			// -> Streak - zero out boundaries
			// -> Cauchy
			// -> Robin
			switch (edge)
			{
			case VFXEpoch::EDGES_2DSIM::TOP:
				if (boundaryType == BOUNDARY::DIRICHLET)
				{
					// TODO: Dirichlet boundary condition processing
				}
				else if (boundaryType == BOUNDARY::NEUMANN_CLOSE) {
					for (int i = 1; i <= m_xCell - 2; i++) {
						data[IDX2D(0, i)] = -data[IDX2D(1, i)];
					}
					boundaryState[0].side = EDGES_2DSIM::TOP;
					boundaryState[0].boundaryType = BOUNDARY::NEUMANN_CLOSE;
				}
				else if (boundaryType == BOUNDARY::NEUMANN_OPEN) {
					for (int i = 1; i <= m_xCell - 2; i++){
						data[IDX2D(0, i)] = data[IDX2D(1, i)];
					}
					boundaryState[0].side = EDGES_2DSIM::TOP;
					boundaryState[0].boundaryType = BOUNDARY::NEUMANN_OPEN;
				}
				else if (boundaryType == BOUNDARY::STREAK) {
					for (int i = 1; i <= m_xCell - 2; i++){
						data[IDX2D(0, i)] = T();
					}
					boundaryState[0].side = EDGES_2DSIM::TOP;
					boundaryState[0].boundaryType = BOUNDARY::STREAK;
				}
				else if (boundaryType == BOUNDARY::CAUCHY) {
					// TODO: Cauchy boundary processing
				}
				else {
					// TODO: Robin boundary processing
				}
				break;
			case VFXEpoch::EDGES_2DSIM::BOTTOM:
				if (boundaryType == BOUNDARY::DIRICHLET)
				{
					// TODO: Dirichlet boundary condition processing
				}
				else if (boundaryType == BOUNDARY::NEUMANN_CLOSE) {
					for (int i = 1; i <= m_xCell - 2; i++) {
						data[IDX2D(m_yCell - 1, i)] = -data[IDX2D(m_yCell - 2, i)];
					}
					boundaryState[1].side = EDGES_2DSIM::BOTTOM;
					boundaryState[1].boundaryType = BOUNDARY::NEUMANN_CLOSE;
				}
				else if (boundaryType == BOUNDARY::NEUMANN_OPEN) {
					for (int i = 1; i <= m_xCell - 2; i++){
						data[IDX2D(m_yCell - 1, i)] = data[IDX2D(m_yCell - 2, i)];
					}
					boundaryState[1].side = EDGES_2DSIM::BOTTOM;
					boundaryState[1].boundaryType = BOUNDARY::NEUMANN_OPEN;
				}
				else if (boundaryType == BOUNDARY::STREAK) {
					for (int i = 1; i <= m_xCell - 2; i++){
						data[IDX2D(m_yCell - 1, i)] = T();
					}
					boundaryState[1].side = EDGES_2DSIM::BOTTOM;
					boundaryState[1].boundaryType = BOUNDARY::STREAK;
				}
				else if (boundaryType == BOUNDARY::CAUCHY) {
					// TODO: Cauchy boundary processing
				}
				else {
					// TODO: Robin boundary processing
				}
				break;
			case VFXEpoch::EDGES_2DSIM::LEFT:
				if (boundaryType == BOUNDARY::DIRICHLET)
				{
					// TODO: Dirichlet boundary condition processing
				}
				else if (boundaryType == BOUNDARY::NEUMANN_CLOSE) {
					for (int i = 1; i <= m_yCell - 2; i++){
						data[IDX2D(i, 0)] = -data[IDX2D(i, 1)];
					}
					boundaryState[2].side = EDGES_2DSIM::LEFT;
					boundaryState[2].boundaryType = BOUNDARY::NEUMANN_CLOSE;
				}
				else if (boundaryType == BOUNDARY::NEUMANN_OPEN) {
					for (int i = 1; i <= m_yCell - 2; i++){
						data[IDX2D(i, 0)] = data[IDX2D(i, 1)];
					}
					boundaryState[2].side = EDGES_2DSIM::LEFT;
					boundaryState[2].boundaryType = BOUNDARY::NEUMANN_OPEN;
				}
				else if (boundaryType == BOUNDARY::STREAK) {
					for (int i = 1; i <= m_xCell - 2; i++){
						data[IDX2D(i, 0)] = T();
					}
					boundaryState[2].side = EDGES_2DSIM::LEFT;
					boundaryState[2].boundaryType = BOUNDARY::STREAK;
				}
				else if (boundaryType == BOUNDARY::CAUCHY) {
					// TODO :Cauchy boundary processing
				}
				else {
					// TODO: Robin boundary processing
				}
				break;
			case VFXEpoch::EDGES_2DSIM::RIGHT:
				if (boundaryType == BOUNDARY::DIRICHLET)
				{
					// TODO: Dirichlet boundary condition processing
				}
				else if (boundaryType == BOUNDARY::NEUMANN_CLOSE) {
					for (int i = 1; i <= m_yCell - 2; i++){
						data[IDX2D(i, m_xCell - 1)] = -data[IDX2D(i, m_xCell - 2)];
					}
					boundaryState[3].side = EDGES_2DSIM::RIGHT;
					boundaryState[3].boundaryType = BOUNDARY::NEUMANN_OPEN;
				}
				else if (boundaryType == BOUNDARY::NEUMANN_OPEN) {
					for (int i = 1; i <= m_yCell - 2; i++){
						data[IDX2D(i, m_xCell - 1)] = data[IDX2D(i, m_xCell - 2)];
					}
					boundaryState[3].side = EDGES_2DSIM::RIGHT;
					boundaryState[3].boundaryType = BOUNDARY::NEUMANN_CLOSE;
				}
				else if (boundaryType == BOUNDARY::STREAK) {
					for (int i = 1; i <= m_xCell - 2; i++){
						data[IDX2D(i, m_xCell - 1)] = T();
					}
					boundaryState[3].side = EDGES_2DSIM::LEFT;
					boundaryState[3].boundaryType = BOUNDARY::STREAK;
				}
				else if (boundaryType == BOUNDARY::CAUCHY) {
					// TODO: Cauchy boundary processing
				}
				else {
					// TODO: Robin boundary processing
				}
				break;
			default:
				break;
			}
		}

		void setBoundariesOnCorners()
		{
			data[IDX2D(0, 0)] = (data[IDX2D(0, 1)] + data[IDX2D(1, 0)]) * 0.5f;
			data[IDX2D(0, m_xCell - 1)] = (data[IDX2D(0, m_xCell - 2)] + data[IDX2D(1, m_xCell - 1)]) * 0.5f;
			data[IDX2D(m_yCell - 1, 0)] = (data[IDX2D(m_yCell - 2, 0)] + data[IDX2D(m_yCell - 1, 1)]) * 0.5f;
			data[IDX2D(m_yCell - 1, m_xCell - 1)] = (data[IDX2D(m_yCell - 2, m_xCell - 1)] + data[IDX2D(m_yCell - 1, m_xCell - 2)]) * 0.5f;
		}

		void clear(){
			m_xCell = 0;
			m_yCell = 0;
			data.clear();
		}
	};

	typedef Grid2D<float> Grid2DfScalarField;
	typedef Grid2D<double> Grid2DdScalarField;
	typedef Grid2D<int> Grid2DiScalarField;
	typedef Grid2D<VFXEpoch::Vector2Df> Grid2DVector2DfField;
	typedef Grid2D<VFXEpoch::Vector2Dd> Grid2DVector2DdField;
	typedef Grid2D<VFXEpoch::Vector2Di> Grid2DVector2DiField;
	typedef Grid2D<VFXEpoch::BOUNDARY_MASK> Grid2DCellTypes;

	template <class T>
	class Grid3D
	{
	public:
		int m_xCell, m_yCell, m_zCell;
		float dx, dy, dz;
		std::vector<T> data;

	private:
		BoundaryState3D boundaryState[6];

	public:
		Grid3D(){ m_xCell = m_yCell = m_zCell = 0; dx = dy = dz = 0.0f; data.clear(); }
		Grid3D(int x, int y, int z){ m_xCell = x; m_yCell = y; m_zCell = z; }
		Grid3D(int x, int y, int z, float _dx, float _dy, float _dz) : m_xCell(x), m_yCell(y), m_zCell(z), dx(_dx), dy(_dy), dz(_dz){ data.clear(); data.resize(m_xCell * m_yCell * m_zCell); }
		Grid3D(const Grid3D& source){ m_xCell = source.m_xCell; m_yCell = source.m_yCell; m_zCell = source.m_zCell; dx = source.dx; dy = source.dx; data = source.data; }
		Grid3D& operator=(const Grid3D& source)
		{
			m_xCell = source.m_xCell;
			m_yCell = source.m_yCell;
			m_zCell = source.m_zCell;
			dx = source.dx;
			dy = source.dy;
			dz = source.dz;

			data = source.data;
			return *this;
		}
		~Grid3D(){ clear(); }

	public:
		void zeroVectors(){
			int size = m_xCell * m_yCell * m_zCell;
			for (int i = 0; i != size; i++) {
				data[i].m_x = data[i].m_y = data[i].m_z = 0.0f;
			}
		}

		void zeroScalars(){
			int size = m_xCell * m_yCell * m_zCell;
			for (int i = 0; i != size; i++) {
				data[i] = 0.0f;
			}
		}

		void ResetDimension(int yCell, int xCell, int zCell) {
			data.clear();
			m_xCell = xCell;
			m_yCell = yCell;
			m_zCell = zCell;
			data.resize(yCell * xCell * zCell);
		}

		void setData(T _data, int i, int j, int k) {
			if (i > (m_yCell - 1) || j > (m_xCell - 1) || k > (m_zCell - 1) || i < 0 || j < 0 || k < 0)
				assert(i <= (m_yCell - 1) && j <= (m_xCell - 1) && k <= (m_zCell - 1));
			else
				data[IDX3D(i, j, k)] = _data;
		}

		T getData(int i, int j, int k)	{
			if (i >(m_yCell - 1) || j >(m_xCell - 1) || k >(m_zCell - 1) || i < 0 || j < 0 || k < 0)
				assert(i <= (m_yCell - 1) && j <= (m_xCell - 1) && k <= (m_zCell - 1));
			else
				return data[IDX3D(i, j, k)];
		}

		Grid3D<T> scale(T f){
			for (int k = 0; k != m_zCell; k++){
				for (int i = 0; i != m_yCell; i++){
					for (int j = 0; j != m_xCell; j++){
						data[k * m_zCell + i * m_yCell + j] *= f;
					}
				}
			}
			return *this;
		}

		const BoundaryState3D* getBoundaryState() const{
			return boundaryState;
		}

		inline int getDimY(){
			return m_yCell;
		}

		inline int getDimX(){
			return m_xCell;
		}

		inline int getDimZ(){
			return m_zCell;
		}

		inline float getDx(){
			return dx;
		}

		inline float getDy(){
			return dy;
		}

		inline float getDz(){
			return dz;
		}

		// TODO: Grid properties
		void setBoundaries(BOUNDARY boundaryType, FACES_3DSIM face)
		{
			switch (face)
			{
			case VFXEpoch::FACES_3DSIM::FRONT:
				if (boundaryType == VFXEpoch::BOUNDARY::DIRICHLET)
				{
					// TODO: Dirichlet boundary condition processing
				}
				else if (boundaryType == VFXEpoch::BOUNDARY::NEUMANN_CLOSE)
				{
					for (int i = 1; i != m_yCell - 1; i++){
						for (int j = 1; j != m_xCell - 1; j++){
							data[IDX3D(0, i, j)] = -data[IDX3D(1, i, j)];
						}
					}
					boundaryState[0].face = FACES_3DSIM::FRONT;
					boundaryState[0].boundaryType = BOUNDARY::NEUMANN_CLOSE;
				}
				else if (boundaryType == VFXEpoch::BOUNDARY::NEUMANN_OPEN)
				{
					for (int i = 1; i != m_yCell - 1; i++){
						for (int j = 1; j != m_xCell - 1; j++){
							data[IDX3D(0, i, j)] = data[IDX3D(1, i, j)];
						}
					}
					boundaryState[0].face = FACES_3DSIM::FRONT;
					boundaryState[0].boundaryType = BOUNDARY::NEUMANN_OPEN;
				}
				else if (boundaryType == VFXEpoch::BOUNDARY::CAUCHY)
				{
					// TODO: Cauchy boundary condition processing
				}
				else
				{
					// TODO: Robin boundary condition processing
				}
				break;
			case VFXEpoch::FACES_3DSIM::REAR:
				if (boundaryType == VFXEpoch::BOUNDARY::DIRICHLET)
				{
					// TODO: Dirichlet boundary condition processing
				}
				else if (boundaryType == VFXEpoch::BOUNDARY::NEUMANN_CLOSE)
				{
					for (int i = 1; i != m_yCell - 1; i++){
						for (int j = 1; j != m_xCell - 1; j++){
							data[IDX3D(m_zCell - 1, i, j)] = -data[IDX3D(m_zCell - 2, i, j)];
						}
					}
					boundaryState[1].face = FACES_3DSIM::REAR;
					boundaryState[1].boundaryType = BOUNDARY::NEUMANN_CLOSE;
				}
				else if (boundaryType == VFXEpoch::BOUNDARY::NEUMANN_OPEN)
				{
					for (int i = 1; i != m_yCell - 1; i++){
						for (int j = 1; j != m_xCell - 1; j++){
							data[IDX3D(m_zCell - 1, i, j)] = data[IDX3D(m_zCell - 2, i, j)];
						}
					}
					boundaryState[1].face = FACES_3DSIM::REAR;
					boundaryState[1].boundaryType = BOUNDARY::NEUMANN_OPEN;
				}
				else if (boundaryType == VFXEpoch::BOUNDARY::CAUCHY)
				{
					// TODO: Cauchy boundary condition processing
				}
				else
				{
					// TODO: Robin boundary condition processing
				}
				break;
			case VFXEpoch::FACES_3DSIM::LEFT:
				if (boundaryType == VFXEpoch::BOUNDARY::DIRICHLET)
				{
					// TODO: Dirichlet boundary condition processing
				}
				else if (boundaryType == VFXEpoch::BOUNDARY::NEUMANN_CLOSE)
				{
					for (int i = 1; i != m_yCell - 1; i++){
						for (int j = 1; j != m_zCell - 1; j++){
							data[IDX3D(j, i, 0)] = -data[IDX3D(j, i, 1)];
						}
					}
					boundaryState[2].face = FACES_3DSIM::LEFT;
					boundaryState[2].boundaryType = BOUNDARY::NEUMANN_CLOSE;
				}
				else if (boundaryType == VFXEpoch::BOUNDARY::NEUMANN_OPEN)
				{
					for (int i = 1; i != m_yCell - 1; i++){
						for (int j = 1; j != m_zCell - 1; j++){
							data[IDX3D(j, i, 0)] = data[IDX3D(j, i, 1)];
						}
					}
					boundaryState[2].face = FACES_3DSIM::LEFT;
					boundaryState[2].boundaryType = BOUNDARY::NEUMANN_OPEN;
				}
				else if (boundaryType == VFXEpoch::BOUNDARY::CAUCHY)
				{
					// TODO: Cauchy boundary condition processing
				}
				else
				{
					// TODO: Robin boundary condition processing
				}
				break;
			case VFXEpoch::FACES_3DSIM::RIGHT:
				if (boundaryType == VFXEpoch::BOUNDARY::DIRICHLET)
				{
					// TODO: Dirichlet boundary condition processing
				}
				else if (boundaryType == VFXEpoch::BOUNDARY::NEUMANN_CLOSE)
				{
					for (int i = 1; i != m_yCell - 1; i++){
						for (int j = 1; j != m_zCell - 1; j++){
							data[IDX3D(j, i, m_xCell - 1)] = -data[IDX3D(j, i, m_xCell - 2)];
						}
					}
					boundaryState[3].face = FACES_3DSIM::RIGHT;
					boundaryState[3].boundaryType = BOUNDARY::NEUMANN_CLOSE;
				}
				else if (boundaryType == VFXEpoch::BOUNDARY::NEUMANN_OPEN)
				{
					for (int i = 1; i != m_yCell - 1; i++){
						for (int j = 1; j != m_zCell - 1; j++){
							data[IDX3D(j, i, m_xCell - 1)] = data[IDX3D(j, i, m_xCell - 2)];
						}
					}
					boundaryState[3].face = FACES_3DSIM::RIGHT;
					boundaryState[3].boundaryType = BOUNDARY::NEUMANN_OPEN;
				}
				else if (boundaryType == VFXEpoch::BOUNDARY::CAUCHY)
				{
					// TODO: Cauchy boundary condition processing
				}
				else
				{
					// TODO: Robin boundary condition processing
				}
				break;
			case VFXEpoch::FACES_3DSIM::TOP:
				if (boundaryType == VFXEpoch::BOUNDARY::DIRICHLET)
				{
					// TODO: Dirichlet boundary condition processing
				}
				else if (boundaryType == VFXEpoch::BOUNDARY::NEUMANN_CLOSE)
				{
					for (int i = 1; i != m_zCell - 1; i++){
						for (int j = 1; j != m_xCell - 1; j++){
							data[IDX3D(i, 0, j)] = -data[IDX3D(i, 1, j)];
						}
					}
					boundaryState[4].face = FACES_3DSIM::TOP;
					boundaryState[4].boundaryType = BOUNDARY::NEUMANN_CLOSE;
				}
				else if (boundaryType == VFXEpoch::BOUNDARY::NEUMANN_OPEN)
				{
					for (int i = 1; i != m_zCell - 1; i++){
						for (int j = 1; j != m_xCell - 1; j++){
							data[IDX3D(i, 0, j)] = data[IDX3D(i, 1, j)];
						}
					}
					boundaryState[4].face = FACES_3DSIM::TOP;
					boundaryState[4].boundaryType = BOUNDARY::NEUMANN_OPEN;
				}
				else if (boundaryType == VFXEpoch::BOUNDARY::CAUCHY)
				{
					// TODO: Cauchy boundary condition processing
				}
				else
				{
					// TODO: Robin boundary condition processing
				}
				break;
			case VFXEpoch::FACES_3DSIM::BOTTOM:
				if (boundaryType == VFXEpoch::BOUNDARY::DIRICHLET)
				{
					// TODO: Dirichlet boundary condition processing
				}
				else if (boundaryType == VFXEpoch::BOUNDARY::NEUMANN_CLOSE)
				{
					for (int i = 1; i != m_zCell - 1; i++){
						for (int j = 1; j != m_xCell - 1; j++){
							data[IDX3D(i, m_yCell - 1, j)] = -data[IDX3D(i, m_yCell - 2, j)];
						}
					}
					boundaryState[5].face = FACES_3DSIM::BOTTOM;
					boundaryState[5].boundaryType = BOUNDARY::NEUMANN_CLOSE;
				}
				else if (boundaryType == VFXEpoch::BOUNDARY::NEUMANN_OPEN)
				{
					for (int i = 1; i != m_zCell - 1; i++){
						for (int j = 1; j != m_xCell - 1; j++){
							data[IDX3D(i, m_yCell - 1, j)] = data[IDX3D(i, m_yCell - 2, j)];
						}
					}
					boundaryState[5].face = FACES_3DSIM::BOTTOM;
					boundaryState[5].boundaryType = BOUNDARY::NEUMANN_OPEN;
				}
				else if (boundaryType == VFXEpoch::BOUNDARY::CAUCHY)
				{
					// TODO: Cauchy boundary condition processing
				}
				else
				{
					// TODO: Robin boundary condition processing
				}
				break;
			default:
				break;
			}
		}

		void setBoundariesOnEdges()
		{
			for (int i = 1; i != m_xCell - 1; i++){
				data[IDX3D(0, 0, i)] = 0.5f * (data[IDX3D(1, 0, i)] + data[IDX3D(0, 1, i)]);
				data[IDX3D(0, m_yCell - 1, i)] = 0.5f * (data[IDX3D(1, m_yCell - 1, i)] + data[IDX3D(0, m_yCell - 2, i)]);
				data[IDX3D(m_zCell - 1, 0, i)] = 0.5f * (data[IDX3D(m_zCell - 2, 0, i)] + data[IDX3D(m_zCell - 1, 1, i)]);
				data[IDX3D(m_zCell - 1, m_yCell - 1, i)] = 0.5f * (data[IDX3D(m_zCell - 1, m_yCell - 2, i)] + data[IDX3D(m_zCell - 2, m_yCell - 1, i)]);
			}

			for (int i = 1; i != m_yCell - 1; i++){
				data[IDX3D(0, i, 0)] = 0.5f * (data[IDX3D(1, i, 0)] + data[IDX3D(0, i, 1)]);
				data[IDX3D(0, i, m_xCell - 1)] = 0.5f * (data[IDX3D(0, i, m_xCell - 2)] + data[IDX3D(1, i, m_xCell - 1)]);
				data[IDX3D(m_zCell - 1, i, 0)] = 0.5f * (data[IDX3D(m_zCell - 2, i, 0)] + data[IDX3D(m_zCell - 1, i, 1)]);
				data[IDX3D(m_zCell - 1, i, m_xCell - 1)] = 0.5f * (data[IDX3D(m_zCell - 2, i, m_xCell - 1)] + data[IDX3D(m_zCell - 1, i, m_xCell - 2)]);
			}

			for (int i = 1; i != m_zCell - 1; i++){
				data[IDX3D(i, 0, 0)] = 0.5f * (data[IDX3D(i, 0, 1)] + data[IDX3D(i, 1, 0)]);
				data[IDX3D(i, m_yCell - 1, 0)] = 0.5f * (data[IDX3D(i, m_yCell - 2, 0)] + data[IDX3D(i, m_yCell - 1, 1)]);
				data[IDX3D(i, 0, m_xCell - 1)] = 0.5f * (data[IDX3D(i, 0, m_xCell - 2)] + data[IDX3D(i, 1, m_xCell - 1)]);
				data[IDX3D(i, m_yCell - 1, m_xCell - 1)] = 0.5f * (data[IDX3D(i, m_yCell - 1, m_xCell - 2)] + data[IDX3D(i, m_yCell - 2, m_xCell - 1)]);
			}
		}

		void setBoundariesOnCorners()
		{
			data[IDX3D(0, 0, 0)] = (data[IDX3D(0, 0, 1)] + data[IDX3D(1, 0, 0)] + data[IDX3D(0, 1, 0)]) / 3.0f;
			data[IDX3D(0, 0, m_xCell - 1)] = (data[IDX3D(0, 0, m_xCell - 2)] + data[IDX3D(1, 0, m_xCell - 1)] + data[IDX3D(0, 1, m_xCell - 1)]) / 3.0f;
			data[IDX3D(0, m_yCell - 1, 0)] = (data[IDX3D(0, m_yCell - 2, 0)] + data[IDX3D(1, m_yCell - 1, 0)] + data[IDX3D(0, m_yCell - 1, 1)]) / 3.0f;
			data[IDX3D(0, m_yCell - 1, m_xCell - 1)] = (data[IDX3D(0, m_yCell - 1, m_xCell - 2)] + data[IDX3D(1, m_yCell - 1, m_xCell - 1)] + data[IDX3D(0, m_yCell - 2, m_xCell - 1)]) / 3.0f;
			data[IDX3D(m_zCell - 1, 0, 0)] = (data[IDX3D(m_zCell - 1, 0, 1)] + data[IDX3D(m_zCell - 1, 1, 0)] + data[IDX3D(m_zCell - 2, 0, 0)]) / 3.0f;
			data[IDX3D(m_zCell - 1, 0, m_xCell - 1)] = (data[IDX3D(m_zCell - 1, 0, m_xCell - 2)] + data[IDX3D(m_zCell - 1, 1, m_xCell - 1)] + data[IDX3D(m_zCell - 2, 0, m_xCell - 1)]) / 3.0f;
			data[IDX3D(m_zCell - 1, m_yCell - 1, 0)] = (data[IDX3D(m_zCell - 1, m_yCell - 2, 0)] + data[IDX3D(m_zCell - 1, m_yCell - 1, 1)] + data[IDX3D(m_zCell - 2, m_zCell - 1, 0)]) / 3.0f;
			data[IDX3D(m_zCell - 1, m_yCell - 1, m_xCell - 1)] = (data[IDX3D(m_zCell - 1, m_yCell - 1, m_xCell - 2)] + data[IDX3D(m_zCell - 2, m_yCell - 1, m_xCell - 1)] + data[IDX3D(m_zCell - 1, m_yCell - 2, m_xCell - 1)]) / 3.0f;
		}

		void clear(){
			m_xCell = 0;
			m_yCell = 0;
			m_zCell = 0;
			data.clear();
			data.resize(0);
		}
	};

	typedef Grid3D<VFXEpoch::Vector3Df> Grid3DVector3DfField;
	typedef Grid3D<VFXEpoch::Vector3Dd> Grid3DVector3DdField;
	typedef Grid3D<VFXEpoch::Vector3Di> Grid3DVector3DiField;
	typedef Grid3D<float> Grid3DfScalarField;
	typedef Grid3D<double> Grid3DdScalarField;
	typedef Grid3D<int> Grid3DiScalarField;
};

#endif
