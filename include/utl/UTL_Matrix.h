/*******************************************************************************
    VFXEPOCH - Physically based simulation VFX

    Copyright (c) 2016 Snow Tsui <trevor.miscellaneous@gmail.com>

    All rights reserved. Use of this source code is governed by
    the MIT license as written in the LICENSE file.
*******************************************************************************/
#ifndef _UTL_MATRIX_H_
#define _UTL_MATRIX_H_

#include "UTL_Vector.h"

#define FOR_2x2_MATRIX for(i=0; i!=2; i++){for(j=0; j!=2; j++){
#define FOR_3x3_MATRIX for(i=0; i!=3; i++){for(j=0; j!=3; j++){
#define FOR_NxN_MATRIX for(i=0; i!=N; i++){for(j=0; j!=N; j++){
#define END_FOR }}

namespace VFXEpoch
{
	template <class T>
	class Matrix3x3
	{
	public:
		T	m00, m01, m02,
			m10, m11, m12,
			m20, m21, m22;
	public:
		Matrix3x3(){}
		Matrix3x3(
			T _v00, T _v01, T _v02,
			T _v10, T _v11, T _v12,
			T _v20, T _v21, T _v22) :
			m00(_v00), m01(_v01), m02(_v02),
			m10(_v10), m11(_v11), m12(_v12),
			m20(_v20), m21(_v21), m22(_v22){}

		Matrix3x3(const Matrix3x3& source){
			m00 = source.m00; m01 = source.m01; m02 = source.m02;
			m10 = source.m10; m11 = source.m11; m12 = source.m12;
			m20 = source.m20; m21 = source.m21; m22 = source.m22;
		}

		~Matrix3x3(){}

		Matrix3x3& operator=(const Matrix3x3& source)	{
			m00 = source.m00; m01 = source.m01; m02 = source.m02;
			m10 = source.m10; m11 = source.m11; m12 = source.m12;
			m20 = source.m20; m21 = source.m21; m22 = source.m22;
			return *this;
		}

		Matrix3x3 operator+(Matrix3x3& source)	{
			return Matrix3x3(m00 + source.m00, m01 + source.m01, m02 + source.m02,
							 m10 + source.m10, m11 + source.m11, m12 + source.m12,
							 m20 + source.m20, m21 + source.m21, m22 + source.m22);
		}

		Matrix3x3 operator-(Matrix3x3& source)	{
			return Matrix3x3(m00 - source.m00, m01 - source.m01, m02 - source.m02,
							 m10 - source.m10, m11 - source.m11, m12 - source.m12,
							 m20 - source.m20, m21 - source.m21, m22 - source.m22);
		}

		Matrix3x3& operator+=(const Matrix3x3& source)	{
			m00 += source.m00; m01 += source.m01; m02 += source.m02;
			m10 += source.m10; m11 += source.m11; m12 += source.m12;
			m20 += source.m20; m21 += source.m20; m22 += source.m22;
			return *this;
		}

		Matrix3x3& operator-=(const Matrix3x3& source)	{
			m00 -= source.m00; m01 -= source.m01; m02 -= source.m02;
			m10 -= source.m10; m11 -= source.m11; m12 -= source.m12;
			m20 -= source.m20; m21 -= source.m20; m22 -= source.m22;
			return *this;
		}

		Matrix3x3 operator*(Matrix3x3& source)	{
			return Matrix3x3(m00 * source.m00 + m01 * source.m10 + m02 * source.m20, m00 * source.m01 + m01 * source.m11 + m02 * source.m21, m00 * source.m02 + m01 * source.m12 + m02 * source.m22,
							 m10 * source.m00 + m11 * source.m10 + m12 * source.m20, m10 * source.m01 + m11 * source.m11 + m12 * source.m21, m10 * source.m02 + m11 * source.m12 + m12 * source.m22,
							 m20 * source.m00 + m21 * source.m10 + m22 * source.m20, m20 * source.m01 + m21 * source.m11 + m22 * source.m21, m20 * source.m02 + m21 * source.m12 + m22 * source.m22);
		}

		Matrix3x3 operator*(T scale)	{
			return Matrix3x3(m00 * scale, m01 * scale, m02 * scale,
							 m10 * scale, m11 * scale, m12 * scale,
							 m20 * scale, m21 * scale, m22 * scale);
		}

		Matrix3x3 operator/(T scale)	{
			return Matrix3x3(m00 / scale, m01 / scale, m02 / scale,
							 m10 / scale, m11 / scale, m12 / scale,
							 m20 / scale, m21 / scale, m22 / scale);
		}

		T Det()	{
			return m00 * (m11 * m22 - m12 * m21) - m01 * (m10 * m22 - m12 * m20) + m02 * (m10 * m21 - m11 * m20);
		}

		void Transpose()	{
			T temp;
			temp = m01; m01 = m10; m10 = temp;
			temp = m02;	m02 = m20; m20 = temp;
			temp = m12;	m12 = m21; m21 = temp;
		}

		void Inverse()	{
			T	i00, i01, i02,
				i10, i11, i12,
				i20, i21, i22;

			T det = this->Det();

			i00 = m11 * m22 - m12 * m21; i01 = m02 * m21 - m01 * m22; i02 = m01 * m12 - m02 * m11;
			i10 = m12 * m20 - m10 * m22; i11 = m00 * m22 - m02 * m20; i12 = m02 * m10 - m00 * m12;
			i20 = m10 * m21 - m11 * m20; i21 = m01 * m20 - m00 * m21; i22 = m00 * m11 - m01 * m10;

			m00 = i00 / det; m01 = i01 / det; m02 = i02 / det;
			m10 = i10 / det; m11 = i11 / det; m12 = i12 / det;
			m20 = i20 / det; m21 = i21 / det; m22 = i22 / det;
		}

		void RoundToEpsilon()	{
			ROUND(m00); ROUND(m01); ROUND(m02);
			ROUND(m10); ROUND(m11); ROUND(m12);
			ROUND(m20); ROUND(m21); ROUND(m22);
		}
	};

	typedef Matrix3x3<float> Matrix3x3f;
	typedef Matrix3x3<double> Matrix3x3d;

	template <class T>
	class MatrixNxN
	{
	private:
		std::vector<T> m_data;
		int size;

	public:
		MatrixNxN(){ size = 0; m_data.clear(); }
		MatrixNxN(int N) : size(N){ m_data.resize(size * size); }
		MatrixNxN(const MatrixNxN& source){ size = source.size; m_data = source.m_data; }
		~MatrixNxN(){ m_data.clear(); }

		MatrixNxN& operator=(const MatrixNxN& source) {
			size = source.size;
			m_data = source.m_data;
			return *this;
		}

		void setData(T data, int i, int j) {
			if (i > (size - 1) || j > (size - 1) || i < 0 || j < 0)
				assert(i <= (size - 1) && j <= (size - 1));
			else
				m_data[i * size + j] = data;
		}

		T getData(int i, int j)	{
			if (i >(size - 1) || j >(size - 1) || i < 0 || j < 0)
				assert(i <= (size - 1) && j <= (size - 1));
			else
				return m_data[i * size + j];
		}

		MatrixNxN operator+(MatrixNxN& source) {
			assert(size == source.size);

			MatrixNxN result(size);
			for (int i = 0; i != size; i++)
			{
				for (int j = 0; j != size; j++)
				{
					result.setData(m_data[i * size + j] + source.m_data[i * size + j], i, j);
				}
			}
			return result;
		}

		MatrixNxN operator-(MatrixNxN& source) {
			assert(size == source.size);

			MatrixNxN result(size);
			for (int i = 0; i != size; i++)
			{
				for (int j = 0; j != size; j++)
				{
					result.setData(m_data[i * size + j] - source.m_data[i * size + j], i, j);
				}
			}

			return result;
		}

		void zeros() {
			for (int i = 0; i != size; i++)
			{
				for (int j = 0; j != size; j++)
				{
					m_data[i * size + j] = 0;
				}
			}
		}

		MatrixNxN getDiagonal()
		{
			MatrixNxN result(size);
			result.zeros();

			for (int i = 0; i != size; i++)
			{
				for (int j = 0; j != size; j++)
				{
					if (i == j)
						result.m_data[i * size + j] = m_data[i * size + j];
				}
			}
			return result;
		}


	};

	typedef MatrixNxN<float> MatrixNxNf;
	typedef MatrixNxN<double> MatrixNxNd;
};

#endif
