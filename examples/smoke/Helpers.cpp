/*******************************************************************************
    VFXEPOCH - Physically based simulation VFX

    Copyright (c) 2016 Snow Tsui <trevor.miscellaneous@gmail.com>

    All rights reserved. Use of this source code is governed by
    the MIT license as written in the LICENSE file.
*******************************************************************************/
#include "Helpers.h"

using namespace Helper;

void
Helper::displayFieldf(int row, int col, float* field){
	cout << std::setprecision(4) << setiosflags(ios::fixed);
	for (int i = 0; i != row; i++){
		for (int j = 0; j != col; j++){
			cout << std::setw(3) << field[i * col + j] << ", ";
			if (j == col - 1)
				cout << endl;
		}
	}
}

void
Helper::displayScalarField(int row, int col, VFXEpoch::Grid2DfScalarField field){
	cout << std::setprecision(4) << setiosflags(ios::fixed);
	for (int i = 0; i != row; i++){
		for (int j = 0; j != col; j++){
			cout << std::setw(7) << field.getData(i, j) << ", ";
			if (j == col - 1)
				cout << endl;
		}
	}
}

void
Helper::displayVectorField(int row, int col, VFXEpoch::Grid2DVector2DfField field){
	cout << std::setprecision(4) << setiosflags(ios::fixed);
	for (int i = 0; i != row; i++){
		for (int j = 0; j != col; j++){
			cout << std::setw(3) << "v(" << field.getData(i, j).m_x << "," << field.getData(i, j).m_y << ") ";
			if (j == col - 1)
				cout << endl;
		}
	}
}

void
Helper::randomInitScalarField(VFXEpoch::Grid2DfScalarField& field, float min, float max)
{
	for (int i = 1; i != field.getDimY() - 1; i++)
	{
		for (int j = 1; j != field.getDimX() - 1; j++)
		{
			float data = VFXEpoch::RandomI(min, max);
			field.setData(data, i, j);
		}
	}
}

void
Helper::randomInitVectorField(VFXEpoch::Grid2DVector2DfField& field, float min, float max)
{
	for (int i = 1; i != field.getDimY() - 1; i++)
	{
		for (int j = 1; j != field.getDimX() - 1; j++)
		{
			Vector2Df data;
			data.m_x = VFXEpoch::RandomI(min, max);
			data.m_y = VFXEpoch::RandomI(min, max);
			field.setData(data, i, j);
		}
	}
}
