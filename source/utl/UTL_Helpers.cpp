/*******************************************************************************
    VFXEPOCH - Physically based simulation VFX

    Copyright (c) 2016 Snow Tsui <trevor.miscellaneous@gmail.com>

    All rights reserved. Use of this source code is governed by
    the MIT license as written in the LICENSE file.
*******************************************************************************/
#include "UTL_Helpers.h"

using namespace Helpers;

void
Helpers::displayFieldf(int row, int col, float* field){
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
Helpers::displayScalarField(VFXEpoch::Grid2DfScalarField field){
	int row = field.getDimY();
	int col = field.getDimX();
	for (int i = 0; i != row + 1; i++){
		for (int j = 0; j != col + 1; j++){
			float m = (float)i;
			float n = (float)j;
			if (0 == i && 0 == j) continue;
			else if(0 == i && 1 < j) {
				cout << std::setprecision(0) << setiosflags(ios::fixed);
				cout << std::setw(7) << n - 1.f << ": ";
				continue;
			}
			else if(0 == i && 1 == j) {
				cout << std::setw(7) << "0: ";
				continue;
			}
			else if(0 == j) {
				cout << std::setprecision(0) << setiosflags(ios::fixed);
				cout << m - 1.f << ": "; continue;
			}
			cout << std::setprecision(4) << setiosflags(ios::fixed);
			cout << std::setw(7) << field(i - 1, j - 1) << ", ";
		}
		cout << endl;
	}
}

void
Helpers::displayVectorField(VFXEpoch::Grid2DVector2DfField field){
	cout << std::setprecision(4) << setiosflags(ios::fixed);
	int row = field.getDimY();
	int col = field.getDimX();

	/*********** DO NOT CHANGE THE RATIO space & index_space ***********/
	int space = 2;
	int index_space = 17;
	/*********** DO NOT CHANGE THE RATIO space & index_space ***********/

	for (int i = 0; i != row + 1; i++){
		for (int j = 0; j != col + 1; j++){
			float m = (float)i;
			float n = (float)j;
			if (0 == i && 0 == j) continue;
			else if(0 == i && 1 < j) {
				cout << std::setprecision(0) << setiosflags(ios::fixed);
				cout << std::setw(index_space) << n - 1.f << ": ";
				continue;
			}
			else if(0 == i && 1 == j) {
				cout << std::setw(6) << "0: ";
				continue;
			}
			else if(0 == j) {
				cout << std::setprecision(0) << setiosflags(ios::fixed);
				cout << m - 1.f << ": "; continue;
			}
			cout << std::setprecision(4) << setiosflags(ios::fixed);
			cout << "v(" << field(i - 1, j - 1).m_x << "," <<	field(i - 1, j - 1).m_y << ") " << std::setw(space);
			0.f > field(i - 1, j - 1).m_x && 0.f > field(i - 1, j - 1).m_y ? cout << std::setw(space) : cout << std::setw(space + 2);
			if((0.f > field(i - 1, j - 1).m_x && 0.f < field(i - 1, j - 1).m_y) ||
			(0.f < field(i - 1, j - 1).m_x && 0.f > field(i - 1, j - 1).m_y)) cout << std::setw(space + 1);
		}
		cout << endl;
		cout << std::setw(0);
	}
}

void
Helpers::displayCellStatus(VFXEpoch::Grid2DCellTypes field){
	int row = field.getDimY();
	int col = field.getDimX();
	for (int i = 0; i != row + 1; i++){
		for (int j = 0; j != col + 1; j++){
			float m = (float)i;
			float n = (float)j;
			if (0 == i && 0 == j) continue;
			else if(0 == i && 1 < j) {
				cout << std::setprecision(0) << setiosflags(ios::fixed);
				cout << std::setw(3) << n - 1.f << ": ";
				continue;
			}
			else if(0 == i && 1 == j) {
				cout << std::setw(8) << "0: ";
				continue;
			}
			else if(0 == j) {
				cout << std::setprecision(0) << setiosflags(ios::fixed);
				cout << m - 1.f << ": "; continue;
			}
			cout << std::setprecision(4) << setiosflags(ios::fixed);
			if(VFXEpoch::BOUNDARY_MASK::SOMETHING == field(i - 1, j - 1))
			cout << std::setw(3) << "B" << ", ";
			else cout << std::setw(3) << "F" << ", ";
		}
		cout << endl;
	}
}

void
Helpers::randomInitScalarField(VFXEpoch::Grid2DfScalarField& field, float min, float max)
{
	for (int i = 1; i != field.getDimY() - 1; i++)
	{
		for (int j = 1; j != field.getDimX() - 1; j++)
		{
			// float data = VFXEpoch::RandomI(min, max);
			// field.setData(data, i, j);

			// TODO: Fix the random number generator which causes core dump
			// field(i, j) = 0.0f;

			field(i, j) = VFXEpoch::RandomF(min, max);
		}
	}
}

void
Helpers::randomInitVectorField(VFXEpoch::Grid2DVector2DfField& field, float min, float max)
{
	for (int i = 1; i != field.getDimY() - 1; i++)
	{
		for (int j = 1; j != field.getDimX() - 1; j++)
		{
			// Vector2Df data(0.0f, 0.0f);
			// data.m_x = VFXEpoch::RandomI(min, max);
			// data.m_y = VFXEpoch::RandomI(min, max);
			// field.setData(data, i, j);

			// TODO: Fix the random number generator which causes core dump
			// field(i, j) = data;

			field(i, j) = VFXEpoch::Vector2Df(VFXEpoch::RandomF(min, max),
																				VFXEpoch::RandomF(min, max));
		}
	}
}

void
Helpers::randomInitCellStatus(VFXEpoch::Grid2DCellTypes &field){
	for (int i = 0; i != field.getDimY() - 1; i++){
		for(int j = 0; j != field.getDimX() - 1; j++){
			int flag = VFXEpoch::RandomI(0, 1);
			if(0 == flag) field(i, j) = VFXEpoch::BOUNDARY_MASK::NOTHING;
			else field(i, j) = VFXEpoch::BOUNDARY_MASK::SOMETHING;
		}
	}
}
