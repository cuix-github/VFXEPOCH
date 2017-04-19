/*******************************************************************************
    VFXEPOCH - Physically based simulation VFX

    Copyright (c) 2016 Snow Tsui <trevor.miscellaneous@gmail.com>

    All rights reserved. Use of this source code is governed by
    the MIT license as written in the LICENSE file.
*******************************************************************************/
#include "Helpers.h"

const int N = 4;
const float source = 1.0f;
const float stride = 1.0 / (N + 1);

using namespace Helper;

int main(int argc, char** argv)
{
	VFXEpoch::LBM2D lbm_solver(N + 2, N + 2);

	system("Pause");
	return 0;
}
