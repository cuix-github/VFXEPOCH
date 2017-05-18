/*******************************************************************************
    VFXEpoch - Physically based simulation VFX

    Copyright (c) 2016 Snow Tsui <trevor.miscellaneous@gmail.com>

    All rights reserved. Use of this source code is governed by
    the MIT license as written in the LICENSE file.
*******************************************************************************/
#ifndef _SIM_BASE_H_
#define _SIM_BASE_H_
#include <vector>
#include "utl/UTL_Grid.h"
#include "utl/UTL_Matrix.h"
#include "utl/UTL_Vector.h"
#include "utl/UTL_General.h"
#include "utl/UTL_LinearSolvers.h"

using namespace std;

namespace VFXEpoch{
  namespace Solvers{
    class Euler_Fluid2D_Base{
    public:
      virtual bool init(){}
      virtual void step(double dt){}
      virtual void advect(){}
      virtual void presure_solve(double dt){}
      virtual void close(){}
      virtual void add_source(){}
    };

    class Euler_Fluid3D_Base{
    public:

      virtual bool init(){}
      virtual void step(double dt){}
      virtual void advect(){}
      virtual void presure_solve(double dt){}
      virtual void close(){}
      virtual void add_source(){}
    };
  }
}

#endif
