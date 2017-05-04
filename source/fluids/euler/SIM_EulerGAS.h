/*******************************************************************************
    VFXEPOCH - Physically based simulation VFX

    Copyright (c) 2016 Snow Tsui <trevor.miscellaneous@gmail.com>

    All rights reserved. Use of this source code is governed by
    the MIT license as written in the LICENSE file.
*******************************************************************************/
#ifndef _SIM_EULER_GAS_H_
#define _SIM_EULER_GAS_H_
#include "utl/UTL_Grid.h"
#include "utl/UTL_Matrix.h"
#include "utl/UTL_Vector.h"
#include "utl/UTL_General.h"
#include "utl/UTL_LinearSolvers.h"
#include "utl/PCGSolver/blas_wrapper.h"
#include "utl/PCGSolver/pcg_solver.h"
#include "utl/PCGSolver/sparse_matrix.h"
#include "utl/PCGSolver/util.h"
#include "fluids/euler/SIM_Base.h"

using namespace VFXEpoch;
using namespace VFXEpoch::Solvers;

namespace VFXEpoch{
  namespace Solvers{

    class EulerGAS2D;
    class EulerGAS3D;

    class EulerGAS2D : public Euler_Fluid2D_Base{
    public:
      EulerGAS2D();
      EulerGAS2D(const EulerGAS2D& src);
      EulerGAS2D& operator=(const EulerGAS2D& rhs);
      ~EulerGAS2D();
    public:
      // TODO: Ovserride the base functions
      bool init() override;
      void step(double dt) override;
      void advect() override;
      void presure_solve(double dt) override;
      void close() override;
    private:
      VFXEpoch::Vector2Di dimension;
      VFXEpoch::Grid2DfScalarField u, u0;
      VFXEpoch::Grid2DfScalarField v, v0;
      VFXEpoch::Grid2DfScalarField d, d0;
      VFXEpoch::Grid2DfScalarField t, t0;
      VFXEpoch::Grid2DfScalarField pressure;
      VFXEpoch::Grid2DfScalarField omega, omega0;
      REAL buoyancy_alpha, buoyancy_beta;
      REAL vort_conf_eps;
      REAL dt;
      REAL tolerance;
      int max_iterations;
      std::vector<VFXEpoch::BOUNDARY_MASK> domain_mask;
    };
  }
}

#endif
