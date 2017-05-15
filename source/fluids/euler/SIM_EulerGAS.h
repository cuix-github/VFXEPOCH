/*******************************************************************************
    VFXEPOCH - Physically based simulation VFX

    Copyright (c) 2016 Snow Tsui <trevor.miscellaneous@gmail.com>

    All rights reserved. Use of this source code is governed by
    the MIT license as written in the LICENSE file.
*******************************************************************************/
#ifndef _SIM_EULER_GAS_H_
#define _SIM_EULER_GAS_H_
#include "fluids/euler/SIM_Base.h"
#include "utl/PCGSolver/util.h"
#include "utl/PCGSolver/sparse_matrix.h"
#include "utl/PCGSolver/blas_wrapper.h"
#include "utl/PCGSolver/pcg_solver.h"

using namespace VFXEpoch;
using namespace VFXEpoch::Solvers;

namespace VFXEpoch{
  namespace Solvers{

    class EulerGAS2D;
    class EulerGAS3D;

    class EulerGAS2D : public Euler_Fluid2D_Base{
    /***************************** User Parameters *****************************/
    public:
      struct Parameters{
      public:
        Parameters(){
          dimension.m_x = 0; dimension.m_y = 0; dt = 0.0;
          dt = 0.0;
          buoyancy_alpha = buoyancy_beta = 0.0;
          tolerance = 0.0;
          max_iterations = 0;
          num_particles = 0;
        }
        Parameters(Vector2Di _dimension, double _dt, double _buoyancy_alpha, double _buoyancy_beta, double _tolerance, int _max_iterations): dimension(_dimension), dt(_dt),
        buoyancy_alpha(_buoyancy_alpha), buoyancy_beta(_buoyancy_beta), tolerance(_tolerance), max_iterations(_max_iterations){}
        Parameters(const Parameters& src){
          dimension = src.dimension;
          dt = src.dt;
          buoyancy_alpha = src.buoyancy_alpha; buoyancy_beta = src.buoyancy_beta;
          tolerance = src.tolerance;
          max_iterations = src.max_iterations;
          num_particles = src.num_particles;
        }
        Parameters& operator=(const Parameters& rhs){
          dimension = rhs.dimension;
          dt = rhs.dt;
          buoyancy_alpha = rhs.buoyancy_alpha; buoyancy_beta = rhs.buoyancy_beta;
          tolerance = rhs.tolerance;
          max_iterations = rhs.max_iterations;
          num_particles = rhs.num_particles;
          return *this;
        }
        ~Parameters(){ this->clear(); }
      public:
        inline void clear(){
          dimension.m_x = dimension.m_y = 0;
          dt = 0.0;
          buoyancy_alpha = buoyancy_beta = 0.0;
          tolerance = 0.0;
          max_iterations = 0;
          num_particles = 0;
        }

        friend inline ostream&
        operator<<(ostream& os, const Parameters& params) {
          os << std::setprecision(4) << setiosflags(ios::fixed);
          os << "Dimension = " << params.dimension.m_x << " x " << params.dimension.m_y << endl;
          os << "Time step = " << params.dt << endl;
          os << "Buoyancy alpha & beta = "
             << params.buoyancy_alpha << ", "
             << params.buoyancy_beta << endl;
          os << "Vorticity Confinement Epsilon = " << params.vort_conf_eps << endl;
          os << "Number of particles = " << params.num_particles << endl;
          os << "Iterations in solver = " << params.max_iterations << endl;
          return os;
        }
      public:
        Vector2Di dimension;
        double dt;
        double buoyancy_alpha, buoyancy_beta;
        double vort_conf_eps;
        double tolerance;
        int max_iterations;
        int num_particles;
      };
    /***************************** User Parameters END *************************/

    public:
      EulerGAS2D();
      EulerGAS2D(Parameters _user_params);
      EulerGAS2D(const EulerGAS2D& src);
      EulerGAS2D& operator=(const EulerGAS2D& rhs);
      ~EulerGAS2D();
    public:
      // TODO: Implement overload functions
      bool init(Parameters params);
      void step();
      void close();
      void add_source(int i, int j);
      void add_particles(VFXEpoch::Particle2D p);
    public:
      void set_boundary(Grid2DCellTypes boundaries);
    public:
      void set_user_params(Parameters params);
      EulerGAS2D::Parameters get_user_params();
    protected:
      void diffuse(Grid2DfScalarField& dest, Grid2DfScalarField ref);
      void advect(Grid2DfScalarField& dest, Grid2DfScalarField ref);
      void advect_particles();
      void trace_rkii();
      void compute_curls();
      void compute_buoyancy();
      void presure_solve();
      VFXEpoch::Vector2Dd get_vel();
    private:
    /*********************** Pressure Solver Parameters ************************/
      struct PressureSolverParams{
        PCGSolver<double> pcg_solver;
        SparseMatrixd sparse_matrix;
        vector<double> rhs;
        vector<double> pressure;
        
        inline void clear(){
          pcg_solver.clear();
          sparse_matrix.clear();
          rhs.clear();
          pressure.clear();
        }
      };
    /*********************** Pressure Solver Parameters END ********************/
    private:
      Grid2DfScalarField u, u0;
      Grid2DfScalarField v, v0;
      Grid2DfScalarField d, d0;
      Grid2DfScalarField t, t0;
      Grid2DfScalarField pressure;
      Grid2DfScalarField omega, omega0;
      Grid2DCellTypes domain_mask;
      vector<VFXEpoch::Particle2D> particles_container;
      Parameters user_params;
      PressureSolverParams pressure_solver_params;
    };
  }
}

#endif
