/*******************************************************************************
    VFXEPOCH - Physically based simulation VFX

    Copyright (c) 2016 Snow Tsui <trevor.miscellaneous@gmail.com>

    All rights reserved. Use of this source code is governed by
    the MIT license as written in the LICENSE file.
*******************************************************************************/
#include "SIM_EulerGAS.h"

// Public
EulerGAS2D::EulerGAS2D(){
  dimension.m_x = dimension.m_y = 0.0;
  buoyancy_alpha = buoyancy_beta = 0.0;
  vort_conf_eps = 0.0;
  dt = 0.0;
  tolerance = 0.0;
  max_iterations = 0;
  this->u.clear(); this->u0.clear();
  this->v.clear(); this->v0.clear();
  this->d.clear(); this->d0.clear();
  this->t.clear(); this->t0.clear();
  this->pressure.clear();
  this->omega.clear(); this->omega0.clear();
  domain_mask.clear();
}

// Public
EulerGAS2D::EulerGAS2D(const EulerGAS2D& src){
  this->dimension = src.dimension;
  this->u = src.u; this->u0 = src.u0;
  this->v = src.v; this->v0 = src.v0;
  this->d = src.v; this->d0 = src.v0;
  this->t = src.v; this->t0 = src.v0;
  this->pressure = src.pressure;
  this->omega = src.v; this->omega0 = src.v0;
  this->buoyancy_alpha = src.buoyancy_alpha;
  this->buoyancy_beta = src.buoyancy_beta;
  this->vort_conf_eps = src.vort_conf_eps;
  this->dt = src.dt;
  this->tolerance = src.tolerance;
  this->max_iterations = src.max_iterations;
  domain_mask = src.domain_mask;
}

// Public
EulerGAS2D::EulerGAS2D(Vector2Di _dimension, REAL _dt, REAL _vort_conf_update, REAL _buoyancy_alpha, REAL _buoyancy_beta, int _max_iterations):dimension(_dimension), dt(_dt), tolerance(0.0), vort_conf_eps(_vort_conf_update), buoyancy_alpha(_buoyancy_alpha), buoyancy_beta(_buoyancy_beta), max_iterations(_max_iterations){
  u.ResetDimension(_dimension.m_x + 1, _dimension.m_y); u0 = u;
  v.ResetDimension(_dimension.m_x, _dimension.m_y + 1); v0 = v;
  d.ResetDimension(_dimension.m_x, _dimension.m_y); d0 = d;
  t.ResetDimension(_dimension.m_x, _dimension.m_y); t0 = t;
  pressure.ResetDimension(_dimension.m_x, _dimension.m_y);
  omega.ResetDimension(_dimension.m_x + 2, _dimension.m_y + 2); omega0 = omega;
  domain_mask.ResetDimension(_dimension.m_x, _dimension.m_y);
}

// Public
EulerGAS2D&
EulerGAS2D::operator=(const EulerGAS2D& rhs){
  this->dimension = rhs.dimension;
  this->u = rhs.u; this->u0 = rhs.u0;
  this->v = rhs.v; this->v0 = rhs.v0;
  this->d = rhs.v; this->d0 = rhs.v0;
  this->t = rhs.v; this->t0 = rhs.v0;
  this->pressure = rhs.pressure;
  this->omega = rhs.v; this->omega0 = rhs.v0;
  this->buoyancy_alpha = rhs.buoyancy_alpha;
  this->buoyancy_beta = rhs.buoyancy_beta;
  this->vort_conf_eps = rhs.vort_conf_eps;
  this->dt = rhs.dt;
  this->tolerance = rhs.tolerance;
  this->max_iterations = rhs.max_iterations;
  domain_mask = rhs.domain_mask;
  return *this;
}

// Public
EulerGAS2D::~EulerGAS2D(){
  this->close();
}

// Public
// Overload from SIM_Base.h -> class Euler_Fluid2D_Base
bool
EulerGAS2D::init(/* arguments */){
  /* TODO: code */
  return true;
}

// Public
// Overload from SIM_Base.h -> class Euler_Fluid2D_Base
void
EulerGAS2D::step(){
  /* TODO: code */
}

// Public
// Overload from SIM_Base.h -> class Euler_Fluid2D_Base
void
EulerGAS2D::close(){
  /* TODO: code */
}

// Public
// Overload from SIM_Base.h -> class Euler_Fluid2D_Base
void
EulerGAS2D::add_source(/* arguments */){
  /* TODO: code */
}

// Public
void
EulerGAS2D::set_boundary(/* arguments */){
  /* TODO: */
}

// Public
void
EulerGAS2D::set_user_params(Parameters params) {
  this->dimension = params.dimension;
  this->buoyancy_alpha = params.buoyancy_alpha;
  this->buoyancy_beta = params.buoyancy_beta;
  this->vort_conf_eps = params.vort_conf_eps;
  this->dt = params.dt;
  this->tolerance = params.tolerance;
  this->max_iterations = params.max_iterations;
}

EulerGAS2D::Parameters
EulerGAS2D::get_user_params(){
  EulerGAS2D::Parameters p;
  p.dimension = this->dimension;
  p.buoyancy_alpha = this->buoyancy_alpha;
  p.buoyancy_beta = this->buoyancy_beta;
  p.vort_conf_eps = this->vort_conf_eps;
  p.dt = this->dt;
  p.tolerance = this->tolerance;
  p.max_iterations = this->max_iterations;
  return p;
}

// Protected
void
EulerGAS2D::diffuse(/* arguments */) {
  /* TODO: code */
}

// Protected
// Overload from SIM_Base.h -> class Euler_Fluid2D_Base
void
EulerGAS2D::advect(Grid2DfScalarField& dest, Grid2DfScalarField ref){
  /* TODO: code */
}

// Protected
// Overload from SIM_Base.h -> class Euler_Fluid2D_Base
void
EulerGAS2D::advect(Grid2DVector2DfField& dest, Grid2DVector2DfField ref){
  /* TODO: code */
}

// Protected
// Overload from SIM_Base.h -> class Euler_Fluid2D_Base
void
EulerGAS2D::presure_solve(double dt){
  /* TODO: code */
}
