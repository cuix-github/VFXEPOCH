/*******************************************************************************
    VFXEPOCH - Physically based simulation VFX

    Copyright (c) 2016 Snow Tsui <trevor.miscellaneous@gmail.com>

    All rights reserved. Use of this source code is governed by
    the MIT license as written in the LICENSE file.
*******************************************************************************/
#include "SIM_EulerGAS.h"

// Public
EulerGAS2D::EulerGAS2D(){
  this->user_params.clear();
  this->u.clear(); this->u0.clear();
  this->v.clear(); this->v0.clear();
  this->d.clear(); this->d0.clear();
  this->t.clear(); this->t0.clear();
  this->pressure.clear();
  this->omega.clear(); this->omega0.clear();
  this->domain_mask.clear();
}

// Public
EulerGAS2D::EulerGAS2D(const EulerGAS2D& src){
  this->u = src.u; this->u0 = src.u0;
  this->v = src.v; this->v0 = src.v0;
  this->d = src.v; this->d0 = src.v0;
  this->t = src.v; this->t0 = src.v0;
  this->pressure = src.pressure;
  this->omega = src.v; this->omega0 = src.v0;
  this->user_params = src.user_params;
  domain_mask = src.domain_mask;
}

// Public
EulerGAS2D::EulerGAS2D(Parameters _user_params):user_params(_user_params){
  u.ResetDimension(_user_params.dimension.m_x + 1, _user_params.dimension.m_y); u0 = u;
  v.ResetDimension(_user_params.dimension.m_x, _user_params.dimension.m_y + 1); v0 = v;
  d.ResetDimension(_user_params.dimension.m_x, _user_params.dimension.m_y); d0 = d;
  t.ResetDimension(_user_params.dimension.m_x, _user_params.dimension.m_y); t0 = t;
  pressure.ResetDimension(_user_params.dimension.m_x, _user_params.dimension.m_y);
  omega.ResetDimension(_user_params.dimension.m_x + 2, _user_params.dimension.m_y + 2);
  omega0.ResetDimension(_user_params.dimension.m_x + 2, _user_params.dimension.m_y + 2);
  domain_mask.ResetDimension(_user_params.dimension.m_x, _user_params.dimension.m_y);
}

// Public
EulerGAS2D&
EulerGAS2D::operator=(const EulerGAS2D& rhs){
  this->u = rhs.u; this->u0 = rhs.u0;
  this->v = rhs.v; this->v0 = rhs.v0;
  this->d = rhs.v; this->d0 = rhs.v0;
  this->t = rhs.v; this->t0 = rhs.v0;
  this->pressure = rhs.pressure;
  this->omega = rhs.v; this->omega0 = rhs.v0;
  this->domain_mask = rhs.domain_mask;
  this->user_params = rhs.user_params;
  return *this;
}

// Public
EulerGAS2D::~EulerGAS2D(){
  this->close();
}

// Public
// Overload from SIM_Base.h -> class Euler_Fluid2D_Base
bool
EulerGAS2D::init(Parameters params){
  this->user_params = params;
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
  u.clear(); u0.clear();
  v.clear(); v0.clear();
  d.clear(); d0.clear();
  t.clear(); t0.clear();
  omega.clear(); omega0.clear();
  pressure.clear();
  domain_mask.clear();
  user_params.clear();
}

// Public
// Overload from SIM_Base.h -> class Euler_Fluid2D_Base
void
EulerGAS2D::add_source(int i, int j){
  /* TODO: code */
}

// Public
void
EulerGAS2D::set_boundary(/* arguments */){
  /* TODO: */
}

// Public
void
EulerGAS2D::set_user_params(Parameters params){
  this->user_params = params;
}

EulerGAS2D::Parameters
EulerGAS2D::get_user_params(){
  return this->user_params;
}

// Protected
void
EulerGAS2D::diffuse(Grid2DfScalarField& dest, Grid2DfScalarField ref){
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
EulerGAS2D::presure_solve(){
  /* TODO: code */
}
