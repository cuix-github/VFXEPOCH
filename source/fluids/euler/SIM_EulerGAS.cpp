/*******************************************************************************
    VFXEPOCH - Physically based simulation VFX

    Copyright (c) 2016 Snow Tsui <trevor.miscellaneous@gmail.com>

    All rights reserved. Use of this source code is governed by
    the MIT license as written in the LICENSE file.
*******************************************************************************/
#include "SIM_EulerGAS.h"

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

EulerGAS2D::EulerGAS2D(VFXEpoch::Vector2Di _dimension, REAL _dt, REAL _vort_conf_update, REAL _buoyancy_alpha, REAL _buoyancy_beta, int _max_iterations):dimension(_dimension), dt(_dt), tolerance(0.0), vort_conf_eps(_vort_conf_update), buoyancy_alpha(_buoyancy_alpha), buoyancy_beta(_buoyancy_beta), max_iterations(_max_iterations){
  u.ResetDimension(_dimension.m_x + 1, _dimension.m_y); u0 = u;
  v.ResetDimension(_dimension.m_x, _dimension.m_y + 1); v0 = v;
  d.ResetDimension(_dimension.m_x, _dimension.m_y); d0 = d;
  t.ResetDimension(_dimension.m_x, _dimension.m_y); t0 = t;
  pressure.ResetDimension(_dimension.m_x, _dimension.m_y);
  omega.ResetDimension(_dimension.m_x + 2, _dimension.m_y + 2); omega0 = omega;
  domain_mask.ResetDimension(_dimension.m_x, _dimension.m_y);
}

EulerGAS2D&
EulerGAS2D::operator=(const EulerGAS2D& rhs){
  return *this;
}

EulerGAS2D::~EulerGAS2D(){
  this->close();
}
