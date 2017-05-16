/*******************************************************************************
    VFXEPOCH - Physically based simulation VFX

    Copyright (c) 2016 Snow Tsui <trevor.miscellaneous@gmail.com>

    All rights reserved. Use of this source code is governed by
    the MIT license as written in the LICENSE file.
*******************************************************************************/
#include "SIM_EulerGAS.h"

// Public
EulerGAS2D::EulerGAS2D(){
  user_params.clear();
  u.clear(); u0.clear();
  v.clear(); v0.clear();
  d.clear(); d0.clear();
  t.clear(); t0.clear();
  pressure.clear();
  omega.clear(); omega0.clear();
  inside_mask.clear();
  pressure_solver_params.clear();
  particles_container.clear();
  domain_boundaries[0].side = EDGES_2DSIM::TOP;
	domain_boundaries[0].boundaryType = BOUNDARY::STREAK;
	domain_boundaries[1].side = EDGES_2DSIM::BOTTOM;
	domain_boundaries[1].boundaryType = BOUNDARY::STREAK;
	domain_boundaries[2].side = EDGES_2DSIM::RIGHT;
	domain_boundaries[2].boundaryType = BOUNDARY::STREAK;
	domain_boundaries[3].side = EDGES_2DSIM::LEFT;
	domain_boundaries[3].boundaryType = BOUNDARY::STREAK;
}

// Public
EulerGAS2D::EulerGAS2D(const EulerGAS2D& src){
  u = src.u; u0 = src.u0;
  v = src.v; v0 = src.v0;
  d = src.v; d0 = src.v0;
  t = src.v; t0 = src.v0;
  pressure = src.pressure;
  omega = src.v; omega0 = src.v0;
  user_params = src.user_params;
  inside_mask = src.inside_mask;
  particles_container = src.particles_container;
}

// Public
EulerGAS2D::EulerGAS2D(Parameters _user_params):user_params(_user_params){
  v.ResetDimension(_user_params.dimension.m_x, _user_params.dimension.m_y + 1); v0 = v;
  u.ResetDimension(_user_params.dimension.m_x + 1, _user_params.dimension.m_y); u0 = u;
  d.ResetDimension(_user_params.dimension.m_x, _user_params.dimension.m_y); d0 = d;
  t.ResetDimension(_user_params.dimension.m_x, _user_params.dimension.m_y); t0 = t;
  pressure.ResetDimension(_user_params.dimension.m_x, _user_params.dimension.m_y);
  omega.ResetDimension(_user_params.dimension.m_x + 2, _user_params.dimension.m_y + 2);
  omega0.ResetDimension(_user_params.dimension.m_x + 2, _user_params.dimension.m_y + 2);
  inside_mask.ResetDimension(_user_params.dimension.m_x, _user_params.dimension.m_y);
  particles_container.resize(_user_params.num_particles);
}

// Public
EulerGAS2D&
EulerGAS2D::operator=(const EulerGAS2D& rhs){
  u = rhs.u; u0 = rhs.u0;
  v = rhs.v; v0 = rhs.v0;
  d = rhs.v; d0 = rhs.v0;
  t = rhs.v; t0 = rhs.v0;
  pressure = rhs.pressure;
  omega = rhs.v; omega0 = rhs.v0;
  inside_mask = rhs.inside_mask;
  user_params = rhs.user_params;
  particles_container = rhs.particles_container;
  return *this;
}

// Public
EulerGAS2D::~EulerGAS2D(){
  close();
}

// Public
// Overload from SIM_Base.h -> class Euler_Fluid2D_Base
bool
EulerGAS2D::init(Parameters params){
  user_params = params;
  v.ResetDimension(user_params.dimension.m_x, user_params.dimension.m_y + 1); v0 = v;
  u.ResetDimension(user_params.dimension.m_x + 1, user_params.dimension.m_y); u0 = u;
  d.ResetDimension(user_params.dimension.m_x, user_params.dimension.m_y); d0 = d;
  t.ResetDimension(user_params.dimension.m_x, user_params.dimension.m_y); t0 = t;
  pressure.ResetDimension(user_params.dimension.m_x, user_params.dimension.m_y);
  omega.ResetDimension(user_params.dimension.m_x + 2, user_params.dimension.m_y + 2); omega0 = omega;
  inside_mask.ResetDimension(user_params.dimension.m_x, user_params.dimension.m_y);
  particles_container.resize(user_params.num_particles);
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
  inside_mask.clear();
  user_params.clear();
  pressure_solver_params.clear();
  particles_container.clear();
}

// Public
// Overload from SIM_Base.h -> class Euler_Fluid2D_Base
void
EulerGAS2D::add_source(int i, int j){
  assert(i >= 0 && i <= d.getDimY() && j >= 0 && j <= d.getDimX());
  d(i, j) = user_params.density_source;
}

// Public
void
EulerGAS2D::add_source(VFXEpoch::VECTOR_COMPONENTS component, int i, int j){
  if(VFXEpoch::VECTOR_COMPONENTS::X == component){
    assert(i >= 0 && i <= u.getDimY() && j >= 0 && j <= u.getDimX());
    u(i, j) = user_params.external_force_strength;
  }
  else if(VFXEpoch::VECTOR_COMPONENTS::Y == component){
    assert(i >= 0 && i <= v.getDimY() && j >= 0 && j <= v.getDimX());
    v(i, j) = user_params.external_force_strength;
  }
  else{
    std::cout << "Component error occured!" << endl;
  }
}

// Public
void
EulerGAS2D::add_particles(VFXEpoch::Particle2D p){
  particles_container.push_back(p);
}

// Public
void 
EulerGAS2D::set_inside_boundary(Grid2DCellTypes boundaries){
  assert(inside_mask.getDimX() == boundaries.getDimX() && inside_mask.getDimY() == boundaries.getDimY());
  inside_mask = boundaries;
}

// Public
void 
EulerGAS2D::set_domain_boundary(const VFXEpoch::BOUNDARY boundary_type, const VFXEpoch::EDGES_2DSIM edge){
  if(VFXEpoch::EDGES_2DSIM::TOP == edge){
    domain_boundaries[0].side = edge;
    domain_boundaries[0].boundaryType = boundary_type;
  } else if(VFXEpoch::EDGES_2DSIM::BOTTOM == edge){
    domain_boundaries[1].side = edge;
    domain_boundaries[1].boundaryType = boundary_type;
  } else if(VFXEpoch::EDGES_2DSIM::LEFT == edge){
    domain_boundaries[2].side = edge;
    domain_boundaries[2].boundaryType = boundary_type;
  } else{
    domain_boundaries[3].side = edge;
    domain_boundaries[3].boundaryType = boundary_type;    
  }
}

// Public
void
EulerGAS2D::set_user_params(Parameters params){
  user_params = params;
}

EulerGAS2D::Parameters
EulerGAS2D::get_user_params(){
  return user_params;
}

// Protected
void 
EulerGAS2D::set_domain_boundary_wrapper(Grid2DfScalarField& field){
  for(int i = 0; i != 4; i++){
    field.setBoundaries(domain_boundaries[i].boundaryType, domain_boundaries[i].side);
  }
  field.setBoundariesOnCorners();
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
void
EulerGAS2D::advect_particles(){
  /* TODO: code */
}

//Protected
void
EulerGAS2D::trace_rkii(const VFXEpoch::Vector2Df& pos, double dt){
  /* TODO: code */
}

// Protected
// Overload from SIM_Base.h -> class Euler_Fluid2D_Base
void
EulerGAS2D::presure_solve(){
  /* TODO: code */
}

// Protected
void
EulerGAS2D::apply_gradients(){
  /* TODO: code */
}

// Protected
VFXEpoch::Vector2Dd
EulerGAS2D::get_vel(){
  VFXEpoch::Vector2Dd result;
  
  /* TODO: code */
  return result;
}
