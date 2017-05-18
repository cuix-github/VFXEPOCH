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
  uw.clear(); vw.clear();
  d.clear(); d0.clear();
  t.clear(); t0.clear();
  omega.clear(); omega0.clear();
  inside_mask.clear(); inside_mask0.clear();
  nodal_solid_phi.clear();
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
  uw = src.uw; vw = src.vw;
  d = src.v; d0 = src.v0;
  t = src.v; t0 = src.v0;
  omega = src.v; omega0 = src.v0;
  user_params = src.user_params;
  inside_mask = src.inside_mask; inside_mask0 = src.inside_mask0;
  nodal_solid_phi = src.nodal_solid_phi;
  particles_container = src.particles_container;
}

// Public
EulerGAS2D::EulerGAS2D(Parameters _user_params):user_params(_user_params){
  v.ResetDimension(_user_params.dimension.m_x, _user_params.dimension.m_y + 1); v0 = v;
  u.ResetDimension(_user_params.dimension.m_x + 1, _user_params.dimension.m_y); u0 = u;
  uw.ResetDimension(_user_params.dimension.m_x, _user_params.dimension.m_y + 1);
  vw.ResetDimension(_user_params.dimension.m_x + 1, _user_params.dimension.m_y + 1);
  d.ResetDimension(_user_params.dimension.m_x, _user_params.dimension.m_y); d0 = d;
  t.ResetDimension(_user_params.dimension.m_x, _user_params.dimension.m_y); t0 = t;
  omega.ResetDimension(_user_params.dimension.m_x + 2, _user_params.dimension.m_y + 2);
  omega0.ResetDimension(_user_params.dimension.m_x + 2, _user_params.dimension.m_y + 2);
  nodal_solid_phi.ResetDimension(_user_params.dimension.m_x + 1, _user_params.dimension.m_y + 1);
  inside_mask.ResetDimension(_user_params.dimension.m_x + 1, _user_params.dimension.m_y + 1); inside_mask0 = inside_mask;
  particles_container.resize(_user_params.num_particles);
}

// Public
EulerGAS2D&
EulerGAS2D::operator=(const EulerGAS2D& rhs){
  u = rhs.u; u0 = rhs.u0;
  v = rhs.v; v0 = rhs.v0;
  uw = rhs.uw; vw = rhs.vw;
  d = rhs.v; d0 = rhs.v0;
  t = rhs.v; t0 = rhs.v0;
  omega = rhs.v; omega0 = rhs.v0;
  inside_mask = rhs.inside_mask; inside_mask0 = rhs.inside_mask0;
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
  uw.ResetDimension(user_params.dimension.m_x, user_params.dimension.m_y + 1);
  vw.ResetDimension(user_params.dimension.m_x + 1, user_params.dimension.m_x);
  d.ResetDimension(user_params.dimension.m_x, user_params.dimension.m_y); d0 = d;
  t.ResetDimension(user_params.dimension.m_x, user_params.dimension.m_y); t0 = t;
  omega.ResetDimension(user_params.dimension.m_x + 2, user_params.dimension.m_y + 2); omega0 = omega;
  inside_mask.ResetDimension(user_params.dimension.m_x, user_params.dimension.m_y);
  nodal_solid_phi.ResetDimension(user_params.dimension.m_x + 1, user_params.dimension.m_y + 1);
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
  inside_mask.clear(); inside_mask0.clear();
  nodal_solid_phi.clear();
  user_params.clear();
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
EulerGAS2D::add_external_force(VFXEpoch::VECTOR_COMPONENTS component, int i, int j){
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
EulerGAS2D::density_diffuse(Grid2DfScalarField& dest, Grid2DfScalarField ref){
  float a = user_params.diff * user_params.dt * user_params.dimension.m_x * user_params.dimension.m_y;
  VFXEpoch::LinearSolver::GSSolve(dest, ref, domain_boundaries, a, 1+4*a, user_params.max_iterations);
}

// Protected
void 
EulerGAS2D::dynamic_diffuse(Grid2DfScalarField& dest, Grid2DfScalarField ref){
  float a = user_params.visc * user_params.dt * user_params.dimension.m_x * user_params.dimension.m_y;
  VFXEpoch::LinearSolver::GSSolve(dest, ref, domain_boundaries, a, 1+4*a, user_params.max_iterations);
}

// Protected
// Overload from SIM_Base.h -> class Euler_Fluid2D_Base
void
EulerGAS2D::advect_vel(){
  // Using RK2 method time integration
  // advect u component of velocity field
  LOOP_GRID2D(u0){
    VFXEpoch::Vector2Df pos((j+0.5f) * user_params.size.m_x, i * user_params.size.m_y);
    pos = trace_rk2(pos, -user_params.dt);
    u0(i, j) = get_vel(pos).m_x;
  }

  // advect v component of velocity field
  LOOP_GRID2D(v0){
    VFXEpoch::Vector2Df pos(j * user_params.size.m_x, (i+0.5f) * user_params.size.m_y);
    pos = trace_rk2(pos, -user_params.dt);
    v0(i, j) = get_vel(pos).m_y;
  }
  u = u0;
  v = v0;
}

// Protected
void 
EulerGAS2D::advect_den(){
  // Using RK2 method time integration
  // advect density field
  // Brutal turning over the boundaries
  assert(d0.getDimX() == inside_mask.getDimX() && d0.getDimY() == inside_mask.getDimY());
  LOOP_GRID2D(d0){
    if(inside_mask(i, j) == VFXEpoch::BOUNDARY_MASK::SOMETHING) continue;
    else{
      VFXEpoch::Vector2Df pos((j+0.5f) * user_params.size.m_x, (i+0.5f) * user_params.size.m_y);
      pos = trace_rk2(pos, -user_params.dt);
      d0(i, j) = get_den(pos);
    }    
  }
  d = d0;
}

// Protected
void
EulerGAS2D::advect_curl(){
  // Using RK2 method time integration
  // advect curl field
  LOOP_GRID2D(omega0){
    VFXEpoch::Vector2Df pos((j+0.5f) * user_params.size.m_x, (i+0.5f) * user_params.size.m_y);
    pos = trace_rk2(pos, -user_params.dt);
    omega0(i, j) = get_curl(pos);
  }
  omega = omega0;
}

// Protected
void
EulerGAS2D::advect_particles(){
  std::vector<VFXEpoch::Particle2D>::iterator ite(0);
  for(ite = particles_container.begin(); ite != particles_container.end(); ite++){
    ite->pos = trace_rk2(ite->pos, user_params.dt);
  }
}

// Public
void
EulerGAS2D::apply_buoyancy(){
  /* TODO: code */
}

//Protected
Vector2Df
EulerGAS2D::trace_rk2(const Vector2Df& pos, float dt){
  Vector2Df vel = get_vel(pos);
  vel = get_vel(pos + 0.5f * dt * vel);
  return Vector2Df(pos + dt * vel);
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
void
EulerGAS2D::extrapolate(Grid2DfScalarField& grid, 
                        const Grid2DfScalarField& weights, 
                        Grid2DCellTypes& mask, 
                        Grid2DCellTypes& mask0){
  LOOP_GRID2D(grid){
    mask(i, j) = weights(i, j) > 0.0f ? BOUNDARY_MASK::NOTHING : BOUNDARY_MASK::SOMETHING;
  }

  for(int i = 0; i != 5; i++){
    mask0 = mask;
    LOOP_GRID2D(grid){
      float sum = 0;
      int count = 0;

      if(BOUNDARY_MASK::SOMETHING == mask0(i, j)){
        if(BOUNDARY_MASK::NOTHING == mask0(i, j+1)){
          sum += grid(i, j+1);
          ++count;
        } else if(BOUNDARY_MASK::NOTHING == mask0(i, j-1)){
          sum += grid(i, j-1);
          ++count;
        } else if(BOUNDARY_MASK::NOTHING == mask0(i+1, j)){
          sum += grid(i+1, j);
          ++count;
        } else if(BOUNDARY_MASK::NOTHING == mask0(i-1, j)){
          sum += grid(i-1, j);
          ++count;
        } else {
          /* TODO: ... */
        }

        if(count > 0){
          grid(i, j) = sum / (float)count;
          mask(i, j) = BOUNDARY_MASK::NOTHING;
        }
      }
    }
  }                       
}

// Protected
void
EulerGAS2D::get_grid_weights(){
  LOOP_GRID2D(uw){
    /* TODO: Compute u weights */
  }

  LOOP_GRID2D(vw){
    /* TODO: Compute v weights */
  }
}

// Protected
void
EulerGAS2D::constraint_vel(){
  /* TODO: code */
}

// Protected
Vector2Df
EulerGAS2D::get_vel(const Vector2Df& pos){
  assert(user_params.size.m_x != 0 && user_params.size.m_y != 0);
  float _u = VFXEpoch::InterpolateGrid(pos / user_params.size.m_x - Vector2Df(0.5f, 0.0f), u);
  float _v = VFXEpoch::InterpolateGrid(pos / user_params.size.m_y - Vector2Df(0.0f, 0.5f), v);
  return Vector2Df(_u, _v);
}

// Protected
float
EulerGAS2D::get_den(const Vector2Df& pos){
  assert(user_params.size.m_x != 0 && user_params.size.m_y != 0 && user_params.size.m_x == user_params.size.m_y);
  float h = user_params.size.m_x;
  return VFXEpoch::InterpolateGrid(pos / h - Vector2Df(0.5f, 0.5f), d);
}

// Protected
float
EulerGAS2D::get_curl(const Vector2Df& pos){
  assert(user_params.size.m_x != 0 && user_params.size.m_y != 0 && user_params.size.m_x == user_params.size.m_y);
  float h = user_params.size.m_x;
  return VFXEpoch::InterpolateGrid(pos / h - Vector2Df(0.5f, 0.5f), omega);
}
