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
  source_locations.clear();
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
  source_locations = src.source_locations;
}

// Public
EulerGAS2D::EulerGAS2D(Parameters _user_params):user_params(_user_params){
  v.Reset(_user_params.dimension.m_x, _user_params.dimension.m_y + 1, _user_params.h, _user_params.h); v0 = v;
  u.Reset(_user_params.dimension.m_x + 1, _user_params.dimension.m_y, _user_params.h, _user_params.h); u0 = u;
  uw.Reset(_user_params.dimension.m_x, _user_params.dimension.m_y + 1, _user_params.h, _user_params.h);
  vw.Reset(_user_params.dimension.m_x + 1, _user_params.dimension.m_y + 1, _user_params.h, _user_params.h);
  d.Reset(_user_params.dimension.m_x, _user_params.dimension.m_y, _user_params.h, _user_params.h); d0 = d;
  t.Reset(_user_params.dimension.m_x, _user_params.dimension.m_y, _user_params.h, _user_params.h); t0 = t;
  omega.Reset(_user_params.dimension.m_x + 2, _user_params.dimension.m_y + 2, _user_params.h, _user_params.h);
  omega0.Reset(_user_params.dimension.m_x + 2, _user_params.dimension.m_y + 2, _user_params.h, _user_params.h);
  nodal_solid_phi.Reset(_user_params.dimension.m_x + 1, _user_params.dimension.m_y + 1, _user_params.h, _user_params.h);
  inside_mask.Reset(_user_params.dimension.m_x + 1, _user_params.dimension.m_y + 1, _user_params.h, _user_params.h); inside_mask0 = inside_mask;
  particles_container.resize(_user_params.num_particles);
  source_locations.resize(0);
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
  source_locations = rhs.source_locations;
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
  v.Reset(user_params.dimension.m_x, user_params.dimension.m_y + 1, user_params.h, user_params.h); v0 = v;
  u.Reset(user_params.dimension.m_x + 1, user_params.dimension.m_y, user_params.h, user_params.h); u0 = u;
  uw.Reset(user_params.dimension.m_x, user_params.dimension.m_y + 1, user_params.h, user_params.h);
  vw.Reset(user_params.dimension.m_x + 1, user_params.dimension.m_x, user_params.h, user_params.h);
  d.Reset(user_params.dimension.m_x, user_params.dimension.m_y, user_params.h, user_params.h); d0 = d;
  t.Reset(user_params.dimension.m_x, user_params.dimension.m_y, user_params.h, user_params.h); t0 = t;
  omega.Reset(user_params.dimension.m_x + 2, user_params.dimension.m_y + 2, user_params.h, user_params.h); omega0 = omega;
  inside_mask.Reset(user_params.dimension.m_x, user_params.dimension.m_y, user_params.h, user_params.h);
  nodal_solid_phi.Reset(user_params.dimension.m_x + 1, user_params.dimension.m_y + 1, user_params.h, user_params.h);
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
  source_locations.clear();
}

// Public
void
EulerGAS2D::set_source_location(int i, int j){
  /* TODO: Try c++11 set function*/
  // Duplicates checking to avoid adding source at the same location
  assert(i >= 0 && j >= 0 && i < user_params.dimension.m_y && j < user_params.dimension.m_x);
  source_locations.push_back(VFXEpoch::Vector2Di(i, j));
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
// TODO: Check fast linear solvercorrectness for dx, dy, dimension in vertical & horizontal
void
EulerGAS2D::density_diffuse(Grid2DfScalarField& dest, Grid2DfScalarField ref){
  float a = user_params.diff * user_params.dt * user_params.dimension.m_x * user_params.dimension.m_y;
  VFXEpoch::LinearSolver::GSSolve(dest, ref, domain_boundaries, a, 1+4*a, user_params.max_iterations);
}

// Protected
// TODO: Check fast linear solver correctness for dx, dy, dimension in vertical & horizontal
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
    VFXEpoch::Vector2Df pos((j+0.5f) * user_params.h, i * user_params.h);
    pos = trace_rk2(pos, -user_params.dt);
    u0(i, j) = get_vel(pos).m_x;
  }

  // advect v component of velocity field
  LOOP_GRID2D(v0){
    VFXEpoch::Vector2Df pos(j * user_params.h, (i+0.5f) * user_params.h);
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
      VFXEpoch::Vector2Df pos((j+0.5f) * user_params.h, (i+0.5f) * user_params.h);
      pos = trace_rk2(pos, -user_params.dt);
      d0(i, j) = get_den(pos);
  }
  d = d0;
}

// Protected
void
EulerGAS2D::advect_tmp(){
  assert(t0.getDimX() == inside_mask.getDimX() && t0.getDimY() == inside_mask.getDimY());
  LOOP_GRID2D(t0){
      VFXEpoch::Vector2Df pos((j+0.5f) * user_params.h, (i+0.5f) * user_params.h);
      pos = trace_rk2(pos, -user_params.dt);
      t0(i, j) = get_tmp(pos);
  }
  t = t0;
}

// Protected
void
EulerGAS2D::advect_curl(){
  // Using RK2 method time integration
  // advect curl field
  assert(omega0.getDimX() == omega.getDimX() && omega0.getDimY() == omega0.getDimY());
  LOOP_GRID2D(omega0){
    VFXEpoch::Vector2Df pos((j+0.5f) * user_params.h, (i+0.5f) * user_params.h);
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

    // Correction particles at the boundaries
    float h = user_params.h;
    float corrections = VFXEpoch::InterpolateGrid(ite->pos / h, nodal_solid_phi);
    if(corrections < 0.0f){
      VFXEpoch::Vector2Df normal;
      VFXEpoch::InterpolateGradient(normal, ite->pos / h, nodal_solid_phi);
      normal.normalize();
      ite->pos -= corrections * normal;
    }
  }
}

// Protected
void
EulerGAS2D::project(){
  pressure_solve();
  apply_gradients();
}

// Protected
void
EulerGAS2D::apply_buoyancy(){
  // Dimension check
  assert(v0.getDimX() == v.getDimX() && v0.getDimY() == v.getDimY());

  float a = user_params.buoyancy_alpha;
  float b = user_params.buoyancy_beta;
  int row = user_params.dimension.m_y;
  int col = user_params.dimension.m_x;
  LOOP_GRID2D(v0){
    if(0 == i || j == 0) continue;
    float average_temperature = (t(i, j) + t(i, j - 1)) * 0.5f;
    float average_density = (d(i, j) + d(i, j - 1)) * 0.5f;
    v0(i, j) = v(i, j) - a * average_density + b * average_temperature;
  }
  v = v0;
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
EulerGAS2D::pressure_solve(){
  // System size check and resize;
  int system_size = user_params.dimension.m_x * user_params.dimension.m_y;
  if(pressure_solver_params.pressure.size() != system_size){
    pressure_solver_params.rhs.resize(system_size);
    pressure_solver_params.pressure.resize(system_size);
    pressure_solver_params.sparse_matrix.resize(system_size);
  }

  VFXEpoch::Grid2DdScalarField div(user_params.dimension.m_x, user_params.dimension.m_x, user_params.h, user_params.h);
  get_grid_weights();
  VFXEpoch::Analysis::computeDivergence_with_weights_mac(div, u, v, uw, vw);
  pressure_solver_params.rhs = div.toVector();
  setup_pressure_coef_matrix();
  bool success = pressure_solver_params.pcg_solver.solve(pressure_solver_params.sparse_matrix,
                                                         pressure_solver_params.rhs,
                                                         pressure_solver_params.pressure,
                                                         user_params.tolerance,
                                                         user_params.max_iterations);
  if(!success){
    std::cout << "Pressure solve failed!" << endl;
  }                                                         
}

// Protected
void
EulerGAS2D::apply_gradients(){
  VFXEpoch::Grid2DfScalarField _pressure(user_params.dimension.m_x, user_params.dimension.m_y);
  std::vector<float> solvedPressure(pressure_solver_params.pressure.begin(), pressure_solver_params.pressure.end());
  VFXEpoch::DataFromVectorToGrid(solvedPressure, _pressure);
  double dt = user_params.dt;
  double dx = user_params.h;
  LOOP_GRID2D(u){
    if(uw(i, j) > 0 && j > 0){
      u(i, j) -= dt * (_pressure(i, j) - _pressure(i, j - 1)) / dx;
    } else {
      u(i, j) = 0.0f;
    }
  }

  LOOP_GRID2D(v){
    if(vw(i, j) > 0 && i > 0){
      v(i, j) -= dt * (_pressure(i, j) - _pressure(i - 1, j)) / dx;
    } else {
      v(i, j) = 0.0f;
    }
  }
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
    uw(i, j) = 1 - VFXEpoch::InteralFrac(nodal_solid_phi(i+1, j), nodal_solid_phi(i, j));
  }

  LOOP_GRID2D(vw){
    vw(i, j) = 1 - VFXEpoch::InteralFrac(nodal_solid_phi(i, j+1), nodal_solid_phi(i, j));
  }
}

// Protected
void
EulerGAS2D::clamp_vel(){
  u0 = u; v0 = v;
  float h = user_params.h;
  LOOP_GRID2D(u){
    if(uw(i, j) == 0.0f){
      VFXEpoch::Vector2Df pos((i+0.5f) * h, j * h);
      VFXEpoch::Vector2Df vel = get_vel(pos);
      VFXEpoch::Vector2Df normal(0.0f, 0.0f);
      normal.normalize();
      VFXEpoch::InterpolateGradient(normal, pos / h, nodal_solid_phi);
      float correction_component = VFXEpoch::Vector2Df::dot(vel, normal);
      vel -= correction_component * normal;
      u0(i, j) = vel.m_x;
    }
  }

  LOOP_GRID2D(v){
    if(vw(i, j) == 0.0f){
      VFXEpoch::Vector2Df pos((i) * h, (j+0.5f) * h);
      VFXEpoch::Vector2Df vel = get_vel(pos);
      VFXEpoch::Vector2Df normal(0.0f, 0.0f);
      VFXEpoch::InterpolateGradient(normal, pos / h, nodal_solid_phi);
      normal.normalize();
      float correction_component = VFXEpoch::Vector2Df::dot(vel, normal);
      vel -= correction_component * normal;
      v0(i, j) = vel.m_y;
    }
  }
  u = u0; v = v0;
}

// Protected
void
EulerGAS2D::setup_pressure_coef_matrix(){
  int row = user_params.dimension.m_y - 1;
  int col = user_params.dimension.m_x - 1;
  int idx = 0;
  float val = 0.0f;
  float dx = user_params.h;
  float dt = user_params.dt;
  int grid_each_row_elements = user_params.dimension.m_x;
  LOOP_GRID2D_WITHOUT_DOMAIN_BOUNDARY(row, col){
    idx = i * user_params.dimension.m_x + j;
    val = uw(i, j+1) * dt / sqrt(dx);
    pressure_solver_params.sparse_matrix.add_to_element(idx, idx,val);
    pressure_solver_params.sparse_matrix.add_to_element(idx, idx + 1, -val);
    val = uw(i, j) * dt / sqrt(dx);
    pressure_solver_params.sparse_matrix.add_to_element(idx, idx, val);
    pressure_solver_params.sparse_matrix.add_to_element(idx, idx - 1, -val);
    val = vw(i+1, j) * dt / sqrt(dx);
    pressure_solver_params.sparse_matrix.add_to_element(idx, idx, val);
    pressure_solver_params.sparse_matrix.add_to_element(idx, idx + grid_each_row_elements, -val);
    val = vw(i, j) * dt / sqrt(dx);
    pressure_solver_params.sparse_matrix.add_to_element(idx, idx, val);
    pressure_solver_params.sparse_matrix.add_to_element(idx, idx - grid_each_row_elements, -val);
  }
}

// Protected
Vector2Df
EulerGAS2D::get_vel(const Vector2Df& pos){
  assert(user_params.h != 0);
  float _u = VFXEpoch::InterpolateGrid(pos / (float)user_params.h - Vector2Df(0.5f, 0.0f), u);
  float _v = VFXEpoch::InterpolateGrid(pos / (float)user_params.h - Vector2Df(0.0f, 0.5f), v);
  return Vector2Df(_u, _v);
}

// Protected
float
EulerGAS2D::get_den(const Vector2Df& pos){
  assert(user_params.h != 0);
  float h = user_params.h;
  return VFXEpoch::InterpolateGrid(pos / h - Vector2Df(0.5f, 0.5f), d);
}

// Protected
float
EulerGAS2D::get_curl(const Vector2Df& pos){
  assert(user_params.h != 0);
  float h = user_params.h;
  return VFXEpoch::InterpolateGrid(pos / h - Vector2Df(0.5f, 0.5f), omega);
}

// Protected
float
EulerGAS2D::get_tmp(const Vector2Df& pos){
  assert(user_params.h != 0);
  float h = user_params.h;
  return VFXEpoch::InterpolateGrid(pos / h - Vector2Df(0.5f, 0.5f), t);
}
