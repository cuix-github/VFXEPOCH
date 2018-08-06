#include "VisualizerHelpers_OpenGL.h"

using namespace VFXEpoch;

#ifdef __APPLE__ 
  #include <GLUT/glut.h> // why does Apple have to put glut.h here...
#else
  #include <GL/glut.h> // ...when everyone else puts it here?
#endif

//#include "vec.h"
#include <cfloat>

//Debug
#include <iostream>
using namespace std;

// TODO: Synchronize with VFXEPOCH library
void 
OpenGL_Utility::draw_grid2d(const VFXEpoch::Vector2Df& origin, float dx, unsigned int nx, unsigned int ny, VFXEpoch::Vector3Df color) {
	glColor3f(color.m_x, color.m_y, color.m_z);
  float width = nx * dx;
	float height = ny * dx;

  glLineWidth(1);
  glBegin(GL_LINES);
  for(int i = 0; i <= nx; i++){
   
    VFXEpoch::Vector2Df a(i * dx, 0);
    VFXEpoch::Vector2Df b(i * dx, height);
    VFXEpoch::Vector2Df oa = origin + a;
    VFXEpoch::Vector2Df ob = origin + b;
   
    double oa_d[2] = {oa[0], oa[1]};
    double ob_d[2] = {ob[0], ob[1]};

    glVertex2dv(oa_d);
    glVertex2dv(ob_d);
  }

  for(int j = 0; j <= ny; j++){

    VFXEpoch::Vector2Df a(0, j * dx);
    VFXEpoch::Vector2Df b(width, j * dx);
    VFXEpoch::Vector2Df oa = origin + a;
    VFXEpoch::Vector2Df ob = origin + b;

    double oa_d[2] = {oa[0], oa[1]};
    double ob_d[2] = {ob[0], ob[1]};
    
    glVertex2dv(oa_d);
    glVertex2dv(ob_d);
  }
	glEnd();
}

void 
OpenGL_Utility::draw_particles2d(const std::vector<VFXEpoch::Vector2Dd>& particles_container, 
  int particle_size, VFXEpoch::Vector3Df color, bool is_round_point) {

  glColor3f(color.m_x, color.m_y, color.m_z);
  glPointSize(particle_size);

  if (is_round_point)
    glEnable( GL_POINT_SMOOTH );

  glBegin(GL_POINTS);
  for(unsigned int i = 0; i < particles_container.size(); ++i) {
    double data[2] = { particles_container[i].m_x, particles_container[i].m_y };
    glVertex2dv(data);
  }
  glEnd();
}

void 
OpenGL_Utility::draw_particles2d(const std::vector<VFXEpoch::Vector2Df>& particles_container, 
  int particle_size, VFXEpoch::Vector3Df color, bool is_round_point) {
  
  glColor3f(color.m_x, color.m_y, color.m_z);
  glPointSize(particle_size);

  if (is_round_point)
    glEnable( GL_POINT_SMOOTH );

  glBegin(GL_POINTS);
  for(unsigned int i = 0; i < particles_container.size(); ++i) {
      double data[2] = { particles_container[i].m_x, particles_container[i].m_y };
      glVertex2dv(data);
  }
  glEnd();
}

void
OpenGL_Utility::draw_particles2d(const std::vector<VFXEpoch::Particle2Dd>& particles_container, 
  int particle_size, VFXEpoch::Vector3Df color, bool is_round_point) {
  
  glColor3f(color.m_x, color.m_y, color.m_z);
  glPointSize(particle_size);

  if (is_round_point)
    glEnable( GL_POINT_SMOOTH );

  glBegin(GL_POINTS);
  for(unsigned int i = 0; i < particles_container.size(); ++i) {
    double data[2] = { particles_container[i].pos.m_x, particles_container[i].pos.m_y };
    glVertex2dv(data);
  }
  glEnd(); 
}

void
OpenGL_Utility::draw_particles2d(const std::vector<VFXEpoch::Particle2Df>& particles_container, 
  int particle_size, VFXEpoch::Vector3Df color, bool is_round_point) {

  glColor3f(color.m_x, color.m_y, color.m_z);
  glPointSize(particle_size);

  if (is_round_point)
    glEnable( GL_POINT_SMOOTH );

  glBegin(GL_POINTS);
  for(unsigned int i = 0; i < particles_container.size(); ++i) {
    double data[2] = { particles_container[i].pos.m_x, particles_container[i].pos.m_y };
    glVertex2dv(data);
  }
  glEnd(); 
}

void 
OpenGL_Utility::draw_arrows(VFXEpoch::Solvers::EulerGAS2D* solver, float arrow_len, VFXEpoch::Vector3Df color){
  if(!solver){
    std::cout << "ERROR: solver is not initialized" << std::endl;
    return;
  }

  glColor3f(color.m_x, color.m_y, color.m_z);
  VFXEpoch::Solvers::EulerGAS2D::Parameters user_params = solver->get_user_params();
  int row = user_params.dimension.m_x;
  int col = user_params.dimension.m_y;

  VFXEpoch::Vector2Df start, end;

  LOOP_GRID2D(row, col) {
    VFXEpoch::Vector2Df pos = VFXEpoch::Vector2Df((j + 0.5) * user_params.h, (i + 0.5) * user_params.h) + user_params.origin;
    start = pos;
    end = pos + arrow_len * solver->get_grid_velocity(pos);
    VFXEpoch::Vector2Df direction = end - start;
    VFXEpoch::Vector2Df dir_norm = direction;    

    //TODO Possibly automatically scale arrowhead length based on vector magnitude
    // if(dir_norm.norm() < 1e-30)
    //   return;

    dir_norm.normalize();
    VFXEpoch::Vector2Df perp(dir_norm[1], -dir_norm[0]);  
    VFXEpoch::Vector2Df tip_left = end + arrow_len / (float)sqrt(2.0) * ( -dir_norm + perp );
    VFXEpoch::Vector2Df tip_right = end + arrow_len / (float)sqrt(2.0) * ( -dir_norm - perp );

    glBegin(GL_LINES);
    glVertex2f(start[0], start[1]);
    glVertex2f(end[0], end[1]);
    glVertex2f(end[0], end[1]);
    glVertex2f(tip_left[0], tip_left[1]);
    glVertex2f(end[0], end[1]);
    glVertex2f(tip_right[0], tip_right[1]);
    glEnd();
  }
}

void 
OpenGL_Utility::draw_circle2d(const VFXEpoch::Vector2Df& center, double rad, int segs, VFXEpoch::Vector3Df color) {
  glColor3f(color.m_x, color.m_y, color.m_z);
  glLineWidth(2);
  glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
  glBegin(GL_POLYGON);
  for(int i = 0; i < segs; i++){
    double cosine = rad * cos(i * 2 * M_PI / (double)(segs));
    double sine = rad * sin(i * 2 * M_PI / (double)(segs));

    VFXEpoch::Vector2Df tmp = VFXEpoch::Vector2Df(cosine, sine) + center;
    double vec[2] = {tmp.m_x, tmp.m_y};
    glVertex2dv(vec);
  }
  glEnd();
}


