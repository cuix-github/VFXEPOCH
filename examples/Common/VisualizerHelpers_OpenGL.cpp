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
OpenGL_Utility::draw_grid2d(const VFXEpoch::Vector2Df& origin, float dx, unsigned int nx, unsigned int ny) {
	float width = nx * dx;
	float height = ny * dx;

  glLineWidth(1);
  glBegin(GL_LINES);
  for(int i = 0; i != nx; i++){
   
    VFXEpoch::Vector2Df a(i * dx, 0);
    VFXEpoch::Vector2Df b(i * dx, height);
    VFXEpoch::Vector2Df oa = origin + a;
    VFXEpoch::Vector2Df ob = origin + b;
   
    double oa_d[2] = {oa[0], oa[1]};
    double ob_d[2] = {ob[0], ob[1]};

    glVertex2dv(oa_d);
    glVertex2dv(ob_d);
  }

  for(int j = 0; j != ny; j++){

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
OpenGL_Utility::draw_particles2d(const std::vector<VFXEpoch::Vector2Dd>& particles_container) {
   glBegin(GL_POINTS);
   for(unsigned int i = 0; i < particles_container.size(); ++i) {
       glVertex2f(particles_container[i].m_x, particles_container[i].m_y);
   }
   glEnd();
}

void 
OpenGL_Utility::draw_particles2d(const std::vector<VFXEpoch::Vector2Df>& particles_container) {
   glBegin(GL_POINTS);
   for(unsigned int i = 0; i < particles_container.size(); ++i) {
       glVertex2f(particles_container[i].m_x, particles_container[i].m_y);
   }
   glEnd();
}

void
OpenGL_Utility::draw_particles2d(const std::vector<VFXEpoch::Particle2Dd>& particles_container) {
  glBegin(GL_POINTS);
  for(unsigned int i = 0; i < particles_container.size(); ++i) {
    glVertex2f(particles_container[i].pos.m_x, particles_container[i].pos.m_y);
  }
  glEnd(); 
}

void
OpenGL_Utility::draw_particles2d(const std::vector<VFXEpoch::Particle2Df>& particles_container) {
  glBegin(GL_POINTS);
  for(unsigned int i = 0; i < particles_container.size(); ++i) {
    glVertex2f(particles_container[i].pos.m_x, particles_container[i].pos.m_y);
  }
  glEnd(); 
}

void 
OpenGL_Utility::draw_circle2d(const VFXEpoch::Vector2Df& center, double rad, int segs) {
  
  glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
  glLineWidth(2);
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