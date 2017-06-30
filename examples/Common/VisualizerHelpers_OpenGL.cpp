#include "VisualizerHelpers_OpenGL.h"

#ifdef __APPLE__ 
  #include <GLUT/glut.h> // why does Apple have to put glut.h here...
#else
  #include <GL/glut.h> // ...when everyone else puts it here?
#endif

//#include "vec.h"
#include <cfloat>

// TODO: Synchronize with VFXEPOCH library
void draw_grid2d(const VFXEpoch::Vector2Df& origin, float dx, float nx, float ny) {
//    float width = nx*dx;
//    float height = ny*dx;
   
//    glBegin(GL_LINES);
//    for(int i = 0; i <= nx; i++) {
//       Vec2f a(i*dx, 0);
//       Vec2f b(i*dx, height);
//       glVertex2fv((origin+a).v); 
//       glVertex2fv((origin+b).v);
//    }
//    for(int j = 0; j <= ny; ++j) {
//       Vec2f a(0,j*dx);
//       Vec2f b(width,j*dx);
//       glVertex2fv((origin + a).v); 
//       glVertex2fv((origin + b).v);
//    }
//    glEnd();
}

void draw_particles2d(const std::vector<VFXEpoch::Vector2Df>& particles_container) {
   glBegin(GL_POINTS);
   for(unsigned int i = 0; i < particles_container.size(); ++i) {
       glVertex2f(particles_container[i].m_x, particles_container[i].m_y);
   }
   glEnd();
}

