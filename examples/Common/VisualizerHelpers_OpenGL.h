
#ifndef VISUALIZER_HELPERS_OPENGL_UTILS_H
#define VISUALIZER_HELPERS_OPENGL_UTILS_H

#include <vector>
//#include "vec.h"
//#include "array2.h"

#include "utl/UTL_Vector.h"
#include "utl/UTL_Grid.h"
#include "utl/UTL_General.h"

namespace VFXEpoch{
    namespace OpenGL_Utility{
        void draw_grid2d(const VFXEpoch::Vector2Df& origin, float dx, unsigned int nx, unsigned int ny);
        void draw_particles2d(const std::vector<VFXEpoch::Vector2Dd>& particles_container);
        void draw_particles2d(const std::vector<VFXEpoch::Vector2Df>& particles_container);
        void draw_particles2d(const std::vector<VFXEpoch::Particle2Dd>& particles_container);
        void draw_particles2d(const std::vector<VFXEpoch::Particle2Df>& particles_container);
        void draw_circle2d(const VFXEpoch::Vector2Df& center, double rad, int segs);
    }
}

#endif