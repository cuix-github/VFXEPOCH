
#ifndef VISUALIZER_HELPERS_OPENGL_UTILS_H
#define VISUALIZER_HELPERS_OPENGL_UTILS_H

#include <vector>
//#include "vec.h"
//#include "array2.h"

#include "utl/UTL_Vector.h"
#include "utl/UTL_Grid.h"

namespace VFXEpoch{
    namespace OpenGL_Utility{
        void draw_grid2d(const VFXEpoch::Vector2Df& origin, float dx, float nx, float ny);
        void draw_particles2d(const std::vector<VFXEpoch::Vector2Df>& particles_container);
    }
}

#endif