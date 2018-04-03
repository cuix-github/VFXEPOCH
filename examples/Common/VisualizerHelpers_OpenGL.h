
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
        void draw_particles2d(const std::vector<VFXEpoch::Vector2Dd>& particles_container, int particle_size = 1, VFXEpoch::Vector3Df color = VFXEpoch::Vector3Df(0.0, 1.0, 0.0), bool is_round_point = false);
        void draw_particles2d(const std::vector<VFXEpoch::Vector2Df>& particles_container, int particle_size = 1, VFXEpoch::Vector3Df color = VFXEpoch::Vector3Df(0.0, 1.0, 0.0), bool is_round_point = false);
        void draw_particles2d(const std::vector<VFXEpoch::Particle2Dd>& particles_container, int particle_size = 1, VFXEpoch::Vector3Df color = VFXEpoch::Vector3Df(0.0, 1.0, 0.0), bool is_round_point = false);
        void draw_particles2d(const std::vector<VFXEpoch::Particle2Df>& particles_container, int particle_size = 1, VFXEpoch::Vector3Df color = VFXEpoch::Vector3Df(0.0, 1.0, 0.0), bool is_round_point = false);
        void draw_arrows(const VFXEpoch::Vector2Df& start, const VFXEpoch::Vector2Df& end, float header_len);
        void draw_arrows(const VFXEpoch::Vector2Dd& start, const VFXEpoch::Vector2Dd& end, double header_len);
        void draw_circle2d(const VFXEpoch::Vector2Df& center, double rad, int segs);
    }
}

#endif