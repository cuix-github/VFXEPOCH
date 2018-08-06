
#ifndef VISUALIZER_HELPERS_OPENGL_UTILS_H
#define VISUALIZER_HELPERS_OPENGL_UTILS_H

#include <vector>
//#include "vec.h"
//#include "array2.h"

#include "utl/UTL_Vector.h"
#include "utl/UTL_Grid.h"
#include "utl/UTL_General.h"
#include "Helpers.h"

namespace VFXEpoch{
    namespace OpenGL_Utility{
        void draw_grid2d(const VFXEpoch::Vector2Df& origin, float dx, unsigned int nx, unsigned int ny, VFXEpoch::Vector3Df color = VFXEpoch::Vector3Df(0.5, 0.5, 0.5));
        void draw_particles2d(const std::vector<VFXEpoch::Vector2Dd>& particles_container, int particle_size = 1, VFXEpoch::Vector3Df color = VFXEpoch::Vector3Df(0.0, 1.0, 0.0), bool is_round_point = false);
        void draw_particles2d(const std::vector<VFXEpoch::Vector2Df>& particles_container, int particle_size = 1, VFXEpoch::Vector3Df color = VFXEpoch::Vector3Df(0.0, 1.0, 0.0), bool is_round_point = false);
        void draw_particles2d(const std::vector<VFXEpoch::Particle2Dd>& particles_container, int particle_size = 1, VFXEpoch::Vector3Df color = VFXEpoch::Vector3Df(0.0, 1.0, 0.0), bool is_round_point = false);
        void draw_particles2d(const std::vector<VFXEpoch::Particle2Df>& particles_container, int particle_size = 1, VFXEpoch::Vector3Df color = VFXEpoch::Vector3Df(0.0, 1.0, 0.0), bool is_round_point = false);
        void draw_arrows(VFXEpoch::Solvers::EulerGAS2D* solver, float header_len, VFXEpoch::Vector3Df color = VFXEpoch::Vector3Df(0.7, 0.7, 0.7));
        void draw_circle2d(const VFXEpoch::Vector2Df& center, double rad, int segs, VFXEpoch::Vector3Df color = VFXEpoch::Vector3Df(0.7, 0.7, 0.7));
    }
}

#endif