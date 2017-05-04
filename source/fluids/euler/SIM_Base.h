/*******************************************************************************
    VFXEpoch - Physically based simulation VFX

    Copyright (c) 2016 Snow Tsui <trevor.miscellaneous@gmail.com>

    All rights reserved. Use of this source code is governed by
    the MIT license as written in the LICENSE file.
*******************************************************************************/
#ifndef _SIM_BASE_H_
#define _SIM_BASE_H_
#include <vector>
#include "utl/UTL_Vector.h"

using namespace std;

namespace VFXEpoch{
  namespace Solvers{
    class Euler_Fluid2D_Base{
    public:
      struct Particle2D{
        Particle2D(){}
        Particle2D(const Particle2D& src):pos(src.pos), vel(src.vel), radius(src.radius){}

        inline
        Particle2D& operator=(const Particle2D& rhs){
          pos = rhs.pos;
          vel = rhs.vel;
          radius = rhs.radius;
          return *this;
        }

        VFXEpoch::Vector2Df pos;
        VFXEpoch::Vector2Df vel;
        static long long instance_counter;
        long long id = instance_counter++;
        int radius = 0.5f;
      };

      virtual bool init(){}
      virtual void step(double dt){}
      virtual void advect(){}
      virtual void presure_solve(double dt){}
      virtual void close(){}

    protected:
      std::vector<Particle2D> particles;
    };

    class Euler_Fluid3D_Base{
    public:
      struct Particle3D{
        Particle3D(){}
        Particle3D(const Particle3D& src):pos(src.pos), vel(src.vel), radius(src.radius){}

        inline
        Particle3D& operator=(const Particle3D& rhs){
          pos = rhs.pos;
          vel = rhs.vel;
          radius = rhs.radius;
          return *this;
        }

        VFXEpoch::Vector3Df pos;
        VFXEpoch::Vector3Df vel;
        static long long instance_counter;
        long long id = instance_counter++;
        int radius = 0.5f;
      };

      virtual bool init(){}
      virtual void step(double dt){}
      virtual void advect(){}
      virtual void presure_solve(double dt){}
      virtual void close(){}

    protected:
      std::vector<Particle3D> particles;
    };
  }
}

#endif
