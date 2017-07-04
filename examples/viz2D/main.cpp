/*******************************************************************************
    VFXEPOCH - Physically based simulation VFX

    Copyright (c) 2016 Snow Tsui <trevor.miscellaneous@gmail.com>

    All rights reserved. Use of this source code is governed by
    the MIT license as written in the LICENSE file.
*******************************************************************************/
#include <cstdio>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <cfloat>

#include "VisualizerHelpers.h"
#include "VisualizerHelpers_OpenGL.h"
#include "Helpers.h"

using namespace VFXEpoch::Gluvi;
using namespace VFXEpoch::OpenGL_Utility;

//PanZoom2D cam(-0.1, -0.35, 1.2);
PanZoom2D cam(0.0, 0.0, 1.0f);
bool load_bin(const char* filename);
void init_data();
void display();
void mouse(int button, int state, int x, int y);
void drag(int x, int y);
void timer(int junk);

std::vector<VFXEpoch::Vector2Df> particles;

int main(int argc, char **argv)
{   
   //Setup viewer stuff
   Gluvi::init("2D Particles Visualizer", &argc, argv, 720, 280);
   init_data();
   Gluvi::camera=&cam;
   Gluvi::userDisplayFunc=display;
   Gluvi::userMouseFunc=mouse;
   Gluvi::userDragFunc=drag;
   glClearColor(0,0,0,1);   
   glutTimerFunc(1000, timer, 0);
   Gluvi::run();

   return 0;
}

void
init_data(){
    // TODO: Initialize data;
    VFXEpoch::Vector2Df v(0.5, 0.5);
    particles.push_back(v);
}

bool
load_bin(const char* filename){
    // TODO: Load binary files
    return true;
}

void
display(){
    // TODO: Draw particles
    glPointSize(5);
    OpenGL_Utility::draw_particles2d(particles);
}

void 
mouse(int button, int state, int x, int y){
    // TODO: Mouse state
}

void
drag(int x, int y){
    // TODO: N/A
}

void 
timer(int junk){
    // TODO: N/A
}


