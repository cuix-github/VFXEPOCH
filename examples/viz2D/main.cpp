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
#include "Helpers.h"

using namespace VFXEpoch::Gluvi;

PanZoom2D cam(-0.1, -0.35, 1.2);
bool loadBin(const char* filename);
void initData();
void display();
void mouse(int button, int state, int x, int y);
void drag(int x, int y);
void timer(int junk);

std::vector<VFXEpoch::Vector2Df> particles;

int main(int argc, char **argv)
{   
   //Setup viewer stuff
   Gluvi::init("2D Particles Visualizer", &argc, argv, 720, 480);
   initData();
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
initData(){
    // TODO: Initialize data;
}

bool
loadBin(const char* filename){
    // TODO: Load binary files
    return true;
}

void
display(){
    // TODO: Draw particles
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


