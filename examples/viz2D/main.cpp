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

// OpenEXR libraries headers
#include <OpenEXR/ImfRgbaFile.h>
#include <OpenEXR/ImfStringAttribute.h>
#include <OpenEXR/ImfMatrixAttribute.h>
#include <OpenEXR/ImfArray.h>
#include <OpenEXR/ImfNamespace.h>

#include "VisualizerHelpers.h"
#include "VisualizerHelpers_OpenGL.h"
#include "Helpers.h"

using namespace VFXEpoch::Gluvi;
using namespace VFXEpoch::OpenGL_Utility;
namespace IMF = OPENEXR_IMF_NAMESPACE;
using namespace IMF;
using namespace IMATH_NAMESPACE;

PanZoom2D cam(-2.5f, 0.0f, 5.0f, false);
bool load_bin(const char* filename);
void write_exrs(const char fileName[], const Rgba *pixels, int width, int height);
void convert_to_exr_rgba(GLubyte* in_pixels, Array2D<Imf::Rgba>& out_pixels, int width, int height);
void init_data();
void display();
void mouse(int button, int state, int x, int y);
void drag(int x, int y);
void timer(int junk);

std::vector<VFXEpoch::Vector2Df> particles;
unsigned int total_frames = 300;
unsigned int frame_counter = 0;
unsigned int width = 960;
unsigned int height = 280;
bool is_write_to_disk = false;

int main(int argc, char **argv)
{   
   //Setup viewer stuff
   Gluvi::init("2D Particles Visualizer", &argc, argv, width, height);
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
    // TODO: Initialize data
    frame_counter = 0;
    char filename[256];
    bool result;
    sprintf(filename, "../../outputs/sims/Particle_data%04d.bin", frame_counter);
    std::string str_filename(filename);
    result = load_bin(filename);
    if(!result){
#ifdef __linux__
        cout << "\033[1;31mERROR:\033[0m" << "Loading simulation file faild!" << endl;
        cout << "File " << "\033[1;33m" << str_filename << "\033[0m" << " does not exist!" << endl;
#endif
        exit(-1);
    }
}

bool
load_bin(const char* filename){
    // Debug: In progress for visualizing particles
    FILE* f;
    f = fopen(filename, "rb");
    if(!f){
        return false;
    }

    int num_of_particles = 0;
    fread(&num_of_particles, sizeof(num_of_particles), 1, f);
    float *data = new float[num_of_particles * 4]; if(!data) return false;
    fread(data, sizeof(data), 4 * num_of_particles, f);
    fclose(f);
    
    for(int i = 0; i != num_of_particles; i++){
        VFXEpoch::Vector2Df v(data[i * 4 + 0], data[i * 4 + 1]);
        particles.push_back(v);
    }
    
    delete [] data;
    return true;
}

void
write_exrs(const char fileName[], const Rgba *pixels, int width, int height){
    RgbaOutputFile file (fileName, width, height, WRITE_RGBA);
    file.setFrameBuffer (pixels, 1, width);
    file.writePixels (height);
}

void 
convert_to_exr_rgba(GLubyte* in_pixels, Array2D<Imf::Rgba>& out_pixels, int width, int height){
    assert(NULL != in_pixels && width >= 0 && height >= 0);
    for(int i = 0; i != height; i++){
        for(int j = 0; j != width; j++){
            int index = 4 * ((height - i - 1) * width + j);
            out_pixels[i][j].r = in_pixels[index];
            out_pixels[i][j].g = in_pixels[index + 1];
            out_pixels[i][j].b = in_pixels[index + 2];
            out_pixels[i][j].a = in_pixels[index + 3];
        }
    }
}

void
display(){
    glPointSize(1);
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
    particles.clear();
    char filename[256];
    bool result;
    frame_counter++;
    sprintf(filename, "../../outputs/sims/Particle_data%04d.bin", frame_counter);
    result = load_bin(filename);
    if(!result){
#ifdef __linux__
        string str_filename(filename);
        cout << "\033[1;31mERROR:\033[0m" << "Loading simulation file faild!" << endl;
        cout << "File " << "\033[1;33m" << str_filename << "\033[0m" << " may not exist!" << endl;
#endif
        exit(0);
    }

    if(is_write_to_disk){
        sprintf(filename, "../../outputs/anims/frame_%04d.exr", frame_counter);
        std::string frame_name(filename);
        cout << "Writing frame " << frame_counter << " to " << frame_name << endl;
        GLubyte* read_pixels = new GLubyte[width * height * 4];
        Array2D<Imf::Rgba> write_pixels(height, width);
        if(!read_pixels){
            cout << "\033[1;31mUnable to write images into disk!\033[0m" << endl;
            exit(-1);
        }
        // --> Get pixels from OpenGL window
        glReadPixels(0, 0, width, height, GL_RGBA, GL_UNSIGNED_BYTE, read_pixels);

        // --> Convert pixels from Glubyte to OpenEXR::Rgba
        convert_to_exr_rgba(read_pixels, write_pixels, width, height);

        // --> Write the frame to a file
        write_exrs(filename, &write_pixels[0][0], width, height);

        delete [] read_pixels;
    }
    glutPostRedisplay();
    glutTimerFunc(30, timer, 0);
}


