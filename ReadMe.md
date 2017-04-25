| **`Linux - ubuntu`** |
|----------------------|
|[![Build Status](https://travis-ci.org/Shakebones/VFXEPOCH.svg?branch=master)]|

### **Note:**
1. gcc-4.8 & ubuntu 14.04 failed
2. gcc-5.4 & ubuntu 14.04 passed
3. gcc-6.2 & ubuntu 14.04 passed

**Commands for checking and switching your gcc version:**
```sh
$ ls -l $(which gcc)
```

*If you have multiple gcc version, then use following commands to switch to
the version that we have successfully tested:*

```sh
$ sudo ln -sf (gcc version) /usr/bin/gcc
```
For example:
```sh
$ sudo ln -sf gcc-5 /usr/bin/gcc
```

### **How to compile:**
1. Except for VFXEPOCH libraries, examples requires only OpenGL GLUT to be installed. Following the link [here](http://kiwwito.com/installing-opengl-glut-libraries-in-ubuntu/) to test your OpenGL GLUT.
2. Go to the root directory of VFXEPOCH. Open the CMakeLists.txt, find the line:
**option(VFXEPOCH_EXAMPLES "Turn ON to build example projects" OFF)** to makre sure the option has been turned off for the very first time compiling the VFXEPOCH library.
3. Use following commands to create a build direcotry and generate makefiles.
```sh
$ mkdir build && cd build && cmake ..
```
4. Then we get it compiled and install to the default public headers and library directory:
```sh
$ make && make install
```
This command will copy all the headers in the source and the generated library file into the **include** and **vfxepoch_libs** folder.
5. Once the library successfully gets compiled and installed ("make install"). Then go back to the root directory:
```sh
$ cd ../../
```
to turn on the option in the CMakeLists.txt to build examples:
```sh
option(VFXEPOCH_EXAMPLES "Turn ON to build example projects" ON)
```
6. Use the following commands collection to delete build folder and create a new one and generate makefiles:
```sh
$ rm -rf build && mkdir build && cd build && cmake ..
```
7. Build it with examples:
```sh
make
```
When it's done, you should see an examples folder under the build directory. Go into it and have fun!

### **Issues for current version:**
**1.** Cannot support non-square dimension /fluids/euler/SIM_Gas
