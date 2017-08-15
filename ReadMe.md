### **STATUS**
| **`Linux - ubuntu`** |
|:----------------------:|
|![Build Status](https://travis-ci.org/Shakebones/VFXEPOCH.svg?branch=master)|

### **Tested Platform**
1. gcc-4.8 & ubuntu 14.04 failed
2. gcc-5.4 & ubuntu 14.04 passed
3. gcc-6.2 & ubuntu 14.04 passed
4. gcc-default & ubuntu 16.04 passed
5. gcc-5.4 & linux-mint 18.2 passed

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

### **Third Party Library Used**
**1. OpenEXR (Including ILMBase)**
&nbsp;&nbsp;&nbsp;&nbsp;Checkout [here](http://www.openexr.com/downloads.html/) for downloading OpenEXR & ILMBase. For "How to correctly install" please refer to the link [here](http://www.openexr.com/documentation.html/)
**2. Alembic**
&nbsp;&nbsp;&nbsp;&nbsp;Checkout [here](https://github.com/alembic/alembic) for downloading and [here](http://docs.alembic.io/#build-alembic) for how to build it from scratch. Don't forget to install into your local path:
```sh
$ sudo make install
```
**3. Eigen**
&nbsp;&nbsp;&nbsp;&nbsp;Checkout [here](http://eigen.tuxfamily.org/index.php?title=Main_Page/) for downloading the library. Note that this library only has header files so there is no need to build.
&nbsp;&nbsp;&nbsp;&nbsp;Additionally, we grab Robert Bridson's Pre-Conditioned Conjugate Gradient (PCG) solver during fluid simulation for pressure solve step. The code has been wrapped up into a specific folder called "/source/util/PCGSolver".


### **How to compile**
1. Except for VFXEPOCH libraries, examples requires only OpenGL GLUT to be installed. Following the link [here](http://kiwwito.com/installing-opengl-glut-libraries-in-ubuntu/) to test your OpenGL GLUT.
2. Go to the root directory of VFXEPOCH folder. Open the CMakeLists.txt, find the line:
**option(VFXEPOCH_EXAMPLES "Turn ON to build example projects" OFF)** to makre sure the option has been turned off for the very first time compiling the VFXEPOCH library.
3. You should find a scripts files called "cleanup.sh" & "setup.sh". Firstly, try "cleanup.sh" to make it naked and ready to build. Then use the bash "setup.sh" to setup and build the projects. Go to the root directory of VFXEPOCH project folder, and following are the scripts to execute "cleanup.sh" & "setup.sh".
```sh
$ ./cleanup.sh && ./setup.sh
```
When it's done, you should see an examples folder under the build directory. Go into it and have fun!

### **Issues for current version:**
