### **STATUS**
| **`Linux - ubuntu`** |
|:----------------------:|
|![Build Status](https://travis-ci.org/Shakebones/VFXEPOCH.svg?branch=master)|

### **Note:**
1. gcc-4.8 & ubuntu 14.04 failed
2. gcc-5.4 & ubuntu 14.04 passed
3. gcc-6.2 & ubuntu 14.04 passed
4. gcc-default & ubuntu 16.04 passed

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
2. Go to the root directory of VFXEPOCH folder. Open the CMakeLists.txt, find the line:
**option(VFXEPOCH_EXAMPLES "Turn ON to build example projects" OFF)** to makre sure the option has been turned off for the very first time compiling the VFXEPOCH library.
3. You should find a scripts files called "cleanup.sh" & "setup.sh". Firstly, try "cleanup.sh" to make it naked and ready to build. Then use the bash "setup.sh" to setup and build the projects. Go to the root directory of VFXEPOCH project folder, and following are the scripts to execute "cleanup.sh" & "setup.sh".
```sh
$ ./cleanup.sh && ./setup.sh
```
When it's done, you should see an examples folder under the build directory. Go into it and have fun!

### **Issues for current version:**
