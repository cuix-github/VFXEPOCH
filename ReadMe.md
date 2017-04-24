### **Note:**
1. gcc-4.8 & ubuntu 14.04 failed
2. gcc-5.4 & ubuntu 14.04 pass
3. gcc-6.2 & ubuntu 14.04 pass

***Commands for checking and switching your gcc version:***
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
**1.** Except for VFXEPOCH libraries, examples requires only OpenGL GLUT to be installed. Following the link [here](http://kiwwito.com/installing-opengl-glut-libraries-in-ubuntu/) to test your OpenGL GLUT.
**2.** Go to the root directory of VFXEPOCH. Open the CMakeLists.txt, find the line option(VFXEPOCH_EXAMPLES "Turn ON to build example projects" OFF) to makre sure the option has been turned off for the very first time compiling the VFXEPOCH library

### **Issues for current version:**
**1.** Cannot support non-square dimension /fluids/euler/SIM_Gas
