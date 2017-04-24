**Note:**
1. gcc-4.8 & ubuntu 14.04 failed
2. gcc-5.4 & ubuntu 14.04 pass
3. gcc-6.2 & ubuntu 14.04 pass

**Commands for checking your gcc/g++ version:**
ls -l $(which g++)
ls -l $(which gcc)

If you have multiple gcc/g++ version, then use following commands to switch to
the version that we have successfully tested:

sudo ln -sf (gcc version) /usr/bin/gcc
e.g: sudo ln -sf gcc-5 /usr/bin/gcc

sudo ln -sf (g++ version) /usr/bin/g++
e.g: sudo ln -sf g++-5 /usr/bin/g++

**How to compile:**
1.Go to the root directory of VFXEPOCH. Open the CMakeLists.txt, find the line option(VFXEPOCH_EXAMPLES "Turn ON to build example projects" OFF) to makre sure the option has been turned off for the very first time compiling the VFXEPOCH library

**2. // TODO: Descriptions**

**Issues for current version:**
1. Cannot support non-square dimension /fluids/euler/SIM_Gas
