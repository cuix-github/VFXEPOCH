### **STATUS**
| **`Linux - ubuntu`** |
|:----------------------:|
|![Build Status](https://travis-ci.org/Shakebones/VFXEPOCH.svg?branch=master)|

### **Tested Platform**
1. gcc/g++-4.8 & ubuntu 14.04 failed
2. gcc/g++-5.4 & ubuntu 14.04 passed
3. gcc/g++-6.2 & ubuntu 14.04 passed
4. gcc/g++-default & ubuntu 16.04 passed
5. gcc/g++-5.4 & linux-mint 18.2 passed
6. gcc/g++-6.3 & linux-mint 18.2 passed

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
Also, check if your g++ is installed as well and switch to the version:
```sh
$ sudo ln -sf (g++ version) /usr/bin/g++
```
Example:
```sh
$ sudo ln -sf g++-5 /usr/bin/g++
```

### **Third Party Library Used**
Gnerally, we include the pre-compiled **Linux** libraries for users, but if they do now work well on your side, we
encourage you compiling one for yourself. They all have clear documents to guide you get the work done.

* **OpenEXR (Including ILMBase)**
<br />&nbsp;&nbsp;&nbsp;&nbsp;Checkout [here](http://www.openexr.com/downloads.html) for downloading OpenEXR & ILMBase. For "How to correctly install" please refer to the link [here](http://www.openexr.com/documentation.html)
We start with building and installing ILMBase first:
	* Extract the ilmbase-2.2.0.tar.gz under VFXEPOCH/external_libs/ILM/ by following commands:
	
	```sh
	$ tar -xvf ilmbase-2.2.0.tar.gz
	```
	to current directory. Then you should get an folder with the same name: ilmbase-2.2.0.
	Now let's dive into to the extracted folder and make a directory called "build", the hierarchy should look like:
	<b>ilmbase-2.2.0/build </b>. Dive into the "build" folder and type <b> cmake .. </b>. Check the options and installing target path as you wish. Then type <b>make -j 8</b>. The option "-j 8" is for 8 threads accelerating compiling process. If you have 12 threads 22 threads (Intel xeon E5 maybe?), change it to "-j 12" or "-j 22". Finally type <b>make install</b> to copy files into the specified location.
	&nbsp;&nbsp;&nbsp;&nbsp;Next step is to install OpenEXR. Please note that, installing OpenEXR would rely on <b>ilmbase</b> library from above we addressed, so don't skip it.
	&nbsp;&nbsp;&nbsp;&nbsp;Here if you simply follow the steps like we did for <b>ilmbase</b>, probably you will get an error said it cannot find <b>half.h</b>. Just for simplifying the step, let's put two more lines in the CMakeLists.txt of OpenEXR (openexr-2.2.0, of course extract like above first!). Using any text editor you like to open CMakeLists.txt. Then add the path of ilmbase header files and libs, for my case, it's:
	```sh
	INCLUDE_DIRECTORIES (
			  ${CMAKE_CURRENT_BINARY_DIR}/config
			  /usr/local/include/OpenEXR/
			  IlmImf
			  IlmImfUtil
			  exrmaketiled
			  exrenvmap
			  exrmakepreview
			  exrmultiview
			  IlmImfFuzzTest
	)

	LINK_DIRECTORIES ( ${ILMBASE_PACKAGE_PREFIX}/lib /usr/local/lib/ )
	```
	We do this because we didn't change any options for the installing of <b>ilmbase</b>, if you are aware of this already, you could definitely fix it before just simply specifying the options for <b>ilmbase</b> CMakeLists.txt. Anyway, it should not be an annoy step, try more goolge!
* **Alembic**
<br />&nbsp;&nbsp;&nbsp;&nbsp;Checkout [here](https://github.com/alembic/alembic) for downloading and [here](http://docs.alembic.io/#build-alembic) for how to build it from scratch. Don't forget to install into your local path:
```sh
$ sudo make install
```
* **Eigen**
<br />&nbsp;&nbsp;&nbsp;&nbsp;Checkout [here](http://eigen.tuxfamily.org/index.php?title=Main_Page) for downloading the library. Note that this library only has header files so there is no need to build.

&nbsp;&nbsp;&nbsp;&nbsp;Additionally, we grab Robert Bridson's Pre-Conditioned Conjugate Gradient (PCG) solver during fluid simulation for pressure solve step. The code has been wrapped up into a specific folder called "/source/util/PCGSolver".

### **How to compile**
1. Except for VFXEPOCH libraries, examples requires only OpenGL GLUT to be installed. Following the link [here](http://kiwwito.com/installing-opengl-glut-libraries-in-ubuntu/) to test your OpenGL GLUT.
2. Go to the root directory of VFXEPOCH folder. Open the CMakeLists.txt, find the line:
**option(VFXEPOCH_EXAMPLES "Turn ON to build example projects" OFF)** to check and  to make sure the option has been turned off for the very first time compiling the VFXEPOCH library.
3. You should find a script files called "cleanup.sh" & "setup.sh". Firstly, try "cleanup.sh" to make it naked and ready to build. Then use the bash "setup.sh" to setup and build the projects. Go to the root directory of VFXEPOCH project folder, and following are the scripts to execute "cleanup.sh" & "setup.sh".
```sh
$ ./cleanup.sh && ./setup.sh
```
When it's done, you should see an examples folder under the build directory. Go into it and have fun!
