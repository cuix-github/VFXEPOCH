# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list

# Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = /usr/bin/ccmake

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/faq/Documents/dev/VFXEPOCH

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/faq/Documents/dev/VFXEPOCH/build

# Include any dependencies generated for this target.
include examples/CMakeFiles/examiner.dir/depend.make

# Include the progress variables for this target.
include examples/CMakeFiles/examiner.dir/progress.make

# Include the compile flags for this target's objects.
include examples/CMakeFiles/examiner.dir/flags.make

examples/CMakeFiles/examiner.dir/Common/Helpers.cpp.o: examples/CMakeFiles/examiner.dir/flags.make
examples/CMakeFiles/examiner.dir/Common/Helpers.cpp.o: ../examples/Common/Helpers.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/faq/Documents/dev/VFXEPOCH/build/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object examples/CMakeFiles/examiner.dir/Common/Helpers.cpp.o"
	cd /home/faq/Documents/dev/VFXEPOCH/build/examples && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/examiner.dir/Common/Helpers.cpp.o -c /home/faq/Documents/dev/VFXEPOCH/examples/Common/Helpers.cpp

examples/CMakeFiles/examiner.dir/Common/Helpers.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/examiner.dir/Common/Helpers.cpp.i"
	cd /home/faq/Documents/dev/VFXEPOCH/build/examples && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/faq/Documents/dev/VFXEPOCH/examples/Common/Helpers.cpp > CMakeFiles/examiner.dir/Common/Helpers.cpp.i

examples/CMakeFiles/examiner.dir/Common/Helpers.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/examiner.dir/Common/Helpers.cpp.s"
	cd /home/faq/Documents/dev/VFXEPOCH/build/examples && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/faq/Documents/dev/VFXEPOCH/examples/Common/Helpers.cpp -o CMakeFiles/examiner.dir/Common/Helpers.cpp.s

examples/CMakeFiles/examiner.dir/Common/Helpers.cpp.o.requires:
.PHONY : examples/CMakeFiles/examiner.dir/Common/Helpers.cpp.o.requires

examples/CMakeFiles/examiner.dir/Common/Helpers.cpp.o.provides: examples/CMakeFiles/examiner.dir/Common/Helpers.cpp.o.requires
	$(MAKE) -f examples/CMakeFiles/examiner.dir/build.make examples/CMakeFiles/examiner.dir/Common/Helpers.cpp.o.provides.build
.PHONY : examples/CMakeFiles/examiner.dir/Common/Helpers.cpp.o.provides

examples/CMakeFiles/examiner.dir/Common/Helpers.cpp.o.provides.build: examples/CMakeFiles/examiner.dir/Common/Helpers.cpp.o

examples/CMakeFiles/examiner.dir/examiner/main.cpp.o: examples/CMakeFiles/examiner.dir/flags.make
examples/CMakeFiles/examiner.dir/examiner/main.cpp.o: ../examples/examiner/main.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/faq/Documents/dev/VFXEPOCH/build/CMakeFiles $(CMAKE_PROGRESS_2)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object examples/CMakeFiles/examiner.dir/examiner/main.cpp.o"
	cd /home/faq/Documents/dev/VFXEPOCH/build/examples && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/examiner.dir/examiner/main.cpp.o -c /home/faq/Documents/dev/VFXEPOCH/examples/examiner/main.cpp

examples/CMakeFiles/examiner.dir/examiner/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/examiner.dir/examiner/main.cpp.i"
	cd /home/faq/Documents/dev/VFXEPOCH/build/examples && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/faq/Documents/dev/VFXEPOCH/examples/examiner/main.cpp > CMakeFiles/examiner.dir/examiner/main.cpp.i

examples/CMakeFiles/examiner.dir/examiner/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/examiner.dir/examiner/main.cpp.s"
	cd /home/faq/Documents/dev/VFXEPOCH/build/examples && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/faq/Documents/dev/VFXEPOCH/examples/examiner/main.cpp -o CMakeFiles/examiner.dir/examiner/main.cpp.s

examples/CMakeFiles/examiner.dir/examiner/main.cpp.o.requires:
.PHONY : examples/CMakeFiles/examiner.dir/examiner/main.cpp.o.requires

examples/CMakeFiles/examiner.dir/examiner/main.cpp.o.provides: examples/CMakeFiles/examiner.dir/examiner/main.cpp.o.requires
	$(MAKE) -f examples/CMakeFiles/examiner.dir/build.make examples/CMakeFiles/examiner.dir/examiner/main.cpp.o.provides.build
.PHONY : examples/CMakeFiles/examiner.dir/examiner/main.cpp.o.provides

examples/CMakeFiles/examiner.dir/examiner/main.cpp.o.provides.build: examples/CMakeFiles/examiner.dir/examiner/main.cpp.o

# Object files for target examiner
examiner_OBJECTS = \
"CMakeFiles/examiner.dir/Common/Helpers.cpp.o" \
"CMakeFiles/examiner.dir/examiner/main.cpp.o"

# External object files for target examiner
examiner_EXTERNAL_OBJECTS =

examples/examiner: examples/CMakeFiles/examiner.dir/Common/Helpers.cpp.o
examples/examiner: examples/CMakeFiles/examiner.dir/examiner/main.cpp.o
examples/examiner: examples/CMakeFiles/examiner.dir/build.make
examples/examiner: source/libVFXEpoch.a
examples/examiner: /usr/lib/x86_64-linux-gnu/libGLU.so
examples/examiner: /usr/lib/x86_64-linux-gnu/libGL.so
examples/examiner: /usr/lib/x86_64-linux-gnu/libSM.so
examples/examiner: /usr/lib/x86_64-linux-gnu/libICE.so
examples/examiner: /usr/lib/x86_64-linux-gnu/libX11.so
examples/examiner: /usr/lib/x86_64-linux-gnu/libXext.so
examples/examiner: /usr/lib/x86_64-linux-gnu/libglut.so
examples/examiner: /usr/lib/x86_64-linux-gnu/libXmu.so
examples/examiner: /usr/lib/x86_64-linux-gnu/libXi.so
examples/examiner: examples/CMakeFiles/examiner.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable examiner"
	cd /home/faq/Documents/dev/VFXEPOCH/build/examples && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/examiner.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
examples/CMakeFiles/examiner.dir/build: examples/examiner
.PHONY : examples/CMakeFiles/examiner.dir/build

examples/CMakeFiles/examiner.dir/requires: examples/CMakeFiles/examiner.dir/Common/Helpers.cpp.o.requires
examples/CMakeFiles/examiner.dir/requires: examples/CMakeFiles/examiner.dir/examiner/main.cpp.o.requires
.PHONY : examples/CMakeFiles/examiner.dir/requires

examples/CMakeFiles/examiner.dir/clean:
	cd /home/faq/Documents/dev/VFXEPOCH/build/examples && $(CMAKE_COMMAND) -P CMakeFiles/examiner.dir/cmake_clean.cmake
.PHONY : examples/CMakeFiles/examiner.dir/clean

examples/CMakeFiles/examiner.dir/depend:
	cd /home/faq/Documents/dev/VFXEPOCH/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/faq/Documents/dev/VFXEPOCH /home/faq/Documents/dev/VFXEPOCH/examples /home/faq/Documents/dev/VFXEPOCH/build /home/faq/Documents/dev/VFXEPOCH/build/examples /home/faq/Documents/dev/VFXEPOCH/build/examples/CMakeFiles/examiner.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : examples/CMakeFiles/examiner.dir/depend
