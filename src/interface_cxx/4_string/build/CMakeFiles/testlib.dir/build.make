# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.18

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Disable VCS-based implicit rules.
% : %,v


# Disable VCS-based implicit rules.
% : RCS/%


# Disable VCS-based implicit rules.
% : RCS/%,v


# Disable VCS-based implicit rules.
% : SCCS/s.%


# Disable VCS-based implicit rules.
% : s.%


.SUFFIXES: .hpux_make_needs_suffix_list


# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/local/Cellar/cmake/3.18.4/bin/cmake

# The command to remove a file.
RM = /usr/local/Cellar/cmake/3.18.4/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/daniel/.julia/dev/SciAlgs/src/interface_cxx/4_string

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/daniel/.julia/dev/SciAlgs/src/interface_cxx/4_string/build

# Include any dependencies generated for this target.
include CMakeFiles/testlib.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/testlib.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/testlib.dir/flags.make

CMakeFiles/testlib.dir/testlib.cpp.o: CMakeFiles/testlib.dir/flags.make
CMakeFiles/testlib.dir/testlib.cpp.o: ../testlib.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/daniel/.julia/dev/SciAlgs/src/interface_cxx/4_string/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/testlib.dir/testlib.cpp.o"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/testlib.dir/testlib.cpp.o -c /Users/daniel/.julia/dev/SciAlgs/src/interface_cxx/4_string/testlib.cpp

CMakeFiles/testlib.dir/testlib.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/testlib.dir/testlib.cpp.i"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/daniel/.julia/dev/SciAlgs/src/interface_cxx/4_string/testlib.cpp > CMakeFiles/testlib.dir/testlib.cpp.i

CMakeFiles/testlib.dir/testlib.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/testlib.dir/testlib.cpp.s"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/daniel/.julia/dev/SciAlgs/src/interface_cxx/4_string/testlib.cpp -o CMakeFiles/testlib.dir/testlib.cpp.s

# Object files for target testlib
testlib_OBJECTS = \
"CMakeFiles/testlib.dir/testlib.cpp.o"

# External object files for target testlib
testlib_EXTERNAL_OBJECTS =

lib/libtestlib.dylib: CMakeFiles/testlib.dir/testlib.cpp.o
lib/libtestlib.dylib: CMakeFiles/testlib.dir/build.make
lib/libtestlib.dylib: /Users/daniel/.julia/artifacts/6017255205dc4fbf4d962903a855a0c631f092dc/lib/libcxxwrap_julia_stl.dylib
lib/libtestlib.dylib: /Users/daniel/.julia/artifacts/6017255205dc4fbf4d962903a855a0c631f092dc/lib/libcxxwrap_julia.0.8.0.dylib
lib/libtestlib.dylib: /Applications/Julia-1.5.app/Contents/Resources/julia/lib/libjulia.1.5.dylib
lib/libtestlib.dylib: CMakeFiles/testlib.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/daniel/.julia/dev/SciAlgs/src/interface_cxx/4_string/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX shared library lib/libtestlib.dylib"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/testlib.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/testlib.dir/build: lib/libtestlib.dylib

.PHONY : CMakeFiles/testlib.dir/build

CMakeFiles/testlib.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/testlib.dir/cmake_clean.cmake
.PHONY : CMakeFiles/testlib.dir/clean

CMakeFiles/testlib.dir/depend:
	cd /Users/daniel/.julia/dev/SciAlgs/src/interface_cxx/4_string/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/daniel/.julia/dev/SciAlgs/src/interface_cxx/4_string /Users/daniel/.julia/dev/SciAlgs/src/interface_cxx/4_string /Users/daniel/.julia/dev/SciAlgs/src/interface_cxx/4_string/build /Users/daniel/.julia/dev/SciAlgs/src/interface_cxx/4_string/build /Users/daniel/.julia/dev/SciAlgs/src/interface_cxx/4_string/build/CMakeFiles/testlib.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/testlib.dir/depend

