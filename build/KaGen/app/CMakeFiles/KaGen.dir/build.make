# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.22

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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = "/mnt/d/informatik/semester 12/masterarbeit/code"

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = "/mnt/d/informatik/semester 12/masterarbeit/code/build"

# Include any dependencies generated for this target.
include KaGen/app/CMakeFiles/KaGen.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include KaGen/app/CMakeFiles/KaGen.dir/compiler_depend.make

# Include the progress variables for this target.
include KaGen/app/CMakeFiles/KaGen.dir/progress.make

# Include the compile flags for this target's objects.
include KaGen/app/CMakeFiles/KaGen.dir/flags.make

KaGen/app/CMakeFiles/KaGen.dir/KaGen.cpp.o: KaGen/app/CMakeFiles/KaGen.dir/flags.make
KaGen/app/CMakeFiles/KaGen.dir/KaGen.cpp.o: ../KaGen/app/KaGen.cpp
KaGen/app/CMakeFiles/KaGen.dir/KaGen.cpp.o: KaGen/app/CMakeFiles/KaGen.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/mnt/d/informatik/semester 12/masterarbeit/code/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object KaGen/app/CMakeFiles/KaGen.dir/KaGen.cpp.o"
	cd "/mnt/d/informatik/semester 12/masterarbeit/code/build/KaGen/app" && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT KaGen/app/CMakeFiles/KaGen.dir/KaGen.cpp.o -MF CMakeFiles/KaGen.dir/KaGen.cpp.o.d -o CMakeFiles/KaGen.dir/KaGen.cpp.o -c "/mnt/d/informatik/semester 12/masterarbeit/code/KaGen/app/KaGen.cpp"

KaGen/app/CMakeFiles/KaGen.dir/KaGen.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/KaGen.dir/KaGen.cpp.i"
	cd "/mnt/d/informatik/semester 12/masterarbeit/code/build/KaGen/app" && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/mnt/d/informatik/semester 12/masterarbeit/code/KaGen/app/KaGen.cpp" > CMakeFiles/KaGen.dir/KaGen.cpp.i

KaGen/app/CMakeFiles/KaGen.dir/KaGen.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/KaGen.dir/KaGen.cpp.s"
	cd "/mnt/d/informatik/semester 12/masterarbeit/code/build/KaGen/app" && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/mnt/d/informatik/semester 12/masterarbeit/code/KaGen/app/KaGen.cpp" -o CMakeFiles/KaGen.dir/KaGen.cpp.s

# Object files for target KaGen
KaGen_OBJECTS = \
"CMakeFiles/KaGen.dir/KaGen.cpp.o"

# External object files for target KaGen
KaGen_EXTERNAL_OBJECTS = \
"/mnt/d/informatik/semester 12/masterarbeit/code/build/KaGen/kagen/CMakeFiles/kagen.dir/ckagen.cpp.o" \
"/mnt/d/informatik/semester 12/masterarbeit/code/build/KaGen/kagen/CMakeFiles/kagen.dir/context.cpp.o" \
"/mnt/d/informatik/semester 12/masterarbeit/code/build/KaGen/kagen/CMakeFiles/kagen.dir/facade.cpp.o" \
"/mnt/d/informatik/semester 12/masterarbeit/code/build/KaGen/kagen/CMakeFiles/kagen.dir/generators/barabassi/barabassi.cpp.o" \
"/mnt/d/informatik/semester 12/masterarbeit/code/build/KaGen/kagen/CMakeFiles/kagen.dir/generators/file/file_graph.cpp.o" \
"/mnt/d/informatik/semester 12/masterarbeit/code/build/KaGen/kagen/CMakeFiles/kagen.dir/generators/generator.cpp.o" \
"/mnt/d/informatik/semester 12/masterarbeit/code/build/KaGen/kagen/CMakeFiles/kagen.dir/generators/geometric/delaunay.cpp.o" \
"/mnt/d/informatik/semester 12/masterarbeit/code/build/KaGen/kagen/CMakeFiles/kagen.dir/generators/geometric/delaunay/delaunay_2d.cpp.o" \
"/mnt/d/informatik/semester 12/masterarbeit/code/build/KaGen/kagen/CMakeFiles/kagen.dir/generators/geometric/delaunay/delaunay_3d.cpp.o" \
"/mnt/d/informatik/semester 12/masterarbeit/code/build/KaGen/kagen/CMakeFiles/kagen.dir/generators/geometric/rgg.cpp.o" \
"/mnt/d/informatik/semester 12/masterarbeit/code/build/KaGen/kagen/CMakeFiles/kagen.dir/generators/geometric/rgg/rgg_2d.cpp.o" \
"/mnt/d/informatik/semester 12/masterarbeit/code/build/KaGen/kagen/CMakeFiles/kagen.dir/generators/geometric/rgg/rgg_3d.cpp.o" \
"/mnt/d/informatik/semester 12/masterarbeit/code/build/KaGen/kagen/CMakeFiles/kagen.dir/generators/gnm/gnm_directed.cpp.o" \
"/mnt/d/informatik/semester 12/masterarbeit/code/build/KaGen/kagen/CMakeFiles/kagen.dir/generators/gnm/gnm_undirected.cpp.o" \
"/mnt/d/informatik/semester 12/masterarbeit/code/build/KaGen/kagen/CMakeFiles/kagen.dir/generators/gnp/gnp_directed.cpp.o" \
"/mnt/d/informatik/semester 12/masterarbeit/code/build/KaGen/kagen/CMakeFiles/kagen.dir/generators/gnp/gnp_undirected.cpp.o" \
"/mnt/d/informatik/semester 12/masterarbeit/code/build/KaGen/kagen/CMakeFiles/kagen.dir/generators/graph500_generator.cpp.o" \
"/mnt/d/informatik/semester 12/masterarbeit/code/build/KaGen/kagen/CMakeFiles/kagen.dir/generators/grid/grid_2d.cpp.o" \
"/mnt/d/informatik/semester 12/masterarbeit/code/build/KaGen/kagen/CMakeFiles/kagen.dir/generators/grid/grid_3d.cpp.o" \
"/mnt/d/informatik/semester 12/masterarbeit/code/build/KaGen/kagen/CMakeFiles/kagen.dir/generators/hyperbolic/hyperbolic.cpp.o" \
"/mnt/d/informatik/semester 12/masterarbeit/code/build/KaGen/kagen/CMakeFiles/kagen.dir/generators/image/image_mesh.cpp.o" \
"/mnt/d/informatik/semester 12/masterarbeit/code/build/KaGen/kagen/CMakeFiles/kagen.dir/generators/image/kargb.cpp.o" \
"/mnt/d/informatik/semester 12/masterarbeit/code/build/KaGen/kagen/CMakeFiles/kagen.dir/generators/kronecker/kronecker.cpp.o" \
"/mnt/d/informatik/semester 12/masterarbeit/code/build/KaGen/kagen/CMakeFiles/kagen.dir/generators/path/path_directed.cpp.o" \
"/mnt/d/informatik/semester 12/masterarbeit/code/build/KaGen/kagen/CMakeFiles/kagen.dir/generators/rmat/generators/dSFMT.cpp.o" \
"/mnt/d/informatik/semester 12/masterarbeit/code/build/KaGen/kagen/CMakeFiles/kagen.dir/generators/rmat/memory.cpp.o" \
"/mnt/d/informatik/semester 12/masterarbeit/code/build/KaGen/kagen/CMakeFiles/kagen.dir/generators/rmat/parallel_do.cpp.o" \
"/mnt/d/informatik/semester 12/masterarbeit/code/build/KaGen/kagen/CMakeFiles/kagen.dir/generators/rmat/rmat.cpp.o" \
"/mnt/d/informatik/semester 12/masterarbeit/code/build/KaGen/kagen/CMakeFiles/kagen.dir/io.cpp.o" \
"/mnt/d/informatik/semester 12/masterarbeit/code/build/KaGen/kagen/CMakeFiles/kagen.dir/io/coordinates.cpp.o" \
"/mnt/d/informatik/semester 12/masterarbeit/code/build/KaGen/kagen/CMakeFiles/kagen.dir/io/dot.cpp.o" \
"/mnt/d/informatik/semester 12/masterarbeit/code/build/KaGen/kagen/CMakeFiles/kagen.dir/io/edgelist.cpp.o" \
"/mnt/d/informatik/semester 12/masterarbeit/code/build/KaGen/kagen/CMakeFiles/kagen.dir/io/graph_format.cpp.o" \
"/mnt/d/informatik/semester 12/masterarbeit/code/build/KaGen/kagen/CMakeFiles/kagen.dir/io/hmetis.cpp.o" \
"/mnt/d/informatik/semester 12/masterarbeit/code/build/KaGen/kagen/CMakeFiles/kagen.dir/io/metis.cpp.o" \
"/mnt/d/informatik/semester 12/masterarbeit/code/build/KaGen/kagen/CMakeFiles/kagen.dir/io/parhip.cpp.o" \
"/mnt/d/informatik/semester 12/masterarbeit/code/build/KaGen/kagen/CMakeFiles/kagen.dir/kagen.cpp.o" \
"/mnt/d/informatik/semester 12/masterarbeit/code/build/KaGen/kagen/CMakeFiles/kagen.dir/sampling/hash.cpp.o" \
"/mnt/d/informatik/semester 12/masterarbeit/code/build/KaGen/kagen/CMakeFiles/kagen.dir/sampling/rng/dSFMT.cpp.o" \
"/mnt/d/informatik/semester 12/masterarbeit/code/build/KaGen/kagen/CMakeFiles/kagen.dir/sampling/spooky/spooky.cpp.o" \
"/mnt/d/informatik/semester 12/masterarbeit/code/build/KaGen/kagen/CMakeFiles/kagen.dir/tlx/thread_pool.cpp.o" \
"/mnt/d/informatik/semester 12/masterarbeit/code/build/KaGen/kagen/CMakeFiles/kagen.dir/tools/postprocessor.cpp.o" \
"/mnt/d/informatik/semester 12/masterarbeit/code/build/KaGen/kagen/CMakeFiles/kagen.dir/tools/statistics.cpp.o" \
"/mnt/d/informatik/semester 12/masterarbeit/code/build/KaGen/kagen/CMakeFiles/kagen.dir/tools/validator.cpp.o"

KaGen/app/KaGen: KaGen/app/CMakeFiles/KaGen.dir/KaGen.cpp.o
KaGen/app/KaGen: KaGen/kagen/CMakeFiles/kagen.dir/ckagen.cpp.o
KaGen/app/KaGen: KaGen/kagen/CMakeFiles/kagen.dir/context.cpp.o
KaGen/app/KaGen: KaGen/kagen/CMakeFiles/kagen.dir/facade.cpp.o
KaGen/app/KaGen: KaGen/kagen/CMakeFiles/kagen.dir/generators/barabassi/barabassi.cpp.o
KaGen/app/KaGen: KaGen/kagen/CMakeFiles/kagen.dir/generators/file/file_graph.cpp.o
KaGen/app/KaGen: KaGen/kagen/CMakeFiles/kagen.dir/generators/generator.cpp.o
KaGen/app/KaGen: KaGen/kagen/CMakeFiles/kagen.dir/generators/geometric/delaunay.cpp.o
KaGen/app/KaGen: KaGen/kagen/CMakeFiles/kagen.dir/generators/geometric/delaunay/delaunay_2d.cpp.o
KaGen/app/KaGen: KaGen/kagen/CMakeFiles/kagen.dir/generators/geometric/delaunay/delaunay_3d.cpp.o
KaGen/app/KaGen: KaGen/kagen/CMakeFiles/kagen.dir/generators/geometric/rgg.cpp.o
KaGen/app/KaGen: KaGen/kagen/CMakeFiles/kagen.dir/generators/geometric/rgg/rgg_2d.cpp.o
KaGen/app/KaGen: KaGen/kagen/CMakeFiles/kagen.dir/generators/geometric/rgg/rgg_3d.cpp.o
KaGen/app/KaGen: KaGen/kagen/CMakeFiles/kagen.dir/generators/gnm/gnm_directed.cpp.o
KaGen/app/KaGen: KaGen/kagen/CMakeFiles/kagen.dir/generators/gnm/gnm_undirected.cpp.o
KaGen/app/KaGen: KaGen/kagen/CMakeFiles/kagen.dir/generators/gnp/gnp_directed.cpp.o
KaGen/app/KaGen: KaGen/kagen/CMakeFiles/kagen.dir/generators/gnp/gnp_undirected.cpp.o
KaGen/app/KaGen: KaGen/kagen/CMakeFiles/kagen.dir/generators/graph500_generator.cpp.o
KaGen/app/KaGen: KaGen/kagen/CMakeFiles/kagen.dir/generators/grid/grid_2d.cpp.o
KaGen/app/KaGen: KaGen/kagen/CMakeFiles/kagen.dir/generators/grid/grid_3d.cpp.o
KaGen/app/KaGen: KaGen/kagen/CMakeFiles/kagen.dir/generators/hyperbolic/hyperbolic.cpp.o
KaGen/app/KaGen: KaGen/kagen/CMakeFiles/kagen.dir/generators/image/image_mesh.cpp.o
KaGen/app/KaGen: KaGen/kagen/CMakeFiles/kagen.dir/generators/image/kargb.cpp.o
KaGen/app/KaGen: KaGen/kagen/CMakeFiles/kagen.dir/generators/kronecker/kronecker.cpp.o
KaGen/app/KaGen: KaGen/kagen/CMakeFiles/kagen.dir/generators/path/path_directed.cpp.o
KaGen/app/KaGen: KaGen/kagen/CMakeFiles/kagen.dir/generators/rmat/generators/dSFMT.cpp.o
KaGen/app/KaGen: KaGen/kagen/CMakeFiles/kagen.dir/generators/rmat/memory.cpp.o
KaGen/app/KaGen: KaGen/kagen/CMakeFiles/kagen.dir/generators/rmat/parallel_do.cpp.o
KaGen/app/KaGen: KaGen/kagen/CMakeFiles/kagen.dir/generators/rmat/rmat.cpp.o
KaGen/app/KaGen: KaGen/kagen/CMakeFiles/kagen.dir/io.cpp.o
KaGen/app/KaGen: KaGen/kagen/CMakeFiles/kagen.dir/io/coordinates.cpp.o
KaGen/app/KaGen: KaGen/kagen/CMakeFiles/kagen.dir/io/dot.cpp.o
KaGen/app/KaGen: KaGen/kagen/CMakeFiles/kagen.dir/io/edgelist.cpp.o
KaGen/app/KaGen: KaGen/kagen/CMakeFiles/kagen.dir/io/graph_format.cpp.o
KaGen/app/KaGen: KaGen/kagen/CMakeFiles/kagen.dir/io/hmetis.cpp.o
KaGen/app/KaGen: KaGen/kagen/CMakeFiles/kagen.dir/io/metis.cpp.o
KaGen/app/KaGen: KaGen/kagen/CMakeFiles/kagen.dir/io/parhip.cpp.o
KaGen/app/KaGen: KaGen/kagen/CMakeFiles/kagen.dir/kagen.cpp.o
KaGen/app/KaGen: KaGen/kagen/CMakeFiles/kagen.dir/sampling/hash.cpp.o
KaGen/app/KaGen: KaGen/kagen/CMakeFiles/kagen.dir/sampling/rng/dSFMT.cpp.o
KaGen/app/KaGen: KaGen/kagen/CMakeFiles/kagen.dir/sampling/spooky/spooky.cpp.o
KaGen/app/KaGen: KaGen/kagen/CMakeFiles/kagen.dir/tlx/thread_pool.cpp.o
KaGen/app/KaGen: KaGen/kagen/CMakeFiles/kagen.dir/tools/postprocessor.cpp.o
KaGen/app/KaGen: KaGen/kagen/CMakeFiles/kagen.dir/tools/statistics.cpp.o
KaGen/app/KaGen: KaGen/kagen/CMakeFiles/kagen.dir/tools/validator.cpp.o
KaGen/app/KaGen: KaGen/app/CMakeFiles/KaGen.dir/build.make
KaGen/app/KaGen: /usr/lib/x86_64-linux-gnu/libgmpxx.so
KaGen/app/KaGen: /usr/lib/x86_64-linux-gnu/libmpfr.so
KaGen/app/KaGen: /usr/lib/x86_64-linux-gnu/libgmp.so
KaGen/app/KaGen: KaGen/extlib/xxHash/cmake_unofficial/libxxhash.a
KaGen/app/KaGen: /usr/lib/x86_64-linux-gnu/openmpi/lib/libmpi_cxx.so
KaGen/app/KaGen: /usr/lib/x86_64-linux-gnu/openmpi/lib/libmpi.so
KaGen/app/KaGen: KaGen/app/CMakeFiles/KaGen.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir="/mnt/d/informatik/semester 12/masterarbeit/code/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable KaGen"
	cd "/mnt/d/informatik/semester 12/masterarbeit/code/build/KaGen/app" && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/KaGen.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
KaGen/app/CMakeFiles/KaGen.dir/build: KaGen/app/KaGen
.PHONY : KaGen/app/CMakeFiles/KaGen.dir/build

KaGen/app/CMakeFiles/KaGen.dir/clean:
	cd "/mnt/d/informatik/semester 12/masterarbeit/code/build/KaGen/app" && $(CMAKE_COMMAND) -P CMakeFiles/KaGen.dir/cmake_clean.cmake
.PHONY : KaGen/app/CMakeFiles/KaGen.dir/clean

KaGen/app/CMakeFiles/KaGen.dir/depend:
	cd "/mnt/d/informatik/semester 12/masterarbeit/code/build" && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" "/mnt/d/informatik/semester 12/masterarbeit/code" "/mnt/d/informatik/semester 12/masterarbeit/code/KaGen/app" "/mnt/d/informatik/semester 12/masterarbeit/code/build" "/mnt/d/informatik/semester 12/masterarbeit/code/build/KaGen/app" "/mnt/d/informatik/semester 12/masterarbeit/code/build/KaGen/app/CMakeFiles/KaGen.dir/DependInfo.cmake" --color=$(COLOR)
.PHONY : KaGen/app/CMakeFiles/KaGen.dir/depend

