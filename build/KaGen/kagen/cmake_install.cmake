# Install script for directory: /mnt/d/informatik/semester 12/masterarbeit/code/KaGen/kagen

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/usr/bin/objdump")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/objects-Release/kagen" TYPE FILE FILES
    "ckagen.cpp.o"
    "context.cpp.o"
    "facade.cpp.o"
    "generators/barabassi/barabassi.cpp.o"
    "generators/file/file_graph.cpp.o"
    "generators/generator.cpp.o"
    "generators/geometric/delaunay.cpp.o"
    "generators/geometric/delaunay/delaunay_2d.cpp.o"
    "generators/geometric/delaunay/delaunay_3d.cpp.o"
    "generators/geometric/rgg.cpp.o"
    "generators/geometric/rgg/rgg_2d.cpp.o"
    "generators/geometric/rgg/rgg_3d.cpp.o"
    "generators/gnm/gnm_directed.cpp.o"
    "generators/gnm/gnm_undirected.cpp.o"
    "generators/gnp/gnp_directed.cpp.o"
    "generators/gnp/gnp_undirected.cpp.o"
    "generators/graph500_generator.cpp.o"
    "generators/grid/grid_2d.cpp.o"
    "generators/grid/grid_3d.cpp.o"
    "generators/hyperbolic/hyperbolic.cpp.o"
    "generators/image/image_mesh.cpp.o"
    "generators/image/kargb.cpp.o"
    "generators/kronecker/kronecker.cpp.o"
    "generators/path/path_directed.cpp.o"
    "generators/rmat/generators/dSFMT.cpp.o"
    "generators/rmat/memory.cpp.o"
    "generators/rmat/parallel_do.cpp.o"
    "generators/rmat/rmat.cpp.o"
    "io.cpp.o"
    "io/coordinates.cpp.o"
    "io/dot.cpp.o"
    "io/edgelist.cpp.o"
    "io/graph_format.cpp.o"
    "io/hmetis.cpp.o"
    "io/metis.cpp.o"
    "io/parhip.cpp.o"
    "kagen.cpp.o"
    "sampling/hash.cpp.o"
    "sampling/rng/dSFMT.cpp.o"
    "sampling/spooky/spooky.cpp.o"
    "tlx/thread_pool.cpp.o"
    "tools/postprocessor.cpp.o"
    "tools/statistics.cpp.o"
    "tools/validator.cpp.o"
    FILES_FROM_DIR "/mnt/d/informatik/semester 12/masterarbeit/code/build/KaGen/kagen/CMakeFiles/kagen.dir/")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include" TYPE FILE FILES "/mnt/d/informatik/semester 12/masterarbeit/code/KaGen/kagen/kagen.h")
endif()

