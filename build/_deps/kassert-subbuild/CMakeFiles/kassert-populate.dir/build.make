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
CMAKE_SOURCE_DIR = "/mnt/d/informatik/semester 12/masterarbeit/code/build/_deps/kassert-subbuild"

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = "/mnt/d/informatik/semester 12/masterarbeit/code/build/_deps/kassert-subbuild"

# Utility rule file for kassert-populate.

# Include any custom commands dependencies for this target.
include CMakeFiles/kassert-populate.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/kassert-populate.dir/progress.make

CMakeFiles/kassert-populate: CMakeFiles/kassert-populate-complete

CMakeFiles/kassert-populate-complete: kassert-populate-prefix/src/kassert-populate-stamp/kassert-populate-install
CMakeFiles/kassert-populate-complete: kassert-populate-prefix/src/kassert-populate-stamp/kassert-populate-mkdir
CMakeFiles/kassert-populate-complete: kassert-populate-prefix/src/kassert-populate-stamp/kassert-populate-download
CMakeFiles/kassert-populate-complete: kassert-populate-prefix/src/kassert-populate-stamp/kassert-populate-update
CMakeFiles/kassert-populate-complete: kassert-populate-prefix/src/kassert-populate-stamp/kassert-populate-patch
CMakeFiles/kassert-populate-complete: kassert-populate-prefix/src/kassert-populate-stamp/kassert-populate-configure
CMakeFiles/kassert-populate-complete: kassert-populate-prefix/src/kassert-populate-stamp/kassert-populate-build
CMakeFiles/kassert-populate-complete: kassert-populate-prefix/src/kassert-populate-stamp/kassert-populate-install
CMakeFiles/kassert-populate-complete: kassert-populate-prefix/src/kassert-populate-stamp/kassert-populate-test
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir="/mnt/d/informatik/semester 12/masterarbeit/code/build/_deps/kassert-subbuild/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_1) "Completed 'kassert-populate'"
	/usr/bin/cmake -E make_directory "/mnt/d/informatik/semester 12/masterarbeit/code/build/_deps/kassert-subbuild/CMakeFiles"
	/usr/bin/cmake -E touch "/mnt/d/informatik/semester 12/masterarbeit/code/build/_deps/kassert-subbuild/CMakeFiles/kassert-populate-complete"
	/usr/bin/cmake -E touch "/mnt/d/informatik/semester 12/masterarbeit/code/build/_deps/kassert-subbuild/kassert-populate-prefix/src/kassert-populate-stamp/kassert-populate-done"

kassert-populate-prefix/src/kassert-populate-stamp/kassert-populate-update:
.PHONY : kassert-populate-prefix/src/kassert-populate-stamp/kassert-populate-update

kassert-populate-prefix/src/kassert-populate-stamp/kassert-populate-build: kassert-populate-prefix/src/kassert-populate-stamp/kassert-populate-configure
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir="/mnt/d/informatik/semester 12/masterarbeit/code/build/_deps/kassert-subbuild/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_2) "No build step for 'kassert-populate'"
	cd "/mnt/d/informatik/semester 12/masterarbeit/code/build/_deps/kassert-build" && /usr/bin/cmake -E echo_append
	cd "/mnt/d/informatik/semester 12/masterarbeit/code/build/_deps/kassert-build" && /usr/bin/cmake -E touch "/mnt/d/informatik/semester 12/masterarbeit/code/build/_deps/kassert-subbuild/kassert-populate-prefix/src/kassert-populate-stamp/kassert-populate-build"

kassert-populate-prefix/src/kassert-populate-stamp/kassert-populate-configure: kassert-populate-prefix/tmp/kassert-populate-cfgcmd.txt
kassert-populate-prefix/src/kassert-populate-stamp/kassert-populate-configure: kassert-populate-prefix/src/kassert-populate-stamp/kassert-populate-patch
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir="/mnt/d/informatik/semester 12/masterarbeit/code/build/_deps/kassert-subbuild/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_3) "No configure step for 'kassert-populate'"
	cd "/mnt/d/informatik/semester 12/masterarbeit/code/build/_deps/kassert-build" && /usr/bin/cmake -E echo_append
	cd "/mnt/d/informatik/semester 12/masterarbeit/code/build/_deps/kassert-build" && /usr/bin/cmake -E touch "/mnt/d/informatik/semester 12/masterarbeit/code/build/_deps/kassert-subbuild/kassert-populate-prefix/src/kassert-populate-stamp/kassert-populate-configure"

kassert-populate-prefix/src/kassert-populate-stamp/kassert-populate-download: kassert-populate-prefix/src/kassert-populate-stamp/kassert-populate-gitinfo.txt
kassert-populate-prefix/src/kassert-populate-stamp/kassert-populate-download: kassert-populate-prefix/src/kassert-populate-stamp/kassert-populate-mkdir
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir="/mnt/d/informatik/semester 12/masterarbeit/code/build/_deps/kassert-subbuild/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_4) "Performing download step (git clone) for 'kassert-populate'"
	cd "/mnt/d/informatik/semester 12/masterarbeit/code/build/_deps" && /usr/bin/cmake -P "/mnt/d/informatik/semester 12/masterarbeit/code/build/_deps/kassert-subbuild/kassert-populate-prefix/tmp/kassert-populate-gitclone.cmake"
	cd "/mnt/d/informatik/semester 12/masterarbeit/code/build/_deps" && /usr/bin/cmake -E touch "/mnt/d/informatik/semester 12/masterarbeit/code/build/_deps/kassert-subbuild/kassert-populate-prefix/src/kassert-populate-stamp/kassert-populate-download"

kassert-populate-prefix/src/kassert-populate-stamp/kassert-populate-install: kassert-populate-prefix/src/kassert-populate-stamp/kassert-populate-build
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir="/mnt/d/informatik/semester 12/masterarbeit/code/build/_deps/kassert-subbuild/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_5) "No install step for 'kassert-populate'"
	cd "/mnt/d/informatik/semester 12/masterarbeit/code/build/_deps/kassert-build" && /usr/bin/cmake -E echo_append
	cd "/mnt/d/informatik/semester 12/masterarbeit/code/build/_deps/kassert-build" && /usr/bin/cmake -E touch "/mnt/d/informatik/semester 12/masterarbeit/code/build/_deps/kassert-subbuild/kassert-populate-prefix/src/kassert-populate-stamp/kassert-populate-install"

kassert-populate-prefix/src/kassert-populate-stamp/kassert-populate-mkdir:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir="/mnt/d/informatik/semester 12/masterarbeit/code/build/_deps/kassert-subbuild/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_6) "Creating directories for 'kassert-populate'"
	/usr/bin/cmake -E make_directory "/mnt/d/informatik/semester 12/masterarbeit/code/build/_deps/kassert-src"
	/usr/bin/cmake -E make_directory "/mnt/d/informatik/semester 12/masterarbeit/code/build/_deps/kassert-build"
	/usr/bin/cmake -E make_directory "/mnt/d/informatik/semester 12/masterarbeit/code/build/_deps/kassert-subbuild/kassert-populate-prefix"
	/usr/bin/cmake -E make_directory "/mnt/d/informatik/semester 12/masterarbeit/code/build/_deps/kassert-subbuild/kassert-populate-prefix/tmp"
	/usr/bin/cmake -E make_directory "/mnt/d/informatik/semester 12/masterarbeit/code/build/_deps/kassert-subbuild/kassert-populate-prefix/src/kassert-populate-stamp"
	/usr/bin/cmake -E make_directory "/mnt/d/informatik/semester 12/masterarbeit/code/build/_deps/kassert-subbuild/kassert-populate-prefix/src"
	/usr/bin/cmake -E make_directory "/mnt/d/informatik/semester 12/masterarbeit/code/build/_deps/kassert-subbuild/kassert-populate-prefix/src/kassert-populate-stamp"
	/usr/bin/cmake -E touch "/mnt/d/informatik/semester 12/masterarbeit/code/build/_deps/kassert-subbuild/kassert-populate-prefix/src/kassert-populate-stamp/kassert-populate-mkdir"

kassert-populate-prefix/src/kassert-populate-stamp/kassert-populate-patch: kassert-populate-prefix/src/kassert-populate-stamp/kassert-populate-update
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir="/mnt/d/informatik/semester 12/masterarbeit/code/build/_deps/kassert-subbuild/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_7) "No patch step for 'kassert-populate'"
	/usr/bin/cmake -E echo_append
	/usr/bin/cmake -E touch "/mnt/d/informatik/semester 12/masterarbeit/code/build/_deps/kassert-subbuild/kassert-populate-prefix/src/kassert-populate-stamp/kassert-populate-patch"

kassert-populate-prefix/src/kassert-populate-stamp/kassert-populate-update:
.PHONY : kassert-populate-prefix/src/kassert-populate-stamp/kassert-populate-update

kassert-populate-prefix/src/kassert-populate-stamp/kassert-populate-test: kassert-populate-prefix/src/kassert-populate-stamp/kassert-populate-install
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir="/mnt/d/informatik/semester 12/masterarbeit/code/build/_deps/kassert-subbuild/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_8) "No test step for 'kassert-populate'"
	cd "/mnt/d/informatik/semester 12/masterarbeit/code/build/_deps/kassert-build" && /usr/bin/cmake -E echo_append
	cd "/mnt/d/informatik/semester 12/masterarbeit/code/build/_deps/kassert-build" && /usr/bin/cmake -E touch "/mnt/d/informatik/semester 12/masterarbeit/code/build/_deps/kassert-subbuild/kassert-populate-prefix/src/kassert-populate-stamp/kassert-populate-test"

kassert-populate-prefix/src/kassert-populate-stamp/kassert-populate-update: kassert-populate-prefix/src/kassert-populate-stamp/kassert-populate-download
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir="/mnt/d/informatik/semester 12/masterarbeit/code/build/_deps/kassert-subbuild/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_9) "Performing update step for 'kassert-populate'"
	cd "/mnt/d/informatik/semester 12/masterarbeit/code/build/_deps/kassert-src" && /usr/bin/cmake -P "/mnt/d/informatik/semester 12/masterarbeit/code/build/_deps/kassert-subbuild/kassert-populate-prefix/tmp/kassert-populate-gitupdate.cmake"

kassert-populate: CMakeFiles/kassert-populate
kassert-populate: CMakeFiles/kassert-populate-complete
kassert-populate: kassert-populate-prefix/src/kassert-populate-stamp/kassert-populate-build
kassert-populate: kassert-populate-prefix/src/kassert-populate-stamp/kassert-populate-configure
kassert-populate: kassert-populate-prefix/src/kassert-populate-stamp/kassert-populate-download
kassert-populate: kassert-populate-prefix/src/kassert-populate-stamp/kassert-populate-install
kassert-populate: kassert-populate-prefix/src/kassert-populate-stamp/kassert-populate-mkdir
kassert-populate: kassert-populate-prefix/src/kassert-populate-stamp/kassert-populate-patch
kassert-populate: kassert-populate-prefix/src/kassert-populate-stamp/kassert-populate-test
kassert-populate: kassert-populate-prefix/src/kassert-populate-stamp/kassert-populate-update
kassert-populate: CMakeFiles/kassert-populate.dir/build.make
.PHONY : kassert-populate

# Rule to build all files generated by this target.
CMakeFiles/kassert-populate.dir/build: kassert-populate
.PHONY : CMakeFiles/kassert-populate.dir/build

CMakeFiles/kassert-populate.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/kassert-populate.dir/cmake_clean.cmake
.PHONY : CMakeFiles/kassert-populate.dir/clean

CMakeFiles/kassert-populate.dir/depend:
	cd "/mnt/d/informatik/semester 12/masterarbeit/code/build/_deps/kassert-subbuild" && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" "/mnt/d/informatik/semester 12/masterarbeit/code/build/_deps/kassert-subbuild" "/mnt/d/informatik/semester 12/masterarbeit/code/build/_deps/kassert-subbuild" "/mnt/d/informatik/semester 12/masterarbeit/code/build/_deps/kassert-subbuild" "/mnt/d/informatik/semester 12/masterarbeit/code/build/_deps/kassert-subbuild" "/mnt/d/informatik/semester 12/masterarbeit/code/build/_deps/kassert-subbuild/CMakeFiles/kassert-populate.dir/DependInfo.cmake" --color=$(COLOR)
.PHONY : CMakeFiles/kassert-populate.dir/depend

