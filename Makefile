# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.11

# Default target executed when no arguments are given to make.
default_target: all

.PHONY : default_target

# Allow only one "make -f Makefile2" at a time, but pass parallelism.
.NOTPARALLEL:


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
CMAKE_COMMAND = /usr/local/bin/cmake

# The command to remove a file.
RM = /usr/local/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /mnt/d/SemII/APT/P2/P2

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /mnt/d/SemII/APT/P2/P2

#=============================================================================
# Targets provided globally by CMake.

# Special rule for the target rebuild_cache
rebuild_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake to regenerate build system..."
	/usr/local/bin/cmake -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : rebuild_cache

# Special rule for the target rebuild_cache
rebuild_cache/fast: rebuild_cache

.PHONY : rebuild_cache/fast

# Special rule for the target edit_cache
edit_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "No interactive CMake dialog available..."
	/usr/local/bin/cmake -E echo No\ interactive\ CMake\ dialog\ available.
.PHONY : edit_cache

# Special rule for the target edit_cache
edit_cache/fast: edit_cache

.PHONY : edit_cache/fast

# The main all target
all: cmake_check_build_system
	$(CMAKE_COMMAND) -E cmake_progress_start /mnt/d/SemII/APT/P2/P2/CMakeFiles /mnt/d/SemII/APT/P2/P2/CMakeFiles/progress.marks
	$(MAKE) -f CMakeFiles/Makefile2 all
	$(CMAKE_COMMAND) -E cmake_progress_start /mnt/d/SemII/APT/P2/P2/CMakeFiles 0
.PHONY : all

# The main clean target
clean:
	$(MAKE) -f CMakeFiles/Makefile2 clean
.PHONY : clean

# The main clean target
clean/fast: clean

.PHONY : clean/fast

# Prepare targets for installation.
preinstall: all
	$(MAKE) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall

# Prepare targets for installation.
preinstall/fast:
	$(MAKE) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall/fast

# clear depends
depend:
	$(CMAKE_COMMAND) -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 1
.PHONY : depend

#=============================================================================
# Target rules for targets named P2

# Build rule for target.
P2: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 P2
.PHONY : P2

# fast build rule for target.
P2/fast:
	$(MAKE) -f CMakeFiles/P2.dir/build.make CMakeFiles/P2.dir/build
.PHONY : P2/fast

Complex.o: Complex.cc.o

.PHONY : Complex.o

# target to build an object file
Complex.cc.o:
	$(MAKE) -f CMakeFiles/P2.dir/build.make CMakeFiles/P2.dir/Complex.cc.o
.PHONY : Complex.cc.o

Complex.i: Complex.cc.i

.PHONY : Complex.i

# target to preprocess a source file
Complex.cc.i:
	$(MAKE) -f CMakeFiles/P2.dir/build.make CMakeFiles/P2.dir/Complex.cc.i
.PHONY : Complex.cc.i

Complex.s: Complex.cc.s

.PHONY : Complex.s

# target to generate assembly for a file
Complex.cc.s:
	$(MAKE) -f CMakeFiles/P2.dir/build.make CMakeFiles/P2.dir/Complex.cc.s
.PHONY : Complex.cc.s

InputImage.o: InputImage.cc.o

.PHONY : InputImage.o

# target to build an object file
InputImage.cc.o:
	$(MAKE) -f CMakeFiles/P2.dir/build.make CMakeFiles/P2.dir/InputImage.cc.o
.PHONY : InputImage.cc.o

InputImage.i: InputImage.cc.i

.PHONY : InputImage.i

# target to preprocess a source file
InputImage.cc.i:
	$(MAKE) -f CMakeFiles/P2.dir/build.make CMakeFiles/P2.dir/InputImage.cc.i
.PHONY : InputImage.cc.i

InputImage.s: InputImage.cc.s

.PHONY : InputImage.s

# target to generate assembly for a file
InputImage.cc.s:
	$(MAKE) -f CMakeFiles/P2.dir/build.make CMakeFiles/P2.dir/InputImage.cc.s
.PHONY : InputImage.cc.s

fft2d.o: fft2d.cc.o

.PHONY : fft2d.o

# target to build an object file
fft2d.cc.o:
	$(MAKE) -f CMakeFiles/P2.dir/build.make CMakeFiles/P2.dir/fft2d.cc.o
.PHONY : fft2d.cc.o

fft2d.i: fft2d.cc.i

.PHONY : fft2d.i

# target to preprocess a source file
fft2d.cc.i:
	$(MAKE) -f CMakeFiles/P2.dir/build.make CMakeFiles/P2.dir/fft2d.cc.i
.PHONY : fft2d.cc.i

fft2d.s: fft2d.cc.s

.PHONY : fft2d.s

# target to generate assembly for a file
fft2d.cc.s:
	$(MAKE) -f CMakeFiles/P2.dir/build.make CMakeFiles/P2.dir/fft2d.cc.s
.PHONY : fft2d.cc.s

# Help Target
help:
	@echo "The following are some of the valid targets for this Makefile:"
	@echo "... all (the default if no target is provided)"
	@echo "... clean"
	@echo "... depend"
	@echo "... rebuild_cache"
	@echo "... edit_cache"
	@echo "... P2"
	@echo "... Complex.o"
	@echo "... Complex.i"
	@echo "... Complex.s"
	@echo "... InputImage.o"
	@echo "... InputImage.i"
	@echo "... InputImage.s"
	@echo "... fft2d.o"
	@echo "... fft2d.i"
	@echo "... fft2d.s"
.PHONY : help



#=============================================================================
# Special targets to cleanup operation of make.

# Special rule to run CMake to check the build system integrity.
# No rule that depends on this can have commands that come from listfiles
# because they might be regenerated.
cmake_check_build_system:
	$(CMAKE_COMMAND) -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 0
.PHONY : cmake_check_build_system
