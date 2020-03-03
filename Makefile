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
CMAKE_COMMAND = /cvmfs/sft.cern.ch/lcg/releases/CMake/3.11.1-daf3a/x86_64-centos7-gcc8-opt/bin/cmake

# The command to remove a file.
RM = /cvmfs/sft.cern.ch/lcg/releases/CMake/3.11.1-daf3a/x86_64-centos7-gcc8-opt/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /user/glegu/testsim/modusim

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /user/glegu/testsim/modusim

#=============================================================================
# Targets provided globally by CMake.

# Special rule for the target install/strip
install/strip: preinstall
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Installing the project stripped..."
	/cvmfs/sft.cern.ch/lcg/releases/CMake/3.11.1-daf3a/x86_64-centos7-gcc8-opt/bin/cmake -DCMAKE_INSTALL_DO_STRIP=1 -P cmake_install.cmake
.PHONY : install/strip

# Special rule for the target install/strip
install/strip/fast: preinstall/fast
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Installing the project stripped..."
	/cvmfs/sft.cern.ch/lcg/releases/CMake/3.11.1-daf3a/x86_64-centos7-gcc8-opt/bin/cmake -DCMAKE_INSTALL_DO_STRIP=1 -P cmake_install.cmake
.PHONY : install/strip/fast

# Special rule for the target edit_cache
edit_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake cache editor..."
	/cvmfs/sft.cern.ch/lcg/releases/CMake/3.11.1-daf3a/x86_64-centos7-gcc8-opt/bin/ccmake -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : edit_cache

# Special rule for the target edit_cache
edit_cache/fast: edit_cache

.PHONY : edit_cache/fast

# Special rule for the target rebuild_cache
rebuild_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake to regenerate build system..."
	/cvmfs/sft.cern.ch/lcg/releases/CMake/3.11.1-daf3a/x86_64-centos7-gcc8-opt/bin/cmake -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : rebuild_cache

# Special rule for the target rebuild_cache
rebuild_cache/fast: rebuild_cache

.PHONY : rebuild_cache/fast

# Special rule for the target list_install_components
list_install_components:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Available install components are: \"Unspecified\""
.PHONY : list_install_components

# Special rule for the target list_install_components
list_install_components/fast: list_install_components

.PHONY : list_install_components/fast

# Special rule for the target install/local
install/local: preinstall
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Installing only the local directory..."
	/cvmfs/sft.cern.ch/lcg/releases/CMake/3.11.1-daf3a/x86_64-centos7-gcc8-opt/bin/cmake -DCMAKE_INSTALL_LOCAL_ONLY=1 -P cmake_install.cmake
.PHONY : install/local

# Special rule for the target install/local
install/local/fast: preinstall/fast
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Installing only the local directory..."
	/cvmfs/sft.cern.ch/lcg/releases/CMake/3.11.1-daf3a/x86_64-centos7-gcc8-opt/bin/cmake -DCMAKE_INSTALL_LOCAL_ONLY=1 -P cmake_install.cmake
.PHONY : install/local/fast

# Special rule for the target install
install: preinstall
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Install the project..."
	/cvmfs/sft.cern.ch/lcg/releases/CMake/3.11.1-daf3a/x86_64-centos7-gcc8-opt/bin/cmake -P cmake_install.cmake
.PHONY : install

# Special rule for the target install
install/fast: preinstall/fast
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Install the project..."
	/cvmfs/sft.cern.ch/lcg/releases/CMake/3.11.1-daf3a/x86_64-centos7-gcc8-opt/bin/cmake -P cmake_install.cmake
.PHONY : install/fast

# The main all target
all: cmake_check_build_system
	$(CMAKE_COMMAND) -E cmake_progress_start /user/glegu/testsim/modusim/CMakeFiles /user/glegu/testsim/modusim/CMakeFiles/progress.marks
	$(MAKE) -f CMakeFiles/Makefile2 all
	$(CMAKE_COMMAND) -E cmake_progress_start /user/glegu/testsim/modusim/CMakeFiles 0
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
# Target rules for targets named G4simu

# Build rule for target.
G4simu: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 G4simu
.PHONY : G4simu

# fast build rule for target.
G4simu/fast:
	$(MAKE) -f CMakeFiles/G4simu.dir/build.make CMakeFiles/G4simu.dir/build
.PHONY : G4simu/fast

G4simu.o: G4simu.cc.o

.PHONY : G4simu.o

# target to build an object file
G4simu.cc.o:
	$(MAKE) -f CMakeFiles/G4simu.dir/build.make CMakeFiles/G4simu.dir/G4simu.cc.o
.PHONY : G4simu.cc.o

G4simu.i: G4simu.cc.i

.PHONY : G4simu.i

# target to preprocess a source file
G4simu.cc.i:
	$(MAKE) -f CMakeFiles/G4simu.dir/build.make CMakeFiles/G4simu.dir/G4simu.cc.i
.PHONY : G4simu.cc.i

G4simu.s: G4simu.cc.s

.PHONY : G4simu.s

# target to generate assembly for a file
G4simu.cc.s:
	$(MAKE) -f CMakeFiles/G4simu.dir/build.make CMakeFiles/G4simu.dir/G4simu.cc.s
.PHONY : G4simu.cc.s

src/AnalysisManager.o: src/AnalysisManager.cc.o

.PHONY : src/AnalysisManager.o

# target to build an object file
src/AnalysisManager.cc.o:
	$(MAKE) -f CMakeFiles/G4simu.dir/build.make CMakeFiles/G4simu.dir/src/AnalysisManager.cc.o
.PHONY : src/AnalysisManager.cc.o

src/AnalysisManager.i: src/AnalysisManager.cc.i

.PHONY : src/AnalysisManager.i

# target to preprocess a source file
src/AnalysisManager.cc.i:
	$(MAKE) -f CMakeFiles/G4simu.dir/build.make CMakeFiles/G4simu.dir/src/AnalysisManager.cc.i
.PHONY : src/AnalysisManager.cc.i

src/AnalysisManager.s: src/AnalysisManager.cc.s

.PHONY : src/AnalysisManager.s

# target to generate assembly for a file
src/AnalysisManager.cc.s:
	$(MAKE) -f CMakeFiles/G4simu.dir/build.make CMakeFiles/G4simu.dir/src/AnalysisManager.cc.s
.PHONY : src/AnalysisManager.cc.s

src/AnalysisMessenger.o: src/AnalysisMessenger.cc.o

.PHONY : src/AnalysisMessenger.o

# target to build an object file
src/AnalysisMessenger.cc.o:
	$(MAKE) -f CMakeFiles/G4simu.dir/build.make CMakeFiles/G4simu.dir/src/AnalysisMessenger.cc.o
.PHONY : src/AnalysisMessenger.cc.o

src/AnalysisMessenger.i: src/AnalysisMessenger.cc.i

.PHONY : src/AnalysisMessenger.i

# target to preprocess a source file
src/AnalysisMessenger.cc.i:
	$(MAKE) -f CMakeFiles/G4simu.dir/build.make CMakeFiles/G4simu.dir/src/AnalysisMessenger.cc.i
.PHONY : src/AnalysisMessenger.cc.i

src/AnalysisMessenger.s: src/AnalysisMessenger.cc.s

.PHONY : src/AnalysisMessenger.s

# target to generate assembly for a file
src/AnalysisMessenger.cc.s:
	$(MAKE) -f CMakeFiles/G4simu.dir/build.make CMakeFiles/G4simu.dir/src/AnalysisMessenger.cc.s
.PHONY : src/AnalysisMessenger.cc.s

src/DetectorConstruction.o: src/DetectorConstruction.cc.o

.PHONY : src/DetectorConstruction.o

# target to build an object file
src/DetectorConstruction.cc.o:
	$(MAKE) -f CMakeFiles/G4simu.dir/build.make CMakeFiles/G4simu.dir/src/DetectorConstruction.cc.o
.PHONY : src/DetectorConstruction.cc.o

src/DetectorConstruction.i: src/DetectorConstruction.cc.i

.PHONY : src/DetectorConstruction.i

# target to preprocess a source file
src/DetectorConstruction.cc.i:
	$(MAKE) -f CMakeFiles/G4simu.dir/build.make CMakeFiles/G4simu.dir/src/DetectorConstruction.cc.i
.PHONY : src/DetectorConstruction.cc.i

src/DetectorConstruction.s: src/DetectorConstruction.cc.s

.PHONY : src/DetectorConstruction.s

# target to generate assembly for a file
src/DetectorConstruction.cc.s:
	$(MAKE) -f CMakeFiles/G4simu.dir/build.make CMakeFiles/G4simu.dir/src/DetectorConstruction.cc.s
.PHONY : src/DetectorConstruction.cc.s

src/DetectorMessenger.o: src/DetectorMessenger.cc.o

.PHONY : src/DetectorMessenger.o

# target to build an object file
src/DetectorMessenger.cc.o:
	$(MAKE) -f CMakeFiles/G4simu.dir/build.make CMakeFiles/G4simu.dir/src/DetectorMessenger.cc.o
.PHONY : src/DetectorMessenger.cc.o

src/DetectorMessenger.i: src/DetectorMessenger.cc.i

.PHONY : src/DetectorMessenger.i

# target to preprocess a source file
src/DetectorMessenger.cc.i:
	$(MAKE) -f CMakeFiles/G4simu.dir/build.make CMakeFiles/G4simu.dir/src/DetectorMessenger.cc.i
.PHONY : src/DetectorMessenger.cc.i

src/DetectorMessenger.s: src/DetectorMessenger.cc.s

.PHONY : src/DetectorMessenger.s

# target to generate assembly for a file
src/DetectorMessenger.cc.s:
	$(MAKE) -f CMakeFiles/G4simu.dir/build.make CMakeFiles/G4simu.dir/src/DetectorMessenger.cc.s
.PHONY : src/DetectorMessenger.cc.s

src/EventAction.o: src/EventAction.cc.o

.PHONY : src/EventAction.o

# target to build an object file
src/EventAction.cc.o:
	$(MAKE) -f CMakeFiles/G4simu.dir/build.make CMakeFiles/G4simu.dir/src/EventAction.cc.o
.PHONY : src/EventAction.cc.o

src/EventAction.i: src/EventAction.cc.i

.PHONY : src/EventAction.i

# target to preprocess a source file
src/EventAction.cc.i:
	$(MAKE) -f CMakeFiles/G4simu.dir/build.make CMakeFiles/G4simu.dir/src/EventAction.cc.i
.PHONY : src/EventAction.cc.i

src/EventAction.s: src/EventAction.cc.s

.PHONY : src/EventAction.s

# target to generate assembly for a file
src/EventAction.cc.s:
	$(MAKE) -f CMakeFiles/G4simu.dir/build.make CMakeFiles/G4simu.dir/src/EventAction.cc.s
.PHONY : src/EventAction.cc.s

src/EventData.o: src/EventData.cc.o

.PHONY : src/EventData.o

# target to build an object file
src/EventData.cc.o:
	$(MAKE) -f CMakeFiles/G4simu.dir/build.make CMakeFiles/G4simu.dir/src/EventData.cc.o
.PHONY : src/EventData.cc.o

src/EventData.i: src/EventData.cc.i

.PHONY : src/EventData.i

# target to preprocess a source file
src/EventData.cc.i:
	$(MAKE) -f CMakeFiles/G4simu.dir/build.make CMakeFiles/G4simu.dir/src/EventData.cc.i
.PHONY : src/EventData.cc.i

src/EventData.s: src/EventData.cc.s

.PHONY : src/EventData.s

# target to generate assembly for a file
src/EventData.cc.s:
	$(MAKE) -f CMakeFiles/G4simu.dir/build.make CMakeFiles/G4simu.dir/src/EventData.cc.s
.PHONY : src/EventData.cc.s

src/G4KleinNishinaModel.o: src/G4KleinNishinaModel.cc.o

.PHONY : src/G4KleinNishinaModel.o

# target to build an object file
src/G4KleinNishinaModel.cc.o:
	$(MAKE) -f CMakeFiles/G4simu.dir/build.make CMakeFiles/G4simu.dir/src/G4KleinNishinaModel.cc.o
.PHONY : src/G4KleinNishinaModel.cc.o

src/G4KleinNishinaModel.i: src/G4KleinNishinaModel.cc.i

.PHONY : src/G4KleinNishinaModel.i

# target to preprocess a source file
src/G4KleinNishinaModel.cc.i:
	$(MAKE) -f CMakeFiles/G4simu.dir/build.make CMakeFiles/G4simu.dir/src/G4KleinNishinaModel.cc.i
.PHONY : src/G4KleinNishinaModel.cc.i

src/G4KleinNishinaModel.s: src/G4KleinNishinaModel.cc.s

.PHONY : src/G4KleinNishinaModel.s

# target to generate assembly for a file
src/G4KleinNishinaModel.cc.s:
	$(MAKE) -f CMakeFiles/G4simu.dir/build.make CMakeFiles/G4simu.dir/src/G4KleinNishinaModel.cc.s
.PHONY : src/G4KleinNishinaModel.cc.s

src/G4LowEPComptonModel.o: src/G4LowEPComptonModel.cc.o

.PHONY : src/G4LowEPComptonModel.o

# target to build an object file
src/G4LowEPComptonModel.cc.o:
	$(MAKE) -f CMakeFiles/G4simu.dir/build.make CMakeFiles/G4simu.dir/src/G4LowEPComptonModel.cc.o
.PHONY : src/G4LowEPComptonModel.cc.o

src/G4LowEPComptonModel.i: src/G4LowEPComptonModel.cc.i

.PHONY : src/G4LowEPComptonModel.i

# target to preprocess a source file
src/G4LowEPComptonModel.cc.i:
	$(MAKE) -f CMakeFiles/G4simu.dir/build.make CMakeFiles/G4simu.dir/src/G4LowEPComptonModel.cc.i
.PHONY : src/G4LowEPComptonModel.cc.i

src/G4LowEPComptonModel.s: src/G4LowEPComptonModel.cc.s

.PHONY : src/G4LowEPComptonModel.s

# target to generate assembly for a file
src/G4LowEPComptonModel.cc.s:
	$(MAKE) -f CMakeFiles/G4simu.dir/build.make CMakeFiles/G4simu.dir/src/G4LowEPComptonModel.cc.s
.PHONY : src/G4LowEPComptonModel.cc.s

src/Particle.o: src/Particle.cc.o

.PHONY : src/Particle.o

# target to build an object file
src/Particle.cc.o:
	$(MAKE) -f CMakeFiles/G4simu.dir/build.make CMakeFiles/G4simu.dir/src/Particle.cc.o
.PHONY : src/Particle.cc.o

src/Particle.i: src/Particle.cc.i

.PHONY : src/Particle.i

# target to preprocess a source file
src/Particle.cc.i:
	$(MAKE) -f CMakeFiles/G4simu.dir/build.make CMakeFiles/G4simu.dir/src/Particle.cc.i
.PHONY : src/Particle.cc.i

src/Particle.s: src/Particle.cc.s

.PHONY : src/Particle.s

# target to generate assembly for a file
src/Particle.cc.s:
	$(MAKE) -f CMakeFiles/G4simu.dir/build.make CMakeFiles/G4simu.dir/src/Particle.cc.s
.PHONY : src/Particle.cc.s

src/PhysicsList.o: src/PhysicsList.cc.o

.PHONY : src/PhysicsList.o

# target to build an object file
src/PhysicsList.cc.o:
	$(MAKE) -f CMakeFiles/G4simu.dir/build.make CMakeFiles/G4simu.dir/src/PhysicsList.cc.o
.PHONY : src/PhysicsList.cc.o

src/PhysicsList.i: src/PhysicsList.cc.i

.PHONY : src/PhysicsList.i

# target to preprocess a source file
src/PhysicsList.cc.i:
	$(MAKE) -f CMakeFiles/G4simu.dir/build.make CMakeFiles/G4simu.dir/src/PhysicsList.cc.i
.PHONY : src/PhysicsList.cc.i

src/PhysicsList.s: src/PhysicsList.cc.s

.PHONY : src/PhysicsList.s

# target to generate assembly for a file
src/PhysicsList.cc.s:
	$(MAKE) -f CMakeFiles/G4simu.dir/build.make CMakeFiles/G4simu.dir/src/PhysicsList.cc.s
.PHONY : src/PhysicsList.cc.s

src/PhysicsMessenger.o: src/PhysicsMessenger.cc.o

.PHONY : src/PhysicsMessenger.o

# target to build an object file
src/PhysicsMessenger.cc.o:
	$(MAKE) -f CMakeFiles/G4simu.dir/build.make CMakeFiles/G4simu.dir/src/PhysicsMessenger.cc.o
.PHONY : src/PhysicsMessenger.cc.o

src/PhysicsMessenger.i: src/PhysicsMessenger.cc.i

.PHONY : src/PhysicsMessenger.i

# target to preprocess a source file
src/PhysicsMessenger.cc.i:
	$(MAKE) -f CMakeFiles/G4simu.dir/build.make CMakeFiles/G4simu.dir/src/PhysicsMessenger.cc.i
.PHONY : src/PhysicsMessenger.cc.i

src/PhysicsMessenger.s: src/PhysicsMessenger.cc.s

.PHONY : src/PhysicsMessenger.s

# target to generate assembly for a file
src/PhysicsMessenger.cc.s:
	$(MAKE) -f CMakeFiles/G4simu.dir/build.make CMakeFiles/G4simu.dir/src/PhysicsMessenger.cc.s
.PHONY : src/PhysicsMessenger.cc.s

src/PrimaryGeneratorAction.o: src/PrimaryGeneratorAction.cc.o

.PHONY : src/PrimaryGeneratorAction.o

# target to build an object file
src/PrimaryGeneratorAction.cc.o:
	$(MAKE) -f CMakeFiles/G4simu.dir/build.make CMakeFiles/G4simu.dir/src/PrimaryGeneratorAction.cc.o
.PHONY : src/PrimaryGeneratorAction.cc.o

src/PrimaryGeneratorAction.i: src/PrimaryGeneratorAction.cc.i

.PHONY : src/PrimaryGeneratorAction.i

# target to preprocess a source file
src/PrimaryGeneratorAction.cc.i:
	$(MAKE) -f CMakeFiles/G4simu.dir/build.make CMakeFiles/G4simu.dir/src/PrimaryGeneratorAction.cc.i
.PHONY : src/PrimaryGeneratorAction.cc.i

src/PrimaryGeneratorAction.s: src/PrimaryGeneratorAction.cc.s

.PHONY : src/PrimaryGeneratorAction.s

# target to generate assembly for a file
src/PrimaryGeneratorAction.cc.s:
	$(MAKE) -f CMakeFiles/G4simu.dir/build.make CMakeFiles/G4simu.dir/src/PrimaryGeneratorAction.cc.s
.PHONY : src/PrimaryGeneratorAction.cc.s

src/RunAction.o: src/RunAction.cc.o

.PHONY : src/RunAction.o

# target to build an object file
src/RunAction.cc.o:
	$(MAKE) -f CMakeFiles/G4simu.dir/build.make CMakeFiles/G4simu.dir/src/RunAction.cc.o
.PHONY : src/RunAction.cc.o

src/RunAction.i: src/RunAction.cc.i

.PHONY : src/RunAction.i

# target to preprocess a source file
src/RunAction.cc.i:
	$(MAKE) -f CMakeFiles/G4simu.dir/build.make CMakeFiles/G4simu.dir/src/RunAction.cc.i
.PHONY : src/RunAction.cc.i

src/RunAction.s: src/RunAction.cc.s

.PHONY : src/RunAction.s

# target to generate assembly for a file
src/RunAction.cc.s:
	$(MAKE) -f CMakeFiles/G4simu.dir/build.make CMakeFiles/G4simu.dir/src/RunAction.cc.s
.PHONY : src/RunAction.cc.s

src/RunActionMessenger.o: src/RunActionMessenger.cc.o

.PHONY : src/RunActionMessenger.o

# target to build an object file
src/RunActionMessenger.cc.o:
	$(MAKE) -f CMakeFiles/G4simu.dir/build.make CMakeFiles/G4simu.dir/src/RunActionMessenger.cc.o
.PHONY : src/RunActionMessenger.cc.o

src/RunActionMessenger.i: src/RunActionMessenger.cc.i

.PHONY : src/RunActionMessenger.i

# target to preprocess a source file
src/RunActionMessenger.cc.i:
	$(MAKE) -f CMakeFiles/G4simu.dir/build.make CMakeFiles/G4simu.dir/src/RunActionMessenger.cc.i
.PHONY : src/RunActionMessenger.cc.i

src/RunActionMessenger.s: src/RunActionMessenger.cc.s

.PHONY : src/RunActionMessenger.s

# target to generate assembly for a file
src/RunActionMessenger.cc.s:
	$(MAKE) -f CMakeFiles/G4simu.dir/build.make CMakeFiles/G4simu.dir/src/RunActionMessenger.cc.s
.PHONY : src/RunActionMessenger.cc.s

src/SensitiveDetector.o: src/SensitiveDetector.cc.o

.PHONY : src/SensitiveDetector.o

# target to build an object file
src/SensitiveDetector.cc.o:
	$(MAKE) -f CMakeFiles/G4simu.dir/build.make CMakeFiles/G4simu.dir/src/SensitiveDetector.cc.o
.PHONY : src/SensitiveDetector.cc.o

src/SensitiveDetector.i: src/SensitiveDetector.cc.i

.PHONY : src/SensitiveDetector.i

# target to preprocess a source file
src/SensitiveDetector.cc.i:
	$(MAKE) -f CMakeFiles/G4simu.dir/build.make CMakeFiles/G4simu.dir/src/SensitiveDetector.cc.i
.PHONY : src/SensitiveDetector.cc.i

src/SensitiveDetector.s: src/SensitiveDetector.cc.s

.PHONY : src/SensitiveDetector.s

# target to generate assembly for a file
src/SensitiveDetector.cc.s:
	$(MAKE) -f CMakeFiles/G4simu.dir/build.make CMakeFiles/G4simu.dir/src/SensitiveDetector.cc.s
.PHONY : src/SensitiveDetector.cc.s

src/StackingAction.o: src/StackingAction.cc.o

.PHONY : src/StackingAction.o

# target to build an object file
src/StackingAction.cc.o:
	$(MAKE) -f CMakeFiles/G4simu.dir/build.make CMakeFiles/G4simu.dir/src/StackingAction.cc.o
.PHONY : src/StackingAction.cc.o

src/StackingAction.i: src/StackingAction.cc.i

.PHONY : src/StackingAction.i

# target to preprocess a source file
src/StackingAction.cc.i:
	$(MAKE) -f CMakeFiles/G4simu.dir/build.make CMakeFiles/G4simu.dir/src/StackingAction.cc.i
.PHONY : src/StackingAction.cc.i

src/StackingAction.s: src/StackingAction.cc.s

.PHONY : src/StackingAction.s

# target to generate assembly for a file
src/StackingAction.cc.s:
	$(MAKE) -f CMakeFiles/G4simu.dir/build.make CMakeFiles/G4simu.dir/src/StackingAction.cc.s
.PHONY : src/StackingAction.cc.s

src/SteppingAction.o: src/SteppingAction.cc.o

.PHONY : src/SteppingAction.o

# target to build an object file
src/SteppingAction.cc.o:
	$(MAKE) -f CMakeFiles/G4simu.dir/build.make CMakeFiles/G4simu.dir/src/SteppingAction.cc.o
.PHONY : src/SteppingAction.cc.o

src/SteppingAction.i: src/SteppingAction.cc.i

.PHONY : src/SteppingAction.i

# target to preprocess a source file
src/SteppingAction.cc.i:
	$(MAKE) -f CMakeFiles/G4simu.dir/build.make CMakeFiles/G4simu.dir/src/SteppingAction.cc.i
.PHONY : src/SteppingAction.cc.i

src/SteppingAction.s: src/SteppingAction.cc.s

.PHONY : src/SteppingAction.s

# target to generate assembly for a file
src/SteppingAction.cc.s:
	$(MAKE) -f CMakeFiles/G4simu.dir/build.make CMakeFiles/G4simu.dir/src/SteppingAction.cc.s
.PHONY : src/SteppingAction.cc.s

src/fileMerger.o: src/fileMerger.cc.o

.PHONY : src/fileMerger.o

# target to build an object file
src/fileMerger.cc.o:
	$(MAKE) -f CMakeFiles/G4simu.dir/build.make CMakeFiles/G4simu.dir/src/fileMerger.cc.o
.PHONY : src/fileMerger.cc.o

src/fileMerger.i: src/fileMerger.cc.i

.PHONY : src/fileMerger.i

# target to preprocess a source file
src/fileMerger.cc.i:
	$(MAKE) -f CMakeFiles/G4simu.dir/build.make CMakeFiles/G4simu.dir/src/fileMerger.cc.i
.PHONY : src/fileMerger.cc.i

src/fileMerger.s: src/fileMerger.cc.s

.PHONY : src/fileMerger.s

# target to generate assembly for a file
src/fileMerger.cc.s:
	$(MAKE) -f CMakeFiles/G4simu.dir/build.make CMakeFiles/G4simu.dir/src/fileMerger.cc.s
.PHONY : src/fileMerger.cc.s

src/stdHit.o: src/stdHit.cc.o

.PHONY : src/stdHit.o

# target to build an object file
src/stdHit.cc.o:
	$(MAKE) -f CMakeFiles/G4simu.dir/build.make CMakeFiles/G4simu.dir/src/stdHit.cc.o
.PHONY : src/stdHit.cc.o

src/stdHit.i: src/stdHit.cc.i

.PHONY : src/stdHit.i

# target to preprocess a source file
src/stdHit.cc.i:
	$(MAKE) -f CMakeFiles/G4simu.dir/build.make CMakeFiles/G4simu.dir/src/stdHit.cc.i
.PHONY : src/stdHit.cc.i

src/stdHit.s: src/stdHit.cc.s

.PHONY : src/stdHit.s

# target to generate assembly for a file
src/stdHit.cc.s:
	$(MAKE) -f CMakeFiles/G4simu.dir/build.make CMakeFiles/G4simu.dir/src/stdHit.cc.s
.PHONY : src/stdHit.cc.s

# Help Target
help:
	@echo "The following are some of the valid targets for this Makefile:"
	@echo "... all (the default if no target is provided)"
	@echo "... clean"
	@echo "... depend"
	@echo "... install/strip"
	@echo "... edit_cache"
	@echo "... G4simu"
	@echo "... rebuild_cache"
	@echo "... list_install_components"
	@echo "... install/local"
	@echo "... install"
	@echo "... G4simu.o"
	@echo "... G4simu.i"
	@echo "... G4simu.s"
	@echo "... src/AnalysisManager.o"
	@echo "... src/AnalysisManager.i"
	@echo "... src/AnalysisManager.s"
	@echo "... src/AnalysisMessenger.o"
	@echo "... src/AnalysisMessenger.i"
	@echo "... src/AnalysisMessenger.s"
	@echo "... src/DetectorConstruction.o"
	@echo "... src/DetectorConstruction.i"
	@echo "... src/DetectorConstruction.s"
	@echo "... src/DetectorMessenger.o"
	@echo "... src/DetectorMessenger.i"
	@echo "... src/DetectorMessenger.s"
	@echo "... src/EventAction.o"
	@echo "... src/EventAction.i"
	@echo "... src/EventAction.s"
	@echo "... src/EventData.o"
	@echo "... src/EventData.i"
	@echo "... src/EventData.s"
	@echo "... src/G4KleinNishinaModel.o"
	@echo "... src/G4KleinNishinaModel.i"
	@echo "... src/G4KleinNishinaModel.s"
	@echo "... src/G4LowEPComptonModel.o"
	@echo "... src/G4LowEPComptonModel.i"
	@echo "... src/G4LowEPComptonModel.s"
	@echo "... src/Particle.o"
	@echo "... src/Particle.i"
	@echo "... src/Particle.s"
	@echo "... src/PhysicsList.o"
	@echo "... src/PhysicsList.i"
	@echo "... src/PhysicsList.s"
	@echo "... src/PhysicsMessenger.o"
	@echo "... src/PhysicsMessenger.i"
	@echo "... src/PhysicsMessenger.s"
	@echo "... src/PrimaryGeneratorAction.o"
	@echo "... src/PrimaryGeneratorAction.i"
	@echo "... src/PrimaryGeneratorAction.s"
	@echo "... src/RunAction.o"
	@echo "... src/RunAction.i"
	@echo "... src/RunAction.s"
	@echo "... src/RunActionMessenger.o"
	@echo "... src/RunActionMessenger.i"
	@echo "... src/RunActionMessenger.s"
	@echo "... src/SensitiveDetector.o"
	@echo "... src/SensitiveDetector.i"
	@echo "... src/SensitiveDetector.s"
	@echo "... src/StackingAction.o"
	@echo "... src/StackingAction.i"
	@echo "... src/StackingAction.s"
	@echo "... src/SteppingAction.o"
	@echo "... src/SteppingAction.i"
	@echo "... src/SteppingAction.s"
	@echo "... src/fileMerger.o"
	@echo "... src/fileMerger.i"
	@echo "... src/fileMerger.s"
	@echo "... src/stdHit.o"
	@echo "... src/stdHit.i"
	@echo "... src/stdHit.s"
.PHONY : help



#=============================================================================
# Special targets to cleanup operation of make.

# Special rule to run CMake to check the build system integrity.
# No rule that depends on this can have commands that come from listfiles
# because they might be regenerated.
cmake_check_build_system:
	$(CMAKE_COMMAND) -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 0
.PHONY : cmake_check_build_system

