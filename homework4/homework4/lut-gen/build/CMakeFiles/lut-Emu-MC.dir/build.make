# CMAKE generated file: DO NOT EDIT!
# Generated by "MinGW Makefiles" Generator, CMake Version 3.27

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

SHELL = cmd.exe

# The CMake executable.
CMAKE_COMMAND = "C:\Program Files\CMake\bin\cmake.exe"

# The command to remove a file.
RM = "C:\Program Files\CMake\bin\cmake.exe" -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = D:\Desktop\GAMES\Games202\GAMES202_HW\homework4\homework4\lut-gen

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = D:\Desktop\GAMES\Games202\GAMES202_HW\homework4\homework4\lut-gen\build

# Include any dependencies generated for this target.
include CMakeFiles/lut-Emu-MC.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/lut-Emu-MC.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/lut-Emu-MC.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/lut-Emu-MC.dir/flags.make

CMakeFiles/lut-Emu-MC.dir/Emu_MC.cpp.obj: CMakeFiles/lut-Emu-MC.dir/flags.make
CMakeFiles/lut-Emu-MC.dir/Emu_MC.cpp.obj: CMakeFiles/lut-Emu-MC.dir/includes_CXX.rsp
CMakeFiles/lut-Emu-MC.dir/Emu_MC.cpp.obj: D:/Desktop/GAMES/Games202/GAMES202_HW/homework4/homework4/lut-gen/Emu_MC.cpp
CMakeFiles/lut-Emu-MC.dir/Emu_MC.cpp.obj: CMakeFiles/lut-Emu-MC.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=D:\Desktop\GAMES\Games202\GAMES202_HW\homework4\homework4\lut-gen\build\CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/lut-Emu-MC.dir/Emu_MC.cpp.obj"
	D:\Program\MinGW\bin\c++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/lut-Emu-MC.dir/Emu_MC.cpp.obj -MF CMakeFiles\lut-Emu-MC.dir\Emu_MC.cpp.obj.d -o CMakeFiles\lut-Emu-MC.dir\Emu_MC.cpp.obj -c D:\Desktop\GAMES\Games202\GAMES202_HW\homework4\homework4\lut-gen\Emu_MC.cpp

CMakeFiles/lut-Emu-MC.dir/Emu_MC.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/lut-Emu-MC.dir/Emu_MC.cpp.i"
	D:\Program\MinGW\bin\c++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E D:\Desktop\GAMES\Games202\GAMES202_HW\homework4\homework4\lut-gen\Emu_MC.cpp > CMakeFiles\lut-Emu-MC.dir\Emu_MC.cpp.i

CMakeFiles/lut-Emu-MC.dir/Emu_MC.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/lut-Emu-MC.dir/Emu_MC.cpp.s"
	D:\Program\MinGW\bin\c++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S D:\Desktop\GAMES\Games202\GAMES202_HW\homework4\homework4\lut-gen\Emu_MC.cpp -o CMakeFiles\lut-Emu-MC.dir\Emu_MC.cpp.s

# Object files for target lut-Emu-MC
lut__Emu__MC_OBJECTS = \
"CMakeFiles/lut-Emu-MC.dir/Emu_MC.cpp.obj"

# External object files for target lut-Emu-MC
lut__Emu__MC_EXTERNAL_OBJECTS =

lut-Emu-MC.exe: CMakeFiles/lut-Emu-MC.dir/Emu_MC.cpp.obj
lut-Emu-MC.exe: CMakeFiles/lut-Emu-MC.dir/build.make
lut-Emu-MC.exe: CMakeFiles/lut-Emu-MC.dir/linkLibs.rsp
lut-Emu-MC.exe: CMakeFiles/lut-Emu-MC.dir/objects1.rsp
lut-Emu-MC.exe: CMakeFiles/lut-Emu-MC.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir=D:\Desktop\GAMES\Games202\GAMES202_HW\homework4\homework4\lut-gen\build\CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable lut-Emu-MC.exe"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles\lut-Emu-MC.dir\link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/lut-Emu-MC.dir/build: lut-Emu-MC.exe
.PHONY : CMakeFiles/lut-Emu-MC.dir/build

CMakeFiles/lut-Emu-MC.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles\lut-Emu-MC.dir\cmake_clean.cmake
.PHONY : CMakeFiles/lut-Emu-MC.dir/clean

CMakeFiles/lut-Emu-MC.dir/depend:
	$(CMAKE_COMMAND) -E cmake_depends "MinGW Makefiles" D:\Desktop\GAMES\Games202\GAMES202_HW\homework4\homework4\lut-gen D:\Desktop\GAMES\Games202\GAMES202_HW\homework4\homework4\lut-gen D:\Desktop\GAMES\Games202\GAMES202_HW\homework4\homework4\lut-gen\build D:\Desktop\GAMES\Games202\GAMES202_HW\homework4\homework4\lut-gen\build D:\Desktop\GAMES\Games202\GAMES202_HW\homework4\homework4\lut-gen\build\CMakeFiles\lut-Emu-MC.dir\DependInfo.cmake "--color=$(COLOR)"
.PHONY : CMakeFiles/lut-Emu-MC.dir/depend

