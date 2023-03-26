cmake_minimum_required(VERSION 3.10.2)
project(fixed_position_Turing_pattern)

set(CMAKE_CXX_STANDARD 14)
add_executable(fixed_position_Turing_pattern main.cpp Brusselator.cpp Brusselator_system.cpp Brusselator.h ran2.cpp)