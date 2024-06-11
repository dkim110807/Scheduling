cmake_minimum_required(VERSION 3.26)
project(Scheduling)

set(CMAKE_CXX_STANDARD 20)

add_executable(Scheduling main.cpp include/debug.h)

target_compile_definitions(Scheduling PRIVATE LOCAL)