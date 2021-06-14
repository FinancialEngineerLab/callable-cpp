cmake_minimum_required(VERSION 3.12)
project(beagle)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "-O2")

include_directories(include /usr/local/include/Eigen)
add_subdirectory(beagle)
add_subdirectory(test)
add_subdirectory(src)
