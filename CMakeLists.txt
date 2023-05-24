cmake_minimum_required(VERSION 3.22)
project(hav_cpp_from_python_2)

set(CMAKE_CXX_STANDARD 14)

include_directories(lib)

add_executable(hav_cpp_from_python_2
        ekf.cpp
        ekf.h
        geo_ned.cpp
        geo_ned.h
        main_utm_2.cpp utm.h utm.cpp )
