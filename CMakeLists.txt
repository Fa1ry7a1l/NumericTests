cmake_minimum_required(VERSION 3.30)
project(tests CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

add_executable(tests main.cpp)

target_compile_options(tests PRIVATE -Wall -Wextra)
target_compile_options(tests PRIVATE $<$<CONFIG:Release>:-O3>)
