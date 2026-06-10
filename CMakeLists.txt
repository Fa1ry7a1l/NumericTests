cmake_minimum_required(VERSION 3.30)
project(tests CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

find_package(Threads REQUIRED)

add_executable(tests main.cpp)

target_include_directories(tests PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/third_party/ntl-install/include"
    "${CMAKE_CURRENT_SOURCE_DIR}/third_party/gmp-install/include"
)

target_link_libraries(tests PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/third_party/ntl-install/lib/libntl.a"
    "${CMAKE_CURRENT_SOURCE_DIR}/third_party/gmp-install/lib/libgmp.a"
    Threads::Threads
)

target_compile_options(tests PRIVATE -Wall -Wextra)
target_compile_options(tests PRIVATE $<$<CONFIG:Release>:-O3>)
