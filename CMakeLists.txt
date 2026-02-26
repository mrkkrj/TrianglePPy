cmake_minimum_required(VERSION 3.15...4.2)
project(triangle_ppy LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)

# binding tools
find_package(Python 3.8 COMPONENTS Interpreter Development REQUIRED)
add_subdirectory(extern/pybind11)

# TrianglePP C++ library
add_subdirectory(extern/TrianglePP)

# Python bindings for it
pybind11_add_module(triangle_ppy triangle_pp_bindings.cpp)

target_include_directories(${PROJECT_NAME} PUBLIC extern/TrianglePP/source)
target_link_libraries(${PROJECT_NAME} PUBLIC TrianglePP)
