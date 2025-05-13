# TrianglePPy
Python bindings for the Triangle++ library. **WIP !!!!!**

## Compilation Instructions

To compile the binding:
 - Ensure *pybind11* is installed (pip install *pybind11* or include it in your project).

 - Compile the C++ code with a command like:
    c++ -O3 -Wall -shared -std=c++17 -fPIC $(python3 -m pybind11 --includes) delaunay_binding.cpp -o delaunay.so

    Adjust the include paths and library links for your C++ Delaunay implementation.

 - Place the resulting delaunay.so in your Python path.


