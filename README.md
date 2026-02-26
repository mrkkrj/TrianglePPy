# Triangle-PPy
Python bindings for the Triangle++ library. **WIP !!!!!**

## Usage

To use the binding:
 
 - compile the C++ code with CMake

 - place the resulting binding .pyd file (Windows) / shared library (Linux) in your Python path.

 - make sure your compiler's redistributable shared libraries can be found by Python!

Use it like this in your Python code:

    from triangle_ppy import Delaunay, DebugOutputLevel, AlgorithmType

    points = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]

    d = Delaunay(points)
    d.set_quality_constraints(angle=20.0, area=0.1)

    d.triangulate(quality=True, trace_level=DebugOutputLevel.Info)

    print(f"Triangle count: {d.triangle_count()}")


## Caution

 - first implementation, only basic functionality supported !!!

 - only tested on Windows with Python 3.12 for now !!!
