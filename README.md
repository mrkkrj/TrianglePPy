# Triangle-PPy
<img src="triangle-PPy-sm.jpg" alt="triangle-PP's logo" width="160"/><br/>Python bindings for the Triangle++ library. **WIP !!!!!**

## Usage

To use the bindings:
 
 - install Python version >= 3.8

 - compile the C++ code with CMake

 - place the resulting bindings .pyd file (Windows) / shared library (Linux) in your Python path.

 - also make sure your compiler's redistributable shared libraries can be found by Python!

Use it like this in your Python code:

    from triangle_ppy import Delaunay, DebugOutputLevel

    points = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]

    d = Delaunay(points)
    d.set_quality_constraints(min_angle=20.0, max_area=0.1)

    d.triangulate(quality=True, trace_level=DebugOutputLevel.Info)

    print(f"Triangle count: {d.triangle_count()}")
    
    faces = d.faces() # faces are oriented triangles!

    for f in faces:
        point1 = f.org()
        point2 = f.dest()
        point3 = f.apex()

        print(f"Triangle: {point1}/{point2}/{point3}")  

For more usage cases have a look into the *usage_example.py* file.

## Caution

 - first implementation, only basic functionality supported !!!

 - only tested on Windows with Python 3.12 for now !!!

 - WIP !!!
