from triangle_ppy import Delaunay, DebugOutputLevel
import sys

output_level = DebugOutputLevel.Info  #.Nothing #.Debug #.Info 

points = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]]
d = Delaunay(points, enable_mesh_indexing=True)

# enable file tracing
#  -- only on Windows!!! --> REMOVE??? -- not working !!!
d.enable_file_io_trace(True)

# constrain Delaunay
d.set_quality_constraints(min_angle=20.0, max_area=0.1)

ok = d.check_constraints()
if not ok:
    print(f"Error - triangulation not guaranteed to succeed with given constraints!")
elif output_level != DebugOutputLevel.Nothing:
    print(f"Quality constraints OK, triangulation guaranteed to succeed")

segment_end_pts = [[0.0, 0.1], [1.0, 1.0]]
d.set_segment_constraint(segment_end_pts)

# Perform triangulation
d.triangulate(quality=True, trace_level=output_level)

triangle_count = d.triangle_count()
edge_count = d.edge_count()
vertice_count = d.vertice_count()

if output_level != DebugOutputLevel.Nothing:
    print(f"Triangle count: {d.triangle_count()}")
    print(f"Edge count: {d.edge_count()}")
    print(f"Vertice count: {d.vertice_count()}")

print(f" >> Triangulation OK ... \n")    



print(f"Finished ---> \n")
