from triangle_ppy import Delaunay, DebugOutputLevel
import sys

output_level = DebugOutputLevel.Info  #.Nothing #.Debug #.Info 

points = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]]
d = Delaunay(points, enable_mesh_indexing=True)

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

# Write segments out
try:
    d.save_segments("./rw_segments_test.poly")
    if output_level != DebugOutputLevel.Nothing:
        print(f"Points & segments saved to a .poly file") 
except:
    print(f"Error - cannot save points & segments to file!") 

# Read segments
duplicate_segments = 0
points = []
segment_endpoints = []
hole_markers = []
region_constr = []

try:
    duplicate_segments = d.read_segments(
                            "./rw_segments_test.poly",
                            points,
                            segment_endpoints,
                            hole_markers,
                            region_constr,
                            trace_level=output_level
                         )
    if output_level != DebugOutputLevel.Nothing:
        print(f"Points & segments read from a .poly file")
except:
    print(f"Error - cannot read points & segments from file!")
    sys.exit(1)

print(f" >> File reading OK ... \n")    

# Show results
print(f"Read points: {points}") 
print(f" -- point count: {len(points)}") 
if len(points) != vertice_count:
    print(f"Error - not enough points read from file!")

print(f"Read segment endpoints: {segment_endpoints}") 
print(f" -- segment count: {len(segment_endpoints)/2:.0f}") 
print(f" -- duplicate segments: {duplicate_segments}") 
if len(segment_endpoints) != segment_end_pts:
    print(f"Error - not enough segments read from file!")

print(f"Read hole markers: {hole_markers}") 
print(f" -- hole count: {len(hole_markers)}") 

print(f"Read region constraints: {region_constr}") 
print(f" -- region constraint count: {len(region_constr)}") 

print(f" >> Read segments OK ... \n")


# --- more ????

# OPEN TODO:: add holes and region constraints !!!


print(f"Finished ---> \n")
