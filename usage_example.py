from triangle_ppy import Delaunay, DebugOutputLevel

# Create a Delaunay object with points
points = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]]

d = Delaunay(points, enable_mesh_indexing=True)

# Set constraints
d.set_quality_constraints(min_angle=20.0, max_area=0.1)

d.set_segment_constraint([[0.0, 0.1], [1.0, 1.0]])
# same as this:
#d.set_segment_constraint(segment_point_indexes = [2, 3])

# Perform triangulation
d.triangulate(quality=True, trace_level=DebugOutputLevel.Info)

# Get results
print(f"Triangle count: {d.triangle_count()}")
print(f"Edge count: {d.edge_count()}")

# Iterate over triangles
faces = d.faces()
for f in faces:
    point1, input_points_idx1 = f.org()
    point2, input_points_idx2 = f.dest()
    point3, input_points_idx3 = f.apex()

    # you can ignore the index!
    pointApx = f.apex()
    pointOrg, _ = f.org()

    # access data
    #  - if index == -1: an added point (i.e. Steiner pt!)
    p1_x = points[input_points_idx1][0] if input_points_idx1 != -1 else point1[0];
    p1_y = points[input_points_idx1][1] if input_points_idx1 != -1 else point1[1];

    p2_x = points[input_points_idx2][0] if input_points_idx2 != -1 else point2[0];
    p2_y = points[input_points_idx2][1] if input_points_idx2 != -1 else point2[1];

    p3_x = points[input_points_idx3][0] if input_points_idx3 != -1 else point3[0];
    p3_y = points[input_points_idx3][1] if input_points_idx3 != -1 else point3[1];

    print(f"Triangle by pts: {point1}/{point2}/{point3} - area: {f.area()}")  
    print(f"Triangle by idx: {[p1_x, p1_y]}/{[p2_x, p2_y]}/{[p3_x, p3_y]}")      
    
# Save to file
d.save_points("./triangulation_result.node")

# Read from file
points = []
if(d.read_points("./triangulation_result.node", points)):
    print(f"Points from file count: {len(points)}")  
    print(f"Points from file: {points}") 
else:
    print(f"Cannot read points from file: {points}") 

